#ifndef KESSLER_FUNCTIONS_IMPL_HPP
#define KESSLER_FUNCTIONS_IMPL_HPP

#include "kessler_functions.hpp"

#include <Kokkos_Core.hpp>

#include <ekat_assert.hpp>

#include <limits>
#include <cmath>
#include <string>
#include <type_traits>

namespace scream {
namespace kessler {

/*
 * kessler_run
 *
 * Implements the Kessler (1969) warm rain scheme as described by
 * Soong & Ogura (1973) and Klemp & Wilhelmson (1978).  The scheme
 * sub-cycles through the moisture processes to satisfy the CFL
 * condition for sedimentation.
 *
 * Strategy:
 *  - Scratch arrays live in a caller-owned KesslerWorkspace, allocated
 *    once and reused across calls (no device allocation per call).
 *  - cpair and rair are 1D (composition-independent) scalar constants,
 *    so per-(col,level) derived quantities that only depend on them
 *    (xk = cp/R, f5 = 4093*lv/cp) collapse to single host-side scalars
 *    captured by value, instead of 2D views read/recomputed every thread.
 *  - Each sub-cycle iteration launches 4 kernels (down from 8): steps
 *    that only differ by a branch on the level index (precl accum +
 *    sedimentation; elapsed-time/mask update + terminal fall speed) are
 *    fused into a single kernel over ncols*nk, and the per-column CFL
 *    min-reduction (previously a TeamPolicy reduction kernel) is folded
 *    into the neighboring per-level kernel via Kokkos::atomic_min.
 *  - Three of those four kernels (precl_sed, micro, dt0_recompute) exit
 *    immediately for a column whose mask is already 0 (converged as of
 *    the previous iteration), instead of doing a full no-op pass over
 *    it -- this makes columns that finish early (e.g. dry columns) cost
 *    almost nothing on later iterations while rain-heavy columns keep
 *    sub-cycling. The 4th (time_velqr) is deliberately NOT guarded this
 *    way: it both writes and would need to read mask for the same
 *    column within one kernel launch, which Kokkos::parallel_for does
 *    not order between indices -- see the note at its call site.
 *  - The sub-cycling while loop runs on the host; convergence is only
 *    checked (device-to-host sync) every CHECK_INTERVAL iterations
 *    instead of every single one. Combined with the mask-guarded early
 *    exits above, the extra iterations this can run past actual
 *    convergence are true no-ops (columns already at mask==0 do
 *    essentially zero work), so this only changes how often the host
 *    peeks at device state, not what gets computed.
 *  - Convention: lyr_surf and lyr_toa are 0-based C++ indices.
 *    lyr_step = +1 if lyr_surf <= lyr_toa, -1 otherwise.
 *    The driver sets lyr_surf=0, lyr_toa=nz-1 (surface at k=0).
 *  - Before sub-cycling, a column whose CFL-limited initial dt0 has
 *    collapsed to (near) zero triggers a fatal EKAT_REQUIRE_MSG (mirrors
 *    the Fortran "bad time splitting" check). This throws a normal C++
 *    exception that unwinds through run_impl and the atmosphere driver,
 *    letting output streams get finalized/closed before the run aborts,
 *    rather than silently spinning through an enormous number of
 *    sub-cycle iterations.
 */
template <typename S, typename D>
void KesslerFunctions<S,D>::kessler_run(
  const int ncols,
  const int nz,
  const Scalar dt,
  const int lyr_surf,
  const int lyr_toa,
  const KesslerData& kd,
  KesslerWorkspace& workspace,
  const view_2d<const Scalar>& rho,
  const view_2d<const Scalar>& z,
  const view_2d<const Scalar>& pk,
  const view_2d<Scalar>& theta,
  const view_2d<Scalar>& qv,
  const view_2d<Scalar>& qc,
  const view_2d<Scalar>& qr,
  const view_1d<Scalar>& precl,
  const view_2d<Scalar>& relhum)
{
  // Physical constants (captured by value in lambdas). cpair/rair are
  // now scalar constants, so xk and f5 -- both derived solely from them
  // -- are single values computed once here rather than per-thread.
  const Scalar lv    = kd.lv;
  const Scalar pref  = kd.pref;    // already in hPa
  const Scalar rhoqr = kd.rhoqr;
  const Scalar cpair = kd.cpair;
  const Scalar xk    = kd.cpair / kd.rair;       // cp/R
  const Scalar f5    = Scalar(4093) * lv / cpair; // condensation-rate constant
  const Scalar f2x   = Scalar(17.27);

  // Vertical direction
  const int lyr_step = (lyr_surf <= lyr_toa) ? 1 : -1;

  // Level range (0-based inclusive)
  const int kmin = (lyr_step > 0) ? lyr_surf : lyr_toa;
  const int kmax = (lyr_step > 0) ? lyr_toa  : lyr_surf;
  const int nk   = kmax - kmin + 1; // nz

  // ---------------------------------------------------------------
  // Scratch views: sized (no-op if already sized) and owned by the
  // caller so no device allocation happens on this call.
  // ---------------------------------------------------------------
  workspace.init(ncols, nz);
  const auto& r            = workspace.r;
  const auto& rhalf        = workspace.rhalf;
  const auto& velqr        = workspace.velqr;
  const auto& sed          = workspace.sed;
  const auto& pc           = workspace.pc;
  const auto& dt0          = workspace.dt0;
  const auto& mask         = workspace.mask;
  const auto& time_counter = workspace.time_counter;
  const auto& precl_acc    = workspace.precl_acc;

  // Layout-aware flattening of the (col,level) loop: pick whichever
  // traversal keeps consecutive thread indices contiguous in memory for
  // the view's actual layout (col-fastest for LayoutLeft, as Kokkos
  // defaults to on GPU; level-fastest for LayoutRight, the CPU default).
  const bool col_fastest =
    std::is_same<typename view_2d<Scalar>::array_layout, Kokkos::LayoutLeft>::value;

  // ---------------------------------------------------------------
  // Kernel 1: reset per-column bookkeeping scalars (ncols-sized)
  // ---------------------------------------------------------------
  Kokkos::parallel_for("kessler_init0",
    Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols),
    KOKKOS_LAMBDA(const int col) {
      dt0(col)          = dt;
      mask(col)         = Scalar(1);
      time_counter(col) = Scalar(0);
      precl_acc(col)    = Scalar(0);
      precl(col)        = Scalar(0);
    });

  // ---------------------------------------------------------------
  // Kernel 2: derived constants + terminal fall speed, with the
  // initial CFL min-reduction for dt0 folded in via atomic_min
  // (replaces a separate TeamPolicy reduction kernel).
  // ---------------------------------------------------------------
  Kokkos::parallel_for("kessler_init_fields",
    Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nk),
    KOKKOS_LAMBDA(const int idx) {
      int col, k;
      if (col_fastest) { col = idx % ncols; k = kmin + idx / ncols; }
      else              { col = idx / nk;    k = kmin + idx % nk; }

      r(col, k)     = Scalar(0.001) * rho(col, k);          // g/cm^3
      rhalf(col, k) = std::sqrt(rho(col, lyr_surf) / rho(col, k));
      pc(col, k)    = Scalar(3.8) / (std::pow(pk(col, k), xk) * pref);
      // Guard against round-off negative qr before computing velqr
      qr(col, k)    = std::max(qr(col, k), Scalar(0));
      velqr(col, k) = Scalar(36.34) * rhalf(col, k) *
                      std::pow(qr(col, k) * r(col, k), Scalar(0.1364));

      if (k != lyr_toa && std::abs(velqr(col, k)) > Scalar(1e-12)) {
        const Scalar dz = z(col, k + lyr_step) - z(col, k);
        Kokkos::atomic_min(&dt0(col), Scalar(0.8) * dz / velqr(col, k));
      }
    });

  // ---------------------------------------------------------------
  // Guard against a pathologically small CFL-limited initial sub-cycle
  // step (mirrors the Fortran "bad time splitting" check, which errors
  // out here rather than entering the sub-cycle loop). Left unchecked,
  // a column with dt0 this small would need an enormous number of
  // sub-cycle iterations to reach dt, effectively hanging the run
  // instead of failing cleanly. This is the one host sync kessler_run
  // needs before the loop, so it's also the last point where a fence
  // would otherwise have been needed.
  // ---------------------------------------------------------------
  {
    using MinLocReducer = Kokkos::MinLoc<Scalar, int>;
    typename MinLocReducer::value_type minloc_result;
    Kokkos::parallel_reduce("kessler_check_dt0",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols),
      KOKKOS_LAMBDA(const int col, typename MinLocReducer::value_type& lminloc) {
        if (dt0(col) < lminloc.val) {
          lminloc.val = dt0(col);
          lminloc.loc = col;
        }
      }, MinLocReducer(minloc_result));

    EKAT_REQUIRE_MSG(minloc_result.val >= Scalar(1e-12),
      "Error! Kessler: bad time splitting -- CFL-limited sub-cycle step "
      "collapsed to (near) zero.\n"
      "  - physics dt:        " + std::to_string(dt) + "\n"
      "  - min sub-cycle dt0: " + std::to_string(minloc_result.val) + "\n"
      "  - local column index: " + std::to_string(minloc_result.loc) +
      " (of " + std::to_string(ncols) + ")\n");
  }

  // ---------------------------------------------------------------
  // Sub-cycling while loop
  // ---------------------------------------------------------------
  constexpr int CHECK_INTERVAL = 4; // convergence sync frequency
  bool all_converged = false;
  int iter = 0;
  while (!all_converged) {
    ++iter;

    // -- Steps 1+2+3 fused: precip accumulation + sedimentation --
    Kokkos::parallel_for("kessler_precl_sed",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nk),
      KOKKOS_LAMBDA(const int idx) {
        int col, k;
        if (col_fastest) { col = idx % ncols; k = kmin + idx / ncols; }
        else              { col = idx / nk;    k = kmin + idx % nk; }

        // Column already converged as of the previous iteration: skip
        // entirely rather than doing a no-op update. mask(col) is only
        // ever written (by kessler_time_velqr) below, so this reflects
        // the state at the end of the previous iteration for the whole
        // of this one.
        if (mask(col) == Scalar(0)) return;

        if (k == lyr_surf) {
          const Scalar p_val = rho(col, lyr_surf) * qr(col, lyr_surf) *
                               velqr(col, lyr_surf) / rhoqr;
          precl(col)     = p_val;
          precl_acc(col) += p_val * dt0(col);
        }

        if (k == lyr_toa) {
          const int kbelow = lyr_toa - lyr_step;
          sed(col, lyr_toa) = -dt0(col) * qr(col, lyr_toa) *
            velqr(col, lyr_toa) /
            (Scalar(0.5) * (z(col, lyr_toa) - z(col, kbelow)));
        } else {
          const int kup = k + lyr_step;
          sed(col, k) = dt0(col) *
            ((r(col, kup) * qr(col, kup) * velqr(col, kup)) -
             (r(col, k)   * qr(col, k)   * velqr(col, k))) /
            (r(col, k) * (z(col, kup) - z(col, k)));
        }
      });

    // -- Step 4: microphysics adjustments (autoconversion, collection,
    //            evaporation, saturation adjustment) --
    Kokkos::parallel_for("kessler_micro",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nk),
      KOKKOS_LAMBDA(const int idx) {
        int col, k;
        if (col_fastest) { col = idx % ncols; k = kmin + idx / ncols; }
        else              { col = idx / nk;    k = kmin + idx % nk; }

        if (mask(col) == Scalar(0)) return;

        const Scalar dt0c = dt0(col);

        // Autoconversion + collection (Klemp & Wilhelmson 1978, Eq 2.13a,b)
        // Semi-implicit treatment of the collection term.
        const Scalar qrprod =
          qc(col, k) -
          (qc(col, k) - dt0c *
            std::max(Scalar(0.001) * (qc(col, k) - Scalar(0.001)),
                     Scalar(0))) /
          (Scalar(1) + dt0c * Scalar(2.2) *
           std::pow(qr(col, k), Scalar(0.875)));
        qc(col, k) = std::max(qc(col, k) - qrprod, Scalar(0));
        qr(col, k) = std::max(qr(col, k) + qrprod + sed(col, k), Scalar(0));

        // Saturation mixing ratio via Teten's formula
        // (Klemp & Wilhelmson 1978, Eq 2.11)
        const Scalar T_pi = pk(col, k) * theta(col, k); // temperature = pi*theta
        const Scalar qvs  = pc(col, k) *
          std::exp(f2x * (T_pi - Scalar(273)) / (T_pi - Scalar(36)));

        // Condensation rate (Durran & Klemp 1983, Eq A13-A14)
        const Scalar dT_pi_36 = T_pi - Scalar(36);
        const Scalar prod = (qv(col, k) - qvs) /
          (Scalar(1) + qvs * f5 / (dT_pi_36 * dT_pi_36));

        // Evaporation rate (Klemp & Wilhelmson 1978, Eq 2.14a,b)
        // dim(qvs, qv) = max(qvs - qv, 0): evaporation only in subsaturated air.
        // rqr^0.2046 and rqr^0.525 share the base rqr; compute log(rqr)
        // once and derive both powers via exp() instead of two
        // independent std::pow calls (each an internal log/exp pair).
        // Guarded against rqr == 0 (common, e.g. dry columns) to avoid a
        // log(0) trap under debug/fpe builds.
        const Scalar rqr = r(col, k) * qr(col, k);
        Scalar rqr_2046 = Scalar(0);
        Scalar rqr_525  = Scalar(0);
        if (rqr > Scalar(0)) {
          const Scalar log_rqr = std::log(rqr);
          rqr_2046 = std::exp(Scalar(0.2046) * log_rqr);
          rqr_525  = std::exp(Scalar(0.525)  * log_rqr);
        }
        const Scalar A = dt0c *
          ((Scalar(1.6) + Scalar(124.9) * rqr_2046) *
           rqr_525 /
           (Scalar(2550000) * pc(col, k) / (Scalar(3.8) * qvs) +
            Scalar(540000))) *
          (std::max(qvs - qv(col, k), Scalar(0)) / (r(col, k) * qvs));
        const Scalar B   = std::max(-prod - qc(col, k), Scalar(0));
        const Scalar ern = std::min(A, std::min(B, qr(col, k)));

        // Saturation adjustment (Durran & Klemp 1983, Eq A1-A4).
        // The early mask(col)==0 return above means every thread that
        // reaches this point is on a still-active column, so the
        // previous mask-blended (msk * new + (1-msk) * old) form
        // simplifies to just applying the update unconditionally.
        const Scalar prod_adj = std::max(prod, -qc(col, k));
        theta(col, k) += lv / (cpair * pk(col, k)) * (prod_adj - ern);
        qv(col, k) = std::max(qv(col, k) - prod_adj + ern, Scalar(0));
        qc(col, k) += prod_adj;
        qr(col, k) = std::max(qr(col, k) - ern, Scalar(0));
      });

    // -- Steps 5+6 fused: elapsed time / mask update (once per column,
    //    guarded on k == kmin) + terminal fall speed recompute (all k) --
    // NOTE: deliberately NOT early-returning on mask(col)==0 here, unlike
    // the other three sub-cycle kernels. This kernel both writes mask
    // (from the k==kmin thread) and would need to read it (for the
    // early-return guard) from every thread for the same column -- with
    // no ordering guarantee between different indices in one
    // Kokkos::parallel_for, that read/write combination would be a data
    // race. The mask(col) multiplier below keeps the update an exact
    // no-op for already-converged columns instead.
    Kokkos::parallel_for("kessler_time_velqr",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nk),
      KOKKOS_LAMBDA(const int idx) {
        int col, k;
        if (col_fastest) { col = idx % ncols; k = kmin + idx / ncols; }
        else              { col = idx / nk;    k = kmin + idx % nk; }

        if (k == kmin) {
          time_counter(col) += mask(col) * dt0(col);
          dt0(col) = std::max(dt - time_counter(col), Scalar(0));
          mask(col) = (std::abs(dt - time_counter(col)) > Scalar(1e-5))
                        ? Scalar(1) : Scalar(0);
        }

        velqr(col, k) = Scalar(36.34) * rhalf(col, k) *
                        std::pow(qr(col, k) * r(col, k), Scalar(0.1364));
      });

    // -- Step 7: recompute dt0 (CFL constraint) via atomic_min --
    // (replaces a separate TeamPolicy reduction kernel; velqr(col,k) was
    // just written by this same thread above, so there is no race).
    Kokkos::parallel_for("kessler_dt0_recompute",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nk),
      KOKKOS_LAMBDA(const int idx) {
        int col, k;
        if (col_fastest) { col = idx % ncols; k = kmin + idx / ncols; }
        else              { col = idx / nk;    k = kmin + idx % nk; }

        if (mask(col) == Scalar(0)) return;

        if (k != lyr_toa && std::abs(velqr(col, k)) > Scalar(1e-12)) {
          const Scalar dz = z(col, k + lyr_step) - z(col, k);
          Kokkos::atomic_min(&dt0(col), Scalar(0.8) * dz / velqr(col, k));
        }
      });

    // -- Step 8: check convergence, only every CHECK_INTERVAL iterations
    //    to reduce host-device sync frequency. Kernels above still run
    //    for every column on every iteration regardless, so this does
    //    not change what is computed -- only how often the host checks. --
    if (iter % CHECK_INTERVAL == 0) {
      int n_active = 0;
      Kokkos::parallel_reduce("kessler_converge",
        Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols),
        KOKKOS_LAMBDA(const int col, int& cnt) {
          if (mask(col) != Scalar(0)) ++cnt;
        }, n_active);
      all_converged = (n_active == 0);
    }

  } // while (!all_converged)

  // ---------------------------------------------------------------
  // Finalise: average precipitation rate over the full time step
  // ---------------------------------------------------------------
  Kokkos::parallel_for("kessler_precl_final",
    Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols),
    KOKKOS_LAMBDA(const int col) {
      precl(col) = precl_acc(col) / dt;
    });

  // ---------------------------------------------------------------
  // Diagnostic: relative humidity
  // ---------------------------------------------------------------
  Kokkos::parallel_for("kessler_relhum",
    Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nk),
    KOKKOS_LAMBDA(const int idx) {
      int col, k;
      if (col_fastest) { col = idx % ncols; k = kmin + idx / ncols; }
      else              { col = idx / nk;    k = kmin + idx % nk; }

      const Scalar T_pi = pk(col, k) * theta(col, k);
      const Scalar qvs  = pc(col, k) *
        std::exp(f2x * (T_pi - Scalar(273)) / (T_pi - Scalar(36)));
      relhum(col, k) = qv(col, k) / qvs * Scalar(100);
    });

  Kokkos::fence();
}

// -----------------------------------------------------------------------

template <typename S, typename D>
void KesslerFunctions<S,D>::kessler_update_timestep_init(
  const int ncols,
  const int nz,
  const view_2d<const Scalar>& temp,
  const view_2d<Scalar>& temp_prev,
  const view_2d<Scalar>& ttend_t)
{
  Kokkos::parallel_for("kessler_ts_init",
    Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nz),
    KOKKOS_LAMBDA(const int idx) {
      const int col = idx / nz;
      const int k   = idx % nz;
      temp_prev(col, k) = temp(col, k);
      ttend_t(col, k)   = Scalar(0);
    });
  Kokkos::fence();
}

// -----------------------------------------------------------------------

template <typename S, typename D>
void KesslerFunctions<S,D>::kessler_update_run(
  const int ncols,
  const int nz,
  const Scalar dt,
  const view_2d<const Scalar>& theta,
  const view_2d<const Scalar>& exner,
  const view_2d<const Scalar>& temp_prev,
  const view_2d<Scalar>& ttend_t)
{
  Kokkos::parallel_for("kessler_update_run",
    Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nz),
    KOKKOS_LAMBDA(const int idx) {
      const int col = idx / nz;
      const int k   = idx % nz;
      ttend_t(col, k) += (theta(col, k) * exner(col, k) - temp_prev(col, k)) / dt;
    });
  Kokkos::fence();
}

// -----------------------------------------------------------------------

template <typename S, typename D>
void KesslerFunctions<S,D>::kessler_update_timestep_final(
  const int ncols,
  const int nz,
  const Scalar gravit,
  const Scalar cpair,
  const view_2d<const Scalar>& temp,
  const view_2d<const Scalar>& zm,
  const view_1d<const Scalar>& phis,
  const view_2d<Scalar>& st_energy)
{
  Kokkos::parallel_for("kessler_ts_final",
    Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nz),
    KOKKOS_LAMBDA(const int idx) {
      const int col = idx / nz;
      const int k   = idx % nz;
      st_energy(col, k) = cpair * temp(col, k) +
                          gravit * zm(col, k) + phis(col);
    });
  Kokkos::fence();
}

} // namespace kessler
} // namespace scream

#endif // KESSLER_FUNCTIONS_IMPL_HPP
