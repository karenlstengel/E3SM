#ifndef KESSLER_FUNCTIONS_IMPL_HPP
#define KESSLER_FUNCTIONS_IMPL_HPP

#include "kessler_functions.hpp"

#include <Kokkos_Core.hpp>

#include <limits>
#include <cmath>

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
 *  - Scratch arrays (r, rhalf, velqr, sed, pc, f5, dt0, mask,
 *    time_counter, precl_acc) are allocated as Kokkos views here
 *    and freed on return.
 *  - The sub-cycling while loop runs on the host; each iteration
 *    dispatches device kernels and then checks convergence via a
 *    device-to-host reduction.
 *  - Convention: lyr_surf and lyr_toa are 0-based C++ indices.
 *    lyr_step = +1 if lyr_surf <= lyr_toa, -1 otherwise.
 *    The driver sets lyr_surf=0, lyr_toa=nz-1 (surface at k=0).
 */
template <typename S, typename D>
void KesslerFunctions<S,D>::kessler_run(
  const int ncols,
  const int nz,
  const Scalar dt,
  const int lyr_surf,
  const int lyr_toa,
  const KesslerData& kd,
  const view_2d<const Scalar>& cpair,
  const view_2d<const Scalar>& rair,
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
  // Physical constants (captured by value in lambdas)
  const Scalar lv    = kd.lv;
  const Scalar pref  = kd.pref;    // already in hPa
  const Scalar rhoqr = kd.rhoqr;
  const Scalar f2x   = Scalar(17.27);

  // Vertical direction
  const int lyr_step = (lyr_surf <= lyr_toa) ? 1 : -1;

  // Level ranges (0-based inclusive)
  const int kmin = (lyr_step > 0) ? lyr_surf : lyr_toa;
  const int kmax = (lyr_step > 0) ? lyr_toa  : lyr_surf;
  const int nk   = kmax - kmin + 1; // nz

  // Range for the sedimentation inner loop: all levels except lyr_toa
  const int kmin_sed = (lyr_step > 0) ? lyr_surf       : lyr_toa + 1;
  const int kmax_sed = (lyr_step > 0) ? lyr_toa - 1    : lyr_surf;
  const int nk_sed   = kmax_sed - kmin_sed + 1; // nz-1

  // ---------------------------------------------------------------
  // Allocate scratch views
  // ---------------------------------------------------------------
  view_2d<Scalar> r            ("kessler_r",            ncols, nz);
  view_2d<Scalar> rhalf        ("kessler_rhalf",        ncols, nz);
  view_2d<Scalar> velqr        ("kessler_velqr",        ncols, nz);
  view_2d<Scalar> sed          ("kessler_sed",          ncols, nz);
  view_2d<Scalar> pc           ("kessler_pc",           ncols, nz);
  view_2d<Scalar> f5           ("kessler_f5",           ncols, nz);
  view_1d<Scalar> dt0          ("kessler_dt0",          ncols);
  view_1d<Scalar> mask         ("kessler_mask",         ncols);
  view_1d<Scalar> time_counter ("kessler_time_counter", ncols);
  view_1d<Scalar> precl_acc    ("kessler_precl_acc",    ncols);

  // ---------------------------------------------------------------
  // Kernel 1: initialise derived constants and terminal fall speed
  // ---------------------------------------------------------------
  Kokkos::parallel_for("kessler_init_fields",
    Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nk),
    KOKKOS_LAMBDA(const int idx) {
      const int col = idx / nk;
      const int k   = kmin + (idx % nk);
      const Scalar xk = cpair(col, k) / rair(col, k);      // cp/R
      f5(col, k)    = Scalar(4093) * lv / cpair(col, k);
      r(col, k)     = Scalar(0.001) * rho(col, k);          // g/cm^3
      rhalf(col, k) = std::sqrt(rho(col, lyr_surf) / rho(col, k));
      pc(col, k)    = Scalar(3.8) / (std::pow(pk(col, k), xk) * pref);
      // Guard against round-off negative qr before computing velqr
      qr(col, k)    = std::max(qr(col, k), Scalar(0));
      velqr(col, k) = Scalar(36.34) * rhalf(col, k) *
                      std::pow(qr(col, k) * r(col, k), Scalar(0.1364));
    });

  // ---------------------------------------------------------------
  // Kernel 2: compute initial sub-cycling time step dt0 via CFL
  //           (min reduction over levels within each column)
  //           and initialise per-column bookkeeping scalars.
  // ---------------------------------------------------------------
  using TeamPol    = Kokkos::TeamPolicy<typename KT::ExeSpace>;
  using MemberType = typename TeamPol::member_type;

  Kokkos::parallel_for("kessler_init_dt0",
    TeamPol(ncols, Kokkos::AUTO),
    KOKKOS_LAMBDA(const MemberType& team) {
      const int col = team.league_rank();

      // Min CFL dt across all interior levels
      Scalar dtmin;
      Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, nk_sed),
        [&](const int idx, Scalar& lmin) {
          const int k = kmin_sed + idx;
          if (std::abs(velqr(col, k)) > Scalar(1e-12)) {
            const Scalar dz = z(col, k + lyr_step) - z(col, k);
            lmin = std::min(lmin, Scalar(0.8) * dz / velqr(col, k));
          }
        },
        Kokkos::Min<Scalar>(dtmin));

      Kokkos::single(Kokkos::PerTeam(team), [&]() {
        dt0(col)          = std::min(dt, dtmin);
        mask(col)         = Scalar(1);
        time_counter(col) = Scalar(0);
        precl_acc(col)    = Scalar(0);
        precl(col)        = Scalar(0);
      });
    });

  Kokkos::fence();

  // ---------------------------------------------------------------
  // Sub-cycling while loop
  // Host drives the loop; convergence checked via parallel_reduce.
  // ---------------------------------------------------------------
  bool all_converged = false;
  while (!all_converged) {

    // -- Step 1: accumulate precipitation (weighted by sub-cycle dt) --
    Kokkos::parallel_for("kessler_precl_accum",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols),
      KOKKOS_LAMBDA(const int col) {
        const Scalar p_val = rho(col, lyr_surf) * qr(col, lyr_surf) *
                             velqr(col, lyr_surf) / rhoqr;
        precl(col)     = p_val;
        precl_acc(col) += mask(col) * p_val * dt0(col);
      });

    // -- Step 2: sedimentation for all levels except lyr_toa --
    Kokkos::parallel_for("kessler_sed_inner",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nk_sed),
      KOKKOS_LAMBDA(const int idx) {
        const int col = idx / nk_sed;
        const int k   = kmin_sed + (idx % nk_sed);
        const int kup = k + lyr_step;
        sed(col, k) = dt0(col) *
          ((r(col, kup) * qr(col, kup) * velqr(col, kup)) -
           (r(col, k)   * qr(col, k)   * velqr(col, k))) /
          (r(col, k) * (z(col, kup) - z(col, k)));
      });

    // -- Step 3: sedimentation at lyr_toa (no flux from above) --
    Kokkos::parallel_for("kessler_sed_toa",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols),
      KOKKOS_LAMBDA(const int col) {
        const int kbelow = lyr_toa - lyr_step;
        sed(col, lyr_toa) = -dt0(col) * qr(col, lyr_toa) *
          velqr(col, lyr_toa) /
          (Scalar(0.5) * (z(col, lyr_toa) - z(col, kbelow)));
      });

    // -- Step 4: microphysics adjustments (autoconversion, collection,
    //            evaporation, saturation adjustment) --
    Kokkos::parallel_for("kessler_micro",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nk),
      KOKKOS_LAMBDA(const int idx) {
        const int col = idx / nk;
        const int k   = kmin + (idx % nk);
        const Scalar msk  = mask(col);
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
          (Scalar(1) + qvs * f5(col, k) / (dT_pi_36 * dT_pi_36));

        // Evaporation rate (Klemp & Wilhelmson 1978, Eq 2.14a,b)
        // dim(qvs, qv) = max(qvs - qv, 0): evaporation only in subsaturated air
        const Scalar rqr = r(col, k) * qr(col, k);
        const Scalar A = dt0c *
          ((Scalar(1.6) + Scalar(124.9) * std::pow(rqr, Scalar(0.2046))) *
           std::pow(rqr, Scalar(0.525)) /
           (Scalar(2550000) * pc(col, k) / (Scalar(3.8) * qvs) +
            Scalar(540000))) *
          (std::max(qvs - qv(col, k), Scalar(0)) / (r(col, k) * qvs));
        const Scalar B   = std::max(-prod - qc(col, k), Scalar(0));
        const Scalar ern = std::min(A, std::min(B, qr(col, k)));

        // Saturation adjustment (Durran & Klemp 1983, Eq A1-A4)
        // Applied only to non-converged columns via mask.
        const Scalar prod_adj = std::max(prod, -qc(col, k));
        theta(col, k) += msk * (lv / (cpair(col, k) * pk(col, k)) *
                                (prod_adj - ern));
        qv(col, k) = msk * std::max(qv(col, k) - prod_adj + ern, Scalar(0))
                   + (Scalar(1) - msk) * qv(col, k);
        qc(col, k) += msk * prod_adj;
        qr(col, k) = msk * std::max(qr(col, k) - ern, Scalar(0))
                   + (Scalar(1) - msk) * qr(col, k);
      });

    // -- Step 5: advance elapsed time and update mask / dt0 --
    Kokkos::parallel_for("kessler_update_time",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols),
      KOKKOS_LAMBDA(const int col) {
        time_counter(col) += mask(col) * dt0(col);
        dt0(col) = std::max(dt - time_counter(col), Scalar(0));
        mask(col) = (std::abs(dt - time_counter(col)) > Scalar(1e-5))
                      ? Scalar(1) : Scalar(0);
      });

    // -- Step 6: recompute terminal fall speed with updated qr --
    Kokkos::parallel_for("kessler_velqr_update",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nk),
      KOKKOS_LAMBDA(const int idx) {
        const int col = idx / nk;
        const int k   = kmin + (idx % nk);
        velqr(col, k) = Scalar(36.34) * rhalf(col, k) *
                        std::pow(qr(col, k) * r(col, k), Scalar(0.1364));
      });

    // -- Step 7: recompute dt0 with updated velqr (CFL constraint) --
    Kokkos::parallel_for("kessler_recompute_dt0",
      TeamPol(ncols, Kokkos::AUTO),
      KOKKOS_LAMBDA(const MemberType& team) {
        const int col = team.league_rank();
        Scalar dtmin;
        Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(team, nk_sed),
          [&](const int idx, Scalar& lmin) {
            const int k = kmin_sed + idx;
            if (std::abs(velqr(col, k)) > Scalar(1e-12)) {
              const Scalar dz = z(col, k + lyr_step) - z(col, k);
              lmin = std::min(lmin, Scalar(0.8) * dz / velqr(col, k));
            }
          },
          Kokkos::Min<Scalar>(dtmin));
        Kokkos::single(Kokkos::PerTeam(team), [&]() {
          dt0(col) = std::min(dt0(col), dtmin);
        });
      });

    // -- Step 8: check convergence (returns to host) --
    int n_active = 0;
    Kokkos::parallel_reduce("kessler_converge",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols),
      KOKKOS_LAMBDA(const int col, int& cnt) {
        if (mask(col) != Scalar(0)) ++cnt;
      }, n_active);
    all_converged = (n_active == 0);

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
      const int col = idx / nk;
      const int k   = kmin + (idx % nk);
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
  const view_2d<const Scalar>& cpair,
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
      st_energy(col, k) = cpair(col, k) * temp(col, k) +
                          gravit * zm(col, k) + phis(col);
    });
  Kokkos::fence();
}

} // namespace kessler
} // namespace scream

#endif // KESSLER_FUNCTIONS_IMPL_HPP
