#ifndef KESSLER_FUNCTIONS_IMPL_HPP
#define KESSLER_FUNCTIONS_IMPL_HPP

#include "kessler_functions.hpp"

#include <ekat_assert.hpp>

#include <Kokkos_Core.hpp>

#include <limits>
#include <cmath>
#include <string>


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
 *    time_counter, precl_acc) live in the caller-owned Workspace `ws`
 *    (see kessler_functions.hpp), persistent across calls instead of
 *    allocated and freed here every time (see Buffer/init_buffers in
 *    eamxx_kessler_process_interface.hpp/.cpp for how the production
 *    caller carves them from one ATMBufferManager allocation, mirroring
 *    the P3/SHOC Buffer pattern).
 *  - Per-(col,level) fields, including the 6 Workspace 2D scratch
 *    arrays, are Pack-typed and operated on Pack-wide using ekat's
 *    Pack-aware math (ekat::exp/pow/sqrt/max/min, not std::, which has
 *    no Pack overloads) -- the same idiom SHOC uses throughout (e.g.
 *    shoc_assumed_pdf_impl.hpp).
 *  - Sedimentation needs the scalar neighbor level (k +/- lyr_step),
 *    which generally falls in an adjacent Pack, not the current one.
 *    That neighbor is assembled via ekat::shift_left/shift_right, the
 *    same mechanism P3's calc_first_order_upwind_step uses for its own
 *    vertical upwind stencil (p3_upwind_impl.hpp).
 *  - A handful of per-level operations are irregular (the CFL-limited
 *    dt0 reductions) or inherently single-level (the lyr_surf/lyr_toa
 *    boundary terms); these use ekat::scalarize() views instead of
 *    fighting Pack lanes for a reduction, matching how P3's own
 *    generalized_sedimentation uses scalarize() for its precip
 *    accumulation.
 *  - Levels beyond nz within the last Pack (padding, when nz is not a
 *    multiple of Pack::n) and, in the sedimentation kernel, the single
 *    level nearest lyr_toa, are excluded via an in-range Mask. Inputs
 *    that would otherwise divide by zero or hit a pow() domain error on
 *    excluded lanes are clamped to a safe placeholder before the op
 *    runs -- the computation itself must not trap under this codebase's
 *    fpe build, even though the result there is never read back out.
 *  - The sub-cycling while loop runs on the host; each iteration
 *    dispatches device kernels and then checks convergence via a
 *    device-to-host reduction. Each column's `dt0` is re-checked
 *    against the CFL condition every sub-cycle (using the just-
 *    updated rain fall speed), exactly as the reference Fortran
 *    `kessler_run` (module `kessler`, atmospheric_physics) and its
 *    JAX translation (`kessler_jax/kessler_run.py`) do -- not
 *    computed once up front. This matters: if the fall speed grows
 *    during the timestep, a `dt0` fixed from the initial fall speed
 *    would violate CFL for later sub-steps.
 *  - Convention: lyr_surf and lyr_toa are 0-based C++ indices.
 *    lyr_step = +1 if lyr_surf <= lyr_toa, -1 otherwise. EAMxx orders
 *    levels top-down, so the production caller passes
 *    lyr_surf = nlevs-1, lyr_toa = 0 (lyr_step = -1).
 */
template <typename S, typename D>
void KesslerMicrophysicsFunctions<S,D>::kessler_run(
  const Int ncols,
    const Int nz,
    const Real dt,
    const Int lyr_surf,
    const Int lyr_toa,
    const Real rhoqr,
    const Real lv,
    const Real pref,
    const Real cpair,
    const Real rair,
    const view_2d<const Pack>& rho,
    const view_2d<const Pack>& z_mid,
    const view_2d<const Pack>& pk,
    const view_2d<Pack>& theta,
    const view_2d<Pack>& qv,
    const view_2d<Pack>& qc,
    const view_2d<Pack>& qr,
    const view_1d<Scalar>& precl,
    const view_2d<Pack>& relhum,
    const Workspace& ws)
{
  // Physical constants (captured by value in lambdas)
  const Real f2x   = Real(17.27);
  const Real xk    = cpair / rair;              // cp/R
  const Real f5_dk = Real(4093) * lv / cpair;    // Durran & Klemp (1983) Eq A13 coefficient

  // Vertical direction
  const int lyr_step = (lyr_surf <= lyr_toa) ? 1 : -1;

  // Level ranges (0-based inclusive)
  const int kmin = (lyr_step > 0) ? lyr_surf : lyr_toa;
  const int kmax = (lyr_step > 0) ? lyr_toa  : lyr_surf;

  // Range for the sedimentation inner loop: all levels except lyr_toa
  const int kmin_sed = (lyr_step > 0) ? lyr_surf       : lyr_toa + 1;
  const int kmax_sed = (lyr_step > 0) ? lyr_toa - 1    : lyr_surf;
  const int nk_sed   = kmax_sed - kmin_sed + 1; // nz-1

  const int nlev_packs = ekat::npack<Pack>(nz);

  // ---------------------------------------------------------------
  // Workspace scratch: persistent, caller-owned (see Workspace in
  // kessler_functions.hpp).
  // ---------------------------------------------------------------
  const view_2d<Pack>& r     = ws.r;
  const view_2d<Pack>& rhalf = ws.rhalf;
  const view_2d<Pack>& velqr = ws.velqr;
  const view_2d<Pack>& sed   = ws.sed;
  const view_2d<Pack>& pc    = ws.pc;
  // ws.f5 is sized/carved as part of the Workspace but unused by the
  // current algorithm -- the f5 Durran & Klemp coefficient below is a
  // scalar constant, not a per-level field.
  const view_1d<Scalar>& dt0          = ws.dt0;
  const view_1d<Scalar>& mask         = ws.mask;
  const view_1d<Scalar>& time_counter = ws.time_counter;
  const view_1d<Scalar>& precl_acc    = ws.precl_acc;

  // Scalarized (zero-copy, flat per-scalar-level) views for the
  // handful of operations that are irregular (CFL reductions) or
  // inherently single-level (surface/TOA boundary terms).
  //
  // ekat::scalarize rebuilds a rank-2 view with the default LayoutRight
  // stride (extent(1)*Pack::n) and drops any stride the input carries, so
  // it is only valid for contiguous views. qr is NOT scalarized: in EAMxx
  // it is a subfield of the bundled "tracers" field (ncol x ntracers x
  // nlev), so its column stride is ntracers*nlev, and a scalarized qr reads
  // other tracers/columns. Single-level qr values are read through the
  // packed view instead: qr(col, k / Pack::n)[k % Pack::n].
  EKAT_REQUIRE_MSG(rho.span_is_contiguous() && z_mid.span_is_contiguous() &&
                   velqr.span_is_contiguous() && sed.span_is_contiguous(),
    "Error! KESSLER: kessler_run scalarizes rho, z_mid, velqr, and sed, "
    "which requires contiguous views.\n");
  const auto rho_s   = ekat::scalarize(rho);
  const auto z_s     = ekat::scalarize(z_mid);
  const auto velqr_s = ekat::scalarize(velqr);
  const auto sed_s   = ekat::scalarize(sed);

  // ---------------------------------------------------------------
  // Kernel 1: initialise derived constants and terminal fall speed
  // ---------------------------------------------------------------
  Kokkos::parallel_for("kessler_init_fields",
    Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nlev_packs),
    KOKKOS_LAMBDA(const int idx) {
      const int col = idx / nlev_packs;
      const int kp  = idx % nlev_packs;

      const auto level = ekat::range<IntPack>(kp * Pack::n);
      const auto valid  = (level >= kmin) && (level <= kmax);

      // Guard denominators/bases so invalid (padding or out-of-range)
      // lanes never hit a division-by-zero or a pow() domain error --
      // their results are never read downstream, but the computation
      // itself must not trap under the fpe build.
      Pack rho_kp = rho(col, kp);
      rho_kp.set(!valid, Real(1));
      Pack pk_kp = pk(col, kp);
      pk_kp.set(!valid, Real(1));

      const Scalar rho_surf = rho_s(col, lyr_surf);

      r(col, kp)     = Real(0.001) * rho_kp;          // g/cm^3
      rhalf(col, kp) = ekat::sqrt(rho_surf / rho_kp);
      pc(col, kp)    = Real(3.8) / (ekat::pow(pk_kp, xk) * pref);
      // Guard against round-off negative qr before computing velqr
      qr(col, kp)    = ekat::max(qr(col, kp), Real(0));
      velqr(col, kp) = Real(36.34) * rhalf(col, kp) *
                       ekat::pow(qr(col, kp) * r(col, kp), Real(0.1364));
    });
  Kokkos::fence();

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
          if (Kokkos::abs(velqr_s(col, k)) > Real(1e-12)) {
            const Scalar dzk = z_s(col, k + lyr_step) - z_s(col, k);
            lmin = Kokkos::min(lmin, Real(0.8) * dzk / velqr_s(col, k));
          }
        },
        Kokkos::Min<Scalar>(dtmin));

      Kokkos::single(Kokkos::PerTeam(team), [&]() {
        dt0(col)          = Kokkos::min(dt, dtmin);
        mask(col)         = Real(1);
        time_counter(col) = Real(0);
        precl_acc(col)    = Real(0);
        precl(col)        = Real(0);
      });
    });

  Kokkos::fence();

  // ---------------------------------------------------------------
  // Guard against a degenerate initial CFL time step. Mirrors the
  // reference Fortran kessler_run's "bad time splitting" abort
  // (kessler.F90: `if (dt0 < 1e-12) then errflg = 1; return`), which
  // the original port of this routine omitted, allowing a column
  // with a pathologically small CFL-limited dt0 to silently subcycle
  // instead of failing loudly.
  // ---------------------------------------------------------------
  Scalar dt0_min = std::numeric_limits<Scalar>::max();
  Kokkos::parallel_reduce("kessler_check_dt0",
    Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols),
    KOKKOS_LAMBDA (const int col, Scalar& lmin) {
      lmin = Kokkos::min(lmin, dt0(col));
    },
    Kokkos::Min<Scalar>(dt0_min));
  EKAT_REQUIRE_MSG(dt0_min >= Real(1e-12),
    "Error! KESSLER: bad time splitting (dt = " + std::to_string(dt) +
    ", dt0 = " + std::to_string(dt0_min) + ").\n");

  // ---------------------------------------------------------------
  // Sub-cycling while loop
  // Host drives the loop; convergence checked via parallel_reduce.
  // ---------------------------------------------------------------
  bool all_converged = false;
  int n_iter = 0;
  while (!all_converged) {

    // -- Step 1: accumulate precipitation (weighted by sub-cycle dt) --
    Kokkos::parallel_for("kessler_precl_accum",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols),
      KOKKOS_LAMBDA(const int col) {
        const Scalar qr_surf = qr(col, lyr_surf / Pack::n)[lyr_surf % Pack::n];
        const Scalar p_val = rho_s(col, lyr_surf) * qr_surf *
                             velqr_s(col, lyr_surf) / rhoqr;
        precl(col)     = p_val;
        precl_acc(col) += mask(col) * p_val * dt0(col);
      });

    // -- Step 2: sedimentation for all levels except lyr_toa --
    Kokkos::parallel_for("kessler_sed_inner",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nlev_packs),
      KOKKOS_LAMBDA(const int idx) {
        const int col = idx / nlev_packs;
        const int kp  = idx % nlev_packs;

        const auto level = ekat::range<IntPack>(kp * Pack::n);
        const auto valid  = (level >= kmin_sed) && (level <= kmax_sed);

        const Pack rqv_k = r(col, kp) * qr(col, kp) * velqr(col, kp);
        const Pack z_k   = z_mid(col, kp);

        // Neighbor pack in the array-index direction lyr_step points.
        // Clamped at the pack-array edge; the lane(s) that would
        // genuinely need an out-of-domain neighbor are always the
        // ones `valid` excludes (kmin_sed/kmax_sed already stop one
        // level short of lyr_toa by construction), so the clamped,
        // otherwise-unused value is never read.
        Pack rqv_up, z_up;
        if (lyr_step > 0) {
          const int kp_nbr = (kp + 1 < nlev_packs) ? kp + 1 : kp;
          const Pack rqv_next = r(col, kp_nbr) * qr(col, kp_nbr) * velqr(col, kp_nbr);
          rqv_up = ekat::shift_left(rqv_next, rqv_k);
          z_up   = ekat::shift_left(z_mid(col, kp_nbr), z_k);
        } else {
          const int kp_nbr = (kp > 0) ? kp - 1 : kp;
          const Pack rqv_prev = r(col, kp_nbr) * qr(col, kp_nbr) * velqr(col, kp_nbr);
          rqv_up = ekat::shift_right(rqv_prev, rqv_k);
          z_up   = ekat::shift_right(z_mid(col, kp_nbr), z_k);
        }

        Pack dz_up = z_up - z_k;
        dz_up.set(!valid, Real(1)); // guard div-by-zero on excluded lanes

        Pack r_k_g = r(col, kp);
        r_k_g.set(!valid, Real(1));

        sed(col, kp) = dt0(col) * (rqv_up - rqv_k) / (r_k_g * dz_up);
      });

    // -- Step 3: sedimentation at lyr_toa (no flux from above) --
    Kokkos::parallel_for("kessler_sed_toa",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols),
      KOKKOS_LAMBDA(const int col) {
        const int kbelow = lyr_toa - lyr_step;
        const Scalar qr_toa = qr(col, lyr_toa / Pack::n)[lyr_toa % Pack::n];
        sed_s(col, lyr_toa) = -dt0(col) * qr_toa *
          velqr_s(col, lyr_toa) /
          (Real(0.5) * (z_s(col, lyr_toa) - z_s(col, kbelow)));
      });

    // -- Step 4: microphysics adjustments (autoconversion, collection,
    //            evaporation, saturation adjustment) --
    Kokkos::parallel_for("kessler_micro",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nlev_packs),
      KOKKOS_LAMBDA(const int idx) {
        const int col = idx / nlev_packs;
        const int kp  = idx % nlev_packs;

        const auto level = ekat::range<IntPack>(kp * Pack::n);
        const auto valid  = (level >= kmin) && (level <= kmax);

        const Pack msk  = mask(col);
        const Pack dt0c = dt0(col);

        // pk's padding lanes are only ever multiplied elsewhere, but
        // here they're also a denominator (theta update below), so
        // guard them the same way kernel 1 does.
        Pack pk_kp = pk(col, kp);
        pk_kp.set(!valid, Real(1));

        // Autoconversion + collection (Klemp & Wilhelmson 1978, Eq 2.13a,b)
        // Semi-implicit treatment of the collection term.
        const Pack qrprod =
          qc(col, kp) -
          (qc(col, kp) - dt0c *
            ekat::max(Real(0.001) * (qc(col, kp) - Real(0.001)),
                      Real(0))) /
          (Real(1) + dt0c * Real(2.2) *
           ekat::pow(qr(col, kp), Real(0.875)));
        // Guard with msk: once a column has converged (elapsed time
        // reached dt), it must not be perturbed further by leftover
        // sub-cycle iterations still running for other columns --
        // matching the Fortran, where each column's do-while loop
        // exits independently and nothing more happens to it.
        qc(col, kp) = msk * ekat::max(qc(col, kp) - qrprod, Real(0))
                    + (Real(1) - msk) * qc(col, kp);
        qr(col, kp) = msk * ekat::max(qr(col, kp) + qrprod + sed(col, kp), Real(0))
                    + (Real(1) - msk) * qr(col, kp);

        // Saturation mixing ratio via Teten's formula
        // (Klemp & Wilhelmson 1978, Eq 2.11)
        const Pack T_pi = pk_kp * theta(col, kp); // temperature = pi*theta
        const Pack qvs  = pc(col, kp) *
          ekat::exp(f2x * (T_pi - Real(273)) / (T_pi - Real(36)));

        // Condensation rate (Durran & Klemp 1983, Eq A13-A14)
        const Pack dT_pi_36 = T_pi - Real(36);
        const Pack prod = (qv(col, kp) - qvs) /
          (Real(1) + qvs * f5_dk / (dT_pi_36 * dT_pi_36));

        // Evaporation rate (Klemp & Wilhelmson 1978, Eq 2.14a,b)
        // dim(qvs, qv) = max(qvs - qv, 0): evaporation only in subsaturated air
        const Pack rqr = r(col, kp) * qr(col, kp);
        const Pack A = dt0c *
          ((Real(1.6) + Real(124.9) * ekat::pow(rqr, Real(0.2046))) *
           ekat::pow(rqr, Real(0.525)) /
           (Real(2550000) * pc(col, kp) / (Real(3.8) * qvs) +
            Real(540000))) *
          (ekat::max(qvs - qv(col, kp), Real(0)) / (r(col, kp) * qvs));
        const Pack B   = ekat::max(-prod - qc(col, kp), Real(0));
        const Pack ern = ekat::min(A, ekat::min(B, qr(col, kp)));

        // Saturation adjustment (Durran & Klemp 1983, Eq A1-A4)
        // Applied only to non-converged columns via mask.
        const Pack prod_adj = ekat::max(prod, -qc(col, kp));
        theta(col, kp) += msk * (lv / (cpair * pk_kp) *
                                (prod_adj - ern));
        qv(col, kp) = msk * ekat::max(qv(col, kp) - prod_adj + ern, Real(0))
                    + (Real(1) - msk) * qv(col, kp);
        qc(col, kp) += msk * prod_adj;
        qr(col, kp) = msk * ekat::max(qr(col, kp) - ern, Real(0))
                    + (Real(1) - msk) * qr(col, kp);
      });

    // -- Step 5: advance elapsed time and update mask / dt0 --
    Kokkos::parallel_for("kessler_update_time",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols),
      KOKKOS_LAMBDA(const int col) {
        time_counter(col) += mask(col) * dt0(col);
        dt0(col) = Kokkos::max(dt - time_counter(col), Real(0));
        mask(col) = (Kokkos::abs(dt - time_counter(col)) > Real(1e-5))
                      ? Real(1) : Real(0);
      });

    // -- Step 6: recompute terminal fall speed with updated qr --
    Kokkos::parallel_for("kessler_velqr_update",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nlev_packs),
      KOKKOS_LAMBDA(const int idx) {
        const int col = idx / nlev_packs;
        const int kp  = idx % nlev_packs;
        velqr(col, kp) = Real(36.34) * rhalf(col, kp) *
                        ekat::pow(qr(col, kp) * r(col, kp), Real(0.1364));
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
            if (Kokkos::abs(velqr_s(col, k)) > Real(1e-12)) {
              const Scalar dzk = z_s(col, k + lyr_step) - z_s(col, k);
              lmin = Kokkos::min(lmin, Real(0.8) * dzk / velqr_s(col, k));
            }
          },
          Kokkos::Min<Scalar>(dtmin));
        Kokkos::single(Kokkos::PerTeam(team), [&]() {
          dt0(col) = Kokkos::min(dt0(col), dtmin);
        });
      });

    // -- Step 8: check convergence (returns to host) --
    // Also count active columns whose sub-cycle step is non-positive:
    // those can never reach dt, and the loop would spin forever.
    int n_active = 0;
    int n_stuck  = 0;
    Kokkos::parallel_reduce("kessler_converge",
      Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols),
      KOKKOS_LAMBDA(const int col, int& cnt, int& stuck) {
        if (mask(col) != Real(0)) {
          ++cnt;
          if (!(dt0(col) > Real(0))) ++stuck;
        }
      }, n_active, n_stuck);
    all_converged = (n_active == 0);

    ++n_iter;
    EKAT_REQUIRE_MSG(n_stuck == 0 && n_iter < 100000,
      "Error! KESSLER: sub-cycling cannot converge (iter = " +
      std::to_string(n_iter) + ", active cols = " + std::to_string(n_active) +
      ", cols with dt0 <= 0 = " + std::to_string(n_stuck) + ").\n");

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
    Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nlev_packs),
    KOKKOS_LAMBDA(const int idx) {
      const int col = idx / nlev_packs;
      const int kp  = idx % nlev_packs;
      const Pack T_pi = pk(col, kp) * theta(col, kp);
      const Pack qvs  = pc(col, kp) *
        ekat::exp(f2x * (T_pi - Real(273)) / (T_pi - Real(36)));
      relhum(col, kp) = qv(col, kp) / qvs * Real(100);
    });
  Kokkos::fence();
}

// -----------------------------------------------------------------------

template <typename S, typename D>
void KesslerMicrophysicsFunctions<S,D>::kessler_update_timestep_init(
  const Int ncols,
    const Int nz,
    const view_2d<const Pack>& temp,
    const view_2d<Pack>& temp_prev,
    const view_2d<Pack>& ttend_t)
{
  const int nlev_packs = ekat::npack<Pack>(nz);
  Kokkos::parallel_for("kessler_ts_init",
    Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nlev_packs),
    KOKKOS_LAMBDA(const int idx) {
      const int col = idx / nlev_packs;
      const int kp  = idx % nlev_packs;
      temp_prev(col, kp) = temp(col, kp);
      ttend_t(col, kp)   = Real(0);
    });
  Kokkos::fence();
}

// -----------------------------------------------------------------------

template <typename S, typename D>
void KesslerMicrophysicsFunctions<S,D>::kessler_update_run(
  const Int ncols,
    const Int nz,
    const Real dt,
    const view_2d<const Pack>& theta,
    const view_2d<const Pack>& exner,
    const view_2d<const Pack>& temp_prev,
    const view_2d<Pack>& ttend_t)
{
  const int nlev_packs = ekat::npack<Pack>(nz);
  Kokkos::parallel_for("kessler_update_run",
    Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nlev_packs),
    KOKKOS_LAMBDA(const int idx) {
      const int col = idx / nlev_packs;
      const int kp  = idx % nlev_packs;
      ttend_t(col, kp) += (theta(col, kp) * exner(col, kp) - temp_prev(col, kp)) / dt;
    });
  Kokkos::fence();
}

// -----------------------------------------------------------------------

template <typename S, typename D>
void KesslerMicrophysicsFunctions<S,D>::kessler_update_timestep_final(
  const Int ncols,
    const Int nz,
    const Real gravit,
    const Real cpair,
    const view_2d<const Pack>& temp,
    const view_2d<const Pack>& z_mid,
    const view_1d<const Scalar>& phis,
    const view_2d<Pack>& st_energy)
{
  const int nlev_packs = ekat::npack<Pack>(nz);
  Kokkos::parallel_for("kessler_ts_final",
    Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nlev_packs),
    KOKKOS_LAMBDA(const int idx) {
      const int col = idx / nlev_packs;
      const int kp  = idx % nlev_packs;
      st_energy(col, kp) = cpair * temp(col, kp) +
                          gravit * z_mid(col, kp) + phis(col);
    });
  Kokkos::fence();
}

} // namespace kessler
} // namespace scream

#endif // KESSLER_FUNCTIONS_IMPL_HPP
