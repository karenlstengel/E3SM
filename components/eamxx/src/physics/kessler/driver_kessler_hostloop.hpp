#ifndef DRIVER_KESSLER_HOSTLOOP_HPP
#define DRIVER_KESSLER_HOSTLOOP_HPP

/*
 * Reference copy of kessler_run as it was before the sub-cycle loop moved
 * into the kernel (kessler_functions_impl.hpp at 49bf6ccdde): the host
 * drives the sub-cycle loop, and every pass sweeps all columns with ~9
 * kernels. Kept next to driver_kessler.cpp so the current kessler_run can be
 * checked against it bit for bit and timed against it
 * (performance_suggestions.md, item D1).
 *
 * Same arguments and Workspace as KMF::kessler_run. n_active, if given,
 * receives the number of still-active columns after each pass, from which
 * the driver derives every column's sub-cycle count: n_active[i] columns
 * needed more than i+1 passes. The arithmetic below is unchanged from the
 * original, apart from that bookkeeping.
 */

#include "kessler_functions.hpp"

#include <ekat_assert.hpp>
#include <Kokkos_Core.hpp>

#include <limits>
#include <string>
#include <vector>

namespace kessler_driver {

template <typename KMF>
void kessler_run_hostloop(
  const int ncols,
  const int nz,
  const typename KMF::Scalar dt,
  const int lyr_surf,
  const int lyr_toa,
  const typename KMF::Scalar rhoqr,
  const typename KMF::Scalar lv,
  const typename KMF::Scalar pref,
  const typename KMF::Scalar cpair,
  const typename KMF::Scalar rair,
  const typename KMF::template view_2d<const typename KMF::Pack>& rho,
  const typename KMF::template view_2d<const typename KMF::Pack>& z_mid,
  const typename KMF::template view_2d<const typename KMF::Pack>& pk,
  const typename KMF::template view_2d<typename KMF::Pack>& theta,
  const typename KMF::template view_2d<typename KMF::Pack>& qv,
  const typename KMF::template view_2d<typename KMF::Pack>& qc,
  const typename KMF::template view_2d<typename KMF::Pack>& qr,
  const typename KMF::template view_1d<typename KMF::Scalar>& precl,
  const typename KMF::template view_2d<typename KMF::Pack>& relhum,
  const typename KMF::Workspace& ws,
  std::vector<int>* n_active_out = nullptr)
{
  using Real     = typename KMF::Scalar;
  using Scalar   = typename KMF::Scalar;
  using Pack     = typename KMF::Pack;
  using IntPack  = typename KMF::IntPack;
  using ExeSpace = typename KMF::KT::ExeSpace;

  const Real f2x   = Real(17.27);
  const Real xk    = cpair / rair;
  const Real f5_dk = Real(4093) * lv / cpair;

  const int lyr_step = (lyr_surf <= lyr_toa) ? 1 : -1;
  const int kmin = (lyr_step > 0) ? lyr_surf : lyr_toa;
  const int kmax = (lyr_step > 0) ? lyr_toa  : lyr_surf;
  const int kmin_sed = (lyr_step > 0) ? lyr_surf       : lyr_toa + 1;
  const int kmax_sed = (lyr_step > 0) ? lyr_toa - 1    : lyr_surf;
  const int nk_sed   = kmax_sed - kmin_sed + 1;

  const int nlev_packs = ekat::npack<Pack>(nz);

  const auto r     = ws.r;
  const auto rhalf = ws.rhalf;
  const auto velqr = ws.velqr;
  const auto sed   = ws.sed;
  const auto pc    = ws.pc;
  const auto dt0          = ws.dt0;
  const auto mask         = ws.mask;
  const auto time_counter = ws.time_counter;
  const auto precl_acc    = ws.precl_acc;

  EKAT_REQUIRE_MSG(rho.span_is_contiguous() && z_mid.span_is_contiguous() &&
                   velqr.span_is_contiguous() && sed.span_is_contiguous(),
    "Error! KESSLER: kessler_run scalarizes rho, z_mid, velqr, and sed, "
    "which requires contiguous views.\n");
  const auto rho_s   = ekat::scalarize(rho);
  const auto z_s     = ekat::scalarize(z_mid);
  const auto velqr_s = ekat::scalarize(velqr);
  const auto sed_s   = ekat::scalarize(sed);

  if (n_active_out) n_active_out->clear();

  Kokkos::parallel_for("kessler_init_fields",
    Kokkos::RangePolicy<ExeSpace>(0, ncols * nlev_packs),
    KOKKOS_LAMBDA(const int idx) {
      const int col = idx / nlev_packs;
      const int kp  = idx % nlev_packs;
      const auto level = ekat::range<IntPack>(kp * Pack::n);
      const auto valid  = (level >= kmin) && (level <= kmax);
      Pack rho_kp = rho(col, kp);
      rho_kp.set(!valid, Real(1));
      Pack pk_kp = pk(col, kp);
      pk_kp.set(!valid, Real(1));
      const Scalar rho_surf = rho_s(col, lyr_surf);
      r(col, kp)     = Real(0.001) * rho_kp;
      rhalf(col, kp) = ekat::sqrt(rho_surf / rho_kp);
      pc(col, kp)    = Real(3.8) / (ekat::pow(pk_kp, xk) * pref);
      qr(col, kp)    = ekat::max(qr(col, kp), Real(0));
      velqr(col, kp) = Real(36.34) * rhalf(col, kp) *
                       ekat::pow(qr(col, kp) * r(col, kp), Real(0.1364));
    });
  Kokkos::fence();

  using TeamPol    = Kokkos::TeamPolicy<ExeSpace>;
  using MemberType = typename TeamPol::member_type;

  Kokkos::parallel_for("kessler_init_dt0",
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
        dt0(col)          = Kokkos::min(dt, dtmin);
        mask(col)         = Real(1);
        time_counter(col) = Real(0);
        precl_acc(col)    = Real(0);
        precl(col)        = Real(0);
      });
    });
  Kokkos::fence();

  Scalar dt0_min = std::numeric_limits<Scalar>::max();
  Kokkos::parallel_reduce("kessler_check_dt0",
    Kokkos::RangePolicy<ExeSpace>(0, ncols),
    KOKKOS_LAMBDA (const int col, Scalar& lmin) {
      lmin = Kokkos::min(lmin, dt0(col));
    },
    Kokkos::Min<Scalar>(dt0_min));
  EKAT_REQUIRE_MSG(dt0_min >= Real(1e-12),
    "Error! KESSLER: bad time splitting (dt = " + std::to_string(dt) +
    ", dt0 = " + std::to_string(dt0_min) + ").\n");

  bool all_converged = false;
  int n_iter = 0;
  while (!all_converged) {

    Kokkos::parallel_for("kessler_precl_accum",
      Kokkos::RangePolicy<ExeSpace>(0, ncols),
      KOKKOS_LAMBDA(const int col) {
        const Scalar qr_surf = qr(col, lyr_surf / Pack::n)[lyr_surf % Pack::n];
        const Scalar p_val = rho_s(col, lyr_surf) * qr_surf *
                             velqr_s(col, lyr_surf) / rhoqr;
        precl(col)     = p_val;
        precl_acc(col) += mask(col) * p_val * dt0(col);
      });

    Kokkos::parallel_for("kessler_sed_inner",
      Kokkos::RangePolicy<ExeSpace>(0, ncols * nlev_packs),
      KOKKOS_LAMBDA(const int idx) {
        const int col = idx / nlev_packs;
        const int kp  = idx % nlev_packs;
        const auto level = ekat::range<IntPack>(kp * Pack::n);
        const auto valid  = (level >= kmin_sed) && (level <= kmax_sed);
        const Pack rqv_k = r(col, kp) * qr(col, kp) * velqr(col, kp);
        const Pack z_k   = z_mid(col, kp);
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
        dz_up.set(!valid, Real(1));
        Pack r_k_g = r(col, kp);
        r_k_g.set(!valid, Real(1));
        sed(col, kp) = dt0(col) * (rqv_up - rqv_k) / (r_k_g * dz_up);
      });

    Kokkos::parallel_for("kessler_sed_toa",
      Kokkos::RangePolicy<ExeSpace>(0, ncols),
      KOKKOS_LAMBDA(const int col) {
        const int kbelow = lyr_toa - lyr_step;
        const Scalar qr_toa = qr(col, lyr_toa / Pack::n)[lyr_toa % Pack::n];
        sed_s(col, lyr_toa) = -dt0(col) * qr_toa *
          velqr_s(col, lyr_toa) /
          (Real(0.5) * (z_s(col, lyr_toa) - z_s(col, kbelow)));
      });

    Kokkos::parallel_for("kessler_micro",
      Kokkos::RangePolicy<ExeSpace>(0, ncols * nlev_packs),
      KOKKOS_LAMBDA(const int idx) {
        const int col = idx / nlev_packs;
        const int kp  = idx % nlev_packs;
        const auto level = ekat::range<IntPack>(kp * Pack::n);
        const auto valid  = (level >= kmin) && (level <= kmax);
        const Pack msk  = mask(col);
        const Pack dt0c = dt0(col);
        Pack pk_kp = pk(col, kp);
        pk_kp.set(!valid, Real(1));
        const Pack qrprod =
          qc(col, kp) -
          (qc(col, kp) - dt0c *
            ekat::max(Real(0.001) * (qc(col, kp) - Real(0.001)),
                      Real(0))) /
          (Real(1) + dt0c * Real(2.2) *
           ekat::pow(qr(col, kp), Real(0.875)));
        qc(col, kp) = msk * ekat::max(qc(col, kp) - qrprod, Real(0))
                    + (Real(1) - msk) * qc(col, kp);
        qr(col, kp) = msk * ekat::max(qr(col, kp) + qrprod + sed(col, kp), Real(0))
                    + (Real(1) - msk) * qr(col, kp);
        const Pack T_pi = pk_kp * theta(col, kp);
        const Pack qvs  = pc(col, kp) *
          ekat::exp(f2x * (T_pi - Real(273)) / (T_pi - Real(36)));
        const Pack dT_pi_36 = T_pi - Real(36);
        const Pack prod = (qv(col, kp) - qvs) /
          (Real(1) + qvs * f5_dk / (dT_pi_36 * dT_pi_36));
        const Pack rqr = r(col, kp) * qr(col, kp);
        const Pack A = dt0c *
          ((Real(1.6) + Real(124.9) * ekat::pow(rqr, Real(0.2046))) *
           ekat::pow(rqr, Real(0.525)) /
           (Real(2550000) * pc(col, kp) / (Real(3.8) * qvs) +
            Real(540000))) *
          (ekat::max(qvs - qv(col, kp), Real(0)) / (r(col, kp) * qvs));
        const Pack B   = ekat::max(-prod - qc(col, kp), Real(0));
        const Pack ern = ekat::min(A, ekat::min(B, qr(col, kp)));
        const Pack prod_adj = ekat::max(prod, -qc(col, kp));
        theta(col, kp) += msk * (lv / (cpair * pk_kp) *
                                (prod_adj - ern));
        qv(col, kp) = msk * ekat::max(qv(col, kp) - prod_adj + ern, Real(0))
                    + (Real(1) - msk) * qv(col, kp);
        qc(col, kp) += msk * prod_adj;
        qr(col, kp) = msk * ekat::max(qr(col, kp) - ern, Real(0))
                    + (Real(1) - msk) * qr(col, kp);
      });

    Kokkos::parallel_for("kessler_update_time",
      Kokkos::RangePolicy<ExeSpace>(0, ncols),
      KOKKOS_LAMBDA(const int col) {
        time_counter(col) += mask(col) * dt0(col);
        dt0(col) = Kokkos::max(dt - time_counter(col), Real(0));
        mask(col) = (Kokkos::abs(dt - time_counter(col)) > Real(1e-5))
                      ? Real(1) : Real(0);
      });

    Kokkos::parallel_for("kessler_velqr_update",
      Kokkos::RangePolicy<ExeSpace>(0, ncols * nlev_packs),
      KOKKOS_LAMBDA(const int idx) {
        const int col = idx / nlev_packs;
        const int kp  = idx % nlev_packs;
        velqr(col, kp) = Real(36.34) * rhalf(col, kp) *
                        ekat::pow(qr(col, kp) * r(col, kp), Real(0.1364));
      });

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

    int n_active = 0;
    int n_stuck  = 0;
    Kokkos::parallel_reduce("kessler_converge",
      Kokkos::RangePolicy<ExeSpace>(0, ncols),
      KOKKOS_LAMBDA(const int col, int& cnt, int& stuck) {
        if (mask(col) != Real(0)) {
          ++cnt;
          if (!(dt0(col) > Real(0))) ++stuck;
        }
      }, n_active, n_stuck);
    all_converged = (n_active == 0);
    if (n_active_out) n_active_out->push_back(n_active);

    ++n_iter;
    EKAT_REQUIRE_MSG(n_stuck == 0 && n_iter < 100000,
      "Error! KESSLER: sub-cycling cannot converge (iter = " +
      std::to_string(n_iter) + ", active cols = " + std::to_string(n_active) +
      ", cols with dt0 <= 0 = " + std::to_string(n_stuck) + ").\n");
  }

  Kokkos::parallel_for("kessler_precl_final",
    Kokkos::RangePolicy<ExeSpace>(0, ncols),
    KOKKOS_LAMBDA(const int col) {
      precl(col) = precl_acc(col) / dt;
    });

  Kokkos::parallel_for("kessler_relhum",
    Kokkos::RangePolicy<ExeSpace>(0, ncols * nlev_packs),
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

} // namespace kessler_driver

#endif // DRIVER_KESSLER_HOSTLOOP_HPP
