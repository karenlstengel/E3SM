#ifndef DRIVER_KESSLER_VARIANTS_HPP
#define DRIVER_KESSLER_VARIANTS_HPP

/*
 * Experimental variants of kessler_run, for timing and checking against the
 * production KMF::kessler_run in driver_kessler.cpp before any of them goes
 * into kessler_functions_impl.hpp (performance_suggestions.md, section I).
 *
 * kessler_run_variant<KMF, SkipDryPow, TeamScratch> is kessler_run with two
 * compile-time switches; every variant must match production bit for bit.
 * SkipDryPow was ported into production kessler_run (as
 * detail::pow_nonneg) on 2026-10-28, so <true, false> is now the
 * production algorithm and <false, false> is production before I2.
 *
 *   SkipDryPow (I2): a pow whose base pack is 0 in every lane is replaced by
 *     0 instead of being called. pow(+-0, y) = +0 for y > 0, so the result is
 *     identical. Applies to pow(qr, 0.875) (autoconversion/collection),
 *     pow(r*qr, 0.2046) and pow(r*qr, 0.525) (evaporation) and
 *     pow(r*qr, 0.1364) (fall speed). qr >= 0 throughout (it is clipped).
 *
 *   TeamScratch (I6): the per-level scratch of one column (r, rhalf, pc,
 *     velqr, sed) lives in Kokkos team scratch (level 0: shared memory on
 *     GPU, a per-team buffer on CPU) instead of the global Workspace arrays.
 *     ws.r/rhalf/pc/velqr/sed are then not written; the per-column outputs
 *     (ws.dt0, mask, time_counter, precl_acc) are.
 *
 * Same arguments, error codes and errors as KMF::kessler_run.
 */

#include "kessler_functions.hpp"

#include <ekat_assert.hpp>
#include <Kokkos_Core.hpp>

#include <string>

namespace kessler_driver {

// pow(x, e), or exactly 0 without calling pow when every lane of x is 0
// (I2, if Skip). Callers pass x >= 0, and pow(+-0, e > 0) = +0, so the
// result is the same either way.
template <bool Skip, typename PackT, typename RealT>
KOKKOS_INLINE_FUNCTION
PackT pow_or_zero (const PackT& x, const RealT e) {
  if constexpr (Skip) {
    if (!(x != RealT(0)).any()) return PackT(RealT(0));
  }
  return ekat::pow(x, e);
}

template <typename KMF, bool SkipDryPow, bool TeamScratch>
void kessler_run_variant(
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
  const typename KMF::Workspace& ws)
{
  using Real       = typename KMF::Scalar;
  using Scalar     = typename KMF::Scalar;
  using Pack       = typename KMF::Pack;
  using IntPack    = typename KMF::IntPack;
  using ExeSpace   = typename KMF::KT::ExeSpace;
  using MemberType = typename KMF::MemberType;
  using TPF        = ekat::TeamPolicyFactory<ExeSpace>;
  using ScratchPacks = Kokkos::View<Pack*, typename ExeSpace::scratch_memory_space,
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  const Real f2x   = Real(17.27);
  const Real xk    = cpair / rair;
  const Real f5_dk = Real(4093) * lv / cpair;

  const int lyr_step = (lyr_surf <= lyr_toa) ? 1 : -1;
  const int kmin     = (lyr_step > 0) ? lyr_surf : lyr_toa;
  const int kmax     = (lyr_step > 0) ? lyr_toa  : lyr_surf;
  const int kmin_sed = (lyr_step > 0) ? lyr_surf     : lyr_toa + 1;
  const int kmax_sed = (lyr_step > 0) ? lyr_toa - 1  : lyr_surf;
  const int kbelow   = lyr_toa - lyr_step;

  const int nlev_packs = ekat::npack<Pack>(nz);
  const int kp_surf = lyr_surf / Pack::n, ks_surf = lyr_surf % Pack::n;
  const int kp_toa  = lyr_toa  / Pack::n, ks_toa  = lyr_toa  % Pack::n;

  const int max_pass      = 100000;
  const int err_bad_split = 2;
  const int err_no_conv   = 3;

  const auto ws_r     = ws.r;
  const auto ws_rhalf = ws.rhalf;
  const auto ws_velqr = ws.velqr;
  const auto ws_sed   = ws.sed;
  const auto ws_pc    = ws.pc;
  const auto dt0_out      = ws.dt0;
  const auto mask         = ws.mask;
  const auto time_counter = ws.time_counter;
  const auto precl_acc    = ws.precl_acc;

  EKAT_REQUIRE_MSG(rho.span_is_contiguous() && z_mid.span_is_contiguous(),
    "Error! KESSLER: kessler_run_variant scalarizes rho and z_mid, "
    "which requires contiguous views.\n");
  const auto rho_s = ekat::scalarize(rho);
  const auto z_s   = ekat::scalarize(z_mid);

  constexpr int n_scratch = 5;  // r, rhalf, pc, velqr, sed
  auto policy = TPF::get_default_team_policy(ncols, nlev_packs);
  if constexpr (TeamScratch) {
    policy.set_scratch_size(0, Kokkos::PerTeam(n_scratch * ScratchPacks::shmem_size(nlev_packs)));
  }

  Kokkos::parallel_for("kessler_run_variant", policy,
    KOKKOS_LAMBDA(const MemberType& team) {
      const int col = team.league_rank();

      // The column body, written once for both storage choices: r_c, rhalf_c,
      // pc_c, velqr_c and sed_c are this column's level arrays, either
      // subviews of the Workspace or team-scratch views.
      auto body = [&](const auto& r_c, const auto& rhalf_c, const auto& pc_c,
                      const auto& velqr_c, const auto& sed_c) {
        auto kp_nbr = [&](const int kp) {
          return (lyr_step > 0) ? ((kp + 1 < nlev_packs) ? kp + 1 : kp)
                                : ((kp > 0) ? kp - 1 : kp);
        };
        auto shift_up = [&](const Pack& nbr, const Pack& p) {
          return (lyr_step > 0) ? ekat::shift_left(nbr, p) : ekat::shift_right(nbr, p);
        };
        auto sed_valid = [&](const int kp) {
          const auto level = ekat::range<IntPack>(kp * Pack::n);
          return (level >= kmin_sed) && (level <= kmax_sed);
        };
        auto cfl_min = [&](const int kp, const Pack& v, Scalar& lmin) {
          const Pack z_k = z_mid(col, kp);
          const Pack dz  = shift_up(z_mid(col, kp_nbr(kp)), z_k) - z_k;
          const auto m   = sed_valid(kp) && (ekat::abs(v) > Real(1e-12));
          ekat_masked_loop_no_vec(m, s) {
            lmin = Kokkos::min(lmin, Real(0.8) * dz[s] / v[s]);
          }
        };

        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_packs), [&](const int kp) {
          const auto level = ekat::range<IntPack>(kp * Pack::n);
          const auto valid = (level >= kmin) && (level <= kmax);
          Pack rho_kp = rho(col, kp);
          rho_kp.set(!valid, Real(1));
          Pack pk_kp = pk(col, kp);
          pk_kp.set(!valid, Real(1));
          const Scalar rho_surf = rho_s(col, lyr_surf);
          r_c(kp)     = Real(0.001) * rho_kp;
          rhalf_c(kp) = ekat::sqrt(rho_surf / rho_kp);
          pc_c(kp)    = Real(3.8) / (ekat::pow(pk_kp, xk) * pref);
          qr(col, kp) = ekat::max(qr(col, kp), Real(0));
          velqr_c(kp) = Real(36.34) * rhalf_c(kp) * pow_or_zero<SkipDryPow>(qr(col, kp) * r_c(kp), Real(0.1364));
        });
        team.team_barrier();

        Scalar dt0;
        {
          Scalar dtmin;
          Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, nlev_packs),
            [&](const int kp, Scalar& lmin) { cfl_min(kp, velqr_c(kp), lmin); },
            Kokkos::Min<Scalar>(dtmin));
          dt0 = Kokkos::min(dt, dtmin);
        }

        int    err    = (dt0 < Real(1e-12)) ? err_bad_split : 0;
        bool   active = (err == 0);
        Scalar tc     = 0;
        Scalar pacc   = 0;
        int    n      = 0;

        while (active) {
          const Scalar qr_surf = qr(col, kp_surf)[ks_surf];
          const Scalar p_val = rho_s(col, lyr_surf) * qr_surf *
                               velqr_c(kp_surf)[ks_surf] / rhoqr;
          pacc += p_val * dt0;

          Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_packs), [&](const int kp) {
            const auto valid = sed_valid(kp);
            const int  kq    = kp_nbr(kp);
            const Pack rqv_k  = r_c(kp) * qr(col, kp) * velqr_c(kp);
            const Pack z_k    = z_mid(col, kp);
            const Pack rqv_up = shift_up(r_c(kq) * qr(col, kq) * velqr_c(kq), rqv_k);
            const Pack z_up   = shift_up(z_mid(col, kq), z_k);
            Pack dz_up = z_up - z_k;
            dz_up.set(!valid, Real(1));
            Pack r_k_g = r_c(kp);
            r_k_g.set(!valid, Real(1));
            Pack s = dt0 * (rqv_up - rqv_k) / (r_k_g * dz_up);
            if (kp == kp_toa) {
              const Scalar qr_toa = qr(col, kp_toa)[ks_toa];
              s[ks_toa] = -dt0 * qr_toa * velqr_c(kp_toa)[ks_toa] /
                          (Real(0.5) * (z_s(col, lyr_toa) - z_s(col, kbelow)));
            }
            sed_c(kp) = s;
          });
          team.team_barrier();

          Scalar dtmin;
          Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, nlev_packs),
            [&](const int kp, Scalar& lmin) {
              const auto level = ekat::range<IntPack>(kp * Pack::n);
              const auto valid = (level >= kmin) && (level <= kmax);
              Pack pk_kp = pk(col, kp);
              pk_kp.set(!valid, Real(1));

              const Pack qrprod =
                qc(col, kp) -
                (qc(col, kp) - dt0 *
                  ekat::max(Real(0.001) * (qc(col, kp) - Real(0.001)), Real(0))) /
                (Real(1) + dt0 * Real(2.2) * pow_or_zero<SkipDryPow>(qr(col, kp), Real(0.875)));
              qc(col, kp) = ekat::max(qc(col, kp) - qrprod, Real(0));
              qr(col, kp) = ekat::max(qr(col, kp) + qrprod + sed_c(kp), Real(0));

              const Pack T_pi = pk_kp * theta(col, kp);
              const Pack qvs  = pc_c(kp) *
                ekat::exp(f2x * (T_pi - Real(273)) / (T_pi - Real(36)));
              const Pack dT_pi_36 = T_pi - Real(36);
              const Pack prod = (qv(col, kp) - qvs) /
                (Real(1) + qvs * f5_dk / (dT_pi_36 * dT_pi_36));

              const Pack rqr = r_c(kp) * qr(col, kp);
              const Pack A = dt0 *
                ((Real(1.6) + Real(124.9) * pow_or_zero<SkipDryPow>(rqr, Real(0.2046))) *
                 pow_or_zero<SkipDryPow>(rqr, Real(0.525)) /
                 (Real(2550000) * pc_c(kp) / (Real(3.8) * qvs) + Real(540000))) *
                (ekat::max(qvs - qv(col, kp), Real(0)) / (r_c(kp) * qvs));
              const Pack B   = ekat::max(-prod - qc(col, kp), Real(0));
              const Pack ern = ekat::min(A, ekat::min(B, qr(col, kp)));

              const Pack prod_adj = ekat::max(prod, -qc(col, kp));
              theta(col, kp) += lv / (cpair * pk_kp) * (prod_adj - ern);
              qv(col, kp) = ekat::max(qv(col, kp) - prod_adj + ern, Real(0));
              qc(col, kp) += prod_adj;
              qr(col, kp) = ekat::max(qr(col, kp) - ern, Real(0));

              const Pack v = Real(36.34) * rhalf_c(kp) *
                             pow_or_zero<SkipDryPow>(qr(col, kp) * r_c(kp), Real(0.1364));
              velqr_c(kp) = v;
              cfl_min(kp, v, lmin);
            },
            Kokkos::Min<Scalar>(dtmin));

          tc += dt0;
          dt0 = Kokkos::max(dt - tc, Real(0));
          active = Kokkos::abs(dt - tc) > Real(1e-5);
          dt0 = Kokkos::min(dt0, dtmin);
          ++n;
          if (active && (!(dt0 > Real(0)) || n >= max_pass)) {
            err = err_no_conv;
            active = false;
          }
          team.team_barrier();
        }

        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_packs), [&](const int kp) {
          const Pack T_pi = pk(col, kp) * theta(col, kp);
          const Pack qvs  = pc_c(kp) *
            ekat::exp(f2x * (T_pi - Real(273)) / (T_pi - Real(36)));
          relhum(col, kp) = qv(col, kp) / qvs * Real(100);
        });

        Kokkos::single(Kokkos::PerTeam(team), [&]() {
          precl(col)        = pacc / dt;
          precl_acc(col)    = pacc;
          dt0_out(col)      = dt0;
          time_counter(col) = tc;
          mask(col)         = Real(err);
        });
      };

      // The Workspace subviews are made outside the if-constexpr: nvcc does
      // not allow an extended lambda to first capture a variable inside one.
      const auto r_w     = Kokkos::subview(ws_r,     col, Kokkos::ALL);
      const auto rhalf_w = Kokkos::subview(ws_rhalf, col, Kokkos::ALL);
      const auto pc_w    = Kokkos::subview(ws_pc,    col, Kokkos::ALL);
      const auto velqr_w = Kokkos::subview(ws_velqr, col, Kokkos::ALL);
      const auto sed_w   = Kokkos::subview(ws_sed,   col, Kokkos::ALL);
      if constexpr (TeamScratch) {
        const ScratchPacks r_c    (team.team_scratch(0), nlev_packs);
        const ScratchPacks rhalf_c(team.team_scratch(0), nlev_packs);
        const ScratchPacks pc_c   (team.team_scratch(0), nlev_packs);
        const ScratchPacks velqr_c(team.team_scratch(0), nlev_packs);
        const ScratchPacks sed_c  (team.team_scratch(0), nlev_packs);
        body(r_c, rhalf_c, pc_c, velqr_c, sed_c);
      } else {
        body(r_w, rhalf_w, pc_w, velqr_w, sed_w);
      }
    });

  int n_bad_split = 0;
  int n_no_conv   = 0;
  Kokkos::parallel_reduce("kessler_run_variant_check",
    Kokkos::RangePolicy<ExeSpace>(0, ncols),
    KOKKOS_LAMBDA(const int col, int& bad_split, int& no_conv) {
      if (mask(col) == Real(err_bad_split)) ++bad_split;
      if (mask(col) == Real(err_no_conv))   ++no_conv;
    }, n_bad_split, n_no_conv);
  EKAT_REQUIRE_MSG(n_bad_split == 0,
    "Error! KESSLER: bad time splitting (dt = " + std::to_string(dt) +
    ", columns with dt0 < 1e-12 = " + std::to_string(n_bad_split) + ").\n");
  EKAT_REQUIRE_MSG(n_no_conv == 0,
    "Error! KESSLER: sub-cycling cannot converge (columns with dt0 <= 0 or >= " +
    std::to_string(max_pass) + " sub-cycles = " + std::to_string(n_no_conv) + ").\n");
}

} // namespace kessler_driver

#endif // DRIVER_KESSLER_VARIANTS_HPP
