#ifndef DRIVER_KESSLER_PREPROCESS_HPP
#define DRIVER_KESSLER_PREPROCESS_HPP

/*
 * The pre-processing KesslerMicrophysics::run_impl does before kessler_run
 * (pk, theta, dz, rho, then z_int and z_mid), in two versions, so
 * driver_kessler.cpp --preprocess can check that they agree and time them
 * (performance_suggestions.md, item I4):
 *
 *   preprocess_old: copied from run_impl (2026-10-28): a RangePolicy over
 *     every (column, level) scalar, then a separate team kernel for the
 *     z_int scan and z_mid.
 *   preprocess_new: one team per column (scan policy): the same
 *     PhysicsFunctions formulas over whole packs, then the z_int scan and
 *     z_mid in the same kernel, as SHOC's preprocess does. Padding lanes get
 *     benign inputs so nothing divides by zero.
 */

#include "kessler_functions.hpp"

#include <ekat_subview_utils.hpp>
#include <Kokkos_Core.hpp>

namespace kessler_driver {

template <typename DeviceT>
struct PreprocessViews {
  using Pack = ekat::Pack<scream::Real, SCREAM_PACK_SIZE>;
  using KT   = ekat::KokkosTypes<DeviceT>;
  template <typename S> using view_2d = typename KT::template view_2d<S>;

  view_2d<const Pack> T_mid, p_mid, pseudo_density, qv;   // inputs
  view_2d<Pack> pk, theta, dz, rho, z_mid, z_int;         // outputs
};

template <typename DeviceT>
void preprocess_old (const int ncols, const int nlevs, const PreprocessViews<DeviceT>& v)
{
  using Real = scream::Real;
  using KT   = ekat::KokkosTypes<DeviceT>;
  using PF   = scream::PhysicsFunctions<DeviceT>;
  using PC   = scream::physics::Constants<Real>;
  using Pack = typename PreprocessViews<DeviceT>::Pack;
  using TPF  = ekat::TeamPolicyFactory<typename KT::ExeSpace>;

  const auto T_mid = v.T_mid; const auto p_mid = v.p_mid;
  const auto pseudo_density = v.pseudo_density; const auto qv = v.qv;
  const auto pk = v.pk; const auto theta = v.theta; const auto dz = v.dz; const auto rho = v.rho;
  const auto z_mid = v.z_mid; const auto z_int = v.z_int;
  const Real inv_ggr = 1/(PC::gravit.value);

  Kokkos::parallel_for(
      "Kessler_preprocess", Kokkos::RangePolicy<typename KT::ExeSpace>(0, ncols * nlevs),
      KOKKOS_LAMBDA(const int i) {
        const int icol = i / nlevs;
        const int klev = i % nlevs;
        pk(icol, klev / Pack::n)[klev % Pack::n]    = PF::exner_function(p_mid(icol, klev / Pack::n)[klev % Pack::n]);
        theta(icol, klev / Pack::n)[klev % Pack::n] = PF::calculate_theta_from_T(T_mid(icol, klev / Pack::n)[klev % Pack::n],p_mid(icol, klev / Pack::n)[klev % Pack::n]);
        dz(icol, klev / Pack::n)[klev % Pack::n]  = PF::calculate_dz(pseudo_density(icol, klev / Pack::n)[klev % Pack::n], p_mid(icol, klev / Pack::n)[klev % Pack::n], T_mid(icol, klev / Pack::n)[klev % Pack::n], qv(icol, klev / Pack::n)[klev % Pack::n]);
        rho(icol, klev / Pack::n)[klev % Pack::n] = inv_ggr*( pseudo_density(icol, klev / Pack::n)[klev % Pack::n]/dz(icol, klev / Pack::n)[klev % Pack::n] );
      });

  const int nlev_packs   = ekat::npack<Pack>(nlevs);
  const auto scan_policy = TPF::get_thread_range_parallel_scan_team_policy(ncols, nlev_packs);
  Kokkos::parallel_for(scan_policy, KOKKOS_LAMBDA (const typename KT::MemberType& team) {
    const int i = team.league_rank();
    auto z_mid_i = ekat::subview(z_mid, i);
    auto dz_i    = ekat::subview(dz, i);
    auto z_int_i = ekat::subview(z_int, i);
    Real zurf  = 0.0;
    PF::calculate_z_int(team, nlevs, dz_i, zurf, z_int_i);
    team.team_barrier();
    PF::calculate_z_mid(team, nlevs, z_int_i, z_mid_i);
    team.team_barrier();
  });
}

template <typename DeviceT>
void preprocess_new (const int ncols, const int nlevs, const PreprocessViews<DeviceT>& v)
{
  using Real    = scream::Real;
  using KT      = ekat::KokkosTypes<DeviceT>;
  using PF      = scream::PhysicsFunctions<DeviceT>;
  using PC      = scream::physics::Constants<Real>;
  using Pack    = typename PreprocessViews<DeviceT>::Pack;
  using IntPack = ekat::Pack<scream::Int, SCREAM_PACK_SIZE>;
  using TPF     = ekat::TeamPolicyFactory<typename KT::ExeSpace>;

  const auto T_mid = v.T_mid; const auto p_mid = v.p_mid;
  const auto pseudo_density = v.pseudo_density; const auto qv = v.qv;
  const auto pk = v.pk; const auto theta = v.theta; const auto dz = v.dz; const auto rho = v.rho;
  const auto z_mid = v.z_mid; const auto z_int = v.z_int;
  const Real inv_ggr = 1/(PC::gravit.value);
  const Real p0      = PC::P0.value;

  const int  nlev_packs  = ekat::npack<Pack>(nlevs);
  const auto scan_policy = TPF::get_thread_range_parallel_scan_team_policy(ncols, nlev_packs);

  Kokkos::parallel_for("Kessler_preprocess", scan_policy,
    KOKKOS_LAMBDA (const typename KT::MemberType& team) {
      const int i = team.league_rank();

      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_packs), [&] (const int k) {
        // Padding lanes (nlevs not a multiple of Pack::n) get benign inputs,
        // so nothing divides by zero under the fpe build; their results
        // are never read.
        const auto valid = ekat::range<IntPack>(k*Pack::n) < nlevs;
        Pack p    = p_mid(i,k);           p.set(!valid, p0);
        Pack pdel = pseudo_density(i,k);  pdel.set(!valid, Real(1));
        Pack T    = T_mid(i,k);           T.set(!valid, Real(1));
        Pack q    = qv(i,k);              q.set(!valid, Real(0));

        // calculate_theta_from_T(T, p) is T / exner_function(p); reuse the
        // Exner value instead of computing its pow a second time.
        const Pack exner = PF::exner_function(p);
        pk(i,k)    = exner;
        theta(i,k) = T / exner;
        dz(i,k)    = PF::calculate_dz(pdel, p, T, q);
        rho(i,k)   = inv_ggr * (pdel / dz(i,k));
      });
      team.team_barrier();

      const auto dz_i    = ekat::subview(dz,    i);
      const auto z_int_i = ekat::subview(z_int, i);
      const auto z_mid_i = ekat::subview(z_mid, i);
      PF::calculate_z_int(team, nlevs, dz_i, Real(0), z_int_i);
      team.team_barrier();
      PF::calculate_z_mid(team, nlevs, z_int_i, z_mid_i);
    });
}

} // namespace kessler_driver

#endif // DRIVER_KESSLER_PREPROCESS_HPP
