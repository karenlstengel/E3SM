#include "eamxx_kessler_process_interface.hpp"
#include "physics/kessler/kessler_functions.hpp"

#include "share/property_checks/field_within_interval_check.hpp"
#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/field/field_utils.hpp"

#include "share/physics/physics_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"

#include <ekat_team_policy_utils.hpp>
#include <ekat_assert.hpp>
#include <ekat_units.hpp>

#include <array>

namespace scream
{
  using namespace kessler;
// =============================================================================
Kessler::Kessler (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Do nothing
}

// =============================================================================
void Kessler::create_requests ()
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  constexpr auto nondim = ekat::units::Units::nondimensional();
  static constexpr auto m2 = (m*m).rename("m2");

  m_grid = m_grids_manager->get_grid("physics");
  const auto& grid_name = m_grid->name();

  m_ncols = m_grid->get_num_local_dofs();
  m_nlevs = m_grid->get_num_vertical_levels();

  // Field layouts
  auto scalar2d     = m_grid->get_2d_scalar_layout();
  auto scalar3d_mid = m_grid->get_3d_scalar_layout(LEV);

  constexpr int ps = Pack::n;

  // Read-only inputs (not modified by Kessler)
  add_field<Required>("p_mid",          scalar3d_mid, Pa, grid_name, ps);
  add_field<Required>("pseudo_density", scalar3d_mid, Pa, grid_name, ps);

  // Updated fields: read at entry and written on exit
  add_field<Updated>("T_mid",         scalar3d_mid, K,    grid_name, ps);
  add_tracer<Updated>("qv", m_grid, kg/kg, ps);
  add_tracer<Updated>("qc", m_grid, kg/kg, ps);
  add_tracer<Updated>("qr", m_grid, kg/kg, ps);
  add_tracer<Updated>("qi", m_grid, kg/kg, ps);

  // Computed outputs
  add_field<Computed>("precl",  scalar2d,     m/s,    grid_name);
  add_field<Computed>("relhum", scalar3d_mid, nondim, grid_name, ps);

  // Initialise Kessler constants from physics constants.
  // cpair/rair are composition-independent dry-air scalar constants
  // (TODO: use composition-dependent fields if/when available).
  using C = physics::Constants<Real>;
  m_kd.lv    = C::LatVap.value;
  m_kd.pref  = C::P0.value / Real(100);  // convert Pa reference pressure to hPa
  m_kd.rhoqr = Real(1000);               // density of fresh liquid water, kg m-3
  m_kd.cpair = C::Cpair.value;
  m_kd.rair  = C::Rair.value;

  // Pulled from SHOC
  // Boundary flux fields for energy and mass conservation checks
  if (has_energy_fixer()) {
    add_field<Computed>("vapor_flux", scalar2d, kg/(m2*s), grid_name);
    add_field<Computed>("water_flux", scalar2d, m/s,       grid_name);
    add_field<Computed>("ice_flux",   scalar2d, m/s,       grid_name);
    add_field<Computed>("heat_flux",  scalar2d, W/m2,      grid_name);
  }
}

// =============================================================================
void Kessler::initialize_impl (const RunType /* run_type */)
{
  // Size the Kessler scratch workspace here (once, at case startup)
  // rather than lazily on the first run_impl call, so its one-time
  // device allocation doesn't land inside the timed simulation loop.
  m_kessler_workspace.init(m_ncols, m_nlevs);

  // Set up energy fixer fields.
  if (has_energy_fixer()) {
    // Set the boundary fluxes to 0.0 at the start of the run
    auto vapor_flux = get_field_out("vapor_flux").get_view<Real*>();
    auto water_flux = get_field_out("water_flux").get_view<Real*>();
    auto ice_flux   = get_field_out("ice_flux").get_view<Real*>();
    auto heat_flux  = get_field_out("heat_flux").get_view<Real*>();

    Kokkos::parallel_for("init_boundary_fluxes", m_ncols, KOKKOS_CLASS_LAMBDA(const int i) {
      vapor_flux(i) = 0.0;
      water_flux(i) = 0.0;
      ice_flux(i)   = 0.0;
      heat_flux(i)  = 0.0;
    });
  }
}

// =============================================================================
void Kessler::run_impl (const double dt)
{
  using PF  = scream::PhysicsFunctions<DefaultDevice>;
  using TPF = ekat::TeamPolicyFactory<KesslerFunc::KT::ExeSpace>;
  using KF  = KesslerFunc;

  const int ncols = m_ncols;
  const int nlevs = m_nlevs;

  // Get Pack views of required fields
  const auto T_mid_pack   = get_field_in("T_mid").get_view<const Pack**>();
  const auto p_mid_pack   = get_field_in("p_mid").get_view<const Pack**>();
  const auto pseudo_dens  = get_field_in("pseudo_density").get_view<const Pack**>();
  const auto qv_pack      = get_field_in("qv").get_view<const Pack**>();

  // Buffer scratch views
  const auto exner  = m_buffer.exner;
  const auto dz     = m_buffer.dz;
  const auto z_int  = m_buffer.z_int;
  const auto z_mid  = m_buffer.z_mid;
  const auto rho    = m_buffer.rho;
  const auto theta  = m_buffer.theta;

  // Physical constants needed inside device kernels
  using C = physics::Constants<Real>;
  const Real rair_val  = C::Rair.value;

  // Compute exner, dz, z, rho, theta using PhysicsFunctions
  const int nlev_packs = ekat::npack<Pack>(nlevs);
  const auto scan_policy = TPF::get_thread_range_parallel_scan_team_policy(ncols, nlev_packs);

  Kokkos::parallel_for(scan_policy,
    KOKKOS_LAMBDA (const KF::KT::MemberType& team) {
      const int i = team.league_rank();

      const auto p_mid_i    = ekat::subview(p_mid_pack, i);
      const auto qv_i       = ekat::subview(qv_pack,    i);
      const auto T_mid_i    = ekat::subview(T_mid_pack, i);
      const auto pseudo_i   = ekat::subview(pseudo_dens, i);
      const auto exner_i    = ekat::subview(exner,  i);
      const auto dz_i       = ekat::subview(dz,     i);
      const auto z_int_i    = ekat::subview(z_int,  i);
      const auto z_mid_i    = ekat::subview(z_mid,  i);
      const auto rho_i      = ekat::subview(rho,    i);
      const auto theta_i    = ekat::subview(theta,  i);

      // Exner function: pi = (p/p0)^(R/cp)
      PF::exner_function<Pack>(team, p_mid_i, exner_i);

      // Layer thicknesses and heights (z_int(i, nlevs) = 0 at surface)
      PF::calculate_dz(team, pseudo_i, p_mid_i, T_mid_i, qv_i, dz_i);
      const Real z_surf = 0;
      team.team_barrier();
      PF::calculate_z_int(team, nlevs, dz_i, z_surf, z_int_i);
      team.team_barrier();
      PF::calculate_z_mid(team, nlevs, z_int_i, z_mid_i);
      team.team_barrier();

      // Dry air density: rho = p / (Rair * T_mid)
      // Note: a more accurate formulation would use virtual temperature
      // to account for moisture effects on density.
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_packs),
        [&](int k) {
          rho_i(k) = p_mid_i(k) / (rair_val * T_mid_i(k));
        });
      team.team_barrier();

      // Potential temperature: theta = T / pi
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_packs),
        [&](int k) {
          theta_i(k) = T_mid_i(k) / exner_i(k);
        });
    });

  // Scalarize Pack views for kessler_run (no SIMD over levels)
  auto T_mid_upd  = get_field_out("T_mid").get_view<Pack**>();
  auto qv_upd     = get_field_out("qv").get_view<Pack**>();
  auto qc_upd     = get_field_out("qc").get_view<Pack**>();
  auto qr_upd     = get_field_out("qr").get_view<Pack**>();
  auto precl_out  = get_field_out("precl").get_view<Real*>();
  auto relhum_out = get_field_out("relhum").get_view<Pack**>();

  // cpair and rair are supplied to kessler_run as scalar constants via
  // m_kd (set once in create_requests), so no per-call fill kernel or
  // 2D constant views are needed here.

  // Scalarize 2D Pack views to Scalar views for the kessler kernel
  auto z_mid_s  = ekat::scalarize(z_mid);
  auto rho_s    = ekat::scalarize(rho);
  auto exner_s  = ekat::scalarize(exner);
  auto theta_s  = ekat::scalarize(theta);
  auto qv_s     = ekat::scalarize(qv_upd);
  auto qc_s     = ekat::scalarize(qc_upd);
  auto qr_s     = ekat::scalarize(qr_upd);
  auto relhum_s = ekat::scalarize(relhum_out);

  // pk in Kessler is the Exner function (pi = T/theta), already computed
  // in exner_s.  z is the geopotential height, in z_mid_s.
  // Convention: lyr_surf=nlevs-1 (surface at last C++ index in EAMxx),
  //             lyr_toa=0.
  const int lyr_surf = nlevs - 1;
  const int lyr_toa  = 0;

  // m_kessler_workspace was already sized in initialize_impl(); the
  // init() call kessler_run makes internally is then just a cheap
  // dimension check, not a device allocation.
  KF::kessler_run(ncols, nlevs, Real(dt),
                  lyr_surf, lyr_toa, m_kd,
                  m_kessler_workspace, rho_s, z_mid_s, exner_s,
                  theta_s, qv_s, qc_s, qr_s,
                  precl_out, relhum_s);

  // Update T_mid from modified theta: T = theta * pi
  Kokkos::parallel_for("kessler_theta_to_T",
    Kokkos::MDRangePolicy<KF::KT::ExeSpace, Kokkos::Rank<2>>(
      {0, 0}, {ncols, nlevs}),
    KOKKOS_LAMBDA(const int col, const int k) {
      ekat::scalarize(T_mid_upd)(col, k) = theta_s(col, k) * exner_s(col, k);
    });

  Kokkos::fence();

  // KF:: update?

  if (has_energy_fixer()) {
    // Set the boundary fluxes to 0.0 at the start of the run
    auto vapor_flux = get_field_out("vapor_flux").get_view<Real*>();
    auto water_flux = get_field_out("water_flux").get_view<Real*>();
    auto ice_flux   = get_field_out("ice_flux").get_view<Real*>();
    auto heat_flux  = get_field_out("heat_flux").get_view<Real*>();

    Kokkos::parallel_for("init_boundary_fluxes", m_ncols, KOKKOS_CLASS_LAMBDA(const int i) {
      vapor_flux(i) = 0.0;
      water_flux(i) = 0.0;
      ice_flux(i)   = 0.0;
      heat_flux(i)  = 0.0;
    });
  }

}

// =============================================================================
void Kessler::finalize_impl ()
{
  // Do nothing
}

// =============================================================================
size_t Kessler::requested_buffer_size_in_bytes () const
{
  const int nlev_packs  = ekat::npack<Pack>(m_nlevs);
  const int nlevi_packs = ekat::npack<Pack>(m_nlevs + 1);
  return Buffer::num_2d_mid * m_ncols * nlev_packs  * sizeof(Pack)
       + Buffer::num_2d_int * m_ncols * nlevi_packs * sizeof(Pack);
}

// =============================================================================
void Kessler::init_buffers (const ATMBufferManager& buffer_manager)
{
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
                   "Error! Buffer size not sufficient for Kessler.\n");

  Pack* mem = reinterpret_cast<Pack*>(buffer_manager.get_memory());
  const int nlev_packs  = ekat::npack<Pack>(m_nlevs);
  const int nlevi_packs = ekat::npack<Pack>(m_nlevs + 1);

  uview_2d* mid_ptrs[Buffer::num_2d_mid] = {
    &m_buffer.exner, &m_buffer.dz, &m_buffer.z_mid, &m_buffer.rho, &m_buffer.theta
  };
  for (int i = 0; i < Buffer::num_2d_mid; ++i) {
    *mid_ptrs[i] = uview_2d(mem, m_ncols, nlev_packs);
    mem += mid_ptrs[i]->size();
  }

  uview_2d_int* int_ptrs[Buffer::num_2d_int] = { &m_buffer.z_int };
  for (int i = 0; i < Buffer::num_2d_int; ++i) {
    *int_ptrs[i] = uview_2d_int(mem, m_ncols, nlevi_packs);
    mem += int_ptrs[i]->size();
  }

  const size_t used =
    (reinterpret_cast<Real*>(mem) - buffer_manager.get_memory()) * sizeof(Real);
  EKAT_REQUIRE_MSG(used == requested_buffer_size_in_bytes(),
                   "Error! Buffer accounting mismatch in Kessler.\n");
}

} // namespace scream
