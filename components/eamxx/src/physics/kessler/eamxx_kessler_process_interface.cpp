#include "eamxx_kessler_process_interface.hpp"
// #include "kessler_eamxx_bridge.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/field/field_utils.hpp"

#include "share/physics/physics_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
// #include "share/physics/eamxx_common_physics_functions_impls.hpp"
#include "share/util/eamxx_timing.hpp"

#include <ekat_assert.hpp>
#include <ekat_units.hpp>

#ifdef EAMXX_HAS_PYTHON
#include "share/atm_process/atmosphere_process_pyhelpers.hpp"
#endif

#include <array>

namespace scream
{
  using namespace kessler;
// =========================================================================================
//  Inputs (these are inherited from AtomoshpereProcess which means we can use the same logger):
//      comm - an EKAT communication group
//      params - a parameter list of options for the process.
//  Outputs:
//      None

KesslerMicrophysics::KesslerMicrophysics (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here usually
  m_atm_logger->info("[EAMxx] Kessler processes constructor");

}

// =========================================================================================
//  Inputs:
//      None
//  Outputs:
//      None

void KesslerMicrophysics::create_requests () //set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using PC = scream::physics::Constants<Real>;

  pref    = PC::P0.value / 100.0;     // Reference pressure; pref_in (Pa -> hPa)
  latvap  = PC::LatVap.value; // Latent heat of vaporization; lv_in
  rhoqr   = PC::RHOW.value;   // rhoqr_in
  gravity = PC::gravit.value; // gravitational acceleration
  Cpair   = PC::Cpair.value; // Specific heat of dry air at constant pressure
  Rair    = PC::Rair.value;  // Gas constant of dry air

  using namespace ekat::units;
  // using namespace ShortFieldTagsNames;

  m_atm_logger->info("[EAMxx] Kessler processes set grids");

  static constexpr auto m2 = (m*m).rename("m2");
  static constexpr auto m3 = (m*m*m).rename("m3");
  static constexpr auto s2 = (s*s).rename("s2");
  constexpr auto nondim    = ekat::units::Units::nondimensional();
  constexpr int packize  = Pack::n;


  // specify which grid to use
  m_grid = m_grids_manager->get_grid("physics");

  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  // Layout for 2D (1d horiz X 1d vertical) variable
  FieldLayout scalar2d_layout { {ShortFieldTagsNames::COL}, {m_num_cols} };

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces
  FieldLayout scalar3d_layout_mid { {ShortFieldTagsNames::COL,ShortFieldTagsNames::LEV}, {m_num_cols,m_num_levs} };
  // Layout for 3D variable defined at interfaces (nlevs+1 values per column)
  FieldLayout scalar3d_layout_int { {ShortFieldTagsNames::COL,ShortFieldTagsNames::ILEV}, {m_num_cols,m_num_levs+1} };
 
  // Fields to use for Kessler microphysics

  // From Field Manager 
  // real(kind_phys),  intent(inout) :: qv(:,:)    ! Water vapor mixing ratio wrt dry air (kg/kg)
  // real(kind_phys),  intent(inout) :: qc(:,:)    ! Cloud water mixing ratio wrt dry air (kg/kg)
  // real(kind_phys),  intent(inout) :: qr(:,:)    ! Rain water mixing ratio wrt dry air (kg/kg)

  add_field<Updated>("T_mid",               scalar3d_layout_mid, K,     grid_name, packize);
  add_field<Computed>("T_mid_prev",          scalar3d_layout_mid, K,     grid_name, packize);
  add_field<Required>("p_mid",              scalar3d_layout_mid, Pa,    grid_name, packize);
  add_field<Required>("pseudo_density",     scalar3d_layout_mid, Pa,    grid_name, packize);
  // add_field<Required>("pseudo_density_dry", scalar3d_layout_mid, Pa,    grid_name, packize);
  add_tracer<Updated>("qv",                 m_grid,              kg/kg,            packize);
  add_tracer<Updated>("qc",                 m_grid,              kg/kg,            packize);
  add_tracer<Updated>("qr",                 m_grid,              kg/kg,            packize);
  add_tracer<Updated>("qi",                 m_grid,              kg/kg,          packize);
  // Locally needed 
  // real(kind_phys),  intent(in)    :: rho(:,:)   ! Dry air density (kg/m^3)
  // real(kind_phys),  intent(in)    :: z(:,:)     ! Heights of thermo. levels (m)
  // real(kind_phys),  intent(in)    :: pk(:,:)    ! Exner function (p/p0)**(R/cp)

  // real(kind_phys),  intent(inout) :: theta(:,:) ! Potential temperature (K)

  // real(kind_phys),  intent(out)   :: precl(:)   ! Precipitation rate (m_water / s)
  // real(kind_phys),  intent(out)   :: relhum(:,:)! Relative humidity in percent

  add_field<Computed>("rho",    scalar3d_layout_mid, kg/m3,  grid_name, packize);
  add_field<Computed>("dz",     scalar3d_layout_mid, m,      grid_name, packize);
  add_field<Computed>("pk",     scalar3d_layout_mid, nondim, grid_name, packize);
  add_field<Computed>("theta",  scalar3d_layout_mid, K,      grid_name, packize);
  add_field<Computed>("precl",  scalar2d_layout,     m/s,    grid_name, packize);
  add_field<Computed>("relhum", scalar3d_layout_mid, nondim, grid_name, packize);

  // Other stuff we need to use
  add_field<Computed>("st_energy", scalar3d_layout_mid, J/kg,  grid_name, packize); // drytatic_energy
  add_field<Updated>("phis",       scalar2d_layout, m2/s2,     grid_name, packize); // geopotential height of surface

  // Tendencies 
  add_field<Computed>("T_mid_tend", scalar3d_layout_mid, K/s, grid_name, packize);
  
  // Pulled from SHOC
  // Boundary flux fields for energy and mass conservation checks
  if (has_energy_fixer()) {
    add_field<Computed>("vapor_flux", scalar2d_layout, kg/(m2*s), grid_name);
    add_field<Computed>("water_flux", scalar2d_layout, m/s,       grid_name);
    add_field<Computed>("ice_flux",   scalar2d_layout, m/s,       grid_name);
    add_field<Computed>("heat_flux",  scalar2d_layout, W/m2,      grid_name);
  }

  add_field<Computed>("z_mid",  scalar3d_layout_mid, m, grid_name, packize);
  add_field<Computed>("z_int",  scalar3d_layout_int, m, grid_name, packize);

}

// =========================================================================================
//  Inputs:
//      run_type - run type, either initial or restart
//  Outputs:
//      None

void KesslerMicrophysics::initialize_impl (const RunType /* run_type */)
{
  m_atm_logger->info("[EAMxx] kessler processes initialize_impl: ");

  // Set post condition checks (these are taken care of by HOMME)
  // <scheme>check_energy_zero_fluxes</scheme>
  // <scheme>check_energycaling</scheme>
  // <scheme>check_energy_chng</scheme>

  add_invariant_check<FieldWithinIntervalCheck>(get_field_out("T_mid"),m_grid,100.0,500.0,false);
  add_invariant_check<FieldWithinIntervalCheck>(get_field_out("qv"),m_grid,1e-13,0.2,true);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("qc"),m_grid,0.0,0.1,true);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("qr"),m_grid,0.0,0.1,true);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("qi"),m_grid,0.0,0.1,true);
  add_postcondition_check<FieldLowerBoundCheck>(get_field_out("precl"),m_grid,0.0,true);

  if (has_energy_fixer()) {
    // Set the boundary fluxes to 0.0 at the start of the run
    auto vapor_flux = get_field_out("vapor_flux").get_view<Real*>();
    auto water_flux = get_field_out("water_flux").get_view<Real*>();
    auto ice_flux   = get_field_out("ice_flux").get_view<Real*>();
    auto heat_flux  = get_field_out("heat_flux").get_view<Real*>();

    Kokkos::parallel_for("init_boundary_fluxes", m_num_cols, KOKKOS_CLASS_LAMBDA(const int i) {
      vapor_flux(i) = 0.0;
      water_flux(i) = 0.0;
      ice_flux(i)   = 0.0;
      heat_flux(i)  = 0.0;
    });
  }

}

// =========================================================================================

// run_impl is called every timestep and where all of the physics happens
// Inputs:
//    - dt - the timestep for the current run step
//  Outputs:
//      None

void KesslerMicrophysics::run_impl (const double dt )
{
  using KMF = kessler::KesslerMicrophysicsFunctions<Real, DefaultDevice>;

  // Pull in variables 
  auto T_mid              = get_field_out("T_mid").get_view<Pack**>();
  auto p_mid              = get_field_in("p_mid").get_view<const Pack**>();
  auto pseudo_density     = get_field_in("pseudo_density").get_view<const Pack**>();
  auto qv                 = get_field_out("qv").get_view<Pack**>();
  auto qc                 = get_field_out("qc").get_view<Pack**>();
  auto qr                 = get_field_out("qr").get_view<Pack**>();
  auto rho                = get_field_out("rho").get_view<Pack**>();
  auto dz                 = get_field_out("dz").get_view<Pack**>();
  auto pk                 = get_field_out("pk").get_view<Pack**>();
  auto theta              = get_field_out("theta").get_view<Pack**>();
  auto precl              = get_field_out("precl").get_view<Real*>();
  auto relhum             = get_field_out("relhum").get_view<Pack**>();
  // Need phis for z_mid computation
  auto phis               = get_field_out("phis").get_view<Real*>();
  // <scheme>kessler_update</scheme> // updates st_energy & temperature related things
  auto T_mid_prev        = get_field_out("T_mid_prev").get_view<Pack**>();
  auto T_mid_tend        = get_field_out("T_mid_tend").get_view<Pack**>();
  auto st_energy         = get_field_out("st_energy").get_view<Pack**>();

  // Geopotential height and interface height for kessler_run
  auto z_mid              = get_field_out("z_mid").get_view<Pack**>();
  auto z_int              = get_field_out("z_int").get_view<Pack**>();

  // EAMxx orders levels top-down (index 0 = TOA, index m_num_levs-1 =
  // surface), so lyr_surf/lyr_toa (see kessler_run's Convention note
  // in kessler_functions_impl.hpp) are the last/first level indices.
  const int lyr_surf = m_num_levs - 1;
  const int lyr_toa  = 0;

  m_atm_logger->info("[EAMxx] kessler run_impl: ");
  
  // NOTE: EAMxx uses a mass-weighted vertical discretization, which is different 
  // than CAM's volume-weighted vertical coordinate. 
  // You will need to derive rho and z from the pseudo_density used in EAMxx.  
  // The change from volume to mass-weighted was to reduce the amount of data movement 
  // and improve numerical accuracy which will help with running reduced precision in the physics parameterizations.
 
  // do the conversion and PF::exner_function, PF::calculate_theta_from_T

  using PF           = scream::PhysicsFunctions<DefaultDevice>;
  using PC           = scream::physics::Constants<Real>;
  const Real inv_ggr = 1/(PC::gravit.value);

  int nlevs = m_num_levs; // local var, to avoid accessing *this.
  const int ncols = m_num_cols;

  start_timer("EAMxx::kessler::run::preprocess");
  Kokkos::parallel_for(
      "Kessler_preprocess", KT::RangePolicy(0, m_num_cols * nlevs),
      KOKKOS_CLASS_LAMBDA(const int i) {
        const int icol = i / nlevs;
        const int klev = i % nlevs;
        pk(icol, klev / Pack::n)[klev % Pack::n]    = PF::exner_function(p_mid(icol, klev / Pack::n)[klev % Pack::n]);
        theta(icol, klev / Pack::n)[klev % Pack::n] = PF::calculate_theta_from_T(T_mid(icol, klev / Pack::n)[klev % Pack::n],p_mid(icol, klev / Pack::n)[klev % Pack::n]);

        // Vertical layer thickness
        dz(icol, klev / Pack::n)[klev % Pack::n]  = PF::calculate_dz(pseudo_density(icol, klev / Pack::n)[klev % Pack::n], p_mid(icol, klev / Pack::n)[klev % Pack::n], T_mid(icol, klev / Pack::n)[klev % Pack::n], qv(icol, klev / Pack::n)[klev % Pack::n]);
        rho(icol, klev / Pack::n)[klev % Pack::n] = inv_ggr*( pseudo_density(icol, klev / Pack::n)[klev % Pack::n]/dz(icol, klev / Pack::n)[klev % Pack::n] );
      }
    );

  // Compute z_mid and z_int into params_in (needed by kessler_run for heights)
  // calculate_z_int() contains a team-level parallelcan, which requires a special policy
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;
  const int nlev_packs   = ekat::npack<Pack>(m_num_levs);
  const auto scan_policy = TPF::get_thread_range_parallel_scan_team_policy(m_num_cols, nlev_packs);

  Kokkos::parallel_for(scan_policy, KOKKOS_CLASS_LAMBDA (const KT::MemberType& team) {
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

  // <!-- MPAS and SE specific scaling of temperature for enforcing energy consistency:
  //       First, calculate the scaling based off cp_or_cv_dycore (from cam_thermo_water_update)
  //       Then, perform the temperature and temperature tendency scaling,
  //       and apply tendencies resulting from such adjustment -->
  // <scheme>check_energycaling</scheme> 
    // real(kind_phys),    intent(out)    :: flx_vap(:)     ! boundary flux of vapor [kg m-2 s-1]
    // real(kind_phys),    intent(out)    :: flx_cnd(:)     ! boundary flux of liquid+ice (precip?) [m s-1]
    // real(kind_phys),    intent(out)    :: flx_ice(:)     ! boundary flux of ice (snow?) [m s-1]
    // real(kind_phys),    intent(out)    :: flxen(:)     ! boundary flux of sensible heat [W m-2]
    // sets all four of the fluxes = 0.0
  if (has_energy_fixer()) {
    auto vapor_flux = get_field_out("vapor_flux").get_view<Real*>();
    auto water_flux = get_field_out("water_flux").get_view<Real*>();
    auto ice_flux   = get_field_out("ice_flux").get_view<Real*>();
    auto heat_flux  = get_field_out("heat_flux").get_view<Real*>();

    Kokkos::parallel_for("check_energy_scaling", m_num_cols, KOKKOS_CLASS_LAMBDA(const int i) {
      vapor_flux(i) = 0.0;
      water_flux(i) = 0.0;
      ice_flux(i)   = 0.0;
      heat_flux(i)  = 0.0;
    });
  }
  stop_timer("EAMxx::kessler::run::preprocess");
  m_atm_logger->info("[EAMxx] kessler run_impl: done with pre-processing");

  // Save the pre-step temperature and zero the tendency accumulator
  // before running the main Kessler update.
  start_timer("EAMxx::kessler::run::kessler_update_timestep_init");
  KMF::kessler_update_timestep_init(ncols, nlevs, T_mid, T_mid_prev, T_mid_tend);
  Kokkos::fence();
  stop_timer("EAMxx::kessler::run::kessler_update_timestep_init");

  start_timer("EAMxx::kessler::run::kessler_run");
  // Workspace views into this call's persistent scratch (see Buffer,
  // populated once in init_buffers()).
  KMF::Workspace ws;
  ws.r            = m_buffer.r;
  ws.rhalf        = m_buffer.rhalf;
  ws.velqr        = m_buffer.velqr;
  ws.sed          = m_buffer.sed;
  ws.pc           = m_buffer.pc;
  ws.f5           = m_buffer.f5;
  ws.dt0          = m_buffer.dt0;
  ws.mask         = m_buffer.mask;
  ws.time_counter = m_buffer.time_counter;
  ws.precl_acc    = m_buffer.precl_acc;

  KMF::kessler_run(ncols, nlevs, Real(dt),
                  lyr_surf, lyr_toa, rhoqr, latvap, pref,
                  Cpair, Rair, rho, z_mid, pk,
                  theta, qv, qc, qr,
                  precl, relhum, ws);

  Kokkos::fence();
  stop_timer("EAMxx::kessler::run::kessler_run");

  // Back out the temperature tendency due to Kessler from theta/exner:
  // ttend = (theta*exner - T_mid_prev) / dt.
  start_timer("EAMxx::kessler::run::kessler_update_run");
  KMF::kessler_update_run(ncols, nlevs, Real(dt), theta, pk, T_mid_prev, T_mid_tend);
  Kokkos::fence();
  stop_timer("EAMxx::kessler::run::kessler_update_run");

  // Apply that tendency to get the updated T_mid: T_mid_prev + ttend*dt is
  // algebraically theta*exner, but going through the tendency (rather than
  // recomputing theta*exner directly) keeps T_mid consistent with the
  // T_mid_tend diagnostic above and is what the rest of the AD reads.
  // start_timer("EAMxx::kessler::run::kessler_apply_T_tend");
  // {
  //   const Real dt_r = Real(dt);
  //   Kokkos::parallel_for("kessler_apply_T_tend",
  //     Kokkos::MDRangePolicy<KMF::KT::ExeSpace, Kokkos::Rank<2>>(
  //       {0, 0}, {ncols, nlevs}),
  //     KOKKOS_LAMBDA(const int col, const int k) {
  //       ekat::scalarize(T_mid_upd)(col, k) = t_prev(col, k) + t_tend(col, k) * dt_r;
  //     });
  //   Kokkos::fence();
  // }
  // stop_timer("EAMxx::kessler::run::kessler_apply_T_tend");

  start_timer("EAMxx::kessler::run::kessler_update_timestep_final");
  KMF::kessler_update_timestep_final(ncols, nlevs, gravity,
                                     Cpair, T_mid, z_mid,
                                     phis, st_energy);
  Kokkos::fence();
  stop_timer("EAMxx::kessler::run::kessler_update_timestep_final");

}

// =========================================================================================
//  Inputs:
//      None
//  Outputs:
//      None

void KesslerMicrophysics::finalize_impl()
{
  m_atm_logger->info("[EAMxx] Kessler processes clean up.");
}


// =========================================================================================
// Buffer/workspace management: one ATMBufferManager allocation for
// kessler_run's persistent scratch (Buffer, see the .hpp), carved into
// unmanaged views the same way P3Microphysics/SHOCMacrophysics carve
// their own Buffer from buffer_manager.get_memory().

size_t KesslerMicrophysics::requested_buffer_size_in_bytes() const
{
  const Int nlev_packs = ekat::npack<Pack>(m_num_levs);

  return Buffer::num_1d_scalar * m_num_cols * sizeof(Real) +
         Buffer::num_2d_vector * m_num_cols * nlev_packs * sizeof(Pack);
}

// =========================================================================================
void KesslerMicrophysics::init_buffers(const ATMBufferManager& buffer_manager)
{
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
                    "Error! Buffers size not sufficient.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  // 1d scalar views (per-column sub-cycle bookkeeping)
  using scalar_1d_view_t = decltype(m_buffer.dt0);
  scalar_1d_view_t* _1d_scalar_view_ptrs[Buffer::num_1d_scalar] = {
    &m_buffer.dt0, &m_buffer.mask, &m_buffer.time_counter, &m_buffer.precl_acc
  };
  for (int i=0; i<Buffer::num_1d_scalar; ++i) {
    *_1d_scalar_view_ptrs[i] = scalar_1d_view_t(mem, m_num_cols);
    mem += _1d_scalar_view_ptrs[i]->size();
  }

  // 2d packed views (per-column, per-level scratch)
  Pack* p_mem = reinterpret_cast<Pack*>(mem);
  const Int nlev_packs = ekat::npack<Pack>(m_num_levs);

  using pack_2d_view_t = decltype(m_buffer.r);
  pack_2d_view_t* _2d_pack_view_ptrs[Buffer::num_2d_vector] = {
    &m_buffer.r, &m_buffer.rhalf, &m_buffer.velqr, &m_buffer.sed, &m_buffer.pc, &m_buffer.f5
  };
  for (int i=0; i<Buffer::num_2d_vector; ++i) {
    *_2d_pack_view_ptrs[i] = pack_2d_view_t(p_mem, m_num_cols, nlev_packs);
    p_mem += _2d_pack_view_ptrs[i]->size();
  }

  size_t used_mem = (reinterpret_cast<Real*>(p_mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(),
                    "Error! Used memory != requested memory for KesslerMicrophysics.");
}

} // namespace scream