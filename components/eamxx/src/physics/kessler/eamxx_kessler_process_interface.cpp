#include "eamxx_kessler_process_interface.hpp"
#include "kessler_eamxx_bridge.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "share/field/field_utils.hpp"

#include "share/physics/physics_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
// #include "share/physics/eamxx_common_physics_functions_impls.hpp"

#include <ekat_assert.hpp>
#include <ekat_units.hpp>

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
//      grids_manager - The grids manager used by the dynamcial core
//  Outputs:
//      None

void KesslerMicrophysics::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  // using namespace ShortFieldTagsNames;

  m_atm_logger->info("[EAMxx] Kessler processes set grids");

  const auto s2           = pow(s,2);
  const auto m2           = pow(m,2);
  const auto m3           = pow(m,3);
  constexpr auto nondim   = ekat::units::Units::nondimensional();
  constexpr int pack_size = Spack::n;


  // specify which grid to use
  m_grid = grids_manager->get_grid("physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  // Layout for 2D (1d horiz X 1d vertical) variable
  FieldLayout scalar2d_layout { {ShortFieldTagsNames::COL}, {m_num_cols} };

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces
  FieldLayout scalar3d_layout_mid { {ShortFieldTagsNames::COL,ShortFieldTagsNames::LEV}, {m_num_cols,m_num_levs} };
 
  // Fields to use for Kessler microphysics

  // From Field Manager 
  // real(kind_phys),  intent(inout) :: qv(:,:)    ! Water vapor mixing ratio wrt dry air (kg/kg)
  // real(kind_phys),  intent(inout) :: qc(:,:)    ! Cloud water mixing ratio wrt dry air (kg/kg)
  // real(kind_phys),  intent(inout) :: qr(:,:)    ! Rain water mixing ratio wrt dry air (kg/kg)

  add_field<Updated>("T_mid",               scalar3d_layout_mid, K,     grid_name, pack_size);
  add_field<Computed>("T_mid_prev",          scalar3d_layout_mid, K,     grid_name, pack_size);
  add_field<Required>("p_mid",              scalar3d_layout_mid, Pa,    grid_name, pack_size);
  add_field<Required>("pseudo_density",     scalar3d_layout_mid, Pa,    grid_name, pack_size);
  // add_field<Required>("pseudo_density_dry", scalar3d_layout_mid, Pa,    grid_name, pack_size);
  add_tracer<Updated>("qv",                 m_grid,              kg/kg,            pack_size);
  add_tracer<Updated>("qc",                 m_grid,              kg/kg,            pack_size);
  add_tracer<Updated>("qr",                 m_grid,              kg/kg,            pack_size);
  add_tracer<Updated>("qi",                 m_grid,              kg/kg,          pack_size);
  // Locally needed 
  // real(kind_phys),  intent(in)    :: rho(:,:)   ! Dry air density (kg/m^3)
  // real(kind_phys),  intent(in)    :: z(:,:)     ! Heights of thermo. levels (m)
  // real(kind_phys),  intent(in)    :: pk(:,:)    ! Exner function (p/p0)**(R/cp)

  // real(kind_phys),  intent(inout) :: theta(:,:) ! Potential temperature (K)

  // real(kind_phys),  intent(out)   :: precl(:)   ! Precipitation rate (m_water / s)
  // real(kind_phys),  intent(out)   :: relhum(:,:)! Relative humidity in percent

  add_field<Computed>("rho",    scalar3d_layout_mid, kg/m3,  grid_name, pack_size);
  add_field<Computed>("dz",     scalar3d_layout_mid, m,      grid_name, pack_size);
  add_field<Computed>("pk",     scalar3d_layout_mid, nondim, grid_name, pack_size);
  add_field<Computed>("theta",  scalar3d_layout_mid, K,      grid_name, pack_size);
  add_field<Computed>("precl",  scalar2d_layout,     m/s,    grid_name, pack_size);
  add_field<Computed>("relhum", scalar3d_layout_mid, nondim, grid_name, pack_size);

  // Other stuff we need to use
  add_field<Computed>("st_energy", scalar3d_layout_mid, J/kg,  grid_name, pack_size); // dry_static_energy
  add_field<Updated>("phis",       scalar2d_layout, m2/s2,     grid_name, pack_size); // geopotential height of surface

  // Tendencies 
  add_field<Computed>("T_mid_tend", scalar3d_layout_mid, K/s, grid_name, pack_size);
  
  // Pulled from SHOC
  // Boundary flux fields for energy and mass conservation checks
  if (has_energy_fixer()) {
    add_field<Computed>("vapor_flux", scalar2d_layout, kg/(m2*s), grid_name);
    add_field<Computed>("water_flux", scalar2d_layout, m/s,       grid_name);
    add_field<Computed>("ice_flux",   scalar2d_layout, m/s,       grid_name);
    add_field<Computed>("heat_flux",  scalar2d_layout, W/m2,      grid_name);
  }
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
  // <scheme>check_energy_scaling</scheme>
  // <scheme>check_energy_chng</scheme>

  add_invariant_check<FieldWithinIntervalCheck>(get_field_out("T_mid"),m_grid,100.0,500.0,false);
  add_invariant_check<FieldWithinIntervalCheck>(get_field_out("qv"),m_grid,1e-13,0.2,true);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("qc"),m_grid,0.0,0.1,true);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("qr"),m_grid,0.0,0.1,true);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("qi"),m_grid,0.0,0.1,true);


  Real P0     = PC::P0;     // Reference pressure; pref_in
  Real latvap = PC::LatVap; // Latent heat of vaporization; lv_in
  Real rhoqr  = PC::RHOW;   // rhoqr_in
  Real gravit = PC::gravit; // gravitational acceleration

  kessler::kessler_eamxx_bridge_init(m_num_cols, m_num_levs, latvap, P0, rhoqr, gravit);

  // Allocate host mirror views for GPU -> CPU Fortran bridge
  params_helpers.h_cpair      = KMF::view_2dh<Real>("kessler.h_cpair",     m_num_cols, m_num_levs);
  params_helpers.h_rair       = KMF::view_2dh<Real>("kessler.h_rair",      m_num_cols, m_num_levs);
  params_helpers.h_rho        = KMF::view_2dh<Real>("kessler.h_rho",       m_num_cols, m_num_levs);
  params_helpers.h_pk         = KMF::view_2dh<Real>("kessler.h_pk",        m_num_cols, m_num_levs);
  params_helpers.h_z_mid      = KMF::view_2dh<Real>("kessler.h_z_mid",     m_num_cols, m_num_levs);
  params_helpers.h_phis       = KMF::view_1dh<Real>("kessler.h_phis",      m_num_cols);

  params_computed.h_theta     = KMF::view_2dh<Real>("kessler.h_theta",     m_num_cols, m_num_levs);
  params_computed.h_qv        = KMF::view_2dh<Real>("kessler.h_qv",        m_num_cols, m_num_levs);
  params_computed.h_qc        = KMF::view_2dh<Real>("kessler.h_qc",        m_num_cols, m_num_levs);
  params_computed.h_qr        = KMF::view_2dh<Real>("kessler.h_qr",        m_num_cols, m_num_levs);
  params_computed.h_precl     = KMF::view_1dh<Real>("kessler.h_precl",     m_num_cols);
  params_computed.h_relhum    = KMF::view_2dh<Real>("kessler.h_relhum",    m_num_cols, m_num_levs);
  params_computed.h_temp_prev = KMF::view_2dh<Real>("kessler.h_temp_prev", m_num_cols, m_num_levs);
  params_computed.h_temp      = KMF::view_2dh<Real>("kessler.h_temp",      m_num_cols, m_num_levs);
  params_computed.h_temp_tend = KMF::view_2dh<Real>("kessler.h_temp_tend", m_num_cols, m_num_levs);
  params_computed.h_st_energy = KMF::view_2dh<Real>("kessler.h_st_energy", m_num_cols, m_num_levs);

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
   
  // Pull in variables 
  auto T_mid              = get_field_out("T_mid").get_view<Spack**>();
  auto p_mid              = get_field_in("p_mid").get_view<const Spack**>();
  auto pseudo_density     = get_field_in("pseudo_density").get_view<const Spack**>();
  auto qv                 = get_field_out("qv").get_view<Spack**>();
  auto qc                 = get_field_out("qc").get_view<Spack**>();
  auto qr                 = get_field_out("qr").get_view<Spack**>();
  auto rho                = get_field_out("rho").get_view<Spack**>();
  auto dz                 = get_field_out("dz").get_view<Spack**>();
  auto pk                 = get_field_out("pk").get_view<Spack**>();
  auto theta              = get_field_out("theta").get_view<Spack**>();
  auto precl              = get_field_out("precl").get_view<Real*>();
  auto relhum             = get_field_out("relhum").get_view<Spack**>();
  // Need phis for z_mid computation
  auto phis               = get_field_out("phis").get_view<Real*>();
  // <scheme>kessler_update</scheme> // updates st_energy & temperature related things
  auto T_mid_prev        = get_field_out("T_mid_prev").get_view<Spack**>();
  auto T_mid_tend        = get_field_out("T_mid_tend").get_view<Spack**>();
  auto st_energy         = get_field_out("st_energy").get_view<Spack**>();

  // Get lyr_surf, lyr_toa
  const int lyr_surf = m_num_levs;
  const int lyr_toa  = 1;

  m_atm_logger->info("[EAMxx] kessler run_impl: ");
  
  // NOTE: EAMxx uses a mass-weighted vertical discretization, which is different 
  // than CAM's volume-weighted vertical coordinate. 
  // You will need to derive rho and z from the pseudo_density used in EAMxx.  
  // The change from volume to mass-weighted was to reduce the amount of data movement 
  // and improve numerical accuracy which will help with running reduced precision in the physics parameterizations.
 
  // do the conversion and PF::exner_function, PF::calculate_theta_from_T

  using PF           = scream::PhysicsFunctions<DefaultDevice>;
  using PC           = scream::physics::Constants<Real>;
  const Real inv_ggr = 1/(PC::gravit);

  int nlevs = m_num_levs; // local var, to avoid accessing *this.

  Kokkos::parallel_for(
      "Kessler_preprocess", KT::RangePolicy(0, m_num_cols * nlevs),
      KOKKOS_CLASS_LAMBDA(const int i) {
        const int icol = i / nlevs;
        const int klev = i % nlevs;
        pk(icol, klev / Spack::n)[klev % Spack::n]    = PF::exner_function(p_mid(icol, klev / Spack::n)[klev % Spack::n]);
        theta(icol, klev / Spack::n)[klev % Spack::n] = PF::calculate_theta_from_T(T_mid(icol, klev / Spack::n)[klev % Spack::n],p_mid(icol, klev / Spack::n)[klev % Spack::n]);

        // Vertical layer thickness
        dz(icol, klev / Spack::n)[klev % Spack::n]  = PF::calculate_dz(pseudo_density(icol, klev / Spack::n)[klev % Spack::n], p_mid(icol, klev / Spack::n)[klev % Spack::n], T_mid(icol, klev / Spack::n)[klev % Spack::n], qv(icol, klev / Spack::n)[klev % Spack::n]);
        rho(icol, klev / Spack::n)[klev % Spack::n] = inv_ggr*( pseudo_density(icol, klev / Spack::n)[klev % Spack::n]/dz(icol, klev / Spack::n)[klev % Spack::n] );
      }
    );

  // set up params_helpers struct
  params_helpers.rho        = rho;
  params_helpers.dz         = dz;
  params_helpers.pk         = pk;
  params_helpers.phis       = phis;

  // setup params_computed struct
  params_computed.theta     = theta;
  params_computed.qv        = qv;
  params_computed.qc        = qc;
  params_computed.qr        = qr;
  params_computed.precl     = precl;
  params_computed.relhum    = relhum;

  params_computed.st_energy = st_energy;
  params_computed.temp      = T_mid;
  params_computed.temp_prev = T_mid_prev;
  params_computed.temp_tend = T_mid_tend;

  // Compute z_mid and z_int into params_in (needed by kessler_run for heights)
  // calculate_z_int() contains a team-level parallel_scan, which requires a special policy
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;
  const int nlev_packs   = ekat::npack<Spack>(m_num_levs);
  const auto scan_policy = TPF::get_thread_range_parallel_scan_team_policy(m_num_cols, nlev_packs);

  printf("nlevs: %d \n", nlevs);

  Kokkos::parallel_for(scan_policy, KOKKOS_CLASS_LAMBDA (const KT::MemberType& team) {
    const int i = team.league_rank();

    auto z_mid_i = ekat::subview(params_helpers.z_mid, i);
    auto dz_i    = ekat::subview(params_helpers.dz, i);
    auto z_int_i = ekat::subview(params_helpers.z_int, i);
    Real z_surf  = 0.0; 

    PF::calculate_z_int(team, nlevs, dz_i, z_surf, z_int_i);
    team.team_barrier();
    PF::calculate_z_mid(team, nlevs, z_int_i, z_mid_i);
    team.team_barrier();
  });

  // Initialize fortran data holders in structs 
  params_helpers.init(m_num_cols, nlevs);
  params_computed.init(m_num_cols, nlevs);

  double dt_timestep = dt;

  // This calls both kessler_rn and kessler_update now
  kessler_eamxx_bridge_run(m_num_cols, nlevs, dt_timestep, lyr_surf, lyr_toa, params_helpers, params_computed);

  // <scheme>qneg</scheme> // this is taken care of by the postcondition checks we have in place for qc, qr, and qi (these checks will set any negative values to 0.0)
  // <scheme>geopotential_temp</scheme> // -> done above when calculating z_mid

  // <scheme>sima_state_diagnostics</scheme> // -> just writes fields to file (taken care of by EAMxx & set in the output_fields.yml file)
  // <>scheme>kessler_diagnostics</scheme> // -> only writes out precl field to file (taken care of by EAMxx & set in the output_fields.yml file)

  // <scheme>thermo_water_update</scheme> // computes enthalpy using cpair (I think computed by the energy fixer provided by homme)

  // <!-- MPAS and SE specific scaling of temperature for enforcing energy consistency:
  //       First, calculate the scaling based off cp_or_cv_dycore (from cam_thermo_water_update)
  //       Then, perform the temperature and temperature tendency scaling,
  //       and apply tendencies resulting from such adjustment -->
  // <scheme>check_energy_scaling</scheme> 
    // real(kind_phys),    intent(out)    :: flx_vap(:)     ! boundary flux of vapor [kg m-2 s-1]
    // real(kind_phys),    intent(out)    :: flx_cnd(:)     ! boundary flux of liquid+ice (precip?) [m s-1]
    // real(kind_phys),    intent(out)    :: flx_ice(:)     ! boundary flux of ice (snow?) [m s-1]
    // real(kind_phys),    intent(out)    :: flx_sen(:)     ! boundary flux of sensible heat [W m-2]
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

  // pretty sure eamxx does this with the energy fixer
  // <scheme>dycore_energy_consistency_adjust</scheme> // TODO ? -> does energy scaling for temperature; 

  // pretty sure eamxx does this with the energy fixer
  // <scheme>apply_tendency_of_air_temperature (nz, t_tend, temp, dtdT_total, dt, errcode, errmsg) TODO ? -> updates dtdT_total

  // <!-- Tendency diagnostics -->
  // <scheme>sima_tend_diagnostics</scheme> // -> only writes out dTdt_total, dudt_total, dvdt_total to file. Can be set in the output_fields.yml file and taken care of by EAMxx

  // Update with the new values from the run
  Kokkos::parallel_for("update_output_1d",m_num_cols, KOKKOS_CLASS_LAMBDA (const int i) {
    precl(i) = params_computed.precl(i);
    }
  );

  Kokkos::parallel_for("update_output_2d",KT::RangePolicy(0, m_num_cols * nlevs), KOKKOS_CLASS_LAMBDA (const int idx) { // clean up
    const int icol = idx/nlevs;
    const int klev = idx%nlevs;

    rho(icol,klev)        = params_helpers.rho(icol,klev);
    dz(icol,klev)         = params_helpers.dz(icol,klev);
    pk(icol,klev)         = params_helpers.pk(icol,klev);

    qv(icol,klev)         = params_computed.qv(icol,klev);
    qc(icol,klev)         = params_computed.qc(icol,klev);
    qr(icol,klev)         = params_computed.qr(icol,klev);
    theta(icol,klev)      = params_computed.theta(icol,klev);
    relhum(icol,klev)     = params_computed.relhum(icol,klev);

    T_mid_prev(icol,klev) = params_computed.temp_prev(icol,klev);
    T_mid_tend(icol,klev) = params_computed.temp_tend(icol,klev);
    T_mid(icol,klev)      = params_computed.temp(icol,klev);
    st_energy(icol,klev)  = params_computed.st_energy(icol,klev);
    // dont need to do anything for z_mid since we recompute it each time :(
    }
  );

  m_atm_logger->info("[EAMxx] kessler run_impl - end ");
}

// =========================================================================================
//  Inputs:
//      None
//  Outputs:
//      None

void KesslerMicrophysics::finalize_impl()
{
  // Do nothing
  m_atm_logger->info("[EAMxx] Kessler processes clean up.");
}
// =========================================================================================
//  Inputs:
//      None
//  Outputs:
//      buffer_size - size in bytes needed for ATMBufferManager buffers

size_t KesslerMicrophysics::requested_buffer_size_in_bytes() const
{
  const int nlevm_packs    = ekat::npack<Spack>(m_num_levs);
  const int nlev_int_packs = ekat::npack<Spack>(m_num_levs+1);

  constexpr auto num_1d_intgr   = KMF::params_helpers::num_1d_intgr   + KMF::params_computed::num_1d_intgr;
  constexpr auto num_1d_scalr   = KMF::params_helpers::num_1d_scalr   + KMF::params_computed::num_1d_scalr;
  constexpr auto num_2d_midlv_c = KMF::params_helpers::num_2d_c       + KMF::params_computed::num_2d_c;
  constexpr auto num_2d_intlv_c = KMF::params_helpers::num_2d_intlv_c + KMF::params_computed::num_2d_intlv_c;
  constexpr auto num_2d_midlv_f = KMF::params_helpers::num_2d_f       + KMF::params_computed::num_2d_f;

  size_t buffer_size = 0;

  buffer_size+= num_1d_intgr   * sizeof(Int)    * m_num_cols;                  // should be 0
  buffer_size+= num_1d_scalr   * sizeof(Scalar) * m_num_cols;                  // should be 2, C++ holders
  buffer_size+= num_1d_scalr   * sizeof(Real)   * m_num_cols;                  // should be 2, for fortran holders
  buffer_size+= num_2d_midlv_c * sizeof(Spack)  * m_num_cols * nlevm_packs;    // should be 14, C++ holders
  buffer_size+= num_2d_intlv_c * sizeof(Spack)  * m_num_cols * nlev_int_packs; // should be 1, C++ holders
  buffer_size+= num_2d_midlv_f * sizeof(Real)   * m_num_cols * m_num_levs;     // should be 14, for fortran holders

  return buffer_size;
}

/*------------------------------------------------------------------------------------------------*/
//  Inputs:
//      buffer_manager - pointer to the ATMBufferManager buffer. This manages memory buffers for the whole system
//  Outputs:
//      None
void KesslerMicrophysics::init_buffers(const ATMBufferManager &buffer_manager)
{
  // TODO - after adding all variables to structs, update this function
  auto buffer_chk = ( buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes() );
  EKAT_REQUIRE_MSG(buffer_chk,"Error! Buffers size not sufficient.\n");
  
  const int nlev_mid_packs = ekat::npack<Spack>(m_num_levs);
  const int nlev_int_packs = ekat::npack<Spack>(m_num_levs+1);

  constexpr auto num_1d_intgr   = KMF::params_helpers::num_1d_intgr   + KMF::params_computed::num_1d_intgr;
  constexpr auto num_1d_scalr   = KMF::params_helpers::num_1d_scalr   + KMF::params_computed::num_1d_scalr;
  constexpr auto num_2d_midlv_c = KMF::params_helpers::num_2d_c       + KMF::params_computed::num_2d_c;
  constexpr auto num_2d_intlv_c = KMF::params_helpers::num_2d_intlv_c + KMF::params_computed::num_2d_intlv_c;
  constexpr auto num_2d_midlv_f = KMF::params_helpers::num_2d_f       + KMF::params_computed::num_2d_f;

  Scalar* scl_mem = reinterpret_cast<Scalar*>(buffer_manager.get_memory());
  //----------------------------------------------------------------------------
  // device 1D integer variables
  KMF::view_1d<Scalar>* ptrs_1d_scalr[num_1d_scalr] = { &params_computed.precl, 
                                                        &params_helpers.phis 
                                                      };
  for (auto& v : ptrs_1d_scalr) {
    *v = KMF::view_1d<Scalar>(scl_mem, m_num_cols);
    scl_mem += v->size();
  }
  //----------------------------------------------------------------------------
  Real* r1_mem = reinterpret_cast<Real*>(scl_mem);
  //----------------------------------------------------------------------------
  // device 1D scalar scalars
  KMF::uview_1d<Real>* ptrs_1d_real[num_1d_scalr] = { &params_computed.f_precl, 
                                                      &params_helpers.f_phis
                                                    };
  for (auto& v : ptrs_1d_real) {
    *v = KMF::uview_1d<Real>(r1_mem, m_num_cols);
    r1_mem += v->size();
  } 
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  Real* r_mem = reinterpret_cast<Real*>(r1_mem);
  //----------------------------------------------------------------------------
  // 2D "f_" views
  KMF::uview_2dl<Real>* midlv_f_ptrs[num_2d_midlv_f] = { &params_helpers.f_cpair,
                                                         &params_helpers.f_rair,
                                                         &params_helpers.f_rho,
                                                         &params_helpers.f_pk,
                                                         &params_helpers.f_z_mid,
                                                         &params_computed.f_theta,
                                                         &params_computed.f_qv,
                                                         &params_computed.f_qc,
                                                         &params_computed.f_qr,
                                                         &params_computed.f_relhum,
                                                         &params_computed.f_temp_prev,
                                                         &params_computed.f_temp,
                                                         &params_computed.f_temp_tend,
                                                         &params_computed.f_st_energy
                                                        };
  for (int i=0; i<num_2d_midlv_f; ++i) {
    *midlv_f_ptrs[i] = KMF::uview_2dl<Real>(r_mem, m_num_cols, m_num_levs);
    r_mem += midlv_f_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  Spack* spk_mem = reinterpret_cast<Spack*>(r_mem);
  //----------------------------------------------------------------------------
  // 2D views 
  KMF::view_2d<Spack>* midlv_c_ptrs[num_2d_midlv_c] = { &params_helpers.rho,
                                                        &params_helpers.dz,
                                                        &params_helpers.pk,
                                                        &params_helpers.z_mid,
                                                        &params_computed.theta,
                                                        &params_computed.qv,
                                                        &params_computed.qc,
                                                        &params_computed.qr,
                                                        &params_computed.relhum,
                                                        &params_computed.temp_prev,
                                                        &params_computed.temp,
                                                        &params_computed.temp_tend,
                                                        &params_computed.st_energy
                                                      };
  for (int i=0; i<num_2d_midlv_c; ++i) {
    *midlv_c_ptrs[i] = KMF::view_2d<Spack>(spk_mem, m_num_cols, nlev_mid_packs);
    spk_mem += midlv_c_ptrs[i]->size();
  }

  // allocate z_int which is our only interface variable. 
  KMF::view_2d<Spack>* intlv_c_ptrs[num_2d_intlv_c] = { &params_helpers.z_int };
  for (int i=0; i<num_2d_intlv_c; ++i) {
     *intlv_c_ptrs[i] = KMF::view_2d<Spack>(spk_mem, m_num_cols, nlev_int_packs);
     spk_mem += intlv_c_ptrs[i]->size();
  }

  //----------------------------------------------------------------------------
  Real* total_mem = reinterpret_cast<Real*>(spk_mem);
  size_t used_mem = (reinterpret_cast<Real*>(total_mem) - buffer_manager.get_memory())*sizeof(Real);
  auto req_mem = requested_buffer_size_in_bytes();
  auto mem_chk = ( used_mem == req_mem );
  EKAT_REQUIRE_MSG(mem_chk,"Error! Used memory ("+ std::to_string(used_mem) + ") != requested memory ("+ std::to_string(req_mem) + ") for Kessler.");
}
} // namespace scream
