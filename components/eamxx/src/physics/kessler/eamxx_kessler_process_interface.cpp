#include "eamxx_kessler_process_interface.hpp"
// #include "kessler_eamxx_bridge.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/field/field_utils.hpp"

#include "share/physics/physics_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
// #include "share/physics/eamxx_common_physics_functions_impls.hpp"

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
  using namespace ekat::units;
  // using namespace ShortFieldTagsNames;

  m_atm_logger->info("[EAMxx] Kessler processes set grids");

  static constexpr auto m2 = (m*m).rename("m2");
  static constexpr auto m3 = (m*m*m).rename("m3");
  static constexpr auto s2 = (s*s).rename("s2");
  constexpr auto nondim    = ekat::units::Units::nondimensional();
  constexpr int pack_size  = Pack::n;


  // specify which grid to use
  m_grid = m_grids_manager->get_grid("physics");

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

  add_field<Computed>("cpair", scalar3d_layout_mid, J/kg/K, grid_name, pack_size);
  add_field<Computed>("rair",  scalar3d_layout_mid, J/kg/K, grid_name, pack_size);
  add_field<Computed>("z_mid",  scalar3d_layout_mid, m, grid_name, pack_size);
  add_field<Computed>("z_int",  scalar3d_layout_mid, m, grid_name, pack_size);

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

  const Real P0     = PC::P0.value;     // Reference pressure; pref_in
  const Real latvap = PC::LatVap.value; // Latent heat of vaporization; lv_in
  const Real rhoqr  = PC::RHOW.value;   // rhoqr_in
  const Real gravit = PC::gravit.value; // gravitational acceleration

  // m_params.set<std::string>("py_module_name","kessler_jax");
  const auto& py_module_name = m_params.get<std::string>("py_module_name");
  const auto& py_module_path = m_params.get<std::string>("py_module_path","./");
  
  m_atm_logger->info("[EAMxx] kessler py_module_name: "+ py_module_name);
  m_atm_logger->info("[EAMxx] kessler py_module_path: "+ py_module_path);


  // The JAX code currently keeps Cpair and Rair as 2D arrays so we need to set them here. 
  #ifdef EAMXX_HAS_PYTHON
    const Real Cpair  = PC::Cpair.value; // Specific heat of dry air at constant pressure
    const Real Rair   = PC::Rair.value;  // Gas constant of dry air
    int nlevs = m_num_levs; // local var, to avoid accessing *this.

    if (has_py_module()) {
      py_module_call("init", latvap, P0, rhoqr, gravit);
      m_atm_logger->info("[EAMxx] kessler called python init");
    }

    auto cpair = get_field_out("cpair").get_view<Real**>();
    auto rair  = get_field_out("rair").get_view<Real**>();

    Kokkos::parallel_for("py_air_const",KT::RangePolicy(0, m_num_cols * nlevs), KOKKOS_CLASS_LAMBDA (const int idx) { 
    const int icol = idx/nlevs;
    const int klev = idx%nlevs;

    cpair(icol,klev)        = Cpair;
    rair(icol,klev)         = Rair;
    }
  );
  #endif
  // kessler::kessler_eamxx_bridge_init(m_num_cols, m_num_levs, latvap, P0, rhoqr, gravit);
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
  auto T_mid_v              = get_field_out("T_mid").get_view<Pack**>();
  auto p_mid_v              = get_field_in("p_mid").get_view<const Pack**>();
  auto pseudo_density_v     = get_field_in("pseudo_density").get_view<const Pack**>();
  auto qv_v                 = get_field_out("qv").get_view<Pack**>();
  auto qc_v                 = get_field_out("qc").get_view<Pack**>();
  auto qr_v                 = get_field_out("qr").get_view<Pack**>();
  auto rho_v                = get_field_out("rho").get_view<Pack**>();
  auto dz_v                 = get_field_out("dz").get_view<Pack**>();
  auto pk_v                 = get_field_out("pk").get_view<Pack**>();
  auto theta_v              = get_field_out("theta").get_view<Pack**>();
  auto precl_v              = get_field_out("precl").get_view<Real*>();
  auto relhum_v             = get_field_out("relhum").get_view<Pack**>();
  // Need phis for z_mid computation
  auto phis_v               = get_field_out("phis").get_view<Real*>();
  // <scheme>kessler_update</scheme> // updates st_energy & temperature related things
  auto T_mid_prev_v        = get_field_out("T_mid_prev").get_view<Pack**>();
  auto T_mid_tend_v        = get_field_out("T_mid_tend").get_view<Pack**>();
  auto st_energy_v         = get_field_out("st_energy").get_view<Pack**>();

  // Geopotential height and interface height for kessler_run
  auto z_mid_v              = get_field_out("z_mid").get_view<Pack**>();
  auto z_int_v              = get_field_out("z_int").get_view<Pack**>();

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
  const Real inv_ggr = 1/(PC::gravit.value);

  int nlevs = m_num_levs; // local var, to avoid accessing *this.

  Kokkos::parallel_for(
      "Kessler_preprocess", KT::RangePolicy(0, m_num_cols * nlevs),
      KOKKOS_CLASS_LAMBDA(const int i) {
        const int icol = i / nlevs;
        const int klev = i % nlevs;
        pk_v(icol, klev / Pack::n)[klev % Pack::n]    = PF::exner_function(p_mid_v(icol, klev / Pack::n)[klev % Pack::n]);
        theta_v(icol, klev / Pack::n)[klev % Pack::n] = PF::calculate_theta_from_T(T_mid_v(icol, klev / Pack::n)[klev % Pack::n],p_mid_v(icol, klev / Pack::n)[klev % Pack::n]);

        // Vertical layer thickness
        dz_v(icol, klev / Pack::n)[klev % Pack::n]  = PF::calculate_dz(pseudo_density_v(icol, klev / Pack::n)[klev % Pack::n], p_mid_v(icol, klev / Pack::n)[klev % Pack::n], T_mid_v(icol, klev / Pack::n)[klev % Pack::n], qv_v(icol, klev / Pack::n)[klev % Pack::n]);
        rho_v(icol, klev / Pack::n)[klev % Pack::n] = inv_ggr*( pseudo_density_v(icol, klev / Pack::n)[klev % Pack::n]/dz_v(icol, klev / Pack::n)[klev % Pack::n] );
      }
    );

  // // set up params_helpers struct
  // params_helpers.rho        = rho;
  // params_helpers.dz         = dz;
  // params_helpers.pk         = pk;
  // params_helpers.phis       = phis;

  // // setup params_computed struct
  // params_computed.theta     = theta;
  // params_computed.qv        = qv;
  // params_computed.qc        = qc;
  // params_computed.qr        = qr;
  // params_computed.precl     = precl;
  // params_computed.relhum    = relhum;

  // params_computed.st_energy = st_energy;
  // params_computed.temp      = T_mid;
  // params_computed.temp_prev = T_mid_prev;
  // params_computed.temp_tend = T_mid_tend;

  // Compute z_mid and z_int into params_in (needed by kessler_run for heights)
  // calculate_z_int() contains a team-level parallel_scan, which requires a special policy
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;
  const int nlev_packs   = ekat::npack<Pack>(m_num_levs);
  const auto scan_policy = TPF::get_thread_range_parallel_scan_team_policy(m_num_cols, nlev_packs);

  Kokkos::parallel_for(scan_policy, KOKKOS_CLASS_LAMBDA (const KT::MemberType& team) {
    const int i = team.league_rank();

    auto z_mid_i = ekat::subview(z_mid_v, i);
    auto dz_i    = ekat::subview(dz_v, i);
    auto z_int_i = ekat::subview(z_int_v, i);
    Real z_surf  = 0.0; 

    PF::calculate_z_int(team, nlevs, dz_i, z_surf, z_int_i);
    team.team_barrier();
    PF::calculate_z_mid(team, nlevs, z_int_i, z_mid_i);
    team.team_barrier();
  });

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
  m_atm_logger->info("[EAMxx] kessler run_impl: done with pre-processing");
  // For python bindings we can't use the field views. Kokkos loops require using views so they are left above.
  auto qv     = get_field_out("qv");
  auto qc     = get_field_out("qc");
  auto qr     = get_field_out("qr");
  auto cpair  = get_field_out("cpair");
  auto rair   = get_field_out("rair");
  auto rho    = get_field_out("rho");
  auto z_mid  = get_field_out("z_mid");
  auto pk     = get_field_out("pk");
  auto theta  = get_field_out("theta");
  auto precl  = get_field_out("precl");
  auto relhum = get_field_out("relhum");
  // Needed for the update functions
  auto phis       = get_field_out("phis");
  auto T_mid      = get_field_out("T_mid");
  auto T_mid_prev = get_field_out("T_mid_prev");
  auto T_mid_tend = get_field_out("T_mid_tend");
  auto st_energy  = get_field_out("st_energy");

  double dt_timestep = dt;
  #ifdef EAMXX_HAS_PYTHON
  m_atm_logger->info("[EAMxx] kessler run_impl: in python check");
    if (has_py_module()) {
      pybind11::array py_qv, py_qc, py_qr,
                      py_cpair, py_rair, py_rho, py_z_mid, py_pk,
                      py_theta, py_precl, py_relhum,
                      py_phis, py_temp, py_temp_prev, py_temp_tend, py_st_energy;

      if (m_params.get<std::string>("py_backend")=="device") {
        py_qv                = get_py_field_dev("qv");
        py_qc                = get_py_field_dev("qc");
        py_qr                = get_py_field_dev("qr");
        py_cpair             = get_py_field_dev("cpair");
        py_rair              = get_py_field_dev("rair");
        py_rho               = get_py_field_dev("rho");
        py_z_mid             = get_py_field_dev("z_mid");
        py_pk                = get_py_field_dev("pk");
        py_theta             = get_py_field_dev("theta");
        py_precl             = get_py_field_dev("precl");
        py_relhum            = get_py_field_dev("relhum");
        py_phis              = get_py_field_dev("phis");
        py_temp              = get_py_field_dev("T_mid");
        py_temp_prev         = get_py_field_dev("T_mid_prev");
        py_temp_tend         = get_py_field_dev("T_mid_tend");
        py_st_energy         = get_py_field_dev("st_energy");

      } else {
        qv.sync_to_host();
        qc.sync_to_host();
        qr.sync_to_host();
        cpair.sync_to_host();
        rair.sync_to_host();
        rho.sync_to_host();
        z_mid.sync_to_host();
        pk.sync_to_host();
        theta.sync_to_host();
        precl.sync_to_host();
        relhum.sync_to_host();
        phis.sync_to_host();
        T_mid.sync_to_host();
        T_mid_prev.sync_to_host();
        T_mid_tend.sync_to_host();
        st_energy.sync_to_host();

        py_qv                = get_py_field_host("qv");
        py_qc                = get_py_field_host("qc");
        py_qr                = get_py_field_host("qr");
        py_cpair             = get_py_field_host("cpair");
        py_rair              = get_py_field_host("rair");
        py_rho               = get_py_field_host("rho");
        py_z_mid             = get_py_field_host("z_mid");
        py_pk                = get_py_field_host("pk");
        py_theta             = get_py_field_host("theta");
        py_precl             = get_py_field_host("precl");
        py_relhum            = get_py_field_host("relhum");
        py_phis              = get_py_field_host("phis");
        py_temp              = get_py_field_host("T_mid");
        py_temp_prev         = get_py_field_host("T_mid_prev");
        py_temp_tend         = get_py_field_host("T_mid_tend");
        py_st_energy         = get_py_field_host("st_energy");
      }

      // NOTE: kessler_run's "z" argument expects heights (z_mid), not layer
      // thickness (dz) 
      py_module_call("run", m_num_cols, nlevs, dt_timestep, lyr_surf, lyr_toa,
                    py_cpair,
                    py_rair,
                    py_rho,
                    py_z_mid,
                    py_pk,
                    py_theta,
                    py_qv,
                    py_qc,
                    py_qr,
                    py_precl,
                    py_relhum);
      m_atm_logger->info("[EAMxx] kessler called python run");

      py_module_call("update", m_num_cols, nlevs, dt_timestep,
                    py_cpair,
                    py_z_mid,
                    py_pk,
                    py_theta,
                    py_phis,
                    py_temp,
                    py_temp_prev,
                    py_temp_tend,
                    py_st_energy);
      m_atm_logger->info("[EAMxx] kessler called python update");

      if (m_params.get<std::string>("py_backend")=="host") {
        qv.sync_to_dev();
        qr.sync_to_dev();
        qc.sync_to_dev();
        cpair.sync_to_dev();
        rair.sync_to_dev();
        pk.sync_to_dev();
        theta.sync_to_dev();
        precl.sync_to_dev();
        relhum.sync_to_dev();
        phis.sync_to_dev();
        T_mid.sync_to_dev();
        T_mid_prev.sync_to_dev();
        T_mid_tend.sync_to_dev();
        st_energy.sync_to_dev();
      }

      m_atm_logger->info("[EAMxx] kessler run_impl - end ");
      return;
    }
  #endif

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

// size_t KesslerMicrophysics::requested_buffer_size_in_bytes() const
// {
//   const int nlevm_packs    = ekat::npack<Pack>(m_num_levs);
//   const int nlev_int_packs = ekat::npack<Pack>(m_num_levs+1);

//   constexpr auto num_1d_intgr   = KMF::params_helpers::num_1d_intgr   + KMF::params_computed::num_1d_intgr;
//   constexpr auto num_1d_scalr   = KMF::params_helpers::num_1d_scalr   + KMF::params_computed::num_1d_scalr;
//   constexpr auto num_2d_midlv_c = KMF::params_helpers::num_2d_c       + KMF::params_computed::num_2d_c;
//   constexpr auto num_2d_intlv_c = KMF::params_helpers::num_2d_intlv_c + KMF::params_computed::num_2d_intlv_c;
//   constexpr auto num_2d_midlv_f = KMF::params_helpers::num_2d_f       + KMF::params_computed::num_2d_f;

//   size_t buffer_size = 0;

//   buffer_size+= num_1d_intgr   * sizeof(Int)    * m_num_cols;                  // should be 0
//   buffer_size+= num_1d_scalr   * sizeof(Scalar) * m_num_cols;                  // should be 2, C++ holders
//   buffer_size+= num_1d_scalr   * sizeof(Real)   * m_num_cols;                  // should be 2, for fortran holders
//   buffer_size+= num_2d_midlv_c * sizeof(Pack)  * m_num_cols * nlevm_packs;    // should be 13, C++ holders
//   buffer_size+= num_2d_intlv_c * sizeof(Pack)  * m_num_cols * nlev_int_packs; // should be 1, C++ holders
//   buffer_size+= num_2d_midlv_f * sizeof(Real)   * m_num_cols * m_num_levs;     // should be 14, for fortran holders

//   return buffer_size;
// }

// /*------------------------------------------------------------------------------------------------*/
// //  Inputs:
// //      buffer_manager - pointer to the ATMBufferManager buffer. This manages memory buffers for the whole system
// //  Outputs:
// //      None
// void KesslerMicrophysics::init_buffers(const ATMBufferManager &buffer_manager)
// {
//   // TODO - after adding all variables to structs, update this function
//   auto buffer_chk = ( buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes() );
//   EKAT_REQUIRE_MSG(buffer_chk,"Error! Buffers size not sufficient.\n");
  
//   const int nlev_mid_packs = ekat::npack<Pack>(m_num_levs);
//   const int nlev_int_packs = ekat::npack<Pack>(m_num_levs+1);

//   constexpr auto num_1d_intgr   = KMF::params_helpers::num_1d_intgr   + KMF::params_computed::num_1d_intgr;
//   constexpr auto num_1d_scalr   = KMF::params_helpers::num_1d_scalr   + KMF::params_computed::num_1d_scalr;
//   constexpr auto num_2d_midlv_c = KMF::params_helpers::num_2d_c       + KMF::params_computed::num_2d_c;
//   constexpr auto num_2d_intlv_c = KMF::params_helpers::num_2d_intlv_c + KMF::params_computed::num_2d_intlv_c;
//   constexpr auto num_2d_midlv_f = KMF::params_helpers::num_2d_f       + KMF::params_computed::num_2d_f;

//   Scalar* scl_mem = reinterpret_cast<Scalar*>(buffer_manager.get_memory());
//   //----------------------------------------------------------------------------
//   // device 1D integer variables
//   KMF::view_1d<Scalar>* ptrs_1d_scalr[num_1d_scalr] = { &params_computed.precl, 
//                                                         &params_helpers.phis 
//                                                       };
//   for (auto& v : ptrs_1d_scalr) {
//     *v = KMF::view_1d<Scalar>(scl_mem, m_num_cols);
//     scl_mem += v->size();
//   }
//   //----------------------------------------------------------------------------
//   Real* r1_mem = reinterpret_cast<Real*>(scl_mem);
//   //----------------------------------------------------------------------------
//   // device 1D scalar scalars
//   KMF::fview_1d<Real>* ptrs_1d_real[num_1d_scalr] = { &params_computed.f_precl, 
//                                                       &params_helpers.f_phis
//                                                     };
//   for (auto& v : ptrs_1d_real) {
//     *v = KMF::fview_1d<Real>(r1_mem, m_num_cols);
//     r1_mem += v->size();
//   } 
//   //----------------------------------------------------------------------------
//   //----------------------------------------------------------------------------
//   Real* r_mem = reinterpret_cast<Real*>(r1_mem);
//   //----------------------------------------------------------------------------
//   // 2D "f_" views
//   KMF::fview_2dl<Real>* midlv_f_ptrs[num_2d_midlv_f] = { &params_helpers.f_cpair,
//                                                          &params_helpers.f_rair,
//                                                          &params_helpers.f_rho,
//                                                          &params_helpers.f_pk,
//                                                          &params_helpers.f_z_mid,
//                                                          &params_computed.f_theta,
//                                                          &params_computed.f_qv,
//                                                          &params_computed.f_qc,
//                                                          &params_computed.f_qr,
//                                                          &params_computed.f_relhum,
//                                                          &params_computed.f_temp_prev,
//                                                          &params_computed.f_temp,
//                                                          &params_computed.f_temp_tend,
//                                                          &params_computed.f_st_energy
//                                                         };
//   for (int i=0; i<num_2d_midlv_f; ++i) {
//     *midlv_f_ptrs[i] = KMF::fview_2dl<Real>(r_mem, m_num_cols, m_num_levs);
//     r_mem += midlv_f_ptrs[i]->size();
//   }
//   //----------------------------------------------------------------------------
//   Pack* spk_mem = reinterpret_cast<Pack*>(r_mem);
//   //----------------------------------------------------------------------------
//   // 2D views 
//   KMF::view_2d<Pack>* midlv_c_ptrs[num_2d_midlv_c] = { &params_helpers.rho,
//                                                         &params_helpers.dz,
//                                                         &params_helpers.pk,
//                                                         &params_helpers.z_mid,
//                                                         &params_computed.theta,
//                                                         &params_computed.qv,
//                                                         &params_computed.qc,
//                                                         &params_computed.qr,
//                                                         &params_computed.relhum,
//                                                         &params_computed.temp_prev,
//                                                         &params_computed.temp,
//                                                         &params_computed.temp_tend,
//                                                         &params_computed.st_energy
//                                                       };
//   for (int i=0; i<num_2d_midlv_c; ++i) {
//     *midlv_c_ptrs[i] = KMF::view_2d<Pack>(spk_mem, m_num_cols, nlev_mid_packs);
//     spk_mem += midlv_c_ptrs[i]->size();
//   }

//   // allocate z_int which is our only interface variable. 
//   KMF::view_2d<Pack>* intlv_c_ptrs[num_2d_intlv_c] = { &params_helpers.z_int };
//   for (int i=0; i<num_2d_intlv_c; ++i) {
//      *intlv_c_ptrs[i] = KMF::view_2d<Pack>(spk_mem, m_num_cols, nlev_int_packs);
//      spk_mem += intlv_c_ptrs[i]->size();
//   }

//   //----------------------------------------------------------------------------
//   Real* total_mem = reinterpret_cast<Real*>(spk_mem);
//   size_t used_mem = (reinterpret_cast<Real*>(total_mem) - buffer_manager.get_memory())*sizeof(Real);
//   auto req_mem = requested_buffer_size_in_bytes();
//   auto mem_chk = ( used_mem == req_mem );
//   EKAT_REQUIRE_MSG(mem_chk,"Error! Used memory ("+ std::to_string(used_mem) + ") != requested memory ("+ std::to_string(req_mem) + ") for Kessler.");
// }
} // namespace scream
