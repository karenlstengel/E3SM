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

  // Set the log filename in the F90 interface
  // const char* logname = m_atm_logger->get_logfile_name().c_str();
  // set_log_file_name_f90(&logname);

}

// =========================================================================================
//  Inputs:
//      grids_manager - The grids manager used by the dynamcial core
//  Outputs:
//      None

void KesslerMicrophysics::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  // using namespace ekat::units;
  // using namespace ShortFieldTagsNames;

  m_atm_logger->info("[EAMxx] Kessler processes set grids");

  constexpr auto K = ekat::units::K;
  constexpr auto Pa = ekat::units::Pa;
  constexpr auto s = ekat::units::s;
  constexpr auto m = ekat::units::m;
  const auto m3    = pow(m,3);
  constexpr auto kg = ekat::units::kg;
  constexpr auto nondim = ekat::units::Units::nondimensional();
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

  add_field<Required>("T_mid",          scalar3d_layout_mid, K,     grid_name, pack_size);
  add_field<Required>("p_mid",          scalar3d_layout_mid, Pa,    grid_name, pack_size);
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa,    grid_name, pack_size);
  add_tracer<Updated>("qv",             m_grid,              kg/kg,            pack_size);
  add_tracer<Updated>("qc",             m_grid,              kg/kg,            pack_size);
  add_tracer<Updated>("qr",             m_grid,              kg/kg,            pack_size);
  // Locally needed 
  // real(kind_phys),  intent(in)    :: rho(:,:)   ! Dry air density (kg/m^3)
  // real(kind_phys),  intent(in)    :: z(:,:)     ! Heights of thermo. levels (m)
  // real(kind_phys),  intent(in)    :: pk(:,:)    ! Exner function (p/p0)**(R/cp)

  // real(kind_phys),  intent(inout) :: theta(:,:) ! Potential temperature (K)

  // real(kind_phys),  intent(out)   :: precl(:)   ! Precipitation rate (m_water / s)
  // real(kind_phys),  intent(out)   :: relhum(:,:)! Relative humidity in percent

  add_field<Computed>("rho",    scalar3d_layout_mid, kg/m3,  grid_name, pack_size);
  add_field<Computed>("z",      scalar3d_layout_mid, m,      grid_name, pack_size);
  add_field<Computed>("pk",     scalar3d_layout_mid, nondim, grid_name, pack_size);
  add_field<Computed>("theta",  scalar3d_layout_mid, K,      grid_name, pack_size);
  add_field<Computed>("precl",  scalar2d_layout,     m/s,    grid_name, pack_size);
  add_field<Computed>("relhum", scalar3d_layout_mid, nondim, grid_name, pack_size);

}

// =========================================================================================
//  Inputs:
//      run_type - run type, either initial or restart
//  Outputs:
//      None

void KesslerMicrophysics::initialize_impl (const RunType /* run_type */)
{
  m_atm_logger->info("[EAMxx] kessler processes initialize_impl: ");

  // Set post condition checks
  // <scheme>check_energy_zero_fluxes</scheme>
  // <scheme>check_energy_scaling</scheme>
  // <scheme>check_energy_chng</scheme>

  // const Real P0 = PC::P0; // Reference pressure; pref_in
  // const Real latvap = PC::LatVap; // Latent heat of vaporization; lv_in
  // const Real rhoqr = PC::RHOW; // rhoqr_in

  // kessler::kessler_eamxx_bridge_init(latvap, P0, rhoqr);
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
  // auto T_mid   = get_field_in("T_mid").get_view<Spack**>();
  // auto p_mid   = get_field_in("p_mid").get_view<Spack**>();
  // auto pseudo_density = get_field_in("pseudo_density").get_view<Spack**>();
  // auto qv      = get_field_in("qv").get_view<Spack**>();
  // auto qc      = get_field_in("qc").get_view<Spack**>();
  // auto qr      = get_field_in("qr").get_view<Spack**>();

  // auto rho     = get_field_out("rho").get_view<Spack**>();
  // auto z       = get_field_out("z").get_view<Spack**>();
  // auto pk      = get_field_out("pk").get_view<Spack**>();
  // auto theta   = get_field_out("theta").get_view<Spack**>();
  // auto precl   = get_field_out("precl").get_view<Real*>();
  // auto relhum  = get_field_out("relhum").get_view<Spack**>();

  // // Get lyr_surf, lyr_toa
  // const int lyr_surf = 0;
  // const int lyr_toa = m_num_levs - 1;

  m_atm_logger->info("[EAMxx] kessler run_impl: ");
  
  // NOTE: EAMxx uses a mass-weighted vertical discretization, which is different 
  // than CAM's volume-weighted vertical coordinate. 
  // You will need to derive rho and z from the pseudo_density used in EAMxx.  
  // The change from volume to mass-weighted was to reduce the amount of data movement 
  // and improve numerical accuracy which will help with running reduced precision in the physics parameterizations.
 
  // do the conversion and PF::exner_function, PF::calculate_theta_from_T
  // KMF::preprocess(T_mid, p_mid, qv, pseudo_density, z, rho, pk);

  // // set up params_in struct
  // // params_in.cpair = cpair;
  // // params_in.rair = rair;
  // params_in.rho = rho;
  // params_in.z = z;
  // params_in.pk = pk;

  // // setup params_out struct
  // params_out.theta = theta;
  // params_out.qv = qv;
  // params_out.qc = qc;
  // params_out.qr = qr;
  // params_out.precl = precl;
  // params_out.relhum = relhum;

  // // Initialize fortran data holders in struct
  // params_in.init(m_num_cols, m_num_levs); 
  // params_out.init(m_num_cols, m_num_levs); 

  // kessler_eamxx_bridge_run(m_num_cols, m_num_levs, dt, lyr_surf, lyr_toa, params_in, params_out); 

  // <scheme>potential_temp_to_temp</scheme>
  // <scheme>dry_to_wet_water_vapor</scheme>
  // <scheme>dry_to_wet_cloud_liquid_water</scheme>
  // <scheme>dry_to_wet_rain</scheme>
  // <scheme>kessler_update</scheme>
  // <scheme>qneg</scheme>
  // <scheme>geopotential_temp</scheme>

  // <scheme>sima_state_diagnostics</scheme>
  // <scheme>kessler_diagnostics</scheme>

  // <scheme>thermo_water_update</scheme>

  // <!-- COMMENTED OUT until qini/liqini/iceini have initialization routines -->
  // <!-- See https://github.com/ESCOMP/atmospheric_physics/issues/222 -->
  // <!-- <scheme>dme_adjust</scheme> -->

  // <!-- MPAS and SE specific scaling of temperature for enforcing energy consistency:
  //       First, calculate the scaling based off cp_or_cv_dycore (from cam_thermo_water_update)
  //       Then, perform the temperature and temperature tendency scaling,
  //       and apply tendencies resulting from such adjustment -->
  // <scheme>check_energy_scaling</scheme>
  // <scheme>dycore_energy_consistency_adjust</scheme>
  // <scheme>apply_tendency_of_air_temperature</scheme>

  // <!-- Tendency diagnostics -->
  // <scheme>sima_tend_diagnostics</scheme>

  // Update with the new values from the run
  // TODO - 
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
  const int nlevm_packs = ekat::npack<Spack>(m_num_levs);
  // const int nlev_int_packs = ekat::npack<Spack>(m_num_levs+1);

  size_t buffer_size = 0;

  buffer_size+= KMF::params_in::num_1d_intgr * sizeof(Int)   * m_num_cols;
  buffer_size+= KMF::params_in::num_1d_scalr * sizeof(Scalar)* m_num_cols;
  buffer_size+= KMF::params_in::num_2d_c * sizeof(Spack) * m_num_cols * nlevm_packs;
  // buffer_size+= KMF::params_in::num_2d_intfc * sizeof(Spack) * m_num_cols * nlev_int_packs;

  buffer_size+= KMF::params_out::num_1d_intgr * sizeof(Int)   * m_num_cols;
  buffer_size+= KMF::params_out::num_1d_scalr * sizeof(Scalar)* m_num_cols;
  buffer_size+= KMF::params_out::num_2d_c * sizeof(Spack) * m_num_cols * nlevm_packs;
  // buffer_size+= KMF::params_out::num_2d_intfc * sizeof(Spack) * m_num_cols * nlev_int_packs;

  // Add for Fortran place holders here 
  constexpr auto num_f_mid = (KMF::params_in::num_2d_f + KMF::params_out::num_2d_f); 
  buffer_size+= num_f_mid * sizeof(Real) * m_num_cols * m_num_levs;
  // buffer_size+= num_f_int * sizeof(Real) * m_num_cols * (m_num_levs+1);

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
  // const int nlev_int_packs = ekat::npack<Spack>(m_num_levs+1);

  // constexpr auto num_1d_intgr = KMF::params_in::num_1d_intgr + KMF::params_out::num_1d_intgr;
  constexpr auto num_1d_scalr = KMF::params_in::num_1d_scalr + KMF::params_out::num_1d_scalr;
  constexpr auto num_2d_midlv_c = KMF::params_in::num_2d_c + KMF::params_out::num_2d_c;
  constexpr auto num_2d_midlv_f = KMF::params_in::num_2d_f + KMF::params_out::num_2d_f;
  // Need to separate out mid level and interface levels?
  
  Scalar* scl_mem = reinterpret_cast<Scalar*>(buffer_manager.get_memory());
  //----------------------------------------------------------------------------
  // device 1D integer variables
  KMF::view_1d<Scalar>* ptrs_1d_scalr[num_1d_scalr]             = { &params_out.precl };
  for (auto& v : ptrs_1d_scalr) {
    *v = KMF::view_1d<Scalar>(scl_mem, m_num_cols);
    scl_mem += v->size();
  }
  //----------------------------------------------------------------------------
  Real* r1_mem = reinterpret_cast<Real*>(scl_mem);
  //----------------------------------------------------------------------------
  // device 1D scalar scalars
  KMF::uview_1d<Real>* ptrs_1d_real[num_1d_scalr]          = { &params_out.f_precl};
  for (auto& v : ptrs_1d_real) {
    *v = KMF::uview_1d<Real>(r1_mem, m_num_cols);
    r1_mem += v->size();
  } 
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  Real* r_mem = reinterpret_cast<Real*>(r1_mem);
  //----------------------------------------------------------------------------
  // 2D "f_" views
  KMF::uview_2dl<Real>* midlv_f_ptrs[num_2d_midlv_f]  = { &params_in.f_cpair, 
                                                      &params_in.f_rair, 
                                                      &params_in.f_rho, 
                                                      &params_in.f_z, 
                                                      &params_in.f_pk,
                                                      &params_out.f_theta,
                                                      &params_out.f_qv,
                                                      &params_out.f_qc,
                                                      &params_out.f_qr,
                                                      &params_out.f_relhum
                                                    };
  for (int i=0; i<num_2d_midlv_f; ++i) {
    *midlv_f_ptrs[i] = KMF::view_2dl<Real>(r_mem, m_num_cols, m_num_levs);
    r_mem += midlv_f_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  Spack* spk_mem = reinterpret_cast<Spack*>(r_mem);
  //----------------------------------------------------------------------------
  // 2D views 
  KMF::view_2d<Spack>* midlv_c_ptrs[num_2d_midlv_c]  = { &params_in.rho, 
                                                      &params_in.z, 
                                                      &params_in.pk,
                                                      &params_out.theta,
                                                      &params_out.qv,
                                                      &params_out.qc,
                                                      &params_out.qr,
                                                      &params_out.relhum
                                                    };
  for (int i=0; i<num_2d_midlv_c; ++i) {
    *midlv_c_ptrs[i] = KMF::view_2d<Spack>(spk_mem, m_num_cols, m_num_levs);
    spk_mem += midlv_c_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  Real* total_mem = reinterpret_cast<Real*>(spk_mem);
  size_t used_mem = (reinterpret_cast<Real*>(total_mem) - buffer_manager.get_memory())*sizeof(Real);
  auto mem_chk = ( used_mem == requested_buffer_size_in_bytes() );
  EKAT_REQUIRE_MSG(mem_chk,"Error! Used memory != requested memory for Kessler.");
}
} // namespace scream