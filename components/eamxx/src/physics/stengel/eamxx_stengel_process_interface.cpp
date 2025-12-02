#include "eamxx_stengel_process_interface.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "share/field/field_utils.hpp"

#include <ekat_assert.hpp>
#include <ekat_units.hpp>
#include <ekat_team_policy_utils.hpp>

#include <array>

namespace scream
{
  using namespace stengel;
// =========================================================================================
//  Inputs (these are inherited from AtomoshpereProcess which means we can use the same logger):
//      comm - an EKAT communication group
//      params - a parameter list of options for the process.

Stengel::Stengel (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here usually
  m_atm_logger->info("[EAMxx] Stengel processes constructor");

  // if we set any starting values for the code we can get them here with:
  // var = params.get<std::type>("var_name");
}

// =========================================================================================
void Stengel::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  // using namespace ekat::units;
  // using namespace ShortFieldTagsNames;

  m_atm_logger->info("[EAMxx] Stengel processes set grids");

  // // The units of mixing ratio Q are technically non-dimensional.
  // // Nevertheless, for output reasons, we like to see 'kg/kg'.
  // auto nondim = Units::nondimensional();
  constexpr auto K = ekat::units::K;
  constexpr auto Pa = ekat::units::Pa;

  // specify which grid to use
  m_grid = grids_manager->get_grid("physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar3d_layout_mid = m_grid->get_3d_scalar_layout(true);

  constexpr int ps = 1;
  add_field<Updated>("T_mid", scalar3d_layout_mid, K, grid_name, ps);
  add_field<Updated>("p_mid", scalar3d_layout_mid, Pa, grid_name, ps);

  // Gather parameters from parameter list:
  // stengel_var1 = m_params.get<double>("stengel_var1",1e-12);  // Default = 1e-12
}

// =========================================================================================
void Stengel::initialize_impl (const RunType /* run_type */)
{
  // NOTE: run_type tells us if this is an initial or restarted run,

  // Set any universal things such as constants or masks (see Pompei example)

  m_atm_logger->info("[EAMxx] stengel processes initialize_impl: ");

  // This is where we setup any starting physics
}

// =========================================================================================

// run_impl is called every timestep and where all of the physics happens
// Inputs:
//    - dt - the timestep for the current run step

void Stengel::run_impl (const double /* dt */)
{

  m_atm_logger->info("[EAMxx] stengel run_impl: ");
  
  // Run with kokkos
  #ifdef EAMXX_ENABLE_GPU

    // Pull in variables .
    auto T_mid   = get_field_out("T_mid").get_view<Spack**>();
    auto p_mid   = get_field_out("p_mid").get_view<Spack**>(); // use this for updated fields too 

    m_atm_logger->info("\t running on GPU! ");
    // TODO - need to open fields as view?
    Real T_mid_max = 0.0;
    Real p_mid_max = 0.0;

    // Kokkos function to get the max value in the field
    StengelFunc::get_max_value(T_mid, T_mid_max);
    m_atm_logger->info("\t\t max value for updated T_mid: "+ std::to_string(T_mid_max));

    // Kokkos function to scale field by scalar value
    StengelFunc::scale_field_2d(T_mid, 2.0);
    StengelFunc::get_max_value(T_mid, T_mid_max);
    m_atm_logger->info("\t\t max value for T_mid scaled by 2x: "+ std::to_string(T_mid_max));

    StengelFunc::scale_field_2d(T_mid, 0.5);
    StengelFunc::get_max_value(T_mid, T_mid_max);
    m_atm_logger->info("\t\t max value for T_mid scaled by 0.5x: "+ std::to_string(T_mid_max));

    //------------
    // Repeat for p_mid
    StengelFunc::get_max_value(p_mid, p_mid_max);
    m_atm_logger->info("\t\t max value for updated p_mid: "+ std::to_string(p_mid_max));

    // Kokkos function to scale field by scalar value
    StengelFunc::scale_field_2d(p_mid, 2.0);
    StengelFunc::get_max_value(p_mid, p_mid_max);
    m_atm_logger->info("\t\t max value for p_mid scaled by 2x: "+ std::to_string(p_mid_max));

    StengelFunc::scale_field_2d(p_mid, 0.5);
    StengelFunc::get_max_value(p_mid, p_mid_max);
    m_atm_logger->info("\t\t max value for p_mid scaled by 0.5x: "+ std::to_string(p_mid_max));

  #else
    m_atm_logger->info("\t running on CPU! ");
    // Pull in variables .
    auto T_mid   = get_field_out("T_mid");
    auto p_mid   = get_field_out("p_mid"); // use this for updated fields too 

    // TODO - change T_mid value here and print like above
    auto T_mid_max = field_max<Real>(T_mid);
    m_atm_logger->info("\t\t max value for updated T_mid: "+ std::to_string(T_mid_max));

    // Can't change value of p_mid directly since p_mid is read-only
    T_mid.scale(2.0);
    T_mid_max = field_max<Real>(T_mid);
    m_atm_logger->info("\t\t max value for T_mid scaled by 2x: "+ std::to_string(T_mid_max));

    T_mid.scale(0.5);
    T_mid_max = field_max<Real>(T_mid);
    m_atm_logger->info("\t\t max value for T_mid scaled by 0.5x: "+ std::to_string(T_mid_max));

    // ------------------ try with p_mid
    auto p_mid_max = field_max<Real>(p_mid);
    m_atm_logger->info("\t\t max value for updated p_mid: "+ std::to_string(p_mid_max));

    // Can't change value of p_mid directly since p_mid is read-only?
    p_mid.scale(2.0);
    p_mid_max = field_max<Real>(p_mid);
    m_atm_logger->info("\t\t max value for p_mid scaled by 2x: "+ std::to_string(p_mid_max));

    p_mid.scale(0.5);
    p_mid_max = field_max<Real>(p_mid);
    m_atm_logger->info("\t\t max value for p_mid scaled by 0.5x: "+ std::to_string(p_mid_max));
  #endif 
}

// =========================================================================================
void Stengel::finalize_impl()
{
  // Do nothing
  m_atm_logger->info("[EAMxx] Stengel processes clean up.");
}
// =========================================================================================

} // namespace scream
