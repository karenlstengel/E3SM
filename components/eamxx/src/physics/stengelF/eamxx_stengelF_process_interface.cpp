#include "eamxx_stengelF_process_interface.hpp"
#include "stengelF_eamxx_bridge.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "share/field/field_utils.hpp"

#include <ekat_assert.hpp>
#include <ekat_units.hpp>

#include <array>

#ifdef EAMXX_HAS_PYTHON
#include "share/atm_process/atmosphere_process_pyhelpers.hpp"
#endif

namespace scream
{
  using namespace stengelF;
// =========================================================================================
//  Inputs (these are inherited from AtomoshpereProcess which means we can use the same logger):
//      comm - an EKAT communication group
//      params - a parameter list of options for the process.

StengelF::StengelF (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here usually
  m_atm_logger->info("[EAMxx] StengelF processes constructor");

}

// =========================================================================================
void StengelF::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  // using namespace ekat::units;
  // using namespace ShortFieldTagsNames;

  m_atm_logger->info("[EAMxx] StengelF processes set grids");

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
}

// =========================================================================================
void StengelF::initialize_impl (const RunType /* run_type */)
{
  m_atm_logger->info("[EAMxx] stengelF processes initialize_impl: ");

  stengelF::stengelF_eamxx_bridge_init(m_num_cols, m_num_levs);
}

// =========================================================================================

// run_impl is called every timestep and where all of the physics happens
// Inputs:
//    - dt - the timestep for the current run step

void StengelF::run_impl (const double /* dt */)
{
  // Pull in variables .
  auto T_mid   = get_field_out("T_mid");
  auto p_mid   = get_field_out("p_mid");

  m_atm_logger->info("[EAMxx] stengelF run_impl: ");
  
  auto T_mid_max = field_max<Real>(T_mid);
  m_atm_logger->info("\t max value for updated T_mid: "+ std::to_string(T_mid_max));
  auto p_mid_max = field_max<Real>(p_mid);
  m_atm_logger->info("\t max value for updated p_mid: "+ std::to_string(p_mid_max));

  // Add the fields we pulled in to the params struct to pass to the bridge:
  params.p_mid = p_mid.get_view<Spack**, Host>();
  params.T_mid = T_mid.get_view<Spack**, Host>();

  stengelF_eamxx_bridge_run(m_num_cols, m_num_levs, params); 

  // Update with the new (should be the same) values from the run
  // T_mid = params.T_mid; // will have to see if this is the correct way to do this
  // p_mid = params.p_mid;

  // T_mid_max = field_max<Real>(T_mid);
  // m_atm_logger->info("\t max value for T_mid after the Fortran bridge: "+ std::to_string(T_mid_max));
  // p_mid_max = field_max<Real>(p_mid);
  // m_atm_logger->info("\t max value for p_mid after the Fortran bridge: "+ std::to_string(p_mid_max));
  
}

// =========================================================================================
void StengelF::finalize_impl()
{
  // Do nothing
  m_atm_logger->info("[EAMxx] StengelF processes clean up.");
}
// =========================================================================================

} // namespace scream
