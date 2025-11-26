#include "eamxx_stengelP_process_interface.hpp"
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
  using namespace stengelP;
// =========================================================================================
//  Inputs (these are inherited from AtomoshpereProcess which means we can use the same logger):
//      comm - an EKAT communication group
//      params - a parameter list of options for the process.

StengelP::StengelP (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here usually
  m_atm_logger->info("[EAMxx] StengelP processes constructor");

  // if we set any starting values for the code we can get them here with:
  // var = params.get<std::type>("var_name");
}

// =========================================================================================
void StengelP::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  // using namespace ekat::units;
  // using namespace ShortFieldTagsNames;

  m_atm_logger->info("[EAMxx] StengelP processes set grids");

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

}

// =========================================================================================
void StengelP::initialize_impl (const RunType /* run_type */)
{
  // NOTE: run_type tells us if this is an initial or restarted run,

  // Set any universal things such as constants or masks (see Pompei example)

  m_atm_logger->info("[EAMxx] stengelP processes initialize_impl: ");

  #ifdef EAMXX_HAS_PYTHON
      if (has_py_module()) {
        try {
          py_module_call("init");
        } catch (const pybind11::error_already_set& e) {
          std::cout << "[stengelP::initialize_impl] Error! Something went wrong while calling the python module's function 'init'.\n"
                      " - module name: " + m_params.get<std::string>("py_module_name") + "\n"
                      " - pybind11 error: " + std::string(e.what()) + "\n";
          throw e;
        }

      }
    #endif

  // This is where we can setup Kokkos functions 
  // This is where we setup any starting physics
}

// =========================================================================================

// run_impl is called every timestep and where all of the physics happens
// Inputs:
//    - dt - the timestep for the current run step

void StengelP::run_impl (const double /* dt */)
{
  // Pull in variables .
  auto T_mid   = get_field_out("T_mid");
  auto p_mid   = get_field_out("p_mid"); // work?

  // TODO - scale and print 
  m_atm_logger->info("[EAMxx] stengelP run_impl: ");

  // Get log file name 
  const char* logname = m_atm_logger->get_logfile_name().c_str();
  
  #ifdef EAMXX_HAS_PYTHON
    if (has_py_module()) {
      // For now, we run Python code only on CPU
      const auto& T_mid          = get_py_field_host("T_mid");
      const auto& p_mid          = get_py_field_host("p_mid");

      // Sync input to host
      //field1_in.sync_to_host();

      // Read in any parameters
      // double param_name = m_params.get<double>("param_name");

      try {
        py_module_call("main", p_mid, T_mid, &logname);
      } catch (const pybind11::error_already_set& e) {
        std::cout << "[stengelP::run_impl] Error! Something went wrong while calling the python module's function 'main'.\n"
                    " - module name: " + m_params.get<std::string>("py_module_name") + "\n"
                    " - pybind11 error: " + std::string(e.what()) + "\n";
        throw e;
      }

      // Sync outputs to dev
      T_mid.sync_to_dev();
      p_mid.sync_to_dev();
    } else
    #endif

}

// =========================================================================================
void StengelP::finalize_impl()
{
  // Do nothing
  m_atm_logger->info("[EAMxx] StengelP processes clean up.");
}
// =========================================================================================

} // namespace scream
