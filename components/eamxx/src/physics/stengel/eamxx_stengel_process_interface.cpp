#include "eamxx_stengel_process_interface.hpp"
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
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  // auto nondim = Units::nondimensional();
  // constexpr auto Pa = ekat::units::Pa;

  // specify which grid to use
  m_grid = grids_manager->get_grid("physics");
  const auto& grid_name = m_grid->name();
  // m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  // m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces
  // FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };

  // Set of fields used strictly as input (Required)
  // constexpr int ps = Pack::n;
  // add_tracer<Required>("tracer1_in", m_grid, kg/kg, ps); // tracers are for advection I think
  // add_field<Required>("qi", scalar3d_layout_mid, nondim, grid_name,ps);

  // Set of fields used strictly as output (Computed)
  // add_field<Computed>("field1_out", scalar3d_layout_mid, Pa, grid_name,ps);
  // add_field<Computed>("field2_out", scalar3d_layout_mid, Pa, grid_name,ps); 

  // Set of fields used as input and output (Updated)
  // add_field<Updated>("field1_updated", scalar3d_layout_mid, Pa,grid_name, ps);

  // FieldLayout scalar3d_layout_mid { {COL,LEV},     {m_num_cols,    m_num_levs  } };
  // constexpr int ps = Spack::n;
  // add_field<Required>("p_mid", scalar3d_layout_mid,  Pa,     grid_name, ps);


  // Gather parameters from parameter list:
  // stengel_var1 = m_params.get<double>("stengel_var1",1e-12);  // Default = 1e-12
}

// =========================================================================================
void Stengel::initialize_impl (const RunType /* run_type */)
{
  // NOTE: run_type tells us if this is an initial or restarted run,

  // Set any universal things such as constants or masks (see Pompei example)

  // Set property checks for fields in this process if needed
  // using Interval = FieldWithinIntervalCheck;
  // add_postcondition_check<Interval>(get_field_out("field1_out"),m_grid,0.0,1.0,false);
  // add_postcondition_check<Interval>(get_field_out("field2_out"),m_grid,0.0,1.0,false);
  // add_postcondition_check<Interval>(get_field_out("field1_updated"),m_grid,0.0,2.0,false);

  m_atm_logger->info("[EAMxx] stengel processes initialize_impl: ");

  // This is where we can setup Kokkos functions 
  // This is where we setup any starting physics
}

// =========================================================================================

// run_impl is called every timestep and where all of the physics happens
// Inputs:
//    - dt - the timestep for the current run step

void Stengel::run_impl (const double /* dt */)
{
  // Pull in variables .
  // auto p_mid   = get_field_in("p_mid");
  // field1_in, field1_out, field2_out, field1_updated
  // auto field1_in = get_field_in("field1_in");
  // auto field1_out = get_field_out("field1_out");
  // auto field2_out = get_field_out("field2_out"); 
  // auto field1_updated = get_field_out("field1_updated"); // use this for updated fields too 

  // // TODO - scale and print 
  m_atm_logger->info("[EAMxx] stengel run_impl: ");
  // auto field1_updated_max = field_max<Real>(field1_updated);
  // m_atm_logger->info("\t max value for updated field 1: "+ std::to_string(field1_updated_max));

  // field1_updated_max.scale(2.0);
  // field1_updated_max = field_max<Real>(field1_updated);
  // m_atm_logger->info("\t max value for updated field 1 scaled by 2x: "+ std::to_string(field1_updated_max));
  
}

// =========================================================================================
void Stengel::finalize_impl()
{
  // Do nothing
  m_atm_logger->info("[EAMxx] Stengel processes clean up.");
}
// =========================================================================================

} // namespace scream
