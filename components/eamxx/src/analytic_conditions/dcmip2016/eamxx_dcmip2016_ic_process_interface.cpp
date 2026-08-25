#include "analytic_conditions/dcmip2016/eamxx_dcmip2016_ic_process_interface.hpp"
#include "analytic_conditions/dcmip2016/dcmip2016_functions.hpp"

#include "share/util/eamxx_units.hpp"

namespace scream {

using namespace ekat::units;

// ==========================================================================
// create_requests
//
// Register all fields this process provides (Computed) so the field manager
// allocates them before initialize_impl runs.
//
// lat, lon, hyam, hybm come from grid geometry data and are fetched
// directly in initialize_impl — they do not go through the field manager.
// ==========================================================================
void DCMIP2016BaroclinicIC::create_requests()
{
  // When HOMME is used with FV physics (PG2 grid), the standard IC workflow reads
  // initial conditions from a NetCDF file on the GLL grid; the dynamics interface
  // then remaps them to PG2 during its own initialize_impl via fv_phys_dyn_to_fv_phys.
  // To match that workflow without requiring an IC file, we declare and fill our
  // analytic ICs on the GLL grid whenever it is distinct from the physics grid.
  // For standalone physics tests (no HOMME) or GLL-physics HOMME runs, we fall
  // back to the physics grid directly.
  auto phys_grid = m_grids_manager->get_grid("physics");
  const bool use_gll = (phys_grid->name() != "physics_gll") &&
                        m_grids_manager->has_grid("physics_gll");
  m_grid     = use_gll ? m_grids_manager->get_grid("physics_gll") : phys_grid;

  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs();
  m_num_levs = m_grid->get_num_vertical_levels();

  // Use the instantiated type to get Pack size consistently with the kernel
  using BWFunctions = dcmip2016::BaroclinicWaveFunctions<Real, DefaultDevice>;
  constexpr int ps = BWFunctions::Pack::n;

  auto scalar2d     = m_grid->get_2d_scalar_layout();
  auto scalar3d_mid = m_grid->get_3d_scalar_layout(LEV);
  auto vector3d_mid = m_grid->get_3d_vector_layout(LEV, 2);

  // Thermodynamic fields filled analytically in initialize_impl
  add_field<Computed>("T_mid",       scalar3d_mid, K,          grid_name, ps);
  add_field<Computed>("horiz_winds", vector3d_mid, m/s,        grid_name, ps);
  add_field<Computed>("ps",          scalar2d,     Pa,         grid_name);
  add_field<Computed>("phis",        scalar2d,     m*m/(s*s),  grid_name);

  // Tracers placed in the "tracers" FieldGroup for Homme and Kessler
  add_tracer<Computed>("qv", m_grid, kg/kg, ps);  // water vapor (analytic)
  add_tracer<Computed>("qc", m_grid, kg/kg, ps);  // cloud liquid (initialized to 0)
  add_tracer<Computed>("qr", m_grid, kg/kg, ps);  // rain water  (initialized to 0)
}

// ==========================================================================
// initialize_impl
//
// Fill all Computed fields by launching the DCMIP2016 analytic-IC kernel on
// the model's device (GPU or CPU) via BaroclinicWaveFunctions::main().
//
// On RunType::Restart the fields are already loaded from the restart file by
// the AtmosphereDriver before this function is called — return early to avoid
// overwriting them.
// ==========================================================================
void DCMIP2016BaroclinicIC::initialize_impl(const RunType run_type)
{
  if (run_type == RunType::Restart) return;

  using BWFunctions = dcmip2016::BaroclinicWaveFunctions<Real, DefaultDevice>;
  using Pack        = BWFunctions::Pack;

  // ---- Geometry: device views of lat/lon (degrees) and hybrid coefficients ----
  // These are plain Real* views — no Pack needed since geometry data is scalar.
  const auto lat_deg = m_grid->get_geometry_data("lat").get_view<const Real*>();
  const auto lon_deg = m_grid->get_geometry_data("lon").get_view<const Real*>();
  const auto hyam    = m_grid->get_geometry_data("hyam").get_view<const Real*>();
  const auto hybm    = m_grid->get_geometry_data("hybm").get_view<const Real*>();

  // ---- Output field views (device, Pack-based to match field manager storage) ----
  auto T_mid_v      = get_field_out("T_mid").get_view<Pack**>();
  auto horiz_winds_v = get_field_out("horiz_winds").get_view<Pack***>();
  auto ps_v          = get_field_out("ps").get_view<Real*>();
  auto phis_v        = get_field_out("phis").get_view<Real*>();
  auto qv_v          = get_field_out("qv").get_view<Pack**>();
  auto qc_v          = get_field_out("qc").get_view<Pack**>();
  auto qr_v          = get_field_out("qr").get_view<Pack**>();

  // ---- DCMIP2016 test parameters ----
  const int  deep  = 0;    // Shallow atmosphere
  const int  moist = 1;    // Include moisture (needed for Kessler microphysics)
  const int  pertt = 0;    // Exponential wind perturbation
  const Real X     = 1.0; // Full Earth (no scaling)

  // Launch kernel; all column-level work runs on DeviceT
  BWFunctions::main(
    m_num_cols, m_num_levs,
    deep, moist, pertt, X,
    lat_deg, lon_deg, hyam, hybm,
    T_mid_v, horiz_winds_v,
    ps_v, phis_v,
    qv_v, qc_v, qr_v);

  // Stamp every output field with t0.  Without this, t=0 output managers
  // would see an invalid timestamp and skip these fields, and NaN precondition
  // checks in subsequent processes would also see them as uninitialised.
  const auto& t0 = start_of_step_ts();
  get_field_out("T_mid"       ).get_header().get_tracking().update_time_stamp(t0);
  get_field_out("horiz_winds" ).get_header().get_tracking().update_time_stamp(t0);
  get_field_out("ps"          ).get_header().get_tracking().update_time_stamp(t0);
  get_field_out("phis"        ).get_header().get_tracking().update_time_stamp(t0);
  get_field_out("qv"          ).get_header().get_tracking().update_time_stamp(t0);
  get_field_out("qc"          ).get_header().get_tracking().update_time_stamp(t0);
  get_field_out("qr"          ).get_header().get_tracking().update_time_stamp(t0);
}

} // namespace scream
