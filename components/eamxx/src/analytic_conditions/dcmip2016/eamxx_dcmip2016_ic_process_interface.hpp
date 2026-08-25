#ifndef EAMXX_DCMIP2016_IC_PROCESS_INTERFACE_HPP
#define EAMXX_DCMIP2016_IC_PROCESS_INTERFACE_HPP

// EAMxx AtmosphereProcess that sets DCMIP2016 Test 1 analytic initial conditions.
//
// This class analytically initializes T_mid, horiz_winds, ps, phis, qv, qc, qr
// using the moist baroclinic wave formulas from Ullrich et al. (2015).  It is
// designed to let EAMxx run DCMIP2016 Test 1 without an IC file.
//
// Grid selection
// --------------
// The process declares its Computed fields on whichever grid is appropriate for
// the IC workflow:
//
//  - Standalone physics (no dynamics):   physics grid ("physics" or "point_grid")
//  - HOMME with GLL physics:             physics grid == "physics_gll"
//  - HOMME with FV physics (PG2 grid):   GLL grid ("physics_gll")
//
// The last case matches the standard HOMME IC workflow: ICs are read on the GLL
// grid and HOMME's fv_phys_dyn_to_fv_phys remapper copies them to PG2 during
// its own initialize_impl.  Because this process declares its output fields on
// the correct grid first, the AtmosphereProcessGroup's sequential-splitting
// logic removes those fields from get_fields_in() at the group level — so
// set_initial_conditions() in the driver never tries to load them from a file.
//
// Usage (no-IC-file YAML)
// -----------------------
//   atm_procs_list: [dcmip2016_baroclinic_wave_ic, homme, kessler]
//
//   initial_conditions:
//     # No 'filename' needed — dcmip2016_baroclinic_wave_ic provides the ICs.
//     # phis = 0 is analytically correct for the flat-surface DCMIP2016 test,
//     # but we still need it listed here so any grid that isn't covered by this
//     # process (e.g., the dynamics GLL grid when using PG2 physics) gets it.
//     phis: 0.0
//     # Surface fields must be given as constants when the coupler is omitted:
//     surf_sens_flux:  0.0
//     surf_evap:       0.0
//     surf_mom_flux:   [0.0, 0.0]
//     sfc_alb_dir_vis: 0.07
//     sfc_alb_dir_nir: 0.07
//     sfc_alb_dif_vis: 0.07
//     sfc_alb_dif_nir: 0.07
//     surf_lw_flux_up: 0.0
//
// This process must be listed FIRST so it fills the fields before dynamics
// or physics try to read them.  It is a no-op during run_impl; all work
// happens once in initialize_impl.
//
// The surface coupler (sc_import / sc_export) should be omitted from
// atm_procs_list.  Provide surface-flux fields as constants in the
// initial_conditions YAML block instead (see example above).

#include "share/atm_process/atmosphere_process.hpp"

#include <ekat_parameter_list.hpp>

#include <string>

namespace scream {

class DCMIP2016BaroclinicIC : public AtmosphereProcess
{
public:

  // Standard EAMxx constructor: comm and parameter list are forwarded to the
  // base class, which stores them.  No DCMIP2016-specific setup needed here.
  DCMIP2016BaroclinicIC(const ekat::Comm& comm, const ekat::ParameterList& params)
    : AtmosphereProcess(comm, params) {}

  // This is a physics-type process (not dynamics, coupler, or group)
  AtmosphereProcessType type() const override {
    return AtmosphereProcessType::Physics;
  }

  // String key used to look up and construct this process via the factory
  // (registered in register_analytic_conditions.hpp under EAMXX_HAS_DCMIP2016_IC).
  // Must match the atm_procs_list entry.
  std::string name() const override {
    return "dcmip2016_baroclinic_wave_ic";
  }

  // Declare which fields this process provides (Computed) and which
  // it needs from other processes or geometry (Required).
  void create_requests() override;

// Cuda/HIP require __device__ lambdas to be in public scope.
// The preprocessor guard matches the pattern used in other EAMxx processes.
#ifndef KOKKOS_ENABLE_CUDA
protected:
#endif

  // Fill Computed fields analytically.  Called once during model init.
  // On RunType::Restart the fields are already loaded from the restart file,
  // so we skip re-initialization to avoid overwriting them.
  void initialize_impl(const RunType run_type) override;

  // IC process: no work per timestep.
  void run_impl(const double dt) override {}

protected:
  void finalize_impl() override {}

  // Physics grid and dimensions, set in create_requests()
  std::shared_ptr<const AbstractGrid> m_grid;
  int m_num_cols = 0;
  int m_num_levs = 0;
};

} // namespace scream

#endif // EAMXX_DCMIP2016_IC_PROCESS_INTERFACE_HPP
