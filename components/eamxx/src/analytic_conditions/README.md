# EAMxx Analytic Conditions

This directory contains EAMxx `AtmosphereProcess` subclasses that analytically
generate initial or boundary conditions for idealized test cases.  They are
**not** physics parameterizations and are not intended for production climate
runs.  Their sole purpose is to let EAMxx run standard dynamical-core or
microphysics test cases without requiring an input NetCDF file.

## Architecture

Each test case lives in its own subdirectory and builds a CMake library.  A
compile-time preprocessor macro (`EAMXX_HAS_<NAME>`) gates the code, and the
process is registered with the `AtmosphereProcessFactory` in
`register_analytic_conditions.hpp`.  That registration function is called by
every EAMxx entry point (MCT coupling, Python bindings) before any processes
are constructed.

```
analytic_conditions/
├── register_analytic_conditions.hpp   # factory registration (include in entry points)
├── CMakeLists.txt                     # calls add_subdirectory for each case
└── dcmip2016/                         # DCMIP2016 Test 1 moist baroclinic wave
    ├── dcmip2016_functions.hpp        # BaroclinicWaveFunctions<S,D> declarations
    ├── dcmip2016_functions_impl.hpp   # template function definitions
    ├── dcmip2016.cpp                  # explicit template instantiation
    ├── eamxx_dcmip2016_ic_process_interface.hpp   # AtmosphereProcess subclass header
    ├── eamxx_dcmip2016_ic_process_interface.cpp   # create_requests / initialize_impl
    └── CMakeLists.txt
```

### Adding a new test case

1. Create a subdirectory under `analytic_conditions/`.
2. Implement an `AtmosphereProcess` subclass that declares its output fields as
   `Computed` in `create_requests()` and fills them in `initialize_impl()`.
3. Build a library and set a `EAMXX_HAS_<NAME>` compile definition in the
   subdirectory `CMakeLists.txt`, and add `add_subdirectory(<name>)` to
   `analytic_conditions/CMakeLists.txt`.
4. Register the process in `register_analytic_conditions.hpp` under the
   corresponding `#ifdef EAMXX_HAS_<NAME>` guard.

---

## DCMIP2016 Test 1 — Moist Baroclinic Wave

**Reference:** Ullrich, Melvin, Staniforth & Jablonowski (2015), *QJRMS*
141(686): "A proposed baroclinic wave test case for deep and shallow-atmosphere
dynamical cores", doi:10.1002/qj.2544.

**Original Fortran:**
`components/homme/src/test_src/dcmip2016-baroclinic.F90`

### What it does

`DCMIP2016BaroclinicIC` (`dcmip2016_baroclinic_wave_ic` in the process
factory) analytically fills the following fields during `initialize_impl`:

| Field         | Units   | Description                         |
|---------------|---------|-------------------------------------|
| `T_mid`       | K       | Temperature at layer midpoints      |
| `horiz_winds` | m/s     | Zonal and meridional winds          |
| `ps`          | Pa      | Surface pressure (uniform p₀)       |
| `phis`        | m²/s²   | Surface geopotential (zero; flat)   |
| `qv`          | kg/kg   | Water-vapor specific humidity       |
| `qc`          | kg/kg   | Cloud liquid (initialized to zero)  |
| `qr`          | kg/kg   | Rain water (initialized to zero)    |

Pressure at each layer midpoint is computed from the hybrid-pressure
coefficients (`hyam`, `hybm`) stored in the grid geometry.  Altitude at each
pressure level is found via a secant-method iteration on the analytic
pressure–temperature profile (`eval_z_from_p`).

All computation runs on the model's device (GPU or CPU) via
`BaroclinicWaveFunctions<Real, DefaultDevice>::main()`, which launches a
`Kokkos::RangePolicy` kernel (one thread per column).

### Grid selection

The process automatically selects the correct grid for its `Computed` field
declarations so that EAMxx's sequential-splitting mechanism eliminates the
need for an IC file:

| Configuration                      | Grid used for ICs  |
|------------------------------------|--------------------|
| Standalone physics (no dynamics)   | `physics`          |
| HOMME with GLL physics             | `physics_gll`      |
| HOMME with FV physics (PG2)        | `physics_gll`      |

For the HOMME+FV physics (PG2) case, this matches the standard IC file
workflow: ICs are declared on the GLL grid, and HOMME's
`fv_phys_dyn_to_fv_phys` remapper copies them to PG2 during its own
`initialize_impl`.  Because this process runs first in `atm_procs_list`, the
`AtmosphereProcessGroup` sequential-splitting logic removes the GLL fields from
the group's `get_fields_in()`, so `set_initial_conditions` in the driver never
requests them from a file.

### Surface flux fields when omitting the coupler

When `sc_import`/`sc_export` are removed from `atm_procs_list`, surface flux
fields have no provider.  They must be constant-initialized in the
`initial_conditions` block — but **only if a process in your list actually
requires them**.  The consumers are:

| Field | Requiring process |
|---|---|
| `surf_sens_flux` | SHOC |
| `surf_mom_flux` | SHOC |
| `surf_evap` | SHOC (Updated) |
| `sfc_alb_dir_vis`, `sfc_alb_dir_nir`, `sfc_alb_dif_vis`, `sfc_alb_dif_nir` | RRTMGP, MAM microphysics |
| `surf_lw_flux_up` | RRTMGP |

For the minimal DCMIP2016+Kessler test none of these processes are present, so
no surface flux constants are needed.  If you later add SHOC or RRTMGP,
add the corresponding entries shown in the extended YAML below.

### Using it (no IC file required)

Place `dcmip2016_baroclinic_wave_ic` **first** in `atm_procs_list` so it fills
all fields before dynamics or physics attempt to read them.

**Minimal configuration (DCMIP2016 + Kessler only):**

```yaml
atmosphere_processes:
  schedule_type: sequential
  atm_procs_list: [dcmip2016_baroclinic_wave_ic, homme, kessler]

initial_conditions:
  # No 'filename' entry — dcmip2016_baroclinic_wave_ic provides the ICs.

  # phis = 0 is analytically correct for this flat-surface test, but it must
  # be listed here so any grid instance not covered by the IC process
  # (e.g., the PG2 phis that HOMME updates) is also initialized.
  phis: 0.0
```

**Extended configuration (adding SHOC and RRTMGP):**

```yaml
atmosphere_processes:
  schedule_type: sequential
  atm_procs_list: [dcmip2016_baroclinic_wave_ic, homme, shoc, rrtmgp, kessler]

initial_conditions:
  phis: 0.0

  # SHOC requires these when sc_import is absent:
  surf_sens_flux:  0.0
  surf_evap:       0.0
  surf_mom_flux:   [0.0, 0.0]

  # RRTMGP requires these when sc_import is absent:
  sfc_alb_dir_vis: 0.07
  sfc_alb_dir_nir: 0.07
  sfc_alb_dif_vis: 0.07
  sfc_alb_dif_nir: 0.07
  surf_lw_flux_up: 0.0
```

### Implementation files

| File | Purpose |
|------|---------|
| `dcmip2016_functions.hpp` | `BaroclinicWaveFunctions<S,D>` struct: type aliases, `Params` constants, `State` struct, and function declarations. |
| `dcmip2016_functions_impl.hpp` | Template definitions for `eval_pressure_temperature`, `eval_z_from_p`, `eval_exponential`, `eval_streamfunction`, `wave_at_point`, and `main`. |
| `dcmip2016.cpp` | Explicit template instantiation: `template struct BaroclinicWaveFunctions<Real, DefaultDevice>;` |
| `eamxx_dcmip2016_ic_process_interface.hpp` | `DCMIP2016BaroclinicIC` class declaration. |
| `eamxx_dcmip2016_ic_process_interface.cpp` | `create_requests()` (field declarations) and `initialize_impl()` (kernel launch + timestamp updates). |

### Key implementation notes

- **Double-precision arithmetic**: All intermediate computations in
  `dcmip2016_functions_impl.hpp` use `double` regardless of `ScalarT`.  This
  ensures the secant iteration in `eval_z_from_p` converges even in
  single-precision (`float`) builds.

- **Pack views, scalar math**: The process interface accepts `Pack**` views
  from the field manager (EAMxx's SIMD convention), but the analytic formulas
  cannot be level-vectorized (each level depends on the altitude found at that
  level by the secant method).  Inside the kernel, column slices are
  scalarized via `ekat::scalarize(ekat::subview(...))` and levels are looped
  sequentially.

- **Timestamps**: `initialize_impl` explicitly stamps every output field with
  `start_of_step_ts()` (i.e., t₀) after the kernel call.  This is required
  because `AtmosphereProcess::initialize` does not call `update_time_stamps`,
  so without the explicit stamp, t=0 output managers and precondition checks in
  subsequent processes would see invalid timestamps.

- **Restart safety**: `initialize_impl` returns immediately on
  `RunType::Restart`; the fields are already loaded from the restart file by
  the driver before this function is called.
