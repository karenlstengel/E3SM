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
eamxx:
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
eamxx:
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

(Both examples nest under the top-level `eamxx:` section, as CIME's buildnml
generates from `namelist_defaults_eamxx.xml` — see any file under
`components/eamxx/tests/*/input.yaml` for other worked examples of this
nesting convention.)

### Configuring the test case

All test-case knobs and reference constants are runtime params, read from the
process's YAML/XML entry (`m_params.get<T>("name", default)` in
`initialize_impl`).  Defaults are set in `namelist_defaults_eamxx.xml` under
`dcmip2016_baroclinic_wave_ic`; the C++ code's own fallback (used only if a
value is absent from *both* the YAML and the namelist defaults, e.g. in a
hand-written standalone YAML) is noted below.

| Param | Meaning | Default | C++ fallback if unset |
|---|---|---|---|
| `deep` | Deep (1) or shallow (0) atmosphere | `0` | `0` |
| `moist` | Moist (1) or dry (0) | `1` | `1` |
| `pertt` | Perturbation type: 0=exponential, 1=stream-function | `0` | `0` |
| `X` | Earth reduced-size scaling factor | `1.0` | `1.0` |
| `rearth` | Earth radius (m) | `6371220.0` | `physics::Constants<Real>::r_earth` (`6.376e6`) |
| `Rd` | Dry-air gas constant (J/kg/K) | `287.0423113650487` | `physics::Constants<Real>::Rair` (`287.042`) |
| `Rvap` | Water-vapor gas constant (J/kg/K) | `461.5046398201599` | `physics::Constants<Real>::RH2O` (`461.505`) |

The `deep`/`moist`/`pertt`/`X` defaults match the DCMIP2016 Test 1 protocol
(Ullrich et al. 2015) and the hardcoded values this process originally used.

The `rearth`/`Rd`/`Rvap` **namelist defaults** intentionally do *not* match
DCMIP2016's canonical literals (`a = 6.376e6`, `Rd = 287.042`,
`Rvap = 461.505` — see "Reference constants" below) — they instead match the
values Storm_SPEED's `moist_baroclinic_wave_dcmip2016` uses by default (CAM's
`physconst` module, itself derived from `shr_const_mod`), so that an
out-of-the-box EAMxx run and an out-of-the-box Storm_SPEED run of this test
use the same background state.

There is no automatic switch between the two conventions (no CIME
compset/testmod distinguishes them yet) — override `rearth`/`Rd`/`Rvap`
explicitly to get the other one. Both variants below are otherwise identical
to the minimal configuration above.

**Storm_SPEED-matching (default — no override needed):**

```yaml
eamxx:
  schedule_type: sequential
  atm_procs_list: [dcmip2016_baroclinic_wave_ic, homme, kessler]

initial_conditions:
  phis: 0.0
```

or, spelled out explicitly instead of relying on the namelist defaults:

```yaml
eamxx:
  schedule_type: sequential
  atm_procs_list: [dcmip2016_baroclinic_wave_ic, homme, kessler]
  dcmip2016_baroclinic_wave_ic:
    rearth: 6371220.0            # shr_const_rearth
    Rd:     287.0423113650487    # shr_const_rdair
    Rvap:   461.5046398201599    # shr_const_rwv

initial_conditions:
  phis: 0.0
```

**DCMIP2016-canonical (Ullrich et al. 2015 / `dcmip2016-baroclinic.F90` literals):**

```yaml
eamxx:
  schedule_type: sequential
  atm_procs_list: [dcmip2016_baroclinic_wave_ic, homme, kessler]
  dcmip2016_baroclinic_wave_ic:
    rearth: 6.376e6
    Rd:     287.042
    Rvap:   461.505

initial_conditions:
  phis: 0.0
```

Equivalently, after case creation, `atmchange` each param instead of editing
the YAML directly:

```bash
./atmchange dcmip2016_baroclinic_wave_ic::rearth=6.376e6
./atmchange dcmip2016_baroclinic_wave_ic::Rd=287.042
./atmchange dcmip2016_baroclinic_wave_ic::Rvap=461.505
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

- **Virtual potential temperature (thetav) is not computed.** The reference
  Fortran (`dcmip2016-baroclinic.F90`'s `baroclinic_wave_test`, and
  Storm_SPEED's `Tv_given_z`-based path) also derives
  `thetav = T*(1+Mvap*qv)*(p0/p)**(Rd/cp)`, but EAMxx has no consumer for a
  `thetav` field on the physics grid (T_mid is the state variable HOMME's
  theta-l dycore itself derives its internal theta from), so this port omits
  it rather than add an orphaned output.  `Params::cp` is kept in
  `dcmip2016_functions.hpp` for exactly this purpose, in case a future
  diagnostic needs it — see the comment above `Params::cp`. If you do want it,
  add a `thetav` `Computed` field in `create_requests()` (scalar3d_mid, K)
  and compute it at the end of `wave_at_point` using the (already-available)
  `Rd`/`p`/`p0` and a `moist_Rd_cp = Rd / P::cp` term, storing the result in a
  new `State::thetav` member.

### Reference constants

DCMIP2016's paper (Ullrich et al. 2015) specifies exact values for the Earth
radius and the dry-air/water-vapor gas constants, which HOMME's standalone
`dcmip2016-baroclinic.F90` hardcodes: `a = 6.376e6` m, `Rd = 287.04` J/kg/K,
`Rvap = 461.50` J/kg/K. This EAMxx port defaults to the same physical
quantities via EAMxx's own `physics::Constants<Real>` (`r_earth = 6.376e6`,
`Rair = 287.042`, `RH2O = 461.505` — matching to the precision each header
happens to use), *when its `rearth`/`Rd`/`Rvap` YAML params are left unset
entirely*. But because `namelist_defaults_eamxx.xml` always supplies a value
for a registered process, an out-of-the-box run instead gets the
**Storm_SPEED-matching defaults** described above in "Configuring the test
case" — Storm_SPEED being a CAM-based implementation of this same test
(`moist_baroclinic_wave_dcmip2016`,
`src/dynamics/tests/initial_conditions/ic_baroclinic.F90`), which pulls
`rearth`/`rair`/`rh2o` from CAM's `physconst` module. `physconst`'s defaults
come from `shr_const_mod`'s `SHR_CONST_REARTH` (`6.37122e6` m — the standard
CESM/E3SM Earth radius, not DCMIP2016's `6.376e6`),  `SHR_CONST_RDAIR`
(`≈287.0423` J/kg/K), and `SHR_CONST_RWV` (`≈461.5046` J/kg/K). Storm_SPEED's
`dctest_baro_kessler.xml` use-case does not override any of these, so an
out-of-the-box Storm_SPEED DCMIP2016 baroclinic-wave run uses the real Earth
radius rather than the DCMIP-specified one; `Rd`/`Rvap` happen to agree with
the DCMIP paper's literals to within their stated precision.  Net effect: the
Earth-radius discrepancy (~0.02%) is the only one with any real physical
significance, and this port lets you pick either convention (or override
independently) via the `rearth`/`Rd`/`Rvap` params documented above.

### Known issue in Storm_SPEED's `ic_baroclinic.F90` (not present in this port)

While comparing this port's math against Storm_SPEED's
`moist_baroclinic_wave_dcmip2016` (`src/dynamics/tests/initial_conditions/ic_baroclinic.F90`),
we found an argument-order mismatch between `evaluate_streamfunction`'s
declaration and its call sites:

- Declaration (`ic_baroclinic.F90:606`): `FUNCTION evaluate_streamfunction(z, lon_local, lat_local)`
  — dummy-argument order is `(z, lon_local, lat_local)`.
- Call sites, inside `uv_given_z` (`ic_baroclinic.F90:563-568`):
  `evaluate_streamfunction(lon, lat ± dxepsilon, z)` — actual-argument order
  is `(lon, lat, z)`.

Fortran binds by position, so the actual longitude lands in the formal `z`
(used for the vertical taper), the actual latitude lands in `lon_local`
(used in `cos(lon_local - pertlon)`), and the actual altitude `z` (meters,
up to ~30 km) lands in `lat_local` (used inside `sin`/`cos` as if it were a
latitude in radians). This would produce physically meaningless wind
perturbations.

The bug is currently **dormant**: Storm_SPEED hardcodes `pertt = 0`
(`ic_baroclinic.F90:56`), so the exponential-perturbation path
(`evaluate_exponential`, which *does* have consistent argument order between
its declaration and call site) is the only one exercised; the
stream-function path (`pertt = 1`) is never called. It would need fixing in
Storm_SPEED before `pertt = 1` could be used there.

This EAMxx port does not have the analogous bug: `eval_streamfunction`'s
declaration and both call sites in `wave_at_point`
(`dcmip2016_functions_impl.hpp`) consistently use `(lon, lat, z)` throughout.
