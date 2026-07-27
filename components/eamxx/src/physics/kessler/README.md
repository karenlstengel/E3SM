# EAMxx Kessler Bridge

This directory contains the EAMxx C++ atmosphere process interface for the Kessler
microphysics scheme.  The Fortran bridge code that connects EAMxx to the Kessler
Fortran scheme is **automatically generated** by the
[xdsl-cpp](https://github.com/xdsl-project/xdsl-cpp) framework (CCPP dialect) rather
than written by hand.  This document describes the full generation process so it can
be reproduced or updated.

## Repository layout

```
kessler/
  eamxx_kessler_process_interface.hpp   C++ class declaration (AtmosphereProcess)
  eamxx_kessler_process_interface.cpp   C++ class implementation (calls generated bridge)
  kessler_functions.hpp                 Kokkos/GPU helper structs (params_helpers, params_computed)
  generated_bridge/                     xdsl-cpp generated files (checked in, do not edit by hand)
    Kessler_ccpp_chost_cap.F90          Fortran BIND(C) entry points
    Kessler_ccpp_cap.F90                Internal CCPP orchestration cap (has OpenACC directives)
    kessler_suite_cap.F90               Suite-level cap (calls scheme subroutines in order)
    ccpp_kinds.F90                      Fortran kind definitions
    Kessler_ccpp_chost_cap.h            C extern "C" declarations (included by .cpp)
    Kessler_chost.hpp                   C++ ergonomics wrapper (struct-based, optional)
    ccpp_kinds.h                        C++ kind/type mappings
  CMakeLists.txt
  README.md                             This file
```

The handwritten bridge files previously in `fortran_bridge/` have been superseded by
the generated files in `generated_bridge/`.  The `fortran_bridge/` directory is kept
for reference but its files are **not compiled**.

The Kessler Fortran scheme and xdsl-cpp metadata live outside this tree:

```
<repo-root>/
  atmospheric_physics/schemes/kessler/
    kessler.F90
    kessler_update.F90
  xdsl-cpp/
    examples/kessler/
      scheme/
        kessler.meta               CCPP metadata for kessler_init / kessler_run
        kessler_update.meta        CCPP metadata for kessler_update_* entry points
        kessler_suite.xml          Suite definition (kessler -> kessler_update)
      host_eamxx/                  EAMxx host metadata (created as part of this work)
        eamxx_kessler_host_mod.meta
        eamxx_kessler_host_sub.meta
      bindc_eamxx_acc/             Raw generator output (copied into generated_bridge/)
```

---

## How the bridge code was generated

### Step 1 — understand the scheme metadata

The Kessler Fortran entry points are described in two `.meta` files that live alongside
the Fortran source in the xdsl-cpp examples directory:

- `xdsl-cpp/examples/kessler/scheme/kessler.meta` — `kessler_init`, `kessler_run`
- `xdsl-cpp/examples/kessler/scheme/kessler_update.meta` — `kessler_update_init`,
  `kessler_update_timestep_init`, `kessler_update_run`, `kessler_update_timestep_final`

Each argument in those files carries a `standard_name`, `units`, `dimensions`, `type`,
`intent`, and (for GPU fields) `memory_space = device`.  These standard names are the
vocabulary the framework uses to match scheme arguments to host-provided variables.

### Step 2 — create the EAMxx host metadata files

The xdsl-cpp framework requires two metadata files describing the *host model* side:

| File | Purpose |
|------|---------|
| `eamxx_kessler_host_mod.meta` | Declares all host module variables (scalars and arrays) with their standard names |
| `eamxx_kessler_host_sub.meta` | Declares the host subroutine entry point and loop-control variables |

These files were written from scratch for EAMxx and are stored at:

```
xdsl-cpp/examples/kessler/host_eamxx/eamxx_kessler_host_mod.meta
xdsl-cpp/examples/kessler/host_eamxx/eamxx_kessler_host_sub.meta
```

#### `eamxx_kessler_host_mod.meta` — key decisions

The file declares every variable that the Kessler schemes need, using:
- `language = c++` so the generator emits C-compatible `extern "C"` bridge code
- C++ variable names that match EAMxx conventions where they differ from the toy
  host example (e.g. `z_mid` not `z`, `temp_tend` not `ttend_t`, `exner` not `pk`)
- `standard_name` values taken verbatim from the scheme `.meta` files so the
  `generate-host-match` pass can resolve every argument

The full variable set and their standard names:

| C++ name | Standard name | Used by |
|----------|---------------|---------|
| `ncol` | `horizontal_dimension` | all |
| `nz` | `vertical_layer_dimension` | all |
| `dt` | `timestep_for_physics` | kessler_run, kessler_update_run |
| `lyr_surf` | `vertical_index_at_surface_adjacent_layer` | kessler_run |
| `lyr_toa` | `vertical_index_at_top_adjacent_layer` | kessler_run |
| `lv` | `latent_heat_of_vaporization_of_water_at_0c` | kessler_init |
| `pref` | `surface_reference_pressure` | kessler_init |
| `rhoqr` | `fresh_liquid_water_density_at_0c` | kessler_init |
| `gravit` | `standard_gravitational_acceleration` | kessler_update_init |
| `scheme_name` | `scheme_name` | kessler_run (output) |
| `cpair` | `composition_dependent_specific_heat_of_dry_air_at_constant_pressure` | kessler_run, kessler_update_timestep_final |
| `rair` | `composition_dependent_gas_constant_of_dry_air` | kessler_run |
| `rho` | `dry_air_density` | kessler_run |
| `z_mid` | `geopotential_height_wrt_surface` | kessler_run, kessler_update_timestep_final |
| `exner` | `dimensionless_exner_function` | kessler_run, kessler_update_run |
| `theta` | `air_potential_temperature` | kessler_run, kessler_update_run |
| `qv` | `water_vapor_mixing_ratio_wrt_dry_air` | kessler_run |
| `qc` | `cloud_liquid_water_mixing_ratio_wrt_dry_air` | kessler_run |
| `qr` | `rain_mixing_ratio_wrt_dry_air` | kessler_run |
| `precl` | `total_precipitation_rate_at_surface` | kessler_run |
| `relhum` | `relative_humidity` | kessler_run |
| `temp` | `air_temperature` | kessler_update_timestep_init/final |
| `temp_prev` | `air_temperature_on_previous_timestep` | kessler_update_* |
| `temp_tend` | `tendency_of_air_temperature_due_to_model_physics` | kessler_update_* |
| `phis` | `surface_geopotential` | kessler_update_timestep_final |
| `st_energy` | `dry_static_energy` | kessler_update_timestep_final |

Note on `horizontal_dimension` vs `horizontal_loop_extent`: the schemes use both
tokens.  `horizontal_dimension` comes from `ncol` in the host mod.
`horizontal_loop_extent` is derived by the framework as `col_end - col_start + 1`
from the host sub variables, and equals `ncol` in EAMxx since all columns are always
processed.

#### `eamxx_kessler_host_sub.meta` — key decisions

```ini
[ccpp-table-properties]
  name = eamxx_kessler_host_sub
  type = host
  language = c++
[ccpp-arg-table]
  name = eamxx_kessler_host_sub
  type = host
[ col_start ]
  standard_name = horizontal_loop_begin
  ...
[ col_end ]
  standard_name = horizontal_loop_end
  ...
[ errmsg ]
  standard_name = ccpp_error_message
  ...
[ errflg ]
  standard_name = ccpp_error_code
  ...
```

`col_start` and `col_end` are passed as `1` and `ncol` by the EAMxx C++ caller since
EAMxx always processes all columns in a single call.

### Step 3 — run the generator

The `ccpp_xdsl` command is installed as part of the xdsl-cpp Python package.  Run it
from the `xdsl-cpp/` directory:

```bash
cd <repo-root>/xdsl-cpp

ccpp_xdsl \
    --suites       examples/kessler/scheme/kessler_suite.xml \
    --scheme-files examples/kessler/scheme/kessler.meta,examples/kessler/scheme/kessler_update.meta \
    --host-files   examples/kessler/host_eamxx/eamxx_kessler_host_mod.meta,examples/kessler/host_eamxx/eamxx_kessler_host_sub.meta \
    --bind-c \
    --directive acc \
    -o examples/kessler/bindc_eamxx_acc
```

Flag summary:

| Flag | Effect |
|------|--------|
| `--suites` | Suite XML defining which schemes run and in what order |
| `--scheme-files` | Comma-separated list of scheme `.meta` files |
| `--host-files` | Comma-separated list of host `.meta` files |
| `--bind-c` | Emit Fortran `BIND(C)` caps and matching C `extern "C"` headers |
| `--directive acc` | Emit `!$acc` OpenACC data movement directives guarded by `#ifdef USE_GPU` |
| `-o` | Output directory |

Omitting `--directive acc` generates the same BIND(C) interface but without any
OpenACC directives (useful for CPU-only builds).

### Step 4 — generated files and what they do

The generator writes eight files:

| File | Description |
|------|-------------|
| `Kessler_ccpp_chost_cap.F90` | **Primary bridge.** Fortran module with six `BIND(C)` subroutines called directly from EAMxx C++: `register`, `initialize`, `finalize`, `timestep_initial`, `timestep_final`, `run` |
| `Kessler_ccpp_chost_cap.h` | Matching C `extern "C"` declarations; included by `eamxx_kessler_process_interface.cpp` |
| `Kessler_chost.hpp` | Optional C++ ergonomics wrapper (struct-based, not used by EAMxx directly) |
| `Kessler_ccpp_cap.F90` | Internal CCPP cap with OpenACC data movement directives; called by the chost cap |
| `kessler_suite_cap.F90` | Suite-level cap that calls `kessler_init/run` and `kessler_update_*` in the correct order |
| `ccpp_kinds.F90` | Fortran kind definitions (`kind_phys`) |
| `ccpp_kinds.h` | C++ kind/type mappings (e.g. `kind_phys` -> `double`) |

The generated BIND(C) entry points and their correspondence to the handwritten bridge:

| Generated entry point | Replaces handwritten call |
|-----------------------|--------------------------|
| `Kessler_chost_physics_register` | (new — no handwritten equivalent) |
| `Kessler_chost_physics_initialize(lv, pref, rhoqr, gravit, ...)` | `kessler_eamxx_bridge_init_c` + `kessler_eamxx_bridge_update_init_c` |
| `Kessler_chost_physics_timestep_initial(ncol, nz, temp, temp_prev, temp_tend, ...)` | First call inside `kessler_eamxx_bridge_update_c` |
| `Kessler_chost_physics_run(ncol, nz, col_start, col_end, dt, ...)` | `kessler_eamxx_bridge_run_c` + `kessler_update_run` call inside `kessler_eamxx_bridge_update_c` |
| `Kessler_chost_physics_timestep_final(ncol, nz, cpair, temp, z_mid, phis, st_energy, ...)` | Final call inside `kessler_eamxx_bridge_update_c` |
| `Kessler_chost_physics_finalize` | (new — no handwritten equivalent) |

After generation, the output directory was copied into this source tree:

```bash
cp -r xdsl-cpp/examples/kessler/bindc_eamxx_acc \
      E3SM/components/eamxx/src/physics/kessler/generated_bridge
```

---

## Changes made to EAMxx source files

### `CMakeLists.txt`

- Added `set(XDSL_GENERATED_PATH ${CMAKE_CURRENT_SOURCE_DIR}/generated_bridge)`
- Replaced the three handwritten bridge sources in `KESSLER_F90_SRCS`:
  ```
  fortran_bridge/kessler_eamxx_bridge.cpp        (removed)
  fortran_bridge/kessler_eamxx_bridge_main.F90   (removed)
  fortran_bridge/kessler_eamxx_bridge_update.F90 (removed)
  ```
  with the four generated Fortran files:
  ```
  ${XDSL_GENERATED_PATH}/Kessler_ccpp_chost_cap.F90
  ${XDSL_GENERATED_PATH}/kessler_suite_cap.F90
  ${XDSL_GENERATED_PATH}/Kessler_ccpp_cap.F90
  ${XDSL_GENERATED_PATH}/ccpp_kinds.F90
  ```
- Replaced `${PATH_TO_LEGACY_CAM_SIMA}/test/include/ccpp_kinds.F90` with the
  generated `${XDSL_GENERATED_PATH}/ccpp_kinds.F90`
- Replaced the `fortran_bridge/` include directory with `${XDSL_GENERATED_PATH}` in
  `target_include_directories`

### `eamxx_kessler_process_interface.cpp`

**Include:** `kessler_eamxx_bridge.hpp` replaced with `Kessler_ccpp_chost_cap.h`.

**`initialize_impl`:** single `kessler_eamxx_bridge_init` call replaced with:
```cpp
Kessler_chost_physics_register(errmsg, &errflg);
Kessler_chost_physics_initialize(latvap, P0, rhoqr, gravit, errmsg, &errflg);
```

**`run_impl`:** single `kessler_eamxx_bridge_run` call (which previously encapsulated
both the transpose logic and the Fortran calls inside `kessler_eamxx_bridge.cpp`)
replaced with the transpose sandwich and three generated cap calls:
```cpp
params_helpers.transpose<c2f>(m_num_cols, nlevs);
params_computed.transpose<c2f>(m_num_cols, nlevs);
Kokkos::fence();

// #if GPU && !OpenACC  ->  use h_* host mirror views
// #else                ->  use f_* Fortran-layout device views
Kessler_chost_physics_timestep_initial(ncol, nz, f_temp, f_temp_prev, f_temp_tend, ...);
Kessler_chost_physics_run(ncol, nz, 1, ncol, dt, lyr_surf, lyr_toa,
    f_cpair, f_rair, f_rho, f_z_mid,
    f_pk,          // pk (Exner function) is passed as the "exner" argument
    f_theta, f_qv, f_qc, f_qr, f_precl, f_relhum,
    f_temp_prev, f_temp_tend, ...);
Kessler_chost_physics_timestep_final(ncol, nz, f_cpair, f_temp, f_z_mid, f_phis, f_st_energy, ...);

params_helpers.transpose<f2c>(m_num_cols, nlevs);
params_computed.transpose<f2c>(m_num_cols, nlevs);
```

The `params_helpers` and `params_computed` structs in `kessler_functions.hpp`, the
`requested_buffer_size_in_bytes()` function, and `init_buffers()` are unchanged.

**`finalize_impl`:** added `Kessler_chost_physics_finalize` call (was previously a
no-op).

---

## Regenerating the bridge code

If the Kessler Fortran scheme or the EAMxx host variables change, regenerate with:

```bash
cd <repo-root>/xdsl-cpp

ccpp_xdsl \
    --suites       examples/kessler/scheme/kessler_suite.xml \
    --scheme-files examples/kessler/scheme/kessler.meta,examples/kessler/scheme/kessler_update.meta \
    --host-files   examples/kessler/host_eamxx/eamxx_kessler_host_mod.meta,examples/kessler/host_eamxx/eamxx_kessler_host_sub.meta \
    --bind-c \
    --directive acc \
    -o examples/kessler/bindc_eamxx_acc

cp -r examples/kessler/bindc_eamxx_acc/* \
      <repo-root>/E3SM/components/eamxx/src/physics/kessler/generated_bridge/
```

Then verify that the signatures in the regenerated `Kessler_ccpp_chost_cap.h` still
match the calls in `eamxx_kessler_process_interface.cpp`.

---

## Background and design notes

### Kessler scheme sources

- CAM-SIMA Fortran scheme: `atmospheric_physics/schemes/kessler/`
- Non-CAM Fortran version with C bindings: `components/eam/src/physics/crm/pam/external/physics/micro/kessler/kessler.f90`
- EAM implementation background: `components/eam/src/physics/crm/pam/external/physics/micro/kessler/Microphysics.h`

Additional initial condition files may need to be downloaded from
`https://web.lcrc.anl.gov/public/e3sm/inputdata/atm/scream/init/`.

### GPU data movement

EAMxx uses Kokkos Pack views (C row-major) internally.  The generated Fortran cap
expects column-major Fortran arrays.  The `params_helpers::transpose` and
`params_computed::transpose` methods in `kessler_functions.hpp` handle this
conversion before and after each Fortran call.

On GPU builds without OpenACC (`EAMXX_ENABLE_GPU && !EAMXX_ENABLE_OPENACC`), Fortran
runs on CPU and the bridge uses host mirror views (`h_*`).  On CPU or GPU+OpenACC
builds it uses Fortran-layout device views (`f_*`) directly.

The `Kessler_ccpp_cap.F90` generated with `--directive acc` contains
`!$acc enter/exit data` and `!$acc data` regions guarded by `#ifdef USE_GPU` that
manage GPU memory for the OpenACC path.

### ATMBufferManager

The Fortran-layout (`f_*`) and C++ Pack (`view_2d<Pack>`) scratch arrays are
allocated through EAMxx's `ATMBufferManager` in `requested_buffer_size_in_bytes()`
and `init_buffers()`.  The buffer layout counts are tracked as `static constexpr int`
members of `params_helpers` and `params_computed` in `kessler_functions.hpp`.
