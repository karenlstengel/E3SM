# EAMxx Kessler Bridge

This directory contains the EAMxx C++ atmosphere process interface for the Kessler
microphysics scheme.  The Fortran bridge code that connects EAMxx to the Kessler
Fortran scheme is **automatically generated** by the
[xdsl-cpp](https://github.com/xdsl-project/xdsl-cpp) framework (CCPP dialect) rather
than written by hand.  This document describes the full generation process so it can
be reproduced or updated.

## Updated list of bugs from xdsl-ccpp bridge generation

1. the cap function calls have `kessler_suite_suite_` names instead of just `kessler_suite_`
2. didn't add in checks to switch between the OpenACC function calls and CPU function calls 
3. xdsl-ccpp missed that `kessler_update_timestep_final` expects the `errflg` before `errmsg` but all of the other kessler & kessler_update calls expect `errmsg` before `errflg`. All of the functions in the kessler & kessler_update have argument orders that match the corresponding meta file entry.

### TODO: fix bug #2 — OpenACC/CPU signature branch missing in `kessler_suite_cap.F90`

Investigated 2026-07-29. Not yet implemented.

**Diagnosis:**

- `kessler/CMakeLists.txt` links `GPU_ports/atmospheric_physics/schemes/kessler/kessler_update.F90`
  when `EAMXX_ENABLE_OPENACC` is ON, and the plain
  `atmospheric_physics/schemes/kessler/kessler_update.F90` otherwise. The two versions have
  **different argument lists** for two subroutines:

  | Subroutine | CPU version | GPU_ports version |
  |---|---|---|
  | `kessler_update_timestep_init` | `(temp, temp_prev, ttend_t, errmsg, errflg)` | `(ncol, nz, temp, temp_prev, ttend_t, errmsg, errflg)` |
  | `kessler_update_timestep_final` | `(nz, cpair, temp, zm, phis, st_energy, errflg, errmsg)` | `(nz, ncol, cpair, temp, zm, phis, st_energy, errflg, errmsg)` |

  (`kessler_run`/`kessler_init` are identical in both trees — no branch needed there.)
- The GPU_ports versions of these subroutines already contain the real `!$acc parallel loop
  collapse(2) deviceptr(...)` directives, which assume the incoming arrays are already device
  pointers (no host staging). Nothing needs to be added to the physics scheme itself.
- `Kessler_ccpp_cap.F90`'s `!$acc enter data copyin(...)` / `USE_GPU` scaffolding is a **different,
  unused architecture** (it assumes a host module `eamxx_kessler_host_mod` supplying host-resident
  module variables; that module was never generated and doesn't exist in this tree). It is not the
  pattern to follow here — do not port those directives.
- `kessler_suite_cap.F90` currently hardcodes only the GPU_ports (`ncol`-included) signature for
  both calls above, unconditionally. This compiles only when `EAMXX_ENABLE_OPENACC=TRUE`; a
  pure-CPU build fails to compile `kessler_suite_cap.F90` because the linked CPU subroutines don't
  take `ncol`.
- The handwritten reference bridge (`fortran_bridge/kessler_eamxx_bridge_update.F90:83-91`, no
  longer compiled but kept for reference) already solves this correctly with an
  `#if defined(EAMXX_ENABLE_GPU) && defined(EAMXX_ENABLE_OPENACC)` / `#else` branch selecting the
  right argument list per subroutine.
- The C++ interface's existing `h_*`/`f_*` selection in `eamxx_kessler_process_interface.cpp`
  (`#if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)` / `#else`) already correctly
  implements all three cases (GPU+OpenACC -> `f_*` device pointers directly; GPU without OpenACC ->
  `h_*` host mirrors; CPU-only -> `f_*` directly). No change needed there.

**To-do list:**

1. Expose `EAMXX_ENABLE_OPENACC` to Fortran compilation. In `kessler/CMakeLists.txt`, after
   `add_library(kessler ...)`, add:
   ```cmake
   if (EAMXX_ENABLE_OPENACC)
     target_compile_definitions(kessler PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:EAMXX_ENABLE_OPENACC>)
   endif()
   ```
   Key this off `EAMXX_ENABLE_OPENACC` alone (not combined with `EAMXX_ENABLE_GPU`) since that's the
   exact condition `kessler/CMakeLists.txt:5` uses to choose which `kessler_update.F90` gets linked.
2. Fix `generated_bridge/kessler_suite_cap.F90`:
   - In `kessler_suite_suite_timestep_initial`, branch the call to `kessler_update_timestep_init`
     on `#ifdef EAMXX_ENABLE_OPENACC` between the `(ncol, nz, ...)` and `(...)` (no `ncol`/`nz`)
     argument lists.
   - In `kessler_suite_suite_timestep_final`, branch the call to `kessler_update_timestep_final`
     the same way (with/without `ncol`).
3. Decide whether to apply the same edit to the raw-generator-output copy at
   `xdsl-ccpp-generated/bindc_eamxx/kessler_suite_cap.F90`, or rely on this README section to
   reapply the fix after future regeneration.
4. Decide the fate of `Kessler_ccpp_cap.F90` / `Kessler_ccpp_cap.h`: it's listed in
   `KESSLER_F90_SRCS` but `use`s `eamxx_kessler_host_mod`, which doesn't exist anywhere in the tree.
   Confirm whether it currently compiles; if not, either stub the missing module or remove
   `Kessler_ccpp_cap.F90`/`.h` from the build since nothing calls into it.
5. Build and verify both configurations:
   - `EAMXX_ENABLE_OPENACC=OFF`: confirm `kessler_suite_cap.F90` compiles against the plain
     `atmospheric_physics/schemes/kessler/kessler_update.F90`.
   - `EAMXX_ENABLE_OPENACC=ON` (derecho GPU): confirm it compiles against
     `GPU_ports/atmospheric_physics/schemes/kessler/kessler_update.F90` and that `-Minfo=accel`
     shows the existing `deviceptr` directives being picked up.
   - Run the standalone eamxx kessler test in both configs and compare `T_mid`, `qv`/`qc`/`qr`,
     `precl` output for correctness, not just successful compilation.

## How xdsl-ccpp works for a C++ host, and where it broke for EAMxx

Investigated 2026-07-29 by reading the `xdsl-ccpp` tool itself
(`multilanguage_plan.md`, `multilanguage_limitations.md`, `xdsl_ccpp/transforms/suite_cap.py`)
and diffing the `.meta`/`.F90` files actually used for generation against the two real Kessler
scheme trees (`atmospheric_physics/` and `GPU_ports/atmospheric_physics/`). This is the first time
xdsl-ccpp has generated a bridge for EAMxx as a C++ host model, so this section records what the
tool assumes, what it documents as unsupported, and the specific root causes behind bugs #1-#3
above. No code changes were made for this investigation.

### The chost-cap mechanism is the tool's real, intended path for C++ hosts

`language = c++` in a host `.meta` file's `[ccpp-table-properties]` block is a fully implemented,
deliberate feature (`multilanguage_plan.md`, "Open Design Questions" -> "`language = c++` for C++
host models [Implemented]"). Setting it auto-activates the chost-cap generation mode
(`Kessler_ccpp_chost_cap.F90` / `Kessler_chost.hpp`) with no extra CLI flags. The tool ships its own
toy reference example at `xdsl-cpp/examples/kessler/host_cpp/`, which `host_eamxx/` was modeled on.
So the overall approach taken here is correct and sanctioned by the tool, not a misuse of it.

### The tool's own docs already flag the GPU/OpenACC gap, and mark it unresolved

`multilanguage_limitations.md` section 2, "GPU Memory Management," states:

> The chost cap is a CPU BIND(C) wrapper. When physics schemes run on a GPU, the C++ host is
> responsible for ensuring arrays are in the correct device memory space before calling the cap.
> **The generated code provides no help with this.**
>
> Kokkos + Fortran OpenACC: Kokkos device allocations (`CudaSpace`) are invisible to the OpenACC
> runtime. The host must either use CUDA Unified Memory (`CudaUVMSpace`) or call `acc_map_data` to
> register already-placed device pointers with the OpenACC runtime before calling the cap.

The doc's own priority table still marks this open ("Blocks real use? Yes, for GPU builds --
Medium-High effort"). All of the tool's shipped examples are CPU-only demos, so the EAMxx Kessler
bridge is the first real exercise of a C++ host driving a chost-cap-generated scheme on GPU with
OpenACC. The GPU_ports `kessler_update.F90` avoids needing `acc_map_data`/UVM by using
`!$acc parallel loop ... deviceptr(...)` clauses, which tell the compiler to trust the incoming
pointer as an already-valid device address rather than going through OpenACC's present-table
bookkeeping -- this works because Kokkos CUDA-space pointers and NVHPC's OpenACC codegen share the
same GPU context, but it's a property of this specific scheme's directives, not something the
xdsl-ccpp-generated cap itself provides or enforces.

### Bug #1 root cause (`kessler_suite_suite_*` naming) -- self-inflicted via the suite XML

`xdsl-cpp/examples/kessler/scheme/kessler_suite.xml` names the suite `kessler_suite`. The generator
(`suite_cap.py`) builds function names as `<suite_name>_suite_<lifecycle>`, so `kessler_suite` +
`_suite_register` -> `kessler_suite_suite_register`. Not a generator defect -- naming the suite
`kessler` instead of `kessler_suite` and regenerating would produce clean `kessler_suite_register`
names.

### Bug #2 root cause (missing OpenACC/CPU signature switch) -- drifted metadata, not a generator defect

Diffing the three copies of `kessler_update.meta` in play:

- `xdsl-cpp/examples/kessler/scheme/kessler_update.meta` -- the copy that actually generated this
  bridge
- `atmospheric_physics/schemes/kessler/kessler_update.meta` -- the CPU scheme's own meta
- `GPU_ports/atmospheric_physics/schemes/kessler/kessler_update.meta` -- the GPU_ports scheme's own
  meta

found:

- The two "real" repo metas (CPU and GPU_ports) are **identical to each other** and both describe
  the CPU-style signature: no `ncol`/`nz` arguments on `kessler_update_timestep_init` /
  `kessler_update_timestep_final`, and no `memory_space = device` tags.
- But `GPU_ports/atmospheric_physics/schemes/kessler/kessler_update.F90`'s actual Fortran already
  has explicit `ncol`/`nz` arguments on those two subroutines (needed for its `!$acc deviceptr`
  clauses) and real `!$acc parallel loop` directives. **Its own committed `.meta` file was never
  updated to match its own Fortran source** -- a pre-existing drift bug in the GPU_ports fork,
  independent of xdsl-ccpp and predating this integration effort.
- The generation-source meta (the one that actually produced `kessler_suite_cap.F90`) is a third,
  hand-customized variant that adds the missing `ncol`/`nz` entries and `memory_space = device`
  annotations so generation would succeed against the GPU_ports signature.

CCPP metadata has no mechanism to express "this scheme's argument list depends on a build flag" --
one `.meta` file describes one signature. Feeding the generator a meta customized for the GPU_ports
signature necessarily produces a suite cap that only compiles against GPU_ports. This is a
consequence of the CPU/GPU_ports Fortran forks having genuinely different signatures with no
metadata mechanism to express both, compounded by the pre-existing GPU_ports metadata drift above --
not a generator bug.

### Bug #3 root cause (errflg/errmsg order) -- same hand-patched meta, likely a copy-paste slip

Across all four `kessler_update_*` subroutines, only `kessler_update_timestep_final` puts `errflg`
before `errmsg` in the real Fortran argument list -- both the CPU and GPU_ports trees agree on this
being the one exception. Both "real" `.meta` files correctly list `errflg` before `errmsg` for this
entry. But the generation-source meta has `errmsg` before `errflg` for `timestep_final` specifically,
while getting the other three subroutines' (already errmsg-first) order right. Consistent with
whoever added the missing `ncol` entry to this meta file applying the common errmsg-first convention
uniformly across all entries without checking that `timestep_final` is the outlier.

### Bottom line

The chost-cap mechanism itself works as designed and is the correct tool for an EAMxx-style C++
host. The actual failure points were: (1) a genuinely unaddressed, tool-documented gap around
GPU+OpenACC pointer interop that no prior use of this tool had exercised, worked around adequately
here by the scheme's own `deviceptr` clauses; and (2) metadata drift/hand-editing around the CPU vs.
GPU_ports argument-signature difference in `kessler_update_timestep_init`/`_final`, which produced a
generation-source meta that could only ever describe one of the two variants.

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
