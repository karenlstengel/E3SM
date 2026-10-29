# EAMxx Kessler Microphysics (Fortran bridge)

The Kessler (1969) warm-rain microphysics scheme as an EAMxx atmosphere
process (`kessler`). On this branch (`ren/kessler-bridge`) the physics is the
original CCPP/CAM-SIMA Fortran from
[ESCOMP/atmospheric_physics](https://github.com/ESCOMP/atmospheric_physics)
(`schemes/kessler/kessler.F90` and `kessler_update.F90`), called from the
C++ process through a C++/Fortran bridge in `fortran_bridge/`. The same
Fortran can run on the CPU, on the host cores of a GPU node, or on the GPU
through OpenACC.

The scheme has three moisture categories: water vapor (`qv`), cloud water
(`qc`) and rain water (`qr`, which falls relative to the air). It does
autoconversion, collection, evaporation and saturation adjustment, plus rain
sedimentation sub-cycled to satisfy the CFL condition. References: Kessler
(1969); Soong & Ogura (1973); Klemp & Wilhelmson (1978), J. Atmos. Sci. 35,
1070-1096; Durran & Klemp (1983), Mon. Wea. Rev. 111, 2341-2361.

Other implementations of the same scheme live on other branches of this
fork, so they can be compared against this one:

| branch               | implementation                                         |
|----------------------|--------------------------------------------------------|
| `ren/kessler-bridge` | this directory: the Fortran scheme through a C++/Fortran bridge |
| `ren/kessler-claude` | a C++/Kokkos port (its README and `performance_suggestions.md` hold the performance history for all three) |
| `ren/kessler-JAX`    | a Python/JAX version called from EAMxx                 |

## Files

| file | contents |
|------|----------|
| `eamxx_kessler_process_interface.hpp/.cpp` | `KesslerMicrophysics`, the EAMxx `AtmosphereProcess`: field requests, buffer setup, pre-processing, the bridge call |
| `kessler_functions.hpp` | `KesslerMicrophysicsFunctions<Scalar, Device>` (KMF): view types and the `params_helpers` / `params_computed` structs that hold each field in C++ layout and in Fortran layout, with their `init()` and `transpose()` |
| `fortran_bridge/kessler_eamxx_bridge.hpp/.cpp` | C++ side of the bridge: `kessler_eamxx_bridge_init/run/finalize`, which transpose the data and call the `bind(C)` Fortran routines |
| `fortran_bridge/kessler_eamxx_bridge_main.F90` | Fortran side for `kessler_init` / `kessler_run` (`bind(C)` wrappers) |
| `fortran_bridge/kessler_eamxx_bridge_update.F90` | Fortran side for `kessler_update_*`, plus applying the temperature tendency |
| `fortran_bridge/kessler_perf_log.F90` | per-routine wall-clock totals of the Fortran calls, written to a CSV at finalize |
| `fortran_bridge/ccpp_kinds.F90` | not compiled; the build uses `ccpp_kinds.F90` from atmospheric_physics (see [Building](#building)) |
| `CMakeLists.txt` | the `kessler` library: bridge sources, the external Fortran sources, include paths; defines `EAMXX_HAS_KESSLER` |
| `bridge/`, `kessler_jax/` | only stale `__pycache__` directories |

## The `kessler` atmosphere process

`KesslerMicrophysics` runs on the `physics` grid (the GLL dynamics grid for
this compset) and is registered in `physics/register_physics.hpp` when
`EAMXX_HAS_KESSLER` is defined. Levels are top-down, as everywhere in EAMxx:
C++ index 0 is the model top; the Fortran gets `lyr_toa = 1` and
`lyr_surf = nlev` (1-based).

Fields:

| field | role | units | notes |
|-------|------|-------|-------|
| `T_mid` | updated | K | overwritten with the Fortran's updated temperature |
| `p_mid`, `pseudo_density` | required | Pa | |
| `qv`, `qc`, `qr`, `qi` (tracers) | updated | kg/kg | mixing ratios wrt dry air; `qi` is only checked, Kessler has no ice |
| `phis` | updated | m2/s2 | surface geopotential (0 for this compset) |
| `rho`, `dz`, `pk`, `theta` | computed | kg/m3, m, -, K | derived in pre-processing |
| `precl` | computed | m/s | surface rain rate |
| `relhum` | computed | % | |
| `T_mid_prev`, `T_mid_tend`, `st_energy` | computed | K, K/s, J/kg | `kessler_update` bookkeeping; `st_energy` is dry static energy |
| `vapor_flux`, `water_flux`, `ice_flux`, `heat_flux` | computed | | only when the energy fixer is on; set to 0 |

`z_mid` and `z_int` are not fields on this branch: they live in the
process's buffer (`params_helpers.z_mid/z_int`) and are recomputed every
call.

What `run_impl` does each call (timers under `EAMxx::kessler::run::`):

1. Point the C++ views in `params_helpers` / `params_computed` at the field
   views (`rho`, `dz`, `pk`, `phis`, `theta`, `qv`, `qc`, `qr`, `precl`,
   `relhum`, `T_mid`, `T_mid_prev`, `T_mid_tend`, `st_energy`), so the
   bridge reads and writes the fields directly.
2. `preprocess`: one kernel, one team per column: `pk` = Exner function,
   `theta = T_mid / pk`, `dz` from `pseudo_density`,
   `rho = pseudo_density / (g * dz)`, then `z_int` (a team scan up from the
   surface) and `z_mid`. Same code as the C++ branch (its
   `performance_suggestions.md`, item I4). The timer is fenced, so on GPU
   it measures the kernel, not only its launch.
3. Zero the Fortran-layout holders (`params_*.init()`).
4. `F90_run`: `kessler_eamxx_bridge_run`, which
   - transposes every field into Fortran layout (`F90_run::transpose_c2f`),
   - on a GPU build without OpenACC, copies them to host mirrors
     (`sync_to_host`),
   - calls `kessler_run`, then `kessler_update_timestep_init`,
     `kessler_update_run`, `T = T_prev + T_tend * dt`, and
     `kessler_update_timestep_final`, all in Fortran,
   - copies back (`sync_to_dev`) and transposes back
     (`F90_run::transpose_f2c`). Because of step 1 this writes the results
     straight into the fields; no copy back is needed.

HOMME takes its forcing from the change in `T_mid` across physics, so the
temperature update in step 4 is what carries the latent heating to the
dynamics. `kessler_update_run` itself only accumulates `T_tend`; in
CAM-SIMA, `apply_tendency_of_air_temperature` applies it, so the bridge
does that step itself (`kessler_eamxx_bridge_update.F90`).

Range checks registered in `initialize_impl` run as pre/postconditions on
every call: `T_mid` in [100, 500] K and `qv` in [1e-13, 0.2] (invariants);
`qc`, `qr`, `qi` in [0, 0.1] and `precl` >= 0 (postconditions). They can be
turned off with `kessler::enable_precondition_checks` /
`enable_postcondition_checks`.

Constants come from `scream::physics::Constants` and are passed to the
Fortran `init` routines once: `latvap`, `P0` (as `pref`), `RHOW` (as
`rhoqr`), `Cpair`, `Rair`, `gravit`.

## Fortran entry points

From atmospheric_physics, called through the `_c` wrappers in
`fortran_bridge/`:

| routine | arguments (in / inout / out) |
|---------|------------------------------|
| `kessler_init` | in: `lv, pref, rhoqr, cpair, rair`; out: `errmsg, errflg` |
| `kessler_run` | in: `ncol, nz, dt, lyr_surf, lyr_toa, rho, z (= z_mid), pk`; inout: `theta, qv, qc, qr`; out: `precl, relhum, scheme_name, errmsg, errflg` |
| `kessler_update_init` | in: `gravit, cpair`; out: `errmsg, errflg` |
| `kessler_update_timestep_init` | in: `temp`; out: `temp_prev, ttend_t, errmsg, errflg` |
| `kessler_update_run` | in: `nz, ncol, dt, theta, exner, temp_prev`; inout: `ttend_t`; out: `errmsg, errflg` |
| `kessler_update_timestep_final` | in: `nz, temp, zm, phis`; out: `st_energy, errflg, errmsg` |

The OpenACC versions (see below) also take `ncol, nz` in
`kessler_update_timestep_init` and `ncol` in
`kessler_update_timestep_final`; `kessler_eamxx_bridge_update.F90` calls
the matching signature for each build. The Fortran `z` argument of
`kessler_run` is the height of each level (`z_mid`), not a layer
thickness.

## Working with a Fortran bridge in an atmosphere process

How this process connects C++ fields to the Fortran scheme, as a pattern
for bridging other Fortran schemes (the structs and transposes follow ZM's
in `physics/zm`).

**1. Fortran side (`fortran_bridge/*.F90`).** One module per group of
routines; each public routine is a thin `bind(C, name="..._c")` wrapper
around the scheme's own routine:

- Scalars are passed by value (`integer(kind=c_int), value`,
  `real(kind=c_real), value`); arrays as explicit-shape
  `real(kind=c_real), dimension(pcols, pver)`, column index first.
- `c_real` is `c_double` or `c_float` from `SCREAM_DOUBLE_PRECISION`
  (`#include "eamxx_config.f"`). The scheme's `kind_phys` comes from
  atmospheric_physics' `test/include/ccpp_kinds.F90`, which is `REAL64`, so
  EAMxx must be built in double precision (the default).
- `pcols`/`pver` are module variables, set in the init wrapper.
- `errmsg`/`errflg` are module variables; the wrappers do not currently
  check `errflg` after a call.

**2. C++ side (`fortran_bridge/kessler_eamxx_bridge.cpp`).** Declare the
`_c` routines in `extern "C"` with matching types (`Int` for `c_int`,
`Real` for `c_real`, `Real*` for arrays) and wrap them in C++ functions the
process calls: `kessler_eamxx_bridge_init` (from `initialize_impl`),
`kessler_eamxx_bridge_run` (from `run_impl`) and
`kessler_eamxx_bridge_finalize` (from `finalize_impl`).

**3. Layout.** EAMxx fields are `view_2d<Pack>` in `LayoutRight`: column
first, levels packed (`Pack::n` levels per element), so one column's levels
are contiguous. The Fortran wants `ncol x nlev` `Real` arrays in
column-major order (column index fastest). Each field therefore has a
second, Fortran-layout view (`f_*`, `LayoutLeft`, `Real`), and
`transpose<c2f>()` / `transpose<f2c>()` in `kessler_functions.hpp` copy
between them, reading single levels as `x(col, k / Pack::n)[k % Pack::n]`.
Tracers (`qv`, `qc`, `qr`) are strided slices of EAMxx's bundled tracer
field; indexing through the packed view handles the stride, so never pass
`.data()` of a tracer view to Fortran directly.

**4. Memory.** The Fortran-layout views, and C++ helpers that are not
fields (`z_mid`, `z_int`), are carved from the process's `ATMBufferManager`
allocation, as P3, SHOC and ZM do:

- `requested_buffer_size_in_bytes()` sums the sizes, using the
  `num_1d_*` / `num_2d_*` counts declared in each struct.
- `init_buffers()` hands consecutive pieces of the buffer to the views, in
  a fixed order.

To add a field: add the C++ and `f_` members to the struct, raise its
counts, add it to the pointer lists in `init_buffers()`, to `init()`, and
to both directions of `transpose()`. The C++ views of fields the process
already owns (`rho`, `theta`, `qv`, ...) are also carved, but `run_impl`
re-points them at the field views every call, so only `z_mid`, `z_int` and
the `f_` views actually use their buffer memory.

**5. Where the Fortran runs.** Three configurations, chosen at build time
by `OPENACC_GPU_OFFLOAD`, which sets `EAMXX_ENABLE_OPENACC` and the
OpenACC Fortran flags in
`cime_config/machines/cmake_macros/derecho_nvidiagpu.cmake`:

| build | Kokkos (HOMME, pre-processing, transposes) | Fortran Kessler | data handed to Fortran |
|-------|--------------------------------------------|-----------------|------------------------|
| CPU (`nvidia`) | Serial, host | host, standard source | `f_` views (host memory) |
| GPU with OpenACC (`nvidiagpu`, `OPENACC_GPU_OFFLOAD=TRUE`) | CUDA | GPU, OpenACC source, `-acc -gpu=cc80` | `f_` views: device pointers, used with `deviceptr` in the OpenACC loops |
| GPU without OpenACC (`nvidiagpu`, `OPENACC_GPU_OFFLOAD=FALSE`; "gpucpu") | CUDA | host, one core per rank, standard source | host mirrors `h_*` (allocated in `initialize_impl`), copied with `deep_copy` every call |

The `#if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)`
blocks in `kessler_functions.hpp`, `kessler_eamxx_bridge.cpp` and
`initialize_impl` select the host-mirror path. With OpenACC, the Fortran
must treat the incoming arrays as device memory (`deviceptr`) and create
its own scratch on the device (`!$acc enter data create`).

**6. Build.** `CMakeLists.txt` puts the bridge `.cpp`/`.F90` files, the
scheme's `.F90` files and the CCPP framework modules they use into one
library, sets `Fortran_MODULE_DIRECTORY` and `LINKER_LANGUAGE Fortran`,
and links it into `eamxx_physics`.

**7. Timing.** GPTL timers cover the C++ steps (`preprocess`, `F90_run`,
the transposes and syncs). Inside the Fortran wrappers each scheme call is
timed with `system_clock` and recorded by `kessler_perf_log` in memory (no
I/O per call). At finalize each rank adds its totals to a CSV holding a
GPTL-style summary (count, total, max, min per routine), at
`$KESSLER_PERF_LOG_PATH` or `./kessler_perf_log.csv` in the run directory,
with the same schema as the JAX branch's `kessler_perf_log.py`.

## The Kessler compset

`F2000-SCREAMv1-KESSLER`
(`2000_SCREAM%KESSLER-noAero_SLND_SICE_SOCN_SROF_SGLC_SWAV`, in
`components/eamxx/cime_config/config_compsets.xml`): EAMxx with the HOMME
dycore and stub land, ice, ocean, river, glacier and wave components.
For `SCREAM%KESSLER`, `namelist_defaults_eamxx.xml` and
`config_component.xml` set:

- `physics` = `mac_aero_mic`, and `mac_aero_mic` = `kessler`: Kessler is the
  only physics. No radiation, turbulence or aerosols.
- `mac_aero_mic::number_of_subcycles`: Kessler calls per physics step. The
  default depends on the grid (e.g. 12 at ne30np4); the perf cases set 1
  or 4 with `atmchange`.
- No topography: `topography_filename` unset, `phis = 0`.
- Zero aerosol fields, and zero values for surface fields that EAMxx's
  surface exporter needs to exist (`precip_*_surf_mass`, `sfc_flux_*`).
- `SCREAM_CMAKE_OPTIONS = "SCREAM_NP 4 SCREAM_NUM_VERTICAL_LEV 58 SCREAM_NUM_TRACERS 10"`.

The user-mods directory
`components/eamxx/cime_config/testdefs/testmods_dirs/eamxx/L58-kessler`
sets 58 levels in `SCREAM_CMAKE_OPTIONS` and aquaplanet-style coupler
settings (fixed orbit with zero eccentricity/obliquity, constant zenith
angle 42.05 deg, zero albedos).

Initial conditions are DCMIP2016 moist baroclinic wave states
(`FKESSLER_<grid>[_dayNN].cam.i.moist_baroclinic_wave_dcmip2016.nc`), set
with `atmchange initial_conditions::filename=...`. The derecho perf cases use
day-10 (raining) and day-1 (no rain) states from
`/glade/derecho/scratch/kstengel/inputdata/atm/scream/init/`. Other EAMxx
initial-condition files can be downloaded from
https://web.lcrc.anl.gov/public/e3sm/inputdata/atm/scream/init/.
Analytic initial conditions (computing the baroclinic wave inside EAMxx,
as StormSPEED does) are not implemented.

## Building

### External Fortran sources

The Fortran scheme is not in this repository. `CMakeLists.txt` takes it
from a checkout of atmospheric_physics at a fixed path:

| build | `PATH_TO_LEGACY_CAM_SIMA` | source |
|-------|---------------------------|--------|
| OpenACC off (CPU, gpucpu) | `/glade/derecho/scratch/kstengel/E3SM/atmospheric_physics` | ESCOMP/atmospheric_physics, tag `atmos_phys0_20_000` |
| OpenACC on (GPU) | `/glade/derecho/scratch/kstengel/E3SM/GPU_ports/atmospheric_physics` | OpenACC port, `johnmauff/atmospheric_physics`, branch `openacc-fixes` |

From the selected checkout it compiles `schemes/kessler/kessler.F90`,
`schemes/kessler/kessler_update.F90`,
`schemes/utilities/geopotential_temp.F90`, `test/include/ccpp_kinds.F90`,
and three modules of the CCPP framework, which must be cloned into
`test/lib/ccpp-framework` first:

```bash
cd <atmospheric_physics>/test && mkdir -p lib && cd lib
git clone -b develop --depth 1 https://github.com/NCAR/ccpp-framework.git
```

To build elsewhere, change the two paths at the top of `CMakeLists.txt`.

### A case

On derecho, from `cime/scripts`:

```bash
./create_newcase --case <case> --compset F2000-SCREAMv1-KESSLER --res ne16_ne16 \
    --machine derecho --compiler nvidia --project <project> \
    --user-mods-dir <E3SM>/components/eamxx/cime_config/testdefs/testmods_dirs/eamxx/L58-kessler
cd <case>
./xmlchange DEBUG=FALSE,NTASKS=128,NTHRDS=1,ROOTPE=0
./case.setup
./xmlchange CAM_TARGET=theta-l_kokkos
./case.build
```

For GPU, use `--compiler nvidiagpu` and, before `case.setup`:

```bash
./xmlchange NTASKS=4,NGPUS_PER_NODE=4,GPU_TYPE=a100
./xmlchange KOKKOS_GPU_OFFLOAD=TRUE,OPENMP_GPU_OFFLOAD=FALSE,OVERSUBSCRIBE_GPU=FALSE
./xmlchange OPENACC_GPU_OFFLOAD=TRUE    # Kessler on the GPU; FALSE runs it on the host cores
```

What the derecho builds come out as:

| | CPU (`nvidia`) | GPU (`nvidiagpu`) |
|-|----------------|-------------------|
| Kokkos backend | Serial (one rank per core, `NTHRDS=1`) | CUDA, sm_80 (A100) |
| `SCREAM_PACK_SIZE` | 16 | 1 |
| ranks per node | 128 | 4 (one per GPU) |
| flags | `-O2 -Mnofma -Mflushz -Kieee` (C, C++, Fortran) | same, plus `--fmad=false` on CUDA kernels; with OpenACC, `-acc -gpu=cc80 -Mnofma` on Fortran |

The `-O2 -Mnofma -Mflushz -Kieee` flags come from
`components/eamxx/cmake/machine-files/derecho.cmake` (currently an
uncommitted change on this branch). They turn off fused multiply-add so the
C++, Fortran and JAX builds and the StormSPEED (CESM) builds can be
compared bit for bit.

Runtime settings used by the perf cases: `./xmlchange ATM_NCPL=48` (30 min
physics step), and with `./atmchange`: `ctl_nl::se_ftype=2`,
`transport_alg=12`, `statefreq=488`, `dt_remap_factor=2`,
`dt_tracer_factor=6`, `hypervis_subcycle=1`, `nu_top=2.5e5`, and per grid:

| grid | `se_tstep` | `se_nsplit` | `hypervis_subcycle_q` |
|------|-----------:|------------:|----------------------:|
| ne16np4 | 300 | 1 | 3 |
| ne30np4 | 300 | 1 | 6 |
| ne60np4 | 150 | 2 | 6 |

### Derecho perf workflow

The C++/Fortran/JAX comparison uses scripts outside this repository, in
`/glade/derecho/scratch/kstengel/GK_perf_utils/workflow/`:

- `build.sh f90 <hw>` builds one executable per hardware variant into
  `e3sm_test/builds/f90_<hw>_<commit>[-dirty<hash>]/`, from this checkout
  (it checks that the branch is `ren/kessler-bridge`). `hw` is `cpu`,
  `gpu` (OpenACC on), `gpucpu` (OpenACC off) or `gpucpumc`. `gpucpumc`
  compiles the OpenACC Fortran for the host cores (`-acc=multicore`)
  through a `KESSLER_ACC_MULTICORE` CMake option that exists only in the
  uncommitted patch of the `f90_gpucpumc_514f63089f-dirtyb61aa4b9` build,
  not in this tree.
- `run.py` creates run cases that point `EXEROOT` at those builds, sets the
  runtime settings above, and records each case in a manifest.

## Status and open items

- The pre-processing kernel and its fenced timer match the C++ branch
  (2026-10-29, uncommitted, not built yet). Next: build, check bit for bit
  against the current `f90_*_514f63089f-*` builds at ne16, then rerun the
  perf matrix.
- The Fortran `errflg` is not checked after the calls.
- Per call, the `params_*.init()` zero-fill before the transpose is an
  extra kernel that could be dropped (the transpose overwrites every
  value).
- The atmospheric_physics paths in `CMakeLists.txt` are hard-coded to one
  user's scratch.
