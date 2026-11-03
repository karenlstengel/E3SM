# EAMxx Kessler Microphysics (Fortran bridge)

The Kessler (1969) warm-rain microphysics scheme as an EAMxx atmosphere
process (`kessler`). On this branch (`ren/kessler-bridge`) the physics is the
CCPP/CAM-SIMA Fortran from
[ESCOMP/atmospheric_physics](https://github.com/ESCOMP/atmospheric_physics)
(`schemes/kessler/kessler.F90` and `kessler_update.F90`), called from the
C++ process through a C++/Fortran bridge in `fortran_bridge/`. It runs on
the CPU or the host cores of a GPU node (the serial scheme), or on the GPU
(an OpenACC port of it). Both are compiled from checkouts outside this
repository, each with local performance changes that keep the results the
same (see [The Fortran Kessler versions](#the-fortran-kessler-versions)).

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
| `kessler_functions.hpp` | `KesslerMicrophysicsFunctions<Scalar, Device>` (KMF): view types and the `params_helpers` / `params_computed` structs that hold each field in C++ layout and in Fortran layout, with their `transpose()` |
| `fortran_bridge/kessler_eamxx_bridge.hpp/.cpp` | C++ side of the bridge: `kessler_eamxx_bridge_init/run/finalize`, which transpose the data and call the `bind(C)` Fortran routines |
| `fortran_bridge/kessler_eamxx_bridge_main.F90` | Fortran side for `kessler_init` / `kessler_run` (`bind(C)` wrappers) |
| `fortran_bridge/kessler_eamxx_bridge_update.F90` | Fortran side for `kessler_update_*`, plus applying the temperature tendency |
| `fortran_bridge/kessler_perf_log.F90` | per-routine wall-clock totals of the Fortran calls, written to a CSV at finalize |
| `fortran_bridge/ccpp_kinds.F90` | not compiled; the build uses `ccpp_kinds.F90` from atmospheric_physics (see [Building](#building)) |
| `driver_bridge_transpose.cpp`, `build_driver_bridge_transpose.sh` | standalone check and timing of the C++ <-> Fortran layout copies against their previous version (`driver_bridge_transpose_ref.hpp`); not built by CMake |
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
3. `F90_run`: `kessler_eamxx_bridge_run`, which
   - transposes into Fortran layout the fields the Fortran reads
     (`F90_run::transpose_c2f`): `rho`, `pk`, `z_mid`, `phis`, `theta`,
     `qv`, `qc`, `qr`, `T_mid`,
   - on a GPU build without OpenACC, copies those to the host
     (`sync_to_host`, inside the c2f timer),
   - calls `kessler_run`, then `kessler_update_timestep_init`,
     `kessler_update_run`, `T = T_prev + T_tend * dt`, and
     `kessler_update_timestep_final`, all in Fortran,
   - copies back (`sync_to_dev`) and transposes back
     (`F90_run::transpose_f2c`) the fields the Fortran writes: `theta`,
     `qv`, `qc`, `qr`, `T_mid`, `relhum`, `precl`, `T_mid_prev`,
     `T_mid_tend`, `st_energy`. Because of step 1 this writes the results
     straight into the fields; no copy back is needed.

HOMME takes its forcing from the change in `T_mid` across physics, so the
temperature update in step 3 is what carries the latent heating to the
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

## The Fortran Kessler versions

The bridge compiles one of two versions of `schemes/kessler/kessler.F90`
(see [External Fortran sources](#external-fortran-sources)). Both have
local changes (2026-11-01, uncommitted in their checkouts) that change
speed only; the versions before them are kept in each checkout's
`schemes/kessler/reference/`, and
`GPU_ports/atmospheric_physics/schemes/kessler/bench/` holds a standalone
benchmark that builds and compares all four (README and results there).

**Serial version** (CPU and gpucpu builds): the scheme's own structure, a
loop over columns with each column's sub-cycling `do while` inside, so
every column does only its own sub-steps. Change: the four powers of rain
water (`(qr*r)**0.1364` in both fall-speed computations, `qr**0.875` in
autoconversion and collection, `(r*qr)**0.2046` and `**0.525` in
evaporation) are skipped where `qr = 0`. `pow(0, e > 0) = 0`, so the
results are the same bit for bit; most levels are dry even in raining
columns. The same change as the C++ branch's item I2.

**OpenACC version** (GPU build). The port this branch started from ran
the sub-cycling as a host loop: every sub-step launched about ten kernels
over all columns, with finished columns masked, plus a reduction read back
on the host to test convergence, until the slowest column was done, so
every column did the slowest column's number of sub-steps (up to ~140 with
rain). It also created and deleted its device work arrays and printed the
iteration count on every call. The C++ branch had the same structure in its
first port. The version now in use (v7 of the tuning steps below):

- **Sub-cycling: one gang per column** (the C++ branch's item D1 in
  OpenACC), its levels spread over the gang's 64 vector lanes; each column
  runs its own sub-steps and stops when it reaches `dt`. The kernel runs at
  most 32 sub-steps per column per launch and is launched again while any
  column has not reached `dt` (one reduction read back per launch: ~5
  launches per call with rain). nvfortran spreads the level loops over the
  vector lanes only inside a counted loop without an `exit` (inside a
  `do while`, or a loop that exits, it runs them one level at a time;
  checked with `-Minfo=accel`), so a column that is done just skips its
  remaining sub-steps of a launch.
- **The column copied into shared memory for its sub-steps.** In the
  `(ncol, nz)` Fortran arrays a column's levels are `ncol` elements apart,
  so every access a warp makes over levels touches 32 cache lines. Each
  running column is copied into gang-private arrays at the start of a
  launch (nvfortran places them in CUDA shared memory; up to
  `nz_max = 128` levels), its sub-steps work there, and the updated
  `theta`, `qv`, `qc`, `qr`, `velqr` are written back at the end.
- **Columns on the vector lanes where no level neighbour is needed:**
  the point-wise initialisation and the relative humidity
  (`collapse(2)`, column index innermost, so a warp reads consecutive
  columns), the first CFL step and the precipitation rate (one thread per
  column).
- **First sub-step one thread per column** when there are at least
  `lanes_first_step_min_cols = 20000` columns (ne60: 48601 per GPU rank):
  every column needs at least one sub-step, and most need only that one
  (all of them without rain). The thread does the levels in order from the
  surface up; each level's sedimentation reads the level above before that
  level is updated, so it is the same computation. Only the columns still
  running go on to the gang-per-column kernel, which is not launched at all
  without rain. With fewer columns there are too few threads to keep the
  GPU busy this way, and the gang kernel does all sub-steps (measured
  crossover ~12000 columns without rain, ~24000 with rain).
- Work arrays and the per-column sub-cycling state stay on the device
  between calls, the print is gone, the constant `f5` and `xk` are scalars,
  `pow` is skipped for dry levels as in the serial version, and a column
  that cannot reach `dt` (`dt0` reaching 0, or 100000 sub-steps) returns an
  error instead of looping forever.

### Performance report (2026-11-01)

Standalone, with `GPU_ports/atmospheric_physics/schemes/kessler/bench/`
(synthetic columns as in the C++ and JAX drivers: 58 levels, top-down, a
given fraction of rainy columns; one A100 or one CPU core; one rank's
columns). Every step was kept only if it stayed bit for bit with the
OpenACC version this branch started from, at ne16 and ne60, rain fractions
0-12%, over 3-step chains, and was faster beyond the ~10% run-to-run spread.

OpenACC `kessler_run`, s/call, every call from the start state (min of 10):

| step | change | ne16 (3457), no rain | ne16, 8% rain | ne60 (48601), no rain | ne60, 8% rain |
|------|--------|------:|------:|------:|------:|
| start | host loop over sub-steps | 0.341e-3 | 24.5e-3 | 3.29e-3 | 266e-3 |
| v1 | one gang per column, sub-cycling inside it | 0.295e-3 | 2.20e-3 | 3.02e-3 | 26.1e-3 |
| v2 | + column in shared memory for the sub-steps | 0.283e-3 | 1.72e-3 | 3.59e-3 | 7.70e-3 |
| v3 | + columns on the lanes for initialisation and relhum | 0.267e-3 | 1.71e-3 | 2.28e-3 | ~7.4e-3 |
| v7 | + first sub-step one thread per column (>= 20000 columns) | **0.267e-3** | **1.73e-3** | **0.848e-3** | **6.46e-3** |
| | final vs start | 1.3x | 14x | 3.9x | 41x |

Tried and not kept (`bench/versions/README.md`): moving the error checks
into the sub-cycling kernel's reduction (v4; 5-10% slower with rain at
ne16) and counting the running columns in a separate thread-per-column
kernel (v5; +3-7% at ne60, -2-3% at ne16). Where the time goes now (nsys,
ne60, per call): without rain 0.73 ms of kernels (thread-per-column first
sub-step 0.41 ms, initialisation 0.21 ms, relhum 0.09 ms); with 8% rain
7.2 ms (gang-per-column sub-cycling 6.0 ms in 5 launches, its reductions
0.44 ms). For comparison, the C++ team-per-column `kessler_run` on the same
columns (before its I2 change, which gave another 1.03-1.14x on GPU):
1.15e-3 (ne16, 8% rain) and 7.56e-3 (ne60, 8% rain).

Serial version, one core: skipping `pow` for dry levels, 1.8x without rain
and 1.5-1.7x with 3-12% rain (ne16: 108 columns, ne60: 1519), bit for bit.

The bridge's C++ <-> Fortran layout copies: see
[Working with a Fortran bridge](#working-with-a-fortran-bridge-in-an-atmosphere-process),
item 3 (2-3.6x, bit for bit).

Results: the new serial version is bit for bit with the old one, and every
OpenACC step is bit for bit with the OpenACC version this branch started
from and repeats itself bit for bit. The OpenACC versions differ from the
serial one by at most 1.7e-9 relative (GPU `pow`/`exp`). (The starting
OpenACC version kept applying two unmasked updates to finished columns with
`dt0 = max(dt - elapsed, 0)`; that `dt0` is exactly 0 in these tests, so
they change nothing.) Model runs with the serial version and OpenACC v7:
see [Status and open items](#status-and-open-items).

Building the benchmark's serial versions on a GPU node needs `-noacc`: the
node's environment (`CRAY_ACCEL_TARGET=nvidia80`) makes nvfortran enable
OpenACC by default. EAMxx's CPU builds run on CPU nodes and are not
affected.

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

They copy only what is needed (2026-11-01): c2f only the arrays the
Fortran reads, f2c only the ones it writes (`params_helpers` holds inputs
only, so it has no f2c). Arrays the Fortran only writes (`relhum`,
`precl`, `T_mid_prev`, `T_mid_tend`, `st_energy`) are not copied in: the
Fortran writes every element. Nothing is zero-filled first (there used to
be an `init()` kernel for that every call). Each direction is one 2-D
`MDRangePolicy` (no integer division per element), with the index the
written array is contiguous in varying fastest (columns for c2f, levels
for f2c), so writes are contiguous and only the reads are strided. Results
are unchanged, bit for bit (`driver_bridge_transpose.cpp`). Per call, old
(zero-fill + c2f + f2c) vs new, one rank's columns:

| build | columns | old | new | speedup |
|-------|--------:|----:|----:|--------:|
| CPU (Serial, one core) | 108 | 0.41 ms | 0.13 ms | 3.1x |
| CPU (Serial, one core) | 1519 | 6.2 ms | 3.1 ms | 2.0x |
| GPU with OpenACC | 3457 | 0.29 ms | 0.09 ms | 3.0x |
| GPU with OpenACC | 48601 | 4.5 ms | 1.2 ms | 3.6x |
| GPU without OpenACC (incl. host copies) | 3457 | 3.9 ms | 2.0 ms | 1.9x |
| GPU without OpenACC (incl. host copies) | 48601 | 46 ms | 30 ms | 1.5x |

In the model the CPU copies cost more than in the driver (ne60: about
13 ms per call per direction, 1.26 s/day), probably because 128 ranks
share the node's memory bandwidth.

**4. Memory.** The Fortran-layout views, and C++ helpers that are not
fields (`z_mid`, `z_int`), are carved from the process's `ATMBufferManager`
allocation, as P3, SHOC and ZM do:

- `requested_buffer_size_in_bytes()` sums the sizes, using the
  `num_1d_*` / `num_2d_*` counts declared in each struct.
- `init_buffers()` hands consecutive pieces of the buffer to the views, in
  a fixed order.

The order of the `f_` views matters: the ones the Fortran reads
(`f_rho` .. `f_temp`) and the ones it writes (`f_theta` .. `f_st_energy`)
are each one contiguous range (checked in `init_buffers()`), so the GPU
build without OpenACC copies each direction to or from the host in one
piece. `initialize_impl` carves the host copies from one host block in the
same order.

To add a field: add the C++ and `f_` members to the struct, raise its
counts, add it to the pointer lists in `init_buffers()` (keeping the read
and write ranges contiguous, and updating the range lengths in
`kessler_eamxx_bridge_run` and the host carving in `initialize_impl`), and
to the directions of `transpose()` it needs. The C++ views of fields the process
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
| GPU without OpenACC (`nvidiagpu`, `OPENACC_GPU_OFFLOAD=FALSE`; "gpucpu") | CUDA | host, one core per rank, standard source | host copies `h_*` (one host block, carved in `initialize_impl`), copied every call: one copy per direction for the 2-D arrays plus one for `phis` / `precl` |

The `#if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)`
blocks in `kessler_functions.hpp`, `kessler_eamxx_bridge.cpp` and
`initialize_impl` select the host-mirror path. With OpenACC, the Fortran
must treat the incoming arrays as device memory (`deviceptr`) and keep its
own scratch on the device: the OpenACC `kessler_run` creates its work
arrays with `!$acc enter data create` on the first call and keeps them
(again only if `ncol` or `nz` change), rather than creating and deleting
them every call.

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

Both checkouts carry uncommitted local changes in `schemes/kessler/`: the
interface changes this branch needs (`cpair`/`rair` passed to
`kessler_init` as scalars instead of per-point arrays to `kessler_run`)
and the performance changes in
[The Fortran Kessler versions](#the-fortran-kessler-versions), with the
earlier versions in `schemes/kessler/reference/`. `build.sh` records only
this repository's commit and diff, so a build's provenance does not show
which state of these checkouts it compiled; note it with the runs, or
commit those checkouts first.

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

- Built 2026-11-01 from this branch at `4b2cebd3f7` (the pre-processing
  kernel and fenced timer of `025012875a`, and the layout copies above)
  with the external Fortran of that time:
  - `f90_cpu_4b2cebd3f7` and `f90_gpucpu_4b2cebd3f7`: the serial version
    described above.
  - `f90_gpu_4b2cebd3f7-dirtycecd5348`: OpenACC v7. Its OpenACC kernels sit
    at the lines of the current v7 source (collapse(2) initialisation,
    thread-per-column first CFL step and first sub-step, the
    `vector_length(64)` gang-per-column loop, precl, relhum). The repository
    diff of this build is only this README.
  - `f90_gpu_4b2cebd3f7` (no suffix) compiled OpenACC v1, before the v7
    source was in place. Do not use it; it has no run cases.
- Model runs (`GK_perf_utils/data/run_manifest.csv`): 24 cases on each of
  the three builds above (72 in all), created 2026-11-02, all
  `case.run success`. They have not yet been compared with the 2026-10-16
  Fortran runs. All the changes keep the results the same, so they should
  be bit for bit in `homme_atm.log`.
- The Fortran `errflg` is not checked after the calls.
- The changes in the external Fortran checkouts are uncommitted. Changing
  only those (not this repository) leaves the build tag the same, so
  `build.sh` then needs `--force` to rebuild.
- The atmospheric_physics paths in `CMakeLists.txt` are hard-coded to one
  user's scratch.
