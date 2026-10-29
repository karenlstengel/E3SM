# EAMxx Kessler Microphysics (C++/Kokkos)

The Kessler (1969) warm-rain microphysics scheme as an EAMxx atmosphere
process (`kessler`), implemented in C++ with Kokkos so the same source runs
on CPU and GPU. It is a port of the CCPP/CAM-SIMA Fortran scheme in
[ESCOMP/atmospheric_physics](https://github.com/ESCOMP/atmospheric_physics),
`schemes/kessler/kessler.F90` and `kessler_update.F90`.

The scheme has three moisture categories: water vapor (`qv`), cloud water
(`qc`, moves with the flow) and rain water (`qr`, falls relative to the
air). It does autoconversion, collection, evaporation and saturation
adjustment, plus rain sedimentation sub-cycled to satisfy the CFL condition.
References: Kessler (1969); Soong & Ogura (1973); Klemp & Wilhelmson (1978),
J. Atmos. Sci. 35, 1070-1096; Durran & Klemp (1983), Mon. Wea. Rev. 111,
2341-2361.

Other implementations of the same scheme live on other branches of this
fork, so they can be compared against this one:

| branch                | implementation                                                    |
|-----------------------|-------------------------------------------------------------------|
| `ren/kessler-claude`  | this directory: C++/Kokkos                                        |
| `ren/kessler-bridge`  | the Fortran scheme called from EAMxx through a C++/Fortran bridge (`fortran_bridge/`) |
| `ren/kessler-JAX`     | a Python/JAX version called from EAMxx                            |

Performance work, results and open items are tracked in
[performance_suggestions.md](performance_suggestions.md).

## Files

| file | contents |
|------|----------|
| `eamxx_kessler_process_interface.hpp/.cpp` | `KesslerMicrophysics`, the EAMxx `AtmosphereProcess`: field requests, buffer setup, pre/post-processing around `kessler_run` |
| `kessler_functions.hpp` | `KesslerMicrophysicsFunctions<Scalar, Device>` (KMF): types, the `Workspace` struct, function declarations |
| `kessler_functions_impl.hpp` | implementations: `kessler_run` and the three `kessler_update_*` functions |
| `eti/kessler.cpp` | explicit template instantiation for CPU builds and GPU builds with relocatable device code (see CMakeLists.txt) |
| `CMakeLists.txt` | the `kessler` library; defines `EAMXX_HAS_KESSLER`, which registers the process in `physics/register_physics.hpp` |
| `driver_kessler.cpp` | standalone test and benchmark driver for `kessler_run` (not built by CMake; see [Standalone driver](#standalone-driver)) |
| `driver_kessler_hostloop.hpp` | the previous, host-loop version of `kessler_run`, kept as a reference for the driver |
| `driver_kessler_variants.hpp` | experimental variants of `kessler_run` (skip-dry-pow, team scratch) for the driver to time and check |
| `driver_kessler_preprocess.hpp` | the old and new versions of `run_impl`'s pre-processing, for the driver's `--preprocess` test |
| `build_driver_kessler.sh` | builds the driver against an existing case build |
| `performance_suggestions.md` | performance history, results, translation guide, open items |
| `fortran_bridge/` | Fortran bridge sources. Not compiled on this branch; used on `ren/kessler-bridge` |

## The `kessler` atmosphere process

`KesslerMicrophysics` runs on the `physics` grid (the GLL dynamics grid for
this compset). Levels are top-down, as everywhere in EAMxx: index 0 is the
model top, `nlev-1` the surface.

Fields:

| field | role | units | notes |
|-------|------|-------|-------|
| `T_mid` | updated | K | latent heating applied through `T_mid_tend` |
| `p_mid`, `pseudo_density` | required | Pa | |
| `qv`, `qc`, `qr`, `qi` (tracers) | updated | kg/kg | mixing ratios wrt dry air; `qi` is only checked, Kessler has no ice |
| `phis` | updated | m2/s2 | surface geopotential (0 for this compset) |
| `rho`, `dz`, `pk`, `theta` | computed | kg/m3, m, -, K | derived in pre-processing |
| `z_mid`, `z_int` | computed | m | heights above the surface |
| `precl` | computed | m/s | surface rain rate, averaged over the step |
| `relhum` | computed | % | |
| `T_mid_prev`, `T_mid_tend`, `st_energy` | computed | K, K/s, J/kg | `kessler_update` bookkeeping; `st_energy` is dry static energy |
| `vapor_flux`, `water_flux`, `ice_flux`, `heat_flux` | computed | | only when the energy fixer is on; set to 0 |

What `run_impl` does each call, in order (each step has an EAMxx timer
under `EAMxx::kessler::run::`):

1. `preprocess`: `pk` = Exner function, `theta = T_mid / pk`, `dz` from
   `pseudo_density`, and `rho = pseudo_density / (g * dz)`; then `z_int`
   and `z_mid` from `dz` (a team-level scan from the surface). One kernel,
   one team per column (the scan policy, as in SHOC's preprocess), working
   on whole packs; padding lanes get benign inputs.
2. `kessler_update_timestep_init`: save `T_mid` in `T_mid_prev`, zero
   `T_mid_tend`.
3. `kessler_run`: the microphysics; updates `theta`, `qv`, `qc`, `qr`,
   and computes `precl` and `relhum`.
4. `kessler_update_run`: `T_mid_tend += (theta * pk - T_mid_prev) / dt`.
5. `kessler_apply_T_tend`: `T_mid = T_mid_prev + T_mid_tend * dt`. HOMME
   takes its forcing from the change in `T_mid` across physics, so without
   this step the latent heating would be lost.
6. `kessler_update_timestep_final`: `st_energy = cp * T_mid + g * z_mid + phis`.

Range checks registered in `initialize_impl` run as pre/postconditions on
every call: `T_mid` in [100, 500] K and `qv` in [1e-13, 0.2] (invariants);
`qc`, `qr`, `qi` in [0, 0.1] and `precl` >= 0 (postconditions). They can be
turned off with `kessler::enable_precondition_checks` /
`enable_postcondition_checks`; at ne16 on GPU they cost more than
`kessler_run` itself (performance_suggestions.md, section I).

Constants come from `scream::physics::Constants`: `pref = P0/100` (hPa),
`latvap`, `rhoqr = RHOW`, `gravit`, `Cpair`, `Rair`.

## C++/Kokkos implementation of `kessler_run`

`kessler_run(ncols, nz, dt, lyr_surf, lyr_toa, rhoqr, lv, pref, cpair, rair,
rho, z_mid, pk, theta, qv, qc, qr, precl, relhum, ws)` in
`kessler_functions_impl.hpp`.

**One kernel, one team per column.** Each column runs its own sub-cycle
loop and does only as many passes as its own CFL limit needs, like the
Fortran, whose `do while` sits inside `do col`. A dry column finishes in one
pass. Inside a team:

- Initialisation over all levels: `r = 0.001*rho`, `rhalf = sqrt(rho_surf/rho)`,
  `pc = 3.8 / (pk^(cp/R) * pref)`, `qr = max(qr, 0)`, terminal fall speed
  `velqr`. Then the first sub-cycle step `dt0 = min(dt, min over levels of
  0.8*dz/velqr)`.
- Each pass: accumulate surface precipitation; sedimentation for all levels
  (no flux into the top level); barrier; microphysics (autoconversion and
  collection, saturation adjustment, evaporation) fused with the new fall
  speed and its CFL limit; advance time; recompute `dt0`; barrier. The CFL
  limit is re-checked every pass, as in the Fortran.
- After the loop: relative humidity, `precl = accumulated / dt`.

Implementation details:

- Level loops are `TeamVectorRange` over packs (`ekat::Pack`, size
  `SCREAM_PACK_SIZE`); math is Pack-aware (`ekat::pow/exp/sqrt/max/min`).
  On the CPU build nvc++ turns `ekat::pow`/`exp` on a 16-wide pack into
  4 calls to its AVX2 vector routines (`__pd_pow_4`, `__pd_exp_4`).
- The four powers of rain water (`qr` in autoconversion/collection, `r*qr`
  in evaporation and in the fall speed) go through `detail::pow_nonneg`,
  which returns 0 without calling `pow` when every lane of the pack is 0.
  `pow(0, x > 0) = 0`, so results are identical; it skips the cost for
  dry levels, which are most of the column even where it rains.
- The vertical neighbour for sedimentation and the CFL limit, which may sit
  in the next pack, comes from `ekat::shift_left/shift_right`.
- Team size from `ekat::TeamPolicyFactory::get_default_team_policy(ncols,
  nlev_packs)`, as in P3/SHOC: 1 thread per team on CPU; on GPU nlev_packs
  rounded up to a multiple of 32, capped at 128 (64 for 58 levels).
- The sub-cycle state (`dt0`, elapsed time, accumulated precipitation, pass
  count, error code) is team-local; every thread keeps an identical copy.
- Padding lanes (when `nz` is not a multiple of the pack size) and the
  excluded boundary level are given safe inputs before any division or pow,
  so floating-point-exception builds do not trap.
- `qr` (and the other tracers) are slices of EAMxx's bundled tracer field,
  so they are strided. They are never `scalarize`d; single levels are read
  as `qr(col, k / Pack::n)[k % Pack::n]`. Only `rho`, `z_mid` and `velqr`
  are scalarized, and their contiguity is checked.
- Errors: a column whose first `dt0` is below 1e-12 ("bad time splitting",
  as in the Fortran), or that cannot finish (`dt0 <= 0`, or 100000 passes),
  stores an error code; after the kernel one reduction on the host throws
  the matching message.
- `lyr_surf`/`lyr_toa` are 0-based; the process passes `nlev-1` and 0.

**Workspace.** Scratch arrays are passed in as `KMF::Workspace`: per level
`r`, `rhalf`, `velqr`, `sed`, `pc`, `f5` (`f5` is unused), and per column
`dt0`, `mask`, `time_counter`, `precl_acc` (written at the end: the final
step, the error code, the elapsed time and the accumulated precipitation).
The process carves them from one `ATMBufferManager` allocation in
`requested_buffer_size_in_bytes()` / `init_buffers()`, like P3 and SHOC.

**History.** The first port ran the sub-cycle loop on the host, sweeping all
columns with about nine kernels per pass until the slowest column was done.
That did `max(passes) * ncols` work instead of `sum(passes)`, which made the
CPU build 3-5.5x slower than Fortran with rain. The current version gives
bit-identical answers; in model runs with rain, `kessler::run` on the
slowest rank is 3-6x faster on CPU and 2.5-14x faster on GPU
(performance_suggestions.md, sections B-D and G-H). Skipping `pow` on dry
packs and the packed pre-processing then gave another 1.15-1.33x on CPU
and 1.03-1.14x on GPU, for 3.7-7.1x (CPU) and 2.6-15x (GPU) over the
first port, and C++ is faster than Fortran in every perf case
(performance_suggestions.md, section J). The old version is kept in
`driver_kessler_hostloop.hpp`.

`kessler_update_timestep_init`, `kessler_update_run` and
`kessler_update_timestep_final` are the `kessler_update.F90` routines
(steps 2, 4 and 6 above).

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

## Building and running

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
./xmlchange KOKKOS_GPU_OFFLOAD=TRUE,OPENACC_GPU_OFFLOAD=FALSE,OPENMP_GPU_OFFLOAD=FALSE,OVERSUBSCRIBE_GPU=FALSE
```

What the derecho builds come out as:

| | CPU (`nvidia`) | GPU (`nvidiagpu`) |
|-|----------------|-------------------|
| Kokkos backend | Serial (one rank per core, `NTHRDS=1`) | CUDA, sm_80 (A100) |
| `SCREAM_PACK_SIZE` | 16 (58 levels = 4 packs) | 1 |
| ranks per node | 128 | 4 (one per GPU) |
| C++ flags | `-O2 -Mnofma -Mflushz -Kieee` | same, plus `--fmad=false` on CUDA kernels |

The flags come from `components/eamxx/cmake/machine-files/derecho.cmake`.
They turn off fused multiply-add, so C++, Fortran and JAX builds can be
compared bit for bit. On CPU and RDC GPU builds, `kessler_functions_impl.hpp`
is compiled once in `eti/kessler.cpp`; on other GPU builds
`kessler_functions.hpp` includes it in every file that uses it.

Runtime settings used by the perf cases: `./xmlchange ATM_NCPL=48` (30 min
physics step), and with `./atmchange`: `ctl_nl::se_ftype=2`, `transport_alg=12`, `statefreq=488`,
`dt_remap_factor=2`, `dt_tracer_factor=6`, `hypervis_subcycle=1`,
`nu_top=2.5e5`, and per grid:

| grid | `se_tstep` | `se_nsplit` | `hypervis_subcycle_q` |
|------|-----------:|------------:|----------------------:|
| ne16np4 | 300 | 1 | 3 |
| ne30np4 | 300 | 1 | 6 |
| ne60np4 | 150 | 2 | 6 |

### Derecho perf workflow

The C++/Fortran/JAX performance comparison uses scripts outside this
repository, in `/glade/derecho/scratch/kstengel/GK_perf_utils/workflow/`:

- `build.sh <impl> <hw>` builds one executable per implementation and
  hardware into `e3sm_test/builds/<impl>_<hw>_<commit>[-dirty<hash>]/`.
  Builds are tagged by commit and uncommitted diff, so older run cases keep
  their executable, and `builds/<impl>_<hw>` points at the newest.
- `run.py` creates run cases that point `EXEROOT` at those builds (no
  compiling), sets the runtime settings above and records each case in a
  manifest.
- `interactive/` holds scripts for building and running small cases inside
  an interactive GPU job, and for comparing old and new builds.

How to use them, and the results, are in performance_suggestions.md
(sections F-H).

### Standalone driver

`driver_kessler.cpp` runs `kessler_run` outside EAMxx on synthetic columns
shaped like one rank of a production case, checks the results and times it.
It is not part of the CMake build; `build_driver_kessler.sh` compiles it with
the Kokkos/EKAT install and flags of an existing case build:

```bash
source <case>/.env_mach_specific.sh
./build_driver_kessler.sh --case-build <EXEROOT> --out ./driver_kessler_cpu
OMP_NUM_THREADS=1 ./driver_kessler_cpu --grid ne16-cpu --rain-frac 0.08 --impl all
./driver_kessler_cpu --perf --grid ne16-cpu --rain-frac 0.08 --kernel-timers
./driver_kessler_cpu --preprocess --grid ne16-cpu --repeat 20
./driver_kessler_cpu --help
```

It checks device vs host results, managed vs carved Workspace, garbage- vs
NaN-filled scratch, strided vs contiguous tracers, latent-heat balance,
`relhum`, and the expected rain pattern. `--impl` takes a list of
implementations (`prod`, `hostloop`, and the variants `base`, `i2`, `i6`,
`i2i6`; `both` = `prod,hostloop`, `all` = every one), and each is compared
with the production `kessler_run` bit for bit. `--preprocess` instead
compares and times the old and new `run_impl` pre-processing. Build and
run it on a compute node (an interactive GPU job runs both the CPU and the
GPU binary).

## Validation

- Standalone: all driver checks pass on CPU and GPU, and the current
  `kessler_run` matches the previous version bit for bit
  (performance_suggestions.md, D1 and F).
- Model: the team-per-column `kessler_run` matches the previous C++ build
  bit for bit in HOMME's global state statistics (`homme_atm.log`), rain
  and no rain: ne16, ne30 and ne60 on CPU and GPU, 1 and 4 subcycles, 1
  and 5 days (48 cases; performance_suggestions.md, G and H). The
  skip-dry-pow and pre-processing changes match it bit for bit in the
  same 48 cases (performance_suggestions.md, J).
- Still open: a direct comparison of `precl`/`qr`/`qc`/`theta` with the
  Fortran reference (`fortran_bridge`, branch `ren/kessler-bridge`).

## Status and open items

- The team-per-column `kessler_run` is committed (`223c6719b2`). The
  skip-dry-pow change (`detail::pow_nonneg`) and the new pre-processing in
  `run_impl` are not committed yet; the full C++ perf matrix with them
  (case suffix `_updated`) is bit-for-bit and faster
  (performance_suggestions.md, J).
- `Workspace::f5` / `Buffer::f5` is allocated but unused; `dt0`,
  `time_counter` and `precl_acc` no longer need to be global arrays.
- Remaining costs and possible improvements: performance_suggestions.md,
  section I. A checklist for future Fortran-to-Kokkos ports is in section 0.
