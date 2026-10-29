# EAMxx Kessler Bridge 

## Background: 

See `components/eam/src/physics/crm/pam/external/physics/micro/kessler/Microphysics.h` for background information on EAM implementation. 

Non-CAM Fortran version of code (with C bindings) found in: `components/eam/src/physics/crm/pam/external/physics/micro/kessler/kessler.f90`

CAM Fortran routines are in: [this repo](https://github.com/ESCOMP/atmospheric_physics)

We are using the Kessler provided in: [atmospheric_physics](git@github.com:ESCOMP/atmospheric_physics.git). 

Note that additional initial conditions files might need to be downloaded from https://web.lcrc.anl.gov/public/e3sm/inputdata/atm/scream/init/. 

## Design & Code Layout

The Kessler suite is in `atmospheric_physics/schemes/kessler`. 

List of parameters: 

Kessler.F90:
kessler_init(real lv_in, real pref_in, real rhoqr_in, char* errmsg, int errflg) {r,r,r,w,w}  
kessler_run(int ncol, int nz, real dt, int lyr_surf, int lyr_toa, real cpair, real rair, real rho, real z, &
        real pk, real theta, real qv, real qc, real qr, real precl, real relhum, char* scheme_name, char* errmsg, int errflg) {r,r,r,r,r,r,r,r,r,r,u,u,u,u,w,w,w,w,w}  

Kessler_update.F90:

kessler_update_init(real gravit_in, char* errmsg, int errflg) {r, w, w}  
kessler_update_timestep_init(real temp, real temp_prev, real ttend_t, char* errmsg, int errflg) {r,w,w,w,w}  
kessler_update_run(int nz, int ncol, real dt, real theta, real exner, real temp_prev, real ttend_t, char* errmsg, int errflg) {r,r,r,r,r,r,u,w,w}  
kessler_update_timestep_final(int nz, real cpair, real temp, real zm, real phis, real st_energy, char* errmsg, int errflg) {r,r,r,r,r,w,w,w}  


## GPU stuff

### Implementation Notes

1. For the Fortran link to work, you must add both the C++ side and Fortran side field variables to the `ATMBufferManager` with `Packagename::init_buffers(const ATMBufferManager &buffer_manager)` and `Packagename::requested_buffer_size_in_bytes()` (see the main interface C++ file).
2. It seems like a good idea to create structs for passing input fields, output fields, and (if needed) parameters/etc. Note that this requires each member of the struct that is a field to be defined for both the C++ side (usually with `Pack` and Kokkos views) and the Fortran side (_unmanaged_ view of type `Real` in the same dimensions as the C++ version UNLESS WE ARE USING OPENACC THEN WE USE MANAGED VIEWS WITH `Real`???). 
3. You should create a transpose function for converting between the Fortran and C++ fields. See `transpose()` in the `packagename_functions.hpp` file. 
4. I _think_ you can only pass to Fortran the `get_field_out("fieldname").get_view<Pack**>().data()` object into the C to Fortran binding. 

### GPU vs. CPU performance 

## Misc notes

Example OpenACC Fortran loop. Recall that it is more efficient to loop over `k` first in fortran 

```fortran
    !$acc parallel deviceptr(p_mid) 
    !$acc loop gang vector collapse(2) reduction(max:p_mid_max)
    do k = 1,pver
        do i = 1,ncol
            p_mid(i,k) = p_mid(i,k) * 2.0
        end do
    end do
    p_mid_max = MAXVAL(p_mid)
    !$acc end parallel
```


# Analytic conditions & test setup

StormSPEED uses https://github.com/ESCOMP/CAM/blob/cam_development/src/dynamics/tests/initial_conditions/ic_baroclinic.F90 c/o Jesse N. This sets: 

requires: vcoord,latvals, lonvals, z_int, 
- U: wind at all levels, u
- V: wind at all levels, v
- T: temperature at all levels (T_mid?)
- PS: surface pressure
- PHIS: phis (surface geopotential)
- Q: [qv, qr, qc, qi] need to check on order here
- can use to save z_mid and p_mid too I think

The other option is: https://github.com/ESCOMP/CAM/tree/cam_development/src/dynamics/tests c/o Jesse N.

Need to implement the baroclinic conditions in EAMxx. Question - what do we do about the aquaplanet? Do we use aquaplanet in StormSPEED? 
Below I added how to change to 58 levels to match StormSPEED but we should be able to swap the vertical levels files for 72 or 58 between the two models directly. 

## TO DO

1. In the namelists.xml we need to add something like `initial_conditions_analytic` to use instead of `initial_conditions_filename`. should contain a string to specify name of analytical function we are using ex:
    ```xml
        <initial_conditions>
            <analytic type="array(string)">UNSET</analytic>
            <analytic>"ic_baroclinic"</analytic>
        </initial_conditions>
    ```
For Kessler, we probably want to make a COMPSET with `nlev=58`. Also need to set the vertical coordinates to match StormSPEED with: 
    ```xml
    <!-- Grids manager specs -->
    <grids_manager>
        <type>homme</type>
        <physics_grid_type>gll</physics_grid_type>
        <physics_grid_type hgrid=".*pg2">pg2</physics_grid_type>
        <physics_grid_rebalance>none</physics_grid_rebalance>
        <dynamics_namelist_file_name>./data/namelist.nl</dynamics_namelist_file_name>
        <vertical_coordinate_filename type="file">UNSET</vertical_coordinate_filename>
        <vertical_coordinate_filename nlev="58">/glade/campaign/cesm/cesmdata/inputdata/atm/cam/inic/cam_vcoords_L58_c250227.nc</vertical_coordinate_filename>
        ...
    </grids_manager>
    ```
and in case script we probably want to have `./create_newcase --case ${CASE_NAME} .... --user-mods-dir ${CCSMROOT}/components/eamxx//cime_config/testdefs/testmods_dirs/eamxx/output/preset/2 ${CCSMROOT}/components/eamxx//cime_config/testdefs/testmods_dirs/eamxx/L58`

where `${CCSMROOT}/components/eamxx//cime_config/testdefs/testmods_dirs/eamxx/L58/shell_comands` contains 
```./xmlchange SCREAM_CMAKE_OPTIONS="`./xmlquery -value SCREAM_CMAKE_OPTIONS | sed 's/SCREAM_NUM_VERTICAL_LEV [0-9][0-9]*/SCREAM_NUM_VERTICAL_LEV 72/'`"```

2. in `/glade/derecho/scratch/kstengel/E3SM/E3SM/components/eamxx/src/control/atmosphere_driver.cpp::create_grids()` line 281 add something like:
    ```cpp
    else if (ic_pl.isParameter("analytic")) {
        // Initial run, if ICs are an analytic function, pass the name.
        auto ic_analytic = ic_pl.get<std::string>("analytic");
        gm_params.set("ic_analytic", ic_analytic);
        m_atm_params.sublist("provenance").set("initial_conditions_analytic",ic_analytic);
    }
    ```
    and line 1255 something like:
    ```cpp
    // If analytic is specified, compute initial conditions on all grids
    if (ic_pl.isParameter("analytic")) {
        // Now loop over all grids, and load from file the needed fields on each grid (if any).
        const auto& ic_analytic = ic_pl.get<std::string>("analytic");
        m_atm_logger->info("    [EAMxx] IC analytic type: " + ic_analytic);

        for (const auto& it : m_grids_manager->get_repo()) {
        const auto& grid = it.second;
        const auto& grid_name = grid->name();

        if (ic_fields_names[grid_name].size()==0) continue;

        std::vector<Field> ic_fields;
        for (const auto& fn : ic_fields_names[grid_name]) {
            ic_fields.push_back(m_field_mgr->get_field(fn,grid_name));
        }
        
        // Call the analytic function. I think this should go to a switch statement to call the correct analytic function but should check with EAMxx team about how they would want this setup
        set_fields_analytic(ic_fields,grid,ic_analytic);

        // I don't think we need to worry about IOP enabled options since we should just do the computation at the lat/lon column in the grid
        
        for (auto& f : ic_fields) {
            f.get_header().get_tracking().update_time_stamp(m_current_ts);
        }
        }
    }
    ```

3. compute the analytical initial conditions. see `/glade/derecho/scratch/kstengel/E3SM/E3SM/components/eamxx/src/share/io/scorpio_input.cpp`. 

## Field perturbation

### StormSPEED: 
???

### EAMxx:

Should be able to do automatically with the `perturbed_fields` and corresponding fields in the namelists.xml options. 

## Performance logging (kessler_perf_log)

Call-level timings for `kessler_run`, logged as CSV so the JAX (Python) and
Fortran (CPU + GPU/OpenACC) implementations can be compared directly.

### Python side: `kessler_perf_log.py`

Used by `kessler.py` (the `writeback` step) and `kessler_jax/kessler_run.py`
(the `kessler_run` JAX call, including first-call JIT compilation time).
Import is optional everywhere it's used (`try/except ImportError`), so
removing/renaming this file never breaks the kessler run path -- it just
stops logging.

`log_call()` only ever touches in-process memory (a dict of running
per-label totals, guarded by a `threading.Lock`) -- there's no file I/O per
call. An earlier version opened the CSV and took a cross-process `flock()`
on every call, which serialized every MPI rank against every other rank on
every timestep and measurably slowed runs down.

```python
from kessler_perf_log import log_call
log_call("kessler_run", ncol, nz, dt, elapsed_seconds)
```

The actual file write happens once per rank, in `flush_log()`, called from
`kessler.py`'s `finalize()` -- wired to
`KesslerMicrophysics::finalize_impl()` in the C++ process interface, so it
runs once per rank at simulation end. Each rank's flush appends its
per-label totals to a small raw/intermediate CSV
(`KESSLER_PERF_LOG_PATH` with a `.raw` suffix) under one `flock()`, then
recomputes the aggregate summary across every rank that has flushed so far
and rewrites it to `KESSLER_PERF_LOG_PATH` -- so that path always holds the
final GPTL-style summary directly, with no separate post-processing step:

```
label,count,walltotal_s,wallmax_s,wallmax_rank,wallmin_s,wallmin_rank,callmax_s,callmax_rank
```

- `label` -- caller-supplied tag for what was timed (e.g. `kessler_run`, `writeback`)
- `count`, `walltotal_s` -- total calls and total wall time for that label, summed across all ranks
- `wallmax_s` / `wallmax_rank` -- the rank with the highest cumulative time for that label, and its time
- `wallmin_s` / `wallmin_rank` -- same, for the lowest cumulative time
- `callmax_s` / `callmax_rank` -- the single largest individual call for that label, and which rank it ran on
  (this is what surfaces first-call JIT compilation outliers, since individual calls are no longer logged as separate rows)

Set `KESSLER_PERF_LOG_PATH` to choose the CSV path, e.g.:

```bash
export KESSLER_PERF_LOG_PATH=/glade/scratch/$USER/kessler_perf_log.csv
```

Falls back to `kessler_perf_log.csv` next to `kessler_perf_log.py` (i.e. in
this source directory) if unset -- set it explicitly for any real run so
timing data doesn't pile up in the checkout, and remove any leftover CSV
(and its `.raw` sibling) from a prior run before starting a new one, since
both are appended-to/aggregated-into rather than overwritten from scratch.

If a rank crashes or is killed before `finalize_impl()` runs, that rank's
timings are simply lost -- an accepted tradeoff for a dev/perf tool, not a
science output.

### Fortran side: `kessler_perf_log.F90`

Same design and same `KESSLER_PERF_LOG_PATH` env var as the Python side, for
the Fortran `kessler_run` in `atmospheric_physics/schemes/kessler/` and
`GPU_ports/atmospheric_physics/schemes/kessler/`. **Not yet wired into
either `kessler_run`** -- copy this file alongside those sources and add the
timing/`call log_call(...)` calls, plus a `call flush_log()` at that
scheme's finalize/timestep_final entry point (mirroring `kessler.py`'s
`finalize()` -> `KesslerMicrophysics::finalize_impl()`), when working on
that branch. Depends on `ccpp_kinds` (present in both of those repos, not in
this EAMxx directory).

```fortran
use kessler_perf_log, only: log_call, flush_log
...
call log_call('kessler_run', ncol, nz, dt, elapsed_seconds)
...
call flush_log()  ! once per rank, at finalize
```

Like the Python side, `log_call()` only updates an in-memory per-label
totals table -- no file I/O per call. `flush_log()` does the file I/O,
meant to run once per rank at finalize: it appends this rank's totals to a
raw intermediate CSV (`KESSLER_PERF_LOG_PATH` with a `.raw` suffix) using
raw `open()`/`write()` via `ISO_C_BINDING` with `O_APPEND` for atomic
appends (Fortran's own `OPEN(POSITION='APPEND')` isn't reliably
kernel-atomic across ranks the way a real `O_APPEND` file descriptor is),
plus `flock()` around that append and the read-back below. It then reads
the raw CSV back with plain Fortran formatted I/O, recomputes the aggregate
summary across every rank that has flushed so far, and rewrites
`KESSLER_PERF_LOG_PATH` with that summary -- same schema as the Python
side's summary CSV, so the two are directly comparable:

```
label,count,walltotal_s,wallmax_s,wallmax_rank,wallmin_s,wallmin_rank,callmax_s,callmax_rank
```
## Standalone pybind driver (`driver_kessler_py.cpp`)

For debugging the cost of moving data between EAMxx's fields and JAX
outside the full model. `driver_kessler_py.cpp` builds Kokkos "fields" laid
out like EAMxx's (padded `(ncol, nlev)` views, `qv`/`qc`/`qr` as strided
slices of one tracers view), and calls `kessler.py` exactly the way
`KesslerMicrophysics` does: EAMxx's `PySession` to start Python and import
the module, NumPy views of the field memory built as `create_py_field()`
builds them, then `init`, `run` and `update` with `run_impl`'s arguments.
It is not part of the CMake build.

| file | contents |
|------|----------|
| `driver_kessler_py.cpp` | the driver |
| `driver_kessler_cai.py` | Python side of the zero-copy device path (`--backend cai`) |
| `driver_kessler_passes.py` | per-column sub-cycle pass counts of the JAX `kessler_run` (`--passes`) |
| `driver_kessler_variants.py` | times and checks the `kessler_jax_experimental/` sub-cycle variants (`--variants`) |
| `build_driver_kessler_py.sh` | builds the driver against a JAX case build (its Kokkos/EKAT, flags, Python and pybind11) |
| `check_kessler_py_bfb.sh` | bit-for-bit check of every backend x `py_run_variant` over grids and rain fractions (run before model tests) |

Backends (`--backend host,device,cai`):

- `host`: `py_backend=host` as `run_impl` does it: `sync_to_host` of the 14
  fields (one fence each; strided tracers through a contiguous helper, as
  `Field::sync_to_host` does), `run`, `update`, `sync_to_dev` of the 12
  outputs. Timed per stage: `sync_to_host`, `py_run`, `py_update`,
  `sync_to_dev`.
- `device`: `py_backend=device` as `run_impl` does it (see
  [py_backend=device](#py_backenddevice)): NumPy arrays over the device
  pointers, contiguous copies of the tracers, no host syncs.
- `cai` (GPU only): each field goes to Python as an object with a
  `__cuda_array_interface__`; `jnp.asarray` imports it on the GPU without a
  host copy, `driver_kessler_cai.py` runs the existing device-resident
  bridges (`*_bridge_device`, "contract 2"), and the driver copies the JAX
  outputs back into the fields on the GPU through their own
  `__cuda_array_interface__`. Strided tracers are packed into contiguous
  device scratch first. No host transfers.

Each backend starts every call from the same state, and each backend's
outputs are compared bit for bit with the first backend's. The driver sets
`JAX_PLATFORMS` (`cuda` on a GPU build), `XLA_PYTHON_CLIENT_PREALLOCATE=false`
(so JAX does not reserve 75% of the GPU that Kokkos also allocates from) and
`KESSLER_PERF_LOG_PATH` unless they are already set.

Build and run on a compute node, with the case's modules and the JAX conda
env active:

```bash
source <jax case>/case/.env_mach_specific.sh
module load conda && conda activate /glade/work/kstengel/conda-envs/jax-kessler
./build_driver_kessler_py.sh --case-build <builds>/jax_gpu/build --out ./driver_kessler_py_gpu
./driver_kessler_py_gpu --grid ne16-gpu --backend host,device,cai --warmup 3 --repeat 10
./driver_kessler_py_gpu --grid ne16-gpu --backend device --py-warmup --warmup 0 --per-call
./driver_kessler_py_gpu --grid ne16-gpu --rain-frac 0.08 --passes
./driver_kessler_py_gpu --grid ne16-gpu --variants auto,split,pallas
./driver_kessler_py_gpu --grid ne16-gpu --backend host,device --run-variant current,pallas,split --py-warmup --warmup 0
./driver_kessler_py_gpu --help
```

`--passes` runs `driver_kessler_passes.py` instead of the backends: the
per-column sub-cycle pass counts of the JAX `kessler_run`, the work its
`vmap(while_loop)` does compared with per-column loops, and timings of the
device-resident bridge on all columns, on all columns without rain (one
pass), and on only the multi-pass columns. The pass counts come from a copy
of `kessler_run_core` with a counter in the loop; `kessler_jax/` itself is not
changed, and the report checks the copy against the original bit for bit.

`--run-variant LIST` runs each backend once per `kessler::py_run_variant`
value through `kessler.set_run_variant()` and `kessler.run()`, as the model
does, and compares all of them with the first.

Checks on every run: outputs bit for bit with the first run's; with
`--nsteps N`, calls are chained N at a time (each from the previous call's
output, as model steps are) and every step is compared; every later chain
must repeat the first chain bit for bit (a race, such as a missing barrier
in a kernel, would show up here); every output finite and `qv`/`qc`/`qr`/
`precl` non-negative. `check_kessler_py_bfb.sh --hw gpu|cpu` runs these over
ne16/ne60, rain fractions 0/0.03/0.08/0.12, every backend and every variant
(GPU: `current,vmap,split,pallas,auto` on host/device/cai; CPU:
`current,vmap,split,map,split_map,auto` on host/device), 3 chained steps
twice, and prints a pass/fail table.

`--variants LIST` runs `driver_kessler_variants.py` instead: the sub-cycle
variants in `kessler_jax_experimental/` (see its README.md) against the
current `kessler_run_bridge_device`, with times and a bit-for-bit comparison.

Results (2026-10-29): `GK_perf_utils/workflow/interactive/jax_driver_2026-10-29/RESULTS.md`.

## py_backend=device

`kessler::py_backend` (`namelist_defaults_eamxx.xml`, default `host`)
chooses how the fields reach Python:

- `host`: `run_impl` copies the 14 fields to host memory
  (`sync_to_host`), `kessler.py` copies them to the GPU and back inside
  every call, and `run_impl` copies the 12 outputs back to the device
  (`sync_to_dev`).
- `device`: no host copies. EAMxx passes NumPy arrays whose data pointers
  are device addresses (`create_py_field<Device>`). Those arrays only carry
  the pointer, shape and strides and must never be read as NumPy (that
  dereferences GPU memory on the host and segfaults); upstream EAMxx's
  `cld_fraction_cupy.py` re-wraps them with cupy for the same reason.
  `kessler.set_backend("device")` (called from `initialize_impl`) turns on
  the device path on a GPU: each array is wrapped in a
  `__cuda_array_interface__` object, which `jnp.asarray` imports on the GPU
  without a copy; `run`/`update` call the device-resident bridges
  (`*_bridge_device`) and write the results back with device-to-device
  `cudaMemcpy` through `ctypes` (`libcudart.so.12` from JAX's CUDA plugin,
  so no cupy is needed). Because those arrays must be contiguous,
  `initialize_impl` makes contiguous clones of the strided tracers
  `qv`/`qc`/`qr` (`m_py_contig`), and `run_impl` copies into them before
  the calls (timer `py_pack`) and out of them after (`py_unpack`). On a CPU
  build device pointers are host pointers, and `kessler.py` uses its host
  path.

`initialize_impl` calls `set_backend` before `warmup` (`kessler::py_warmup`),
and on the device path `warmup` passes its dummy fields the way EAMxx does
(NumPy arrays over `cudaMalloc`'d device copies), through the same import
and write-back. That matters for JAX's cache: arrays imported from device
pointers are committed to the GPU, NumPy arrays given to `jnp.asarray` are
not, and `jit` compiles separately for each. Warmed with host arrays, the
device path recompiled all four bridges and `update()`'s `add`/`multiply`
on the first call (~1.3 s at ne16) and `kessler_run`/`timestep_final` again
on the second (~1 s). Check with
`driver_kessler_py_gpu --backend device --py-warmup --warmup 0 --per-call`,
which prints the JIT cache misses (JAX lowerings) of every call; it should
print 0 for every call on both paths.

In the driver, `device` matches `host` bit for bit on GPU (ne16, ne60) and
CPU (ne16), and takes 13.4e-3 vs 22.2e-3 s per call at ne16 and 0.101 vs
0.240 s at ne60 on one A100. To use it in a case:

```bash
./atmchange kessler::py_backend=device
```

and set `XLA_PYTHON_CLIENT_PREALLOCATE=false` in the run environment, as
`JAX_PLATFORMS` is set: otherwise JAX reserves 75% of the GPU next to
Kokkos (not yet tried in a model run).

## Sub-cycle loop: restructuring prototypes (`kessler_jax_experimental/`)

Almost all of the JAX `kessler_run` time is the rain sub-cycle.
`kessler_jax/kessler_run.py` runs each column's sub-cycle as a
`lax.while_loop` and `jax.vmap`s it over the columns. Under `vmap` JAX runs
one loop until the slowest column is done and masks the finished columns,
so every column does the slowest column's number of passes. On the
driver's synthetic columns (`--passes`), every column needs either 1 pass
(dry) or 17-144 (rainy), and the current loop does 12x to 40x the
necessary work (17.7x at 8% rainy columns). On the GPU every loop iteration
also reads the loop flag back on the host and launches several kernels.

`kessler_jax_experimental/` (see its README.md) keeps `kessler_run_core`'s
physics, copied unchanged, and changes only how the loop runs.
`kessler_jax/` itself is not changed, so the LAFT translation stays as
generated. Variants:

- `split`: one pass on all columns (dry columns finish there), then the
  loop on only the columns still running, gathered into a fixed-size
  buffer and scattered back.
- `map`: `lax.map` over columns, so each column runs its own loop, as in
  the Fortran. Meant for the CPU.
- `pallas`: the whole sub-cycle as one Pallas kernel (GPU). Each program
  is one column, its levels are spread over the program's threads, and the
  column state stays in registers: the C++ D1 structure (one team per
  column). Each column stops after its own passes, and the loop's exit test
  never leaves the GPU.
- `auto`: `pallas` on a GPU, `map` otherwise.

Pallas on the A100 compiles through Triton (JAX 0.6.2). It cannot slice
values, so the level k+1 neighbour goes through memory (store,
`plgpu.debug_barrier()`, gather load, barrier). Blocks must be powers of
two, so the 58 levels are padded to 64. float64 `exp`/`pow`/`sqrt` match
XLA's bit for bit, and float64 `dot` crashes Triton, so it is not used.

Results with `--variants` (2026-10-29, deg0053, synthetic columns with 8%
rainy columns, `kessler_run_bridge_device` only, s/call). Every variant
matches the current code bit for bit (GPU at 3%, 8% and 12% rain, CPU at
8%):

| grid     | current  | split    | map      | pallas   | best speedup   | C++ D1   |
|----------|---------:|---------:|---------:|---------:|---------------:|---------:|
| GPU ne16 | 10.9e-3  | 5.8e-3   | -        | 1.5e-3   | 7.3x (pallas)  | 1.15e-3  |
| GPU ne60 | 97.9e-3  | 27.6e-3  | -        | 7.4e-3   | 13.2x (pallas) | 7.56e-3  |
| CPU ne16 | 33.9e-3  | 9.4e-3   | 3.4e-3   | -        | 9.8x (map)     | 3.70e-3  |
| CPU ne60 | 0.511    | 0.138    | 48.6e-3  | -        | 10.5x (map)    | 52.3e-3  |

GPU times vary by about 20% from run to run. At 3% rain Pallas is 8.6x
(ne16) and 20.7x (ne60) faster.

### Using a variant in a case: `kessler::py_run_variant`

`kessler::py_run_variant` (`namelist_defaults_eamxx.xml`, default `current`)
selects the variant at run time, so one build covers all of them:

```bash
./atmchange kessler::py_run_variant=pallas    # or auto, split, map, split_map, vmap
```

`initialize_impl` passes it to `kessler.set_run_variant()` after
`set_backend()` and before `warmup()`, so `kessler::py_warmup` compiles the
selected variant. `kessler.run()` then calls the variant's bridge on both
`py_backend` paths (`update()` is unchanged). With `run.py`:
`--jax-variant current pallas split --py-backend host device`.

Checked through `kessler.py` in the driver (`check_kessler_py_bfb.sh`,
2026-10-29, deg0027): every variant on every backend is bit for bit with
`host/current` at every step of 3-step chains, repeats itself bit for bit,
passes the sanity checks, and has 0 JIT cache misses after warmup, on GPU
and CPU, ne16 and ne60, rain fractions 0, 0.03, 0.08 and 0.12 (16 of 16
cases, 504 run/step comparisons). Per call (all of
`kessler.run` + `update` + data movement), ne16, 8% rain:

| s/call           | current | split   | map     | pallas  |
|------------------|--------:|--------:|--------:|--------:|
| GPU, host path   | 23.7e-3 | 18.9e-3 | -       | 13.1e-3 |
| GPU, device path | 13.5e-3 | 9.7e-3  | -       | 4.1e-3  |
| CPU              | 36.2e-3 | 11.5e-3 | 5.7e-3  | -       |

Not done yet: a model run with a variant, and checking the rain against real
model state (only synthetic columns so far). Details:
`GK_perf_utils/workflow/interactive/jax_driver_2026-10-29/RESULTS.md`,
sections 7 and 8.
