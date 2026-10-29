# Kessler Kokkos Performance Suggestions

Last updated 2026-10-29. Code state: branch `ren/kessler-claude` at
`223c6719b2` (which commits the team-per-column `kessler_run`, D1) plus the
uncommitted skip-dry-pow (I2) and pre-processing (I4) changes; the original
perf numbers in B and C refer to the pre-port code at `49bf6ccdde`. Model
performance numbers come from `GK_perf_utils/kessler_perf_trends.md` and
the run cases in sections G and H; standalone numbers come from
`driver_kessler.cpp` (section F).

**Summary (2026-10-28):** the cause in section B is confirmed. A
team-per-column version of `kessler_run` (fix D1) gives bit-identical
results and, with rain, is 10-33x faster on CPU and 9-22x faster on GPU
(standalone, one rank's worth of columns, ne16 and ne60). It is now
`kessler_run` in `kessler_functions_impl.hpp` (committed in `223c6719b2`);
the old host-loop version is kept for comparison in
`driver_kessler_hostloop.hpp`.

Full model, ne16 1 day, run on an interactive GPU node (section G): the
port is bit-for-bit identical to the previous C++ build (HOMME global
state stats, rain and noRain, CPU and GPU). With rain, `kessler::run` on
the slowest rank drops 3.4-4.0x on CPU (8.35 -> 2.11 s at 64 ranks) and
4.9x on GPU (0.863 -> 0.176 s); noRain is unchanged. Projected onto the
original 128-rank queue run, C++ Kessler would be ~1.0-1.2 s on the
slowest rank vs Fortran's 1.38 s.

Queue runs of the full perf matrix (section H, 48 cases, all
bit-for-bit with the original runs): C++ CPU Kessler is now faster than
Fortran in every case. At ne60 CPU with rain, the case the problem was
reported on, `kessler::run` went from 65.5 to 9.5 s/day (steady state,
slowest rank), against Fortran's 11.8.

After that, two smaller bit-for-bit changes were made and checked in ne16
model runs (section I): skipping `pow` for dry levels (I2, `kessler_run`
1.14x faster on CPU with rain) and faster pre-processing in `run_impl`
(I4, `preprocess` 1.2x on CPU, 1.8-2.6x on GPU). The full C++ matrix was
rerun with them under the same suffix `_updated` (section J): all 48
cases bit-for-bit, and `kessler::run` on the slowest rank another
1.15-1.33x faster than D1 alone on CPU and 1.03-1.14x on GPU. ne60 CPU
with rain is now 7.4 s/day (steady state), against 65.5 originally and
11.8 for Fortran. The D1-only results are kept in section H and
`results_updated_2026-10-28.txt`.

## Contents

- [0. Translation guide: Fortran -> C++/Kokkos (for LLM-assisted ports)](#0-translation-guide-fortran---ckokkos-for-llm-assisted-ports)
- [A. Current state: what has been done, and corrections to earlier notes](#a-current-state)
- [B. The CPU performance problem (C++-CPU was 3-5.5x slower than Fortran-CPU with rain)](#b-the-cpu-performance-problem)
- [C. Confirming the cause (done)](#c-confirming-the-cause-done-2026-10-28)
- [D. Possible fixes, in priority order](#d-possible-fixes-in-priority-order)
- [E. Other open items](#e-other-open-items)
- [F. Standalone driver: how to run it, and results](#f-standalone-driver-how-to-run-it-and-results)
- [G. Full-model check of the port on an interactive node (2026-10-28)](#g-full-model-check-of-the-port-on-an-interactive-node-2026-10-28)
- [H. Queue runs of the port (suffix `_updated`, 2026-10-28)](#h-queue-runs-of-the-port-suffix-_updated-2026-10-28)
- [I. Review after the port: remaining costs and possible improvements](#i-review-after-the-port-remaining-costs-and-possible-improvements)
- [J. Queue runs with I2 + I4 (suffix `_updated`, rerun 2026-10-28)](#j-queue-runs-with-i2--i4-suffix-_updated-rerun-2026-10-28)

---

## 0. Translation guide: Fortran -> C++/Kokkos (for LLM-assisted ports)

Lessons from porting Kessler, written as instructions for the next
port of a Fortran physics scheme into EAMxx C++/Kokkos with an LLM.
Each rule names the Kessler problem it would have prevented. Give
this section to the LLM at the start, together with the Fortran
source and one existing EAMxx scheme to copy patterns from (P3 or
SHOC).

### 0.1 Before writing any code

- Read the Fortran driver and the scheme together and write down, for
  every argument: units, which index is the surface, whether it is
  per-column or per-level, and what it really is. Do not trust names:
  Kessler's Fortran `dz` argument is `z_mid` (heights), not a layer
  thickness, and the first port got it wrong.
- Write down the Fortran loop structure, especially any loop that runs
  a data-dependent number of times (sub-cycling, iteration to
  convergence) and where it sits relative to `do col`. That placement
  decides the Kokkos design (0.2). Kessler's `do while` is inside
  `do col`; the first port moved it outside all columns, which made
  C++-CPU 3-5.5x slower than Fortran with rain (sections B-D).
- List every error exit, abort and limit in the Fortran (e.g.
  "bad time splitting", `dt0 < 1e-12`). Each must exist in the port.
- Record the index convention: Fortran 1-based, C++ 0-based; EAMxx
  levels are top-down (0 = TOA). State `lyr_surf`/`lyr_toa` explicitly.
- Note the compiler flags the comparison builds use (derecho:
  `-O2 -Mnofma -Mflushz -Kieee`, and `--fmad=false` on GPU). Bit-for-bit
  comparisons are only meaningful with matching flags.

### 0.2 Requirements for the translated code

**Structure**

- Keep the Fortran's per-column independence. A column loop with a
  per-column while/do-while inside becomes ONE `TeamPolicy` kernel, one
  team per column, with the while loop inside the team. Never turn it
  into a host loop over passes that sweeps all columns: the work
  becomes `max(passes) * ncols` instead of `sum(passes)`, and every pass
  adds launches and a host sync.
- Do not "mask out" finished columns by computing and multiplying by
  0/1 (`msk * new + (1 - msk) * old`). That still does all the work. If
  a column is finished, it should do nothing.
- Use `ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(ncols, nlev_packs)`
  like P3/SHOC, not `Kokkos::AUTO`. On A100 with Pack 1 and 58 levels it
  gives 64 threads; AUTO gave 128 and was ~40% slower without rain at
  ne60 (F).
- Level loops: `TeamVectorRange` over packs; reductions over levels:
  team `parallel_reduce`. Vertical neighbours across pack boundaries:
  `ekat::shift_left`/`shift_right` (as P3's upwind code does).
- Fuse loops that only touch the same pack index (e.g. micro + fall
  speed + CFL limit) and put `team_barrier()` only where a loop reads
  another pack's data written by the previous loop.
- Per-column scalars that change during the loop (`dt0`, elapsed time,
  accumulators, counters) are team-local variables, not global views.
- Errors inside kernels: set a per-column error code, reduce once on
  the host after the kernel, then throw the same message the Fortran
  would. No host sync inside the loop.

**Data and memory**

- Scratch arrays come from a caller-owned Workspace carved from the
  ATMBufferManager (P3/SHOC Buffer pattern), not allocated per call.
  Do not size or carve arrays the algorithm does not use (Kessler's `f5`).
- Initialise every scratch value the kernel reads. Kessler's first GPU
  build read uninitialised `mask`/`dt0` (`mask = 6966.7`) and never
  converged.
- Never `ekat::scalarize` a view that may be strided. EAMxx tracers (qv,
  qc, qr, ...) are subviews of one bundled tracers field, so their
  column stride is `ntracers*nlev`. `scalarize` drops the stride and reads
  other tracers/columns; that produced wrong `precl` in a 10-day run
  while the state looked fine. Read single levels through the packed
  view: `q(col, k / Pack::n)[k % Pack::n]`. Add `EKAT_REQUIRE_MSG`
  contiguity checks for anything that is scalarized.
- Pack padding lanes (nlev not a multiple of `Pack::n`) and excluded
  boundary levels must not divide by zero or call pow on a bad value:
  set guarded inputs (`x.set(!valid, 1)`) before the operation; the fpe
  build traps even if the result is never used.
- Use `ekat::` math on Packs and `Kokkos::` math on scalars.
  `Kokkos::max/min/pow/exp` only take plain numbers, and `Kokkos::max`
  compares with `<`, which on a Pack gives a per-lane Mask, so
  `Kokkos::max(pack, 0.0)` does not compile. The ekat versions are
  inlined per-lane loops: `ekat::exp/sqrt/abs/log` call `Kokkos::exp` etc.
  (in an `omp simd` loop), `ekat::pow` calls `std::pow`, and
  `ekat::max/min` call `ekat::impl::max/min`. So there is no speed or
  answer difference between the two; on GPU (Pack 1) `ekat::pow` is one
  pow call. For scalars (e.g. per-column `dt0`, CFL minimum) use
  `Kokkos::min/max/abs`, not `std::`, so the code also compiles for GPU.
  One catch: with a NaN input `ekat::max` differs between CPU and GPU.
  On CPU it is `std::max`, `(a < b) ? b : a`, so `max(NaN, 0) = NaN`; on
  GPU ekat uses `(a > b) ? a : b`, so `max(NaN, 0) = 0`. A clip such as
  `qr = ekat::max(qr, 0)` therefore passes a NaN through on CPU but
  silently turns it into 0 on GPU, which can hide an upstream problem in
  GPU runs only. When chasing NaNs, compare a CPU run, or check for NaN
  before the clip.
- Physical constants are passed as `Real` scalars; derived constants
  (cp/R etc.) are computed once on the host.

**Faithfulness**

- Port every per-iteration step exactly where the Fortran has it
  (Kessler re-checks `dt0` against CFL every sub-cycle; computing it
  once up front violates CFL when the fall speed grows).
- No "TEMPORARY" limits or shortcuts (the first port had an
  `n_iter < 10` cap). If one is needed, put it in a tracked TODO here
  and in the code, and remove it before any perf run.

### 0.3 Validation requirements (build these with the port, not after)

- A standalone driver (like `driver_kessler.cpp`) that runs the routine
  on synthetic columns at production size per rank, with:
  - device vs host instantiation (to round-off),
  - Workspace managed vs carved (bit for bit),
  - Workspace poisoned with garbage vs NaN (bit for bit: catches
    reads of uninitialised scratch),
  - tracers as strided subviews vs contiguous (bit for bit: catches
    dropped strides),
  - physical checks (conservation such as the latent-heat balance,
    non-negative water, known no-op cases),
  - a realistic fraction of "active" columns (here: rain in 3-12% of
    columns, varying amounts), not every column or every other one.
- When restructuring for performance, keep the previous version in
  the driver (as `driver_kessler_hostloop.hpp`) and require bit-for-bit
  agreement with it on CPU and GPU before any model run.
- Model check: same rank layout and node type, compare the HOMME
  global state statistics (`homme_atm.log`) between old and new builds;
  they must match exactly. Then compare against the Fortran reference
  (`fortran_bridge`) for the scientific check.

### 0.4 Performance requirements

- Write down the work model of the Fortran and of the port (e.g. "sum
  over columns of passes" vs "max * ncols") before measuring. If they
  differ, the port is wrong regardless of how fast each kernel is.
- Benchmark early, at one rank's production size (ncols per rank,
  nlev, dt), on both CPU and GPU, with and without the expensive regime
  (here rain). A port that ties Fortran only in the cheap regime is not
  done.
- Time with a per-kernel breakdown (driver `--kernel-timers`) so you see
  whether time is in real work, launches/syncs, or wasted work.
- In the model, look at the slowest rank (timer max), not only the
  mean: physics load is spatially clustered.
- Remember the EAMxx overhead around the kernel (pre/postcondition
  checks, preprocessing) when judging speedups (section I).

### 0.5 Working process with an LLM

- Give the LLM: the Fortran source, its driver/caller, the EAMxx
  process interface it must fit, this section, and one reference
  scheme (P3/SHOC) to copy idioms from.
- Ask for the Kokkos design (kernel/policy structure, which loops are
  per column, what is team-local) and the work model BEFORE code, and
  review it against 0.1-0.2.
- Translate in small steps, each followed by the driver checks in 0.3.
  Do not accept "done" for an item without a test or a grep that shows
  it in the code: earlier versions of this file listed four
  optimisations as done that were never committed (see A,
  "Corrections").
- Keep this file in step with the code: record the commit/build each
  result comes from, and mark items done only with evidence.
- Keep production APIs clean: debugging/benchmark hooks (pass counts,
  reference implementations) live in the driver, not in production
  signatures.
- Use builds tagged by commit + diff (`workflow/build.sh`) so older run
  cases keep their executable, and compare old vs new builds on the
  same hardware and layout.

### 0.6 Checklist before calling a port done

- [ ] Loop structure matches the Fortran (per-column iteration inside a
  per-column kernel); work model written down and equal to Fortran's.
- [ ] Every Fortran error exit and limit ported; no temporary limits.
- [ ] Argument meanings, units, index direction verified against the
  Fortran caller.
- [ ] No scalarize on strided views; contiguity asserted where used.
- [ ] All scratch initialised; poison and NaN runs match bit for bit.
- [ ] Device vs host, managed vs carved, strided vs contiguous all pass.
- [ ] Physical checks pass (conservation, non-negativity, known no-ops).
- [ ] Standalone timing at production size, CPU and GPU, both regimes,
  within expectations of the work model.
- [ ] Model run bit-for-bit vs the previous build (if refactoring) and
  compared with the Fortran reference.
- [ ] This file updated with commit, build and results.

### 0.7 Kokkos lessons and best practices

General Kokkos/EKAT lessons from the Kessler work (D1, I2-I6, and the
timing work in F-J), for any EAMxx physics kernel. The rules in 0.2
say what the port must do; this section says why, and what was learned
about Kokkos along the way.

**Parallel structure**

- Match the parallel structure to where the work is. Data-dependent
  work per column (sub-cycling, iteration) belongs inside one team per
  column, so each column does only its own passes. Moving that loop to
  the host cost 10-33x on CPU and 9-22x on GPU in Kessler (F), more
  than any kernel-level tuning could win back.
- Kernel launches and host syncs are not free. The old host loop did
  ~9 launches and a reduction to the host per pass; on GPU the gain from
  I4 was mostly one launch fewer per call (J). Fuse kernels that run over
  the same columns, and keep reductions to the host to one per call,
  after the main kernel.
- A team-level `parallel_scan` (e.g. `PF::calculate_z_int`) needs the
  policy from `TPF::get_thread_range_parallel_scan_team_policy`, not the
  default team policy; use that policy for the whole kernel that
  contains the scan, as SHOC's preprocess does (I4).
- Team size: use `ekat::TeamPolicyFactory`, not `Kokkos::AUTO`. For
  short columns (58 levels, Pack 1 on GPU) AUTO picked 128 threads,
  half of them idle; 64 was up to 1.5x faster (F). Measure team size
  when the level count changes a lot.
- Values every thread in a team computes from team-wide inputs (the
  result of a team reduction, a scalar read from a view) can be
  team-local variables, one identical copy per thread, with no
  broadcast. Write per-column results once, inside
  `Kokkos::single(Kokkos::PerTeam(team), ...)`.
- `team_barrier()` only where a loop reads data another thread wrote
  in the previous loop (e.g. a neighbouring pack). Fusing loops that
  touch only their own pack removes barriers as well as passes.

**CPU and GPU from one source**

- On the derecho CPU build the backend is Serial with Pack 16: a team is
  one thread and `TeamVectorRange` is a plain loop, so all of the CPU
  speed comes from the Pack (SIMD) and from not doing wasted work. On
  GPU, Pack 1, the parallelism is the threads. Code written over packs
  with `ekat::` math runs well on both.
- nvc++ vectorizes `ekat::pow`/`exp` on a 16-wide pack into four AVX2
  vector calls (`__pd_pow_4`, `__pd_exp_4`), even with `-Mnofma -Kieee`
  (I3). Check this in the binary (`nm`/`objdump`) rather than from
  `-Minfo` alone: its line numbers are ambiguous across inlined headers.
- Skipping work per pack (I2: no `pow` when a whole pack has no rain)
  pays on CPU, where one test skips 16 lanes (1.15-1.33x in the model).
  On GPU each thread tests its own level and threads in a warp wait for
  the slowest one, so the same change gave only 1.03-1.14x (J). Test a
  whole pack (`(x != 0).any()`), never branch per lane inside a pack.
- Transcendental functions dominate the cost of physics kernels once
  the structure is right. Reusing a value instead of recomputing it
  (I4: `theta = T / exner` with the Exner value just computed, instead
  of a second `pow`) gave 1.85x on CPU where packing and fusing alone
  gave 1.1x. Look for repeated `pow`/`exp`/`log` first; keep the
  arithmetic identical if the change must be bit-for-bit.
- Team scratch (level-0 `team_scratch`, shared memory on GPU) did not
  help (I6): a column's per-level arrays already stay in cache on CPU
  and L1/L2 on GPU. Do not expect it to pay unless the data is reused
  across many threads or does not fit in cache.
- nvcc restrictions on extended (`KOKKOS_LAMBDA`) lambdas: a variable
  cannot be captured for the first time inside an `if constexpr` branch
  (make it outside, as `driver_kessler_variants.hpp` does). In a class
  member, use `KOKKOS_CLASS_LAMBDA` and copy members to local variables
  first, so the device code does not read through `this`.

**Measuring**

- Fence before stopping a timer. Kokkos kernels are asynchronous on GPU,
  so a timer stopped without `Kokkos::fence()` measures the launch, and
  the kernel is charged to the next fenced timer (`run_impl`'s
  `preprocess` timer did this until 2026-10-29, I4). The driver's `--kernel-timers` fences after
  every kernel, so use it for the breakdown and an unfenced run for the
  total.
- Time a standalone driver at one rank's production size (columns per
  rank, levels, Pack size, team size) on both backends, with warm-up
  runs and best-of-N, before running the model. Then check the model's
  slowest rank: the driver numbers carried over (I2 on CPU: 1.3x standalone,
  1.15-1.29x in the model with rain), but the model adds EAMxx overhead
  (property checks, preprocessing) the driver does not see.
- Keep every experimental variant in the driver behind a compile-time
  switch and compare it bit for bit with production
  (`driver_kessler_variants.hpp`, `--impl all`). This made it cheap to
  test I2 and I6 and drop the one that did not help.

---

## A. Current state

Done and present in the code:

1. **Scratch view pre-allocation** (Plan A from the 2026-09-22 notes).
   The 10 scratch arrays are in `KMF::Workspace` (`kessler_functions.hpp`)
   and are passed to `kessler_run` as `const Workspace& ws`. The process
   interface carves them out of one ATMBufferManager allocation in
   `requested_buffer_size_in_bytes()`/`init_buffers()`
   (`eamxx_kessler_process_interface.cpp`), the same way P3 and SHOC
   do. `driver_kessler.cpp` can build the workspace either way
   ("managed" or "carved") and checks that both give the same answer.
   The 09-22 call-site and z-vs-dz problems are fixed: `kessler_run`
   now takes `z_mid`, and `run_impl` passes `pref`, `ncols` and `pk`.
2. **cpair/rair are scalars.**
   `cpair` and `rair` are passed to `kessler_run` as `Real` scalars (there is
   no KesslerData struct). `xk = cp/R` and `f5_dk = 4093*lv/cp` are computed
   once on the host.
3. **dt0 is re-checked against CFL every sub-cycle.**
   `kessler_recompute_dt0` does this, matching the Fortran reference. The
   "bad time splitting" abort (`dt0 < 1e-12`) has also been restored.
4. **Standalone driver extended for performance work (2026-10-28).**
   `driver_kessler.cpp` can time `kessler_run` against the previous
   host-loop version (`driver_kessler_hostloop.hpp`), check that they
   agree bit for bit, and report per-column sub-cycle counts and
   per-kernel times. See section F.
5. **D1 ported (2026-10-28, committed in `223c6719b2`):** `kessler_run` in
   `kessler_functions_impl.hpp` is now one team per column. See D1. The
   Workspace is unchanged so far.
6. **New builds and interactive-node cases (2026-10-28):** the port is built
   as `builds/cpp_{cpu,gpu}_7ecd5d94bb-dirty9022632f` (`workflow/build.sh`;
   `builds/cpp_cpu` and `builds/cpp_gpu` now point at them). The original
   builds (`cpp_{cpu,gpu}_49bf6ccdde-dirty349ee6c7`) are untouched. 12
   ne16 1-day cases (`_port_*` and `_base_*`) were run on interactive node
   deg0053. See section G.
7. **I2 and I4 (2026-10-28, uncommitted):** `kessler_run` skips `pow` on
   all-dry packs (`detail::pow_nonneg`), and `run_impl`'s pre-processing
   is one packed team-per-column kernel. Built as
   `builds/cpp_{cpu,gpu}_223c6719b2-dirtyf15739b6`; bit-for-bit in ne16
   interactive runs and in the full 48-case queue matrix (section J). I3
   (vectorization) and the I6 team-scratch variant were checked and
   left out (no gain). See section I.

**Corrections** (written against the pre-port code; the D1 port replaced
the host loop, so the first three items below no longer apply). The
previous version of this file said the items below were done, but they
were NOT in the code:

- Old item 2, `CHECK_INTERVAL=4`. The convergence reduction
  (`kessler_converge`) still runs every iteration, and each run is a
  device-to-host sync.
- Old item 7, mask-guarded early exit. No sub-cycle kernel skips
  converged columns. `kessler_precl_accum`, `kessler_sed_inner`,
  `kessler_sed_toa`, `kessler_micro`, `kessler_velqr_update` and
  `kessler_recompute_dt0` all do the full computation for every
  column on every iteration, including every pow/exp call.
  `kessler_micro` only uses `msk` to throw the result away for
  converged columns.
- Old item 8, fence removed before the sub-cycle loop. Both fences
  are still there: one after `kessler_init_fields` and one after
  `kessler_init_dt0`. `kessler_check_dt0` (added later for the
  bad-time-splitting abort) also reduces into a host scalar, which
  syncs anyway.
- "f5 was dropped." `ws.f5` / `Buffer::f5` is still in the Workspace and is
  still sized and carved (`num_2d_vector = 6`), but `kessler_run` never
  uses it. That wastes `ncols*nlev_packs` Packs of buffer memory.
- `driver_kessler.cpp`'s `--help` text still talked about a "TEMPORARY
  `n_iter < 10`" limit. Fixed 2026-10-28: it now says 100000.

---

## B. The CPU performance problem

Symptoms (ne60, s per simulated day, subcycles = 1, `kessler::run`):

|             | Fortran | C++  |                            |
|-------------|--------:|-----:|----------------------------|
| cpu, noRain |     3.5 |  3.1 |                            |
| cpu, rain   |    11.8 | 65.5 | 5.5x slower                |
| gpu, rain   |    20.0 |  6.8 | C++ is the fastest option  |

- The gap shows up only with rain, and it grows with grid size:
  about 3-4x at ne16/ne30 and about 5x at ne60.
- Without rain, the C++ and Fortran CPU kernels take the same time
  (e.g. 0.13 vs 0.13 s/day at ne30 sc1).
- Rain cost is set by the number of CFL sub-steps, not by how many
  times Kessler is called (trend 7 in `kessler_perf_trends.md`).

**Cause** (CONFIRMED 2026-10-28, see C1/C3/C5; FIXED by the D1 port, see
G/H): sub-cycling ran as one loop over all columns instead of one loop
per column. The description below is of the pre-port code.

The Fortran reference (`atmospheric_physics/schemes/kessler/kessler.F90`)
runs its sub-cycle `do while` loop INSIDE `do col = 1, ncol`. Each
column does only its own sub-cycles: a dry column has `velqr = 0`, so
`dt0 = dt` and it finishes in 1 pass.

    Fortran work = SUM over columns of each column's sub-cycle count

The pre-port C++ code ran the `while (!all_converged)` loop on the
host, outside every kernel (kept for reference in
`driver_kessler_hostloop.hpp`). Each pass sweeps all `ncols * nlev_packs`, and
converged columns are not skipped (see the corrections in A).

    C++ work = MAX over columns of the sub-cycle count * ncols

With rain in 3-12% of columns, `max * ncols` is several times the sum.
This fits every symptom above:

- noRain: max = sum = 1, so C++ ties with Fortran.
- A bigger grid has more columns per rank, so the worst column on
  each rank needs more sub-cycles.
- On GPU all columns run at once, so the time is set by the slowest
  column either way. The extra work lands on idle threads, which
  is why the GPU number looks good. (The D1 prototype shows GPU
  also gains 9-22x, from dropping ~9 launches and a host sync per
  pass, see F.)

Secondary CPU costs of the pre-port code. The measurements in C3 show
these were small next to the wasted work itself:

- Memory. Each pass runs about 7 separate kernels, and each one
  streams the full (ncols, nlev) fields. Fortran keeps one column
  (nlev levels) in L1 cache for all of its sub-cycles.
- A host sync every pass (`kessler_converge`), plus kernel launch
  overhead. Measured: < 1% of CPU time. The CPU build uses Kokkos
  Serial, so there is no OpenMP fork/join.
- The RangePolicy kernels compute `col = idx / nlev_packs` and
  `kp = idx % nlev_packs` on every iteration.
- The pow/exp calls in `kessler_micro` and `kessler_velqr_update`
  run on dry columns too, where qr = 0. This is where the wasted
  time goes: `kessler_micro` is ~75% and `kessler_velqr_update` ~14%.

---

## C. Confirming the cause (done 2026-10-28)

All measured with `driver_kessler.cpp` (section F), not in the model.
The per-column counts come from the D1 prototype, which is
bit-identical to production, so production's `n_iter = max(npass)`.

- [x] **C1. Count the extra column passes.**
  Done in the driver instead of behind a `KESSLER_SUBCYCLE_STATS`
  flag in `kessler_functions_impl.hpp`: the D1 prototype writes each
  column's sub-cycle count (`npass`), and the driver prints max,
  mean, prod vs team column-passes, and a histogram.
  Result: most columns need 1 pass; rainy ones need 17-256. The
  worst column needs ~140 passes at both ne16 and ne60.

  | rain frac | prod/team column-passes, ne16 cpu | ne60 cpu |
  |----------:|----------------------------------:|---------:|
  | 0.03      | 33.1                              | 41.2     |
  | 0.08      | 17.5                              | 17.8     |
  | 0.12      | 12.1                              | 12.2     |

  Expected waste >> 1 with rain and = 1 without: confirmed.

- [x] **C2. Make `driver_kessler.cpp`'s rain look like the real runs.**
  - `--rain-frac F`: a fraction F of the columns is rainy, spread
    evenly, each with its own rain/cloud amount (log-uniform
    0.25-2x). `--rain` is now `--rain-frac 0.5`.
  - `Kokkos::Timer` around every `kessler_run` call; `--repeat N` and
    `--warmup N`, with the input state restored (untimed) before
    each repeat.
  - `--grid ne16-cpu|ne16-gpu|ne60-cpu|ne60-gpu` sets ncols/nz/dt
    to one rank of the production cases (58 levels, dt 1800 s;
    13826 or 194402 columns over 128 CPU or 4 GPU ranks).
  - Built with `build_driver_kessler.sh --case-build` pointing at
    `builds/cpp_cpu_*` and `builds/cpp_gpu_*` (section F).

  Not yet done: checking the synthetic rain against real model
  state (sub-cycle counts from an actual dcmip2016 day 10-15).

- [x] **C3. Per-kernel CPU breakdown.**
  `--kernel-timers`, built into the driver with Kokkos tool
  callbacks (no Kokkos Tools install on derecho). ne60-cpu,
  rain frac 0.08, production `kessler_run`:

  | kernel                  | share                             |
  |-------------------------|-----------------------------------|
  | `kessler_micro`         | 76.1%                             |
  | `kessler_velqr_update`  | 13.9%                             |
  | `kessler_sed_inner`     | 7.6%                              |
  | `kessler_recompute_dt0` | 1.4%                              |
  | everything else         | < 1% (incl. `kessler_converge`)   |

  ne16-cpu is the same within a few percent. Time is spread
  over the passes (wasted work), not in sync/launch overhead.

- [x] **C4. CPU build configuration** (`builds/cpp_cpu_49bf6ccdde-dirty*`):
  - Compiler nvc++ (nvhpc 26.1), `-O2 -Mnofma -Mflushz -Kieee`.
  - `SCREAM_PACK_SIZE = 16` with 58 levels: 4 packs, 6 padding
    lanes (~10% wasted lanes).
  - Kokkos backend is Serial (the CMakeCache lists
    `Kokkos_ENABLE_OPENMP`, but the installed Kokkos is the
    "nothreads" one), `NTHRDS_ATM = 1`, 128 ranks per node.
  - Not done: compiler vectorization report for `kessler_micro`
    (`-Minfo=vect` for nvc++). Lower priority now.

  GPU build: CUDA sm_80, `SCREAM_PACK_SIZE = 1`, 4 ranks.

- [x] **C5. Bound.** The D1 prototype measures this directly instead of
  estimating it: its speedup follows the C1 waste ratio closely
  (e.g. ne16-cpu, 8% rain: waste 17.5, speedup 14.7x). The loop
  structure is nearly the whole story. Projected per day
  (48 calls/day, 8% rain):
  - ne60 cpu: prod 0.755 s/call = 36 s/day; team 2.5 s/day
    (Fortran in the model: 11.8 s/day)
  - ne16 cpu: prod 0.054 s/call = 2.6 s/day; team 0.18 s/day
    (model: C++ 3.6, Fortran 0.9 s/day)

  (Checked in the model, H: ne60 cpu rain is 9.5 s/day steady state
  after the port, not 2.5. The standalone projection assumed evenly
  spread rain; in the model the rain is clustered, so the slowest rank
  does much more real work.)
  The standalone ne60 prod number (36 s/day) is below the model's
  65.5. The driver runs one rank alone on a node; in production
  128 ranks share cache and memory bandwidth, which hurts prod's
  repeated full sweeps most. Check in the model after the port.

---

## D. Possible fixes, in priority order

- [x] **D1. Move the sub-cycle loop into the kernel, one team per column, like
  Fortran.** (PORTED 2026-10-28, committed in `223c6719b2`; ne16 full-model check done on
  an interactive node, G; full queue matrix done and bit-for-bit, H.)
  Now `kessler_run` itself (`kessler_functions_impl.hpp`); it was
  prototyped as `kessler_run_team` in a driver header first. Same
  arguments and Workspace as before. Structure:
  - One kernel (`"kessler_run"`) with
    `ekat::TeamPolicyFactory::get_default_team_policy(ncols, nlev_packs)`,
    as P3/SHOC do. Each team runs the
    init, its own `while (active)` sub-cycle loop, and relhum.
  - Level loops are `TeamVectorRange` over packs; the CFL min is a
    team `parallel_reduce`; the `shift_left`/`shift_right` stencil is
    unchanged.
  - Per pass: precl accumulate (every thread, replicated);
    sed (inner + TOA, one loop); `team_barrier`; micro fused with
    the velqr update and the CFL min (`velqr(kp)` and its CFL limit
    only depend on `qr(kp)`); time update; `team_barrier`.
  - `dt0`, `time_counter`, `precl_acc`, the pass count and the error
    code are team-local scalars, kept identically by every thread
    (all their inputs are team-wide), so no broadcast is needed.
    Final values are written to `ws` at the end.
  - `ws.mask` returns 0 (converged), 2 (bad time splitting) or 3
    (cannot converge: `dt0 <= 0` or >= 100000 passes). One
    `parallel_reduce` after the kernel throws the same errors as
    `kessler_run`.

  Validation: matches the old host loop bit for bit (rtol = atol = 0)
  on CPU (Serial, Pack 16) and GPU (CUDA, Pack 1), for the
  managed, carved, NaN-poisoned and strided-tracer runs, and the
  host run. All driver checks pass.

  Results (s per call, best of 10; full tables in F):
  - ne60 cpu, 8% rain: 0.755 -> 0.0523 (14.4x)
  - ne60 gpu, 8% rain: 0.0897 -> 0.0060 (14.8x, team size 64)
  - no rain: about the same on CPU (0.98x); 1.10x faster on GPU
    with team size 64.

  GPU team size: use 64, not AUTO. With Pack 1 and 58 levels, 64
  covers the column in one sweep. AUTO picks 128 and made the
  no-rain case 0.71x the speed of prod at ne60; 256 is ~2x
  slower again. On CPU, AUTO gives a team size of 1. The
  TeamPolicyFactory default gives exactly this: min(128,
  nlev_packs rounded up to 32) on GPU (64 for 58 levels), 1 on
  CPU.

  After the port (ported `kessler_run` vs old host loop, ne16,
  bit for bit for 3 steps at rain frac 0, 0.08 and 0.5, all
  driver checks pass on CPU and GPU):
  - cpu, 8% rain: 5.40e-2 -> 3.59e-3 s/call (15.0x)
  - gpu, no rain: 1.90e-4 -> 1.14e-4 s/call (1.67x)
  - gpu, 8% rain: 1.45e-2 -> 1.37e-3 s/call (10.6x)

  Full model (section G): bit-for-bit identical to the previous
  C++ build at ne16 on CPU (64 and 128 ranks) and GPU (4 ranks),
  rain and noRain. `kessler::run` on the slowest rank with rain:
  CPU 3.4-4.0x, GPU 4.9x faster.

  Still to do:
  - [done, H] Queue runs matching the original perf cases, ne16/ne30/ne60.
  - Compare against the `fortran_bridge` reference (the port only
    matches the previous C++ build, which had its own validation).
  - Shrink the Workspace: `dt0`, `time_counter` and `precl_acc` are
    team-local and only written at the end (keep `mask`, or
    replace it with a small error view); drop `f5` (see D6).
    `velqr` and `sed` could move to team scratch; `r`, `rhalf` and `pc`
    could be computed on the fly. (Re-run `driver_kessler.cpp` and
    the G cases after this; it changes the process interface.)

- [ ] **D2. (Stopgap) Make converged columns cheap.** Return early when
  `mask(col) == 0` in the sub-cycle kernels. NOT NEEDED if D1 is
  ported: D1 already does only each column's own passes.

- [ ] **D3. (Stopgap) Compact the active columns.** NOT NEEDED if D1 is
  ported.

- [ ] **D4. Reduce the host syncs** (check convergence every N passes).
  NOT NEEDED: the per-pass sync is < 1% of CPU time (C3), and D1
  removes it.

- [ ] **D5. Fuse kernels within a pass.** Done as part of D1 (micro +
  velqr update + dt0 recompute, and sed inner + TOA, are fused).

- [ ] **D6. Small cleanups:**
  - [x] Remove the fences after `kessler_init_fields` and
    `kessler_init_dt0` (done by the D1 port).
  - [ ] Remove `ws.f5` / `Buffer::f5` (set `num_2d_vector = 5`; update
    `init_buffers`, `requested_buffer_size_in_bytes` and
    `driver_kessler.cpp`'s `make_workspace` / `poison_workspace`).
  - [ ] Use MDRangePolicy or TeamPolicy instead of the `idx /` and `%`
    in the remaining RangePolicy kernels
    (`kessler_update_timestep_init`, `kessler_update_run`,
    `kessler_update_timestep_final`).

---

## E. Other open items

- Old item 6, TeamPolicy for the dt0 recompute. Superseded by D1.
- D1 is ported, passes `driver_kessler.cpp` on CPU and GPU, and is
  bit-for-bit with the previous C++ build in ne16 model runs (G) and in
  every queue run (H), including ne60 CPU (65.5 -> 9.5 s/day, below
  Fortran's 11.8). Still to do: compare against the `fortran_bridge`
  reference.
- I2 and I4 (section I) are in the code, uncommitted, and bit-for-bit in
  the full C++ queue matrix (section J; 1.15-1.33x more on CPU, 1.03-1.14x
  on GPU over D1 alone). The `preprocess` timer now fences (2026-10-29,
  see I4). Still to do: commit them, and port I4 and the timer fence to
  the bridge and JAX branches so all three implementations have the same
  pre-processing cost and timers.
- Further improvements after the port: see section I.
- Load imbalance across ranks (G3): with the port, the rainiest CPU
  rank's Kessler time is ~7x the mean. That is real work (most of that
  rank's columns rain), so Kessler itself has little left to gain there;
  Fortran has the same imbalance (same slowest rank, 45). The queue runs
  (H) show C++ ahead of Fortran everywhere, so this is not needed to
  close the gap.
- `driver_kessler.cpp`'s `--help` text: fixed (2026-10-28).
- Check sub-cycle counts from real model state against the synthetic
  `--rain-frac` profile (C2).
- Compiler vectorization report for the CPU build (C4): done as item I3
  (pow/exp are vectorized).

---

## F. Standalone driver: how to run it, and results

Files:

| file                          | what it is                                                                                    |
|-------------------------------|-----------------------------------------------------------------------------------------------|
| `driver_kessler.cpp`          | test driver (checks, comparisons, timing)                                                     |
| `driver_kessler_hostloop.hpp` | the pre-D1 `kessler_run` (host-driven sub-cycle loop), `kessler_driver::kessler_run_hostloop`, for comparison |
| `driver_kessler_variants.hpp` | experimental variants of `kessler_run` (I2, I6), `kessler_driver::kessler_run_variant`         |
| `driver_kessler_preprocess.hpp` | old and new `run_impl` pre-processing (I4), for `--preprocess`                                 |
| `build_driver_kessler.sh`     | builds the driver against a case build (`--src-from-case` also copies the `.hpp`; `--cxxflags` adds flags) |

Build (from this directory) on the interactive compute node, not a
login node, with the matching case's module environment:

```bash
E=/glade/derecho/scratch/kstengel/E3SM/e3sm_test
source $E/kessler_langOp_perf/CPP_ne16np4_1day_cpu_sc1_rain/case/.env_mach_specific.sh
./build_driver_kessler.sh --case-build $E/builds/cpp_cpu_49bf6ccdde-dirty349ee6c7/build \
    --out $PWD/driver_kessler_cpu
source $E/kessler_langOp_perf/CPP_ne16np4_1day_gpu_sc1_rain/case/.env_mach_specific.sh
./build_driver_kessler.sh --case-build $E/builds/cpp_gpu_49bf6ccdde-dirty349ee6c7/build \
    --out $PWD/driver_kessler_gpu
```

Run on a GPU node (e.g. an interactive gpudev job; the CPU binary runs
on the same node's CPU):

```bash
OMP_NUM_THREADS=1 ./driver_kessler_cpu --perf --grid ne16-cpu --rain-frac 0.08
CUDA_VISIBLE_DEVICES=0 ./driver_kessler_gpu --perf --grid ne16-gpu --rain-frac 0.08
# correctness, both implementations, all checks:
./driver_kessler_cpu --grid ne16-cpu --rain-frac 0.08 --impl both --nsteps 3
# per-kernel breakdown (fences after every kernel; on GPU take
# totals from a run without it):
./driver_kessler_cpu --perf --grid ne60-cpu --rain-frac 0.08 --kernel-timers
```

New options: `--impl LIST` (comma-separated: `prod`, `hostloop`, and the I2/I6 variants
in `driver_kessler_variants.hpp`: `base`, `i2`, `i6`, `i2i6`; `both`, `all`), `--rain-frac F`, `--grid G`,
`--repeat N`, `--warmup N`, `--kernel-timers`, `--host-ref`, `--preprocess`
(compare and time the old and new pre-processing instead of `kessler_run`),
`--perf` (= `--ws managed --no-poison --no-strided --no-host-ref
--impl both --repeat 5 --warmup 1`). `--help` lists them all.

Results, 2026-10-28, deg0053 (A100 node; CPU runs single-threaded on
its AMD CPU). Seconds per `kessler_run` call, best of 10 repeats after 1
warmup, one rank's worth of columns. These were measured before the
port: "prod" is the old host loop (now `kessler_run_hostloop`) and "team"
is the D1 prototype (now `kessler_run`), with AUTO team size unless noted.
The prototype's `--team-size` option was removed with the port.

**ne16** (cpu: 108 cols, gpu: 3457 cols, 58 levels, dt 1800 s)

| rain | cpu prod | cpu team | speedup | gpu prod | gpu team | speedup |
|-----:|---------:|---------:|--------:|---------:|---------:|--------:|
| 0    | 6.08e-4  | 6.04e-4  | 1.01x   | 1.79e-4  | 1.42e-4  | 1.26x   |
| 0.03 | 4.22e-2  | 1.68e-3  | 25.1x   | 1.33e-2  | 1.10e-3  | 12.1x   |
| 0.08 | 5.43e-2  | 3.70e-3  | 14.7x   | 1.18e-2  | 1.15e-3  | 10.3x   |
| 0.12 | 5.65e-2  | 5.27e-3  | 10.7x   | 1.17e-2  | 1.32e-3  | 8.9x    |

GPU team size at 0.08: 32 -> 1.63e-3, 64 -> 1.16e-3, 128 -> 1.15e-3.

**ne60** (cpu: 1519 cols, gpu: 48601 cols, 58 levels, dt 1800 s)

| rain | cpu prod | cpu team | speedup | gpu prod | gpu team | speedup |
|-----:|---------:|---------:|--------:|---------:|---------:|--------:|
| 0    | 8.30e-3  | 8.47e-3  | 0.98x   | 1.08e-3  | 1.53e-3  | 0.71x   |
| 0.03 | 7.74e-1  | 2.35e-2  | 33.0x   | 8.92e-2  | 4.10e-3  | 21.8x   |
| 0.08 | 7.55e-1  | 5.23e-2  | 14.4x   | 8.97e-2  | 7.56e-3  | 11.9x   |
| 0.12 | 7.52e-1  | 7.63e-2  | 9.9x    | 8.98e-2  | 1.02e-2  | 8.8x    |

GPU team size (speedup vs prod in the same run):

| rain | 32    | 64    | AUTO (128) | 256  |
|-----:|------:|------:|-----------:|-----:|
| 0    | 1.23x | 1.10x | 0.71x      |      |
| 0.08 | 12.7x | 14.8x | 11.9x      | ~6x  |

---

## G. Full-model check of the port on an interactive node (2026-10-28)

Everything here was built and run inside interactive job 7991092 on
deg0053 (1 GPU node: 64 AMD cores / 128 hardware threads, 4 A100s),
reached by ssh from a login node; nothing went through the queue.
Scripts and the results table are kept in
`GK_perf_utils/workflow/interactive/`:

| file                      | what it does                                                                                 |
|---------------------------|----------------------------------------------------------------------------------------------|
| `make_deg_cases.sh`       | create the 12 cases (`run.py --no-submit`, then NTASKS and MPI_RUN_COMMAND changes); arg: `cpu`, `gpu` or both |
| `run_deg_cases.sh`        | run cases one at a time with `case.submit --no-batch` (needs `PBS_JOBID=<interactive job id>`) |
| `compare_deg_cases.py`    | timing table + bit-for-bit check                                                             |
| `results_2026-10-28.txt`  | output of `compare_deg_cases.py` for these runs                                              |

### G1. Builds

`./build.sh cpp cpu` and `./build.sh cpp gpu` (`workflow/`), run on the
node: `builds/cpp_{cpu,gpu}_7ecd5d94bb-dirty9022632f`, ~12 min (CPU) and
~21 min (GPU). HEAD `7ecd5d94bb` vs the original builds' `49bf6ccdde`
differs only in docs, the GPTL timer print width, and `derecho.cmake`
moving `--fmad=false` onto `Kokkos::kokkoscore` (same flags as before);
the code difference between the builds is the `kessler_run` port.

### G2. Cases

`run.py --impl cpp --res ne16np4 --days 1 --ic rain noRain --no-submit
--manifest GK_perf_utils/data/run_manifest_interactive.csv` (kept out of
`data/run_manifest.csv`, since these are not queue runs), with suffixes:

- `_port_r64` / `_port_r128` (cpu) and `_port` (gpu): new builds
- `_base_r64` / `_base_r128` (cpu) and `_base` (gpu): original builds,
  through a build root whose `cpp_*` links point at
  `cpp_*_49bf6ccdde-dirty349ee6c7` (`run.py --build-root ... --allow-src-drift`)

CPU r64 cases: `NTASKS=64`, `case.setup --reset`, `BUILD_COMPLETE=TRUE`.

Launch: mpibind takes ranks from the PBS select line, which the
interactive job does not have (no mpiprocs), so every case sets
`MPI_RUN_COMMAND`:

```bash
# cpu (N=64: one rank per core; N=128: ranks 64-127 on the hyperthread siblings of 0-63)
mpiexec --label --line-buffer -n N --ppn N --cpu-bind depth -d 1
# gpu
mpiexec --label --line-buffer -n 4 --ppn 4 --cpu-bind depth -d 16 set_gpu_rank
```

`set_gpu_rank` needs `PBS_NODEFILE`, which an ssh session into the job
does not set: `export PBS_NODEFILE=/var/spool/pbs/aux/<jobid> NGPUS=4`.

### G3. Results

ne16 1 day, rain = day-10 IC, sc1; s per simulated day for the run, s per
day summed over calls for the Kessler timers; "max" = slowest rank, mean
over ranks in parentheses.

| layout     | IC     | run base | run port | `kessler::run` max (mean), base | port          | `kessler_run` max, base | port   | state |
|------------|--------|---------:|---------:|--------------------------------:|--------------:|------------------------:|-------:|-------|
| gpu 4      | rain   |    3.910 |    3.222 | 0.863 (0.651)                   | 0.176 (0.147) | 0.7680                  | 0.0791 | BFB   |
| gpu 4      | noRain |    3.155 |    3.118 | 0.106 (0.105)                   | 0.102 (0.101) | 0.0088                  | 0.0052 | BFB   |
| cpu 64     | rain   |   27.427 |   21.067 | 8.354 (1.378)                   | 2.109 (0.292) | 8.2564                  | 2.0208 | BFB   |
| cpu 64     | noRain |   18.737 |   19.072 | 0.184 (0.161)                   | 0.179 (0.159) | 0.0740                  | 0.0745 | BFB   |
| cpu 128 HT | rain   |   25.682 |   19.984 | 7.305 (0.936)                   | 2.146 (0.272) | 7.2088                  | 2.0516 | BFB   |
| cpu 128 HT | noRain |   19.554 |   19.621 | 0.177 (0.161)                   | 0.176 (0.159) | 0.0692                  | 0.0714 | BFB   |

Speedups with rain, slowest rank: `kessler::run` gpu 4.9x, cpu64 4.0x,
cpu128 3.4x; `kessler_run` gpu 9.7x, cpu64 4.1x, cpu128 3.5x. Whole
model with rain: 1.2-1.3x. noRain: unchanged (within noise).
BFB = the HOMME global state lines (u, v, T, ps, M, Q1-Q4 at every
statefreq step) in `homme_atm.log` are identical between port and base.

Interpretation:

- Smaller than the standalone 10-30x because, in the model, rain is
  clustered: the slowest rank is mostly rainy columns, so the old
  loop wasted little there (`max*ncols ~ sum`). Its remaining cost is
  real work. The mean-rank Kessler time dropped 3.5-5x, and the
  slowest rank is now ~7x the mean (load imbalance).
- Reference, original queue runs (128 ranks on 128 physical cores,
  ne16 rain sc1): `kessler::run` max C++ 4.16 s, Fortran 1.38 s (both on
  rank 45); run time C++ 13.1, Fortran 12.1 s/day. Scaling the C++
  number by the 3.4-4.0x above gives ~1.0-1.2 s, i.e. at or below
  Fortran. To be confirmed by G4.
- Whole-model times here (19-27 s/day) are much higher than the
  queue runs (9-13 s/day) because the rest of the model (dycore etc.)
  runs on 64 cores, or 128 hyperthreads, of a GPU node instead of 128
  cores of a CPU node. Compare across these cases only, not with the
  queue runs.
- Single runs, no repeats; noRain differences of a few percent are
  run-to-run noise.

### G4. Next: queue runs matching the original cases

From `workflow/`, with `builds/cpp_{cpu,gpu}` now pointing at the port:

```bash
./run.py --impl cpp --hw cpu gpu --res ne16np4 --days 1 \
         --ic rain noRain --suffix port
```

(add `--subcycles 1 4`, and ne30np4/ne60np4, to cover the full matrix).
`run.py` checks that the E3SM checkout matches the build's
`git_diff.patch`. This file was edited after the
build, so it will refuse; this file is not compiled, so either pass
`--allow-src-drift` or rebuild (`./build.sh cpp cpu|gpu`, ~12-21 min) to
get a build tag that matches the checkout.

---

## H. Queue runs of the port (suffix `_updated`, 2026-10-28)

**Note (2026-10-28, evening):** the C++ matrix was resubmitted with I2 and
I4 (build `cpp_{cpu,gpu}_223c6719b2-dirtyf15739b6`) under the same suffix
`_updated`, using `workflow/cpp_batch_{11,14,51,54}.sh` (`run.py --force`).
`--force` deletes and recreates each case, so the `CPP_*_updated` case
directories now hold the I2+I4 runs, and the numbers in this section (D1
only) can no longer be regenerated from them. They are kept here and in
`GK_perf_utils/workflow/interactive/results_updated_2026-10-28.txt`. The
I2+I4 results are in section J.

The full original perf matrix (`CPP_<res>_<1|5>day_<cpu|gpu>_sc<1|4>_<rain|noRain>`)
rerun through the queue with the port build
(`cpp_{cpu,gpu}_7ecd5d94bb-dirty9022632f`), same layouts as the originals
(128 ranks on a CPU node, 4 ranks on a GPU node). Every `_updated` case
was checked to use the new build and every original the old one.
Table: `GK_perf_utils/workflow/interactive/compare_updated.py`, saved in
`results_updated_2026-10-28.txt`. All 48 have finished (the ne60 CPU
runs last, 2026-10-28 15:34).

Times per simulated day. `kessler::run` max = slowest rank; F90 max = the
Fortran case with the same name (`ne*np4_...`, no prefix), slowest rank.

| case                   | run orig | run new | speedup | `kessler::run` max orig | new   | speedup | F90 max |
|------------------------|---------:|--------:|--------:|------------------------:|------:|--------:|--------:|
| ne16 1day cpu sc1 rain |    13.09 |    9.80 | 1.34x   |                   4.163 | 1.204 | 3.5x    |   1.380 |
| ne16 1day cpu sc4 rain |    13.14 |   10.00 | 1.31x   |                   4.424 | 1.420 | 3.1x    |   1.727 |
| ne16 5day cpu sc1 rain |    13.58 |   10.01 | 1.36x   |                   3.746 | 0.883 | 4.2x    |   1.022 |
| ne16 5day cpu sc4 rain |    13.78 |    9.98 | 1.38x   |                   4.062 | 1.079 | 3.8x    |   1.335 |
| ne30 1day cpu sc1 rain |    51.81 |   41.48 | 1.25x   |                  14.469 | 3.883 | 3.7x    |   4.524 |
| ne30 1day cpu sc4 rain |    52.92 |   42.16 | 1.26x   |                  15.337 | 4.566 | 3.4x    |   5.589 |
| ne30 5day cpu sc1 rain |    53.61 |   41.42 | 1.29x   |                  13.460 | 2.847 | 4.7x    |   3.463 |
| ne30 5day cpu sc4 rain |    54.22 |   41.97 | 1.29x   |                  14.889 | 3.667 | 4.1x    |   4.586 |
| ne16 1day gpu sc1 rain |     3.88 |    3.21 | 1.21x   |                   0.858 | 0.172 | 5.0x    |   1.761 |
| ne16 1day gpu sc4 rain |     4.26 |    3.50 | 1.22x   |                   1.167 | 0.464 | 2.5x    |   2.115 |
| ne30 1day gpu sc1 rain |     9.96 |    8.30 | 1.20x   |                   1.897 | 0.231 | 8.2x    |   3.633 |
| ne30 1day gpu sc4 rain |    10.33 |    8.73 | 1.18x   |                   2.261 | 0.588 | 3.8x    |   4.192 |
| ne60 1day gpu sc1 rain |    58.10 |   52.21 | 1.11x   |                   6.362 | 0.496 | 12.8x   |  18.572 |
| ne60 1day gpu sc4 rain |    58.71 |   52.72 | 1.11x   |                   6.961 | 1.058 | 6.6x    |  20.105 |
| ne60 5day gpu sc1 rain |    57.07 |   50.65 | 1.13x   |                   6.736 | 0.486 | 13.9x   |  19.726 |
| ne60 5day gpu sc4 rain |    57.71 |   51.25 | 1.13x   |                   7.361 | 1.054 | 7.0x    |  21.287 |
| ne60 1day cpu sc1 rain |   398.58 |  353.08 | 1.13x   |                  65.510 | 19.650 | 3.3x   |  23.298 |
| ne60 1day cpu sc4 rain |   403.68 |  354.58 | 1.14x   |                  71.487 | 22.324 | 3.2x   |  27.998 |
| ne60 5day cpu sc1 rain |   404.12 |  353.07 | 1.14x   |                  65.515 | 11.550 | 5.7x   |  14.134 |
| ne60 5day cpu sc4 rain |   409.83 |  351.46 | 1.17x   |                  67.008 | 14.813 | 4.5x   |  18.847 |

noRain (all 24): `kessler::run` within 0.9-1.1x of the original, except
the ne60 5-day runs (gpu sc1 1.7x, sc4 1.2x; cpu sc1 1.9x, sc4 1.4x),
where rain develops during the 5 days. Run times within noise (one
outlier, ne30 1day cpu sc4 noRain 0.93x, with Kessler unchanged).

**ne60 CPU, steady state.** `kessler_perf_trends.md` and section B
quote the steady-state rate, (5-day total - 1-day total)/4 of
`kessler::run` on the slowest rank, which leaves out the first day; the
table above is the plain average per day. On the steady-state measure
(s per simulated day):

| ne60 cpu      | C++ orig | C++ new | Fortran |
|---------------|---------:|--------:|--------:|
| sc1, rain     |    65.52 |    9.53 |   11.84 |
| sc4, rain     |    65.89 |   12.94 |   16.56 |
| sc1, noRain*  |     3.10 |    1.48 |    3.47 |
| sc4, noRain*  |     7.96 |    5.64 |   13.93 |

(* the noRain IC develops rain after day 1, so days 2-5 are not dry.)
The case the problem was reported on, ne60 CPU with rain, went from
65.5 to 9.5 s/day (6.9x), now below Fortran's 11.8.

What this shows:

- Answers: all 48 pairs are bit-for-bit identical (HOMME global state
  statistics) to the original runs.
- CPU with rain: Kessler on the slowest rank is 3.1-5.7x faster, and
  C++ is now faster than Fortran in every CPU case (e.g. ne16 1day sc1:
  1.20 vs 1.38 s/day; ne30 1day sc1: 3.88 vs 4.52; ne60 1day sc1: 19.65
  vs 23.30). Whole model 1.13-1.38x faster (ne60: 1.13-1.17x, where the
  dycore dominates); ne16 CPU rain is now ~10 s/day vs ~9 noRain.
- GPU with rain: Kessler 5-14x faster at sc1, growing with grid size
  (more columns per GPU, so more waste removed), and 2.5-7x at sc4.
  C++ GPU was already faster than Fortran GPU; the gap is now 10-40x.
- sc4 gains are smaller than sc1: each call covers 1/4 of the physics
  step, so there are fewer sub-cycles to waste per call, and the fixed
  per-call costs inside `kessler::run` (preprocess, checks, tendencies)
  are a larger share.
- After the port the slowest CPU rank is still ~7-8x the mean
  (load imbalance from clustered rain; see E).
- Job times (case.run start to end) of the ne60 CPU runs: `_updated`
  up to 372 s (1 day) and 1783 s (5 days); the Fortran ne60 CPU cases
  up to 380 s and 1805 s. `run.py`'s `WALLCLOCK_PERF` for ne60 CPU was
  lowered from 20/80 min to 1.5x the longest job, rounded up to 5 min:
  cpp 10/45 min, f90 10/50 min (1 day/5 day), 2026-10-28.

---

## I. Review after the port: remaining costs and possible improvements

Review of `kessler_functions_impl.hpp` (ported `kessler_run`) and
`eamxx_kessler_process_interface.cpp` `run_impl`, 2026-10-28. The review
itself changed no code; I2 and I4 below were implemented afterwards, so
the timer table is from before them. Timer breakdown from the `_updated` queue runs (ne16, 1 day,
s per simulated day; mean over ranks, slowest rank in [ ]):

| timer                       | cpu sc1 rain  | cpu sc1 noRain | gpu sc1 rain  | gpu sc1 noRain |
|-----------------------------|---------------|----------------|---------------|----------------|
| `kessler::run` (total)      | 0.153 [1.204] | 0.089 [0.101]  | 0.147 [0.172] | 0.102 [0.103]  |
| `kessler_run` (the kernel)  | 0.098 [1.148] | 0.034 [0.039]  | 0.051 [0.078] | 0.005 [0.005]  |
| run-postcondition-checks    | 0.027         | 0.026          | 0.061         | 0.061          |
| run::preprocess             | 0.021         | 0.021          | 0.012         | 0.012          |
| run-precondition-checks     | 0.002         | 0.002          | 0.007         | 0.007          |
| 4 small update kernels      | 0.003         | 0.003          | 0.003         | 0.003          |

(ne30 gpu sc4 rain: total 0.523, kernel 0.104, postcondition checks
0.269.)

So on GPU, and on CPU without rain, the kernel is no longer the main
cost: most of `kessler::run` is EAMxx's field checks and preprocessing.
On CPU with rain, the remaining cost is real sub-cycling work on the
rainiest rank.

Possible improvements, roughly in order of payoff. "BFB" = should keep
answers bit-for-bit identical, so the driver and G/H checks still apply.

- [ ] **I1. Postcondition/precondition checks** (EAMxx infrastructure, not
  Kessler code). The FieldWithinInterval/LowerBound checks on
  T_mid, qv, qc, qr, qi and precl (`initialize_impl`) each do a
  reduction and a host sync every call. Now the largest single
  cost on GPU: 60% of `kessler::run` without rain, ~50% at ne30 sc4;
  ~30% on CPU without rain. The Fortran cases pay the same (0.059
  s/day on GPU), so it does not change C++ vs Fortran, but it caps
  what any Kessler change can show.
  - Runtime switch: `atmchange kessler::enable_postcondition_checks=false`
    (and `enable_precondition_checks`). Only as a labelled
    experiment, applied to C++, Fortran and JAX alike.
  - Whether the checks are slower than necessary (~0.2 ms per
    check per call on GPU at ne16) is a question for EAMxx's
    property-check code.

- [x] **I2. Skip pow on all-dry packs (BFB).** Each pass does 4 pow and 1 exp
  per level: `pow(qr, 0.875)` in qrprod, `pow(rqr, 0.2046)` and
  `pow(rqr, 0.525)` in the evaporation term A, and `pow(qr*r, 0.1364)`
  for velqr. Where qr == 0 these are exactly 0 (`pow(0, x>0) = +0`),
  so a pack whose qr lanes are all 0 can use 0 without calling
  pow. Check qr at the point each pow uses it: qrprod uses qr
  before the update; A and velqr use qr after qrprod and sed are
  added. The exp in qvs is needed everywhere. On CPU (Pack 16)
  this skips the pows above the rain in rainy columns; on GPU
  (Pack 1) it helps where a whole warp of levels is dry. Targets
  the 1.15 s/day slowest-rank kernel time; rough guess 20-40% off,
  depending on how deep the rain is. Validate with
  `driver_kessler.cpp --impl both` (bit for bit vs hostloop).

  **Ported 2026-10-28 (uncommitted; checked in the model with I4, below).** Production
  `kessler_run` now calls `detail::pow_nonneg` (in
  `kessler_functions_impl.hpp`) for the four rain powers. After the port
  the driver checks all pass and production still matches `hostloop` and
  every variant bit for bit (CPU and GPU, 0/8/50% rain, 3 steps).
  Production vs production-before-I2 (`--impl prod,base`), 8% rain: CPU
  1.30x (ne16), 1.34x (ne60); GPU 1.11x (ne16), ~1.0x (ne60, noisy run:
  1.18x in the variant test below).

  **Model check 2026-10-28 (with I4, interactive node deg0022):** build
  `cpp_{cpu,gpu}_223c6719b2-dirtyf15739b6` (I2 + I4) vs
  `cpp_{cpu,gpu}_7ecd5d94bb-dirty9022632f` (before both), ne16 1 day, CPU 64
  ranks and GPU 4 ranks, rain and noRain (cases `*_i2i4*` vs `*_pre*`;
  `GK_perf_utils/workflow/interactive/results_i2i4_2026-10-28.txt`). HOMME
  global state statistics bit-for-bit in all four pairs. Slowest rank,
  s/day, old -> new:

  | ne16                     | rain                    | noRain                  |
  |--------------------------|-------------------------|-------------------------|
  | cpu `kessler_run` (I2)   | 1.953 -> 1.718 (1.14x)  | 0.080 -> 0.056 (1.43x)  |
  | cpu `preprocess` (I4)    | 0.039 -> 0.032 (1.22x)  | 0.039 -> 0.033 (1.18x)  |
  | cpu `kessler::run`       | 2.041 -> 1.796 (1.14x)  | 0.179 -> 0.170 (1.05x)  |
  | gpu `kessler_run` (I2)   | 0.078 -> 0.075 (1.04x)  | 0.0052 -> 0.0043 (1.2x) |
  | gpu `preprocess` (I4)    | 0.016 -> 0.006 (2.62x)  | 0.011 -> 0.006 (1.83x)  |
  | gpu `kessler::run`       | 0.172 -> 0.163 (1.06x)  | 0.101 -> 0.095 (1.06x)  |

  I2 gains less in the model than standalone (1.14x vs 1.3x on CPU with
  rain): the slowest rank is the rainiest, so fewer of its packs are dry.
  Over all ranks (mean) CPU `kessler::run` with rain went 0.291 -> 0.234
  (1.24x). Whole-model run times did not change beyond noise.

  Variant test before the port (standalone driver, deg0022):
  `kessler_run_variant<KMF, SkipDryPow=true, ...>` in
  `driver_kessler_variants.hpp` (`--impl i2`). Bit-for-bit with production
  on CPU and GPU (0, 8%, 50% rain, 3 steps, all driver checks). Speedup vs
  production, s per call, best of 10:

  | grid / rain    | cpu i2 | gpu i2 |
  |----------------|-------:|-------:|
  | ne16, 0        | 1.64x  | 1.12x  |
  | ne16, 0.08     | 1.34x  | 1.08x  |
  | ne16, 0.12     | 1.33x  | 1.14x  |
  | ne60, 0        | 1.60x  | 1.22x  |
  | ne60, 0.08     | 1.38x  | 1.18x  |
  | ne60, 0.12     | 1.31x  | 1.19x  |

  (The synthetic no-rain profile has `qr = 0` everywhere, so every pow is
  skipped; the `base` variant, the same code with both switches off, ran
  at 0.94-1.06x of production, which bounds the noise.) A small change,
  bit-for-bit, ~25-30% off the CPU kernel with rain, so it was ported.

- [x] **I3. CPU vectorization report** (was C4). The micro math is
  the remaining CPU cost, and it is not known whether nvc++
  vectorizes `ekat::pow`/`exp` over 16-wide packs under `-Mnofma
  -Kieee`. Compile the driver with `-Minfo=vect` (no code change) and
  check the `kessler_run` micro loop. If pow/exp do not vectorize,
  this could be worth more than I2.

  **Done 2026-10-28: pow and exp are vectorized; nothing to gain.**
  `build_driver_kessler.sh --cxxflags "-Minfo=vect,inline"` gave the report
  (`GK_perf_utils/workflow/interactive/variants_2026-10-28/vect_report.log`),
  but its line numbers are ambiguous across inlined headers, so the CPU
  binary was checked directly (`nm`/`objdump`). In production
  `kessler_run` every `ekat::exp` is a call to `__pd_exp_4` (nvc++'s vector
  exp, 4 doubles, AVX2: the widest the AMD Milan CPUs support), and
  `ekat::pow` calls `__pd_pow_4` with a scalar `__pd_pow_1` fallback that
  only runs if nvc++'s run-time aliasing check fails. A 16-wide pack takes
  4 vector calls. `ekat::pow` is called out of line rather than inlined,
  which is cheap next to 16 pows.

- [x] **I4. Preprocess is unpacked scalar work** (BFB if the formulas stay the
  same). `run_impl`'s `Kessler_preprocess` indexes
  `(icol, klev / Pack::n)[klev % Pack::n]` for every level in a
  RangePolicy, so it does not use the packs. 0.021 s/day on CPU,
  ~60% of the kernel's own time without rain. Write it over packs
  (`PF::exner_function` etc. have Pack overloads), or fold it, and
  the z_int/z_mid scan, into the per-column kernel.

  **Tested 2026-10-28 (standalone driver, deg0022), before it went into
  `run_impl`: bit-for-bit, ~1.85x on CPU.** `driver_kessler.cpp --preprocess` runs the
  current pre-processing (`preprocess_old`, copied from `run_impl`) and the
  new one (`preprocess_new`), both in `driver_kessler_preprocess.hpp`, on
  synthetic columns with `qv` as a strided tracer slice, and compares `pk`,
  `theta`, `dz`, `rho`, `z_mid`, `z_int` on the real levels. The new
  version is one team per column (scan policy, as SHOC), works on whole
  packs, guards padding lanes, and computes `theta = T / pk` from the Exner
  value it just made: `calculate_theta_from_T(T, p)` is `T /
  exner_function(p)`, so this is the same arithmetic with one `pow` per
  level instead of two. Results, all six outputs identical on CPU (Pack
  16), GPU and host, ne16 and ne60:

  | per call   | old       | new       | old/new |
  |------------|----------:|----------:|--------:|
  | CPU ne16   | 2.75e-4 s | 1.47e-4 s | 1.87x   |
  | CPU ne60   | 3.90e-3 s | 2.11e-3 s | 1.85x   |
  | GPU ne16   | 4.29e-5 s | 3.64e-5 s | 1.18x   |
  | GPU ne60   | 3.17e-4 s | 2.15e-4 s | 1.48x   |

  Without the Exner reuse (packing and fusing only) the CPU gain was just
  1.09-1.12x: the cost is the pows, not the indexing. In the model this is
  roughly 0.006 s/day per rank off `kessler::run` at ne16 CPU sc1 (about
  7% without rain), 4x that at sc4, plus one kernel launch fewer per call
  on GPU. The `preprocess` block is the same in the bridge and JAX
  branches (variable names differ), so the change carries over; apply it
  to all three to keep the comparison fair.

  **Put into `run_impl` 2026-10-28 (uncommitted) and checked in the
  model** together with I2: bit-for-bit, `preprocess` 1.2x faster on CPU
  and 1.8-2.6x on GPU at ne16 (table under I2). On GPU, read that number
  with care: in these runs `run_impl` stopped the `preprocess` timer
  without a `Kokkos::fence()` (as it did before I4), so the timer mostly
  measured kernel launches, and the kernel itself was charged to the next
  fenced timer (`kessler_update_timestep_init`). The gain shown is mostly
  the launch saved. **Fixed 2026-10-29 (uncommitted, not yet in any
  build):** `run_impl` now fences before stopping the `preprocess` timer,
  like every other sub-timer. Runs from then on report the real GPU
  `preprocess` time, and their `preprocess` and
  `kessler_update_timestep_init` sub-timers are not comparable with
  earlier runs (totals are). CPU (Serial) numbers are unaffected. Queue runs of the full
  C++ matrix with I2 + I4 run 2026-10-28 (suffix `_updated`): all 48
  bit-for-bit, results in section J. Still to do: the same change in the bridge and
  JAX branches (same two kernels; variables `params_helpers.*` and `*_v`;
  keep `theta = T / exner` a single division by the same Exner value),
  including the `Kokkos::fence()` before `stop_timer(...preprocess)`, so
  the three implementations' sub-timers measure the same thing.

- [ ] **I5. Fuse the small update kernels** (BFB if `T_mid = T_mid_prev +
  T_mid_tend*dt` stays as written). `kessler_update_timestep_init`,
  `kessler_update_run`, `kessler_apply_T_tend` and
  `kessler_update_timestep_final` are separate sweeps, each
  followed by a `Kokkos::fence`, 0.003 s/day together. Folding them
  into the per-column kernel saves launches and syncs, mainly on
  GPU. Small gain.

- [ ] **I6. GPU memory traffic and Workspace size.**
  - `sed` and `velqr` could move to team scratch (2 x nlev_packs
    Packs, ~1 KB per team on GPU) instead of global memory; `r`,
    `rhalf` and `pc` could be computed on the fly.
  - Shrink the Workspace (D6): drop `f5`; `dt0`, `time_counter` and
    `precl_acc` are now only written once at the end. Saves memory,
    not time.

  **Team scratch tested 2026-10-28 (standalone driver): no gain.**
  `kessler_run_variant<KMF, ..., TeamScratch=true>` (`--impl i6`) keeps
  `r`, `rhalf`, `pc`, `velqr` and `sed` in Kokkos level-0 team scratch
  (~2.3 KB per team; shared memory on GPU). Bit-for-bit with production.
  Speedup vs production: CPU 0.98-1.05x, GPU 0.96-1.02x; combined with
  I2 (`--impl i2i6`) CPU matches I2 alone and GPU is 0-7% better than I2
  alone (1.12-1.28x), within the run-to-run spread. Each column's arrays
  already stay in cache (CPU) or L1/L2 (GPU), so moving them gains
  nothing measurable. Not worth porting for speed; the Workspace shrink
  (memory only) can still be done separately.

- [ ] **I7. Load imbalance.** The rainiest rank does 7-8x the mean Kessler
  work (H). That is domain decomposition, not Kessler code;
  Fortran has the same imbalance.

- [ ] **I8. Per-call logging.** `run_impl` writes two info-level lines per call
  (96 per day) and the perf runs set `atm_log_level=info`. Use debug
  level for these, or run perf cases at a higher log level.

Not BFB (would change answers; need a validation decision first; not
worth doing before I1-I4):

- `pc = 3.8 / (pk^xk * pref)` is mathematically `3.8 * p0 / (p * pref)`,
  which drops one pow per level per call.
- `pow(rqr, 0.2046)` and `pow(rqr, 0.525)` could share one `log(rqr)`
  (`exp(a*L)` each), fewer transcendental calls per pass.

Checked, no problem found: the per-column sub-cycle state, barriers
(one after sed, one at the end of each pass; relhum runs after the last
pass's barrier), the error codes and their single host check, and the
`max_pass` limit (same 100000 as before).

---

## J. Queue runs with I2 + I4 (suffix `_updated`, rerun 2026-10-28)

The full C++ perf matrix rerun with I2 (skip `pow` on dry packs) and I4
(packed pre-processing) on top of D1: build
`cpp_{cpu,gpu}_223c6719b2-dirtyf15739b6`, same 48 cases, layouts and
suffix as H (`workflow/cpp_batch_{11,14,51,54}.sh`, `run.py --force`), so
the `CPP_*_updated` case directories now hold these runs. All finished
2026-10-28. Files in `GK_perf_utils/workflow/interactive/`:

- `results_updated_i2i4_2026-10-29.txt`: `compare_updated.py` output,
  these runs vs the original C++ cases and Fortran.
- `results_d1_vs_i2i4_2026-10-29.txt`: these runs vs the D1-only runs
  (`results_updated_2026-10-28.txt`), plus the `kessler_run` and
  `preprocess` sub-timers.

**Answers:** all 48 cases are bit-for-bit with the original C++ runs
(HOMME global state statistics), so also with the D1-only runs.

**`kessler::run`, slowest rank, s per simulated day, with rain** ("orig"
= the original host-loop C++ cases, "D1" = section H, "vs" columns are
speedups, F90 / C++ = how much faster C++ is than Fortran):

| case                   |   orig |     D1 | D1 + I2 + I4 | vs D1 | vs orig |    F90 | F90 / C++ |
|------------------------|-------:|-------:|-------------:|------:|--------:|-------:|----------:|
| ne16 1day cpu sc1 rain |  4.163 |  1.204 |        1.049 | 1.15x |    4.0x |  1.380 |      1.3x |
| ne16 1day cpu sc4 rain |  4.424 |  1.420 |        1.191 | 1.19x |    3.7x |  1.727 |      1.5x |
| ne16 5day cpu sc1 rain |  3.746 |  0.883 |        0.685 | 1.29x |    5.5x |  1.022 |      1.5x |
| ne16 5day cpu sc4 rain |  4.062 |  1.079 |        0.899 | 1.20x |    4.5x |  1.335 |      1.5x |
| ne30 1day cpu sc1 rain | 14.469 |  3.883 |        3.189 | 1.22x |    4.5x |  4.524 |      1.4x |
| ne30 1day cpu sc4 rain | 15.337 |  4.566 |        3.754 | 1.22x |    4.1x |  5.589 |      1.5x |
| ne30 5day cpu sc1 rain | 13.460 |  2.847 |        2.431 | 1.17x |    5.5x |  3.463 |      1.4x |
| ne30 5day cpu sc4 rain | 14.889 |  3.667 |        2.967 | 1.24x |    5.0x |  4.586 |      1.5x |
| ne60 1day cpu sc1 rain | 65.510 | 19.650 |       16.663 | 1.18x |    3.9x | 23.298 |      1.4x |
| ne60 1day cpu sc4 rain | 71.487 | 22.324 |       17.746 | 1.26x |    4.0x | 27.998 |      1.6x |
| ne60 5day cpu sc1 rain | 65.515 | 11.550 |        9.237 | 1.25x |    7.1x | 14.134 |      1.5x |
| ne60 5day cpu sc4 rain | 67.008 | 14.813 |       12.127 | 1.22x |    5.5x | 18.847 |      1.6x |
| ne16 1day gpu sc1 rain |  0.858 |  0.172 |        0.163 | 1.06x |    5.3x |  1.761 |     10.8x |
| ne16 1day gpu sc4 rain |  1.167 |  0.464 |        0.445 | 1.04x |    2.6x |  2.115 |      4.8x |
| ne16 5day gpu sc1 rain |  0.890 |  0.175 |        0.162 | 1.08x |    5.5x |  1.801 |     11.1x |
| ne16 5day gpu sc4 rain |  1.222 |  0.472 |        0.432 | 1.09x |    2.8x |  2.195 |      5.1x |
| ne30 1day gpu sc1 rain |  1.897 |  0.231 |        0.215 | 1.07x |    8.8x |  3.633 |     16.9x |
| ne30 1day gpu sc4 rain |  2.261 |  0.588 |        0.555 | 1.06x |    4.1x |  4.192 |      7.6x |
| ne30 5day gpu sc1 rain |  1.959 |  0.231 |        0.214 | 1.08x |    9.2x |  3.732 |     17.4x |
| ne30 5day gpu sc4 rain |  2.343 |  0.573 |        0.536 | 1.07x |    4.4x |  4.315 |      8.1x |
| ne60 1day gpu sc1 rain |  6.362 |  0.496 |        0.456 | 1.09x |   14.0x | 18.572 |     40.7x |
| ne60 1day gpu sc4 rain |  6.961 |  1.058 |        0.989 | 1.07x |    7.0x | 20.105 |     20.3x |
| ne60 5day gpu sc1 rain |  6.736 |  0.486 |        0.454 | 1.07x |   14.8x | 19.726 |     43.4x |
| ne60 5day gpu sc4 rain |  7.361 |  1.054 |        0.974 | 1.08x |    7.6x | 21.287 |     21.9x |

Ranges over the 12 cases of each kind, I2 + I4 vs D1 only:

|             | `kessler::run` slowest rank | mean over ranks | model run time | C++ vs Fortran (slowest rank) |
|-------------|----------------------------:|----------------:|---------------:|------------------------------:|
| cpu rain    | 1.15-1.29x                  | 1.17-1.29x      | 0.98-1.06x     | 1.3-1.6x faster               |
| cpu noRain  | 1.12-1.33x                  | 1.15-1.31x      | 0.99-1.08x     | 1.9-2.9x faster               |
| gpu rain    | 1.04-1.09x                  | 1.06-1.08x      | 0.99-1.01x     | 4.8-43x faster                |
| gpu noRain  | 1.03-1.14x                  | 1.05-1.13x      | 1.00-1.01x     | 1.4-5.2x faster               |

**ne60 CPU, steady state** ((5-day - 1-day)/4, as in H), s per
simulated day:

| ne60 cpu      | C++ orig |    D1 | D1 + I2 + I4 | Fortran |
|---------------|---------:|------:|-------------:|--------:|
| sc1, rain     |    65.52 |  9.53 |         7.38 |   11.84 |
| sc4, rain     |    65.89 | 12.94 |        10.72 |   16.56 |
| sc1, noRain*  |     3.10 |  1.48 |         1.26 |    3.47 |
| sc4, noRain*  |     7.96 |  5.64 |         4.79 |   13.93 |

(* the noRain IC develops rain after day 1.) The case the problem was
reported on, ne60 CPU with rain, is now 7.4 s/day against 65.5 originally
(8.9x) and Fortran's 11.8 (1.6x faster than Fortran).

**Sub-timers, ne16 1 day sc1** (D1 numbers from the table at the top of
section I; `kessler_run` = slowest rank, `preprocess` = mean over ranks):

|                        | cpu rain       | cpu noRain     | gpu rain       | gpu noRain     |
|------------------------|----------------|----------------|----------------|----------------|
| `kessler_run` (I2)     | 1.148 -> 1.002 | 0.039 -> 0.026 | 0.078 -> 0.074 | 0.005 -> 0.004 |
| `preprocess` (I4)      | 0.021 -> 0.015 | 0.021 -> 0.015 | 0.012 -> 0.006 | 0.012 -> 0.006 |

What this shows:

- The queue runs agree with the ne16 interactive check (I2, I4): CPU
  `kessler::run` is 1.15-1.33x faster than D1 alone, rain or not; GPU
  1.03-1.14x. The CPU gain is larger at sc4 without rain (1.3x), where
  `preprocess` runs 4 times per physics step and `kessler_run` itself is
  cheap. (In these runs the GPU `preprocess` number mostly measures kernel launches; see
  the timer note under I4.)
- Against the original C++ build, CPU with rain is now 3.7-7.1x faster
  and GPU with rain 2.6-14.8x; against Fortran, C++ is faster in all 48
  cases (CPU with rain 1.3-1.6x).
- Whole-model run times are unchanged within noise: Kessler is a small
  part of the step once D1 is in.
- The remaining CPU cost with rain is the sub-cycling on the rainiest
  rank (I7); on GPU it is mostly the property checks (I1).
