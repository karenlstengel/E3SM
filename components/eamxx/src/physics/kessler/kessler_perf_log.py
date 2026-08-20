"""
Lightweight in-memory call-timer for the kessler JAX path, flushed to a
summary CSV once per rank at simulation finalize.

GPTL (eamxx_timing.txt) is only reachable from C++ -- there is no pybind11
binding exposing scream::start_timer/stop_timer to the embedded Python
interpreter -- so JAX-side timings are logged here instead. Import is
optional everywhere it's used (`try/except ImportError`) so removing/
renaming this file never breaks the kessler run path.

log_call() only ever touches in-process memory (protected by a
threading.Lock, since EAMxx's embedded Python interpreter can be called
from multiple host threads) -- there is no file I/O per call. Previously
log_call() opened the CSV and took a cross-process flock() on every call;
with hundreds of calls/rank that flock() serialized every MPI rank against
every other rank on every timestep, which measurably slowed the run down.

flush_log() does the file I/O instead, and is meant to be called exactly
once per rank, at simulation finalize (see kessler.py's finalize(), wired
to KesslerMicrophysics::finalize_impl() in the C++ process interface) -- so
the one-time flock() there costs nothing. If a rank never reaches finalize
(crash, MPI_Abort), that rank's timings are simply lost -- an accepted
tradeoff for a dev/perf-instrumentation tool, not a science output.

flush_log() appends this rank's per-label totals to a small raw
intermediate CSV (KESSLER_PERF_LOG_PATH with a ".raw" suffix), then
recomputes the GPTL-style summary (count/walltotal/wallmax/wallmin per
label, across all ranks that have flushed so far) and rewrites
KESSLER_PERF_LOG_PATH itself with that summary -- so the CSV at
KESSLER_PERF_LOG_PATH is always the final summary directly; no separate
post-processing script needed. Because every rank's append +
summary-rewrite happens inside one flock() critical section, whichever
rank's critical section runs last (in real time) is guaranteed to see
every other rank's already-committed row, so the summary is exactly
correct once every rank has finalized, regardless of what order ranks
arrive in.

Set KESSLER_PERF_LOG_PATH to override the default CSV location (which
otherwise sits next to this file, in the source tree).
"""
import csv
import fcntl
import os
import threading

_DEFAULT_LOG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "kessler_perf_log.csv")
_LOG_PATH = os.environ.get("KESSLER_PERF_LOG_PATH", _DEFAULT_LOG_PATH)
_LOG_ROOT, _LOG_EXT = os.path.splitext(_LOG_PATH)
_RAW_PATH = _LOG_ROOT + ".raw" + (_LOG_EXT or ".csv")

_RAW_HEADER = ["label", "rank", "count", "walltotal_s", "callmax_s", "callmin_s"]
_SUMMARY_HEADER = ["label", "count", "walltotal_s", "wallmax_s", "wallmax_rank",
                    "wallmin_s", "wallmin_rank", "callmax_s", "callmax_rank"]

_lock = threading.Lock()
# label -> [count, walltotal_s, callmax_s, callmin_s]
_totals = {}

# Rank env vars checked in order, covering the common MPI launchers
# (OpenMPI, MPICH/Intel MPI, Slurm srun).
_RANK_ENV_VARS = ("OMPI_COMM_WORLD_RANK", "PMI_RANK", "SLURM_PROCID", "MPI_LOCALRANKID", "PALS_RANKID")


def _detect_rank():
    for var in _RANK_ENV_VARS:
        val = os.environ.get(var)
        if val is not None:
            try:
                return int(val)
            except ValueError:
                pass
    return 0


_RANK = _detect_rank()


def log_call(label, ncol, nz, dt, elapsed):
    full_label = "a:EAMxx::kessler::run::py_run::" + label
    with _lock:
        entry = _totals.setdefault(full_label, [0, 0.0, elapsed, elapsed])
        entry[0] += 1
        entry[1] += elapsed
        if elapsed > entry[2]:
            entry[2] = elapsed
        if elapsed < entry[3]:
            entry[3] = elapsed


def flush_log():
    with _lock:
        totals = {label: tuple(v) for label, v in _totals.items()}
    if not totals:
        return

    with open(_RAW_PATH, "a+", newline="") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)  # one-time cost: only taken once per rank, at finalize
        try:
            writer = csv.writer(f)
            if os.fstat(f.fileno()).st_size == 0:
                writer.writerow(_RAW_HEADER)
            for label, (count, walltotal, callmax, callmin) in totals.items():
                writer.writerow([label, _RANK, count, walltotal, callmax, callmin])
            f.flush()
            os.fsync(f.fileno())

            f.seek(0)
            rows = list(csv.reader(f))[1:]  # skip header

            summary = {}
            for row_label, row_rank, row_count, row_total, row_cmax, row_cmin in rows:
                row_rank = int(row_rank)
                row_count = int(row_count)
                row_total = float(row_total)
                row_cmax = float(row_cmax)
                s = summary.setdefault(row_label, {
                    "count": 0, "walltotal": 0.0,
                    "wallmax": (None, float("-inf")), "wallmin": (None, float("inf")),
                    "callmax": (None, float("-inf")),
                })
                s["count"] += row_count
                s["walltotal"] += row_total
                if row_total > s["wallmax"][1]:
                    s["wallmax"] = (row_rank, row_total)
                if row_total < s["wallmin"][1]:
                    s["wallmin"] = (row_rank, row_total)
                if row_cmax > s["callmax"][1]:
                    s["callmax"] = (row_rank, row_cmax)

            with open(_LOG_PATH, "w", newline="") as sf:
                swriter = csv.writer(sf)
                swriter.writerow(_SUMMARY_HEADER)
                for label, s in summary.items():
                    swriter.writerow([
                        label, s["count"], s["walltotal"],
                        s["wallmax"][1], s["wallmax"][0],
                        s["wallmin"][1], s["wallmin"][0],
                        s["callmax"][1], s["callmax"][0],
                    ])
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)
