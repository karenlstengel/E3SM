import os
import time
import numpy as np
os.environ["JAX_ENABLE_X64"] = "1"     # MUST precede any jax import
os.environ["JAX_PLATFORMS"] = os.environ.get("JAX_PLATFORMS", "cpu")  # MUST precede any jax import; honors the driver script's export, else forces cpu
import jax
import jax.numpy as jnp
jax.config.update("jax_enable_x64", True)
# recommended: persistent compilation cache 
_cache_dir = os.environ.get(
    "JAX_COMPILATION_CACHE_DIR",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),".jax_cache"),
)
jax.config.update("jax_compilation_cache_dir", _cache_dir)
jax.config.update("jax_persistent_cache_min_compile_time_secs", 0)
jax.config.update("jax_persistent_cache_min_entry_size_bytes", -1)

import sys
sys.path.insert(0, "/glade/derecho/scratch/kstengel/E3SM/E3SM/components/eamxx/src/physics/kessler/kessler_jax") 

# kessler_update is a low-arithmetic-intensity scheme: the host-facing
# bridge (contract 1) is dominated by the PCIe transfer -- measured 0.33x
# vs. serial Fortran, worse than not porting it at all. The device-resident
# bridge (contract 2) chains all three per-step phases on the GPU with no
# transfers between them and is the pattern this scheme strictly requires
# (measured 181x for the full 3-phase step). See USING_THE_BRIDGE.md and
# BRIDGE_CONTRACT2_notes.md in the laft-kessler-update project.
from bridge.kessler_update_init_bridge import kessler_update_init_bridge
from bridge.kessler_update_timestep_init_bridge import kessler_update_timestep_init_bridge_device
from bridge.kessler_update_run_bridge import kessler_update_run_bridge_device, to_device
from bridge.kessler_update_timestep_final_bridge import kessler_update_timestep_final_bridge_device

# kessler_jax.kessler_run_core expects arrays already in JAX (nz, ncol)
# layout -- kessler_run_bridge is what does the (ncol, nz) <-> (nz, ncol)
# layout conversion (in-jit, PATH C: pure H2D/D2H with the axis reversal
# fused into the jitted kernel by XLA), so EAMxx's native (ncol, nz)
# arrays go through the bridge, not kessler_jax.kessler_run directly.
from bridge.kessler_init_bridge import kessler_init_bridge
from bridge.kessler_run_bridge import kessler_run_bridge

try:
    from kessler_perf_log import log_call as _log_perf_call, flush_log as _flush_perf_log
except ImportError:
    _log_perf_call = None
    _flush_perf_log = None

# Setup a few global variables that we set values for with init() and then use in run()
latvap = None
pref = None
rhoqr = None
gravit = None 
errmsg = ""
errflg = 0
scheme_name = "kessler"

# latvap, P0, rhoqr, gravit
def init(lv_in, pref_in, rhoqr_in, gravit_in):
    global latvap, pref, rhoqr, gravit
    global errmsg, errflg

    errmsg, errflg, latvap, pref, rhoqr = kessler_init_bridge(lv_in, pref_in, rhoqr_in, errmsg, errflg, latvap, pref, rhoqr)

    # kessler_update carries its own MODULE variable `gravit` (INOUT,
    # threaded from here through every kessler_update_timestep_final call
    # in update() below) -- separate from kessler's own MODULE vars
    # (latvap/pref/rhoqr) above.
    errmsg, errflg, gravit = kessler_update_init_bridge(gravit_in, "", 0, gravit)

# calls the kessler_run function from https://github.com/NCAR/llm-fortran-modernization/tree/main/fortran2jax-kessler/_officialJAX
# these arrays should be automatically updated back in EAMxx if everything is setup correctly.

def run(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, zm, pk, theta, qv, qc, qr, precl, relhum):

    global latvap, pref, rhoqr
    global errmsg, errflg, scheme_name

    # Compute — kessler_run_bridge (contract 1, host-facing) handles the
    # (ncol, nz) <-> (nz, ncol) layout conversion in-jit and does its own
    # H2D/D2H; called directly on EAMxx's native (ncol, nz) arrays.
    #
    # Timed the same way as fortran_bridge/kessler_eamxx_bridge_main.F90's
    # system_clock wrap around its `kessler_run` call, and logged under the
    # same label, so py_run and F90_run rows line up in the perf CSV.
    _t0_run = time.perf_counter()
    theta_out, qv_out, qc_out, qr_out, precl_out, relhum_out, scheme_name, errmsg, errflg, latvap, pref, rhoqr = kessler_run_bridge(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, zm, pk, theta, qv, qc, qr, precl, relhum, scheme_name, errmsg, errflg, latvap, pref, rhoqr)
    if _log_perf_call is not None:
        try:
            _log_perf_call("kessler_run", ncol, nz, dt, time.perf_counter() - _t0_run)
        except Exception:
            pass

    # theta/qv/qc/qr/precl/relhum are zero-copy views into EAMxx's field
    # buffers (created in create_py_field()); the C++ caller (py_module_call)
    # discards this function's return value, so results must be written
    # back in place for them to reach EAMxx at all.
    #
    # KNOWN PERF ISSUE (tracked, not yet fixed): this write-back block
    # ("writeback" in kessler_perf_log.csv) costs about as much as the
    # kessler_run JAX call itself on the GPU production run (~0.0087s/call
    # vs ~0.0099s/call, over 576 calls/rank) -- a real, measured chunk of
    # a:EAMxx::kessler::run's total that isn't the JAX computation at all.
    #
    # Ruled out: GPU-resident destination memory. Re-ran the identical GPU
    # case with py_backend=host (theta/qv/qc/qr/precl/relhum become plain
    # host arrays instead of pybind11 views over Kokkos device memory) --
    # writeback cost was unchanged (5.024s vs 5.026s total over the run),
    # and so was a:EAMxx::kessler::run's overall total (87.768s vs 87.001s).
    # So this isn't about host vs. device memory on either side.
    #
    # Leading hypothesis instead: fixed per-call overhead from doing SIX
    # separate JAX-array -> NumPy conversions/assignments, rather than a
    # data-volume/bandwidth cost. ~28MB across the six arrays should take
    # ~2ms at ordinary memory bandwidth; we're seeing ~8.7ms/call (~1.45ms
    # per array), which is in the range of per-call dispatch overhead
    # (JAX array materialization, buffer-protocol negotiation), not copy
    # time -- and that would explain why it doesn't move with py_backend.
    #
    # Proposed fix (not yet implemented -- holding off for now): theta, qv,
    # qc, qr, and relhum are all the same (ncol, nz) shape. Stack them into
    # one array before returning from kessler_run_core (e.g.
    # jnp.stack([...], axis=0)), materialize that ONE stacked array here
    # instead of five separate ones, then slice it apart for the five
    # assignments below. That cuts six separate JAX->NumPy crossings down
    # to two (the stacked block + precl), which should cut this cost
    # roughly proportionally if the fixed-per-call-overhead hypothesis is
    # right. Would need a matching change in kessler_jax/kessler_run.py's
    # return values (and kessler_run's unpacking) to un-stack before this
    # runs -- not just a change here.

    _t0_writeback = time.perf_counter()
    theta[...]  = theta_out
    qv[...]     = qv_out
    qc[...]     = qc_out
    qr[...]     = qr_out
    precl[...]  = precl_out
    relhum[...] = relhum_out
    if _log_perf_call is not None:
        try:
            _log_perf_call("writeback", ncol, nz, dt, time.perf_counter() - _t0_writeback)
        except Exception:
            pass

def update(ncol, nz, dt, cpair, zm, pk, theta, phis, temp, temp_prev, temp_tend, st_energy):
    global gravit
    global errmsg, errflg

    # Contract 2 (device-resident): upload each array ONCE (to_device is a
    # pure H2D copy -- no host-side permute; the bridge_device wrappers
    # reverse axes inside their jitted region), chain timestep_init -> run ->
    # timestep_final entirely on the GPU with no transfers between phases,
    # then fetch ONCE at the end. Per-phase errflg is left on the device
    # (contract 2's own convention) and fetched together with the arrays,
    # not checked per call.
    temp_d  = to_device(temp)
    theta_d = to_device(theta)
    exner_d = to_device(pk)
    cpair_d = to_device(cpair)
    zm_d    = to_device(zm)
    phis_d  = to_device(phis)
    zeros_d = jnp.zeros((ncol, nz), dtype=jnp.float64)

    # Per-phase timing below mirrors fortran_bridge/kessler_eamxx_bridge_update.F90's
    # system_clock wrap around kessler_update_timestep_init/_run/_timestep_final
    # (same three labels, so py_run and F90_run rows line up in the perf CSV).
    # Unlike the Fortran calls, these three are intentionally left un-synced
    # on the device -- contract 2's whole point is to chain them with no D2H
    # between phases and fetch once at the end (see the batched device_get
    # below) -- so perf_counter() here times host-side dispatch, not device
    # execution. A true per-phase device time would need jax.block_until_ready()
    # after each call, which would force a sync point between phases and
    # defeat that chaining.
    _t0 = time.perf_counter()
    temp_prev_d, ttend_t_d, e1 = kessler_update_timestep_init_bridge_device(
        temp=temp_d, temp_prev=zeros_d, ttend_t=zeros_d, errflg=0)
    if _log_perf_call is not None:
        try:
            _log_perf_call("kessler_update_timestep_init", ncol, nz, dt, time.perf_counter() - _t0)
        except Exception:
            pass

    _t0 = time.perf_counter()
    ttend_t_d, e2 = kessler_update_run_bridge_device(
        nz=nz, ncol=ncol, dt=dt, theta=theta_d, exner=exner_d,
        temp_prev=temp_prev_d, ttend_t=ttend_t_d, errflg=0)
    if _log_perf_call is not None:
        try:
            _log_perf_call("kessler_update_run", ncol, nz, dt, time.perf_counter() - _t0)
        except Exception:
            pass

    _t0 = time.perf_counter()
    st_energy_d, e3, gravit = kessler_update_timestep_final_bridge_device(
        nz=nz, cpair=cpair_d, temp=temp_d, zm=zm_d, phis=phis_d,
        st_energy=zeros_d, errflg=0, gravit=gravit)
    if _log_perf_call is not None:
        try:
            _log_perf_call("kessler_update_timestep_final", ncol, nz, dt, time.perf_counter() - _t0)
        except Exception:
            pass

    # ONE batched D2H fetch for all three phases' outputs and errflg scalars.
    temp_prev_out, temp_tend_out, st_energy_out, e1, e2, e3 = jax.device_get(
        (temp_prev_d, ttend_t_d, st_energy_d, e1, e2, e3))

    errflg = max(int(e1), int(e2), int(e3))
    errmsg = "" if errflg == 0 else "kessler_update: bad time splitting"

    # temp_prev/temp_tend/st_energy are zero-copy views into EAMxx's field
    # buffers (same as run()'s writeback) -- results must be written back
    # in place for them to reach EAMxx.
    temp_prev[...] = temp_prev_out
    temp_tend[...] = temp_tend_out
    st_energy[...] = st_energy_out

# Called once per rank from KesslerMicrophysics::finalize_impl() -- flushes
# this rank's in-memory perf totals (accumulated by log_call() in run()/
# kessler_run() above) out to the perf-log CSV.
def finalize():
    if _flush_perf_log is not None:
        try:
            _flush_perf_log()
        except Exception:
            pass