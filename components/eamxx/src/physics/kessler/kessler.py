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

# from kessler_jax_T.kessler_run import kessler_run_vec
from kessler_jax.kessler_update_timestep_init import kessler_update_timestep_init
from kessler_jax.kessler_update_run import kessler_update_run
from kessler_jax.kessler_update_timestep_final import kessler_update_timestep_final

# kessler_jax.kessler_run now does its own (ncol, nz) <-> (nz, ncol) transpose
# on-device, inside its @jax.jit core (see the NATIVE-LAYOUT ADAPTER comments
# in kessler_jax/kessler_run.py) -- so it can be called directly on EAMxx's
# native-layout arrays, without going through bridge.kessler_run_bridge's
# NumPy-level (host-side, non-jit) transpose.
from kessler_jax.kessler_run import kessler_run

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
    
    latvap = lv_in
    pref = pref_in
    rhoqr = rhoqr_in
    gravit = gravit_in

# calls the kessler_run function from https://github.com/NCAR/llm-fortran-modernization/tree/main/fortran2jax-kessler/_officialJAX
# these arrays should be automatically updated back in EAMxx if everything is setup correctly.
def run(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, zm, pk, theta, qv, qc, qr, precl, relhum):

    global latvap, pref, rhoqr
    global errmsg, errflg, scheme_name

    # Compute — {proc_name}_core is @jax.jit decorated; JIT fires on first call.
    # Called directly on native (ncol, nz) arrays -- no bridge/NumPy transpose.
    theta_out, qv_out, qc_out, qr_out, precl_out, relhum_out, scheme_name, errmsg, errflg, latvap, pref, rhoqr = kessler_run(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, zm, pk, theta, qv, qc, qr, precl, relhum, scheme_name, errmsg, errflg, latvap, pref, rhoqr)

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

    theta[...]  = theta_out
    qv[...]     = qv_out
    qc[...]     = qc_out
    qr[...]     = qr_out
    precl[...]  = precl_out
    relhum[...] = relhum_out

def update(ncol, nz, dt, cpair, zm, pk, theta, phis, temp, temp_prev, temp_tend, st_energy):
    global gravit
    global errmsg, errflg, scheme_name


    # jnp.transpose (not `.T`) so these dispatch through JAX, not NumPy --
    # `.T` on a plain pybind11/NumPy array always resolves to NumPy's own
    # implementation, regardless of what's imported here; jnp.transpose
    # accepts the raw array-like directly and does the JAX conversion +
    # transpose together. phis is 1D, so it's passed through untransposed
    # below (a transpose would be a no-op there anyway).
    temp_compute = jnp.transpose(temp)
    temp_prev_compute = jnp.transpose(temp_prev)
    temp_tend_compute = jnp.transpose(temp_tend)
    theta_compute = jnp.transpose(theta)
    pk_compute = jnp.transpose(pk)
    cpair_compute = jnp.transpose(cpair)
    zm_compute = jnp.transpose(zm)
    st_energy_compute = jnp.transpose(st_energy)

    # Call the kessler_update functions

    print("Calling kessler_update_* kessler_jax/kessler_update_*.py")
    # Compute — {proc_name}_core is @jax.jit decorated; JIT fires on first call
    temp_prev_out, temp_tend_out, errmsg, errflg = kessler_update_timestep_init(temp_compute, temp_prev_compute, temp_tend_compute, errmsg, errflg)

    # Compute — {proc_name}_core is @jax.jit decorated; JIT fires on first call
    temp_tend_out, errmsg, errflg = kessler_update_run(nz, ncol, dt, theta_compute, pk_compute, temp_prev_out, temp_tend_out, errmsg, errflg)

    # Compute — {proc_name}_core is @jax.jit decorated; JIT fires on first call
    st_energy_out, errflg, errmsg, gravit = kessler_update_timestep_final(nz, cpair_compute, temp_compute, zm_compute, phis, st_energy_compute, errflg, errmsg, gravit)

    # These are already JAX arrays so we can just use .T to transpose them back to EAMxx's native layout (ncol, nz) and write them back in place.
    temp_prev[...] = temp_prev_out.T
    temp_tend[...] = temp_tend_out.T
    st_energy[...] = st_energy_out.T