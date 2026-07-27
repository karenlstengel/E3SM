import os
import numpy as np
os.environ["JAX_ENABLE_X64"] = "1"     # MUST precede any jax import
os.environ["JAX_PLATFORMS"] = os.environ.get("JAX_PLATFORMS", "cpu")  # MUST precede any jax import; honors the driver script's export, else forces cpu
import jax
jax.config.update("jax_enable_x64", True)
# recommended: persistent compilation cache (the driver sets this too)
_cache_dir = os.environ.get(
    "JAX_COMPILATION_CACHE_DIR",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),".jax_cache"),
)
jax.config.update("jax_compilation_cache_dir", _cache_dir)
jax.config.update("jax_persistent_cache_min_compile_time_secs", 0)
jax.config.update("jax_persistent_cache_min_entry_size_bytes", -1)

import sys
sys.path.insert(0, "/glade/derecho/scratch/kstengel/E3SM/E3SM/components/eamxx/src/physics/kessler/kessler_jax")   # this is a transposed layout 
sys.path.insert(0, "/glade/derecho/scratch/kstengel/E3SM/E3SM/components/eamxx/src/physics/kessler/kessler_jax_T")   # this is the layout that matches E3SM
sys.path.insert(0, "/glade/derecho/scratch/kstengel/E3SM/E3SM/components/eamxx/src/physics/kessler/bridge")   # exposes bridge/ which does the transpose to use the kessler_jax/ files
from kessler_jax_T.kessler_run import kessler_run
from bridge.kessler_run_bridge import kessler_run_bridge

# Setup a few global variables that we set values for with init() and then use in run()
latvap = None
pref = None
rhoqr = None
gravit = None # Currently not used in the JAX code but adding in case we translate the kessler_update calls. 
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
def run(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, pk, theta, qv, qc, qr, precl, relhum):

    print("Calling kessler_run from kessler_JAX_T")
    global latvap, pref, rhoqr
    global errmsg, errflg, scheme_name

    # Compute — {proc_name}_core is @jax.jit decorated; JIT fires on first call
    # theta_out, qv_out, qc_out, qr_out, precl_out, relhum_out, scheme_name, errmsg, errflg, latvap, pref, rhoqr = kessler_run_bridge(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, pk, theta, qv, qc, qr, precl, relhum, scheme_name, errmsg, errflg, latvap, pref, rhoqr)

    theta_out, qv_out, qc_out, qr_out, precl_out, relhum_out, scheme_name, errmsg, errflg, latvap, pref, rhoqr = kessler_run(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, pk, theta, qv, qc, qr, precl, relhum, scheme_name, errmsg, errflg, latvap, pref, rhoqr)

    # theta/qv/qc/qr/precl/relhum are zero-copy views into EAMxx's field
    # buffers (created in create_py_field()); the C++ caller (py_module_call)
    # discards this function's return value, so results must be written
    # back in place for them to reach EAMxx at all.
    theta[...]  = theta_out
    qv[...]     = qv_out
    qc[...]     = qc_out
    qr[...]     = qr_out
    precl[...]  = precl_out
    relhum[...] = relhum_out

    print("after kessler_run_T")
    print(f"precl max: {np.max(precl)}")
    # Run with no transposes below:
    # theta, qv, qc, qr, precl, relhum, scheme_name, errmsg, errflg, _, _, _ = kessler_run(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, pk, theta, qv, qc, qr, precl, relhum, scheme_name, errmsg, errflg, latvap, pref, rhoqr)
