import os
os.environ["JAX_ENABLE_X64"] = "1"     # MUST precede any jax import
import jax
jax.config.update("jax_enable_x64", True)
# recommended: persistent compilation cache (the driver sets this too)
jax.config.update("jax_compilation_cache_dir", "JAX_cache/")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 0)

import sys
sys.path.insert(0, "/path/to/_officialJAX")   # exposes kessler_jax/ and bridge/
from kessler_jax.kessler_init import kessler_init
from kessler_jax.kessler_run import kessler_run

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

    global latvap, pref, rhoqr
    global errmsg, errflg, scheme_name

    # Compute — {proc_name}_core is @jax.jit decorated; JIT fires on first call
    theta, qv, qc, qr, precl, relhum, scheme_name, errmsg, errflg, _, _, _ = kessler_run(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, pk, theta, qv, qc, qr, precl, relhum, scheme_name, errmsg, errflg, latvap, pref, rhoqr)
