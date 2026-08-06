# LLM: Gemini 3.1 PRO via Gemini CLI code agent — generated translation (copied from translations_2sdExp_JDdata/gemini/, header added on promotion to _officialJAX)
"""
Bridge: kessler_run
Strategy: ALWAYS transpose + contiguous; JIT handled by _core decorator

MODULE: kessler
MODULE vars: lv, pref, rhoqr
Procedure type: COMPUTE

MODULE_VARIABLE_POLICY (INOUT Pattern):
- MODULE vars passed as function parameters (INOUT)
- Driver manages MODULE state explicitly
- No hidden state dictionary
- INIT: Returns updated MODULE vars to driver
- COMPUTE: Returns updated MODULE vars to driver
"""

"""
Layout Conversion Functions

STRATEGY:
- ALWAYS transpose (guaranteed correct physics)
- ALWAYS copy to contiguous (predictable performance)
- JIT handled by _core decorator, not by the bridge

Fortran (column-major) ↔ JAX (row-major)

ON-DEVICE TRANSPOSE (CHANGED): the axis-swap itself now happens via
jnp.transpose, on-device, instead of NumPy's arr.T + ascontiguousarray on
the host. The old NumPy version read the (possibly device-resident, e.g.
py_backend=device) input array from the host to do the reshuffle there,
then re-uploaded it -- a real host-side memory copy (and possible
host<->device migration) paid on every call, every timestep. The
kessler_update_*_core functions this bridge calls are unchanged (Level-1
whole-array ops, no vmap/while_loop/gather), so there's no expectation of
the compile-time regression seen when moving kessler_run's transpose
inside its jit core -- this is a pure "do the reshuffle on the device
instead of the host" change, nothing moves inside a jit boundary.
"""

import numpy as np
import jax.numpy as jnp

def to_row_major_1d(arr):
    """Fortran (col-major) → JAX (row-major) for 1D arrays.

    CHANGED: a transpose of a 1D array is a no-op, so the old
    np.ascontiguousarray(arr.T) here was just a wasted host-side copy.
    Upload as-is.
    """
    return jnp.asarray(arr, dtype=jnp.float64)

def to_col_major_1d(arr):
    """JAX (row-major) → Fortran (col-major) for 1D arrays.

    CHANGED: same reasoning as to_row_major_1d -- 1D transpose is a no-op,
    so only the host-side Fortran-contiguity fixup is needed.
    """
    return np.asfortranarray(np.asarray(arr), dtype=np.float64)

def to_row_major_2d(arr):
    """Fortran (col-major) → JAX (row-major) for 2D arrays.

    CHANGED: upload first (jnp.asarray, native layout, no reshuffling),
    then transpose on-device (jnp.transpose) instead of transposing on the
    host with NumPy before upload.
    """
    return jnp.transpose(jnp.asarray(arr, dtype=jnp.float64))

def to_col_major_2d(arr):
    """JAX (row-major) → Fortran (col-major) for 2D arrays.

    CHANGED: transpose on-device first (jnp.transpose, cheap while the
    data is still GPU-resident), then bring the already-correctly-ordered
    result to the host. Only the Fortran-contiguity fixup (a host-only
    concept -- EAMxx/Fortran reads column-major memory) still happens on
    the host; the actual axis-swap no longer does.
    """
    return np.asfortranarray(np.asarray(jnp.transpose(arr)), dtype=np.float64)

# copied from llm-fortran-modernization/laft-kessler-update/_officialJAX/optimized/jax/
from kessler_jax.kessler_update_timestep_init import kessler_update_timestep_init
from kessler_jax.kessler_update_run import kessler_update_run
from kessler_jax.kessler_update_timestep_final import kessler_update_timestep_final

def kessler_update_bridge(ncol, nz, dt, cpair, zm, exner, theta, phis, temp, temp_prev, ttend_t, st_energy, scheme_name, errmsg, errflg, gravit):
    """
    Production bridge for kessler_update.

    Strategy: transpose inputs/outputs (CPU), call wrapper directly.
    kessler_update_*_core are @jax.jit decorated — JIT and GPU execution
    are handled there, not here.

    MODULE VARIABLES (from kessler):
      - lv (INOUT)
      - pref (INOUT)
      - rhoqr (INOUT)
      - gravit (INOUT)

    Procedure type: COMPUTE
    Pattern: MODULE vars passed as INOUT parameters
    """
        # Convert inputs
    temp_compute = to_row_major_2d(temp)
    temp_prev_compute = to_row_major_2d(temp_prev)
    ttend_t_compute = to_row_major_2d(ttend_t)
    theta_compute = to_row_major_2d(theta)
    exner_compute = to_row_major_2d(exner)
    cpair_compute = to_row_major_2d(cpair)
    zm_compute = to_row_major_2d(zm)
    phis_compute = to_row_major_1d(phis)
    st_energy_compute = to_row_major_2d(st_energy)

    # Compute — {proc_name}_core is @jax.jit decorated; JIT fires on first call
    temp_prev_out, ttend_t_out, errmsg_out, errflg_out = kessler_update_timestep_init(temp_compute, temp_prev_compute, ttend_t_compute, errmsg, errflg)

    # Compute — {proc_name}_core is @jax.jit decorated; JIT fires on first call
    ttend_t_out, errmsg_out, errflg_out = kessler_update_run(nz, ncol, dt, theta_compute, exner_compute, temp_prev_out, ttend_t_out, errmsg, errflg)

    # Compute — {proc_name}_core is @jax.jit decorated; JIT fires on first call
    st_energy_out, errflg_out, errmsg_out, gravit = kessler_update_timestep_final(nz, cpair_compute, temp_compute, zm_compute, phis_compute, st_energy_compute, errflg, errmsg, gravit)

    # Convert outputs and return MODULE vars (INOUT pattern)
    temp_prev_fortran = to_col_major_2d(temp_prev_out)
    ttend_t_fortran = to_col_major_2d(ttend_t_out)
    st_energy_fortran = to_col_major_2d(st_energy_out)

    return temp_prev_fortran, ttend_t_fortran, st_energy_fortran, scheme_name, errmsg_out, errflg_out, gravit
