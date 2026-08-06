# ---------------------------------------------------------------------------
# Translated by: Claude Fable 5 (claude-fable-5)
# Pass: final
# ---------------------------------------------------------------------------

"""
JAX translation of kessler_update_timestep_init (kessler_update.F90).

Snapshots the pre-physics temperature and zeroes its tendency:
    temp_prev = temp
    ttend_t   = 0

Bridge Mode layout: 2-D arrays arrive row-major (nz, ncol). Both Fortran loop
bounds come from SIZE(temp) — the loops cover the whole arrays — so this is
Level 1 of the VECTORIZATION PRIORITY LADDER: whole-array ops, no loops, no
slices. There are no integer shape arguments, so the jit decorator needs no
static_argnames.
"""

import os
os.environ["JAX_ENABLE_X64"] = "1"

import functools
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
from jax import lax


@jax.jit
def kessler_update_timestep_init_core(temp, temp_prev, ttend_t, errflg):
    """
    Core computation for kessler_update_timestep_init. Pure JAX, jit-compiled.

    Operates on row-major (nz, ncol) arrays; bridge handles layout.
    Strings (errmsg) stay in the wrapper per the JIT boundary rules.
    """
    errflg = jnp.asarray(0, dtype=jnp.int32)             # [JAX]  Fortran: errflg = 0

    # Fortran: temp_prev(i,k) = temp(i,k)  over the full SIZE(temp) extent
    temp_prev = temp                                     # [JAX-VEC]
    # Fortran: ttend_t(i,k) = 0._kind_phys  over the full extent
    ttend_t = jnp.zeros_like(temp)                       # [JAX-VEC]  float64 via temp

    return temp_prev, ttend_t, errflg


def kessler_update_timestep_init(temp, temp_prev, ttend_t, errmsg, errflg):
    """
    Wrapper for kessler_update_timestep_init (Bridge Mode).

    Handles the errmsg string (JAX cannot); passes everything else to the
    jitted core.

    Returns: temp_prev, ttend_t, errmsg, errflg
    """
    errmsg = ""                                          # [PY]  Fortran: errmsg = ''
    temp_prev, ttend_t, errflg = kessler_update_timestep_init_core(
        temp, temp_prev, ttend_t, errflg)
    return temp_prev, ttend_t, errmsg, errflg

# Pass-4 self-check: dtype — errflg is explicit int32; zeros_like(temp)
# inherits float64 from the bridge-supplied array under JAX_ENABLE_X64.
# Annotation tags present. Signatures and return order match the bridge
# contract exactly.
