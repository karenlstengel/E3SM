# ---------------------------------------------------------------------------
# Translated by: Claude Fable 5 (claude-fable-5)
# Pass: final
# ---------------------------------------------------------------------------
import os
os.environ["JAX_ENABLE_X64"] = "1"

import functools
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
from jax import lax

"""
JAX translation of the Fortran subroutine `kessler_update_timestep_init`.

Start-of-step bookkeeping over the FULL array extent (n1 = SIZE(temp,1),
n2 = SIZE(temp,2)): snapshot the pre-physics temperature and zero the
tendency accumulator:
    temp_prev(i,k) = temp(i,k);  ttend_t(i,k) = 0
Arrays arrive in the core as (nz, ncol) float64; no static ints (the bounds
are the whole arrays).
"""


@jax.jit
def kessler_update_timestep_init_core(temp, temp_prev, ttend_t, errflg):
    """
    Pure JAX compute on (nz, ncol) float64 arrays; no static ints. The
    intent(out) arguments temp_prev / ttend_t are fully overwritten over the
    whole extent, so their incoming values are unused. Returns: temp_prev, ttend_t, errflg
    """
    temp_prev = temp.astype(jnp.float64)                                   # [JAX-VEC]
    ttend_t = jnp.zeros_like(temp, dtype=jnp.float64)                      # [JAX-VEC]
    errflg = jnp.asarray(0, dtype=jnp.int32)                               # [JAX]
    return temp_prev, ttend_t, errflg


def kessler_update_timestep_init(temp, temp_prev, ttend_t, errmsg, errflg):
    """
    Host wrapper: calls the jitted core, owns `errmsg` (never enters the
    core), converts errflg to a Python int. Returns: temp_prev, ttend_t, errmsg, errflg
    """
    temp_prev, ttend_t, errflg = kessler_update_timestep_init_core(
        temp, temp_prev, ttend_t, errflg)                                   # [PY]
    errflg = int(errflg)                                                    # [PY]
    errmsg = ""                                                             # [PY]
    return temp_prev, ttend_t, errmsg, errflg
