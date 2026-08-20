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
JAX translation of the Fortran subroutine `kessler_update_run`.

Accumulates the total temperature tendency over the physics time step:
    ttend_t(i,klev) += (theta(i,klev)*exner(i,klev) - temp_prev(i,klev)) / dt
for klev = 1..nz, i = 1..ncol. Arrays arrive in the core as (nz, ncol) float64;
`nz`/`ncol` are static ints so the Fortran loop bounds are static slices.
"""


@functools.partial(jax.jit, static_argnames=("nz", "ncol"))
def kessler_update_run_core(nz, ncol, dt, theta, exner, temp_prev, ttend_t,
                            errflg):
    """
    Pure JAX compute on (nz, ncol) float64 arrays. nz / ncol are static
    Python ints (shape bounds); dt is a traced float64 scalar; errflg is
    returned as a JAX int32 (always 0 — the Fortran has no error branch).
    Returns: ttend_t, errflg
    """
    # Fortran: do klev=1,nz / do i=1,ncol — one fused elementwise update
    # over the [:nz, :ncol] block; elements outside the bounds pass through.
    ttend_t = ttend_t.at[:nz, :ncol].add(                                 # [JAX-VEC]
        (theta[:nz, :ncol] * exner[:nz, :ncol] - temp_prev[:nz, :ncol]) / dt)
    errflg = jnp.asarray(0, dtype=jnp.int32)                              # [JAX] CHANGED: dtype explicit (already) — tagged
    return ttend_t, errflg


def kessler_update_run(nz, ncol, dt, theta, exner, temp_prev, ttend_t,
                       errmsg, errflg):
    """
    Host wrapper: calls the jitted core, owns the CHARACTER output `errmsg`
    (never enters the core) and converts errflg to a Python int.
    Returns: ttend_t, errmsg, errflg
    """
    ttend_t, errflg = kessler_update_run_core(
        nz, ncol, dt, theta, exner, temp_prev, ttend_t, errflg)             # [PY]
    errflg = int(errflg)                                                   # [PY]
    errmsg = ""                                                            # [PY]
    return ttend_t, errmsg, errflg
