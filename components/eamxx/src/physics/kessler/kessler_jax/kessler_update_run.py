# ---------------------------------------------------------------------------
# Translated by: Claude Fable 5 (claude-fable-5)
# Pass: final
# ---------------------------------------------------------------------------

"""
JAX translation of kessler_update_run (kessler_update.F90).

Backs out the total air-temperature tendency from the updated potential
temperature:  ttend_t += (theta * exner - temp_prev) / dt

Bridge Mode layout: all 2-D arrays arrive row-major (nz, ncol) — Fortran
a(i, klev) reads as a[klev-1, i-1] here. The Fortran loop nest
(do klev=1,nz / do i=1,ncol) is fully independent in both dimensions →
Level 1 of the VECTORIZATION PRIORITY LADDER: one fused jnp expression,
no loops. The klev/i bounds are honored with static slices [:nz, :ncol]
so elements outside the loop range pass through unchanged, exactly as in
Fortran.
"""

import os
os.environ["JAX_ENABLE_X64"] = "1"

import functools
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
from jax import lax


@functools.partial(jax.jit, static_argnames=("nz", "ncol"))
def kessler_update_run_core(nz, ncol, dt, theta, exner, temp_prev, ttend_t,
                            errflg):
    """
    Core computation for kessler_update_run. Pure JAX, jit-compiled.

    Operates on row-major (nz, ncol) arrays; bridge handles layout.
    Strings (errmsg) stay in the wrapper per the JIT boundary rules.
    """
    errflg = jnp.asarray(0, dtype=jnp.int32)                     # [JAX]  Fortran: errflg = 0

    # Fortran: ttend_t(i,klev) = ttend_t(i,klev)
    #            + (theta(i,klev)*exner(i,klev) - temp_prev(i,klev)) / dt
    # over klev=1..nz (axis 0), i=1..ncol (axis 1) — slices are [STATIC-INT]
    tend_inc = (theta[:nz, :ncol] * exner[:nz, :ncol]
                - temp_prev[:nz, :ncol]) / dt                     # [JAX-VEC]
    ttend_t = ttend_t.at[:nz, :ncol].add(tend_inc)               # [JAX-VEC]

    return ttend_t, errflg


def kessler_update_run(nz, ncol, dt, theta, exner, temp_prev, ttend_t,
                       errmsg, errflg):
    """
    Wrapper for kessler_update_run (Bridge Mode).

    Handles the errmsg string (JAX cannot); passes everything else to the
    jitted core.

    Returns: ttend_t, errmsg, errflg
    """
    errmsg = ""                                                  # [PY]  Fortran: errmsg = ''
    ttend_t, errflg = kessler_update_run_core(
        nz, ncol, dt, theta, exner, temp_prev, ttend_t, errflg)
    return ttend_t, errmsg, errflg

# Pass-4 self-check: dtype — the only allocation is errflg (explicit int32);
# all float math inherits float64 from the bridge-supplied arrays under
# JAX_ENABLE_X64. Annotation tags present on every non-trivial statement.
# Signatures and return order match the bridge contract exactly.
