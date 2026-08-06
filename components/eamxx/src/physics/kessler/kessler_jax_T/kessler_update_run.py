# ---------------------------------------------------------------------------
# Translated by: Claude Fable 5 (claude-fable-5)
# Pass: final (LITERAL baseline — Fortran-faithful transliteration)
# JIT boundary check: PASSED
# ---------------------------------------------------------------------------

"""
LITERAL JAX transliteration of kessler_update_run (kessler_update.F90).

Backs out the total air-temperature tendency:
    ttend_t(i,klev) += (theta(i,klev)*exner(i,klev) - temp_prev(i,klev)) / dt

Fortran-native layout (NO bridge): all 2-D arrays are (ncol, nz) — Fortran
a(i, klev) reads as a[i-1, klev-1] here, SAME index order, 0-based. The
Fortran loop nest (do klev=1,nz outer / do i=1,ncol inner) is preserved 1:1
as nested serial lax.fori_loop — NOT vectorized; that is the idiomatic
build's job. Signature, dtypes, x64 header, and errmsg handling are held
identical to the idiomatic build (only loops + layout differ).
"""

import os
os.environ["JAX_ENABLE_X64"] = "1"

import functools
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
from jax import lax


@functools.partial(jax.jit, static_argnames=("nz", "ncol"))
def kessler_update_run_core(nz, ncol, dt, theta, exner, temp_prev, ttend_t,
                            errflg):
    """
    Core computation for kessler_update_run. Pure JAX, jit-compiled.

    Operates on Fortran-native (ncol, nz) arrays — no bridge, no transpose.
    Strings (errmsg) stay in the wrapper per the JIT boundary rules.
    """
    errflg = jnp.asarray(0, dtype=jnp.int32)                     # [JAX]  Fortran: errflg = 0

    # Fortran: do klev = 1, nz  /  do i = 1, ncol   — preserved as nested
    # serial loops, same nesting depth and order (transliteration build).
    def klev_body(klev, ttend_t_c):                              # [JAX-FORI]  do klev=1,nz
        def i_body(i, ttend_t_i):                                # [JAX-FORI]  do i=1,ncol
            # Fortran: ttend_t(i,klev) = ttend_t(i,klev)
            #            + ((theta(i,klev)*exner(i,klev) - temp_prev(i,klev)) / dt)
            return ttend_t_i.at[i, klev].set(                    # [JAX]  single-index write
                ttend_t_i[i, klev]
                + ((theta[i, klev] * exner[i, klev]
                    - temp_prev[i, klev]) / dt))
        return lax.fori_loop(0, ncol, i_body, ttend_t_c)

    ttend_t = lax.fori_loop(0, nz, klev_body, ttend_t)           # [JAX-FORI]

    return ttend_t, errflg


def kessler_update_run(nz, ncol, dt, theta, exner, temp_prev, ttend_t,
                       errmsg, errflg):
    """
    Wrapper for kessler_update_run (Fortran-native — no bridge).

    Handles the errmsg string (JAX cannot); passes everything else to the
    jitted core.

    Returns: ttend_t, errmsg, errflg
    """
    errmsg = ""                                                  # [PY]  Fortran: errmsg = ''
    ttend_t, errflg = kessler_update_run_core(
        nz, ncol, dt, theta, exner, temp_prev, ttend_t, errflg)
    return ttend_t, errmsg, errflg

# Pass-3 self-check: every Fortran do loop is a lax.fori_loop of the same
# nesting depth and order; no vmap, no broadcast collapse; single-index
# .at[i, klev].set writes; dtypes/x64/signature identical to the idiomatic
# build; layout Fortran-native (ncol, nz).


# ---------------------------------------------------------------------------
# Vectorized variant (NOT wired up / not called by kessler.py).
#
# kessler_update_run is a pure elementwise update — no dependency between
# grid points — so the nested serial lax.fori_loop above lowers to nz*ncol
# sequential scalar scatter-writes. That's the likely cause of the GPU hang:
# each of those iterations is a real kernel-dispatch round trip on GPU, vs.
# negligible overhead on CPU. This version computes the identical result as
# one elementwise expression on the native (ncol, nz) layout — no transpose.
# ---------------------------------------------------------------------------

@jax.jit
def kessler_update_run_core_vec(nz, ncol, dt, theta, exner, temp_prev, ttend_t,
                                 errflg):
    """
    Vectorized core computation for kessler_update_run.

    Same math as kessler_update_run_core, expressed as a single elementwise
    expression instead of a nested lax.fori_loop over (ncol, nz). nz/ncol are
    accepted for signature compatibility but unused (shapes are static).
    """
    errflg = jnp.asarray(0, dtype=jnp.int32)                     # [JAX]  Fortran: errflg = 0

    # Fortran: ttend_t(i,klev) = ttend_t(i,klev)
    #            + ((theta(i,klev)*exner(i,klev) - temp_prev(i,klev)) / dt)
    ttend_t = ttend_t + ((theta * exner - temp_prev) / dt)        # [JAX]

    return ttend_t, errflg


def kessler_update_run_vec(nz, ncol, dt, theta, exner, temp_prev, ttend_t,
                            errmsg, errflg):
    """
    Vectorized wrapper for kessler_update_run (Fortran-native — no bridge).

    Drop-in alternative to kessler_update_run: same signature and return
    values, computed without the nested serial loops.
    """
    errmsg = ""                                                  # [PY]  Fortran: errmsg = ''
    ttend_t, errflg = kessler_update_run_core_vec(
        nz, ncol, dt, theta, exner, temp_prev, ttend_t, errflg)
    return ttend_t, errmsg, errflg


# ---------------------------------------------------------------------------
# NumPy variant (NOT wired up / not called by kessler.py).
#
# Another pure elementwise update with no control flow — nothing here needs
# jax.jit/tracing. Plain NumPy does the identical computation with no
# compilation step and no device dispatch overhead.
# ---------------------------------------------------------------------------

def kessler_update_run_np(nz, ncol, dt, theta, exner, temp_prev, ttend_t,
                           errmsg, errflg):
    """
    NumPy (no JAX) alternative to kessler_update_run.

    Drop-in alternative: same signature and return values, using plain
    NumPy instead of JAX. nz/ncol are accepted for signature compatibility
    but unused. Inputs are expected to be NumPy arrays (e.g. from
    py_backend=host).
    """
    errmsg = ""                                                  # [PY]  Fortran: errmsg = ''
    errflg = 0                                                    # [PY]  Fortran: errflg = 0

    # Fortran: ttend_t(i,klev) = ttend_t(i,klev)
    #            + ((theta(i,klev)*exner(i,klev) - temp_prev(i,klev)) / dt)
    ttend_t = ttend_t + ((theta * exner - temp_prev) / dt)        # [NUMPY]

    return ttend_t, errmsg, errflg
