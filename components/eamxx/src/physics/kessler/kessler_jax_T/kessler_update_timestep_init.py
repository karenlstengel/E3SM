# ---------------------------------------------------------------------------
# Translated by: Claude Fable 5 (claude-fable-5)
# Pass: final (LITERAL baseline — Fortran-faithful transliteration)
# JIT boundary check: PASSED
# ---------------------------------------------------------------------------

"""
LITERAL JAX transliteration of kessler_update_timestep_init (kessler_update.F90).

Snapshots the pre-physics temperature and zeroes its tendency:
    temp_prev(i,k) = temp(i,k)
    ttend_t(i,k)   = 0

Fortran-native layout (NO bridge): 2-D arrays are (ncol, nz) — Fortran
a(i, k) reads as a[i-1, k-1], SAME index order, 0-based. Both Fortran loop
bounds come from SIZE(temp) (n1 = SIZE(temp,1), n2 = SIZE(temp,2)); shapes
are static under jit, so the bounds are the array's own shape. The loop nest
(do k=1,n2 outer / do i=1,n1 inner) is preserved 1:1 as nested serial
lax.fori_loop — NOT vectorized. Signature, dtypes, x64 header, and errmsg
handling are held identical to the idiomatic build.
"""

import os
os.environ["JAX_ENABLE_X64"] = "1"

import functools
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
from jax import lax


@jax.jit
def kessler_update_timestep_init_core(temp, temp_prev, ttend_t, errflg):
    """
    Core computation for kessler_update_timestep_init. Pure JAX, jit-compiled.

    Operates on Fortran-native (ncol, nz) arrays — no bridge, no transpose.
    Strings (errmsg) stay in the wrapper per the JIT boundary rules.
    """
    errflg = jnp.asarray(0, dtype=jnp.int32)             # [JAX]  Fortran: errflg = 0

    n1 = temp.shape[0]                                   # [STATIC-INT]  Fortran: SIZE(temp,dim=1)
    n2 = temp.shape[1]                                   # [STATIC-INT]  Fortran: SIZE(temp,dim=2)

    # Fortran: do k = 1, n2  /  do i = 1, n1   — preserved as nested serial
    # loops, same nesting depth and order (transliteration build).
    def k_body(k, carry):                                # [JAX-FORI]  do k=1,n2
        def i_body(i, carry_i):                          # [JAX-FORI]  do i=1,n1
            temp_prev_i, ttend_t_i = carry_i
            # Fortran: temp_prev(i,k) = temp(i,k)
            temp_prev_i = temp_prev_i.at[i, k].set(temp[i, k])   # [JAX]
            # Fortran: ttend_t(i,k) = 0._kind_phys
            ttend_t_i = ttend_t_i.at[i, k].set(0.0)              # [JAX]
            return temp_prev_i, ttend_t_i
        return lax.fori_loop(0, n1, i_body, carry)

    temp_prev, ttend_t = lax.fori_loop(0, n2, k_body, (temp_prev, ttend_t))  # [JAX-FORI]

    return temp_prev, ttend_t, errflg


def kessler_update_timestep_init(temp, temp_prev, ttend_t, errmsg, errflg):
    """
    Wrapper for kessler_update_timestep_init (Fortran-native — no bridge).

    Handles the errmsg string (JAX cannot); passes everything else to the
    jitted core.

    Returns: temp_prev, ttend_t, errmsg, errflg
    """
    errmsg = ""                                          # [PY]  Fortran: errmsg = ''
    temp_prev, ttend_t, errflg = kessler_update_timestep_init_core(
        temp, temp_prev, ttend_t, errflg)
    return temp_prev, ttend_t, errmsg, errflg

# Pass-3 self-check: both Fortran do loops preserved as lax.fori_loop, same
# nesting depth/order; single-index .at[i, k].set writes; no vmap, no
# broadcast collapse; dtypes/x64/signature identical to the idiomatic build;
# layout Fortran-native (ncol, nz).


# ---------------------------------------------------------------------------
# Vectorized variant (NOT wired up / not called by kessler.py).
#
# kessler_update_timestep_init is a pure elementwise assignment — there is no
# dependency between grid points — so the nested serial lax.fori_loop above
# lowers to n1*n2 sequential scalar scatter-writes. On CPU that's tolerable;
# on GPU each loop iteration costs a real kernel-dispatch round trip, so for
# ne30-size grids this is a likely contributor to the GPU "hang" (it never
# deadlocks, it's just catastrophically slow). This version computes the same
# result as one elementwise expression, still on the native (ncol, nz)
# layout — no transpose needed.
# ---------------------------------------------------------------------------

@jax.jit
def kessler_update_timestep_init_core_vec(temp, temp_prev, ttend_t, errflg):
    """
    Vectorized core computation for kessler_update_timestep_init.

    Same math as kessler_update_timestep_init_core, expressed as a single
    elementwise assignment instead of a nested lax.fori_loop over (ncol, nz).
    """
    errflg = jnp.asarray(0, dtype=jnp.int32)             # [JAX]  Fortran: errflg = 0

    temp_prev = temp                                     # [JAX]  temp_prev(i,k) = temp(i,k)
    ttend_t   = jnp.zeros_like(ttend_t)                  # [JAX]  ttend_t(i,k) = 0._kind_phys

    return temp_prev, ttend_t, errflg


def kessler_update_timestep_init_vec(temp, temp_prev, ttend_t, errmsg, errflg):
    """
    Vectorized wrapper for kessler_update_timestep_init (Fortran-native — no bridge).

    Drop-in alternative to kessler_update_timestep_init: same signature and
    return values, computed without the nested serial loops.
    """
    errmsg = ""                                          # [PY]  Fortran: errmsg = ''
    temp_prev, ttend_t, errflg = kessler_update_timestep_init_core_vec(
        temp, temp_prev, ttend_t, errflg)
    return temp_prev, ttend_t, errmsg, errflg


# ---------------------------------------------------------------------------
# NumPy variant (NOT wired up / not called by kessler.py).
#
# This function has no control flow and nothing to differentiate — it's a
# value copy and a zero-fill. There's no need for jax.jit/tracing at all;
# plain NumPy does the same elementwise work with no compilation step and no
# device dispatch overhead. This is a candidate for computing directly on
# the host (or even folding into the C++ Kokkos preprocessing kernel) rather
# than routing through JAX at all.
# ---------------------------------------------------------------------------

def kessler_update_timestep_init_np(temp, temp_prev, ttend_t, errmsg, errflg):
    """
    NumPy (no JAX) alternative to kessler_update_timestep_init.

    Drop-in alternative: same signature and return values, using plain
    NumPy instead of JAX. temp/temp_prev/ttend_t are expected to be NumPy
    arrays (e.g. from py_backend=host).
    """
    errmsg = ""                                          # [PY]  Fortran: errmsg = ''
    errflg = 0                                            # [PY]  Fortran: errflg = 0

    temp_prev = np.array(temp, copy=True)                 # [NUMPY]  temp_prev(i,k) = temp(i,k)
    ttend_t   = np.zeros_like(ttend_t)                    # [NUMPY]  ttend_t(i,k) = 0._kind_phys

    return temp_prev, ttend_t, errmsg, errflg
