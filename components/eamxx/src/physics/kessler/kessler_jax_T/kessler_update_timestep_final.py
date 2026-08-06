# ---------------------------------------------------------------------------
# Translated by: Claude Fable 5 (claude-fable-5)
# Pass: final (LITERAL baseline — Fortran-faithful transliteration)
# JIT boundary check: PASSED
# ---------------------------------------------------------------------------

"""
LITERAL JAX transliteration of kessler_update_timestep_final (kessler_update.F90).

Computes the dry static energy:
    st_energy(i,klev) = temp(i,klev)*cpair(i,klev) + gravit*zm(i,klev) + phis(i)

Fortran-native layout (NO bridge): 2-D arrays are (ncol, nz) — Fortran
a(i, klev) reads as a[i-1, klev-1], SAME index order, 0-based; phis(i) is
1-D (ncol,). The i-loop bound comes from SIZE(cpair,dim=1) (static shape
under jit); the klev bound nz is an argument. The loop nest (do klev=1,nz
outer / do i=1,n1 inner) is preserved 1:1 as nested serial lax.fori_loop —
NOT vectorized.

NOTE on argument order: this subroutine takes errflg BEFORE errmsg — the
signature below preserves that Fortran quirk exactly (same as the idiomatic
build).

MODULE VARIABLES (from kessler_update):
    gravit: MODULE variable (INOUT) — read here, threaded per
    MODULE_VARIABLE_POLICY and returned unchanged.
"""

import os
os.environ["JAX_ENABLE_X64"] = "1"

import functools
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
from jax import lax


@functools.partial(jax.jit, static_argnames=("nz",))
def kessler_update_timestep_final_core(nz, cpair, temp, zm, phis, st_energy,
                                       errflg, gravit):
    """
    Core computation for kessler_update_timestep_final. Pure JAX, jit-compiled.

    Operates on Fortran-native (ncol, nz) arrays — no bridge, no transpose.
    Strings (errmsg) stay in the wrapper per the JIT boundary rules.
    """
    errflg = jnp.asarray(0, dtype=jnp.int32)                     # [JAX]  Fortran: errflg = 0

    n1 = cpair.shape[0]                                          # [STATIC-INT]  Fortran: SIZE(cpair,dim=1)

    # Fortran: do klev = 1, nz  /  do i = 1, n1   — preserved as nested
    # serial loops, same nesting depth and order (transliteration build).
    def klev_body(klev, st_energy_c):                            # [JAX-FORI]  do klev=1,nz
        def i_body(i, st_energy_i):                              # [JAX-FORI]  do i=1,n1
            # Fortran: st_energy(i,klev) = (temp(i,klev)*cpair(i,klev))
            #            + (gravit*zm(i,klev)) + phis(i)
            return st_energy_i.at[i, klev].set(                  # [JAX]  single-index write
                (temp[i, klev] * cpair[i, klev])
                + (gravit * zm[i, klev]) + phis[i])
        return lax.fori_loop(0, n1, i_body, st_energy_c)

    st_energy = lax.fori_loop(0, nz, klev_body, st_energy)       # [JAX-FORI]

    return st_energy, errflg, gravit


def kessler_update_timestep_final(nz, cpair, temp, zm, phis, st_energy,
                                  errflg, errmsg, gravit):
    """
    Wrapper for kessler_update_timestep_final (Fortran-native — no bridge).

    Handles the errmsg string (JAX cannot); passes everything else to the
    jitted core. Preserves the Fortran errflg-before-errmsg argument order.

    Returns: st_energy, errflg, errmsg, gravit
    """
    errmsg = ""                                                  # [PY]  Fortran: errmsg = ''
    st_energy, errflg, gravit = kessler_update_timestep_final_core(
        nz, cpair, temp, zm, phis, st_energy, errflg, gravit)
    return st_energy, errflg, errmsg, gravit

# Pass-3 self-check: both Fortran do loops preserved as lax.fori_loop, same
# nesting depth/order; single-index .at[i, klev].set writes; no vmap, no
# broadcast collapse; errflg-before-errmsg order and gravit INOUT threading
# identical to the idiomatic build; layout Fortran-native (ncol, nz).


# ---------------------------------------------------------------------------
# Vectorized variant (NOT wired up / not called by kessler.py).
#
# kessler_update_timestep_final is a pure elementwise update — no dependency
# between grid points — so the nested serial lax.fori_loop above lowers to
# nz*n1 sequential scalar scatter-writes, each a real kernel-dispatch round
# trip on GPU. This version computes the identical result as one elementwise
# expression (with phis broadcast over the level axis) on the native
# (ncol, nz) layout — no transpose needed.
# ---------------------------------------------------------------------------

@functools.partial(jax.jit, static_argnames=("nz",))
def kessler_update_timestep_final_core_vec(nz, cpair, temp, zm, phis, st_energy,
                                            errflg, gravit):
    """
    Vectorized core computation for kessler_update_timestep_final.

    Same math as kessler_update_timestep_final_core, expressed as a single
    elementwise expression instead of a nested lax.fori_loop over (ncol, nz).
    nz is accepted for signature compatibility but unused (shapes are static).
    """
    errflg = jnp.asarray(0, dtype=jnp.int32)                     # [JAX]  Fortran: errflg = 0

    # Fortran: st_energy(i,klev) = (temp(i,klev)*cpair(i,klev))
    #            + (gravit*zm(i,klev)) + phis(i)
    st_energy = (temp * cpair) + (gravit * zm) + phis[:, None]   # [JAX]  phis(i) broadcast over klev

    return st_energy, errflg, gravit


def kessler_update_timestep_final_vec(nz, cpair, temp, zm, phis, st_energy,
                                       errflg, errmsg, gravit):
    """
    Vectorized wrapper for kessler_update_timestep_final (Fortran-native — no bridge).

    Drop-in alternative to kessler_update_timestep_final: same signature
    (including the errflg-before-errmsg quirk) and return values, computed
    without the nested serial loops.
    """
    errmsg = ""                                                  # [PY]  Fortran: errmsg = ''
    st_energy, errflg, gravit = kessler_update_timestep_final_core_vec(
        nz, cpair, temp, zm, phis, st_energy, errflg, gravit)
    return st_energy, errflg, errmsg, gravit


# ---------------------------------------------------------------------------
# NumPy variant (NOT wired up / not called by kessler.py).
#
# Another pure elementwise update with no control flow — nothing here needs
# jax.jit/tracing. Plain NumPy does the identical computation with no
# compilation step and no device dispatch overhead.
# ---------------------------------------------------------------------------

def kessler_update_timestep_final_np(nz, cpair, temp, zm, phis, st_energy,
                                      errflg, errmsg, gravit):
    """
    NumPy (no JAX) alternative to kessler_update_timestep_final.

    Drop-in alternative: same signature (including the errflg-before-errmsg
    quirk) and return values, using plain NumPy instead of JAX. nz is
    accepted for signature compatibility but unused. Inputs are expected to
    be NumPy arrays (e.g. from py_backend=host).
    """
    errmsg = ""                                                  # [PY]  Fortran: errmsg = ''
    errflg = 0                                                    # [PY]  Fortran: errflg = 0

    # Fortran: st_energy(i,klev) = (temp(i,klev)*cpair(i,klev))
    #            + (gravit*zm(i,klev)) + phis(i)
    st_energy = (temp * cpair) + (gravit * zm) + phis[:, None]    # [NUMPY]  phis(i) broadcast over klev

    return st_energy, errflg, errmsg, gravit
