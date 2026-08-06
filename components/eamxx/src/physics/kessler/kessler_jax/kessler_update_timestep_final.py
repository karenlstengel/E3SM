# ---------------------------------------------------------------------------
# Translated by: Claude Fable 5 (claude-fable-5)
# Pass: final
# ---------------------------------------------------------------------------

"""
JAX translation of kessler_update_timestep_final (kessler_update.F90).

Computes the dry static energy:
    st_energy = temp * cpair + gravit * zm + phis

Bridge Mode layout: 2-D arrays arrive row-major (nz, ncol) — Fortran
a(i, klev) reads as a[klev-1, i-1] here; phis(i) is 1-D (ncol,), which
broadcasts across the level axis. The Fortran loop nest
(do klev=1,nz / do i=1,n1 with n1=SIZE(cpair,dim=1)) is fully independent →
Level 1 of the VECTORIZATION PRIORITY LADDER. The i loop covers ALL columns
(bound comes from SIZE, not an argument) so axis 1 is used in full; the klev
bound nz is an argument, honored with a static slice [:nz, :].

NOTE on argument order: this subroutine takes errflg BEFORE errmsg — the
signature below preserves that Fortran quirk exactly (bridge contract).

MODULE VARIABLES (from kessler_update):
    gravit: MODULE variable (INOUT) — read here, threaded per
    MODULE_VARIABLE_POLICY and returned unchanged.
"""

import os
os.environ["JAX_ENABLE_X64"] = "1"

import functools
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
from jax import lax


@functools.partial(jax.jit, static_argnames=("nz",))
def kessler_update_timestep_final_core(nz, cpair, temp, zm, phis, st_energy,
                                       errflg, gravit):
    """
    Core computation for kessler_update_timestep_final. Pure JAX, jit-compiled.

    Operates on row-major (nz, ncol) arrays; bridge handles layout.
    Strings (errmsg) stay in the wrapper per the JIT boundary rules.
    """
    errflg = jnp.asarray(0, dtype=jnp.int32)                     # [JAX]  Fortran: errflg = 0

    # Fortran: st_energy(i,klev) = temp(i,klev)*cpair(i,klev)
    #                            + gravit*zm(i,klev) + phis(i)
    # phis (ncol,) broadcasts across the level axis of the (nz, ncol) arrays.
    st_full = temp * cpair + gravit * zm + phis                  # [JAX-VEC]
    # klev loop bound is the nz argument — write rows [:nz] only; [STATIC-INT] slice
    st_energy = st_energy.at[:nz, :].set(st_full[:nz, :])        # [JAX-VEC]

    return st_energy, errflg, gravit


def kessler_update_timestep_final(nz, cpair, temp, zm, phis, st_energy,
                                  errflg, errmsg, gravit):
    """
    Wrapper for kessler_update_timestep_final (Bridge Mode).

    Handles the errmsg string (JAX cannot); passes everything else to the
    jitted core. Preserves the Fortran errflg-before-errmsg argument order.

    Returns: st_energy, errflg, errmsg, gravit
    """
    errmsg = ""                                                  # [PY]  Fortran: errmsg = ''
    st_energy, errflg, gravit = kessler_update_timestep_final_core(
        nz, cpair, temp, zm, phis, st_energy, errflg, gravit)
    return st_energy, errflg, errmsg, gravit

# Pass-4 self-check: dtype — errflg is explicit int32; all float math inherits
# float64 from the bridge-supplied arrays under JAX_ENABLE_X64. Annotation
# tags present. Signature keeps the Fortran errflg-before-errmsg order and the
# gravit INOUT slot; return order matches the bridge contract exactly.
