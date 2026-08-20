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
JAX translation of the Fortran subroutine `kessler_update_timestep_final`.

Dry static energy diagnostic at the end of the step:
    st_energy(i,klev) = temp*cpair + gravit*zm + phis(i)
for klev = 1..nz over all columns (n1 = SIZE(cpair,1)). Arrays arrive in the
core as (nz, ncol) float64; phis is (ncol,); `nz` is a static int. The
MODULE variable `gravit` (kessler_update) follows the INOUT pattern: passed
in (set by kessler_update_init), used, returned unchanged.
"""


@functools.partial(jax.jit, static_argnames=("nz",))
def kessler_update_timestep_final_core(nz, cpair, temp, zm, phis, st_energy,
                                       errflg, gravit):
    """
    Pure JAX compute. nz static; cpair/temp/zm/st_energy (nz, ncol) float64;
    phis (ncol,) float64; gravit traced float64 scalar (MODULE var, INOUT).
    Returns: st_energy, errflg, gravit
    """
    # Fortran: do klev=1,nz / do i=1,n1 — one fused elementwise expression on
    # the [:nz, :] block; phis(i) broadcast along the level axis.
    st_energy = st_energy.at[:nz, :].set(                                  # [JAX-VEC]
        temp[:nz, :] * cpair[:nz, :] + gravit * zm[:nz, :] + phis[None, :])
    errflg = jnp.asarray(0, dtype=jnp.int32)                               # [JAX]
    return st_energy, errflg, gravit


def kessler_update_timestep_final(nz, cpair, temp, zm, phis, st_energy,
                                  errflg, errmsg, gravit):
    """
    Host wrapper: calls the jitted core, owns `errmsg` (never enters the
    core), converts errflg to a Python int. Returns: st_energy, errflg, errmsg, gravit
    """
    st_energy, errflg, gravit = kessler_update_timestep_final_core(
        nz, cpair, temp, zm, phis, st_energy, errflg, gravit)              # [PY]
    errflg = int(errflg)                                                    # [PY]
    errmsg = ""                                                             # [PY]
    return st_energy, errflg, errmsg, gravit
