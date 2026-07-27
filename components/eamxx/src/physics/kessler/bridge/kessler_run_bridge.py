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
"""

import numpy as np
import jax.numpy as jnp

def to_row_major_1d(arr):
    """Fortran (col-major) → JAX (row-major) for 1D arrays"""
    return jnp.asarray(np.ascontiguousarray(arr.T), dtype=jnp.float64)

def to_col_major_1d(arr):
    """JAX (row-major) → Fortran (col-major) for 1D arrays"""
    return np.asfortranarray(np.array(arr).T, dtype=np.float64)

def to_row_major_2d(arr):
    """Fortran (col-major) → JAX (row-major) for 2D arrays"""
    return jnp.asarray(np.ascontiguousarray(arr.T), dtype=jnp.float64)

def to_col_major_2d(arr):
    """JAX (row-major) → Fortran (col-major) for 2D arrays"""
    return np.asfortranarray(np.array(arr).T, dtype=np.float64)

from kessler_jax.kessler_run import kessler_run

def kessler_run_bridge(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, pk, theta, qv, qc, qr, precl, relhum, scheme_name, errmsg, errflg, lv, pref, rhoqr):
    """
    Production bridge for kessler_run.

    Strategy: transpose inputs/outputs (CPU), call wrapper directly.
    kessler_run_core is @jax.jit decorated — JIT and GPU execution
    are handled there, not here.

    MODULE VARIABLES (from kessler):
      - lv (INOUT)
      - pref (INOUT)
      - rhoqr (INOUT)
    
    Procedure type: COMPUTE
    Pattern: MODULE vars passed as INOUT parameters
    """
    # Convert inputs
    cpair_compute = to_row_major_2d(cpair)
    rair_compute = to_row_major_2d(rair)
    rho_compute = to_row_major_2d(rho)
    z_compute = to_row_major_2d(z)
    pk_compute = to_row_major_2d(pk)
    theta_compute = to_row_major_2d(theta)
    qv_compute = to_row_major_2d(qv)
    qc_compute = to_row_major_2d(qc)
    qr_compute = to_row_major_2d(qr)
    precl_compute = to_row_major_1d(precl)
    relhum_compute = to_row_major_2d(relhum)

    print(f"ncol: {ncol}, nz: {nz}, dt: {dt}, lyr_surf: {lyr_surf}, lyr_toa: {lyr_toa}, latvap: {lv}, pref: {pref}, rhoqr: {rhoqr}")

    print(f"cpair shape: {cpair_compute.shape}")
    print(f"rair shape: {rair_compute.shape}")
    print(f"rho shape: {rho_compute.shape}")
    print(f"z shape: {z_compute.shape}")
    print(f"pk shape: {pk_compute.shape}")
    print(f"theta shape: {theta_compute.shape}")
    print(f"qv shape: {qv_compute.shape}")
    print(f"qc shape: {qc_compute.shape}")
    print(f"qr shape: {qr_compute.shape}")
    print(f"precl shape: {precl_compute.shape}")
    print(f"relhum shape: {relhum_compute.shape}")

    # Compute — {proc_name}_core is @jax.jit decorated; JIT fires on first call
    theta_out, qv_out, qc_out, qr_out, precl_out, relhum_out, scheme_name_out, errmsg_out, errflg_out, lv, pref, rhoqr = kessler_run(ncol, nz, dt, lyr_surf, lyr_toa, cpair_compute, rair_compute, rho_compute, z_compute, pk_compute, theta_compute, qv_compute, qc_compute, qr_compute, precl_compute, relhum_compute, scheme_name, errmsg, errflg, lv, pref, rhoqr)

    # Convert outputs and return MODULE vars (INOUT pattern)
    theta_fortran = to_col_major_2d(theta_out)
    qv_fortran = to_col_major_2d(qv_out)
    qc_fortran = to_col_major_2d(qc_out)
    qr_fortran = to_col_major_2d(qr_out)
    precl_fortran = to_col_major_1d(precl_out)
    relhum_fortran = to_col_major_2d(relhum_out)
    print("kessler_run_bridge: after kessler_run and transposes.")
    return theta_fortran, qv_fortran, qc_fortran, qr_fortran, precl_fortran, relhum_fortran, scheme_name_out, errmsg_out, errflg_out, lv, pref, rhoqr
