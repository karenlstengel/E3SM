# ---------------------------------------------------------------------------
# Translated by: Claude Opus 4.8 (claude-opus-4-8[1m])
# Build: LITERAL transliteration baseline (_literal_baseline)
# Pass: final
# ---------------------------------------------------------------------------
# JIT boundary check: PASSED - static integers stay outside traced control,
# strings remain in the wrapper, and the core returns only JAX-compatible data.
import os
os.environ["JAX_ENABLE_X64"] = "1"

import functools
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
from jax import lax

"""
Kessler warm-rain microphysics — FORTRAN-FAITHFUL TRANSLITERATION.

This is the loop-for-loop, layout-for-layout baseline (NOT the idiomatic build):
  * Layout is Fortran-native (ncol, nz) — NO bridge, NO transpose. A 2-D array
    is indexed exactly like the Fortran: a(col, klev) -> a[col, klev] (0-based).
  * The Fortran `do col` loop is a SERIAL lax.fori_loop over the column axis
    (axis 0); every `do klev` loop is a serial lax.fori_loop; the `do while`
    subcycle is a lax.while_loop. No jax.vmap, no whole-array broadcast.
  * f5/xk are scalars recomputed inside the level loop exactly as in Fortran;
    f5 is carried out with its final (top-of-atmosphere) value, matching the
    Fortran scalar-carryover into the adjustment loop.

Held identical to the idiomatic build: signature, module-var INOUT threading,
dtypes, x64 header, string handling, 1->0 index-base conversion.
"""


@functools.partial(jax.jit, static_argnames=("ncol", "nz", "lyr_surf", "lyr_toa"))
def kessler_run_core(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, pk,
                     theta, qv, qc, qr, precl, relhum, errflg, lv, pref, rhoqr):
    """
    Pure JAX compute core for Kessler (literal transliteration).

    MODULE VARIABLES (from kessler):
        lv: MODULE variable (INOUT)
        pref: MODULE variable (INOUT)
        rhoqr: MODULE variable (INOUT)
    """
    s_idx = lyr_surf - 1                                   # [STATIC-INT]
    t_idx = lyr_toa - 1                                    # [STATIC-INT]
    lyr_step = -1 if lyr_surf > lyr_toa else 1             # [STATIC-INT]
    nlev = abs(t_idx - s_idx) + 1                          # [STATIC-INT] levels surf..toa inclusive

    f2x = 17.27                                            # [PY] constant

    # Fortran initializes precl and errflg at entry (lines 146, 148).
    precl = jnp.zeros(ncol, dtype=jnp.float64)             # [JAX] precl = 0
    errflg = jnp.asarray(0, dtype=jnp.int32)               # [JAX] errflg = 0

    def col_body(col, state):                              # [JAX-FORI] do col = 1, ncol  (SERIAL)
        theta, qv, qc, qr, precl, relhum, errflg = state

        # Per-column 1-D profiles (shape (nz,)) — Fortran-native: row slice on axis 0.
        cpair_c = cpair[col, :]                            # [JAX] a(col,:) -> a[col,:]
        rair_c  = rair[col, :]
        rho_c   = rho[col, :]
        z_c     = z[col, :]
        pk_c    = pk[col, :]
        theta_c = theta[col, :]
        qv_c    = qv[col, :]
        qc_c    = qc[col, :]
        qr_c    = qr[col, :]

        rho_surf = rho_c[s_idx]                            # [JAX] rho(col, lyr_surf)

        # Per-column scratch (declared inside the body — never threaded through).
        r     = jnp.zeros(nz, dtype=jnp.float64)           # [JAX]
        rhalf = jnp.zeros(nz, dtype=jnp.float64)           # [JAX]
        velqr = jnp.zeros(nz, dtype=jnp.float64)           # [JAX]
        pc    = jnp.zeros(nz, dtype=jnp.float64)           # [JAX]

        # ---- Loop 1: constants + terminal velocity (do klev = surf, toa, step) ----
        def loop1(i, carry):                               # [JAX-FORI]
            r, rhalf, pc, velqr, qr_c, f5 = carry
            klev = s_idx + i * lyr_step                    # [JAX] running level index
            f5 = 4093.0 * lv / cpair_c[klev]               # [JAX] scalar (carried out at TOA value)
            xk = cpair_c[klev] / rair_c[klev]              # [JAX] 1/kappa = cp/R
            r     = r.at[klev].set(0.001 * rho_c[klev])    # [JAX] r(klev)
            rhalf = rhalf.at[klev].set(jnp.sqrt(rho_surf / rho_c[klev]))  # [JAX] rhalf(klev)
            pc    = pc.at[klev].set(3.8 / ((pk_c[klev] ** xk) * pref))    # [JAX] pc(klev)
            qr_v  = jnp.maximum(qr_c[klev], 0.0)           # [JAX] MAX(qr,0)
            qr_c  = qr_c.at[klev].set(qr_v)
            velqr = velqr.at[klev].set(36.34 * rhalf[klev] * (qr_v * r[klev]) ** 0.1364)  # [JAX] velqr(klev)
            return r, rhalf, pc, velqr, qr_c, f5

        r, rhalf, pc, velqr, qr_c, f5 = lax.fori_loop(
            0, nlev, loop1, (r, rhalf, pc, velqr, qr_c, jnp.asarray(0.0, dtype=jnp.float64))
        )

        # ---- CFL time step (do klev = surf, toa - step, step) ----
        def cfl_loop(i, dt0):                              # [JAX-FORI]
            klev = s_idx + i * lyr_step                    # [JAX]
            v = velqr[klev]                                # [JAX]
            active = jnp.abs(v) > 1.0e-12                  # [JAX-WHERE] guard (was PY-IF in Fortran)
            cand = 0.8 * (z_c[klev + lyr_step] - z_c[klev]) / jnp.where(active, v, 1.0)  # [JAX]
            return jnp.where(active, jnp.minimum(dt0, cand), dt0)  # [JAX-WHERE]

        dt0 = lax.fori_loop(0, nlev - 1, cfl_loop, dt)     # [JAX] dt0 seeded with dt

        # Bad time-splitting: Fortran returns errflg=1. Here: flag the column and
        # skip its subcycle (the while condition below gates on err_c == 0).
        err_c = jnp.where(dt0 < 1.0e-12,                   # [JAX-WHERE]
                          jnp.asarray(1, dtype=jnp.int32),
                          jnp.asarray(0, dtype=jnp.int32))

        # ---- Subcycle (do while abs(dt - time_counter) > 1e-5) ----
        def cond(carry):                                   # [JAX-WHILE] condition
            _th, _qv, _qc, _qr, _v, _pacc, _d0, tc = carry
            return (jnp.abs(dt - tc) > 1.0e-5) & (err_c == 0)  # [JAX]

        def body(carry):                                   # [JAX-WHILE] body
            theta_c, qv_c, qc_c, qr_c, velqr, precl_acc, dt0, tc = carry

            # Precipitation rate over the subcycled step, accumulate (time-weighted).
            precl_now = rho_surf * qr_c[s_idx] * velqr[s_idx] / rhoqr  # [JAX]
            precl_acc = precl_acc + precl_now * dt0        # [JAX]

            # Mass-weighted sedimentation (do klev = surf, toa - step, step) + TOA term.
            def sed_loop(i, sed):                          # [JAX-FORI]
                klev = s_idx + i * lyr_step                # [JAX]
                kn = klev + lyr_step                       # [JAX]
                num = (r[kn] * qr_c[kn] * velqr[kn]) - (r[klev] * qr_c[klev] * velqr[klev])  # [JAX]
                den = r[klev] * (z_c[kn] - z_c[klev])      # [JAX]
                return sed.at[klev].set(dt0 * num / den)   # [JAX] sed(klev)

            sed = lax.fori_loop(0, nlev - 1, sed_loop, jnp.zeros(nz, dtype=jnp.float64))
            sed = sed.at[t_idx].set(                       # [JAX] sed(lyr_toa) boundary statement
                -dt0 * qr_c[t_idx] * velqr[t_idx]
                / (0.5 * (z_c[t_idx] - z_c[t_idx - lyr_step]))
            )

            # ---- Adjustment terms (do klev = surf, toa, step) ----
            def adj_loop(i, carry):                        # [JAX-FORI]
                theta_c, qv_c, qc_c, qr_c = carry
                klev = s_idx + i * lyr_step                # [JAX]

                # Autoconversion / collection (uses OLD qc, qr).
                qrprod = qc_c[klev] - (qc_c[klev] - dt0 * jnp.maximum(
                    0.001 * (qc_c[klev] - 0.001), 0.0)) / (1.0 + dt0 * 2.2 * qr_c[klev] ** 0.875)  # [JAX]
                qc_c = qc_c.at[klev].set(jnp.maximum(qc_c[klev] - qrprod, 0.0))          # [JAX] qc updated
                qr_c = qr_c.at[klev].set(jnp.maximum(qr_c[klev] + qrprod + sed[klev], 0.0))  # [JAX] qr updated

                # Teten's saturation mixing ratio; condensation rate (uses updated qr).
                pkth = pk_c[klev] * theta_c[klev]          # [JAX]
                qvs = pc[klev] * jnp.exp(f2x * (pkth - 273.0) / (pkth - 36.0))           # [JAX] qvs
                prod = (qv_c[klev] - qvs) / (1.0 + qvs * f5 / (pkth - 36.0) ** 2)        # [JAX] prod

                ern = jnp.minimum(                          # [JAX] evaporation rate
                    dt0 * (((1.6 + 124.9 * (r[klev] * qr_c[klev]) ** 0.2046)
                            * (r[klev] * qr_c[klev]) ** 0.525)
                           / (2550000.0 * pc[klev] / (3.8 * qvs) + 540000.0))
                    * (jnp.maximum(qvs - qv_c[klev], 0.0) / (r[klev] * qvs)),
                    jnp.minimum(jnp.maximum(-prod - qc_c[klev], 0.0), qr_c[klev]),
                )

                cnd = jnp.maximum(prod, -qc_c[klev])       # [JAX] max(prod, -qc) (qc = post-line-246)
                theta_c = theta_c.at[klev].set(            # [JAX] saturation adjustment
                    theta_c[klev] + lv / (cpair_c[klev] * pk_c[klev]) * (cnd - ern))
                qv_c = qv_c.at[klev].set(jnp.maximum(qv_c[klev] - cnd + ern, 0.0))       # [JAX]
                qc_c = qc_c.at[klev].set(qc_c[klev] + cnd)  # [JAX]
                qr_c = qr_c.at[klev].set(jnp.maximum(qr_c[klev] - ern, 0.0))             # [JAX]
                return theta_c, qv_c, qc_c, qr_c

            theta_c, qv_c, qc_c, qr_c = lax.fori_loop(
                0, nlev, adj_loop, (theta_c, qv_c, qc_c, qr_c))

            # Elapsed time.
            tc = tc + dt0                                  # [JAX]

            # Recalculate terminal velocity (do klev = surf, toa, step).
            def vel_loop(i, velqr):                        # [JAX-FORI]
                klev = s_idx + i * lyr_step                # [JAX]
                return velqr.at[klev].set(
                    36.34 * rhalf[klev] * (qr_c[klev] * r[klev]) ** 0.1364)  # [JAX]
            velqr = lax.fori_loop(0, nlev, vel_loop, velqr)

            # Recompute the time step (do klev = surf, toa - step, step).
            dt0 = jnp.maximum(dt - tc, 0.0)                # [JAX]
            def cfl_loop2(i, dt0):                         # [JAX-FORI]
                klev = s_idx + i * lyr_step                # [JAX]
                v = velqr[klev]                            # [JAX]
                active = jnp.abs(v) > 1.0e-12              # [JAX-WHERE]
                cand = 0.8 * (z_c[klev + lyr_step] - z_c[klev]) / jnp.where(active, v, 1.0)  # [JAX]
                return jnp.where(active, jnp.minimum(dt0, cand), dt0)  # [JAX-WHERE]
            dt0 = lax.fori_loop(0, nlev - 1, cfl_loop2, dt0)

            return theta_c, qv_c, qc_c, qr_c, velqr, precl_acc, dt0, tc

        init = (theta_c, qv_c, qc_c, qr_c, velqr,
                jnp.asarray(0.0, dtype=jnp.float64),       # precl_acc
                dt0,
                jnp.asarray(0.0, dtype=jnp.float64))       # time_counter
        (theta_c, qv_c, qc_c, qr_c, velqr,
         precl_acc, dt0, tc) = lax.while_loop(cond, body, init)  # [JAX-WHILE]

        # Average precipitation rate over the physics time step.
        precl_c = precl_acc / dt                           # [JAX]

        # Diagnostic relative humidity (do klev = surf, toa, step).
        relhum_c = relhum[col, :]                          # [JAX] start from input slice
        def rh_loop(i, relhum_c):                          # [JAX-FORI]
            klev = s_idx + i * lyr_step                    # [JAX]
            pkth = pk_c[klev] * theta_c[klev]              # [JAX]
            qvs = pc[klev] * jnp.exp(f2x * (pkth - 273.0) / (pkth - 36.0))  # [JAX]
            return relhum_c.at[klev].set(qv_c[klev] / qvs * 100.0)          # [JAX] relhum(col,klev)
        relhum_c = lax.fori_loop(0, nlev, rh_loop, relhum_c)

        # Write per-column results back into the (ncol, nz) arrays (row slice, axis 0).
        theta  = theta.at[col, :].set(theta_c)             # [JAX]
        qv     = qv.at[col, :].set(qv_c)                   # [JAX]
        qc     = qc.at[col, :].set(qc_c)                   # [JAX]
        qr     = qr.at[col, :].set(qr_c)                   # [JAX]
        precl  = precl.at[col].set(precl_c)                # [JAX]
        relhum = relhum.at[col, :].set(relhum_c)           # [JAX]
        errflg = jnp.maximum(errflg, err_c)                # [JAX] accumulate error flag
        return theta, qv, qc, qr, precl, relhum, errflg

    theta, qv, qc, qr, precl, relhum, errflg = lax.fori_loop(  # [JAX-FORI] SERIAL column loop
        0, ncol, col_body,
        (theta, qv, qc, qr, precl, relhum, errflg),
    )

    return theta, qv, qc, qr, precl, relhum, errflg, lv, pref, rhoqr


def kessler_run(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, pk,
                theta, qv, qc, qr, precl, relhum, scheme_name, errmsg, errflg,
                lv, pref, rhoqr):
    """
    Host wrapper for Kessler (literal transliteration).

    MODULE VARIABLES (from kessler):
        lv: MODULE variable (INOUT)
        pref: MODULE variable (INOUT)
        rhoqr: MODULE variable (INOUT)
    """
    scheme_name = "KESSLER"                                # [PY]
    errmsg = ""                                            # [PY]
    errflg = 0                                             # [PY]

    if dt <= 0.0:                                          # [PY-IF]
        precl = jnp.zeros(ncol, dtype=jnp.float64)
        return theta, qv, qc, qr, precl, relhum, scheme_name, \
            "KESSLER called with nonpositive dt", 1, lv, pref, rhoqr

    theta, qv, qc, qr, precl, relhum, errflg, lv, pref, rhoqr = kessler_run_core(
        ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, pk,
        theta, qv, qc, qr, precl, relhum, errflg, lv, pref, rhoqr,
    )

    if int(errflg) != 0:                                  # [PY-IF]
        errmsg = f"KESSLER: bad time splitting dt={dt}"

    return theta, qv, qc, qr, precl, relhum, scheme_name, errmsg, errflg, lv, pref, rhoqr
