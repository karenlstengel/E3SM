# LLM: Claude Fable 5 (claude-fable-5) — generated translation (copied from translations/claude/jax/, header added on promotion to _officialJAX)
# ---------------------------------------------------------------------------
# Translated by: Claude Fable 5 (claude-fable-5)
# Pass: final
# ---------------------------------------------------------------------------
"""JAX translation of the Fortran subroutine `kessler_run` (module `kessler`).

Kessler warm-rain microphysics: autoconversion/accretion of cloud water to
rain, saturation adjustment, rain evaporation, and sub-cycled rain
sedimentation with a CFL-limited inner timestep.

Design:
- Arrays arrive in JAX layout (nz, ncol) — the bridge already transposed
  from Fortran (ncol, nz). Fortran a(col, klev) -> a[klev-1, col-1].
- lyr_surf / lyr_toa / lyr_step are trace-time statics: the level axis is
  gathered into canonical surface-first order (row 0 = surface, row -1 =
  TOA) so the neighbour `klev+lyr_step` is always row j+1, then scattered
  back at the end. When the native order is already surface-first the
  gather/scatter is skipped entirely at trace time.
- The data-dependent sub-cycle (`do while`) is a lax.while_loop per column,
  vmapped over the column axis; all per-level loops inside the body are
  whole-array vectorized ops.

Fortran fidelity notes:
- `f5` is a SCALAR in the Fortran, overwritten on every iteration of the
  first per-level loop; the sub-cycle then reads the surviving value, i.e.
  the one computed at klev = lyr_toa. `cpair`/`rair` are spatially-uniform
  physical constants (scalars here, not per-level fields), so
  `f5 = 4093 * lv / cpair` is exactly that surviving value already.
- Error semantics are per column, without Fortran's early RETURN: a column
  with a bad time split (initial dt0 < 1e-12) skips its sub-cycle and keeps
  its post-floor state (qr floored, theta/qv/qc unchanged, precl 0), and
  errflg comes back 1 if any column is bad. For nonpositive dt all fields
  pass through unchanged (precl/relhum zero) with errflg = 1. In both cases
  the caller aborts on errflg, exactly as with the Fortran.
- relhum is intent(out); levels outside the lyr_surf..lyr_toa range (none
  in practice) are returned as 0.

MODULE VARIABLES (from kessler), INOUT pattern: lv, pref, rhoqr — passed in,
used read-only, returned unchanged.
"""

import os
os.environ["JAX_ENABLE_X64"] = "1"

import functools
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
from jax import lax


@functools.partial(jax.jit, static_argnames=('ncol', 'nz', 'lyr_surf', 'lyr_toa'))
def kessler_run_core(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, pk,
                     theta, qv, qc, qr, precl, relhum, errflg, lv, pref, rhoqr):
    """Pure JAX kernel for kessler_run.

    cpair/rair are spatially-uniform scalars. Other 2-D arrays are
    (nz, ncol); precl is (ncol,). Returns
    theta, qv, qc, qr, precl, relhum, errflg, lv, pref, rhoqr.
    """
    # --- canonical surface-first level order (statics -> trace-time Python) ---
    s = lyr_surf - 1                                   # [STATIC-INT] 1-based -> 0-based
    t = lyr_toa - 1                                    # [STATIC-INT]
    step = 1 if t >= s else -1                         # [STATIC-INT] [PY-IF] on statics: Fortran lyr_step
    order = tuple(range(s, t + step, step))            # [STATIC-INT] canonical level order
    nlev = len(order)                                  # [STATIC-INT]
    # Fortran fidelity: sed(lyr_toa) reads z(col, lyr_toa - lyr_step) — for a
    # one-level range that neighbour lies OUTSIDE the loop range but must be
    # inside the array, exactly as in the Fortran (else it is UB there too).
    assert 0 <= t - step <= nz - 1, \
        "kessler_run: z neighbour of lyr_toa out of bounds"  # [PY] trace-time check on statics
    identity = order == tuple(range(nz))               # [STATIC-INT]

    if identity:                                       # [PY-IF] on [STATIC-INT] — resolved at trace time
        rho_o, z_o, pk_o = rho, z, pk                  # [PY] aliasing only
        theta_o, qv_o, qc_o, qr_o = theta, qv, qc, qr  # [PY] aliasing only
        idx = None                                     # [PY]
    else:
        idx = jnp.asarray(order, dtype=jnp.int32)      # [JAX] static gather indices
        rho_o = rho[idx, :]                            # [JAX-VEC]
        z_o = z[idx, :]                                # [JAX-VEC]
        pk_o = pk[idx, :]                              # [JAX-VEC]
        theta_o = theta[idx, :]                        # [JAX-VEC]
        qv_o = qv[idx, :]                              # [JAX-VEC]
        qc_o = qc[idx, :]                              # [JAX-VEC]
        qr_o = qr[idx, :]                              # [JAX-VEC]

    # --- first per-level loop: r, rhalf, pc, qr floor, initial velqr ---
    f2x = 17.27                                        # [PY] float64 literal (inherits from operands)
    r_o = 0.001 * rho_o                                # [JAX-VEC] r(klev)
    rhalf_o = jnp.sqrt(rho_o[0:1, :] / rho_o)          # [JAX-VEC] sqrt(rho_surf/rho)
    xk = cpair / rair                                  # [PY] cpair/rair are spatially-uniform scalars
    pc_o = 3.8 / (pk_o ** xk * pref)                   # [JAX-VEC] pc(klev)
    qr0_o = jnp.maximum(qr_o, 0.0)                     # [JAX-VEC] Fortran qr floor
    velqr0_o = 36.34 * rhalf_o * (qr0_o * r_o) ** 0.1364  # [JAX-VEC] initial fallspeed
    # Fortran fidelity: f5 is a SCALAR overwritten every klev iteration; the
    # sub-cycle below uses its value from the LAST iteration (klev = lyr_toa).
    # cpair is spatially uniform, so this is a single scalar shared by all columns.
    f5 = 4093.0 * lv / cpair                           # [PY] scalar carry-over
    rho_surf_col = rho_o[0, :]                         # [JAX-VEC] rho at the surface, per column
    # z neighbour of the TOA level, from the FULL native array: equals the
    # canonical zc[-2] whenever nlev >= 2, and is the Fortran
    # z(col, lyr_toa - lyr_step) when the range is a single level.
    z_below_top_col = z[t - step, :]                   # [JAX-VEC] static row index

    def _column(dt_c, rho_surf, z_below_top, r, rhalf, pc, zc, pkc,
                th0, qv0, qc0, qr0, velqr_init):       # [JAX-VMAP] per-column sub-cycle
        """Sub-cycle for one column; all 1-D args are (nlev,) canonical order."""
        dz = zc[1:] - zc[:-1]                          # [JAX-VEC] z(klev+step) - z(klev)

        def _cfl(dt0, velqr):
            if nlev < 2:                               # [PY-IF] on [STATIC-INT]: Fortran CFL loop is empty
                return dt0
            v = velqr[:-1]                             # [JAX-VEC] body levels only
            # JIT safety: mask BEFORE dividing so v == 0 never produces
            # inf/nan inside the jitted graph (both jnp.where branches are
            # always evaluated under jit).
            mask = jnp.abs(v) > 1.0e-12                # [JAX-VEC]
            cand = 0.8 * dz / jnp.where(mask, v, 1.0)  # [JAX-WHERE] guarded division
            cand = jnp.where(mask, cand, jnp.inf)      # [JAX-WHERE] inert where |velqr| <= 1e-12
            return jnp.minimum(dt0, jnp.min(cand))     # [JAX-VEC] sequential min == min-reduce

        dt0_init = _cfl(dt_c, velqr_init)              # [JAX]
        bad_split = (dt_c > 0.0) & (dt0_init < 1.0e-12)  # [JAX] Fortran 'bad time splitting'
        # JIT safety: under jit there is no early `return` — a column with
        # nonpositive dt or a bad time split must SKIP the sub-cycle, or the
        # while_loop condition can never become false (dt0 <= 0 stalls
        # time_counter; dt < 0 diverges it) and the kernel hangs.
        col_ok = (dt_c > 0.0) & jnp.logical_not(bad_split)  # [JAX]

        def cond(state):
            tc, dt0, pacc, th, qvv, qcc, qrr, velqr = state
            return (jnp.abs(dt_c - tc) > 1.0e-5) & col_ok  # [JAX] Fortran do-while condition

        def body(state):
            tc, dt0, pacc, th, qvv, qcc, qrr, velqr = state
            # precipitation rate at the surface, accumulated over the sub-step
            precl_inst = rho_surf * qrr[0] * velqr[0] / rhoqr  # [JAX]
            pacc = pacc + precl_inst * dt0             # [JAX]
            # sedimentation term
            flux = r * qrr * velqr                     # [JAX-VEC]
            sed_body = dt0 * (flux[1:] - flux[:-1]) / (r[:-1] * dz)  # [JAX-VEC]
            sed_top = -dt0 * qrr[-1] * velqr[-1] / (0.5 * (zc[-1] - z_below_top))  # [JAX] z(lyr_toa) - z(lyr_toa-lyr_step)
            sed = jnp.concatenate([sed_body, sed_top[None]])  # [JAX-VEC]
            # autoconversion + accretion (qrprod), then updates in Fortran order
            qrprod = qcc - (qcc - dt0 * jnp.maximum(0.001 * (qcc - 0.001), 0.0)) \
                / (1.0 + dt0 * 2.2 * qrr ** 0.875)     # [JAX-VEC]
            qcc = jnp.maximum(qcc - qrprod, 0.0)       # [JAX-VEC]
            qrr = jnp.maximum(qrr + qrprod + sed, 0.0)  # [JAX-VEC]
            # saturation adjustment and rain evaporation
            pt = pkc * th                              # [JAX-VEC]
            qvs = pc * jnp.exp(f2x * (pt - 273.0) / (pt - 36.0))  # [JAX-VEC]
            prod = (qvv - qvs) / (1.0 + qvs * f5 / (pt - 36.0) ** 2)  # [JAX-VEC] f5 = scalar carry-over
            rqr = r * qrr                              # [JAX-VEC]
            ern = jnp.minimum(                         # [JAX-VEC] Fortran 3-arg MIN
                jnp.minimum(
                    dt0 * (((1.6 + 124.9 * rqr ** 0.2046) * rqr ** 0.525)
                           / (2550000.0 * pc / (3.8 * qvs) + 540000.0))
                    * (jnp.maximum(qvs - qvv, 0.0) / (r * qvs)),  # [JAX-VEC] DIM(qvs,qv) = max(qvs-qv,0)
                    jnp.maximum(-prod - qcc, 0.0)),
                qrr)
            th = th + lv / (cpair * pkc) * (jnp.maximum(prod, -qcc) - ern)  # [JAX-VEC] cpair: closed-over scalar
            qvv = jnp.maximum(qvv - jnp.maximum(prod, -qcc) + ern, 0.0)  # [JAX-VEC]
            qcc = qcc + jnp.maximum(prod, -qcc)        # [JAX-VEC]
            qrr = jnp.maximum(qrr - ern, 0.0)          # [JAX-VEC]
            # advance the sub-cycle clock, recompute fallspeed and next dt0
            tc = tc + dt0                              # [JAX]
            velqr = 36.34 * rhalf * (qrr * r) ** 0.1364  # [JAX-VEC]
            dt0 = _cfl(jnp.maximum(dt_c - tc, 0.0), velqr)  # [JAX]
            return (tc, dt0, pacc, th, qvv, qcc, qrr, velqr)

        # JIT safety: pin the scalar carry entries to strong float64 —
        # weak-typed Python 0.0 initializers change dtype after the first
        # body iteration and lax.while_loop rejects the carry mismatch.
        state0 = (jnp.asarray(0.0, dtype=jnp.float64),  # [JAX] time_counter
                  dt0_init,
                  jnp.asarray(0.0, dtype=jnp.float64),  # [JAX] precl_acc
                  th0, qv0, qc0, qr0, velqr_init)
        tc, dt0_f, pacc, th, qvv, qcc, qrr, velqr = lax.while_loop(  # [JAX-WHILE] Fortran do-while sub-cycle
            cond, body, state0)

        precl_col = pacc / dt_c                        # [JAX] precl_acc / dt
        # final relative humidity from the updated state
        pt = pkc * th                                  # [JAX-VEC]
        qvs = pc * jnp.exp(f2x * (pt - 273.0) / (pt - 36.0))  # [JAX-VEC]
        relhum_col = qvv / qvs * 100.0                 # [JAX-VEC]
        return th, qvv, qcc, qrr, precl_col, relhum_col, bad_split

    _vcolumn = jax.vmap(                               # [JAX-VMAP] columns are independent
        _column,
        in_axes=(None, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
        out_axes=(1, 1, 1, 1, 0, 1, 0),
    )
    (theta_n, qv_n, qc_n, qr_n, precl_n, relhum_n, bad_split) = _vcolumn(
        dt, rho_surf_col, z_below_top_col, r_o, rhalf_o, pc_o, z_o,
        pk_o, theta_o, qv_o, qc_o, qr0_o, velqr0_o)

    # --- scatter back to the native level order ---
    if identity:                                       # [PY-IF] on [STATIC-INT] — resolved at trace time
        theta_out, qv_out, qc_out, qr_out, relhum_out = (
            theta_n, qv_n, qc_n, qr_n, relhum_n)       # [PY] aliasing only
    else:
        theta_out = theta.at[idx, :].set(theta_n)      # [JAX-VEC] inverse of the entry gather
        qv_out = qv.at[idx, :].set(qv_n)               # [JAX-VEC]
        qc_out = qc.at[idx, :].set(qc_n)               # [JAX-VEC]
        qr_out = qr.at[idx, :].set(qr_n)               # [JAX-VEC]
        relhum_out = jnp.zeros_like(relhum).at[idx, :].set(relhum_n)  # [JAX-VEC] intent(out): base 0

    # --- error flags: nonpositive dt, bad time splitting ---
    dt_pos = dt > 0.0                                  # [JAX]
    theta_out = jnp.where(dt_pos, theta_out, theta)    # [JAX-WHERE] dt<=0: inputs pass through
    qv_out = jnp.where(dt_pos, qv_out, qv)             # [JAX-WHERE]
    qc_out = jnp.where(dt_pos, qc_out, qc)             # [JAX-WHERE]
    qr_out = jnp.where(dt_pos, qr_out, qr)             # [JAX-WHERE] original (unfloored) qr
    precl_out = jnp.where(dt_pos, precl_n, jnp.zeros_like(precl_n))  # [JAX-WHERE]
    relhum_out = jnp.where(dt_pos, relhum_out, jnp.zeros_like(relhum_out))  # [JAX-WHERE]
    errflg = jnp.where(dt_pos,                         # [JAX-WHERE] Fortran errflg = 1 branches
                       jnp.where(jnp.any(bad_split), 1, 0),
                       1).astype(jnp.int32)            # [JAX] pin dtype: Fortran INTEGER; fresh value, incoming errflg never read

    return (theta_out, qv_out, qc_out, qr_out, precl_out, relhum_out,
            errflg, lv, pref, rhoqr)


def kessler_run(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, pk,
                theta, qv, qc, qr, precl, relhum, scheme_name, errmsg, errflg,
                lv, pref, rhoqr):
    """Wrapper for kessler_run (Bridge Mode): strings handled here, numerics
    in kessler_run_core. Returns the full 12-tuple in Fortran argument order.
    """
    (theta, qv, qc, qr, precl, relhum, errflg,
     lv, pref, rhoqr) = kessler_run_core(
        ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, pk,
        theta, qv, qc, qr, precl, relhum, errflg, lv, pref, rhoqr)

    scheme_name = "KESSLER"                            # [PY] CHARACTER outputs live in the wrapper
    errmsg = ""                                        # [PY]
    if int(errflg) != 0:                               # [PY-IF] host-side concretization (wrapper only)
        if float(dt) <= 0.0:                           # [PY-IF] pick the matching Fortran message
            errmsg = "KESSLER called with nonpositive dt"
        else:
            errmsg = f"KESSLER: bad time splitting {float(dt)}"
    return (theta, qv, qc, qr, precl, relhum, scheme_name, errmsg, errflg,
            lv, pref, rhoqr)
