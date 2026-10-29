"""
Sub-cycle pass counts for the JAX kessler_run (driver_kessler_py.cpp --passes).

kessler_jax/kessler_run.py runs each column's sub-cycle as a lax.while_loop
and vmaps it over the columns. Under vmap, JAX runs one loop until every
column is done and masks the finished ones, so the work is
max(passes) * ncol instead of sum(passes) (C++ branch,
performance_suggestions.md, C1 and D1). This module measures how big that
gap is on the driver's columns, and what the restructurings would buy:

  - per-column pass counts (how many sub-steps each column needs);
  - work of the current vmap(while_loop) vs per-column loops, and vs
    "split": one pass over all columns, then the loop only on the columns
    that need more than one;
  - times of the device-resident bridge on all columns, on all columns with
    no rain (one pass each: what split's first pass costs), and on only the
    columns that need more than one pass (what split's loop costs).

The pass counts come from kessler_run_core_npass below: a COPY of
kessler_jax.kessler_run.kessler_run_core with a counter added to the
while_loop carry (lines marked "# NPASS"). kessler_jax/ is not changed, so the
LAFT translation stays as generated; report() checks that the copy's outputs
match kessler_run_core's bit for bit, which catches the two drifting apart.
"""

import functools
import time

import numpy as np
import jax
import jax.numpy as jnp
from jax import lax

import kessler as K
from bridge.kessler_run_bridge import kessler_run_bridge_device
from kessler_jax.kessler_run import kessler_run_core


# ---------------------------------------------------------------------------
# COPY of kessler_run_core (kessler_jax/kessler_run.py), plus a per-column
# pass counter. Same arguments; returns kessler_run_core's outputs and npass.
# Comments trimmed; see the original for the Fortran fidelity notes.
# ---------------------------------------------------------------------------
@functools.partial(jax.jit, static_argnames=('ncol', 'nz', 'lyr_surf', 'lyr_toa'))
def kessler_run_core_npass(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, pk,
                           theta, qv, qc, qr, precl, relhum, errflg, lv, pref, rhoqr):
    s = lyr_surf - 1
    t = lyr_toa - 1
    step = 1 if t >= s else -1
    order = tuple(range(s, t + step, step))
    nlev = len(order)
    assert 0 <= t - step <= nz - 1, \
        "kessler_run: z neighbour of lyr_toa out of bounds"
    identity = order == tuple(range(nz))

    if identity:
        rho_o, z_o, pk_o = rho, z, pk
        theta_o, qv_o, qc_o, qr_o = theta, qv, qc, qr
        idx = None
    else:
        idx = jnp.asarray(order, dtype=jnp.int32)
        rho_o = rho[idx, :]
        z_o = z[idx, :]
        pk_o = pk[idx, :]
        theta_o = theta[idx, :]
        qv_o = qv[idx, :]
        qc_o = qc[idx, :]
        qr_o = qr[idx, :]

    f2x = 17.27
    r_o = 0.001 * rho_o
    rhalf_o = jnp.sqrt(rho_o[0:1, :] / rho_o)
    xk = cpair / rair
    pc_o = 3.8 / (pk_o ** xk * pref)
    qr0_o = jnp.maximum(qr_o, 0.0)
    velqr0_o = 36.34 * rhalf_o * (qr0_o * r_o) ** 0.1364
    f5 = 4093.0 * lv / cpair
    rho_surf_col = rho_o[0, :]
    z_below_top_col = z[t - step, :]

    def _column(dt_c, rho_surf, z_below_top, r, rhalf, pc, zc, pkc,
                th0, qv0, qc0, qr0, velqr_init):
        dz = zc[1:] - zc[:-1]

        def _cfl(dt0, velqr):
            if nlev < 2:
                return dt0
            v = velqr[:-1]
            mask = jnp.abs(v) > 1.0e-12
            cand = 0.8 * dz / jnp.where(mask, v, 1.0)
            cand = jnp.where(mask, cand, jnp.inf)
            return jnp.minimum(dt0, jnp.min(cand))

        dt0_init = _cfl(dt_c, velqr_init)
        bad_split = (dt_c > 0.0) & (dt0_init < 1.0e-12)
        col_ok = (dt_c > 0.0) & jnp.logical_not(bad_split)

        def cond(state):
            tc, dt0, pacc, th, qvv, qcc, qrr, velqr, n = state           # NPASS
            return (jnp.abs(dt_c - tc) > 1.0e-5) & col_ok

        def body(state):
            tc, dt0, pacc, th, qvv, qcc, qrr, velqr, n = state           # NPASS
            precl_inst = rho_surf * qrr[0] * velqr[0] / rhoqr
            pacc = pacc + precl_inst * dt0
            flux = r * qrr * velqr
            sed_body = dt0 * (flux[1:] - flux[:-1]) / (r[:-1] * dz)
            sed_top = -dt0 * qrr[-1] * velqr[-1] / (0.5 * (zc[-1] - z_below_top))
            sed = jnp.concatenate([sed_body, sed_top[None]])
            qrprod = qcc - (qcc - dt0 * jnp.maximum(0.001 * (qcc - 0.001), 0.0)) \
                / (1.0 + dt0 * 2.2 * qrr ** 0.875)
            qcc = jnp.maximum(qcc - qrprod, 0.0)
            qrr = jnp.maximum(qrr + qrprod + sed, 0.0)
            pt = pkc * th
            qvs = pc * jnp.exp(f2x * (pt - 273.0) / (pt - 36.0))
            prod = (qvv - qvs) / (1.0 + qvs * f5 / (pt - 36.0) ** 2)
            rqr = r * qrr
            ern = jnp.minimum(
                jnp.minimum(
                    dt0 * (((1.6 + 124.9 * rqr ** 0.2046) * rqr ** 0.525)
                           / (2550000.0 * pc / (3.8 * qvs) + 540000.0))
                    * (jnp.maximum(qvs - qvv, 0.0) / (r * qvs)),
                    jnp.maximum(-prod - qcc, 0.0)),
                qrr)
            th = th + lv / (cpair * pkc) * (jnp.maximum(prod, -qcc) - ern)
            qvv = jnp.maximum(qvv - jnp.maximum(prod, -qcc) + ern, 0.0)
            qcc = qcc + jnp.maximum(prod, -qcc)
            qrr = jnp.maximum(qrr - ern, 0.0)
            tc = tc + dt0
            velqr = 36.34 * rhalf * (qrr * r) ** 0.1364
            dt0 = _cfl(jnp.maximum(dt_c - tc, 0.0), velqr)
            return (tc, dt0, pacc, th, qvv, qcc, qrr, velqr, n + 1)      # NPASS

        state0 = (jnp.asarray(0.0, dtype=jnp.float64),
                  dt0_init,
                  jnp.asarray(0.0, dtype=jnp.float64),
                  th0, qv0, qc0, qr0, velqr_init,
                  jnp.asarray(0, dtype=jnp.int32))                       # NPASS
        tc, dt0_f, pacc, th, qvv, qcc, qrr, velqr, n = lax.while_loop(  # NPASS
            cond, body, state0)

        precl_col = pacc / dt_c
        pt = pkc * th
        qvs = pc * jnp.exp(f2x * (pt - 273.0) / (pt - 36.0))
        relhum_col = qvv / qvs * 100.0
        return th, qvv, qcc, qrr, precl_col, relhum_col, bad_split, n   # NPASS

    _vcolumn = jax.vmap(
        _column,
        in_axes=(None, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
        out_axes=(1, 1, 1, 1, 0, 1, 0, 0),                              # NPASS
    )
    (theta_n, qv_n, qc_n, qr_n, precl_n, relhum_n, bad_split, npass) = _vcolumn(  # NPASS
        dt, rho_surf_col, z_below_top_col, r_o, rhalf_o, pc_o, z_o,
        pk_o, theta_o, qv_o, qc_o, qr0_o, velqr0_o)

    if identity:
        theta_out, qv_out, qc_out, qr_out, relhum_out = (
            theta_n, qv_n, qc_n, qr_n, relhum_n)
    else:
        theta_out = theta.at[idx, :].set(theta_n)
        qv_out = qv.at[idx, :].set(qv_n)
        qc_out = qc.at[idx, :].set(qc_n)
        qr_out = qr.at[idx, :].set(qr_n)
        relhum_out = jnp.zeros_like(relhum).at[idx, :].set(relhum_n)

    dt_pos = dt > 0.0
    theta_out = jnp.where(dt_pos, theta_out, theta)
    qv_out = jnp.where(dt_pos, qv_out, qv)
    qc_out = jnp.where(dt_pos, qc_out, qc)
    qr_out = jnp.where(dt_pos, qr_out, qr)
    precl_out = jnp.where(dt_pos, precl_n, jnp.zeros_like(precl_n))
    relhum_out = jnp.where(dt_pos, relhum_out, jnp.zeros_like(relhum_out))
    errflg = jnp.where(dt_pos,
                       jnp.where(jnp.any(bad_split), 1, 0),
                       1).astype(jnp.int32)

    return (theta_out, qv_out, qc_out, qr_out, precl_out, relhum_out,
            errflg, lv, pref, rhoqr, npass)                              # NPASS


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def _time(f, repeat):
    """Min seconds per call of f() (which blocks), after one untimed call."""
    f()
    best = float("inf")
    for _ in range(repeat):
        t0 = time.perf_counter()
        f()
        best = min(best, time.perf_counter() - t0)
    return best


def report(ncol, nz, dt, lyr_surf, lyr_toa, rho, zm, pk, theta, qv, qc, qr, precl, relhum,
           repeat=10):
    """Print pass-count statistics and timings for one input state.

    The arrays are host NumPy arrays in EAMxx's (ncol, nz) order, as
    kessler.run() gets them with py_backend=host; kessler.init() must have run.
    Returns True if the counting copy matches kessler_run_core bit for bit.
    """
    host = [np.ascontiguousarray(a, dtype=np.float64)
            for a in (rho, zm, pk, theta, qv, qc, qr, precl, relhum)]
    dev = [jnp.asarray(a) for a in host]
    scal = dict(cpair=K.cpair, rair=K.rair, errflg=0, lv=K.latvap, pref=K.pref, rhoqr=K.rhoqr)

    # --- the copy vs the original, on the same (nz, ncol) inputs ---------------
    T = [jnp.transpose(a) if a.ndim == 2 else a for a in dev]
    args = dict(ncol=ncol, nz=nz, dt=dt, lyr_surf=lyr_surf, lyr_toa=lyr_toa,
                rho=T[0], z=T[1], pk=T[2], theta=T[3], qv=T[4], qc=T[5], qr=T[6],
                precl=T[7], relhum=T[8], **scal)
    ref = jax.device_get(kessler_run_core(**args)[:7])
    out = jax.device_get(kessler_run_core_npass(**args))
    npass = np.asarray(out[-1], dtype=np.int64)
    same = all(np.array_equal(np.asarray(a), np.asarray(b)) for a, b in zip(ref, out[:7]))

    # --- pass statistics ------------------------------------------------------------
    total = int(npass.sum())
    mx = int(npass.max())
    multi = npass > 1
    n_multi = int(multi.sum())
    mx_multi = int(npass[multi].max()) if n_multi else 0
    work_vmap = mx * ncol                          # what vmap(while_loop) runs
    work_split = ncol + mx_multi * n_multi         # one pass on all, then the loop on the rest
    print(f"\n--- kessler_run sub-cycle passes (ncol={ncol}, nz={nz}, dt={dt:g}, "
          f"jax backend {jax.default_backend()})")
    print(f"  copy of kessler_run_core with a pass counter: "
          f"{'matches the original bit for bit' if same else 'DIFFERS from the original'}")
    print(f"  passes per column: min {int(npass.min())}, mean {npass.mean():.2f}, max {mx}")
    print(f"  columns needing 1 pass: {ncol - n_multi} ({100.0 * (ncol - n_multi) / ncol:.1f}%), "
          f"more than 1: {n_multi} (max {mx_multi})")
    edges = [1, 2, 17, 33, 65, 129, 257]
    hist = []
    for lo, hi in zip(edges, edges[1:] + [None]):
        n = int(((npass >= lo) & (npass < hi)).sum()) if hi else int((npass >= lo).sum())
        label = f"{lo}" if hi == lo + 1 else (f"{lo}-{hi - 1}" if hi else f">={lo}")
        hist.append(f"{label}: {n}")
    print("  histogram: " + ", ".join(hist))
    print("  column-passes (work):")
    print(f"    per-column loops (sum)        {total:>12d}   1.00x")
    print(f"    split (ncol + max_multi*n)    {work_split:>12d}   {work_split / total:.2f}x")
    print(f"    vmap(while_loop) (max*ncol)   {work_vmap:>12d}   {work_vmap / total:.2f}x")

    # --- times of the device-resident bridge ---------------------------------------
    def bridge(arrs, n):
        o = kessler_run_bridge_device(
            ncol=n, nz=nz, dt=dt, lyr_surf=lyr_surf, lyr_toa=lyr_toa,
            rho=arrs[0], z=arrs[1], pk=arrs[2], theta=arrs[3], qv=arrs[4], qc=arrs[5],
            qr=arrs[6], precl=arrs[7], relhum=arrs[8], **scal)
        jax.block_until_ready(o[:6])

    t_all = _time(lambda: bridge(dev, ncol), repeat)
    dry = list(dev)
    dry[5] = jnp.zeros_like(dev[5])                # qc = 0
    dry[6] = jnp.zeros_like(dev[6])                # qr = 0: one pass per column
    t_dry = _time(lambda: bridge(dry, ncol), repeat)
    print(f"  kessler_run_bridge_device, s/call (min of {repeat}):")
    print(f"    all columns                   {t_all:.4e}   ({mx} loop iterations, "
          f"{t_all / max(mx, 1):.2e} s each)")
    print(f"    all columns, no rain or cloud {t_dry:.4e}   (1 iteration: split's first pass)")
    if n_multi:
        sel = np.nonzero(multi)[0]
        sub = [jnp.asarray(a[sel]) for a in host]
        t_sub = _time(lambda: bridge(sub, n_multi), repeat)
        print(f"    only the {n_multi} multi-pass columns {t_sub:.4e}   ({mx_multi} loop iterations, "
              f"{t_sub / max(mx_multi, 1):.2e} s each)")
        print(f"    split estimate (sum of the two above) {t_dry + t_sub:.4e}   "
              f"({t_all / (t_dry + t_sub):.2f}x vs all columns)")
    return same
