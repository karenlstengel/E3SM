"""
Restructured sub-cycle loops for kessler_run (prototypes).

kessler_jax.kessler_run.kessler_run_core runs each column's sub-cycle as a
lax.while_loop and vmaps it over the columns. Under vmap JAX runs one loop
until the slowest column is done, so every column does max(passes) passes
(most need 1, rainy ones 17-144). The variants here keep the physics of
kessler_run_core and change only how the loop is run:

  auto       pallas on a GPU, map otherwise: the fastest measured on each
             (driver_kessler_py --variants).
  vmap       the current structure, rebuilt from the pieces below. Used to
             check that the pieces reproduce kessler_run_core.
  split      one pass over all columns (dry columns finish there), then the
             vmapped loop on only the columns still running, gathered into a
             fixed-size buffer and scattered back. The buffer size is chosen
             at run time from a few static sizes with lax.switch, so there is
             one compile.
  map        lax.map over the columns: each column runs its own loop, as in
             the Fortran (do col / do while). Sequential; meant for the CPU.
  split_map  split, with lax.map over the gathered columns.
  pallas     the sub-cycle as one Pallas kernel, one column per program
             (kessler_run_pallas.py); GPU.
  pallas_lanes  the same with one column per thread (kessler_run_pallas.py).

kessler_run_core_variant has kessler_run_core's arguments and outputs, and
kessler_run_bridge_device_variant has kessler_run_bridge_device's (the
layout transposes included), each with a static `variant` argument.

Setup and finish (everything outside the sub-cycle) are copied from
kessler_run_core unchanged; the per-column body is copied from its _column.
Assumes at least two levels (EAMxx always has more).
"""

import functools
import os

import jax
import jax.numpy as jnp
from jax import lax

VARIANTS = ("auto", "vmap", "split", "map", "split_map", "pallas", "pallas_lanes")
F2X = 17.27


# ---------------------------------------------------------------------------
# Per-column pieces, from kessler_run_core's _column. A column's parameters
# p = (rho_surf, z_below_top, r, rhalf, pc, zc, pkc, dz, ok) and loop state
# s = (tc, dt0, pacc, th, qv, qc, qr, velqr); level arrays are (nlev,).
# ---------------------------------------------------------------------------
def _cfl(dt0, velqr, dz):
    v = velqr[:-1]
    mask = jnp.abs(v) > 1.0e-12
    cand = 0.8 * dz / jnp.where(mask, v, 1.0)
    cand = jnp.where(mask, cand, jnp.inf)
    return jnp.minimum(dt0, jnp.min(cand))


def _cond(s, p, dt):
    return (jnp.abs(dt - s[0]) > 1.0e-5) & p[8]


def _body(s, p, dt, lv, cpair, rhoqr, f5):
    tc, dt0, pacc, th, qvv, qcc, qrr, velqr = s
    rho_surf, z_below_top, r, rhalf, pc, zc, pkc, dz, ok = p
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
    qvs = pc * jnp.exp(F2X * (pt - 273.0) / (pt - 36.0))
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
    dt0 = _cfl(jnp.maximum(dt - tc, 0.0), velqr, dz)
    return (tc, dt0, pacc, th, qvv, qcc, qrr, velqr)


def _column_loop(s, p, consts):
    dt = consts[0]
    return lax.while_loop(lambda s: _cond(s, p, dt), lambda s: _body(s, p, *consts), s)


def _column_step(s, p, consts):
    """One iteration of the loop, if the column is still running (what the
    first iteration of the vmapped while_loop does to each column)."""
    new = _body(s, p, *consts)
    run = _cond(s, p, consts[0])
    return jax.tree.map(lambda a, b: jnp.where(run, a, b), new, s)


# ---------------------------------------------------------------------------
# Sub-cycle drivers. s and p are batched over axis 0 (columns).
# ---------------------------------------------------------------------------
def _loop_vmap(s, p, consts):
    return jax.vmap(_column_loop, in_axes=(0, 0, None))(s, p, consts)


def _loop_map(s, p, consts):
    return lax.map(lambda sp: _column_loop(sp[0], sp[1], consts), (s, p))


def _buffer_sizes(ncol):
    """Static sizes of the gathered buffer: ncol/64, /16, /4 and ncol by
    default, or ncol/d for each d in KESSLER_SPLIT_DIVISORS (comma-separated;
    read at trace time, for tuning). Each size is one branch of the
    lax.switch, compiled once: more sizes fit the running columns more
    tightly but take longer to compile (halving steps down to ncol/256: CPU
    split 3.6x -> 6x faster than the current code, GPU no faster, first
    compile ~1 s -> 7-8 s)."""
    env = os.environ.get("KESSLER_SPLIT_DIVISORS")
    divisors = [int(d) for d in env.split(",")] if env else [64, 16, 4, 1]
    return sorted({max(1, -(-ncol // d)) for d in divisors + [1]})


def _split(s, p, consts, ncol, loop):
    s1 = jax.vmap(_column_step, in_axes=(0, 0, None))(s, p, consts)
    running = jax.vmap(_cond, in_axes=(0, 0, None))(s1, p, consts[0])
    n = jnp.sum(running)
    sizes = _buffer_sizes(ncol)

    def branch(size):
        def f(s1):
            # Running columns first; the rest of the buffer points past the end
            idx = jnp.nonzero(running, size=size, fill_value=ncol)[0]
            take = lambda x: jnp.take(x, idx, axis=0, mode="clip")
            sg = jax.tree.map(take, s1)
            pg = jax.tree.map(take, p)
            pg = pg[:8] + (pg[8] & (idx < ncol),)        # padding entries do not run
            sf = loop(sg, pg, consts)
            return jax.tree.map(lambda full, part: full.at[idx].set(part, mode="drop"), s1, sf)
        return f

    which = sum((n > c).astype(jnp.int32) for c in sizes[:-1])  # smallest size >= n
    return lax.switch(which, [branch(c) for c in sizes], s1)


def _subcycle(variant, ncol, s, p, consts):
    if variant == "auto":
        variant = "pallas" if jax.default_backend() == "gpu" else "map"
    if variant == "vmap":
        return _loop_vmap(s, p, consts)
    if variant == "map":
        return _loop_map(s, p, consts)
    if variant == "split":
        return _split(s, p, consts, ncol, _loop_vmap)
    if variant == "split_map":
        return _split(s, p, consts, ncol, _loop_map)
    if variant == "pallas":
        from kessler_jax_experimental.kessler_run_pallas import subcycle_pallas
        return subcycle_pallas(s, p, consts)
    if variant == "pallas_lanes":
        from kessler_jax_experimental.kessler_run_pallas import subcycle_pallas_lanes
        return subcycle_pallas_lanes(s, p, consts)
    raise ValueError(f"unknown variant {variant!r}; one of {VARIANTS}")


# ---------------------------------------------------------------------------
# kessler_run_core with a choice of sub-cycle
# ---------------------------------------------------------------------------
@functools.partial(jax.jit, static_argnames=('variant', 'ncol', 'nz', 'lyr_surf', 'lyr_toa'))
def kessler_run_core_variant(variant, ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, pk,
                             theta, qv, qc, qr, precl, relhum, errflg, lv, pref, rhoqr):
    """kessler_run_core (kessler_jax/kessler_run.py), with the sub-cycle run
    as `variant`. Same arguments ((nz, ncol) arrays) and outputs."""
    # --- setup: as kessler_run_core ---
    s_ = lyr_surf - 1
    t = lyr_toa - 1
    step = 1 if t >= s_ else -1
    order = tuple(range(s_, t + step, step))
    nlev = len(order)
    assert nlev >= 2, "kessler_jax_experimental: needs at least two levels"
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

    r_o = 0.001 * rho_o
    rhalf_o = jnp.sqrt(rho_o[0:1, :] / rho_o)
    xk = cpair / rair
    pc_o = 3.8 / (pk_o ** xk * pref)
    qr0_o = jnp.maximum(qr_o, 0.0)
    velqr0_o = 36.34 * rhalf_o * (qr0_o * r_o) ** 0.1364
    f5 = 4093.0 * lv / cpair
    rho_surf_col = rho_o[0, :]
    z_below_top_col = z[t - step, :]

    # --- per-column start of the sub-cycle (kessler_run_core's _column) ---
    # Columns on axis 0 from here: (ncol, nlev)
    zc = z_o.T
    dz = zc[:, 1:] - zc[:, :-1]
    velqr_init = velqr0_o.T
    dt0_init = jax.vmap(_cfl, in_axes=(None, 0, 0))(dt, velqr_init, dz)
    bad_split = (dt > 0.0) & (dt0_init < 1.0e-12)
    col_ok = (dt > 0.0) & jnp.logical_not(bad_split)
    zeros = jnp.zeros((ncol,), dtype=jnp.float64)
    s0 = (zeros, dt0_init, zeros, theta_o.T, qv_o.T, qc_o.T, qr0_o.T, velqr_init)
    p = (rho_surf_col, z_below_top_col, r_o.T, rhalf_o.T, pc_o.T, zc, pk_o.T, dz, col_ok)
    consts = (dt, lv, cpair, rhoqr, f5)

    tc, dt0_f, pacc, th, qvv, qcc, qrr, velqr = _subcycle(variant, ncol, s0, p, consts)

    precl_n = pacc / dt
    pt = pk_o.T * th
    qvs = pc_o.T * jnp.exp(F2X * (pt - 273.0) / (pt - 36.0))
    relhum_n = (qvv / qvs * 100.0).T
    theta_n, qv_n, qc_n, qr_n = th.T, qvv.T, qcc.T, qrr.T

    # --- finish: as kessler_run_core ---
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
            errflg, lv, pref, rhoqr)


@functools.partial(jax.jit, static_argnames=('variant', 'ncol', 'nz', 'lyr_surf', 'lyr_toa', 'errflg'))
def kessler_run_bridge_device_variant(variant, ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, pk,
                                      theta, qv, qc, qr, precl, relhum, errflg, lv, pref, rhoqr):
    """bridge/kessler_run_bridge.py's kessler_run_bridge_device (device arrays
    in EAMxx's (ncol, nz) order in and out), on kessler_run_core_variant."""
    T = jnp.transpose
    out = kessler_run_core_variant(variant, ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair,
                                   T(rho), T(z), T(pk), T(theta), T(qv), T(qc), T(qr),
                                   precl, T(relhum), errflg, lv, pref, rhoqr)
    theta, qv, qc, qr, precl, relhum, errflg, lv, pref, rhoqr = out
    return T(theta), T(qv), T(qc), T(qr), precl, T(relhum), errflg, lv, pref, rhoqr
