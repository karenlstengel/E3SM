"""
kessler_run's sub-cycle as one Pallas kernel (GPU; prototypes).

Called from kessler_run_variants.py with the same loop state and per-column
parameters as the other variants. The whole sub-cycle of every column runs
in one kernel launch, and each loop iteration's exit test stays on the GPU.
Two layouts:

subcycle_pallas (variant "pallas"): one column per program, its levels
spread over the program's threads, the whole column state in registers.
This is the C++ D1 structure (one team per column): each column stops after
its own passes, and the passes of one column run level-parallel. Reading
level k+1 needs a shift, and the Triton backend cannot slice values. By
default the shift goes through memory: store, barrier (debug_barrier, a
__syncthreads), gather-load at k+1, barrier. The alternative (shift=
"select") is a masked sum over a (levels x levels) selector, one nonzero
term plus zeros, which is exact but costs a cross-thread reduction of
levels^2 values per pass (about 1.6 ms of 2 ms at ne16).

subcycle_pallas_lanes (variant "pallas_lanes"): a block of columns per
program, one column per lane, each pass looping over the levels in order
as the Fortran does. No shift needed (each thread reads its own column's
level k+1 from memory, and updating in ascending order keeps it the old
value), but there is only one thread per column: at ne16 too few threads to
hide the latency of the float64 pow/exp chain.

Restrictions of Pallas's Triton backend (A100, JAX 0.6.2) these work
around: no slicing of values, block shapes in powers of two (levels padded
to 64), and no float64 dot.
"""
import functools
import os

import jax
import jax.numpy as jnp
from jax import lax
from jax.experimental import pallas as pl
from jax.experimental.pallas import triton as plgpu

F2X = 17.27


def _cparams(num_warps):
    return plgpu.CompilerParams(num_warps=num_warps, num_stages=1)


# ---------------------------------------------------------------------------
# Layout A: one column per program, state in registers
# ---------------------------------------------------------------------------
NSC = 16   # per-column scalars: dt, lv, cpair, rhoqr, f5, rho_surf, z_below_top, tc, dt0, pacc, ok


def _kernel_cols(nlev, nlp, shift, r_ref, rhalf_ref, pc_ref, zc_ref, pk_ref, dz_ref,
                 th0_ref, qv0_ref, qc0_ref, qr0_ref, vq0_ref, cs_ref,
                 th_ref, qv_ref, qc_ref, qr_ref, vq_ref, co_ref, fl_ref):
    r, rhalf, pc, zc, pkc, dz = (x[0, :] for x in (r_ref, rhalf_ref, pc_ref, zc_ref, pk_ref, dz_ref))
    cs = cs_ref[0, :]
    j = jnp.arange(NSC)

    def scalar(i):
        return jnp.sum(jnp.where(j == i, cs, 0.0))

    dt, lv, cpair, rhoqr, f5, rho_surf, z_below_top, tc0, dt00, pacc0, okf = (
        scalar(i) for i in range(11))
    ok = okf > 0.5
    k = jnp.arange(nlp)
    below_top = k < nlev - 1
    top = k == nlev - 1
    surf = k == 0
    if shift == "memory":
        # x[k + 1] at k through memory: store the column, wait for every
        # thread of the program, read it back one level up, wait again
        # (so the next store cannot overtake a slow reader). fl_ref is an
        # output used as scratch.
        k_up = jnp.minimum(k + 1, nlp - 1)

        def shift_up(x):
            fl_ref[0, :] = x
            plgpu.debug_barrier()
            y = pl.load(fl_ref, (0, k_up))
            plgpu.debug_barrier()
            return y
    else:
        up = (k[:, None] + 1) == k[None, :]       # up[k, i]: level i is k + 1

        def shift_up(x):                           # masked sum: one term plus zeros (exact)
            return jnp.sum(jnp.where(up, x[None, :], 0.0), axis=1)

    def cond(c):
        return (jnp.abs(dt - c[0]) > 1.0e-5) & ok

    def body(c):
        # kessler_run_core's _column body, whole column at once
        tc, dt0, pacc, th, qvv, qcc, qrr, velqr = c
        precl_inst = jnp.sum(jnp.where(surf, rho_surf * qrr * velqr / rhoqr, 0.0))
        pacc = pacc + precl_inst * dt0
        flux = r * qrr * velqr
        sed_body = dt0 * (shift_up(flux) - flux) / (r * dz)
        sed_top = -dt0 * qrr * velqr / (0.5 * (zc - z_below_top))
        sed = jnp.where(top, sed_top, sed_body)
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
        mask = below_top & (jnp.abs(velqr) > 1.0e-12)
        cand = 0.8 * dz / jnp.where(mask, velqr, 1.0)
        dt0 = jnp.minimum(jnp.maximum(dt - tc, 0.0), jnp.min(jnp.where(mask, cand, jnp.inf)))
        return (tc, dt0, pacc, th, qvv, qcc, qrr, velqr)

    c0 = (tc0, dt00, pacc0, th0_ref[0, :], qv0_ref[0, :], qc0_ref[0, :], qr0_ref[0, :], vq0_ref[0, :])
    tc, dt0, pacc, th, qvv, qcc, qrr, velqr = lax.while_loop(cond, body, c0)
    th_ref[0, :] = th
    qv_ref[0, :] = qvv
    qc_ref[0, :] = qcc
    qr_ref[0, :] = qrr
    vq_ref[0, :] = velqr
    co_ref[0, :] = jnp.where(j == 0, tc, jnp.where(j == 1, dt0, jnp.where(j == 2, pacc, 0.0)))


def subcycle_pallas(s, p, consts, num_warps=None, shift=None):
    """The sub-cycle of every column, as kessler_run_variants' loops take and
    return it: s = (tc, dt0, pacc, th, qv, qc, qr, velqr), p = (rho_surf,
    z_below_top, r, rhalf, pc, zc, pkc, dz, ok), batched over columns
    (axis 0); consts = (dt, lv, cpair, rhoqr, f5). Layout A."""
    # For tuning (read at trace time): KESSLER_PALLAS_NUM_WARPS, and
    # KESSLER_PALLAS_SHIFT = memory (store, barrier, load) or select (masked sum)
    if num_warps is None:
        num_warps = int(os.environ.get("KESSLER_PALLAS_NUM_WARPS", "2"))
    if shift is None:
        shift = os.environ.get("KESSLER_PALLAS_SHIFT", "memory")
    tc, dt0, pacc, th, qv, qc, qr, vq = s
    rho_surf, z_below_top, r, rhalf, pc, zc, pk, dz, ok = p
    ncol, nlev = th.shape
    f64 = jnp.float64
    nlp = 1 << (nlev - 1).bit_length()

    def lev(x):           # (ncol, n) -> (ncol, nlp), padding 1.0 (never used)
        return jnp.pad(x, ((0, 0), (0, nlp - x.shape[1])), constant_values=1.0)

    per_col = [jnp.broadcast_to(jnp.asarray(c, f64), (ncol,)) for c in consts] + [
        rho_surf, z_below_top, tc, dt0, pacc, ok.astype(f64)]
    cs = jnp.stack(per_col + [jnp.zeros((ncol,), f64)] * (NSC - len(per_col)), axis=1)

    lvl = pl.BlockSpec((1, nlp), lambda i: (i, 0))
    scs = pl.BlockSpec((1, NSC), lambda i: (i, 0))
    outs = pl.pallas_call(
        functools.partial(_kernel_cols, nlev, nlp, shift),
        grid=(ncol,),
        in_specs=[lvl] * 11 + [scs],
        out_specs=[lvl] * 5 + [scs, lvl],
        out_shape=[jax.ShapeDtypeStruct((ncol, nlp), f64)] * 5
                  + [jax.ShapeDtypeStruct((ncol, NSC), f64), jax.ShapeDtypeStruct((ncol, nlp), f64)],
        compiler_params=_cparams(num_warps),
    )(lev(r), lev(rhalf), lev(pc), lev(zc), lev(pk), lev(dz),
      lev(th), lev(qv), lev(qc), lev(qr), lev(vq), cs)
    th_o, qv_o, qc_o, qr_o, vq_o, co, _ = outs
    back = lambda x: x[:, :nlev]
    return (co[:, 0], co[:, 1], co[:, 2],
            back(th_o), back(qv_o), back(qc_o), back(qr_o), back(vq_o))


# ---------------------------------------------------------------------------
# Layout B: columns on lanes, loop over levels
# ---------------------------------------------------------------------------


def _kernel_lanes(nlev, r_ref, rhalf_ref, pc_ref, zc_ref, pk_ref, dz_ref,
            th0_ref, qv0_ref, qc0_ref, qr0_ref, vq0_ref, cs_ref, cst_ref,
            th_ref, qv_ref, qc_ref, qr_ref, vq_ref, co_ref):
    def rd(ref, k):
        return ref[pl.ds(k, 1), :]

    def wr(ref, k, v):
        ref[pl.ds(k, 1), :] = v

    # The outputs hold the state: start them from the inputs
    def init(k, c):
        wr(th_ref, k, rd(th0_ref, k))
        wr(qv_ref, k, rd(qv0_ref, k))
        wr(qc_ref, k, rd(qc0_ref, k))
        wr(qr_ref, k, rd(qr0_ref, k))
        wr(vq_ref, k, rd(vq0_ref, k))
        return c
    lax.fori_loop(0, nlev, init, 0)

    dt, lv, cpair, rhoqr, f5 = (rd(cst_ref, i) for i in range(5))
    rho_surf = rd(cs_ref, 0)
    z_below_top = rd(cs_ref, 1)
    ok = rd(cs_ref, 5) > 0.5

    def running(tc):
        return (jnp.abs(dt - tc) > 1.0e-5) & ok

    def cond(c):
        return jnp.max(jnp.where(running(c[0]), 1.0, 0.0)) > 0.5

    def body(c):
        tc, dt0, pacc = c
        act = running(tc)
        precl_inst = rho_surf * rd(qr_ref, 0) * rd(vq_ref, 0) / rhoqr
        pacc_n = pacc + precl_inst * dt0

        def level(k, sed):
            # kessler_run_core's per-level updates, for level k of each column
            r = rd(r_ref, k)
            th = rd(th_ref, k)
            qvv = rd(qv_ref, k)
            qcc = rd(qc_ref, k)
            qrr = rd(qr_ref, k)
            pc = rd(pc_ref, k)
            pkc = rd(pk_ref, k)
            qrprod = qcc - (qcc - dt0 * jnp.maximum(0.001 * (qcc - 0.001), 0.0)) \
                / (1.0 + dt0 * 2.2 * qrr ** 0.875)
            qcc2 = jnp.maximum(qcc - qrprod, 0.0)
            qrr2 = jnp.maximum(qrr + qrprod + sed, 0.0)
            pt = pkc * th
            qvs = pc * jnp.exp(F2X * (pt - 273.0) / (pt - 36.0))
            prod = (qvv - qvs) / (1.0 + qvs * f5 / (pt - 36.0) ** 2)
            rqr = r * qrr2
            ern = jnp.minimum(
                jnp.minimum(
                    dt0 * (((1.6 + 124.9 * rqr ** 0.2046) * rqr ** 0.525)
                           / (2550000.0 * pc / (3.8 * qvs) + 540000.0))
                    * (jnp.maximum(qvs - qvv, 0.0) / (r * qvs)),
                    jnp.maximum(-prod - qcc2, 0.0)),
                qrr2)
            th2 = th + lv / (cpair * pkc) * (jnp.maximum(prod, -qcc2) - ern)
            qvv2 = jnp.maximum(qvv - jnp.maximum(prod, -qcc2) + ern, 0.0)
            qcc3 = qcc2 + jnp.maximum(prod, -qcc2)
            qrr3 = jnp.maximum(qrr2 - ern, 0.0)
            vq2 = 36.34 * rd(rhalf_ref, k) * (qrr3 * r) ** 0.1364
            wr(th_ref, k, jnp.where(act, th2, th))
            wr(qv_ref, k, jnp.where(act, qvv2, qvv))
            wr(qc_ref, k, jnp.where(act, qcc3, qcc))
            wr(qr_ref, k, jnp.where(act, qrr3, qrr))
            wr(vq_ref, k, jnp.where(act, vq2, rd(vq_ref, k)))
            return vq2

        def below_top(k, cmin):
            # sedimentation from the old level k+1 (not yet updated this pass)
            r = rd(r_ref, k)
            dz = rd(dz_ref, k)
            flux = r * rd(qr_ref, k) * rd(vq_ref, k)
            flux_up = rd(r_ref, k + 1) * rd(qr_ref, k + 1) * rd(vq_ref, k + 1)
            sed = dt0 * (flux_up - flux) / (r * dz)
            vq2 = level(k, sed)
            # CFL candidate from the new fall speed (kessler_run_core's _cfl)
            mask = jnp.abs(vq2) > 1.0e-12
            cand = 0.8 * dz / jnp.where(mask, vq2, 1.0)
            return jnp.minimum(cmin, jnp.where(mask, cand, jnp.inf))

        cmin = lax.fori_loop(0, nlev - 1, below_top, jnp.full_like(dt, jnp.inf))
        kt = nlev - 1
        sed_top = -dt0 * rd(qr_ref, kt) * rd(vq_ref, kt) / (0.5 * (rd(zc_ref, kt) - z_below_top))
        level(kt, sed_top)
        tc_n = tc + dt0
        dt0_n = jnp.minimum(jnp.maximum(dt - tc_n, 0.0), cmin)
        return (jnp.where(act, tc_n, tc), jnp.where(act, dt0_n, dt0), jnp.where(act, pacc_n, pacc))

    tc, dt0, pacc = lax.while_loop(cond, body, (rd(cs_ref, 2), rd(cs_ref, 3), rd(cs_ref, 4)))
    wr(co_ref, 0, tc)
    wr(co_ref, 1, dt0)
    wr(co_ref, 2, pacc)


def subcycle_pallas_lanes(s, p, consts, block_cols=32, num_warps=1, pad_levels=True):
    """Layout B (see the top of this file); arguments as subcycle_pallas."""
    tc, dt0, pacc, th, qv, qc, qr, vq = s
    rho_surf, z_below_top, r, rhalf, pc, zc, pk, dz, ok = p
    ncol, nlev = th.shape
    f64 = jnp.float64
    bc = block_cols
    ncolp = -(-ncol // bc) * bc
    nlp = (1 << (nlev - 1).bit_length()) if pad_levels else nlev

    def lev(x):           # (ncol, n) -> (nlp, ncolp), padding 1.0 (never used)
        x = x.T
        return jnp.pad(x, ((0, nlp - x.shape[0]), (0, ncolp - ncol)), constant_values=1.0)

    def row(x):           # (ncol,) -> (1, ncolp); padding columns get ok = 0
        return jnp.pad(x.astype(f64)[None, :], ((0, 0), (0, ncolp - ncol)))

    zero = jnp.zeros((1, ncolp), f64)
    cs = jnp.concatenate([row(rho_surf), row(z_below_top), row(tc), row(dt0), row(pacc),
                          row(ok), zero, zero])
    cst = jnp.broadcast_to(jnp.stack([jnp.asarray(c, f64) for c in consts]
                                     + [jnp.zeros((), f64)] * 3)[:, None], (8, ncolp))

    lvl = pl.BlockSpec((nlp, bc), lambda i: (0, i))
    r8 = pl.BlockSpec((8, bc), lambda i: (0, i))
    outs = pl.pallas_call(
        functools.partial(_kernel_lanes, nlev),
        grid=(ncolp // bc,),
        in_specs=[lvl] * 11 + [r8, r8],
        out_specs=[lvl] * 5 + [r8],
        out_shape=[jax.ShapeDtypeStruct((nlp, ncolp), f64)] * 5
                  + [jax.ShapeDtypeStruct((8, ncolp), f64)],
        compiler_params=_cparams(num_warps),
    )(lev(r), lev(rhalf), lev(pc), lev(zc), lev(pk), lev(dz),
      lev(th), lev(qv), lev(qc), lev(qr), lev(vq), cs, cst)
    th_o, qv_o, qc_o, qr_o, vq_o, co = outs
    back = lambda x: x[:nlev, :ncol].T
    return (co[0, :ncol], co[1, :ncol], co[2, :ncol],
            back(th_o), back(qv_o), back(qc_o), back(qr_o), back(vq_o))
