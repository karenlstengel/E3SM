"""
Compare the sub-cycle variants of kessler_jax_experimental with the current
kessler_run (driver_kessler_py.cpp --variants).

For one input state, runs bridge/kessler_run_bridge.py's
kessler_run_bridge_device (the current code) and
kessler_jax_experimental's kessler_run_bridge_device_variant for each
variant, on the same device arrays, and prints for each variant the time
per call (min, after a compile call; waiting on the device, without
copying the outputs to the host), the same for the columns with no rain or
cloud (one pass each: the fixed cost around the loop), and how its outputs
differ from the current code's: bit for bit, or the number of differing
values and the largest relative difference per field.
"""

import time
import traceback

import numpy as np
import jax
import jax.numpy as jnp

import kessler as K
from bridge.kessler_run_bridge import kessler_run_bridge_device
from kessler_jax_experimental.kessler_run_variants import kessler_run_bridge_device_variant

FIELDS = ("theta", "qv", "qc", "qr", "precl", "relhum", "errflg")


def _time(f, repeat):
    t0 = time.perf_counter()
    out = f()                                      # compiles
    first = time.perf_counter() - t0
    best = float("inf")
    for _ in range(repeat):
        t0 = time.perf_counter()
        f()
        best = min(best, time.perf_counter() - t0)
    return out, first, best


def _diff(ref, out):
    """'BFB', or per-field 'n differ (max rel)' for the fields that differ."""
    parts = []
    for name, a, b in zip(FIELDS, ref, out):
        a, b = np.asarray(a), np.asarray(b)
        d = a != b
        if d.any():
            rel = np.abs(a[d] - b[d]) / np.maximum(np.abs(a[d]), 1e-300)
            parts.append(f"{name} {int(d.sum())} differ (max rel {rel.max():.1e})")
    return "BFB" if not parts else "; ".join(parts)


def compare(ncol, nz, dt, lyr_surf, lyr_toa, rho, zm, pk, theta, qv, qc, qr, precl, relhum,
            variants=("vmap", "split", "map", "split_map", "pallas", "pallas_lanes"), repeat=10):
    """Arrays: host NumPy, EAMxx's (ncol, nz) order (as kessler.run() gets
    them with py_backend=host); kessler.init() must have run. Returns False
    if a variant failed to run."""
    dev = [jnp.asarray(np.ascontiguousarray(a, dtype=np.float64))
           for a in (rho, zm, pk, theta, qv, qc, qr, precl, relhum)]
    kw = dict(ncol=ncol, nz=nz, dt=dt, lyr_surf=lyr_surf, lyr_toa=lyr_toa, cpair=K.cpair, rair=K.rair,
              rho=dev[0], z=dev[1], pk=dev[2], theta=dev[3], qv=dev[4], qc=dev[5], qr=dev[6],
              precl=dev[7], relhum=dev[8], errflg=0, lv=K.latvap, pref=K.pref, rhoqr=K.rhoqr)

    dry = dict(kw, qc=jnp.zeros_like(dev[5]), qr=jnp.zeros_like(dev[6]))   # one pass per column

    def ready(o):                     # waits on the GPU; no copy to the host
        return jax.block_until_ready(o[:7])

    def timed(f):
        out, first, t = _time(lambda: ready(f(**kw)), repeat)
        _, _, t_dry = _time(lambda: ready(f(**dry)), repeat)
        return jax.device_get(out), first, t, t_dry

    ref, first, t_ref, d_ref = timed(kessler_run_bridge_device)
    print(f"\n--- kessler_run sub-cycle variants (ncol={ncol}, nz={nz}, jax backend "
          f"{jax.default_backend()}), s/call, min of {repeat}; 'no rain': same columns with "
          f"qc = qr = 0 (one pass each)")
    print(f"  {'':<12} {'s/call':>10}  {'speedup':>7}  {'no rain':>10}")
    print(f"  {'current':<12} {t_ref:.4e}  {1.0:6.2f}x  {d_ref:.4e}  (first call {first:.2f} s)")
    ok = True
    for v in variants:
        if v.startswith("pallas") and jax.default_backend() != "gpu":
            print(f"  {v:<12} skipped (GPU only)")
            continue
        try:
            f = lambda **a: kessler_run_bridge_device_variant(v, **a)
            out, first, t, t_dry = timed(f)
        except Exception as e:
            ok = False
            print(f"  {v:<12} FAILED: {type(e).__name__}: {str(e).splitlines()[0][:300]}")
            traceback.print_exc()
            continue
        print(f"  {v:<12} {t:.4e}  {t_ref / t:6.2f}x  {t_dry:.4e}  (first call {first:.2f} s)  {_diff(ref, out)}")
    return ok
