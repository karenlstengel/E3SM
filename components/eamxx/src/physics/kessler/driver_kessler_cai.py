"""
Zero-copy device path for driver_kessler_py.cpp (--backend cai).

EAMxx's py_backend=device wraps device pointers in plain NumPy arrays
(create_py_field<Device>), which NumPy and JAX then read as host memory.
This module tries the alternative: the driver passes each field as an
object exposing __cuda_array_interface__ (CAI, version 3), jnp.asarray
imports it on the GPU without a host copy, the device-resident bridges
(contract 2) run on it, and the outputs go back to the driver as JAX
arrays, whose own __cuda_array_interface__ gives the driver the device
pointer to copy from. No data crosses to the host.

run() and update() do the same computation as kessler.run() and
kessler.update(), with the same jitted bridge functions and the same
module state (kessler.latvap/pref/rhoqr/gravit/errflg), so results
should match the host path bit for bit. They return their outputs
instead of writing them back in place: the inputs are CAI objects, not
writable arrays.

Inputs must be C-contiguous (the driver packs strided tracers into
contiguous scratch first). Every call blocks until JAX is done, so the
driver can copy the outputs as soon as the call returns.
"""

import time

import jax
import jax.numpy as jnp

import kessler as K
from bridge.kessler_run_bridge import kessler_run_bridge_device
from bridge.kessler_update_timestep_init_bridge import kessler_update_timestep_init_bridge_device
from bridge.kessler_update_run_bridge import kessler_update_run_bridge_device
from bridge.kessler_update_timestep_final_bridge import kessler_update_timestep_final_bridge_device


# Seconds spent in the last run()/update() call, by part; the driver adds
# these to its per-stage times (import = jnp.asarray of the CAI inputs,
# compute = the jitted bridges until block_until_ready, rest = scalar
# fetches and bookkeeping).
last_times = {}


def _dev(x):
    # Same as the bridges' to_device(); on a CAI object jnp.asarray imports
    # the device buffer instead of reading host memory.
    return jnp.asarray(x, dtype=jnp.float64)


def run(ncol, nz, dt, lyr_surf, lyr_toa, rho, zm, pk, theta, qv, qc, qr, precl, relhum):
    """kessler.run() on device buffers. Returns (theta, qv, qc, qr, precl, relhum)."""
    t0 = time.perf_counter()
    rho, zm, pk, theta, qv, qc, qr, precl, relhum = (
        _dev(x) for x in (rho, zm, pk, theta, qv, qc, qr, precl, relhum))
    t1 = time.perf_counter()
    out = kessler_run_bridge_device(
        ncol=ncol, nz=nz, dt=dt, lyr_surf=lyr_surf, lyr_toa=lyr_toa,
        cpair=K.cpair, rair=K.rair,
        rho=rho, z=zm, pk=pk, theta=theta, qv=qv, qc=qc, qr=qr, precl=precl, relhum=relhum,
        errflg=K.errflg, lv=K.latvap, pref=K.pref, rhoqr=K.rhoqr)
    theta_o, qv_o, qc_o, qr_o, precl_o, relhum_o, errflg_o, K.latvap, K.pref, K.rhoqr = out
    arrays = jax.block_until_ready((theta_o, qv_o, qc_o, qr_o, precl_o, relhum_o))
    t2 = time.perf_counter()
    K.errflg = int(errflg_o)   # as kessler_run_bridge does (.item())
    last_times.clear()
    last_times.update({"import": t1 - t0, "compute": t2 - t1, "rest": time.perf_counter() - t2})
    return arrays


def update(ncol, nz, dt, zm, pk, theta, phis, temp, temp_prev, temp_tend, st_energy):
    """kessler.update() on device buffers. Returns (temp, temp_prev, temp_tend, st_energy)."""
    t0 = time.perf_counter()
    zm, pk, theta, phis, temp = (_dev(x) for x in (zm, pk, theta, phis, temp))
    t1 = time.perf_counter()
    zeros_d = jnp.zeros((ncol, nz), dtype=jnp.float64)
    temp_prev_d, ttend_t_d, e1 = kessler_update_timestep_init_bridge_device(
        temp=temp, temp_prev=zeros_d, ttend_t=zeros_d, errflg=0)
    ttend_t_d, e2 = kessler_update_run_bridge_device(
        nz=nz, ncol=ncol, dt=dt, theta=theta, exner=pk,
        temp_prev=temp_prev_d, ttend_t=ttend_t_d, errflg=0)
    temp_new_d = temp_prev_d + ttend_t_d * dt
    st_energy_d, e3, K.gravit = kessler_update_timestep_final_bridge_device(
        nz=nz, cpair=K.cpair, temp=temp_new_d, zm=zm, phis=phis,
        st_energy=zeros_d, errflg=0, gravit=K.gravit)
    arrays = jax.block_until_ready((temp_new_d, temp_prev_d, ttend_t_d, st_energy_d))
    t2 = time.perf_counter()
    K.errflg = max(int(e1), int(e2), int(e3))
    K.errmsg = "" if K.errflg == 0 else "kessler_update: bad time splitting"
    last_times.clear()
    last_times.update({"import": t1 - t0, "compute": t2 - t1, "rest": time.perf_counter() - t2})
    return arrays
