"""
Bridge: kessler_update_timestep_init
Strategy: PATH C device-side layout conversion — pure H2D/D2H,
in-jit axis reversal fused by XLA; public contract unchanged

MODULE_VARIABLE_POLICY (INOUT Pattern):
- MODULE vars passed as function parameters (INOUT)
- Driver manages MODULE state explicitly
- No hidden state dictionary
- INIT: Returns updated MODULE vars to driver
- COMPUTE: Returns updated MODULE vars to driver
"""

"""
Layout Conversion (PATH C: device-side)

STRATEGY:
- Ship host arrays to the device AS-IS (pure H2D, no host permute)
- Reverse axes INSIDE the jitted device wrapper — XLA fuses the
  permute into the kernel (HBM bandwidth, not host memcpy)
- Pure D2H on the way back: outputs already (ncol, nz) C-order,
  fetched in ONE batched jax.device_get (single stream sync)

Two bridge contracts per array-bearing procedure:
  1. {proc}_bridge        — host-facing: standard C-order NumPy arrays
                            in Fortran (ncol, nz) index order in/out
                            (H2D + jitted device wrapper + ONE batched D2H)
  2. {proc}_bridge_device — device-resident: the jitted wrapper itself;
                            jax.Arrays already on the device in/out, no
                            transfers. For callers that keep the model
                            state on the device across the per-step call
                            sequence (fetch once at the end).
"""

import functools

import numpy as np
import jax
import jax.numpy as jnp

def to_device(arr):
    """Host → device unchanged (pure H2D; no host-side permute)"""
    return jnp.asarray(arr, dtype=jnp.float64)

def to_host(arr):
    """Device → host NumPy for a SINGLE array (blocking D2H).

    The bridge body does NOT use this per-array path: outputs are
    fetched together in one batched jax.device_get (one sync total
    instead of one per array). Kept for tests and interactive use."""
    return np.array(arr, dtype=np.float64)

from kessler_jax.kessler_update_timestep_init import kessler_update_timestep_init_core

# Device-side layout wrapper (PATH C): arrays arrive in Fortran
# (ncol, nz) index order; axes are reversed INSIDE the jitted region
# so the permute runs on the device and XLA can fuse it into the
# kernel. CHARACTER args never cross this boundary (not JAX types).
# PUBLIC device-resident entry — see its docstring.
@functools.partial(jax.jit, static_argnames=('errflg',))
def kessler_update_timestep_init_bridge_device(temp, temp_prev, ttend_t, errflg):
    """
    DEVICE-RESIDENT bridge entry for kessler_update_timestep_init (contract 2, no transfers).

    Same arguments as kessler_update_timestep_init_bridge minus CHARACTER params, in the
    same order; array arguments are jax.Arrays ALREADY ON THE DEVICE in
    Fortran (ncol, nz) index order (the layout to_device(host) produces);
    integer/logical scalars are static Python values; MODULE vars are
    passed and returned (INOUT pattern). Returns the written non-CHARACTER
    params in argument order, then MODULE vars — as device arrays, left
    on the device (no device_get). Chain these entries across the
    per-step procedure sequence and fetch once at the end; the outputs
    are bit-identical to kessler_update_timestep_init_bridge's (same jitted code — the
    host bridge is this function wrapped in to_device / device_get).
    """
    # In-jit input reversal: reverse ALL axes on rank>=2
    # (covers rank 3+; rank 0/1 pass through)
    temp = jnp.transpose(temp)
    temp_prev = jnp.transpose(temp_prev)
    ttend_t = jnp.transpose(ttend_t)
    temp_prev, ttend_t, errflg = kessler_update_timestep_init_core(temp, temp_prev, ttend_t, errflg)
    # In-jit output reversal: back to Fortran (ncol, nz) order
    temp_prev = jnp.transpose(temp_prev)
    ttend_t = jnp.transpose(ttend_t)
    return temp_prev, ttend_t, errflg

_kessler_update_timestep_init_device = kessler_update_timestep_init_bridge_device  # alias (private name kept for compatibility)

def kessler_update_timestep_init_bridge(temp, temp_prev, ttend_t, errmsg, errflg):
    """
    Production HOST-FACING bridge for kessler_update_timestep_init (contract 1).

    Device-resident callers use kessler_update_timestep_init_bridge_device instead
    (contract 2: device arrays in/out, no transfers); this function is
    exactly to_device -> kessler_update_timestep_init_bridge_device -> batched device_get.

    Strategy (PATH C): ship arrays to the device unchanged (pure H2D),
    reverse axes inside the jitted device wrapper (fused by XLA), call
    kessler_update_timestep_init_core, reverse back in-jit, pure D2H. Inputs and outputs
    stay standard C-order NumPy in Fortran (ncol, nz) index order.

    CHARACTER args stay host-side: this path calls the jitted core
    directly, so wrapper-level string handling does not run here;
    CHARACTER inputs/outputs pass through the bridge unchanged.
    """
    # Ship arrays to the device unchanged (pure H2D, no host permute)
    temp_dev = to_device(temp)
    temp_prev_dev = to_device(temp_prev)
    ttend_t_dev = to_device(ttend_t)

    # Compute — in-jit layout conversion + jitted core
    temp_prev_out, ttend_t_out, errflg_out = kessler_update_timestep_init_bridge_device(temp=temp_dev, temp_prev=temp_prev_dev, ttend_t=ttend_t_dev, errflg=errflg)

    # Convert outputs — ONE batched D2H (single sync); layout already
    # Fortran (ncol, nz) C-order
    temp_prev_host, ttend_t_host, errflg_host = jax.device_get((temp_prev_out, ttend_t_out, errflg_out))
    temp_prev_fortran = np.asarray(temp_prev_host, dtype=np.float64)
    ttend_t_fortran = np.asarray(ttend_t_host, dtype=np.float64)
    errflg_fortran = np.asarray(errflg_host).item()

    return temp_prev_fortran, ttend_t_fortran, errmsg, errflg_fortran
