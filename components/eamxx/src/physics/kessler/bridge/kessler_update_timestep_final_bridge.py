"""
Bridge: kessler_update_timestep_final
Strategy: PATH C device-side layout conversion — pure H2D/D2H,
in-jit axis reversal fused by XLA; public contract unchanged

MODULE: kessler_update
MODULE vars: gravit
Procedure type: COMPUTE

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

from kessler_jax.kessler_update_timestep_final import kessler_update_timestep_final_core

# Device-side layout wrapper (PATH C): arrays arrive in Fortran
# (ncol, nz) index order; axes are reversed INSIDE the jitted region
# so the permute runs on the device and XLA can fuse it into the
# kernel. CHARACTER args never cross this boundary (not JAX types).
# PUBLIC device-resident entry — see its docstring.
@functools.partial(jax.jit, static_argnames=('nz', 'errflg',))
def kessler_update_timestep_final_bridge_device(nz, cpair, temp, zm, phis, st_energy, errflg, gravit):
    """
    DEVICE-RESIDENT bridge entry for kessler_update_timestep_final (contract 2, no transfers).

    Same arguments as kessler_update_timestep_final_bridge minus CHARACTER params, in the
    same order; array arguments are jax.Arrays ALREADY ON THE DEVICE in
    Fortran (ncol, nz) index order (the layout to_device(host) produces);
    integer/logical scalars are static Python values; MODULE vars are
    passed and returned (INOUT pattern). Returns the written non-CHARACTER
    params in argument order, then MODULE vars — as device arrays, left
    on the device (no device_get). Chain these entries across the
    per-step procedure sequence and fetch once at the end; the outputs
    are bit-identical to kessler_update_timestep_final_bridge's (same jitted code — the
    host bridge is this function wrapped in to_device / device_get).
    """
    # In-jit input reversal: reverse ALL axes on rank>=2
    # (covers rank 3+; rank 0/1 pass through). cpair is a scalar,
    # so it passes through unchanged.
    temp = jnp.transpose(temp)
    zm = jnp.transpose(zm)
    st_energy = jnp.transpose(st_energy)
    st_energy, errflg, gravit = kessler_update_timestep_final_core(nz, cpair, temp, zm, phis, st_energy, errflg, gravit)
    # In-jit output reversal: back to Fortran (ncol, nz) order
    st_energy = jnp.transpose(st_energy)
    return st_energy, errflg, gravit

_kessler_update_timestep_final_device = kessler_update_timestep_final_bridge_device  # alias (private name kept for compatibility)

def kessler_update_timestep_final_bridge(nz, cpair, temp, zm, phis, st_energy, errflg, errmsg, gravit):
    """
    Production HOST-FACING bridge for kessler_update_timestep_final (contract 1).

    Device-resident callers use kessler_update_timestep_final_bridge_device instead
    (contract 2: device arrays in/out, no transfers); this function is
    exactly to_device -> kessler_update_timestep_final_bridge_device -> batched device_get.

    Strategy (PATH C): ship arrays to the device unchanged (pure H2D),
    reverse axes inside the jitted device wrapper (fused by XLA), call
    kessler_update_timestep_final_core, reverse back in-jit, pure D2H. Inputs and outputs
    stay standard C-order NumPy in Fortran (ncol, nz) index order.

    CHARACTER args stay host-side: this path calls the jitted core
    directly, so wrapper-level string handling does not run here;
    CHARACTER inputs/outputs pass through the bridge unchanged.

    MODULE VARIABLES (from kessler_update):
      - gravit (INOUT)
    
    Procedure type: COMPUTE
    Pattern: MODULE vars passed as INOUT parameters
    """
    # Ship arrays to the device unchanged (pure H2D, no host permute).
    # cpair is a spatially-uniform scalar, passed through as-is.
    temp_dev = to_device(temp)
    zm_dev = to_device(zm)
    phis_dev = to_device(phis)
    st_energy_dev = to_device(st_energy)

    # Compute — in-jit layout conversion + jitted core
    st_energy_out, errflg_out, gravit = kessler_update_timestep_final_bridge_device(nz=nz, cpair=cpair, temp=temp_dev, zm=zm_dev, phis=phis_dev, st_energy=st_energy_dev, errflg=errflg, gravit=gravit)

    # Convert outputs — ONE batched D2H (single sync); layout already
    # Fortran (ncol, nz) C-order
    st_energy_host, errflg_host = jax.device_get((st_energy_out, errflg_out))
    st_energy_fortran = np.asarray(st_energy_host, dtype=np.float64)
    errflg_fortran = np.asarray(errflg_host).item()

    return st_energy_fortran, errflg_fortran, errmsg, gravit
