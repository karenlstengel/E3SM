"""
Bridge: kessler_update_init
Strategy: PATH C device-side layout conversion — pure H2D/D2H,
in-jit axis reversal fused by XLA; public contract unchanged

MODULE: kessler_update
MODULE vars: gravit
Procedure type: INIT

MODULE_VARIABLE_POLICY (INOUT Pattern):
- MODULE vars passed as function parameters (INOUT)
- Driver manages MODULE state explicitly
- No hidden state dictionary
- INIT: Returns updated MODULE vars to driver
- COMPUTE: Returns updated MODULE vars to driver
"""

# No array conversions needed (all parameters are scalars)

from kessler_jax.kessler_update_init import kessler_update_init

def kessler_update_init_bridge(gravit_in, errmsg, errflg, gravit):
    """
    Production bridge for kessler_update_init (scalar-only).

    No array arguments — no layout conversion exists; calls the
    translated wrapper directly (strings and Python control flow
    stay host-side).

    MODULE VARIABLES (from kessler_update):
      - gravit (INOUT)
    
    Procedure type: INIT
    Pattern: MODULE vars passed as INOUT parameters
    """
    errmsg_out, errflg_out, gravit = kessler_update_init(gravit_in, errmsg, errflg, gravit)

    return errmsg_out, errflg_out, gravit
