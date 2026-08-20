# LLM: Gemini 3.1 PRO via Gemini CLI code agent — generated translation (copied from _2sd_exp_JDdata/gemini/, header added on promotion to _officialJAX)
"""
Bridge: kessler_init
Strategy: ALWAYS transpose + contiguous; JIT handled by _core decorator

MODULE: kessler
MODULE vars: lv, pref, rhoqr
Procedure type: INIT

MODULE_VARIABLE_POLICY (INOUT Pattern):
- MODULE vars passed as function parameters (INOUT)
- Driver manages MODULE state explicitly
- No hidden state dictionary
- INIT: Returns updated MODULE vars to driver
- COMPUTE: Returns updated MODULE vars to driver
"""

# No array conversions needed (all parameters are scalars)

from kessler_jax.kessler_init import kessler_init

def kessler_init_bridge(lv_in, pref_in, rhoqr_in, errmsg, errflg, lv, pref, rhoqr):
    """
    Production bridge for kessler_init.

    Strategy: transpose inputs/outputs (CPU), call wrapper directly.
    kessler_init_core is @jax.jit decorated — JIT and GPU execution
    are handled there, not here.

    MODULE VARIABLES (from kessler):
      - lv (INOUT)
      - pref (INOUT)
      - rhoqr (INOUT)
    
    Procedure type: INIT
    Pattern: MODULE vars passed as INOUT parameters
    """
    # Convert inputs

    # Compute — {proc_name}_core is @jax.jit decorated; JIT fires on first call
    errmsg_out, errflg_out, lv, pref, rhoqr = kessler_init(lv_in, pref_in, rhoqr_in, errmsg, errflg, lv, pref, rhoqr)

    # Convert outputs and return MODULE vars (INOUT pattern)

    return errmsg_out, errflg_out, lv, pref, rhoqr
