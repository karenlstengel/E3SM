# LLM: Claude Fable 5 (claude-fable-5) — generated translation (copied from translations/claude/jax/, header added on promotion to _officialJAX)
# ---------------------------------------------------------------------------
# Translated by: Claude Fable 5 (claude-fable-5)
# Pass: final
# ---------------------------------------------------------------------------
"""SCALAR-ONLY translation of the Fortran subroutine `kessler_init` (module
`kessler`).

This procedure has no array arguments and is not called from any JAX/JIT
context, so it is plain Python: no JAX, no `_core` function (per the
scalar-only translation policy).

MODULE VARIABLES (from kessler), handled via the INOUT pattern — passed in
and returned, no hidden state:
    lv    : latent heat of vaporization [J/kg]        (INOUT)
    pref  : reference pressure, stored in hPa          (INOUT)
    rhoqr : density of liquid water [kg/m^3]           (INOUT)
"""


def kessler_init(lv_in, pref_in, rhoqr_in, errmsg, errflg, lv, pref, rhoqr):
    """Initialize the kessler module constants.

    Mirrors the Fortran:
        errmsg = ''
        errflg = 0
        lv     = lv_in
        pref   = pref_in / 100.  (Pa -> hPa)
        rhoqr  = rhoqr_in

    Returns: errmsg, errflg, lv, pref, rhoqr
    """
    errmsg = ""
    errflg = 0
    lv = float(lv_in)
    pref = float(pref_in) / 100.0
    rhoqr = float(rhoqr_in)
    return errmsg, errflg, lv, pref, rhoqr
