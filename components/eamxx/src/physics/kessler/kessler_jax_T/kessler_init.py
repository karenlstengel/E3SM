# ---------------------------------------------------------------------------
# Translated by: Claude Opus 4.8 (claude-opus-4-8[1m])
# Build: LITERAL transliteration baseline (_literal_baseline)
# Pass: final
# ---------------------------------------------------------------------------

"""
Plain Python translation of the scalar-only Kessler initialization procedure.
This procedure uses only scalar values and explicitly threaded module variables.

Scalar-only: no arrays → no layout and no loops, so this is identical to the
idiomatic build (nothing to transliterate).
"""


def kessler_init(lv_in, pref_in, rhoqr_in, errmsg, errflg, lv, pref, rhoqr):
    """
    Initialize scalar module parameters for the Kessler microphysics scheme.

    MODULE VARIABLES (from kessler):
        lv: MODULE variable (INOUT)
        pref: MODULE variable (INOUT)
        rhoqr: MODULE variable (INOUT)
    """
    errmsg = ""
    errflg = 0
    lv = float(lv_in)
    pref = float(pref_in) / 100.0
    rhoqr = float(rhoqr_in)

    return errmsg, errflg, lv, pref, rhoqr
