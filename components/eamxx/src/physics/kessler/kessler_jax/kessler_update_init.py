# ---------------------------------------------------------------------------
# Translated by: Claude Fable 5 (claude-fable-5)
# Pass: final
# ---------------------------------------------------------------------------
"""SCALAR-ONLY translation of the Fortran subroutine `kessler_update_init`.

Plain Python, no JAX: the procedure has no array arguments and is not called
from any jitted context. The MODULE variable `gravit` (from `kessler_update`)
follows the INOUT pattern — passed in, set from `gravit_in`, returned.
"""


def kessler_update_init(gravit_in, errmsg, errflg, gravit):
    """
    Store the gravitational acceleration in the module variable.

    MODULE VARIABLES (from kessler_update):
        gravit: MODULE variable (INOUT)

    Returns: errmsg, errflg, gravit
    """
    errmsg = ""
    errflg = 0
    gravit = float(gravit_in)
    return errmsg, errflg, gravit
