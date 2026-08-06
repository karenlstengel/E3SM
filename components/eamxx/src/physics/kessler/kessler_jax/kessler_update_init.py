# ---------------------------------------------------------------------------
# Translated by: Claude Fable 5 (claude-fable-5)
# Pass: final
# ---------------------------------------------------------------------------

"""
Plain Python translation of the scalar-only kessler_update initialization
procedure. Sets the gravitational-acceleration module constant.

Fortran source: kessler_update.F90, subroutine kessler_update_init.
SCALAR-ONLY: no arrays, no JAX (jax_required=False — validated on CPU).
"""


def kessler_update_init(gravit_in, errmsg, errflg, gravit):
    """
    Initialize the kessler_update module constant.

    MODULE VARIABLES (from kessler_update):
        gravit: MODULE variable (INOUT)

    Fortran:
        errmsg = ''
        errflg = 0
        gravit = gravit_in
    """
    errmsg = ""                    # [PY]
    errflg = 0                     # [PY]
    gravit = float(gravit_in)      # [PY]  module var set from intent(in) arg

    return errmsg, errflg, gravit
