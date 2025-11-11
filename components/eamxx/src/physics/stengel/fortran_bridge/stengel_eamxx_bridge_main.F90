module stengel_eamxx_bridge_main

  use iso_c_binding
  use cam_logfile,   only: iulog
  use shr_sys_mod,   only: shr_sys_flush
  use stengel_eamxx_bridge_params, only: masterproc, r8, pcols, pver, pverp, top_lev ! TODO - required_arguments = p_mid, T_mid
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! public methods
  public :: stengel_eamxx_bridge_init_c
  public :: stengel_eamxx_bridge_run_c

!===================================================================================================
#include "eamxx_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif
!===================================================================================================
contains
!===================================================================================================

subroutine stengel_eamxx_bridge_init_c( pcol_in, pver_in ) bind(C)
  ! Do stuff here
end subroutine stengel_eamxx_bridge_init_c

subroutine stengel_eamxx_bridge_run_c( pcol_in, pver_in ) bind(C)
  ! Do stuff here
end subroutine stengel_eamxx_bridge_run_c

end module stengel_eamxx_bridge_main