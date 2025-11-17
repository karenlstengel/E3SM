module stengel_eamxx_bridge_main

  use iso_c_binding
  use cam_logfile,   only: iulog
  use shr_sys_mod,   only: shr_sys_flush
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! public methods
  public :: stengel_eamxx_bridge_init_c
  public :: stengel_eamxx_bridge_run_c

  ! Public variables?
  integer, public :: pcols
  integer, public :: pver

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

subroutine stengel_eamxx_bridge_init_c( pcol, pver ) bind(C)
  ! Define uses here
  !-----------------------------------------------------------------------------
  ! Arguments
  integer(kind=c_int), value, intent(in) :: pcol
  integer(kind=c_int), value, intent(in) :: pver

  ! Set dimensions of fields
  pcols = pcol
  pver  = pver

  ! Do stuff here - initialize fields to 0?

  return
end subroutine stengel_eamxx_bridge_init_c

!===================================================================================================

subroutine stengel_eamxx_bridge_run_c( ncol, p_mid, T_mid ) bind(C)
  ! Define uses here
  !-----------------------------------------------------------------------------
  ! Arguments
  integer(kind=c_int),                value, intent(in   ) :: ncol
  real(kind=c_real),  dimension(pcols,pver), intent(inout) :: p_mid
  real(kind=c_real),  dimension(pcols,pver), intent(inout) :: T_mid

  ! Vars to store max values 
  real(kind=c_real) :: p_mid_max, T_mid_max

  ! Do stuff here - scale both p_mid and T_mid by 0.5, print, and then by 2.0, print
  p_mid = p_mid * 0.5 
  T_mid = T_mid * 0.5

  p_mid_max = MAXVAL(p_mid) 
  T_mid_max = MAXVAL(T_mid) 

  ! TODO - print but how? need to see how the logger works
  write(iulog,*) 'Fortran, max value of 0.5 p_mid ', p_mid_max
  write(iulog,*) 'Fortran, max value of 0.5 T_mid ', T_mid_max

  p_mid = p_mid * 2.0 
  T_mid = T_mid * 2.0

  p_mid_max = MAXVAL(p_mid) 
  T_mid_max = MAXVAL(T_mid)

  ! TODO - print but how? need to see how the logger works
  write(iulog,*) 'Fortran, max value of 2 p_mid ', p_mid_max
  write(iulog,*) 'Fortran, max value of 2 T_mid ', T_mid_max
  
  return
end subroutine stengel_eamxx_bridge_run_c

!===================================================================================================

end module stengel_eamxx_bridge_main