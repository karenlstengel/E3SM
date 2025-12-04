module kessler_eamxx_bridge_main

  use iso_c_binding
  use openacc_utils
  use cam_logfile,   only: iulog ! kinds instead of cam_logfile?
  use shr_sys_mod,   only: shr_sys_flush
  use spmd_utils,      only: masterproc
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! public methods
  public :: kessler_eamxx_bridge_init_c
  public :: kessler_eamxx_bridge_run_c
  public :: set_log_file_name_f90_c

  ! Public variables?
  integer, public            :: pcols
  integer, public            :: pver
  character(len=256), public :: log_fname = ""

!===================================================================================================
#include "eamxx_config.f"
# define c_real c_double
!===================================================================================================
contains
!===================================================================================================

subroutine kessler_eamxx_bridge_init_c( pcol_in, pver_in ) bind(C, name="kessler_eamxx_bridge_init_c")
  ! Define uses here
  !-----------------------------------------------------------------------------
  ! Arguments
  integer(kind=c_int), value, intent(in) :: pcol_in
  integer(kind=c_int), value, intent(in) :: pver_in

  ! Set dimensions of fields
  pcols = pcol_in
  pver  = pver_in

  ! TODO - 

  return
end subroutine kessler_eamxx_bridge_init_c

!===================================================================================================

subroutine kessler_eamxx_bridge_run_c( ncol, TODO) bind(C, name="kessler_eamxx_bridge_run_c")
  ! Define uses here
  !-----------------------------------------------------------------------------
  ! Arguments
  integer(kind=c_int),                value, intent(in   ) :: ncol
  ! TODO - add in all agruments and types.

  ! TODO - add in the calls to do the physics. 
  return
end subroutine kessler_eamxx_bridge_run_c

!===================================================================================================

! TODO - fix the IO for running in parallel (I.e. do the masterproc thing)
subroutine set_log_file_name_f90_c(c_str) bind(C, name="set_log_file_name_f90_c")
  type (c_ptr), intent(in) :: c_str
  !
  ! Local(s)
  !
  character(len=256), pointer :: full_name
  character(len=256) :: path, fname
  integer :: len, slash, ierr

  call c_f_pointer(c_str,full_name)
  len = index(full_name, C_NULL_CHAR) -1
  if (len>0) then
    ! Search last slash in the (trimmed) full name
    slash = index(full_name(1:len),'/',back=.true.)

    ! Note: if there's no slash (relative filename),
    ! then slash=0, and path is the empty string.
    ! Otherwise, path ends with the slash
    path = full_name(1:slash)
    fname = full_name(slash+1:len)

    log_fname = trim(path)//fname

    ! Create the log file on root rank...
    open (unit=iulog,file=trim(log_fname), &
          action='WRITE', access='SEQUENTIAL', position="append")
    
    write(iulog,*) " ---- KESSLER TEST ----"
    flush(iulog)

  endif
end subroutine set_log_file_name_f90_c
!===================================================================================================

end module kessler_eamxx_bridge_main