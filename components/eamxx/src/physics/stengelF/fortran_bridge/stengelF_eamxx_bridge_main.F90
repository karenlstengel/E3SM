module stengelF_eamxx_bridge_main

  use iso_c_binding
  use openacc_utils
  use mpi
  use cam_logfile,   only: iulog
  use shr_sys_mod,   only: shr_sys_flush
  use shr_file_mod
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! public methods
  public :: stengelF_eamxx_bridge_init_c
  public :: stengelF_eamxx_bridge_run_c
  public :: set_log_file_name_f90_c

  ! Public variables?
  integer, public            :: pcols
  integer, public            :: pver
  character(len=256), public :: log_fname = ""
  logical, public            :: masterproc

!===================================================================================================
#include "eamxx_config.f"
# define c_real c_double
!===================================================================================================
contains
!===================================================================================================

subroutine stengelF_eamxx_bridge_init_c( pcol_in, pver_in ) bind(C, name="stengelF_eamxx_bridge_init_c")
  ! Define uses here
  !-----------------------------------------------------------------------------
  ! Arguments
  integer(kind=c_int), value, intent(in) :: pcol_in
  integer(kind=c_int), value, intent(in) :: pver_in

  ! Local variables
  integer :: mpi_rank, ierror

  ! Set dimensions of fields
  pcols = pcol_in
  pver  = pver_in

  call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, ierror)
  masterproc = .false.
  if (mpi_rank==0) masterproc = .true.

  return
end subroutine stengelF_eamxx_bridge_init_c

!===================================================================================================

subroutine stengelF_eamxx_bridge_run_c( ncol, p_mid, T_mid ) bind(C, name="stengelF_eamxx_bridge_run_c")
  ! Define uses here
  !-----------------------------------------------------------------------------
  ! Arguments
  integer(kind=c_int),                value, intent(in   ) :: ncol
  real(kind=c_real),  dimension(pcols,pver), intent(inout) :: p_mid
  real(kind=c_real),  dimension(pcols,pver), intent(inout) :: T_mid

  ! Vars to store max values 
  real(kind=c_real) :: p_mid_max, T_mid_max

  ! Local variables for looping
  integer :: i,k

  ! initialize reduction scalars
  p_mid_max = 0.0
  T_mid_max = 0.0

  ! Scale and get max value for each field
  !$acc parallel deviceptr(p_mid, T_mid) reduction(max: p_mid_max, T_mid_max)
  !$acc loop gang vector collapse(2)
  do k = 1,pver
    do i = 1,ncol
      p_mid(i,k) = p_mid(i,k) * 0.5
      T_mid(i,k) = T_mid(i,k) * 0.5
      p_mid_max = max(p_mid_max,p_mid(i,k))
      t_mid_max = max(t_mid_max,t_mid(i,k))
    end do
  end do
  !$acc end parallel 

  ! TODO - need to fix this in terms of running in parallel. (if (masterproc) write...)
  if (masterproc) then
    write(*,*) "Fortran, max value of 0.5 p_mid ", p_mid_max
    write(*,*) "Fortran, max value of 0.5 T_mid ", T_mid_max
  end if

  ! Scale and get max value for each field
  !$acc parallel deviceptr(p_mid, T_mid) reduction(max: p_mid_max, T_mid_max)
  !$acc loop gang vector collapse(2)
  do k = 1,pver
    do i = 1,ncol
      p_mid(i,k) = p_mid(i,k) * 2.0
      T_mid(i,k) = T_mid(i,k) * 2.0
      p_mid_max = max(p_mid_max,p_mid(i,k))
      t_mid_max = max(t_mid_max,t_mid(i,k))
    end do
  end do
  !$acc end parallel

  ! TODO - need to fix this in terms of running in parallel. also newline?
  if (masterproc) then
    write(*,*) "Fortran, max value of 2 p_mid ", p_mid_max
    write(*,*) "Fortran, max value of 2 T_mid ", T_mid_max
  end if

  return
end subroutine stengelF_eamxx_bridge_run_c

!===================================================================================================

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

    iulog = shr_file_getunit()
    if (masterproc) then
      ! Set the log file on root rank...
      open (unit=iulog,file=trim(log_fname), &
            action='WRITE', access='SEQUENTIAL', position="append")
      write(iulog,*) " ---- STENGELF ----"
      flush(iulog)
    endif
  endif
end subroutine set_log_file_name_f90_c
!===================================================================================================

end module stengelF_eamxx_bridge_main