module kessler_eamxx_bridge_main

  use iso_c_binding
  ! use openacc_utils
  use cam_logfile,   only: iulog ! kinds instead of cam_logfile?
  use shr_sys_mod,   only: shr_sys_flush
  ! use spmd_utils,      only: masterproc

  ! Kessler code from CAM-SIMA
  use kessler
  use kessler_update 
  use ccpp_kinds, only:  kind_phys
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! public methods
  public :: kessler_eamxx_bridge_init_c
  public :: kessler_eamxx_bridge_run_c
  public :: set_log_file_name_f90_c ! Might remove this

  ! Public variables
  integer, public            :: pcols
  integer, public            :: pver
  character(len=256), public :: log_fname = ""
  character(len=64),  public :: scheme_name = ""
  character(len=512), public :: errmsg = ""
  integer, public            :: errflg = 0

!===================================================================================================
#include "eamxx_config.f"
# define c_real c_double
!===================================================================================================
contains
!===================================================================================================

subroutine kessler_eamxx_bridge_init_c( pcol_in, pver_in, lv_in, pref_in, rhoqr_in) bind(C, name="kessler_eamxx_bridge_init_c")
  ! Define uses here
  !-----------------------------------------------------------------------------
  ! Arguments
  integer(kind=c_int), value, intent(in) :: pcol_in
  integer(kind=c_int), value, intent(in) :: pver_in

  ! Things to pass along to the Kessler base code
  real(kind_phys), value,    intent(in)  :: lv_in    ! latent heat of vaporization, J/kg
  real(kind_phys), value,    intent(in)  :: pref_in  ! reference pressure, Pa
  real(kind_phys), value,    intent(in)  :: rhoqr_in ! density of fresh liquid water, kg/m^3

  ! Set dimensions of fields
  pcols = pcol_in
  pver  = pver_in

  errmsg = "temp"
  errflg = 0
  scheme_name = "KESSLER"

  ! Call the Kessler init function 
  call kessler_init(lv_in, pref_in, rhoqr_in, errmsg, errflg)

end subroutine kessler_eamxx_bridge_init_c

!===================================================================================================

subroutine kessler_eamxx_bridge_run_c( ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, &
        pk, theta, qv, qc, qr, precl, relhum) bind(C, name="kessler_eamxx_bridge_run_c")
  ! Define uses here
  !-----------------------------------------------------------------------------
  ! Arguments
  integer(kind=c_int), value,                intent(in)    :: ncol     ! Number of columns
  integer(kind=c_int), value,                intent(in)    :: nz       ! Number of vertical levels
  real(kind_phys),     value,                intent(in)    :: dt       ! Physics time step (s)
  integer(kind=c_int), value,                intent(in)    :: lyr_surf ! Index of surface layer in the vertical coordinate
  integer(kind=c_int), value,                intent(in)    :: lyr_toa  ! Index of top of the atmosphere in the vertical coordinate
  real(kind=c_real),  dimension(pcols,pver), intent(in)    :: cpair    ! Specific_heat_of_dry_air_at_constant_pressure (J/kg/K)
  real(kind_phys),    dimension(pcols,pver), intent(in)    :: rair     ! Gas constant of dry air (J/kg/K)
  real(kind_phys),    dimension(pcols,pver), intent(in)    :: rho      ! Dry air density (kg/m^3)
  real(kind_phys),    dimension(pcols,pver), intent(in)    :: z        ! Heights of thermo. levels (m)
  real(kind_phys),    dimension(pcols,pver), intent(in)    :: pk       ! Exner function (p/p0)**(R/cp)

  real(kind_phys),    dimension(pcols,pver), intent(inout) :: theta    ! Potential temperature (K)
  real(kind_phys),    dimension(pcols,pver), intent(inout) :: qv       ! Water vapor mixing ratio wrt dry air (kg/kg)
  real(kind_phys),    dimension(pcols,pver), intent(inout) :: qc       ! Cloud water mixing ratio wrt dry air (kg/kg)
  real(kind_phys),    dimension(pcols,pver), intent(inout) :: qr       ! Rain water mixing ratio wrt dry air (kg/kg)

  real(kind_phys),    dimension(pcols),      intent(out)   :: precl    ! Precipitation rate (m_water / s)
  real(kind_phys),    dimension(pcols,pver), intent(out)   :: relhum   ! Relative humidity in percent

  ! Call the Kessler run function
  call kessler_run(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, &
        pk, theta, qv, qc, qr, precl, relhum, scheme_name, errmsg, errflg)
  
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