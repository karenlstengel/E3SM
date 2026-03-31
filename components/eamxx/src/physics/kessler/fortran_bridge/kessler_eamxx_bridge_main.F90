module kessler_eamxx_bridge_main

  use iso_c_binding
  use mpi
  use openacc_utils
  use cam_logfile,   only: iulog ! kinds instead of cam_logfile?
  use shr_sys_mod,   only: shr_sys_flush
  ! use spmd_utils,      only: masterproc

  ! Kessler code from CAM-SIMA
  use kessler
  use ccpp_kinds, only:  kind_phys
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! public methods
  public :: kessler_eamxx_bridge_init_c
  public :: kessler_eamxx_bridge_run_c

  ! Public variables
  integer, public            :: pcols
  integer, public            :: pver
  character(len=64),  public :: scheme_name = ""
  character(len=512), public :: errmsg = ""
  integer, public            :: errflg = 0
  logical, public            :: masterproc

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

  integer :: mpi_rank, ierror

  ! Set dimensions of fields
  pcols = pcol_in
  pver  = pver_in

  errmsg = "temp"
  errflg = 0
  scheme_name = "KESSLER"

  ! Call the Kessler init function 
  call kessler_init(lv_in, pref_in, rhoqr_in, errmsg, errflg)

  call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, ierror)
  masterproc = .false.
  if (mpi_rank==0) masterproc = .true.

end subroutine kessler_eamxx_bridge_init_c

!===================================================================================================

subroutine kessler_eamxx_bridge_run_c( ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z_mid, &
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
  real(kind_phys),    dimension(pcols,pver), intent(in)    :: z_mid    ! Heights of thermo. levels (m)
  real(kind_phys),    dimension(pcols,pver), intent(in)    :: pk       ! Exner function (p/p0)**(R/cp)

  real(kind_phys),    dimension(pcols,pver), intent(inout) :: theta    ! Potential temperature (K)
  real(kind_phys),    dimension(pcols,pver), intent(inout) :: qv       ! Water vapor mixing ratio wrt dry air (kg/kg)
  real(kind_phys),    dimension(pcols,pver), intent(inout) :: qc       ! Cloud water mixing ratio wrt dry air (kg/kg)
  real(kind_phys),    dimension(pcols,pver), intent(inout) :: qr       ! Rain water mixing ratio wrt dry air (kg/kg)

  real(kind_phys),    dimension(pcols),      intent(out)   :: precl    ! Precipitation rate (m_water / s)
  real(kind_phys),    dimension(pcols,pver), intent(out)   :: relhum   ! Relative humidity in percent

  integer :: i,k
  ! real(kind=c_real) :: relhum_max, pk_max, theta_max, qv_max

  ! real(kind=c_real) :: qv_sum, qc_sum, qr_sum
  if (masterproc) then
    write(*,*) "ncol: ", ncol
    write(*,*) "nz: ", nz
    write(*,*) "cpair: ", cpair(1,1)
    write(*,*) "z_mid: ", z_mid(1,1)
    write(*,*) "precl: ", precl(1)
  end if

  ! Call the Kessler run function
  call kessler_run(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z_mid, &
        pk, theta, qv, qc, qr, precl, relhum, scheme_name, errmsg, errflg)
  
  
  ! relhum_max = 0.0
  ! pk_max = 0.0
  ! theta_max = 0.0
  ! qv_max = 0.0

  ! do k = 1,pver
  !   do i = 1,ncol
  !     relhum_max = max(relhum_max,relhum(i,k))
  !     pk_max = max(pk_max,pk(i,k))
  !     theta_max = max(theta_max,theta(i,k))
  !     qv_max = max(qv_max,qv(i,k))
  !   end do
  ! end do

end subroutine kessler_eamxx_bridge_run_c

!===================================================================================================

end module kessler_eamxx_bridge_main