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
  use kessler_perf_log, only: log_call, flush_log
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! public methods
  public :: kessler_eamxx_bridge_init_c
  public :: kessler_eamxx_bridge_run_c
  public :: kessler_eamxx_bridge_finalize_c

  ! Public variables
  integer, public            :: pcols
  integer, public            :: pver
  character(len=64),  public :: scheme_name = ""
  character(len=512), public :: errmsg = ""
  integer, public            :: errflg = 0
  logical, public            :: masterproc

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

subroutine kessler_eamxx_bridge_init_c( pcol_in, pver_in, lv_in, pref_in, rhoqr_in) bind(C, name="kessler_eamxx_bridge_init_c")
  ! Define uses here
  !-----------------------------------------------------------------------------
  ! Arguments
  integer(kind=c_int), value, intent(in) :: pcol_in
  integer(kind=c_int), value, intent(in) :: pver_in

  ! Things to pass along to the Kessler base code
  real(kind=c_real), value,    intent(in)  :: lv_in    ! latent heat of vaporization, J/kg
  real(kind=c_real), value,    intent(in)  :: pref_in  ! reference pressure, Pa
  real(kind=c_real), value,    intent(in)  :: rhoqr_in ! density of fresh liquid water, kg/m^3

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
  real(kind=c_real),     value,                intent(in)    :: dt       ! Physics time step (s)
  integer(kind=c_int), value,                intent(in)    :: lyr_surf ! Index of surface layer in the vertical coordinate
  integer(kind=c_int), value,                intent(in)    :: lyr_toa  ! Index of top of the atmosphere in the vertical coordinate
  real(kind=c_real),  dimension(pcols,pver), intent(in)    :: cpair    ! Specific_heat_of_dry_air_at_constant_pressure (J/kg/K)
  real(kind=c_real),    dimension(pcols,pver), intent(in)    :: rair     ! Gas constant of dry air (J/kg/K)
  real(kind=c_real),    dimension(pcols,pver), intent(in)    :: rho      ! Dry air density (kg/m^3)
  real(kind=c_real),    dimension(pcols,pver), intent(in)    :: z_mid    ! Heights of thermo. levels (m)
  real(kind=c_real),    dimension(pcols,pver), intent(in)    :: pk       ! Exner function (p/p0)**(R/cp)

  real(kind=c_real),    dimension(pcols,pver), intent(inout) :: theta    ! Potential temperature (K)
  real(kind=c_real),    dimension(pcols,pver), intent(inout) :: qv       ! Water vapor mixing ratio wrt dry air (kg/kg)
  real(kind=c_real),    dimension(pcols,pver), intent(inout) :: qc       ! Cloud water mixing ratio wrt dry air (kg/kg)
  real(kind=c_real),    dimension(pcols,pver), intent(inout) :: qr       ! Rain water mixing ratio wrt dry air (kg/kg)

  real(kind=c_real),    dimension(pcols),      intent(out)   :: precl    ! Precipitation rate (m_water / s)
  real(kind=c_real),    dimension(pcols,pver), intent(out)   :: relhum   ! Relative humidity in percent

  integer :: i,k
  ! real(kind=c_real) :: relhum_max, pk_max, theta_max, qv_max

  integer(kind=8) :: count_start, count_end, count_rate
  real(kind_phys) :: elapsed

  ! Call the Kessler run function
  call system_clock(count_start, count_rate)
  call kessler_run(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z_mid, &
        pk, theta, qv, qc, qr, precl, relhum, scheme_name, errmsg, errflg)
  call system_clock(count_end)
  elapsed = real(count_end - count_start, kind_phys) / real(count_rate, kind_phys)
  call log_call('kessler_run', ncol, nz, dt, elapsed)

end subroutine kessler_eamxx_bridge_run_c

!===================================================================================================

subroutine kessler_eamxx_bridge_finalize_c() bind(C, name="kessler_eamxx_bridge_finalize_c")
  ! Flush this rank's in-memory kessler_perf_log totals to CSV. Meant to be
  ! called exactly once per rank, at simulation finalize -- mirrors the
  ! JAX-side kessler.py's finalize() -> _flush_perf_log().
  call flush_log()

end subroutine kessler_eamxx_bridge_finalize_c

!===================================================================================================

end module kessler_eamxx_bridge_main