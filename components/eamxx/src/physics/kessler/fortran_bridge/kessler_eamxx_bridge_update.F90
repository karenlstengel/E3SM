module kessler_eamxx_bridge_update

  use iso_c_binding
  use openacc_utils
  use cam_logfile,   only: iulog ! kinds instead of cam_logfile?
  use shr_sys_mod,   only: shr_sys_flush
  ! use spmd_utils,      only: masterproc

  ! Kessler code from CAM-SIMA
  use kessler_update
  use ccpp_kinds, only:  kind_phys
  use kessler_perf_log, only: log_call
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! public methods
  public :: kessler_eamxx_bridge_update_init_c
  public :: kessler_eamxx_bridge_update_c

  ! Public variables
  integer, public            :: pcols
  integer, public            :: pver
  character(len=64),  public :: scheme_name = ""
  character(len=512), public :: errmsg = ""
  integer, public            :: errflg = 0

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

subroutine kessler_eamxx_bridge_update_init_c(pcol_in, pver_in, gravit_in, cpair_in) bind(C, name="kessler_eamxx_bridge_update_init_c")
  ! Define uses here
  !-----------------------------------------------------------------------------
  ! Arguments
  integer(kind=c_int), value, intent(in) :: pcol_in
  integer(kind=c_int), value, intent(in) :: pver_in

  ! Things to pass along to the Kessler base code
  real(kind=c_real), value,    intent(in)  :: gravit_in    ! gravity acceleration m/s^2
  real(kind=c_real), value,    intent(in)  :: cpair_in     ! specific heat of dry air at constant pressure, J/kg/K

  ! Set dimensions of fields
  pcols = pcol_in
  pver  = pver_in

  errmsg = "temp"
  errflg = 0
  scheme_name = "KESSLER"

  ! Call the Kessler update init function
  call kessler_update_init(gravit_in, cpair_in, errmsg, errflg)

end subroutine kessler_eamxx_bridge_update_init_c

!===================================================================================================

subroutine kessler_eamxx_bridge_update_c( ncol, nz, dt, pk, theta, temp_prev, temp, temp_tend, z_mid, phis, st_energy) bind(C, name="kessler_eamxx_bridge_update_c")
  ! Define uses here
  !-----------------------------------------------------------------------------
  ! Arguments
  integer(kind=c_int),  value,                intent(in)    :: ncol      ! Number of columns
  integer(kind=c_int),  value,                intent(in)    :: nz        ! Number of vertical levels
  real(kind=c_real),    value,                intent(in)    :: dt        ! Physics time step (s)

  real(kind=c_real),    dimension(pcols,pver), intent(in)    :: pk        ! Exner function (p/p0)**(R/cp)
  real(kind=c_real),    dimension(pcols,pver), intent(in)    :: theta     ! Potential temperature (K)

  real(kind=c_real),    dimension(pcols,pver), intent(inout) :: temp_prev ! Temperature at previous timestep (K)
  real(kind=c_real),    dimension(pcols,pver), intent(inout) :: temp      ! Temperature updated (K)
  real(kind=c_real),    dimension(pcols,pver), intent(out)   :: temp_tend ! Temperature tendency (K)
  real(kind=c_real),    dimension(pcols,pver), intent(in)    :: z_mid     ! Geopotential height at each level (m)
  real(kind=c_real),    dimension(pcols),      intent(in)    :: phis      ! Geopotential height of surface (m2/s2)
  real(kind=c_real),    dimension(pcols,pver), intent(out)   :: st_energy ! Dry static energy J/kg

  integer(kind=8) :: count_start, count_end, count_rate
  real(kind_phys) :: elapsed

  ! Call the Kessler update Functions
  #if defined(EAMXX_ENABLE_GPU) && defined(EAMXX_ENABLE_OPENACC)
    call system_clock(count_start, count_rate)
    call kessler_update_timestep_init(ncol, nz, temp, temp_prev, temp_tend, errmsg, errflg)
    call system_clock(count_end)
    elapsed = real(count_end - count_start, kind_phys) / real(count_rate, kind_phys)
    call log_call('kessler_update_timestep_init', ncol, nz, dt, elapsed)

    call system_clock(count_start, count_rate)
    call kessler_update_run(nz, ncol, dt, theta, pk, temp_prev, temp_tend, errmsg, errflg)
    call system_clock(count_end)
    elapsed = real(count_end - count_start, kind_phys) / real(count_rate, kind_phys)
    call log_call('kessler_update_run', ncol, nz, dt, elapsed)

    call system_clock(count_start, count_rate)
    call kessler_update_timestep_final(nz, ncol, temp, z_mid, phis, st_energy, errflg, errmsg)
    call system_clock(count_end)
    elapsed = real(count_end - count_start, kind_phys) / real(count_rate, kind_phys)
    call log_call('kessler_update_timestep_final', ncol, nz, dt, elapsed)
  #else
    call system_clock(count_start, count_rate)
    call kessler_update_timestep_init(temp, temp_prev, temp_tend, errmsg, errflg)
    call system_clock(count_end)
    elapsed = real(count_end - count_start, kind_phys) / real(count_rate, kind_phys)
    call log_call('kessler_update_timestep_init', ncol, nz, dt, elapsed)

    call system_clock(count_start, count_rate)
    call kessler_update_run(nz, ncol, dt, theta, pk, temp_prev, temp_tend, errmsg, errflg)
    call system_clock(count_end)
    elapsed = real(count_end - count_start, kind_phys) / real(count_rate, kind_phys)
    call log_call('kessler_update_run', ncol, nz, dt, elapsed)

    call system_clock(count_start, count_rate)
    call kessler_update_timestep_final(nz, temp, z_mid, phis, st_energy, errflg, errmsg)
    call system_clock(count_end)
    elapsed = real(count_end - count_start, kind_phys) / real(count_rate, kind_phys)
    call log_call('kessler_update_timestep_final', ncol, nz, dt, elapsed)
  #endif
end subroutine kessler_eamxx_bridge_update_c

!===================================================================================================

end module kessler_eamxx_bridge_update