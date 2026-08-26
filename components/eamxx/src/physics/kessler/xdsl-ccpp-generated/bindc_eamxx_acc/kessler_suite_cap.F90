module kessler_suite_cap
  
  use ccpp_kinds
  use kessler, only: kessler_init
  use kessler, only: kessler_run
  use kessler_update, only: kessler_update_init
  use kessler_update, only: kessler_update_run
  use kessler_update, only: kessler_update_timestep_final
  use kessler_update, only: kessler_update_timestep_init
  
  implicit none
  private

  character(len=16) :: ccpp_suite_state = 'uninitialized'
  character(len=16), parameter :: const_in_time_step = 'in_time_step'
  character(len=16), parameter :: const_initialized = 'initialized'
  character(len=16), parameter :: const_uninitialized = 'uninitialized'
  public :: kessler_suite_register
  public :: kessler_suite_initialize
  public :: kessler_suite_finalize
  public :: kessler_suite_init_physics
  public :: kessler_suite_physics
  public :: kessler_suite_timestep_init_physics
  public :: kessler_suite_timestep_final_physics
  public :: kessler_suite_final_physics

CONTAINS
  
  subroutine kessler_suite_register(errflg, errmsg) 
    integer, intent(out) :: errflg
    character(len=512), intent(out) :: errmsg
    
    errflg = 0    
    errmsg = ''
  end subroutine kessler_suite_register 
  
  subroutine kessler_suite_initialize(errflg, errmsg) 
    integer, intent(out) :: errflg
    character(len=512), intent(out) :: errmsg
    
    errflg = 0    
    errmsg = ''
    if (.NOT. (const_uninitialized .eq. ccpp_suite_state)) then
      write(errmsg, '(3a)') "Invalid initial CCPP state, '", trim(ccpp_suite_state),              &
        "' in kessler_suite_initialize"
      errflg = 1      
    end if
    ccpp_suite_state = const_initialized
  end subroutine kessler_suite_initialize 
  
  subroutine kessler_suite_finalize(errflg, errmsg) 
    integer, intent(out) :: errflg
    character(len=512), intent(out) :: errmsg
    
    errflg = 0    
    errmsg = ''
    if (.NOT. (const_initialized .eq. ccpp_suite_state)) then
      write(errmsg, '(3a)') "Invalid initial CCPP state, '", trim(ccpp_suite_state),              &
        "' in kessler_suite_finalize"
      errflg = 1      
    end if
    ccpp_suite_state = const_uninitialized
  end subroutine kessler_suite_finalize 
  
  subroutine kessler_suite_init_physics(lv_in, pref_in, rhoqr_in, gravit_in, errmsg, errflg) 
    real(kind=kind_phys), intent(in) :: lv_in
    real(kind=kind_phys), intent(in) :: pref_in
    real(kind=kind_phys), intent(in) :: rhoqr_in
    real(kind=kind_phys), intent(in) :: gravit_in
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg
    
    errflg = 0    
    errmsg = ''
    if (.NOT. (const_initialized .eq. ccpp_suite_state)) then
      write(errmsg, '(3a)') "Invalid initial CCPP state, '", trim(ccpp_suite_state),              &
        "' in kessler_suite_init_physics"
      errflg = 1      
    end if
    if (errflg .eq. 0) then
      call kessler_init(lv_in=lv_in, pref_in=pref_in, rhoqr_in=rhoqr_in, errmsg=errmsg,           &
        errflg=errflg)
    end if
    if (errflg .eq. 0) then
      call kessler_update_init(gravit_in=gravit_in, errmsg=errmsg, errflg=errflg)
    end if
  end subroutine kessler_suite_init_physics 
  
  subroutine kessler_suite_physics(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, pk,      &
    theta, qv, qc, qr, precl, relhum, temp_prev, ttend_t, scheme_name, errmsg, errflg) 
    integer, intent(in) :: ncol
    integer, intent(in) :: nz
    real(kind=kind_phys), intent(in) :: dt
    integer, intent(in) :: lyr_surf
    integer, intent(in) :: lyr_toa
    real(kind=kind_phys), target, intent(in) :: cpair(:, :)
    real(kind=kind_phys), target, intent(in) :: rair(:, :)
    real(kind=kind_phys), target, intent(in) :: rho(:, :)
    real(kind=kind_phys), target, intent(in) :: z(:, :)
    real(kind=kind_phys), target, intent(in) :: pk(:, :)
    real(kind=kind_phys), target, intent(inout) :: theta(:, :)
    real(kind=kind_phys), target, intent(inout) :: qv(:, :)
    real(kind=kind_phys), target, intent(inout) :: qc(:, :)
    real(kind=kind_phys), target, intent(inout) :: qr(:, :)
    real(kind=kind_phys), target, intent(inout) :: precl(:)
    real(kind=kind_phys), target, intent(inout) :: relhum(:, :)
    real(kind=kind_phys), target, intent(in) :: temp_prev(:, :)
    real(kind=kind_phys), target, intent(inout) :: ttend_t(:, :)
    character(len=64), intent(out) :: scheme_name
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg
    
    errflg = 0    
    errmsg = ''
    if (.NOT. (const_in_time_step .eq. ccpp_suite_state)) then
      write(errmsg, '(3a)') "Invalid initial CCPP state, '", trim(ccpp_suite_state),              &
        "' in kessler_suite_physics"
      errflg = 1      
    end if
    if (errflg .eq. 0) then
      call kessler_run(ncol=ncol, nz=nz, dt=dt, lyr_surf=lyr_surf, lyr_toa=lyr_toa, cpair=cpair,  &
        rair=rair, rho=rho, z=z, pk=pk, theta=theta, qv=qv, qc=qc, qr=qr, precl=precl,            &
        relhum=relhum, scheme_name=scheme_name, errmsg=errmsg, errflg=errflg)
    end if
    if (errflg .eq. 0) then
      call kessler_update_run(nz=nz, ncol=ncol, dt=dt, theta=theta, exner=pk,                     &
        temp_prev=temp_prev, ttend_t=ttend_t, errmsg=errmsg, errflg=errflg)
    end if
  end subroutine kessler_suite_physics 
  
  subroutine kessler_suite_timestep_init_physics(ncol, nz, temp, temp_prev, ttend_t, errmsg,      &
    errflg) 
    integer, intent(in) :: ncol
    integer, intent(in) :: nz
    real(kind=kind_phys), target, intent(in) :: temp(:, :)
    real(kind=kind_phys), target, intent(inout) :: temp_prev(:, :)
    real(kind=kind_phys), target, intent(inout) :: ttend_t(:, :)
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg
    
    errflg = 0    
    errmsg = ''
    if (errflg .eq. 0) then
      call kessler_update_timestep_init(ncol=ncol, nz=nz, temp=temp, temp_prev=temp_prev,         &
        ttend_t=ttend_t, errmsg=errmsg, errflg=errflg)
    end if
    ccpp_suite_state = const_in_time_step
  end subroutine kessler_suite_timestep_init_physics 
  
  subroutine kessler_suite_timestep_final_physics(nz, ncol, cpair, temp, zm, phis, st_energy,     &
    errmsg, errflg) 
    integer, intent(in) :: nz
    integer, intent(in) :: ncol
    real(kind=kind_phys), target, intent(in) :: cpair(:, :)
    real(kind=kind_phys), target, intent(in) :: temp(:, :)
    real(kind=kind_phys), target, intent(in) :: zm(:, :)
    real(kind=kind_phys), target, intent(in) :: phis(:)
    real(kind=kind_phys), target, intent(inout) :: st_energy(:, :)
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg
    
    errflg = 0    
    errmsg = ''
    if (errflg .eq. 0) then
      call kessler_update_timestep_final(nz=nz, ncol=ncol, cpair=cpair, temp=temp, zm=zm,         &
        phis=phis, st_energy=st_energy, errmsg=errmsg, errflg=errflg)
    end if
    ccpp_suite_state = const_initialized
  end subroutine kessler_suite_timestep_final_physics 
  
  subroutine kessler_suite_final_physics(errflg, errmsg) 
    integer, intent(out) :: errflg
    character(len=512), intent(out) :: errmsg
    
    errflg = 0    
    errmsg = ''
    if (.NOT. (const_initialized .eq. ccpp_suite_state)) then
      write(errmsg, '(3a)') "Invalid initial CCPP state, '", trim(ccpp_suite_state),              &
        "' in kessler_suite_final_physics"
      errflg = 1      
    end if
  end subroutine kessler_suite_final_physics 
end module kessler_suite_cap
