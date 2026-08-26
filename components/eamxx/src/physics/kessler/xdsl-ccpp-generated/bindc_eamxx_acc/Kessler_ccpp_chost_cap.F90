module Kessler_ccpp_chost_cap

  use ccpp_kinds, only: kind_phys
  use iso_c_binding
  use kessler_suite_cap, only: kessler_suite_register
  use kessler_suite_cap, only: kessler_suite_initialize
  use kessler_suite_cap, only: kessler_suite_finalize
  use kessler_suite_cap, only: kessler_suite_physics
  use kessler_suite_cap, only: kessler_suite_timestep_init_physics
  use kessler_suite_cap, only: kessler_suite_timestep_final_physics
  use kessler_suite_cap, only: kessler_suite_init_physics
  use kessler_suite_cap, only: kessler_suite_final_physics

  implicit none
  private

  public :: Kessler_chost_physics_register
  public :: Kessler_chost_physics_initialize
  public :: Kessler_chost_physics_finalize
  public :: Kessler_chost_physics_run
  public :: Kessler_chost_physics_timestep_initial
  public :: Kessler_chost_physics_timestep_final
  public :: Kessler_chost_physics_physics_initial
  public :: Kessler_chost_physics_physics_final

contains

  subroutine Kessler_chost_physics_register(errmsg, errflg) &
      bind(C, name='Kessler_chost_physics_register')
    character(kind=c_char, len=1), intent(out) :: errmsg(*)
    integer(c_int),               intent(out) :: errflg
    integer :: i
    character(len=512) :: errmsg_f

    errmsg_f = ' '
    errflg = 0
    call kessler_suite_register(errflg, errmsg_f)
    do i = 1, len_trim(errmsg_f)
      errmsg(i) = errmsg_f(i:i)
    end do
    errmsg(len_trim(errmsg_f)+1) = c_null_char
  end subroutine Kessler_chost_physics_register

  subroutine Kessler_chost_physics_initialize(errmsg, errflg) &
      bind(C, name='Kessler_chost_physics_initialize')
    character(kind=c_char, len=1), intent(out) :: errmsg(*)
    integer(c_int),               intent(out) :: errflg
    integer :: i
    character(len=512) :: errmsg_f

    errmsg_f = ' '
    errflg = 0
    call kessler_suite_initialize(errflg, errmsg_f)
    do i = 1, len_trim(errmsg_f)
      errmsg(i) = errmsg_f(i:i)
    end do
    errmsg(len_trim(errmsg_f)+1) = c_null_char
  end subroutine Kessler_chost_physics_initialize

  subroutine Kessler_chost_physics_finalize(errmsg, errflg) &
      bind(C, name='Kessler_chost_physics_finalize')
    character(kind=c_char, len=1), intent(out) :: errmsg(*)
    integer(c_int),               intent(out) :: errflg
    integer :: i
    character(len=512) :: errmsg_f

    errmsg_f = ' '
    errflg = 0
    call kessler_suite_finalize(errflg, errmsg_f)
    do i = 1, len_trim(errmsg_f)
      errmsg(i) = errmsg_f(i:i)
    end do
    errmsg(len_trim(errmsg_f)+1) = c_null_char
  end subroutine Kessler_chost_physics_finalize

  subroutine Kessler_chost_physics_run( &
      ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z_mid, exner, theta, qv, qc,  &
      qr, precl, relhum, temp_prev, temp_tend, scheme_name, errmsg, errflg) &
      bind(C, name='Kessler_chost_physics_run')
    integer(c_int), value, intent(in) :: ncol
    integer(c_int), value, intent(in) :: nz
    real(c_double), value, intent(in) :: dt
    integer(c_int), value, intent(in) :: lyr_surf
    integer(c_int), value, intent(in) :: lyr_toa
    real(c_double), target, intent(in) :: cpair(ncol, nz)
    real(c_double), target, intent(in) :: rair(ncol, nz)
    real(c_double), target, intent(in) :: rho(ncol, nz)
    real(c_double), target, intent(in) :: z_mid(ncol, nz)
    real(c_double), target, intent(in) :: exner(ncol, nz)
    real(c_double), target, intent(inout) :: theta(ncol, nz)
    real(c_double), target, intent(inout) :: qv(ncol, nz)
    real(c_double), target, intent(inout) :: qc(ncol, nz)
    real(c_double), target, intent(inout) :: qr(ncol, nz)
    real(c_double), target, intent(inout) :: precl(ncol)
    real(c_double), target, intent(inout) :: relhum(ncol, nz)
    real(c_double), target, intent(in) :: temp_prev(ncol, nz)
    real(c_double), target, intent(inout) :: temp_tend(ncol, nz)
    character(kind=c_char, len=1), intent(out) :: scheme_name(*)
    character(kind=c_char, len=1), intent(out) :: errmsg(*)
    integer(c_int),               intent(out) :: errflg
    integer :: i
    character(len=64)  :: scheme_name_f
    character(len=512) :: errmsg_f

    errmsg_f = ' '
    scheme_name_f = ' '
    errflg = 0
    call kessler_suite_physics( &
        ncol, nz, real(dt, kind_phys), lyr_surf, lyr_toa, cpair, rair, rho,  &
        z_mid, exner, theta, qv, qc, qr, precl, relhum, temp_prev,  &
        temp_tend, scheme_name_f, errmsg_f, errflg)
    do i = 1, len_trim(scheme_name_f)
      scheme_name(i) = scheme_name_f(i:i)
    end do
    scheme_name(len_trim(scheme_name_f)+1) = c_null_char
    do i = 1, len_trim(errmsg_f)
      errmsg(i) = errmsg_f(i:i)
    end do
    errmsg(len_trim(errmsg_f)+1) = c_null_char
  end subroutine Kessler_chost_physics_run

  subroutine Kessler_chost_physics_timestep_initial( &
      ncol, nz, temp, temp_prev, temp_tend, errmsg, errflg) &
      bind(C, name='Kessler_chost_physics_timestep_initial')
    integer(c_int), value, intent(in) :: ncol
    integer(c_int), value, intent(in) :: nz
    real(c_double), target, intent(in) :: temp(ncol, nz)
    real(c_double), target, intent(inout) :: temp_prev(ncol, nz)
    real(c_double), target, intent(inout) :: temp_tend(ncol, nz)
    character(kind=c_char, len=1), intent(out) :: errmsg(*)
    integer(c_int),               intent(out) :: errflg
    integer :: i
    character(len=512) :: errmsg_f

    errmsg_f = ' '
    errflg = 0
    call kessler_suite_timestep_init_physics( &
        ncol, nz, temp, temp_prev, temp_tend, errmsg_f, errflg)
    do i = 1, len_trim(errmsg_f)
      errmsg(i) = errmsg_f(i:i)
    end do
    errmsg(len_trim(errmsg_f)+1) = c_null_char
  end subroutine Kessler_chost_physics_timestep_initial

  subroutine Kessler_chost_physics_timestep_final( &
      ncol, nz, cpair, temp, z_mid, phis, st_energy, errmsg, errflg) &
      bind(C, name='Kessler_chost_physics_timestep_final')
    integer(c_int), value, intent(in) :: ncol
    integer(c_int), value, intent(in) :: nz
    real(c_double), target, intent(in) :: cpair(ncol, nz)
    real(c_double), target, intent(in) :: temp(ncol, nz)
    real(c_double), target, intent(in) :: z_mid(ncol, nz)
    real(c_double), target, intent(in) :: phis(ncol)
    real(c_double), target, intent(inout) :: st_energy(ncol, nz)
    character(kind=c_char, len=1), intent(out) :: errmsg(*)
    integer(c_int),               intent(out) :: errflg
    integer :: i
    character(len=512) :: errmsg_f

    errmsg_f = ' '
    errflg = 0
    call kessler_suite_timestep_final_physics( &
        nz, ncol, cpair, temp, z_mid, phis, st_energy, errmsg_f, errflg)
    do i = 1, len_trim(errmsg_f)
      errmsg(i) = errmsg_f(i:i)
    end do
    errmsg(len_trim(errmsg_f)+1) = c_null_char
  end subroutine Kessler_chost_physics_timestep_final

  subroutine Kessler_chost_physics_physics_initial( &
      lv, pref, rhoqr, gravit, errmsg, errflg) &
      bind(C, name='Kessler_chost_physics_physics_initial')
    real(c_double), value, intent(in) :: lv
    real(c_double), value, intent(in) :: pref
    real(c_double), value, intent(in) :: rhoqr
    real(c_double), value, intent(in) :: gravit
    character(kind=c_char, len=1), intent(out) :: errmsg(*)
    integer(c_int),               intent(out) :: errflg
    integer :: i
    character(len=512) :: errmsg_f

    errmsg_f = ' '
    errflg = 0
    call kessler_suite_init_physics( &
        real(lv, kind_phys), real(pref, kind_phys), real(rhoqr, kind_phys),  &
        real(gravit, kind_phys), errmsg_f, errflg)
    do i = 1, len_trim(errmsg_f)
      errmsg(i) = errmsg_f(i:i)
    end do
    errmsg(len_trim(errmsg_f)+1) = c_null_char
  end subroutine Kessler_chost_physics_physics_initial

  subroutine Kessler_chost_physics_physics_final(errmsg, errflg) &
      bind(C, name='Kessler_chost_physics_physics_final')
    character(kind=c_char, len=1), intent(out) :: errmsg(*)
    integer(c_int),               intent(out) :: errflg
    integer :: i
    character(len=512) :: errmsg_f

    errmsg_f = ' '
    errflg = 0
    call kessler_suite_final_physics(errflg, errmsg_f)
    do i = 1, len_trim(errmsg_f)
      errmsg(i) = errmsg_f(i:i)
    end do
    errmsg(len_trim(errmsg_f)+1) = c_null_char
  end subroutine Kessler_chost_physics_physics_final

end module Kessler_ccpp_chost_cap
