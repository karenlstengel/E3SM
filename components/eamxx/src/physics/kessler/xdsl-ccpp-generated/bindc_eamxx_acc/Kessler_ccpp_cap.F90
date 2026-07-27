module Kessler_ccpp_cap
  
  use ccpp_kinds
  use iso_c_binding
  use eamxx_kessler_host_mod, only: cpair
  use eamxx_kessler_host_mod, only: dt
  use eamxx_kessler_host_mod, only: exner
  use eamxx_kessler_host_mod, only: gravit
  use eamxx_kessler_host_mod, only: lv
  use eamxx_kessler_host_mod, only: lyr_surf
  use eamxx_kessler_host_mod, only: lyr_toa
  use eamxx_kessler_host_mod, only: ncol
  use eamxx_kessler_host_mod, only: nz
  use eamxx_kessler_host_mod, only: phis
  use eamxx_kessler_host_mod, only: precl
  use eamxx_kessler_host_mod, only: pref
  use eamxx_kessler_host_mod, only: qc
  use eamxx_kessler_host_mod, only: qr
  use eamxx_kessler_host_mod, only: qv
  use eamxx_kessler_host_mod, only: rair
  use eamxx_kessler_host_mod, only: relhum
  use eamxx_kessler_host_mod, only: rho
  use eamxx_kessler_host_mod, only: rhoqr
  use eamxx_kessler_host_mod, only: scheme_name
  use eamxx_kessler_host_mod, only: st_energy
  use eamxx_kessler_host_mod, only: temp
  use eamxx_kessler_host_mod, only: temp_prev
  use eamxx_kessler_host_mod, only: temp_tend
  use eamxx_kessler_host_mod, only: theta
  use eamxx_kessler_host_mod, only: z_mid
  use kessler_suite_cap, only: kessler_suite_suite_finalize
  use kessler_suite_cap, only: kessler_suite_suite_initialize
  use kessler_suite_cap, only: kessler_suite_suite_physics
  use kessler_suite_cap, only: kessler_suite_suite_register
  use kessler_suite_cap, only: kessler_suite_suite_timestep_final
  use kessler_suite_cap, only: kessler_suite_suite_timestep_initial
  
  implicit none
  private

  character(len=13), parameter :: str_kessler_suite = 'kessler_suite'
  character(len=7), parameter :: str_physics = 'physics'
  public :: Kessler_ccpp_physics_register
  public :: Kessler_ccpp_physics_initialize
  public :: Kessler_ccpp_physics_finalize
  public :: Kessler_ccpp_physics_timestep_initial
  public :: Kessler_ccpp_physics_timestep_final
  public :: Kessler_ccpp_physics_run
  public :: ccpp_physics_suite_list
  public :: ccpp_physics_suite_part_list
  public :: ccpp_physics_suite_variables

CONTAINS
  
  subroutine Kessler_ccpp_physics_register(suite_name, errmsg, errflg) BIND(C,                    &
    name='Kessler_ccpp_physics_register') 
    character(kind=c_char, len=1), intent(in) :: suite_name(*)
    character(kind=c_char, len=1), intent(out) :: errmsg(*)
    integer(c_int), intent(out) :: errflg
    integer :: ccpp_c2f_i
    character(len=512) :: suite_name_f
    character(len=512) :: errmsg_f
    
    suite_name_f = ' '
    do ccpp_c2f_i = 1, len(suite_name_f) 
      if (suite_name(ccpp_c2f_i) == c_null_char) exit
      suite_name_f(ccpp_c2f_i:ccpp_c2f_i) = suite_name(ccpp_c2f_i)
    end do 
    errmsg_f = ' '
    errflg = 0    
    if (trim(suite_name_f) .eq. 'kessler_suite') then
      call kessler_suite_suite_register(errflg, errmsg_f)
    else
      write(errmsg_f, '(3a)') "No suite named ", trim(suite_name_f), " found"
      errflg = 1      
    end if
    do ccpp_c2f_i = 1, len_trim(errmsg_f) 
      errmsg(ccpp_c2f_i) = errmsg_f(ccpp_c2f_i:ccpp_c2f_i)
    end do 
    errmsg(len_trim(errmsg_f)+1) = c_null_char
  end subroutine Kessler_ccpp_physics_register 
  
  subroutine Kessler_ccpp_physics_initialize(suite_name, errmsg, errflg) BIND(C,                  &
    name='Kessler_ccpp_physics_initialize') 
    character(kind=c_char, len=1), intent(in) :: suite_name(*)
    character(kind=c_char, len=1), intent(out) :: errmsg(*)
    integer(c_int), intent(out) :: errflg
    integer :: ccpp_c2f_i
    character(len=512) :: suite_name_f
    character(len=512) :: errmsg_f
    
    suite_name_f = ' '
    do ccpp_c2f_i = 1, len(suite_name_f) 
      if (suite_name(ccpp_c2f_i) == c_null_char) exit
      suite_name_f(ccpp_c2f_i:ccpp_c2f_i) = suite_name(ccpp_c2f_i)
    end do 
    errmsg_f = ' '
    errflg = 0    
    if (trim(suite_name_f) .eq. 'kessler_suite') then
      call kessler_suite_suite_initialize(lv, pref, rhoqr, gravit, errmsg_f, errflg)
    else
      write(errmsg_f, '(3a)') "No suite named ", trim(suite_name_f), " found"
      errflg = 1      
    end if
    do ccpp_c2f_i = 1, len_trim(errmsg_f) 
      errmsg(ccpp_c2f_i) = errmsg_f(ccpp_c2f_i:ccpp_c2f_i)
    end do 
    errmsg(len_trim(errmsg_f)+1) = c_null_char
  end subroutine Kessler_ccpp_physics_initialize 
  
  subroutine Kessler_ccpp_physics_finalize(suite_name, errmsg, errflg) BIND(C,                    &
    name='Kessler_ccpp_physics_finalize') 
    character(kind=c_char, len=1), intent(in) :: suite_name(*)
    character(kind=c_char, len=1), intent(out) :: errmsg(*)
    integer(c_int), intent(out) :: errflg
    integer :: ccpp_c2f_i
    character(len=512) :: suite_name_f
    character(len=512) :: errmsg_f
    
    suite_name_f = ' '
    do ccpp_c2f_i = 1, len(suite_name_f) 
      if (suite_name(ccpp_c2f_i) == c_null_char) exit
      suite_name_f(ccpp_c2f_i:ccpp_c2f_i) = suite_name(ccpp_c2f_i)
    end do 
    errmsg_f = ' '
    errflg = 0    
    if (trim(suite_name_f) .eq. 'kessler_suite') then
      call kessler_suite_suite_finalize(errflg, errmsg_f)
    else
      write(errmsg_f, '(3a)') "No suite named ", trim(suite_name_f), " found"
      errflg = 1      
    end if
    do ccpp_c2f_i = 1, len_trim(errmsg_f) 
      errmsg(ccpp_c2f_i) = errmsg_f(ccpp_c2f_i:ccpp_c2f_i)
    end do 
    errmsg(len_trim(errmsg_f)+1) = c_null_char
  end subroutine Kessler_ccpp_physics_finalize 
  
  subroutine Kessler_ccpp_physics_timestep_initial(suite_name, errmsg, errflg) BIND(C,            &
    name='Kessler_ccpp_physics_timestep_initial') 
    character(kind=c_char, len=1), intent(in) :: suite_name(*)
    character(kind=c_char, len=1), intent(out) :: errmsg(*)
    integer(c_int), intent(out) :: errflg
    integer :: ccpp_c2f_i
    character(len=512) :: suite_name_f
    character(len=512) :: errmsg_f
    
    suite_name_f = ' '
    do ccpp_c2f_i = 1, len(suite_name_f) 
      if (suite_name(ccpp_c2f_i) == c_null_char) exit
      suite_name_f(ccpp_c2f_i:ccpp_c2f_i) = suite_name(ccpp_c2f_i)
    end do 
    errmsg_f = ' '
    errflg = 0    
    if (trim(suite_name_f) .eq. 'kessler_suite') then
#ifdef USE_GPU
      !$acc enter data copyin(temp, temp_prev, temp_tend)
#endif
      call kessler_suite_suite_timestep_initial(ncol, nz, temp, temp_prev, temp_tend, errmsg_f,   &
        errflg)
    else
      write(errmsg_f, '(3a)') "No suite named ", trim(suite_name_f), " found"
      errflg = 1      
    end if
    do ccpp_c2f_i = 1, len_trim(errmsg_f) 
      errmsg(ccpp_c2f_i) = errmsg_f(ccpp_c2f_i:ccpp_c2f_i)
    end do 
    errmsg(len_trim(errmsg_f)+1) = c_null_char
  end subroutine Kessler_ccpp_physics_timestep_initial 
  
  subroutine Kessler_ccpp_physics_timestep_final(suite_name, errmsg, errflg) BIND(C,              &
    name='Kessler_ccpp_physics_timestep_final') 
    character(kind=c_char, len=1), intent(in) :: suite_name(*)
    character(kind=c_char, len=1), intent(out) :: errmsg(*)
    integer(c_int), intent(out) :: errflg
    integer :: ccpp_c2f_i
    character(len=512) :: suite_name_f
    character(len=512) :: errmsg_f
    
    suite_name_f = ' '
    do ccpp_c2f_i = 1, len(suite_name_f) 
      if (suite_name(ccpp_c2f_i) == c_null_char) exit
      suite_name_f(ccpp_c2f_i:ccpp_c2f_i) = suite_name(ccpp_c2f_i)
    end do 
    errmsg_f = ' '
    errflg = 0    
    if (trim(suite_name_f) .eq. 'kessler_suite') then
#ifdef USE_GPU
      !$acc data copyin(phis) copyout(st_energy)
#endif
      call kessler_suite_suite_timestep_final(nz, ncol, cpair, temp, z_mid, phis, st_energy,      &
        errmsg_f, errflg)
#ifdef USE_GPU
      !$acc end data
#endif
#ifdef USE_GPU
      !$acc exit data delete(cpair, temp, z_mid)
#endif
    else
      write(errmsg_f, '(3a)') "No suite named ", trim(suite_name_f), " found"
      errflg = 1      
    end if
    do ccpp_c2f_i = 1, len_trim(errmsg_f) 
      errmsg(ccpp_c2f_i) = errmsg_f(ccpp_c2f_i:ccpp_c2f_i)
    end do 
    errmsg(len_trim(errmsg_f)+1) = c_null_char
  end subroutine Kessler_ccpp_physics_timestep_final 
  
  subroutine Kessler_ccpp_physics_run(suite_name, suite_part, col_start, col_end, errmsg,         &
    errflg) BIND(C, name='Kessler_ccpp_physics_run') 
    character(kind=c_char, len=1), intent(in) :: suite_name(*)
    character(kind=c_char, len=1), intent(in) :: suite_part(*)
    integer(c_int), value, intent(in) :: col_start
    integer(c_int), value, intent(in) :: col_end
    character(kind=c_char, len=1), intent(inout) :: errmsg(*)
    integer(c_int), intent(inout) :: errflg
    integer :: ccpp_c2f_i
    character(len=512) :: suite_name_f
    character(len=512) :: suite_part_f
    character(len=512) :: errmsg_f
    
    suite_name_f = ' '
    do ccpp_c2f_i = 1, len(suite_name_f) 
      if (suite_name(ccpp_c2f_i) == c_null_char) exit
      suite_name_f(ccpp_c2f_i:ccpp_c2f_i) = suite_name(ccpp_c2f_i)
    end do 
    suite_part_f = ' '
    do ccpp_c2f_i = 1, len(suite_part_f) 
      if (suite_part(ccpp_c2f_i) == c_null_char) exit
      suite_part_f(ccpp_c2f_i:ccpp_c2f_i) = suite_part(ccpp_c2f_i)
    end do 
    errmsg_f = ' '
    do ccpp_c2f_i = 1, len(errmsg_f) 
      if (errmsg(ccpp_c2f_i) == c_null_char) exit
      errmsg_f(ccpp_c2f_i:ccpp_c2f_i) = errmsg(ccpp_c2f_i)
    end do 
    errflg = 0    
    if (trim(suite_name_f) .eq. 'kessler_suite') then
      if (trim(suite_part_f) .eq. 'physics') then
#ifdef USE_GPU
        !$acc enter data copyin(cpair, z_mid)
#endif
#ifdef USE_GPU
        !$acc data copy(theta(col_start:col_end, 1:nz), qv(col_start:col_end, 1:nz), &
        !$acc      qc(col_start:col_end, 1:nz), qr(col_start:col_end, 1:nz)) &
        !$acc      copyin(rair(col_start:col_end, 1:nz), rho(col_start:col_end, 1:nz), &
        !$acc      exner(col_start:col_end, 1:nz)) copyout(precl, relhum(col_start:col_end, 1:nz))
#endif
        call kessler_suite_suite_physics(col_start, col_end, nz, dt, lyr_surf, lyr_toa,           &
          cpair(col_start:col_end, 1:nz), rair(col_start:col_end, 1:nz), rho(col_start:col_end,   &
          1:nz), z_mid(col_start:col_end, 1:nz), exner(col_start:col_end, 1:nz),                  &
          theta(col_start:col_end, 1:nz), qv(col_start:col_end, 1:nz), qc(col_start:col_end,      &
          1:nz), qr(col_start:col_end, 1:nz), precl, relhum(col_start:col_end, 1:nz),             &
          temp_prev(col_start:col_end, 1:nz), temp_tend(col_start:col_end, 1:nz), scheme_name,    &
          errmsg_f, errflg)
#ifdef USE_GPU
        !$acc end data
#endif
#ifdef USE_GPU
        !$acc exit data copyout(temp_prev, temp_tend)
#endif
      else
        write(errmsg_f, '(3a)') "No suite part named ", trim(suite_part_f),                       &
          " found in suite kessler_suite"
        errflg = 1        
      end if
    else
      write(errmsg_f, '(3a)') "No suite named ", trim(suite_name_f), " found"
      errflg = 1      
    end if
    do ccpp_c2f_i = 1, len_trim(errmsg_f) 
      errmsg(ccpp_c2f_i) = errmsg_f(ccpp_c2f_i:ccpp_c2f_i)
    end do 
    errmsg(len_trim(errmsg_f)+1) = c_null_char
  end subroutine Kessler_ccpp_physics_run 
  
  subroutine ccpp_physics_suite_list(suites) 
    character(len=*), allocatable, intent(out) :: suites(:)
    
    allocate(suites(1))
    suites(1) = str_kessler_suite
  end subroutine ccpp_physics_suite_list 
  
  subroutine ccpp_physics_suite_part_list(suite_name, part_list, errmsg, errflg) 
    character(len=*), intent(in) :: suite_name
    character(len=*), allocatable, intent(out) :: part_list(:)
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg
    
    errflg = 0    
    if (trim(suite_name) .eq. 'kessler_suite') then
      allocate(part_list(1))
      part_list(1) = str_physics
    else
      write(errmsg, '(3a)') "No suite named ", trim(suite_name), " found"
      errflg = 1      
    end if
  end subroutine ccpp_physics_suite_part_list 
  subroutine ccpp_physics_suite_variables(suite_name, var_list, errmsg, errflg, input_vars,       &
    output_vars)
    character(len=*), intent(in) :: suite_name
    character(len=*), allocatable, intent(out) :: var_list(:)
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg
    logical, optional, intent(in) :: input_vars
    logical, optional, intent(in) :: output_vars
    logical :: do_input, do_output
    errmsg = ''
    errflg = 0
    do_input = .true.
    do_output = .true.
    if (present(input_vars)) do_input = input_vars
    if (present(output_vars)) do_output = output_vars
    if (trim(suite_name) .eq. 'kessler_suite') then
      if (do_input .and. .not. do_output) then
        allocate(var_list(22))
        var_list(1) = 'air_potential_temperature           '
        var_list(2) = 'air_temperature                     '
        var_list(3) = 'air_temperature_on_previous_timestep'
        var_list(4) = 'cloud_liquid_water_mixing_ratio_wrt_dry_air'
        var_list(5) = 'composition_dependent_gas_constant_of_dry_air'
        var_list(6) = 'composition_dependent_specific_heat_of_dry_air_at_constant_pressure'
        var_list(7) = 'dimensionless_exner_function        '
        var_list(8) = 'dry_air_density                     '
        var_list(9) = 'fresh_liquid_water_density_at_0c    '
        var_list(10) = 'geopotential_height_wrt_surface     '
        var_list(11) = 'horizontal_dimension                '
        var_list(12) = 'latent_heat_of_vaporization_of_water_at_0c'
        var_list(13) = 'rain_mixing_ratio_wrt_dry_air       '
        var_list(14) = 'standard_gravitational_acceleration '
        var_list(15) = 'surface_geopotential                '
        var_list(16) = 'surface_reference_pressure          '
        var_list(17) = 'tendency_of_air_temperature_due_to_model_physics'
        var_list(18) = 'timestep_for_physics                '
        var_list(19) = 'vertical_index_at_surface_adjacent_layer'
        var_list(20) = 'vertical_index_at_top_adjacent_layer'
        var_list(21) = 'vertical_layer_dimension            '
        var_list(22) = 'water_vapor_mixing_ratio_wrt_dry_air'
      else if (.not. do_input .and. do_output) then
        allocate(var_list(12))
        var_list(1) = 'air_potential_temperature           '
        var_list(2) = 'air_temperature_on_previous_timestep'
        var_list(3) = 'ccpp_error_code                     '
        var_list(4) = 'ccpp_error_message                  '
        var_list(5) = 'cloud_liquid_water_mixing_ratio_wrt_dry_air'
        var_list(6) = 'dry_static_energy                   '
        var_list(7) = 'rain_mixing_ratio_wrt_dry_air       '
        var_list(8) = 'relative_humidity                   '
        var_list(9) = 'scheme_name                         '
        var_list(10) = 'tendency_of_air_temperature_due_to_model_physics'
        var_list(11) = 'total_precipitation_rate_at_surface '
        var_list(12) = 'water_vapor_mixing_ratio_wrt_dry_air'
      else
        allocate(var_list(28))
        var_list(1) = 'air_potential_temperature           '
        var_list(2) = 'air_temperature                     '
        var_list(3) = 'air_temperature_on_previous_timestep'
        var_list(4) = 'ccpp_error_code                     '
        var_list(5) = 'ccpp_error_message                  '
        var_list(6) = 'cloud_liquid_water_mixing_ratio_wrt_dry_air'
        var_list(7) = 'composition_dependent_gas_constant_of_dry_air'
        var_list(8) = 'composition_dependent_specific_heat_of_dry_air_at_constant_pressure'
        var_list(9) = 'dimensionless_exner_function        '
        var_list(10) = 'dry_air_density                     '
        var_list(11) = 'dry_static_energy                   '
        var_list(12) = 'fresh_liquid_water_density_at_0c    '
        var_list(13) = 'geopotential_height_wrt_surface     '
        var_list(14) = 'horizontal_dimension                '
        var_list(15) = 'latent_heat_of_vaporization_of_water_at_0c'
        var_list(16) = 'rain_mixing_ratio_wrt_dry_air       '
        var_list(17) = 'relative_humidity                   '
        var_list(18) = 'scheme_name                         '
        var_list(19) = 'standard_gravitational_acceleration '
        var_list(20) = 'surface_geopotential                '
        var_list(21) = 'surface_reference_pressure          '
        var_list(22) = 'tendency_of_air_temperature_due_to_model_physics'
        var_list(23) = 'timestep_for_physics                '
        var_list(24) = 'total_precipitation_rate_at_surface '
        var_list(25) = 'vertical_index_at_surface_adjacent_layer'
        var_list(26) = 'vertical_index_at_top_adjacent_layer'
        var_list(27) = 'vertical_layer_dimension            '
        var_list(28) = 'water_vapor_mixing_ratio_wrt_dry_air'
      end if
    else
      write(errmsg, '(3a)') "No suite named ", trim(suite_name), " found"
      errflg = 1
    end if
  end subroutine ccpp_physics_suite_variables
end module Kessler_ccpp_cap
