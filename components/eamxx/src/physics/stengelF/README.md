# Adding new package to EAMxx

### Things of note

1. For the Fortran link to work, you must add both the C++ side and Fortran side field variables to the `ATMBufferManager` with `Packagename::init_buffers(const ATMBufferManager &buffer_manager)` and `Packagename::requested_buffer_size_in_bytes()` (see the main interface C++ file).
2. It seems like a good idea to create structs for passing input fields, output fields, and (if needed) parameters/etc. Note that this requires each member of the struct that is a field to be defined for both the C++ side (usually with `Spack` and Kokkos views) and the Fortran side (_unmanaged_ view of type `Real` in the same dimensions as the C++ version). 
3. You should create a transpose function for converting between the Fortran and C++ fields. See `transpose()` in the `packagename_functions.hpp` file. 
4. I _think_ you can only pass to Fortran the `get_field_out("fieldname").get_view<Spack**>().data()` object into the C to Fortran binding. 

### The main C++ code 

Create a new folder in `components/eamxx/src/physics/` with the name of your package: `components/eamxx/src/physics/PACKAGENAME`.

In this folder, we must add the following files: 

1. `eamxx_packagename_process_interface.cpp`
  ```cpp
  #include "eamxx_packagename_process_interface.hpp"
  #include "packagename_eamxx_bridge.hpp"
  #include "share/property_checks/field_within_interval_check.hpp"
  #include "share/field/field_utils.hpp"

  #include <ekat_assert.hpp>
  #include <ekat_units.hpp>

  #include <array>

  #ifdef EAMXX_HAS_PYTHON
  #include "share/atm_process/atmosphere_process_pyhelpers.hpp"
  #endif

  namespace scream
  {
    using namespace packagename;
  // =========================================================================================
  //  Inputs (these are inherited from AtomoshpereProcess which means we can use the same logger):
  //      comm - an EKAT communication group
  //      params - a parameter list of options for the process.

  Packagename::Packagename (const ekat::Comm& comm, const ekat::ParameterList& params)
    : AtmosphereProcess(comm, params)
  {
    // Nothing to do here usually
    m_atm_logger->info("[EAMxx] Packagename processes constructor");

    // Set the log filename in the F90 interface
    const char* logname = m_atm_logger->get_logfile_name().c_str();
    set_log_file_name_f90(&logname);

  }

  // =========================================================================================
  void Packagename::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
  {
    // using namespace ekat::units;
    // using namespace ShortFieldTagsNames;

    m_atm_logger->info("[EAMxx] Packagename processes set grids");

    constexpr auto K = ekat::units::K;
    constexpr auto Pa = ekat::units::Pa;

    // specify which grid to use
    m_grid = grids_manager->get_grid("physics");
    const auto& grid_name = m_grid->name();
    m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
    m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

    FieldLayout scalar3d_layout_mid = m_grid->get_3d_scalar_layout(true);

    constexpr int ps = 1;
    add_field<Updated>("T_mid", scalar3d_layout_mid, K, grid_name, ps);
    add_field<Updated>("p_mid", scalar3d_layout_mid, Pa, grid_name, ps);
  }

  // =========================================================================================
  void Packagename::initialize_impl (const RunType /* run_type */)
  {
    m_atm_logger->info("[EAMxx] packagename processes initialize_impl: ");

    packagename::packagename_eamxx_bridge_init(m_num_cols, m_num_levs);
  }

  // =========================================================================================

  // run_impl is called every timestep and where all of the physics happens
  // Inputs:
  //    - dt - the timestep for the current run step

  void Packagename::run_impl (const double /* dt */)
  {
    // Pull in variables .
    auto T_mid   = get_field_out("T_mid");
    auto p_mid   = get_field_out("p_mid");

    m_atm_logger->info("[EAMxx] packagename run_impl: ");
    
    auto T_mid_max = field_max<Real>(T_mid);
    m_atm_logger->info("\t max value for updated T_mid: "+ std::to_string(T_mid_max));
    auto p_mid_max = field_max<Real>(p_mid);
    m_atm_logger->info("\t max value for updated p_mid: "+ std::to_string(p_mid_max));

    // Add the fields we pulled in to the params struct to pass to the bridge:
    params.p_mid = p_mid.get_view<Spack**, Host>();
    params.T_mid = T_mid.get_view<Spack**, Host>();

    // Initialize fortran data holders in struct
    params.init(m_num_cols, m_num_levs);

    packagename_eamxx_bridge_run(m_num_cols, m_num_levs, params); 

    // Update with the new (should be the same) values from the run
    // TODO - figure out how to read the data from the params struct back into the fields on the C++ side.
    // Look at ZM but also probably others

    // for (int i=0; i<m_num_cols; ++i) {
    //   for (int j=0; j<m_num_levs; ++j) {
    //     p_mid(i,j) = params.p_mid(i,j);
    //     T_mid(i,j) = params.T_mid(i,j);
    //   }
    // }

    // T_mid_max = field_max<Real>(T_mid);
    // m_atm_logger->info("\t max value for T_mid after the Fortran bridge: "+ std::to_string(T_mid_max));
    // p_mid_max = field_max<Real>(p_mid);
    // m_atm_logger->info("\t max value for p_mid after the Fortran bridge: "+ std::to_string(p_mid_max));
    
  }

  // =========================================================================================
  void Packagename::finalize_impl()
  {
    // Do nothing
    m_atm_logger->info("[EAMxx] Packagename processes clean up.");
  }
  // =========================================================================================

  size_t Packagename::requested_buffer_size_in_bytes() const
  {
    const int nlevm_packs = ekat::npack<Spack>(m_num_levs);
    size_t buffer_size = 0;

    buffer_size+= PackagenameFunc::params::num_2d_midlv_c_views * sizeof(Spack) * m_num_cols * nlevm_packs;
    buffer_size+= PackagenameFunc::params::num_2d_midlv_f_views * sizeof(Real)  * m_num_cols * m_num_levs;

    return buffer_size;
  }

  /*------------------------------------------------------------------------------------------------*/

  void Packagename::init_buffers(const ATMBufferManager &buffer_manager)
  {
    auto buffer_chk = ( buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes() );
    EKAT_REQUIRE_MSG(buffer_chk,"Error! Buffers size not sufficient.\n");

    const int nlevm_packs = ekat::npack<Spack>(m_num_levs);

    constexpr auto num_2d_midlv_c_views = PackagenameFunc::params::num_2d_midlv_c_views;
    constexpr auto num_2d_midlv_f_views = PackagenameFunc::params::num_2d_midlv_f_views;
    
    //----------------------------------------------------------------------------
    Real* r_mem = reinterpret_cast<Real*>(buffer_manager.get_memory());
    //----------------------------------------------------------------------------
    // 2D "f_" views on mid-point levels
    PackagenameFunc::uview_2dl<Real>* midlv_f_ptrs[num_2d_midlv_f_views]  = { &params.f_p_mid, &params.f_T_mid};
    for (int i=0; i<num_2d_midlv_f_views; ++i) {
      *midlv_f_ptrs[i] = PackagenameFunc::uview_2dl<Real>(r_mem, m_num_cols, m_num_levs);
      r_mem += midlv_f_ptrs[i]->size();
    }
    //----------------------------------------------------------------------------
    Spack* spk_mem = reinterpret_cast<Spack*>(r_mem);
    //----------------------------------------------------------------------------
    // 2D views on mid-point levels
    PackagenameFunc::view_2d<Spack>* midlv_c_ptrs[num_2d_midlv_c_views]  = { &params.p_mid, &params.T_mid};
    for (int i=0; i<num_2d_midlv_c_views; ++i) {
      *midlv_c_ptrs[i] = PackagenameFunc::view_2d<Spack>(spk_mem, m_num_cols, nlevm_packs);
      spk_mem += midlv_c_ptrs[i]->size();
    }
    //----------------------------------------------------------------------------
    Real* total_mem = reinterpret_cast<Real*>(spk_mem);
    size_t used_mem = (reinterpret_cast<Real*>(total_mem) - buffer_manager.get_memory())*sizeof(Real);
    auto mem_chk = ( used_mem == requested_buffer_size_in_bytes() );
    EKAT_REQUIRE_MSG(mem_chk,"Error! Used memory != requested memory for Packagename.");
  }
  } // namespace scream
  ```

2. `eamxx_packagename_process_interface.hpp`
  ```cpp
  #ifndef SCREAM_PACKAGENAME_HPP
  #define SCREAM_PACKAGENAME_HPP

  #include "physics/packagename/packagename_functions.hpp"
  #include "share/atm_process/atmosphere_process.hpp"
  #include "share/atm_process/ATMBufferManager.hpp"

  #include <ekat_parameter_list.hpp>

  #include <string>

  namespace scream
  {

  /*
  * The class responsible to do packagename physics
  *
  * The AD should store exactly ONE instance of this class stored
  * in its list of subcomponents (the AD should make sure of this).
  */

  class Packagename : public AtmosphereProcess
  {
  public:
    using PackagenameFunc = packagename::PackagenameFunctions<Real, DefaultDevice>;
    using Spack        = PackagenameFunc::Spack;
    using Pack         = ekat::Pack<Real,Spack::n>;

    // Constructors
    Packagename (const ekat::Comm& comm, const ekat::ParameterList& params);

    // The type of subcomponent
    AtmosphereProcessType type () const override { return AtmosphereProcessType::Physics; }

    // The name of the subcomponent
    std::string name () const override { return "packagename"; }

    void set_grids(
      const std::shared_ptr<const GridsManager> grids_manager) override;
    
    // Define the protected functions, usually at least initialize_impl, run_impl
    // and finalize_impl, but others could be included.  See
    // eamxx_template_process_interface.cpp for definitions of each of these.
    #ifndef KOKKOS_ENABLE_CUDA
    protected:
    #endif
    void initialize_impl(const RunType run_type) override;
    void run_impl(const double dt) override;
    void finalize_impl() override;

    // Computes bytes needed in buffers
    size_t requested_buffer_size_in_bytes() const;

    // Set the variables using memory provided by the ATMBufferManager. Needed for Fortran?
    void init_buffers(const ATMBufferManager &buffer_manager);

    // Keep track of field dimensions
    std::shared_ptr<const AbstractGrid> m_grid;
    Int m_num_cols;
    Int m_num_levs;

    // Parameters struct to pass through fortran bridge
    PackagenameFunc::params params;

  }; // class Packagename

  } // namespace scream

  #endif // SCREAM_PACKAGENAME_HPP
  ```

3. `packagename_functions.hpp`
  ```cpp
  #ifndef PACKAGENAME_FUNCTIONS_HPP
  #define PACKAGENAME_FUNCTIONS_HPP

  #include "share/core/eamxx_types.hpp"

  #include <ekat_pack_kokkos.hpp>
  #include <ekat_workspace.hpp>

  namespace scream {
  namespace packagename {

  template <typename ScalarT, typename DeviceT>
  struct PackagenameFunctions
  {

    //
    // ------- Types --------
    //

    using Scalar = ScalarT;
    using Device = DeviceT;

    template <typename S>
    using BigPack = ekat::Pack<S,SCREAM_PACK_SIZE>;
    template <typename S>
    using SmallPack = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;

    using Pack = BigPack<Scalar>;
    using Spack = SmallPack<Scalar>;

    using KT = KokkosTypes<Device>;
    using MemberType = typename KT::MemberType;

    template <typename S> using view_1d   = typename KT::template view_1d<S>;
    template <typename S> using view_2d   = typename KT::template view_2d<S>;
    template <typename S> using view_2dl  = typename KT::template lview<S**>;
    template <typename S> using uview_1d  = typename ekat::template Unmanaged<view_1d<S> >;
    template <typename S> using uview_2d  = typename ekat::template Unmanaged<view_2d<S> >;
    template <typename S> using uview_2dl = typename ekat::template Unmanaged<view_2dl<S> >;
    // ----------------------------------------
    // Structs
    struct params {
      // From field manager
      view_2d<Spack>  p_mid;
      view_2d<Spack>  T_mid;

      // For fortran
      uview_2dl<Real>  f_p_mid;
      uview_2dl<Real>  f_T_mid;

      // Set number of variables for ATMBufferManager
      static constexpr int num_2d_midlv_c_views = 2;
      static constexpr int num_2d_midlv_f_views = 2;

      // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
      void init(int ncol_in, int pver_in) {
        Real init_fill_value = -999;

        for (int i=0; i<ncol_in; ++i) {
          for (int j=0; j<pver_in; ++j) {
            f_p_mid(i,j) = init_fill_value;
            f_T_mid(i,j) = init_fill_value;
          }
        }
      }; // End init

      // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
      template <ekat::TransposeDirection::Enum D>
      void transpose(int ncol_in, int pver_in) {
        auto pverp = pver_in+1;
        if (D == ekat::TransposeDirection::c2f) {
          for (int i=0; i<ncol_in; ++i) {
            for (int j=0; j<pver_in; ++j) {
              f_p_mid(i,j) = p_mid(i,j/Spack::n)[j%Spack::n];
              f_T_mid(i,j) = T_mid(i,j/Spack::n)[j%Spack::n];
            }
          }
        }
        if (D == ekat::TransposeDirection::f2c) {
          for (int i=0; i<ncol_in; ++i) {
            // mid-point level variables
            for (int j=0; j<pver_in; ++j) {
              p_mid(i,j/Spack::n)[j%Spack::n] = f_p_mid(i,j);
              T_mid(i,j/Spack::n)[j%Spack::n] = f_T_mid(i,j);
            }
          }
        }
      }; // End transpose

    }; // end Struct params

  }; // struct PackagenameFunctions

  } // namespace packagename
  } // namespace scream

  #endif // PACKAGENAME_FUNCTIONS_HPP
  ```

4. `CMakeLists.txt`
  ```cmake
    set(PACKAGENAME_F90_SRCS
      # ----------------------------------------------------------------------------
      # EAMxx side C++ bridge methods
      ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge/packagename_eamxx_bridge.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge/packagename_eamxx_bridge_main.F90
      # list others as needed

      # ----------------------------------------------------------------------------
      # Fortran PACKAGENAME code - legacy (DNE)
      # ${PATH_TO_LEGACY_PACKAGENAME}/packagename_main.F90
      # List all needed 

      # ----------------------------------------------------------------------------
      # shared fotran modules
      ${SCREAM_BASE_DIR}/../../share/util/shr_sys_mod.F90
      ${SCREAM_BASE_DIR}/../../share/util/shr_kind_mod.F90
      ${SCREAM_BASE_DIR}/../../share/util/shr_assert_mod.F90.in

      # ----------------------------------------------------------------------------
      # misc other fortran modules if needed
      ${SCREAM_BASE_DIR}/../eam/src/utils/cam_abortutils.F90
      ${SCREAM_BASE_DIR}/../eam/src/control/cam_logfile.F90
      # ----------------------------------------------------------------------------
    )


    # List of all cpp source files 
    set(PACKAGENAME_SRCS
      eamxx_packagename_process_interface.cpp
      # packagename.cpp
    )

    # List of all hpp header files
    set(PACKAGENAME_HEADERS
      eamxx_packagename_process_interface.hpp
      packagename_functions.hpp
    )

    # Adds the library to eamxx_physics 
    add_library(packagename ${PACKAGENAME_F90_SRCS} ${PACKAGENAME_SRCS})

    SET(USE_OPENACC TRUE)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OPENACC_Fortran_FLAGS}")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OPENACC_Linker_FLAGS}")

    set_target_properties(packagename PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
    )

    target_compile_definitions(packagename PUBLIC EAMXX_HAS_PACKAGENAME)

    target_include_directories(packagename PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge
      ${CMAKE_CURRENT_BINARY_DIR}/modules
      # ${PATH_TO_LEGACY_PACKAGENAME}
    )

    set_target_properties(packagename PROPERTIES LINKER_LANGUAGE Fortran)
    target_link_libraries(packagename Kokkos::kokkos)

    target_link_libraries(packagename eamxx_physics_share scream_share)
    target_compile_options(packagename PUBLIC)

    if (TARGET eamxx_physics)
      # Add this library to eamxx_physics
      target_link_libraries(eamxx_physics INTERFACE packagename)
    endif()
  ```

### Fortran bridge

In `components/eamxx/src/physics/PACKAGENAME` add a subdirectory `fortran_bridge`. In this subdirectory add:

1. `packagename_eamxx_bridge.cpp`
  ```cpp
  #include "packagename_eamxx_bridge.hpp"
  #include "packagename_functions.hpp"
  #include "share/core/eamxx_types.hpp"

  #include <ekat_pack_kokkos.hpp>
  #include <ekat_workspace.hpp>

  using scream::Real;
  using scream::Int;

  // A C++ interface to packagename fortran calls and vice versa

  extern "C" {
      void set_log_file_name_f90_c(const char** fname);

      void packagename_eamxx_bridge_init_c(Int pcols, Int pver );

      void packagename_eamxx_bridge_run_c(Int pcols, Real* p_mid, Real* T_mid);
  } // extern "C" : end _c decls

  namespace scream {
      namespace packagename {

      void set_log_file_name_f90(const char** logname){
          set_log_file_name_f90_c(logname);
      }

      void packagename_eamxx_bridge_init( Int pcols, Int pver ){
          packagename_eamxx_bridge_init_c( pcols, pver );
      }

      void packagename_eamxx_bridge_run( Int pcols, Int pver, PackagenameFunc::params &params){ 
          //----------------------------------------------------------------------------
          // Need to transpose to match how Fortran handles things
          params.transpose<ekat::TransposeDirection::c2f>(pcols,pver);

          packagename_eamxx_bridge_run_c(pcols, params.f_p_mid.data(), params.f_T_mid.data()); 

          // Transpose back to C++ convention
          params.transpose<ekat::TransposeDirection::f2c>(pcols,pver);

          //----------------------------------------------------------------------------
          }

      // end _c impls

      } // namespace packagename
  } // namespace scream
  ```

2. `packagename_eamxx_bridge.hpp`
  ```cpp
  #include "share/core/eamxx_types.hpp"

  #include <array>
  #include <utility>
  #include <memory>   // for shared_ptr

  #include "packagename_functions.hpp"

  // Bridge functions to call fortran version of Packagename functions from C++

  namespace scream {
  namespace packagename {

      using PackagenameFunc = packagename::PackagenameFunctions<Real, DefaultDevice>;

      // Glue functions to call fortran from from C++ with the Data struct
      void packagename_eamxx_bridge_init( Int pcols, Int pver );
      void packagename_eamxx_bridge_run( Int pcols, Int pver, PackagenameFunc::params &params); 
      void set_log_file_name_f90(const char** logname);

      extern "C" { // _f function decls
      }

      }  // namespace packagename
  }  // namespace scream
  ```

3. `packagename_eamxx_bridge_main.F90` and any other `*.F90` modules you need:
  ```fortran
    module packagename_eamxx_bridge_main

    use iso_c_binding
    use openacc
    use cam_logfile,   only: iulog ! kinds instead of cam_logfile?
    use shr_sys_mod,   only: shr_sys_flush
    !-----------------------------------------------------------------------------
    implicit none
    private
    !-----------------------------------------------------------------------------
    ! public methods
    public :: packagename_eamxx_bridge_init_c
    public :: packagename_eamxx_bridge_run_c
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

  subroutine packagename_eamxx_bridge_init_c( pcol_in, pver_in ) bind(C, name="packagename_eamxx_bridge_init_c")
    ! Define uses here
    !-----------------------------------------------------------------------------
    ! Arguments
    integer(kind=c_int), value, intent(in) :: pcol_in
    integer(kind=c_int), value, intent(in) :: pver_in

    ! Set dimensions of fields
    pcols = pcol_in
    pver  = pver_in

    return
  end subroutine packagename_eamxx_bridge_init_c

  !===================================================================================================

  subroutine packagename_eamxx_bridge_run_c( ncol, p_mid, T_mid ) bind(C, name="packagename_eamxx_bridge_run_c")
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

    ! Do stuff here - scale both p_mid and T_mid by 0.5
    !$acc parallel deviceptr(p_mid, T_mid)
    !$acc loop gang vector collapse(2)
    do i = 1,pcols
      do k = 1,pver
        p_mid(i,k) = p_mid(i,k) * 0.5
        T_mid(i,k) = T_mid(i,k) * 0.5
      end do
    end do
    !$acc end parallel

    ! get max value to print out
    p_mid_max = MAXVAL(p_mid) 
    T_mid_max = MAXVAL(T_mid) 

    ! TODO - need to fix this in terms of running in parallel. also newline?
    write(iulog,*) "Fortran, max value of 0.5 p_mid ", p_mid_max
    write(iulog,*) "Fortran, max value of 0.5 T_mid ", T_mid_max

    !$acc parallel deviceptr(p_mid, T_mid)
    !$acc loop gang vector collapse(2)
    do i = 1,pcols
      do k = 1,pver
        p_mid(i,k) = p_mid(i,k) * 2.0
        T_mid(i,k) = T_mid(i,k) * 2.0
      end do
    end do
    !$acc end parallel

    ! get max value to print out
    p_mid_max = MAXVAL(p_mid) 
    T_mid_max = MAXVAL(T_mid) 

    ! TODO - need to fix this in terms of running in parallel. also newline?
    write(iulog,*) "Fortran, max value of 2 p_mid ", p_mid_max
    write(iulog,*) "Fortran, max value of 2 T_mid ", T_mid_max
    
    return
  end subroutine packagename_eamxx_bridge_run_c

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

      ! Create the log file on root rank...
      open (unit=iulog,file=trim(log_fname), &
            action='WRITE', access='SEQUENTIAL', position="append")
      
      write(iulog,*) " ---- PACKAGENAME TEST ----"
      flush(iulog)

    endif
  end subroutine set_log_file_name_f90_c
  !===================================================================================================

  end module packagename_eamxx_bridge_main
  ```

Note: I think this is the bare bones for what would be needed to set this bridge up. Will update as needed once I actually start setting this up.

## Hooking up the noodles 

In `components/eamxx/src/physics/CMakeLists.txt`:

```cmake
add_subdirectory(packagename)
```

In `components/eamxx/src/mct_coupling/CMakeLists.txt`:

```cmake
# add new physics package to list here
set (SCREAM_LIBS
     scream_share
     scream_control
     ${dynLibName}
     p3
     shoc
     scream_rrtmgp
     eamxx_cosp
     cld_fraction
     spa
     iop_forcing
     nudging
     tms
     packagename 
     )
```

In `components/src/physics/register_physics.hpp`: 

```cpp
#ifndef EAMXX_HAS_PACKAGENAME
#include "eamxx/src/physics/packagename/eamxx_packagename_process_interface.hpp"
#endif

// later in the file

#ifndef EAMXX_HAS_PACKAGENAME
  proc_factory.register_product("packagename", &create_atmosphere_process<PACKAGENAME>);
```

In `components/eamxx/cime_config/namelist_defaults_scream.xml`, add:

```xml
    <!-- Packagename -->
    <packagename inherit="atm_proc_base">
      <!-- Add in any required variables for the package here -->
    </packagename>
    <!-- later if needed  -->
    <initial_conditions>
    <variable               >default_value</variable>
    <!-- Later -->

    <mac_aero_mic inherit="atm_proc_group">
      <atm_procs_list>shoc,cld_fraction,spa,p3,packagename</atm_procs_list>
      <atm_procs_list hgrid=".*pg2">tms,shoc,cld_fraction,spa,p3,packagename</atm_procs_list> <!-- TMS only available for PG2 -->
      <atm_procs_list COMPSET=".*SCREAM%MAM4xx.*" hgrid=".*pg2">tms,shoc,cld_fraction,mam4_aci,p3,packagename</atm_procs_list> <!-- TMS only available for PG2 -->
      <atm_procs_list COMPSET=".*SCREAM.*noAero">shoc,cld_fraction,p3,packagename</atm_procs_list>
      <atm_procs_list COMPSET=".*SCREAM.*noAero" hgrid=".*pg2">tms,shoc,cld_fraction,p3,packagename</atm_procs_list> <!-- TMS only available for PG2 -->
      ...
    </mac_aero_mic>
```

An alternative to adding the new package to the `atm_procs_list` is to add it at runtime as: 

```bash
    ./atmchange mac_aero_mic::atm_procs_list+=packagename
```

--- 

## Kokkos and GPU performance 

TBD