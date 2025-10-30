# Adding new package to EAMxx

### The code

Create a new folder in `components/eamxx/src/physics/` with the name of your package: `components/eamxx/src/physics/PACKAGENAME`.

In this folder, we must add the following files: 

1. `eamxx_packagename_process_interface.cpp`
  ```cpp
  #include "eamxx_packagename_process_interface.hpp"
  #include "share/property_checks/field_within_interval_check.hpp"
  #include "share/field/field_utils.hpp"

  #include <ekat_assert.hpp>
  #include <ekat_units.hpp>

  #include <array>

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

    // if we set any starting values for the code we can get them here with:
    // var = params.get<std::type>("var_name");
  }

  // =========================================================================================
  void Packagename::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
  {
    // Note that we inherit the logger used by the overall atmospheric process
    m_atm_logger->info("[EAMxx] Packagename processes set grids");

    // Set units
    // auto nondim = Units::nondimensional();
    constexpr auto K = ekat::units::K;
    constexpr auto Pa = ekat::units::Pa;

    // specify which grid to use
    m_grid = grids_manager->get_grid("physics");
    const auto& grid_name = m_grid->name();
    m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
    m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

    FieldLayout scalar3d_layout_mid = m_grid->get_3d_scalar_layout(true);

    constexpr int ps = 1; // constexpr int ps = Pack::n; for multilevel fields
    add_field<Updated>("T_mid", scalar3d_layout_mid, K, grid_name, ps);
    add_field<Updated>("p_mid", scalar3d_layout_mid, Pa, grid_name, ps);

    // Set of fields used strictly as input (Required)

    add_field<Required>("field1_in", m_grid, kg/kg, ps);

    // Set of fields used strictly as output (Computed)
    add_field<Computed>("field1_out", scalar3d_layout_mid, Pa, grid_name,ps);

    // Set of fields used as input and output (Updated)
    add_field<Updated>("field1_updated", scalar3d_layout_mid, Pa,grid_name, ps);

    // Gather parameters from parameter list:
    packagename_var1 = m_params.get<double>("packagename_var1",1e-12);  // Default = 1e-12
  }

  // =========================================================================================
  void Packagename::initialize_impl (const RunType /* run_type */)
  {
    // NOTE: run_type tells us if this is an initial or restarted run,

    // Set any universal things such as constants or masks (see Pompei example)

    m_atm_logger->info("[EAMxx] packagename processes initialize_impl: ");

    // This is where we can setup Kokkos functions 
    // This is where we setup any starting physics
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
    auto field1_in = get_field_in("field1_in"); // get_field_in() sets the field as read-only regardless of if it is required/computed/updated
    auto field1_out = get_field_out("field1_out");
    auto field1_updated = get_field_out("field1_updated"); // use this for updated fields too 

    // Logger call 
    m_atm_logger->info("[EAMxx] packagename run_impl: ");

    // Scale and print p_mid as an example
    auto p_mid_max = field_max<Real>(p_mid);
    m_atm_logger->info("\t max value for updated p_mid: "+ std::to_string(p_mid_max));

    p_mid.scale(2.0);
    p_mid_max = field_max<Real>(p_mid);
    m_atm_logger->info("\t max value for p_mid scaled by 2x: "+ std::to_string(p_mid_max));

    p_mid.scale(0.5);
    p_mid_max = field_max<Real>(p_mid);
    m_atm_logger->info("\t max value for p_mid scaled by 0.5x: "+ std::to_string(p_mid_max));

  }

  // =========================================================================================
  void Packagename::finalize_impl()
  {
    // Do nothing usually
  }
  // =========================================================================================

  } // namespace scream
  ```

2. `eamxx_packagename_process_interface.hpp`
  ```cpp
  #ifndef SCREAM_PACKAGENAME_HPP
  #define SCREAM_PACKAGENAME_HPP

  #include "physics/packagename/packagename_functions.hpp"
  #include "share/atm_process/atmosphere_process.hpp"

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
    using Spack           = PackagenameFunc::Spack;
    using Smask           = PackagenameFunc::Smask;
    using Pack            = ekat::Pack<Real,Spack::n>;

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

    // Keep track of field dimensions and the iteration count
    Int m_num_cols;
    Int m_num_levs;

    // Parameters here if needed. 
    // These can be set in namelist_defaults_scream.xml (see below) for default values and etc
    Real packagename_var1;
    Real packagename_var2;

    std::shared_ptr<const AbstractGrid> m_grid;
  }; // class Packagename

  } // namespace scream

  #endif // SCREAM_PACKAGENAME_HPP
  ```

3. `CMakeLists.txt`
  ```cmake
    # List of all cpp source files needed for the package to work
    set(PACKAGENAME_SRCS
      eamxx_packagename_process_interface.cpp
      # packagename.cpp
    )

    # List of all hpp header files needed for the package to work
    set(PACKAGENAME_HEADERS
      eamxx_packagename_process_interface.hpp
      packagename_functions.hpp
    )

    # Adds the library to eamxx_physics 
    add_library(packagename ${PACKAGENAME_SRCS})
    target_compile_definitions(packagename PUBLIC EAMXX_HAS_PACKAGENAME)
    target_link_libraries(packagename eamxx_physics_share scream_share)
    target_compile_options(packagename PUBLIC)

    if (TARGET eamxx_physics)
      # Add this library to eamxx_physics
      target_link_libraries(eamxx_physics INTERFACE packagename)
    endif()
  ```

### Hooking up the noodles 

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

## Bridges

### Fortran

In `components/eamxx/src/physics/PACKAGENAME` add a subdirectory `fortran_bridge`. In this subdirectory add:

1. `packagename_eamxx_bridge.cpp`
    ```cpp
    #include "packagename_eamxx_bridge.hpp"

    using scream::Real;
    using scream::Int;

    // A C++ interface to packagename fortran calls and vice versa

    extern "C" {packagename_eamxx_bridge_run_c( required_args);
      } // extern "C" : end _c decls

      namespace scream {
      namespace packagename {

      void packagename_eamxx_bridge_init( Int pcol, Int pver ){
        packagename_eamxx_bridge_init_c( pcol, pver );
      }

      void packagename_eamxx_bridge_run( Int ncol, Int pver,
                                ZMF::packagename_input_state& packagename_input,
                                ZMF::packagename_output_tend& packagename_output,
                                ZMF::packagename_runtime_opt& packagename_opts
                              ){
        //----------------------------------------------------------------------------
        // Need to transpose to match how Fortran handles things
        // packagename_input.transpose<ekat::TransposeDirection::c2f>(ncol,pver);
        // packagename_output.transpose<ekat::TransposeDirection::c2f>(ncol,pver);

        packagename_eamxx_bridge_run_c( required_arguments);

        // Transpose back to C++ convention
        // packagename_input.transpose<ekat::TransposeDirection::f2c>(ncol,pver);
        // packagename_output.transpose<ekat::TransposeDirection::f2c>(ncol,pver);

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

    // Bridge functions to call fortran version of ZM functions from C++

    namespace scream {
    namespace packagename {

    using ZMF = packagename::Functions<Real, DefaultDevice>;

    // Glue functions to call fortran from from C++ with the Data struct
    void packagename_eamxx_bridge_init( Int pcols, Int pver );
    void packagename_eamxx_bridge_run( Int ncol, Int pver, required_arguments);

    extern "C" { // _f function decls
    }

    }  // namespace packagename
    }  // namespace scream
    ```

3. `packagename_eamxx_bridge_main.F90` and any other `*.F90` modules you need:
    ```fortran
    module packagename_eamxx_bridge_main

      use iso_c_binding
      use cam_logfile,   only: iulog
      use shr_sys_mod,   only: shr_sys_flush
      use packagename_eamxx_bridge_params, only: masterproc, r8, pcols, pver, pverp, top_lev
      !-----------------------------------------------------------------------------
      implicit none
      private
      !-----------------------------------------------------------------------------
      ! public methods
      public :: packagename_eamxx_bridge_init_c
      public :: packagename_eamxx_bridge_run_c

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

    subroutine packagename_eamxx_bridge_init_c( pcol_in, pver_in ) bind(C)
      ! Do stuff here
    end subroutine packagename_eamxx_bridge_init_c

    subroutine packagename_eamxx_bridge_run_c( pcol_in, pver_in ) bind(C)
      ! Do stuff here
    end subroutine packagename_eamxx_bridge_run_c
  end module packagename_eamxx_bridge_main
  ```

Update `components/eamxx/src/physics/PACKAGENAME/CMakeLists.txt`:

```cmake
# Set a legacy path if needed 
set(PATH_TO_LEGACY_PACKAGENAME ${SCREAM_BASE_DIR}/path/to/code/) 

set(PACKAGENAME_F90_SRCS
  # ----------------------------------------------------------------------------
  # EAMxx side C++ bridge methods
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge/packagename_eamxx_bridge.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge/packagename_eamxx_bridge_main.F90
  # list others as needed

  # ----------------------------------------------------------------------------
  # Fortran PACKAGENAME code
  ${PATH_TO_LEGACY_PACKAGENAME}/packagename_*.F90
  # List all needed 

  # ----------------------------------------------------------------------------
  # shared fotran modules
  ${SCREAM_BASE_DIR}/../../share/util/shr_sys_mod.F90
  ${SCREAM_BASE_DIR}/../../share/util/shr_kind_mod.F90
  ${SCREAM_BASE_DIR}/../../share/util/shr_assert_mod.F90.in

  # ----------------------------------------------------------------------------
  # misc other fotran modules if needed
  
  # ----------------------------------------------------------------------------
)

# List of all cpp source files needed for the package to work
set(PACKAGENAME_CXX_SRCS
  eamxx_packagename_process_interface.cpp
)

# List of all hpp header files needed for the package to work
set(PACKAGENAME_CXX_HEADERS
  eamxx_packagename_process_interface.hpp
  packagename_functions.hpp
)

# Adds the library to eamxx_physics 
# -----------------------------------------

add_library(packagename ${PACKAGENAME_F90_SRCS} ${PACKAGENAME_CXX_SRCS})

set_target_properties(packagename PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)
target_compile_definitions(packagename PUBLIC EAMXX_HAS_PACKAGENAME)

target_include_directories(packagename PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  ${PATH_TO_LEGACY_PACKAGENAME}
)

target_link_libraries(packagename eamxx_physics_share scream_share)
target_compile_options(packagename PUBLIC)

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE packagename)
endif()
```

To finish linking everything together, we update: 

```cpp
  // =========================================================================================
  void Packagename::initialize_impl (const RunType /* run_type */)
  {
    // NOTE: run_type tells us if this is an initial or restarted run,

    // Set any universal things such as constants or masks (see Pompei example)

    m_atm_logger->info("[EAMxx] packagename processes initialize_impl: ");

    // This is where we can setup Kokkos functions 
    // This is where we setup any starting physics
    packagename::packagename_eamxx_bridge_init(m_pcol, m_nlev);
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
    auto field1_in = get_field_in("field1_in"); // get_field_in() sets the field as read-only regardless of if it is required/computed/updated
    auto field1_out = get_field_out("field1_out");
    auto field1_updated = get_field_out("field1_updated"); // use this for updated fields too 

    // Logger call 
    m_atm_logger->info("[EAMxx] packagename run_impl: ");

    packagename::packagename_eamxx_bridge_run(m_pcol, m_nlev, required_arguments);

    // Do anything else we need to do

  }

  // =========================================================================================
  void Packagename::finalize_impl()
  {
    // Do nothing usually
  }
  // =========================================================================================
```

Note: I think this is the bare bones for what would be needed to set this bridge up. Will update as needed once I actually start setting this up.

### Python

In `components/eamxx/src/physics/PACKAGENAME` add `packagename.py`:

```python
# Imports 
import numpy as np

# Any initialization step can be done here
# This method is called during Packagename::initialize_impl
def init ():
    pass

# Probably what we would call during run_impl
def main (required_arguments):
    # Do whatever physics needed 

# Define any other needed functions

```

Update `eamxx_packagename_process_interface.cpp` as follows:

```cpp
  // All other imports here

  #ifdef EAMXX_HAS_PYTHON
  #include "share/atm_process/atmosphere_process_pyhelpers.hpp"
  #endif

  void Packagename::initialize_impl (const RunType /* run_type */)
  {
    // Previously defined initialization steps 

    #ifdef EAMXX_HAS_PYTHON
      if (has_py_module()) {
        try {
          py_module_call("init");
        } catch (const pybind11::error_already_set& e) {
          std::cout << "[Packagename::initialize_impl] Error! Something went wrong while calling the python module's function 'init'.\n"
                      " - module name: " + m_params.get<std::string>("py_module_name") + "\n"
                      " - pybind11 error: " + std::string(e.what()) + "\n";
          throw e;
        }

      }
    #endif
  }

  void Packagename::run_impl (const double /* dt */)
{
  // Pull in variables.
    auto T_mid          = get_field_out("T_mid");
    auto p_mid          = get_field_out("p_mid");
    auto field1_in      = get_field_in("field1_in"); // get_field_in() sets the field as read-only regardless of if it is required/computed/updated
    auto field1_out     = get_field_out("field1_out");
    auto field1_updated = get_field_out("field1_updated"); // use this for updated fields too 

  #ifdef EAMXX_HAS_PYTHON
    if (has_py_module()) {
      // For now, we run Python code only on CPU
      const auto& T_mid          = get_py_field_host("T_mid");
      const auto& p_mid          = get_py_field_host("p_mid");
      const auto& field1_in      = get_py_field_host("field1_in");
      const auto& field1_out     = get_py_field_host("field1_out");
      const auto& field1_updated = get_py_field_host("field1_updated");

      // Sync input to host
      field1_in.sync_to_host();

      // Read in any parameters
      // double param_name = m_params.get<double>("param_name");

      try {
        py_module_call("main", required_arguments);
      } catch (const pybind11::error_already_set& e) {
        std::cout << "[Packagename::run_impl] Error! Something went wrong while calling the python module's function 'main'.\n"
                    " - module name: " + m_params.get<std::string>("py_module_name") + "\n"
                    " - pybind11 error: " + std::string(e.what()) + "\n";
        throw e;
      }

      // Sync outputs to dev
      T_mid.sync_to_dev();
      p_mid.sync_to_dev();
      field1_in.sync_to_dev();
      field1_out.sync_to_dev();
      field1_updated.sync_to_dev();
    } else
    #endif
    // Here is where the normal calls w/out python would go
  }
```

WIP - how to run Python code on the GPU?

## Kokkos and GPU performance 

TBD