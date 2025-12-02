# Adding new package to EAMxx

### The code with Python bridge

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

    // if we set any starting values for the code we can get them here with:
    // var = params.get<std::type>("var_name");
  }

  // =========================================================================================
  void Packagename::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
  {
    // using namespace ekat::units;
    // using namespace ShortFieldTagsNames;

    m_atm_logger->info("[EAMxx] Packagename processes set grids");

    // // The units of mixing ratio Q are technically non-dimensional.
    // // Nevertheless, for output reasons, we like to see 'kg/kg'.
    // auto nondim = Units::nondimensional();
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
    // NOTE: run_type tells us if this is an initial or restarted run,

    // Set any universal things such as constants or masks (see Pompei example)

    m_atm_logger->info("[EAMxx] packagename processes initialize_impl: ");

    #ifdef EAMXX_HAS_PYTHON
        if (has_py_module()) {
          try {
            py_module_call("init");
          } catch (const pybind11::error_already_set& e) {
            std::cout << "[packagename::initialize_impl] Error! Something went wrong while calling the python module's function 'init'.\n"
                        " - module name: " + m_params.get<std::string>("py_module_name") + "\n"
                        " - pybind11 error: " + std::string(e.what()) + "\n";
            throw e;
          }

        }
      #endif

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
    auto p_mid   = get_field_out("p_mid"); // work?

    // TODO - scale and print 
    m_atm_logger->info("[EAMxx] packagename run_impl: ");

    // Get log file name 
    const char* logname = m_atm_logger->get_logfile_name().c_str();
    
    #ifdef EAMXX_HAS_PYTHON
      if (has_py_module()) {
        // For now, we run Python code only on CPU
        const auto& T_mid          = get_py_field_host("T_mid");
        const auto& p_mid          = get_py_field_host("p_mid");

        // Sync input to host
        //field1_in.sync_to_host();

        // Read in any parameters
        // double param_name = m_params.get<double>("param_name");

        try {
          py_module_call("main", p_mid, T_mid, &logname);
        } catch (const pybind11::error_already_set& e) {
          std::cout << "[packagename::run_impl] Error! Something went wrong while calling the python module's function 'main'.\n"
                      " - module name: " + m_params.get<std::string>("py_module_name") + "\n"
                      " - pybind11 error: " + std::string(e.what()) + "\n";
          throw e;
        }

        // Sync outputs to dev
        T_mid.sync_to_dev();
        p_mid.sync_to_dev();
      } else
      #endif

  }

  // =========================================================================================
  void Packagename::finalize_impl()
  {
    // Do nothing
    m_atm_logger->info("[EAMxx] Packagename processes clean up.");
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

    // Parameters for fun
    Real packagename_var1;
    Real packagename_var2;

    std::shared_ptr<const AbstractGrid> m_grid;
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

    using Mask = ekat::Mask<BigPack<Scalar>::n>;
    using Smask = ekat::Mask<SmallPack<Scalar>::n>;

    using KT = KokkosTypes<Device>;
    using MemberType = typename KT::MemberType;

    template <typename S>
    using view_1d = typename KT::template view_1d<S>;
    template <typename S>
    using view_2d = typename KT::template view_2d<S>;
    template <typename S>
    using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;

  }; // struct Functions

  } // namespace packagename
  } // namespace scream

  #endif // PACKAGENAME_FUNCTIONS_HPP
  ```

4. `CMakeLists.txt`
  ```cmake
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

  # Python setup - from the ML correction package
  include(ScreamUtils)
      if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.11.0")
      message(STATUS "Downloading Pybind11")
      include(FetchContent)

      FetchContent_Declare(pybind11 GIT_REPOSITORY https://github.com/pybind/pybind11.git GIT_TAG v3.0.1)
      FetchContent_MakeAvailable(pybind11)
  else()
      message(FATAL_ERROR "pybind11 is missing. Use CMake >= 3.11 or download it")
  endif()
  find_package(Python REQUIRED COMPONENTS Interpreter Development)

  # Adds the library to eamxx_physics 
  add_library(packagename ${PACKAGENAME_SRCS})
  target_compile_definitions(packagename PUBLIC EAMXX_HAS_PACKAGENAME)

  target_include_directories(packagename PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  target_include_directories(packagename SYSTEM PUBLIC ${PYTHON_INCLUDE_DIRS})

  target_link_libraries(packagename eamxx_physics_share scream_share pybind11::pybind11 Python::Python)
  target_compile_options(packagename PUBLIC)

  if (TARGET eamxx_physics)
    # Add this library to eamxx_physics
    target_link_libraries(eamxx_physics INTERFACE packagename)
  endif()
  ```
---

## Python

In `components/eamxx/src/physics/PACKAGENAME` add `packagename.py`:

```python
  # Imports 
  import numpy as np

  # Any initialization step can be done here
  # This method is called during Packagename::initialize_impl
  def init ():
      pass

  # Probably what we would call during run_impl
  def main (p_mid, T_mid, fname=None):
      # Print out max value of p_mid and T_mid
      if fname is not None:
          with open(fname, 'a') as log_file:
              log_file.write("\t Python - Max value for p_mid: ", np.max(p_mid))
              log_file.write("\t Python - Max value for T_mid: ", np.max(T_mid))

              # scale both by 0.5
              p_mid = p_mid * 0.5
              T_mid = T_mid * 0.5

              # log max values
              log_file.write("\t Python - Max value for p_mid scaled by 0.5: ", np.max(p_mid))
              log_file.write("\t Python - Max value for T_mid scaled by 0.5: ", np.max(T_mid))

              #scale both by 2.0
              p_mid = p_mid * 2.0
              T_mid = T_mid * 2.0

              # log max values
              log_file.write("\t Python - Max value for p_mid scaled by 2.0: ", np.max(p_mid))
              log_file.write("\t Python - Max value for T_mid scaled by 2.0: ", np.max(T_mid))

      else:
          # scale both by 0.5
          p_mid = p_mid * 0.5
          T_mid = T_mid * 0.5

          #scale both by 2.0
          p_mid = p_mid * 2.0
          T_mid = T_mid * 2.0

  # Define any other needed functions
```

To check: we can use `get_py_field_dev()` to get fields on the device and also use `sync_to_host()`. I am not sure why the examples specify to run python implementations only on the CPU.

TODO - figure out how to pass the logfile name into python function

---

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

Finally, add the Python binary location in 'components/eamxx/cmake/machine-files/derecho.cmake':

```cmake
# Set Python info
# need to have pybind11 and mpi4py installed and Python >= 3.9.2
OPTION(EAMXX_ENABLE_PYTHON "" ON)
set(Python_EXECUTABLE /path/to/python)
```
