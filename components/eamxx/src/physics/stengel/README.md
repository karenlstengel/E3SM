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
    using PNF = packagename::PNFtions<Real, DefaultDevice>;
    using Spack           = PNF::Spack;
    using Smask           = PNF::Smask;
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

## Kokkos and GPU performance 

At runtime we need to add/update the following to our submission script:

```bash
  MYCOMPILER=gnugpu # previously intel

  # other setup calls

  ./xmlchange NTASKS=4
  ./xmlchange NTHRDS=1
  ./xmlchange NGPUS_PER_NODE=4
  ./xmlchange GPU_TYPE=a100 # NVIDIA A100 GPUs in Derecho
  ./xmlchange OPENACC_GPU_OFFLOAD=FALSE
  ./xmlchange OPENMP_GPU_OFFLOAD=FALSE
  ./xmlchange KOKKOS_GPU_OFFLOAD=TRUE
  ./xmlchange OVERSUBSCRIBE_GPU=FALSE
  ./xmlchange ROOTPE='0'
  ./xmlchange DOUT_S=false

  # Previously:
  # ./xmlchange NTASKS=128
  # ./xmlchange NTHRDS=1
  # ./xmlchange ROOTPE='0'  
```

In `components/eamxx/src/physics/PACKAGENAME/`, we then update the C++ as follows:
1. `eamxx_packagename_process_interface.cpp`
  ```cpp
  #include "eamxx_packagename_process_interface.hpp"
  #include "share/property_checks/field_within_interval_check.hpp"
  #include "share/field/field_utils.hpp"

  #include <ekat_assert.hpp>
  #include <ekat_units.hpp>
  #include <ekat_team_policy_utils.hpp>

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
    m_atm_logger->info("[EAMxx] Packagename processes constructor");

    // if we set any starting values for the code we can get them here with:
    // var = params.get<std::type>("var_name");
  }

  // =========================================================================================
  void Packagename::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
  {

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
    // NOTE: run_type tells us if this is an initial or restarted run,

    // Set any universal things such as constants or masks (see Pompei example)

    m_atm_logger->info("[EAMxx] packagename processes initialize_impl: ");

    // This is where we setup any starting physics
  }

  // =========================================================================================

  // run_impl is called every timestep and where all of the physics happens
  // Inputs:
  //    - dt - the timestep for the current run step

  void Packagename::run_impl (const double /* dt */)
  {

    m_atm_logger->info("[EAMxx] packagename run_impl: ");
    
    // Run with kokkos
    #ifdef EAMXX_ENABLE_GPU

      // Pull in variables .
      auto T_mid   = get_field_out("T_mid").get_view<Spack**>();
      auto p_mid   = get_field_out("p_mid").get_view<Spack**>(); // use this for updated fields too 

      m_atm_logger->info("\t running on GPU! ");
      // TODO - need to open fields as view?
      Real T_mid_max = 0.0;
      Real p_mid_max = 0.0;

      // Kokkos function to get the max value in the field
      PNF::get_max_value(T_mid, T_mid_max);
      m_atm_logger->info("\t\t max value for updated T_mid: "+ std::to_string(T_mid_max));

      // Kokkos function to scale field by scalar value
      PNF::scale_field_2d(T_mid, 2.0);
      PNF::get_max_value(T_mid, T_mid_max);
      m_atm_logger->info("\t\t max value for T_mid scaled by 2x: "+ std::to_string(T_mid_max));

      PNF::scale_field_2d(T_mid, 0.5);
      PNF::get_max_value(T_mid, T_mid_max);
      m_atm_logger->info("\t\t max value for T_mid scaled by 0.5x: "+ std::to_string(T_mid_max));

      //------------
      // Repeat for p_mid
      PNF::get_max_value(p_mid, p_mid_max);
      m_atm_logger->info("\t\t max value for updated p_mid: "+ std::to_string(p_mid_max));

      // Kokkos function to scale field by scalar value
      PNF::scale_field_2d(p_mid, 2.0);
      PNF::get_max_value(p_mid, p_mid_max);
      m_atm_logger->info("\t\t max value for p_mid scaled by 2x: "+ std::to_string(p_mid_max));

      PNF::scale_field_2d(p_mid, 0.5);
      PNF::get_max_value(p_mid, p_mid_max);
      m_atm_logger->info("\t\t max value for p_mid scaled by 0.5x: "+ std::to_string(p_mid_max));
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
    using KT    = ekat::KokkosTypes<DefaultDevice>;
    using PNF   = packagename::PNFtions<Real, DefaultDevice>;
    using Spack = PNF::Spack;
    using Smask = PNF::Smask;
    using Pack  = ekat::Pack<Real,Spack::n>;

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

    // Keep track of field dimensions
    Int m_num_cols;
    Int m_num_levs;

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
  #include <ekat_team_policy_utils.hpp>

  namespace scream {
  namespace packagename {

  template <typename ScalarT, typename DeviceT>
  struct PNFtions
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

    // GPU/Kokkos related things
    using KT = ekat::KokkosTypes<Device>;
    using TeamPolicy  = typename KokkosTypes<Device>::TeamPolicy;

    template <typename S>
    using view_1d = typename KT::template view_1d<S>;
    template <typename S>
    using view_2d = typename KT::template view_2d<S>;

    // KOKKOS_FUNCTION 
    static void scale_field_2d(view_2d<Spack> &field, Real alpha);
    // KOKKOS_FUNCTION 
    static void get_max_value(view_2d<Spack> &field, Real &max_value);

  }; // struct Functions

  template<typename S, typename D>
  // KOKKOS_FUNCTION
  void PNFtions<S,D>::scale_field_2d(view_2d<Spack> &field, Real alpha) {

    const int ni = static_cast<int>(field.extent(0));
    const int nj = static_cast<int>(field.extent(1));

    using MDPolicy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
    MDPolicy mdp({0,0}, {ni,nj});

    // create a device mirror and copy host->device, run kernel, copy device->host
    auto field_dev = Kokkos::create_mirror_view_and_copy(DefaultDevice(), field); // or try HostDevice()
    Kokkos::parallel_for(mdp, KOKKOS_LAMBDA(const int i, const int j) {
        field_dev(i,j) *= alpha;
    });
    Kokkos::deep_copy(field, field_dev);
  };

  template<typename S, typename D>
  // KOKKOS_FUNCTION
  void PNFtions<S,D>::get_max_value(view_2d<Spack> &field, Real &max_value) {
    // compute global max over all i,j in a single parallel_reduce (no nested lambdas)
    const int ni = static_cast<int>(field.extent(0));
    const int nj = static_cast<int>(field.extent(1));
    Real global_max = -std::numeric_limits<Real>::infinity();

    using MDPolicy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
    MDPolicy mdp({0,0}, {ni,nj});

    Kokkos::parallel_reduce("GetMaxValueMD", mdp,
      KOKKOS_LAMBDA(const int i, const int j, Real &local_max) {
        const Real v = field(i,j)[0];
        if (v > local_max) local_max = v;
      },
      Kokkos::Max<Real>(global_max)
    );

    // store result on host
    max_value = global_max;
  };

  } // namespace packagename
  } // namespace scream

  #endif // PACKAGENAME_FUNCTIONS_HPP
  ```