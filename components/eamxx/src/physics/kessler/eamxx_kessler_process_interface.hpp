#ifndef SCREAM_KESSLER_HPP
#define SCREAM_KESSLER_HPP

#include "physics/kessler/kessler_functions.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "share/atm_process/ATMBufferManager.hpp"

#include "share/physics/physics_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
// #include "share/physics/eamxx_common_physics_functions_impls.hpp"

#include <ekat_parameter_list.hpp>

#include <string>

namespace scream
{

/*
 * The class responsible to do Kessler microphysics
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
*/

class KesslerMicrophysics : public AtmosphereProcess
{
public:
  using KT  = ekat::KokkosTypes<DefaultDevice>;
  using KMF = kessler::KesslerMicrophysicsFunctions<Real, DefaultDevice>;
  using PF  = scream::PhysicsFunctions<DefaultDevice>;
  using PC  = scream::physics::Constants<Real>;

  using Scalar = KMF::Scalar;
  using Pack = KMF::Pack;

  // Constructors
  KesslerMicrophysics (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const override { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const override { return "kessler"; }

  // Create grid-dependent field requests
  void create_requests() override;
  // Old method 
  // void set_grids(
  //   const std::shared_ptr<const GridsManager> grids_manager) override;
  
  // Define the protected functions, usually at least initialize_impl, run_impl
  // and finalize_impl, but others could be included.  See
  // eamxx_template_process_interface.cpp for definitions of each of these.
  #ifndef KOKKOS_ENABLE_CUDA
    protected:
  #endif
    void initialize_impl(const RunType run_type) override;
    void run_impl(const double dt) override;
  protected:
    void finalize_impl() override;

    // Computes bytes needed in buffers
    size_t requested_buffer_size_in_bytes() const;

    // Set the variables using memory provided by the ATMBufferManager. Needed for Fortran?
    void init_buffers(const ATMBufferManager &buffer_manager);

    // Keep track of field dimensions
    std::shared_ptr<const AbstractGrid> m_grid;
    int m_num_cols;
    int m_num_levs;

    // Parameters structs to pass through fortran bridge
    KMF::params_helpers params_helpers; // Helper variables 
    KMF::params_computed params_computed; // Variables computed in the run_impl function

}; // class Kessler

} // namespace scream

#endif // SCREAM_KESSLER_HPP