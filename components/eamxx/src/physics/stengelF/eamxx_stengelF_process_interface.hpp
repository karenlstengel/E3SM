#ifndef SCREAM_STENGELF_HPP
#define SCREAM_STENGELF_HPP

#include "physics/stengelF/stengelF_functions.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "share/atm_process/ATMBufferManager.hpp"

#include <ekat_parameter_list.hpp>

#include <string>

namespace scream
{

/*
 * The class responsible to do stengelF physics
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
*/

class StengelF : public AtmosphereProcess
{
public:
  using StengelFFunc = stengelF::StengelFFunctions<Real, DefaultDevice>;
  using Spack        = StengelFFunc::Spack;
  using Pack         = ekat::Pack<Real,Spack::n>;

  // Constructors
  StengelF (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const override { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const override { return "stengelF"; }

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
  StengelFFunc::params params;

}; // class StengelF

} // namespace scream

#endif // SCREAM_STENGELF_HPP
