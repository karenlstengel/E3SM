#ifndef SCREAM_STENGELP_HPP
#define SCREAM_STENGELP_HPP

#include "physics/stengelP/stengelP_functions.hpp"
#include "share/atm_process/atmosphere_process.hpp"

#include <ekat_parameter_list.hpp>

#include <string>

namespace scream
{

/*
 * The class responsible to do stengelP physics
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
*/

class StengelP : public AtmosphereProcess
{
public:
  using StengelPFunc = stengelP::StengelPFunctions<Real, DefaultDevice>;
  using Spack           = StengelPFunc::Spack;
  using Smask           = StengelPFunc::Smask;
  using Pack            = ekat::Pack<Real,Spack::n>;

  // Constructors
  StengelP (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const override { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const override { return "stengelP"; }

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
  Real stengelP_var1;
  Real stengelP_var2;

  std::shared_ptr<const AbstractGrid> m_grid;
}; // class StengelP

} // namespace scream

#endif // SCREAM_STENGELP_HPP
