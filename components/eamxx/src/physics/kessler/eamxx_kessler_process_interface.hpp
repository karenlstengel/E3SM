#ifndef SCREAM_KESSLER_PROCESS_HPP
#define SCREAM_KESSLER_PROCESS_HPP

#include "physics/kessler/kessler_functions.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "share/atm_process/ATMBufferManager.hpp"

#include <ekat_parameter_list.hpp>

#include <string>

namespace scream
{

/*
 * The class responsible for running the Kessler (1969) warm rain
 * microphysics parameterization as an EAMxx AtmosphereProcess.
 *
 * Required fields (input):
 *   T_mid         - air temperature at layer midpoints (K)
 *   p_mid         - air pressure at layer midpoints (Pa)
 *   pseudo_density- layer thickness in pressure units (Pa)
 *   phis          - surface geopotential (m^2 s^-2)
 *   qv, qc, qr    - moisture tracers (kg kg^-1, wrt dry air)
 *
 * Computed / updated fields (output):
 *   T_mid         - updated by Kessler latent heating
 *   qv, qc, qr    - updated moisture fields
 *   precl         - total precipitation rate at surface (m s^-1)
 *   relhum        - relative humidity (%)
 *
 * The AD should store exactly ONE instance of this class in its list
 * of subcomponents.
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

  template <typename S> using uview_1d = typename ekat::template Unmanaged<KT::view_1d<S>>;
  template <typename S> using uview_2d = typename ekat::template Unmanaged<KT::view_2d<S>>;

  // Constructors
  KesslerMicrophysics (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const override { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const override { return "kessler"; }

  // Create grid-dependent field requests
  void create_requests() override;

  // Buffer/workspace management: request one ATMBufferManager allocation
  // for kessler_run's persistent scratch instead of allocating and
  // freeing it on every call (mirrors P3Microphysics/SHOCMacrophysics'
  // own Buffer + requested_buffer_size_in_bytes()/init_buffers()).
  size_t requested_buffer_size_in_bytes() const override;
  void init_buffers(const ATMBufferManager& buffer_manager) override;
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

    // Keep track of field dimensions
    std::shared_ptr<const AbstractGrid> m_grid;
    Int m_num_cols;
    Int m_num_levs;

    // Parameters
    Real Cpair;
    Real Rair;
    Real latvap;
    Real pref;
    Real rhoqr;
    Real gravity;

    // Persistent scratch for kessler_run, carved from one
    // ATMBufferManager allocation in init_buffers() (see
    // KMF::Workspace in kessler_functions.hpp, which these fields
    // populate at each run_impl call). Same 6x Pack view_2d + 4x Real
    // view_1d split as KMF::Workspace, unmanaged since the memory is
    // owned by the ATMBufferManager, not this struct.
    struct Buffer {
      static constexpr int num_2d_vector = KMF::Workspace::num_2d_vector;
      static constexpr int num_1d_scalar = KMF::Workspace::num_1d_scalar;

      uview_2d<Pack> r, rhalf, velqr, sed, pc, f5;
      uview_1d<Real> dt0, mask, time_counter, precl_acc;
    };
    Buffer m_buffer;

}; // class Kessler

} // namespace scream

#endif // SCREAM_KESSLER_PROCESS_HPP
