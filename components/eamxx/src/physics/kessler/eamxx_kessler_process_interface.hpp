#ifndef SCREAM_KESSLER_PROCESS_HPP
#define SCREAM_KESSLER_PROCESS_HPP

#include "physics/kessler/kessler_functions.hpp"
#include "share/atm_process/atmosphere_process.hpp"

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

class Kessler : public AtmosphereProcess
{
  using KesslerFunc = kessler::KesslerFunctions<Real, DefaultDevice>;
  using KesslerData = KesslerFunc::KesslerData;
  using Pack        = ekat::Pack<Real, SCREAM_PACK_SIZE>;
  using view_2d     = KesslerFunc::view_2d<Pack>;
  using uview_2d    = ekat::Unmanaged<view_2d>;
  using view_2d_int = KesslerFunc::view_2d<Pack>;
  using uview_2d_int = ekat::Unmanaged<view_2d_int>;

public:

  Kessler (const ekat::Comm& comm, const ekat::ParameterList& params);

  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }
  std::string name () const { return "kessler"; }

  void create_requests ();

  // Buffer for intermediate / scratch Pack views
  struct Buffer {
    // 2D midpoint scratch arrays
    static constexpr int num_2d_mid = 5; // exner, dz, z_mid, rho, theta
    uview_2d exner, dz, z_mid, rho, theta;

    // 2D interface scratch array
    static constexpr int num_2d_int = 1; // z_int
    uview_2d_int z_int;
  };

#ifndef KOKKOS_ENABLE_CUDA
protected:
#endif

  void run_impl (const double dt);

protected:

  void initialize_impl (const RunType run_type);
  void finalize_impl   ();

  size_t requested_buffer_size_in_bytes () const;
  void   init_buffers (const ATMBufferManager& buffer_manager);

  Buffer m_buffer;

  int  m_ncols;
  int  m_nlevs;

  KesslerData m_kd;

  std::shared_ptr<const AbstractGrid> m_grid;

}; // class Kessler

} // namespace scream

#endif // SCREAM_KESSLER_PROCESS_HPP
