#ifndef KESSLER_FUNCTIONS_HPP
#define KESSLER_FUNCTIONS_HPP

#include "share/core/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>

namespace scream {
namespace kessler {

/*
 * KesslerFunctions is a stateless struct used to encapsulate the
 * Kessler (1969) warm rain microphysics parameterization.
 *
 * The Kessler scheme contains three moisture categories: water vapor,
 * cloud water (liquid water that moves with the flow), and rain water
 * (liquid water that falls relative to the surrounding air).
 *
 * References:
 *   Klemp & Wilhelmson (1978), J. Atmos. Sci., 35, 1070-1096
 *   Durran & Klemp (1983), Mon. Wea. Rev., 111, 2341-2361
 */
template <typename ScalarT, typename DeviceT>
struct KesslerFunctions
{
  //
  // ------- Types --------
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  using KT = ekat::KokkosTypes<Device>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;

  //
  // ------- Constants struct --------
  //

  // Physical constants required by the Kessler parameterization.
  // pref must be supplied in hPa (the Fortran init divides the Pa input by 100).
  struct KesslerData {
    Scalar lv;    // latent heat of vaporization (J kg-1)
    Scalar pref;  // reference pressure (hPa)
    Scalar rhoqr; // density of fresh liquid water (kg m-3)
  };

  //
  // --------- Functions ---------
  //

  // Main Kessler warm rain microphysics.
  // All 2D arrays have layout (ncols, nz).
  // lyr_surf and lyr_toa are 0-based level indices; the sign of their
  // difference determines lyr_step (+1 or -1).
  // cpair and rair may vary by column and level (composition-dependent).
  static void kessler_run(
    const int ncols,
    const int nz,
    const Scalar dt,
    const int lyr_surf,
    const int lyr_toa,
    const KesslerData& kd,
    const view_2d<const Scalar>& cpair,
    const view_2d<const Scalar>& rair,
    const view_2d<const Scalar>& rho,
    const view_2d<const Scalar>& z,
    const view_2d<const Scalar>& pk,
    const view_2d<Scalar>& theta,
    const view_2d<Scalar>& qv,
    const view_2d<Scalar>& qc,
    const view_2d<Scalar>& qr,
    const view_1d<Scalar>& precl,
    const view_2d<Scalar>& relhum);

  // Save the current temperature and zero the temperature tendency accumulator.
  // Called once per physics time step before kessler_run.
  static void kessler_update_timestep_init(
    const int ncols,
    const int nz,
    const view_2d<const Scalar>& temp,
    const view_2d<Scalar>& temp_prev,
    const view_2d<Scalar>& ttend_t);

  // Back out the temperature tendency due to Kessler from theta and exner.
  // ttend_t += (theta * exner - temp_prev) / dt
  static void kessler_update_run(
    const int ncols,
    const int nz,
    const Scalar dt,
    const view_2d<const Scalar>& theta,
    const view_2d<const Scalar>& exner,
    const view_2d<const Scalar>& temp_prev,
    const view_2d<Scalar>& ttend_t);

  // Compute dry static energy: st_energy = cpair*temp + gravit*zm + phis
  static void kessler_update_timestep_final(
    const int ncols,
    const int nz,
    const Scalar gravit,
    const view_2d<const Scalar>& cpair,
    const view_2d<const Scalar>& temp,
    const view_2d<const Scalar>& zm,
    const view_1d<const Scalar>& phis,
    const view_2d<Scalar>& st_energy);

}; // struct KesslerFunctions

} // namespace kessler
} // namespace scream

// If a GPU build without relocatable device code, include the full
// implementation in every translation unit; otherwise ETI is used.
#if defined(EAMXX_ENABLE_GPU) && !defined(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE) \
                                && !defined(KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
# include "kessler_functions_impl.hpp"
#endif

#endif // KESSLER_FUNCTIONS_HPP
