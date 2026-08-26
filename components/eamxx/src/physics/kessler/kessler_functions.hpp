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
  // cpair and rair are now 1D (composition-independent) scalar constants,
  // matching the updated Fortran interface -- previously these were 2D
  // (ncols, nz) arrays even though every call filled them with the same
  // uniform dry-air value.
  struct KesslerData {
    Scalar lv;    // latent heat of vaporization (J kg-1)
    Scalar pref;  // reference pressure (hPa)
    Scalar rhoqr; // density of fresh liquid water (kg m-3)
    Scalar cpair; // specific heat of dry air at constant pressure (J kg-1 K-1)
    Scalar rair;  // gas constant of dry air (J kg-1 K-1)
  };

  //
  // ------- Workspace struct --------
  //

  // Scratch views used internally by kessler_run. Owned by the caller and
  // passed in by reference so they can be allocated once (e.g. at process
  // initialization) and reused across every call, instead of being
  // allocated and freed on the device on every single call.
  struct KesslerWorkspace {
    view_2d<Scalar> r, rhalf, velqr, sed, pc;
    view_1d<Scalar> dt0, mask, time_counter, precl_acc;
    int ncols = 0;
    int nz    = 0;

    // (Re)allocates the scratch views if the requested dimensions differ
    // from what is currently allocated. Safe to call every time step --
    // a no-op once ncols/nz stop changing (the common case).
    void init(const int ncols_in, const int nz_in) {
      if (ncols_in == ncols && nz_in == nz) return;
      ncols = ncols_in;
      nz    = nz_in;
      r            = view_2d<Scalar>("kessler_r",            ncols, nz);
      rhalf        = view_2d<Scalar>("kessler_rhalf",        ncols, nz);
      velqr        = view_2d<Scalar>("kessler_velqr",        ncols, nz);
      sed          = view_2d<Scalar>("kessler_sed",          ncols, nz);
      pc           = view_2d<Scalar>("kessler_pc",            ncols, nz);
      dt0          = view_1d<Scalar>("kessler_dt0",           ncols);
      mask         = view_1d<Scalar>("kessler_mask",          ncols);
      time_counter = view_1d<Scalar>("kessler_time_counter",  ncols);
      precl_acc    = view_1d<Scalar>("kessler_precl_acc",     ncols);
    }
  };

  //
  // --------- Functions ---------
  //

  // Main Kessler warm rain microphysics.
  // All 2D arrays have layout (ncols, nz).
  // lyr_surf and lyr_toa are 0-based level indices; the sign of their
  // difference determines lyr_step (+1 or -1).
  // workspace must have been sized via KesslerWorkspace::init(ncols, nz)
  // (kessler_run also calls init(), which is a no-op if already sized).
  static void kessler_run(
    const int ncols,
    const int nz,
    const Scalar dt,
    const int lyr_surf,
    const int lyr_toa,
    const KesslerData& kd,
    KesslerWorkspace& workspace,
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
  // cpair is a 1D (composition-independent) scalar constant.
  static void kessler_update_timestep_final(
    const int ncols,
    const int nz,
    const Scalar gravit,
    const Scalar cpair,
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
