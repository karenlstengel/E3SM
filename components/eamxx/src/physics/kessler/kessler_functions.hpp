#ifndef KESSLER_FUNCTIONS_HPP
#define KESSLER_FUNCTIONS_HPP

#include "share/physics/physics_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/core/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_pack_math.hpp>
#include <ekat_team_policy_utils.hpp>
#include <ekat_workspace.hpp>

/*
 * KesslerMicrophysicsFunctions encapsulates the Kessler (1969) warm
 * rain microphysics parameterization. kessler_run additionally takes
 * a caller-owned Workspace of persistent scratch views (see Workspace
 * below), so the struct itself carries no state, but a given call's
 * behavior depends on the externally-owned scratch passed in.
 *
 * The Kessler scheme contains three moisture categories: water vapor,
 * cloud water (liquid water that moves with the flow), and rain water
 * (liquid water that falls relative to the surrounding air).
 *
 * References:
 *   Klemp & Wilhelmson (1978), J. Atmos. Sci., 35, 1070-1096
 *   Durran & Klemp (1983), Mon. Wea. Rev., 111, 2341-2361
 */

namespace scream {
namespace kessler {

template <typename ScalarT, typename DeviceT>
struct KesslerMicrophysicsFunctions
{

  //
  // ------- Types --------
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  using Pack    = ekat::Pack<Scalar,SCREAM_PACK_SIZE>;
  using IntPack = ekat::Pack<Int,SCREAM_PACK_SIZE>;

  using KT      = ekat::KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;
  using TeamPolicy = typename KokkosTypes<Device>::TeamPolicy;

  template <typename S> using view_1d   = typename KT::template view_1d<S>;
  template <typename S> using view_2d   = typename KT::template view_2d<S>;
  template <typename S> using view_2dl  = typename KT::template lview<S**>;

  //
  // --------- Workspace ---------
  //

  // Persistent scratch used internally by kessler_run, owned by the
  // caller so it can be allocated once (e.g. carved from one
  // ATMBufferManager allocation by KesslerMicrophysics::init_buffers,
  // mirroring the P3/SHOC Buffer pattern) instead of allocated and
  // freed on every call. The 2D fields are per-(col,level), packed
  // like the rest of the scheme's field views; the 1D fields are one
  // value per column (sub-cycle bookkeeping), so they stay unpacked
  // -- matching how P3's own Buffer keeps its per-column fields
  // (e.g. precip_liq_surf_flux) as plain Real while its per-level
  // fields are Pack.
  struct Workspace {
    static constexpr int num_2d_vector = 6;
    static constexpr int num_1d_scalar = 4;

    // Per-(col,level) scratch, sized (ncols, nlev_packs).
    view_2d<Pack> r, rhalf, velqr, sed, pc, f5;
    // Per-column scratch, sized (ncols).
    view_1d<Scalar> dt0, mask, time_counter, precl_acc;
  };

  //
  // --------- Functions ---------
  //

  // Main Kessler warm rain microphysics.
  // All 2D arrays have layout (ncols, nz), packed along the level dimension.
  // lyr_surf and lyr_toa are 0-based level indices; the sign of their
  // difference determines lyr_step (+1 or -1).
  // cpair and rair may vary by column and level (composition-dependent).
  // z_mid is level height at layer midpoints; sedimentation uses the
  // height difference between adjacent levels (note: the reference
  // Fortran's "dz" argument to kessler_run is this same z_mid, not a
  // layer thickness -- that naming carried over confusingly).
  // ws is caller-owned persistent scratch (see Workspace above).
  static void kessler_run(
    const Int ncols,
    const Int nz,
    const Real dt,
    const Int lyr_surf,
    const Int lyr_toa,
    const Real rhoqr,
    const Real lv,
    const Real pref,
    const Real cpair,
    const Real rair,
    const view_2d<const Pack>& rho,
    const view_2d<const Pack>& z_mid,
    const view_2d<const Pack>& pk,
    const view_2d<Pack>& theta,
    const view_2d<Pack>& qv,
    const view_2d<Pack>& qc,
    const view_2d<Pack>& qr,
    const view_1d<Scalar>& precl,
    const view_2d<Pack>& relhum,
    const Workspace& ws);

  // Save the current temperature and zero the temperature tendency accumulator.
  // Called once per physics time step before kessler_run.
  static void kessler_update_timestep_init(
    const Int ncols,
    const Int nz,
    const view_2d<const Pack>& temp,
    const view_2d<Pack>& temp_prev,
    const view_2d<Pack>& ttend_t);

  // Back out the temperature tendency due to Kessler from theta and exner.
  // ttend_t += (theta * exner - temp_prev) / dt
  static void kessler_update_run(
    const Int ncols,
    const Int nz,
    const Real dt,
    const view_2d<const Pack>& theta,
    const view_2d<const Pack>& exner,
    const view_2d<const Pack>& temp_prev,
    const view_2d<Pack>& ttend_t);

  // Compute dry static energy: st_energy = cpair*temp + gravit*zm + phis
  // phis (surface geopotential) is one value per column, not per-level,
  // so it stays unpacked even though the rest of this function's fields
  // are packed along the level dimension.
  static void kessler_update_timestep_final(
    const Int ncols,
    const Int nz,
    const Real gravit,
    const Real cpair,
    const view_2d<const Pack>& temp,
    const view_2d<const Pack>& z_mid,
    const view_1d<const Scalar>& phis,
    const view_2d<Pack>& st_energy);

}; // struct KesslerMicrophysicsFunctions

} // namespace kessler
} // namespace scream

// If a GPU build without relocatable device code, include the full
// implementation in every translation unit; otherwise ETI is used.
#if defined(EAMXX_ENABLE_GPU) && !defined(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE) \
                                && !defined(KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
# include "kessler_functions_impl.hpp"
#endif

#endif // KESSLER_FUNCTIONS_HPP
