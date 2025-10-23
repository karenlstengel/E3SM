#ifndef STENGEL_FUNCTIONS_HPP
#define STENGEL_FUNCTIONS_HPP

#include "share/core/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_workspace.hpp>

namespace scream {
namespace stengel {

template <typename ScalarT, typename DeviceT>
struct StengelFunctions
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

  // static void main(
  //   const Int nj, 
  //   const Int nk,
  //   const Real ice_threshold,
  //   const Real ice_4out_threshold,
  //   const view_2d<const Pack>& qi, 
  //   const view_2d<const Pack>& liq_cld_frac, 
  //   const view_2d<Pack>& ice_cld_frac, 
  //   const view_2d<Pack>& tot_cld_frac,
  //   const view_2d<Pack>& ice_cld_frac_4out, 
  //   const view_2d<Pack>& tot_cld_frac_4out);

}; // struct Functions

} // namespace stengel
} // namespace scream

#endif // STENGEL_FUNCTIONS_HPP
