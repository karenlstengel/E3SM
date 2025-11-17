#ifndef STENGELF_FUNCTIONS_HPP
#define STENGELF_FUNCTIONS_HPP

#include "share/core/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_workspace.hpp>

namespace scream {
namespace stengelF {

template <typename ScalarT, typename DeviceT>
struct StengelFFunctions
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
  template <typename S> 
  using uview_2dl = typename ekat::template Unmanaged<view_2dl<S> >;

  // ----------------------------------------
  // Structs
  struct params {
    // From field manager
    view_2d<const Spack>  p_mid;
    view_2d<      Spack>  T_mid;

    // For fortran
    uview_2dl<Real>  f_p_mid;
    uview_2dl<Real>  f_T_mid;
  } // end Struct params

}; // struct StengelFFunctions

} // namespace stengelF
} // namespace scream

#endif // STENGELF_FUNCTIONS_HPP
