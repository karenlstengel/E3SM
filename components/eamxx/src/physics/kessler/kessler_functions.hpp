#ifndef KESSLER_FUNCTIONS_HPP
#define KESSLER_FUNCTIONS_HPP

#include "share/physics/physics_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/core/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_workspace.hpp>

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

  template <typename S>
  using BigPack = ekat::Pack<S,SCREAM_PACK_SIZE>;
  template <typename S>
  using SmallPack = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;

  using Pack = BigPack<Scalar>;
  using Spack = SmallPack<Scalar>;

  using KT = KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;

  template <typename S> using view_1d   = typename KT::template view_1d<S>;
  template <typename S> using view_2d   = typename KT::template view_2d<S>;
  template <typename S> using view_2dl  = typename KT::template lview<S**>;
  template <typename S> using uview_1d  = typename ekat::template Unmanaged<view_1d<S> >;
  template <typename S> using uview_2d  = typename ekat::template Unmanaged<view_2d<S> >;
  template <typename S> using uview_2dl = typename ekat::template Unmanaged<view_2dl<S> >;
  // ----------------------------------------
  // Structs
  struct params {
    // From field manager
    view_2d<Spack>  p_mid;
    view_2d<Spack>  T_mid;

    // For fortran
    view_2dl<Real>  f_p_mid;
    view_2dl<Real>  f_T_mid;

    // Set number of variables for ATMBufferManager
    static constexpr int num_2d_midlv_c_views = 2;
    static constexpr int num_2d_midlv_f_views = 2;

    // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
    void init(int ncol_in, int pver_in) {
      Real init_fill_value = -999;

      for (int i=0; i<ncol_in; ++i) {
        for (int j=0; j<pver_in; ++j) {
          f_p_mid(i,j) = init_fill_value;
          f_T_mid(i,j) = init_fill_value;
        }
      }
    }; // End init

    // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
    template <ekat::TransposeDirection::Enum D>
    void transpose(int ncol_in, int pver_in) {
      auto pverp = pver_in+1;
      if (D == ekat::TransposeDirection::c2f) {
        for (int i=0; i<ncol_in; ++i) {
          for (int j=0; j<pver_in; ++j) {
            f_p_mid(i,j) = p_mid(i,j/Spack::n)[j%Spack::n];
            f_T_mid(i,j) = T_mid(i,j/Spack::n)[j%Spack::n];
          }
        }
      }
      if (D == ekat::TransposeDirection::f2c) {
        for (int i=0; i<ncol_in; ++i) {
          // mid-point level variables
          for (int j=0; j<pver_in; ++j) {
            p_mid(i,j/Spack::n)[j%Spack::n] = f_p_mid(i,j);
            T_mid(i,j/Spack::n)[j%Spack::n] = f_T_mid(i,j);
          }
        }
      }
    }; // End transpose

  }; // end Struct params

}; // struct KesslerMicrophysicsFunctions

} // namespace kessler
} // namespace scream

#endif // KESSLER_FUNCTIONS_HPP