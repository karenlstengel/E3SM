#ifndef STENGELF_FUNCTIONS_HPP
#define STENGELF_FUNCTIONS_HPP

#include "share/core/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_team_policy_utils.hpp>
#include <ekat_workspace.hpp>

namespace scream
{
namespace stengelF
{

template <typename ScalarT, typename DeviceT> struct StengelFFunctions {

  //
  // ------- Types --------
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  template <typename S> using BigPack   = ekat::Pack<S, SCREAM_PACK_SIZE>;
  template <typename S> using SmallPack = ekat::Pack<S, SCREAM_SMALL_PACK_SIZE>;

  using Pack  = BigPack<Scalar>;
  using Spack = SmallPack<Scalar>;

  using Mask  = ekat::Mask<BigPack<Scalar>::n>;
  using Smask = ekat::Mask<SmallPack<Scalar>::n>;

  // GPU/Kokkos related things
  using KT         = ekat::KokkosTypes<Device>;
  using TeamPolicy = typename KokkosTypes<Device>::TeamPolicy;

  template <typename S> using view_1d  = typename KT::template view_1d<S>;
  template <typename S> using view_2d  = typename KT::template view_2d<S>;
  template <typename S> using view_2dl = typename KT::template lview<S **>;

  #if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
    template <typename S> using fview_1d  = typename ekat::template Unmanaged<view_1d<S> >;
    template <typename S> using fview_2d  = typename ekat::template Unmanaged<view_2d<S> >;
    template <typename S> using fview_2dl = typename ekat::template Unmanaged<view_2dl<S> >;
    template <typename S> using view_1dh  = typename view_1d<S>::HostMirror;
    template <typename S> using view_2dh  = typename view_2dl<S>::HostMirror;
  #else 
    template <typename S> using fview_1d  = typename KT::template view_1d<S>;
    template <typename S> using fview_2d  = typename KT::template view_2d<S>;
    template <typename S> using fview_2dl = typename KT::template lview<S **>;
  #endif
  // ----------------------------------------
  // Structs
  struct params {
    params() = default;

    // From field manager
    view_2d<Spack> p_mid;
    view_2d<Spack> T_mid;

    // For fortran
    fview_2dl<Real> f_p_mid;
    fview_2dl<Real> f_T_mid;
    
    #if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
      // Host mirrors
      view_2dh<Real> h_p_mid;
      view_2dh<Real> h_T_mid;
    #endif
    

    // Set number of variables for ATMBufferManager
    static constexpr int num_2d_midlv_c_views = 2;
    static constexpr int num_2d_midlv_f_views = 2;

    // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
    void
    init(int ncol_in, int pver_in)
    {
      Real init_fill_value = -999;

      Kokkos::parallel_for(
          Kokkos::MDRangePolicy<typename KT::ExeSpace, Kokkos::Rank<2>>({0, 0}, {ncol_in, pver_in}),
          KOKKOS_CLASS_LAMBDA(const int i, const int j) {
            f_p_mid(i, j) = init_fill_value;
            f_T_mid(i, j) = init_fill_value;
          });
    }; // End init

    // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
    template <ekat::TransposeDirection::Enum D>
    void
    transpose(int ncol_in, int pver_in)
    {
      auto pver_in_packs = ekat::npack<Spack>(pver_in);
      if (D == ekat::TransposeDirection::c2f) {

        Kokkos::parallel_for(
            "transpose c2f", KT::RangePolicy(0, ncol_in * pver_in_packs),
            KOKKOS_CLASS_LAMBDA(const int i) {
              const int icol            = i / pver_in_packs;
              const int klev            = i % pver_in_packs;
              f_p_mid(icol, klev) = p_mid(icol, klev / Spack::n)[klev % Spack::n];
              f_T_mid(icol, klev) = T_mid(icol, klev / Spack::n)[klev % Spack::n];
            }
          );
          #if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
            // Copy from device to host mirrors for Fortran
            Kokkos::deep_copy(h_p_mid, f_p_mid);
            Kokkos::deep_copy(h_T_mid, f_T_mid);
          #endif
      }
      if (D == ekat::TransposeDirection::f2c) {
        #if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
          // Copy from host mirrors back to device
          Kokkos::deep_copy(f_p_mid,   h_p_mid);
          Kokkos::deep_copy(f_T_mid,   h_T_mid);
        #endif

        Kokkos::parallel_for(
            "transpose f2c", KT::RangePolicy(0, ncol_in * pver_in_packs),
            KOKKOS_CLASS_LAMBDA(const int i) {
              const int icol                                      = i / pver_in_packs;
              const int klev                                      = i % pver_in_packs;
              p_mid(icol, klev / Spack::n)[klev % Spack::n] = f_p_mid(icol, klev);
              T_mid(icol, klev / Spack::n)[klev % Spack::n] = f_T_mid(icol, klev);
            }
          );
      }
    }; // End transpose

  }; // end Struct params

}; // struct StengelFFunctions

} // namespace stengelF
} // namespace scream

#endif // STENGELF_FUNCTIONS_HPP
