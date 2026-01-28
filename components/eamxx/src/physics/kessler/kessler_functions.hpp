#ifndef KESSLER_FUNCTIONS_HPP
#define KESSLER_FUNCTIONS_HPP

#include "share/physics/physics_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/core/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_team_policy_utils.hpp>
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

  using KT         = ekat::KokkosTypes<Device>;
  // using MemberType = typename KT::MemberType;
  using TeamPolicy = typename KokkosTypes<Device>::TeamPolicy;

  template <typename S> using view_1d   = typename KT::template view_1d<S>;
  template <typename S> using view_2d   = typename KT::template view_2d<S>;
  template <typename S> using view_2dl  = typename KT::template lview<S**>;
  template <typename S> using uview_1d  = typename ekat::template Unmanaged<view_1d<S> >;
  template <typename S> using uview_2d  = typename ekat::template Unmanaged<view_2d<S> >;
  template <typename S> using uview_2dl = typename ekat::template Unmanaged<view_2dl<S> >;

  // ----------------------------------------
  // Structs
  struct params_in {
    params_in() = default;
    // Needed inputs to kessler microphysics

    // integer,          intent(in)    :: lyr_surf   ! Index of surface layer in the vertical coordinate
    // integer,          intent(in)    :: lyr_toa    ! Index of top of the atmosphere in the vertical coordinate
    // real(kind_phys),  intent(in)    :: cpair(:,:) ! Specific_heat_of_dry_air_at_constant_pressure (J/kg/K)
    // real(kind_phys),  intent(in)    :: rair(:,:)  ! Gas constant of dry air (J/kg/K)
    // real(kind_phys),  intent(in)    :: rho(:,:)   ! Dry air density (kg/m^3)
    // real(kind_phys),  intent(in)    :: z(:,:)     ! Heights of thermo. levels (m)
    // real(kind_phys),  intent(in)    :: pk(:,:)    ! Exner function (p/p0)**(R/cp)

    // view_2d<Spack>  cpair;
    // view_2d<Spack>  rair;
    view_2d<Spack>  rho;
    view_2d<Spack>  dz;
    view_2d<Spack>  pk;

    // Fortran holders/in Fortran format
    uview_2dl<Real>  f_cpair;
    uview_2dl<Real>  f_rair;
    uview_2dl<Real>  f_rho;
    uview_2dl<Real>  f_dz;
    uview_2dl<Real>  f_pk;
    // Set number of variables for ATMBufferManager
    static constexpr int num_1d_intgr = 0;  // number of 1D integer views
    static constexpr int num_1d_scalr = 0;  // number of 1D scalar views
    static constexpr int num_2d_c     = 3;  // number of 2D field views for C++
    static constexpr int num_2d_f     = 5;  // number of 2D field views for fortran

    // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
    void init(int ncol_in, int pver_in) { // TODO - Kokko-ize this
      using PC  = scream::physics::Constants<Real>;

      const Real cpair  = PC::Cpair; // Specific heat of dry air at constant pressure
      const Real Rair   = PC::Rair;  // Gas constant of dry air

      Real init_fill_value = 0;

      using MDPolicy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
      MDPolicy mdp({0,0}, {ncol_in, pver_in});

      Kokkos::parallel_for("init_params_in", mdp,
        KOKKOS_CLASS_LAMBDA(const int i, const int j) {
          f_cpair(i,j) = cpair;
          f_rair(i,j) = Rair;
          f_rho(i,j) = init_fill_value;
          f_dz(i,j) = init_fill_value;
          f_pk(i,j) = init_fill_value;
          }
        );
        
    }; // End init

    // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
    template <ekat::TransposeDirection::Enum D>
    void transpose(int ncol_in, int pver_in) { // TODO - Kokko-ize this
      auto pver_in_packs = ekat::npack<Spack>(pver_in);

      // using MDPolicy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
      // MDPolicy mdp({0,0}, {ncol_in, pver_in});

      if (D == ekat::TransposeDirection::c2f) {

        Kokkos::parallel_for(
            "transpose c2f", KT::RangePolicy(0, ncol_in * pver_in_packs),
            KOKKOS_CLASS_LAMBDA(const int i) {
              const int icol            = i / pver_in_packs;
              const int klev            = i % pver_in_packs;
              f_rho(icol, klev) = rho(icol, klev / Spack::n)[klev % Spack::n];
              f_dz(icol, klev) = dz(icol, klev / Spack::n)[klev % Spack::n];
              f_pk(icol, klev) = pk(icol, klev / Spack::n)[klev % Spack::n];
            }
          );

        // for (int i=0; i<ncol_in; ++i) {
        //   for (int j=0; j<pver_in; ++j) {
        //     // f_cpair(i,j) = cpair(i,j/Spack::n)[j%Spack::n];
        //     // f_rair(i,j) = rair(i,j/Spack::n)[j%Spack::n];
        //     f_rho(i,j) = rho(i,j/Spack::n)[j%Spack::n];
        //     f_dz(i,j) = dz(i,j/Spack::n)[j%Spack::n];
        //     f_pk(i,j) = pk(i,j/Spack::n)[j%Spack::n];
        //   }
        // }
      }
      if (D == ekat::TransposeDirection::f2c) {  // Not needed but leaving in just in case/temporary

        Kokkos::parallel_for(
            "transpose f2c", KT::RangePolicy(0, ncol_in * pver_in_packs),
            KOKKOS_CLASS_LAMBDA(const int i) {
              const int icol                                      = i / pver_in_packs;
              const int klev                                      = i % pver_in_packs;
              rho(icol, klev / Spack::n)[klev % Spack::n] = f_rho(icol, klev);
              dz(icol, klev / Spack::n)[klev % Spack::n] = f_dz(icol, klev);
              pk(icol, klev / Spack::n)[klev % Spack::n] = f_pk(icol, klev);
            }
          );

        // for (int i=0; i<ncol_in; ++i) {
        //   // mid-point level variables
        //   for (int j=0; j<pver_in; ++j) {
        //     // cpair(i,j/Spack::n)[j%Spack::n] = f_cpair(i,j);
        //     // rair(i,j/Spack::n)[j%Spack::n] = f_rair(i,j);
        //     rho(i,j/Spack::n)[j%Spack::n] = f_rho(i,j);
        //     dz(i,j/Spack::n)[j%Spack::n] = f_dz(i,j);
        //     pk(i,j/Spack::n)[j%Spack::n] = f_pk(i,j);
        //   }
        // }
      }
    }; // End transpose

  }; // end Struct params_in

  struct params_out {

    params_out() = default;
    // Needed outputs to kessler microphysics

    // real(kind_phys),  intent(inout) :: theta(:,:) ! Potential temperature (K)
    // real(kind_phys),  intent(inout) :: qv(:,:)    ! Water vapor mixing ratio wrt dry air (kg/kg)
    // real(kind_phys),  intent(inout) :: qc(:,:)    ! Cloud water mixing ratio wrt dry air (kg/kg)
    // real(kind_phys),  intent(inout) :: qr(:,:)    ! Rain water mixing ratio wrt dry air (kg/kg)
    // real(kind_phys),  intent(out)   :: precl(:)   ! Precipitation rate (m_water / s)
    // real(kind_phys),  intent(out)   :: relhum(:,:)! Relative humidity in percent

    view_2d<Spack>  theta;
    view_2d<Spack>  qv;
    view_2d<Spack>  qc;
    view_2d<Spack>  qr;
    view_1d<Scalar> precl;
    view_2d<Spack>  relhum;

    // For fortran (left layout) versions
    uview_2dl<Real>  f_theta;
    uview_2dl<Real>  f_qv;
    uview_2dl<Real>  f_qc;
    uview_2dl<Real>  f_qr;
    uview_1d<Real>   f_precl;  
    uview_2dl<Real>  f_relhum;

    // Set number of variables for ATMBufferManager
    static constexpr int num_1d_intgr = 0;  // number of 1D integer views
    static constexpr int num_1d_scalr = 1;  // number of 1D scalar views
    static constexpr int num_2d_c     = 5;  // number of 2D fields
    static constexpr int num_2d_f     = 5;  // number of 2D fields

    // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
    void init(int ncol_in, int pver_in) {
      Real init_fill_value = 0;

      using MDPolicy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
      MDPolicy mdp({0,0}, {ncol_in, pver_in});
      
      Kokkos::parallel_for("init_params_out", mdp,
        KOKKOS_CLASS_LAMBDA(const int i, const int j) {
          f_theta(i,j) = init_fill_value;
          f_qv(i,j) = init_fill_value;
          f_qc(i,j) = init_fill_value;
          f_qr(i,j) = init_fill_value;
          f_relhum(i,j) = init_fill_value;
          f_precl(i) = init_fill_value;
          }
        );

      // for (int i=0; i<ncol_in; ++i) {
      //   for (int j=0; j<pver_in; ++j) {
      //     f_theta(i,j) = init_fill_value;
      //     f_qv(i,j) = init_fill_value;
      //     f_qc(i,j) = init_fill_value;
      //     f_qr(i,j) = init_fill_value;
      //     f_relhum(i,j) = init_fill_value;
      //   }
      //   f_precl(i) = init_fill_value;
      // }
    }; // End init

    // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
    template <ekat::TransposeDirection::Enum D>
    void transpose(int ncol_in, int pver_in) { // TODO - Kokko-ize this
      auto pver_in_packs = ekat::npack<Spack>(pver_in);

      if (D == ekat::TransposeDirection::c2f) {

        Kokkos::parallel_for(
            "transpose c2f", KT::RangePolicy(0, ncol_in * pver_in_packs),
            KOKKOS_CLASS_LAMBDA(const int i) {
              const int icol            = i / pver_in_packs;
              const int klev            = i % pver_in_packs;
              f_theta(icol, klev) = theta(icol, klev / Spack::n)[klev % Spack::n];
              f_qv(icol, klev) = qv(icol, klev / Spack::n)[klev % Spack::n];
              f_qc(icol, klev) = qc(icol, klev / Spack::n)[klev % Spack::n];
              f_qr(icol, klev) = qr(icol, klev / Spack::n)[klev % Spack::n];
              f_relhum(icol, klev) = relhum(icol, klev / Spack::n)[klev % Spack::n];
              f_precl(icol) = precl(icol);
            }
          );
        // for (int i=0; i<ncol_in; ++i) {
        //   for (int j=0; j<pver_in; ++j) {
        //     f_theta(i,j) = theta(i,j/Spack::n)[j%Spack::n];
        //     f_qv(i,j) = qv(i,j/Spack::n)[j%Spack::n];
        //     f_qc(i,j) = qc(i,j/Spack::n)[j%Spack::n];
        //     f_qr(i,j) = qr(i,j/Spack::n)[j%Spack::n];
        //     f_relhum(i,j) = relhum(i,j/Spack::n)[j%Spack::n];
        //   }
        //   f_precl(i) = precl(i);
        // }
      }
      if (D == ekat::TransposeDirection::f2c) {

        Kokkos::parallel_for(
            "transpose f2c", KT::RangePolicy(0, ncol_in * pver_in_packs),
            KOKKOS_CLASS_LAMBDA(const int i) {
              const int icol                                      = i / pver_in_packs;
              const int klev                                      = i % pver_in_packs;
              theta(icol, klev / Spack::n)[klev % Spack::n] = f_theta(icol, klev);
              qv(icol, klev / Spack::n)[klev % Spack::n] = f_qv(icol, klev);
              qc(icol, klev / Spack::n)[klev % Spack::n] = f_qc(icol, klev);
              qr(icol, klev / Spack::n)[klev % Spack::n] = f_qr(icol, klev);
              relhum(icol, klev / Spack::n)[klev % Spack::n] = f_relhum(icol, klev);
            }
          );
        // for (int i=0; i<ncol_in; ++i) {
        //   // mid-point level variables
        //   for (int j=0; j<pver_in; ++j) {
        //     theta(i,j/Spack::n)[j%Spack::n] = f_theta(i,j);
        //     qv(i,j/Spack::n)[j%Spack::n] = f_qv(i,j);
        //     qc(i,j/Spack::n)[j%Spack::n] = f_qc(i,j);
        //     qr(i,j/Spack::n)[j%Spack::n] = f_qr(i,j);
        //     relhum(i,j/Spack::n)[j%Spack::n] = f_relhum(i,j);
        //   }
        //   precl(i) = f_precl(i);
        // }
      }
    }; // End transpose

  }; // end Struct params_out

}; // struct KesslerMicrophysicsFunctions
  
} // namespace kessler
} // namespace scream

#endif // KESSLER_FUNCTIONS_HPP