#ifndef KESSLER_FUNCTIONS_HPP
#define KESSLER_FUNCTIONS_HPP

#include "share/physics/physics_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/physics/eamxx_common_physics_functions_impls.hpp"
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
  struct params_in {
    // Needed inputs to kessler microphysics

    // integer,          intent(in)    :: lyr_surf   ! Index of surface layer in the vertical coordinate
    // integer,          intent(in)    :: lyr_toa    ! Index of top of the atmosphere in the vertical coordinate
    // real(kind_phys),  intent(in)    :: cpair(:,:) ! Specific_heat_of_dry_air_at_constant_pressure (J/kg/K)
    // real(kind_phys),  intent(in)    :: rair(:,:)  ! Gas constant of dry air (J/kg/K)
    // real(kind_phys),  intent(in)    :: rho(:,:)   ! Dry air density (kg/m^3)
    // real(kind_phys),  intent(in)    :: z(:,:)     ! Heights of thermo. levels (m)
    // real(kind_phys),  intent(in)    :: pk(:,:)    ! Exner function (p/p0)**(R/cp)

    view_2d<Spack>  cpair;
    view_2d<Spack>  rair;
    view_2d<Spack>  rho;
    view_2d<Spack>  z;
    view_2d<Spack>  pk;

    // Fortran holders/in Fortran format
    uview_2dl<Real>  f_cpair;
    uview_2dl<Real>  f_rair;
    uview_2dl<Real>  f_rho;
    uview_2dl<Real>  f_z;
    uview_2dl<Real>  f_pk;
    // Set number of variables for ATMBufferManager
    static constexpr int num_1d_intgr = 0;  // number of 1D integer views
    static constexpr int num_1d_scalr = 0;  // number of 1D scalar views
    static constexpr int num_2d       = 5;  // number of 2D field views

    // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
    void init(int ncol_in, int pver_in) {

      const Real cpair  = PC::Cpair; // Specific heat of dry air at constant pressure
      const Real Rair   = PC::Rair;  // Gas constant of dry air

      Real init_fill_value = -999;

      for (int i=0; i<ncol_in; ++i) {
        for (int j=0; j<pver_in; ++j) {
          f_cpair(i,j) = cpair;
          f_rair(i,j) = Rair;
          f_rho(i,j) = init_fill_value;
          f_z(i,j) = init_fill_value;
          f_pk(i,j) = init_fill_value;
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
            f_cpair(i,j) = cpair(i,j/Spack::n)[j%Spack::n];
            f_rair(i,j) = rair(i,j/Spack::n)[j%Spack::n];
            f_rho(i,j) = rho(i,j/Spack::n)[j%Spack::n];
            f_z(i,j) = z(i,j/Spack::n)[j%Spack::n];
            f_pk(i,j) = pk(i,j/Spack::n)[j%Spack::n];
          }
        }
      }
      if (D == ekat::TransposeDirection::f2c) { // Not needed but leaving in just in case/temporary
        for (int i=0; i<ncol_in; ++i) {
          // mid-point level variables
          for (int j=0; j<pver_in; ++j) {
            cpair(i,j/Spack::n)[j%Spack::n] = f_cpair(i,j);
            rair(i,j/Spack::n)[j%Spack::n] = f_rair(i,j);
            rho(i,j/Spack::n)[j%Spack::n] = f_rho(i,j);
            z(i,j/Spack::n)[j%Spack::n] = f_z(i,j);
            pk(i,j/Spack::n)[j%Spack::n] = f_pk(i,j);
          }
        }
      }
    }; // End transpose

  }; // end Struct params_in

  struct params_out {
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
    uview_2dl<Real>  f_relhum;

    // Set number of variables for ATMBufferManager
    static constexpr int num_1d_intgr = 0;  // number of 1D integer views
    static constexpr int num_1d_scalr = 1;  // number of 1D scalar views
    static constexpr int num_2d       = 5;  // number of 2D fields

    // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
    void init(int ncol_in, int pver_in) {
      Real init_fill_value = -999;

      for (int i=0; i<ncol_in; ++i) {
        for (int j=0; j<pver_in; ++j) {
          f_theta(i,j) = init_fill_value;
          f_qv(i,j) = init_fill_value;
          f_qc(i,j) = init_fill_value;
          f_qr(i,j) = init_fill_value;
          f_relhum(i,j) = init_fill_value;
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
            f_theta(i,j) = theta(i,j/Spack::n)[j%Spack::n];
            f_qv(i,j) = qv(i,j/Spack::n)[j%Spack::n];
            f_qc(i,j) = qc(i,j/Spack::n)[j%Spack::n];
            f_qr(i,j) = qr(i,j/Spack::n)[j%Spack::n];
            f_relhum(i,j) = relhum(i,j/Spack::n)[j%Spack::n];
          }
        }
      }
      if (D == ekat::TransposeDirection::f2c) {
        for (int i=0; i<ncol_in; ++i) {
          // mid-point level variables
          for (int j=0; j<pver_in; ++j) {
            theta(i,j/Spack::n)[j%Spack::n] = f_theta(i,j);
            qv(i,j/Spack::n)[j%Spack::n] = f_qv(i,j);
            qc(i,j/Spack::n)[j%Spack::n] = f_qc(i,j);
            qr(i,j/Spack::n)[j%Spack::n] = f_qr(i,j);
            relhum(i,j/Spack::n)[j%Spack::n] = f_relhum(i,j);
          }
        }
      }
    }; // End transpose

  }; // end Struct params_out

}; // struct KesslerMicrophysicsFunctions

} // namespace kessler
} // namespace scream

#endif // KESSLER_FUNCTIONS_HPP