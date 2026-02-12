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
  // template <typename S> using uview_1d  = typename ekat::template Unmanaged<view_1d<S> >;
  // template <typename S> using uview_2d  = typename ekat::template Unmanaged<view_2d<S> >;
  // template <typename S> using uview_2dl = typename ekat::template Unmanaged<view_2dl<S> >;

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
    view_2dl<Real>  f_cpair;
    view_2dl<Real>  f_rair;
    view_2dl<Real>  f_rho;
    view_2dl<Real>  f_dz;
    view_2dl<Real>  f_pk;
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
    void transpose(int ncol_in, int pver_in) { 
      // auto pver_in_packs = ekat::npack<Spack>(pver_in);

      if (D == ekat::TransposeDirection::c2f) {

        Kokkos::parallel_for(
            "transpose c2f", KT::RangePolicy(0, ncol_in * pver_in),
            KOKKOS_CLASS_LAMBDA(const int i) {
              const int icol            = i / pver_in;
              const int klev            = i % pver_in;
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
            "transpose f2c", KT::RangePolicy(0, ncol_in * pver_in),
            KOKKOS_CLASS_LAMBDA(const int i) {
              const int icol = i / pver_in;
              const int klev = i % pver_in;
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
    view_2dl<Real>  f_theta;
    view_2dl<Real>  f_qv;
    view_2dl<Real>  f_qc;
    view_2dl<Real>  f_qr;
    view_1d<Real>   f_precl;  
    view_2dl<Real>  f_relhum;

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
    }; // End init

    // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
    template <ekat::TransposeDirection::Enum D>
    void transpose(int ncol_in, int pver_in) { // TODO - Kokko-ize this
      // auto pver_in_packs = ekat::npack<Spack>(pver_in);

      if (D == ekat::TransposeDirection::c2f) {

        Kokkos::parallel_for(
            "transpose c2f", KT::RangePolicy(0, ncol_in * pver_in),
            KOKKOS_CLASS_LAMBDA(const int i) {
              const int icol = i / pver_in;
              const int klev = i % pver_in;
              f_theta(icol, klev) = theta(icol, klev / Spack::n)[klev % Spack::n];
              f_qv(icol, klev) = qv(icol, klev / Spack::n)[klev % Spack::n];
              f_qc(icol, klev) = qc(icol, klev / Spack::n)[klev % Spack::n];
              f_qr(icol, klev) = qr(icol, klev / Spack::n)[klev % Spack::n];
              f_relhum(icol, klev) = relhum(icol, klev / Spack::n)[klev % Spack::n];
              f_precl(icol) = precl(icol);
            }
          );
      }
      if (D == ekat::TransposeDirection::f2c) {

        Kokkos::parallel_for(
            "transpose f2c", KT::RangePolicy(0, ncol_in * pver_in),
            KOKKOS_CLASS_LAMBDA(const int i) {
              const int icol = i / pver_in;
              const int klev  = i % pver_in;
              theta(icol, klev / Spack::n)[klev % Spack::n] = f_theta(icol, klev);
              qv(icol, klev / Spack::n)[klev % Spack::n] = f_qv(icol, klev);
              qc(icol, klev / Spack::n)[klev % Spack::n] = f_qc(icol, klev);
              qr(icol, klev / Spack::n)[klev % Spack::n] = f_qr(icol, klev);
              relhum(icol, klev / Spack::n)[klev % Spack::n] = f_relhum(icol, klev);
              precl(icol) = f_precl(icol);
            }
          );
      }
    }; // End transpose

  }; // end Struct params_out

  struct params_update {

    params_update() = default;
    // Needed variables to update after kessler microphysics

    // real(kind_phys),  intent(inout) :: temp_prev(:,:)  ! Previous temperature (K)
    // real(kind_phys),  intent(inout) :: temp(:,:)       ! Current temperature (K)
    // real(kind_phys),  intent(inout) :: temp_tend(:,:)  ! Temperature tendency (K/s)
    // real(kind_phys),  intent(inout) :: zm(:,:)         ! Mass of dry air in layer (kg/m^2)
    // real(kind_phys),  intent(inout) :: phis(:)         ! Surface geopotential (m^2/s^2)
    // real(kind_phys),  intent(inout) :: st_energy(:,:)  ! Surface energy (J/m^2)

    view_2d<Spack>   temp_prev;
    view_2d<Spack>   temp;
    view_2d<Spack>   temp_tend;
    view_2d<Spack>   z_mid;
    view_2d<Spack>   z_int; // Helper, doesn't need to get copied to fortran
    view_1d<Scalar>  phis;
    view_2d<Spack>   st_energy;

    view_2dl<Real>  f_temp_prev;
    view_2dl<Real>  f_temp;
    view_2dl<Real>  f_temp_tend;
    view_2dl<Real>  f_z_mid;
    view_1d<Real>   f_phis;
    view_2dl<Real>  f_st_energy;

    // Set number of variables for ATMBufferManager
    static constexpr int num_1d_intgr = 0;  // number of 1D integer views
    static constexpr int num_1d_scalr = 1;  // number of 1D scalar views
    static constexpr int num_2d_c     = 6;  // number of 2D fields
    static constexpr int num_2d_f     = 5;  // number of 2D fields

    // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
    void init(int ncol_in, int pver_in) {
      Real init_fill_value = 0;

      using MDPolicy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
      MDPolicy mdp({0,0}, {ncol_in, pver_in});
      
      Kokkos::parallel_for("init_params_update", mdp,
        KOKKOS_CLASS_LAMBDA(const int i, const int j) {
          f_temp_prev(i,j) = init_fill_value;
          f_temp(i,j) = init_fill_value;
          f_temp_tend(i,j) = init_fill_value;
          f_z_mid(i,j) = init_fill_value;
          f_st_energy(i,j) = init_fill_value;

          f_phis(i) = init_fill_value;
          }
        );

    }; // End init

    // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
    template <ekat::TransposeDirection::Enum D>
    void transpose(int ncol_in, int pver_in) { // TODO - Kokko-ize this
      // auto pver_in_packs = ekat::npack<Spack>(pver_in);

      if (D == ekat::TransposeDirection::c2f) {

        Kokkos::parallel_for(
            "transpose c2f", KT::RangePolicy(0, ncol_in * pver_in),
            KOKKOS_CLASS_LAMBDA(const int i) {
              const int icol = i / pver_in;
              const int klev = i % pver_in;

              f_temp_prev(icol, klev) = temp_prev(icol, klev / Spack::n)[klev % Spack::n];
              f_temp(icol, klev) = temp(icol, klev / Spack::n)[klev % Spack::n];
              f_temp_tend(icol, klev) = temp_tend(icol, klev / Spack::n)[klev % Spack::n];
              f_z_mid(icol, klev) = z_mid(icol, klev / Spack::n)[klev % Spack::n];
              f_st_energy(icol, klev) = st_energy(icol, klev / Spack::n)[klev % Spack::n];
              f_phis(icol) = phis(icol);
            }
          );
      }
      if (D == ekat::TransposeDirection::f2c) {

        Kokkos::parallel_for(
            "transpose f2c", KT::RangePolicy(0, ncol_in * pver_in),
            KOKKOS_CLASS_LAMBDA(const int i) {
              const int icol = i / pver_in;
              const int klev = i % pver_in;

              temp_prev(icol, klev / Spack::n)[klev % Spack::n] = f_temp_prev(icol, klev);
              temp(icol, klev / Spack::n)[klev % Spack::n] = f_temp(icol, klev);
              temp_tend(icol, klev / Spack::n)[klev % Spack::n] = f_temp_tend(icol, klev);
              z_mid(icol, klev / Spack::n)[klev % Spack::n] = f_z_mid(icol, klev); // could skip this
              st_energy(icol, klev / Spack::n)[klev % Spack::n] = f_st_energy(icol, klev);
              phis(icol) = f_phis(icol); // could skip this
            }
          );
      }
    }; // End transpose

  }; // end Struct params_update

}; // struct KesslerMicrophysicsFunctions
  
} // namespace kessler
} // namespace scream

#endif // KESSLER_FUNCTIONS_HPP