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

  using Pack    = ekat::Pack<Scalar,SCREAM_PACK_SIZE>;
  // using IntPack = ekat::Pack<Int,SCREAM_PACK_SIZE>;

  using KT      = ekat::KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;
  using TeamPolicy = typename KokkosTypes<Device>::TeamPolicy;

  template <typename S> using view_1d   = typename KT::template view_1d<S>;
  template <typename S> using view_2d   = typename KT::template view_2d<S>;
  template <typename S> using view_2dl  = typename KT::template lview<S**>;
  
//   #if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
//     // Needed for running the Fortran without OpenACC
//     template <typename S> using fview_1d  = typename ekat::template Unmanaged<view_1d<S> >;
//     template <typename S> using fview_2d  = typename ekat::template Unmanaged<view_2d<S> >;
//     template <typename S> using fview_2dl = typename ekat::template Unmanaged<view_2dl<S> >;
//     template <typename S> using view_1dh  = typename view_1d<S>::HostMirror;
//     template <typename S> using view_2dh  = typename view_2dl<S>::HostMirror;
//   #else 
//     template <typename S> using fview_1d  = typename KT::template view_1d<S>;
//     template <typename S> using fview_2d  = typename KT::template view_2d<S>;
//     template <typename S> using fview_2dl = typename KT::template lview<S **>;
//   #endif

//   // ----------------------------------------
//   // Structs
// struct params_helpers {
//     params_helpers() = default;

//     // Needed for kessler_run or to compute things to pass to kessler
//     // real(kind_phys),  intent(in)    :: cpair(:,:) ! Specific_heat_of_dry_air_at_constant_pressure (J/kg/K)
//     // real(kind_phys),  intent(in)    :: rair(:,:)  ! Gas constant of dry air (J/kg/K)
//     // real(kind_phys),  intent(in)    :: rho(:,:)   ! Dry air density (kg/m^3)
//     // real(kind_phys),  intent(in)    :: z_mid(:,:) ! Heights of thermo. levels (m)
//     // real(kind_phys),  intent(in)    :: pk(:,:)    ! Exner function (p/p0)**(R/cp)

//     // real(kind_phys),  intent(inout) :: phis(:)    ! Surface geopotential (m^2/s^2)

//     // Needed for kessler_run
//     view_2d<Pack>  rho;
//     view_2d<Pack>  dz; // Helper, doesn't need a fortran view
//     view_2d<Pack>  pk;
//     view_2d<Pack>  z_mid;
//     view_2d<Pack>  z_int; // Helper, doesn't need a fortran view and is an interface variable. 
//     // Needed for the kessler_update
//     view_1d<Scalar>  phis;

//     // kessler_run Fortran holders/in Fortran format
//     fview_2dl<Real>  f_cpair;
//     fview_2dl<Real>  f_rair;
//     fview_2dl<Real>  f_rho;
//     fview_2dl<Real>  f_pk;
//     fview_2dl<Real>  f_z_mid;
//     // kessler_update Fortran holders/in Fortran format
//     fview_1d<Real>   f_phis;

//     #if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
//       // Host mirror views for passing to Fortran code on CPU
      
//       view_2dh<Real>   h_cpair;
//       view_2dh<Real>   h_rair;
//       view_2dh<Real>   h_rho;
//       view_2dh<Real>   h_pk;
//       view_2dh<Real>   h_z_mid;
//       view_1dh<Real>   h_phis;
//     #endif

//     // Set number of variables for ATMBufferManager
//     static constexpr int num_1d_intgr   = 0;  // number of 1D integer views
//     static constexpr int num_1d_scalr   = 1;  // number of 1D scalar views (phis or f_phis)
//     static constexpr int num_2d_c       = 4;  // number of 2D fields (dz, rho, pk, z_mid)
//     static constexpr int num_2d_f       = 5;  // number of 2D fields (f_cpair, f_rair, f_rho, f_pk, f_z_mid)
//     static constexpr int num_2d_intlv_c = 1;  // for z_int, which is an interface variable.
    
//     // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
//     void init(int ncol_in, int pver_in) { // TODO - Kokko-ize this
//       using PC  = scream::physics::Constants<Real>;

//       const Real cpair  = PC::Cpair.value; // Specific heat of dry air at constant pressure
//       const Real Rair   = PC::Rair.value;  // Gas constant of dry air

//       Real init_fill_value = 0;

//       // using MDPolicy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
//       // MDPolicy mdp({0,0}, {ncol_in, pver_in});

//       // Kokkos::parallel_for(Kokkos::MDRangePolicy<typename KT::ExeSpace, Kokkos::Rank<2>>({0, 0}, {ncol_in, pver_in}),
//       //   KOKKOS_CLASS_LAMBDA(const int i, const int j) {
//       Kokkos::parallel_for(
//         "init params_computed", KT::RangePolicy(0, ncol_in * pver_in),
//         KOKKOS_CLASS_LAMBDA(const int i) {
//           const int icol = i / pver_in;
//           const int klev = i % pver_in;
//           f_cpair(icol, klev) = cpair;
//           f_rair(icol, klev) = Rair;
//           f_rho(icol, klev) = init_fill_value;
//           f_pk(icol, klev) = init_fill_value;
//           f_z_mid(icol, klev) = init_fill_value;

//           f_phis(icol) = init_fill_value;
//           }
//         );

//     }; // End init

//     // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
//     template <ekat::TransposeDirection::Enum D>
//     void transpose(int ncol_in, int pver_in) { 
//       // auto pver_in_packs = ekat::npack<Pack>(pver_in);

//       if (D == ekat::TransposeDirection::c2f) {

//         Kokkos::parallel_for(
//           "transpose c2f", KT::RangePolicy(0, ncol_in * pver_in),
//           KOKKOS_CLASS_LAMBDA(const int i) {
//             const int icol            = i / pver_in;
//             const int klev            = i % pver_in;
//             // Don't need to transpose cpair, rair

//             f_rho(icol, klev) = rho(icol, klev / Pack::n)[klev % Pack::n];
//             f_pk(icol, klev) = pk(icol, klev / Pack::n)[klev % Pack::n];
//             f_z_mid(icol, klev) = z_mid(icol, klev / Pack::n)[klev % Pack::n];

//             f_phis(icol) = phis(icol);
//           }
//         );
//         #if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
//           // Copy from device to host mirrors for Fortran
          
//           Kokkos::deep_copy(h_cpair, f_cpair);
//           Kokkos::deep_copy(h_rair,  f_rair);
//           Kokkos::deep_copy(h_rho,   f_rho);
//           Kokkos::deep_copy(h_pk,    f_pk);
//           Kokkos::deep_copy(h_z_mid, f_z_mid);
//           Kokkos::deep_copy(h_phis,  f_phis);
//         #endif
//       }
//       if (D == ekat::TransposeDirection::f2c) {  // Not needed but leaving in just in case/temporary
//         #if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
//           // Copy from host mirrors back to device
          
//           Kokkos::deep_copy(f_rho,   h_rho);
//           Kokkos::deep_copy(f_pk,    h_pk);
//           Kokkos::deep_copy(f_z_mid, h_z_mid);
//           Kokkos::deep_copy(f_phis,  h_phis);
//         #endif

//         Kokkos::parallel_for(
//             "transpose f2c", KT::RangePolicy(0, ncol_in * pver_in),
//             KOKKOS_CLASS_LAMBDA(const int i) {
//               const int icol = i / pver_in;
//               const int klev = i % pver_in;
//               // Don't need to transpose cpair, rair

//               rho(icol, klev / Pack::n)[klev % Pack::n] = f_rho(icol, klev);
//               pk(icol, klev / Pack::n)[klev % Pack::n] = f_pk(icol, klev);
//               z_mid(icol, klev / Pack::n)[klev % Pack::n] = f_z_mid(icol, klev);

//               phis(icol) = f_phis(icol); // could skip this
//             }
//           );
//       }
//     }; // End transpose
// }; // End Struct params_helpers

// struct params_computed {
//     params_computed() = default;
//     // Parameters that are computed/updated by kessler 

//     // real(kind_phys),  intent(inout) :: theta(:,:) ! Potential temperature (K)
//     // real(kind_phys),  intent(inout) :: qv(:,:)    ! Water vapor mixing ratio wrt dry air (kg/kg)
//     // real(kind_phys),  intent(inout) :: qc(:,:)    ! Cloud water mixing ratio wrt dry air (kg/kg)
//     // real(kind_phys),  intent(inout) :: qr(:,:)    ! Rain water mixing ratio wrt dry air (kg/kg)
//     // real(kind_phys),  intent(out)   :: precl(:)   ! Precipitation rate (m_water / s)
//     // real(kind_phys),  intent(out)   :: relhum(:,:)! Relative humidity in percent

//     // real(kind_phys),  intent(inout) :: temp_prev(:,:)  ! Previous temperature (K)
//     // real(kind_phys),  intent(inout) :: temp(:,:)       ! Current temperature (K)
//     // real(kind_phys),  intent(inout) :: temp_tend(:,:)  ! Temperature tendency (K/s)
//     // real(kind_phys),  intent(inout) :: st_energy(:,:)  ! Surface energy (J/m^2)

//     // kessler_run C++ view
//     view_2d<Pack>  theta;
//     view_2d<Pack>  qv;
//     view_2d<Pack>  qc;
//     view_2d<Pack>  qr;
//     view_1d<Scalar> precl;
//     view_2d<Pack>  relhum;
//     // kessler_update C++ view
//     view_2d<Pack>   temp_prev;
//     view_2d<Pack>   temp;
//     view_2d<Pack>   temp_tend;
//     view_2d<Pack>   st_energy;
  
//     // kessler_run fortran (left layout) versions
//     fview_2dl<Real>  f_theta;
//     fview_2dl<Real>  f_qv;
//     fview_2dl<Real>  f_qc;
//     fview_2dl<Real>  f_qr;
//     fview_1d<Real>   f_precl;  
//     fview_2dl<Real>  f_relhum;
//     // kessler_update fortran (left layout) versions
//     fview_2dl<Real>  f_temp_prev;
//     fview_2dl<Real>  f_temp;
//     fview_2dl<Real>  f_temp_tend;
//     fview_2dl<Real>  f_st_energy;

//     #if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
//       // Host mirror views for passing to Fortran code on CPU
      
//       view_2dh<Real>   h_theta;
//       view_2dh<Real>   h_qv;
//       view_2dh<Real>   h_qc;
//       view_2dh<Real>   h_qr;
//       view_1dh<Real>   h_precl;
//       view_2dh<Real>   h_relhum;
//       view_2dh<Real>   h_temp_prev;
//       view_2dh<Real>   h_temp;
//       view_2dh<Real>   h_temp_tend;
//       view_2dh<Real>   h_st_energy;
//     #endif

//     // Set number of variables for ATMBufferManager
//     static constexpr int num_1d_intgr   = 0; // number of 1D integer views
//     static constexpr int num_1d_scalr   = 1; // number of 1D scalar views (precl or f_precl)
//     static constexpr int num_2d_c       = 9; // number of 2D fields (theta, qv, qc, qr, relhum, temp_prev, temp, temp_tend, st_energy)
//     static constexpr int num_2d_f       = 9; // number of 2D fields (f_theta, f_qv, f_qc, f_qr, f_relhum, f_temp_prev, f_temp, f_temp_tend, f_st_energy)
//     static constexpr int num_2d_intlv_c = 0; // no interface variables here.

//     // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
//     void init(int ncol_in, int pver_in) {
//       Real init_fill_value = 0.0;

//       // using MDPolicy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
//       // MDPolicy mdp({0,0}, {ncol_in, pver_in});
      
//       // Kokkos::parallel_for(Kokkos::MDRangePolicy<typename KT::ExeSpace, Kokkos::Rank<2>>({0, 0}, {ncol_in, pver_in}),
//       //   KOKKOS_CLASS_LAMBDA(const int i, const int j) {
//       Kokkos::parallel_for(
//         "init params_computed", KT::RangePolicy(0, ncol_in * pver_in),
//         KOKKOS_CLASS_LAMBDA(const int i) {
//           const int icol = i / pver_in;
//           const int klev = i % pver_in;
//           f_theta(icol, klev) = init_fill_value;
//           f_qv(icol, klev) = init_fill_value;
//           f_qc(icol, klev) = init_fill_value;
//           f_qr(icol, klev) = init_fill_value;
//           f_relhum(icol, klev) = init_fill_value;

//           f_precl(icol) = init_fill_value;

//           f_temp_prev(icol, klev) = init_fill_value;
//           f_temp(icol, klev) = init_fill_value;
//           f_temp_tend(icol, klev) = init_fill_value;
//           f_st_energy(icol, klev) = init_fill_value;
//           }
//         );
//     }; // End init

//     // Modified from the ZM implementation in components/eamxx/src/physics/zm/zm_functions.hpp
//     template <ekat::TransposeDirection::Enum D>
//     void transpose(int ncol_in, int pver_in) { // TODO - Kokko-ize this
//       // auto pver_in_packs = ekat::npack<Pack>(pver_in);

//       if (D == ekat::TransposeDirection::c2f) {

//         Kokkos::parallel_for(
//           "transpose c2f", KT::RangePolicy(0, ncol_in * pver_in),
//           KOKKOS_CLASS_LAMBDA(const int i) {
//             const int icol = i / pver_in;
//             const int klev = i % pver_in;
//             f_theta(icol, klev) = theta(icol, klev / Pack::n)[klev % Pack::n];
//             f_qv(icol, klev) = qv(icol, klev / Pack::n)[klev % Pack::n];
//             f_qc(icol, klev) = qc(icol, klev / Pack::n)[klev % Pack::n];
//             f_qr(icol, klev) = qr(icol, klev / Pack::n)[klev % Pack::n];
//             f_relhum(icol, klev) = relhum(icol, klev / Pack::n)[klev % Pack::n];
//             f_precl(icol) = precl(icol);

//             f_temp_prev(icol, klev) = temp_prev(icol, klev / Pack::n)[klev % Pack::n];
//             f_temp(icol, klev) = temp(icol, klev / Pack::n)[klev % Pack::n];
//             f_temp_tend(icol, klev) = temp_tend(icol, klev / Pack::n)[klev % Pack::n];
//             f_st_energy(icol, klev) = st_energy(icol, klev / Pack::n)[klev % Pack::n];
            
//           }
//         );
//         #if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
//           // Copy from device to host mirrors for Fortran
          
//           Kokkos::deep_copy(h_theta, f_theta);
//           Kokkos::deep_copy(h_qv, f_qv);
//           Kokkos::deep_copy(h_qc, f_qc);
//           Kokkos::deep_copy(h_qr, f_qr);
//           Kokkos::deep_copy(h_relhum, f_relhum);
//           Kokkos::deep_copy(h_precl, f_precl);
//           Kokkos::deep_copy(h_temp_prev, f_temp_prev);
//           Kokkos::deep_copy(h_temp, f_temp);
//           Kokkos::deep_copy(h_temp_tend, f_temp_tend);
//           Kokkos::deep_copy(h_st_energy, f_st_energy);
//         #endif
//       }
//       if (D == ekat::TransposeDirection::f2c) {
//         #if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
//           // Copy from host mirrors back to device
          
//           Kokkos::deep_copy(f_theta, h_theta);
//           Kokkos::deep_copy(f_qv, h_qv);
//           Kokkos::deep_copy(f_qc, h_qc);
//           Kokkos::deep_copy(f_qr, h_qr);
//           Kokkos::deep_copy(f_relhum, h_relhum);
//           Kokkos::deep_copy(f_precl, h_precl);
//           Kokkos::deep_copy(f_temp_prev, h_temp_prev);
//           Kokkos::deep_copy(f_temp, h_temp);
//           Kokkos::deep_copy(f_temp_tend, h_temp_tend);
//           Kokkos::deep_copy(f_st_energy, h_st_energy);
//         #endif

//         Kokkos::parallel_for(
//             "transpose f2c", KT::RangePolicy(0, ncol_in * pver_in),
//             KOKKOS_CLASS_LAMBDA(const int i) {
//               const int icol = i / pver_in;
//               const int klev  = i % pver_in;
//               theta(icol, klev / Pack::n)[klev % Pack::n] = f_theta(icol, klev);
//               qv(icol, klev / Pack::n)[klev % Pack::n] = f_qv(icol, klev);
//               qc(icol, klev / Pack::n)[klev % Pack::n] = f_qc(icol, klev);
//               qr(icol, klev / Pack::n)[klev % Pack::n] = f_qr(icol, klev);
//               relhum(icol, klev / Pack::n)[klev % Pack::n] = f_relhum(icol, klev);

//               precl(icol) = f_precl(icol);

//               temp_prev(icol, klev / Pack::n)[klev % Pack::n] = f_temp_prev(icol, klev);
//               temp(icol, klev / Pack::n)[klev % Pack::n] = f_temp(icol, klev);
//               temp_tend(icol, klev / Pack::n)[klev % Pack::n] = f_temp_tend(icol, klev);
//               st_energy(icol, klev / Pack::n)[klev % Pack::n] = f_st_energy(icol, klev);
//             }
//           );
//       }
//     }; // End transpose
// }; // End params_computed 

}; // struct KesslerMicrophysicsFunctions
  
} // namespace kessler
} // namespace scream

#endif // KESSLER_FUNCTIONS_HPP