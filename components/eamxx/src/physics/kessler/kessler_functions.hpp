#ifndef KESSLER_FUNCTIONS_HPP
#define KESSLER_FUNCTIONS_HPP

#include "share/physics/physics_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/core/eamxx_types.hpp"
#include "share/util/eamxx_timing.hpp"

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
  using IntPack = ekat::Pack<Int,SCREAM_PACK_SIZE>;

  using KT      = ekat::KokkosTypes<Device>;
  // using MemberType = typename KT::MemberType;
  using TeamPolicy = typename KokkosTypes<Device>::TeamPolicy;

  template <typename S> using view_1d   = typename KT::template view_1d<S>;
  template <typename S> using view_2d   = typename KT::template view_2d<S>;
  template <typename S> using view_2dl  = typename KT::template lview<S**>;
  
  #if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
    // Needed for running the Fortran without OpenACC
    template <typename S> using fview_1d  = typename ekat::template Unmanaged<view_1d<S> >;
    template <typename S> using fview_2d  = typename ekat::template Unmanaged<view_2d<S> >;
    template <typename S> using fview_2dl = typename ekat::template Unmanaged<view_2dl<S> >;
    // Host copies of the Fortran-layout arrays, carved from one host block
    // (see KesslerMicrophysics::initialize_impl), so a whole group of them
    // moves between host and device in one copy
    template <typename S> using view_1dh  = Kokkos::View<S*,  Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
    template <typename S> using view_2dh  = Kokkos::View<S**, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
  #else 
    template <typename S> using fview_1d  = typename KT::template view_1d<S>;
    template <typename S> using fview_2d  = typename KT::template view_2d<S>;
    template <typename S> using fview_2dl = typename KT::template lview<S **>;
  #endif

  // C++ <-> Fortran layout copies (transposes): one work item per (column,
  // level), with no integer division per item. The index that varies
  // fastest (across a GPU warp, or in the innermost CPU loop) is the one the
  // written array is contiguous in: columns when writing the Fortran
  // (LayoutLeft) arrays, levels when writing EAMxx's (column, level-pack)
  // arrays. So writes are contiguous and only the reads are strided, the
  // cheaper side (c2f, which wrote strided, cost ~4x f2c on the GPU).
  using C2FPolicy = Kokkos::MDRangePolicy<typename KT::ExeSpace,
                      Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>;
  using F2CPolicy = Kokkos::MDRangePolicy<typename KT::ExeSpace,
                      Kokkos::Rank<2, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>;

  // ----------------------------------------
  // Structs
struct params_helpers {
    params_helpers() = default;

    // Needed for kessler_run or to compute things to pass to kessler
    // real(kind_phys),  intent(in)    :: rho(:,:)   ! Dry air density (kg/m^3)
    // real(kind_phys),  intent(in)    :: z_mid(:,:) ! Heights of thermo. levels (m)
    // real(kind_phys),  intent(in)    :: pk(:,:)    ! Exner function (p/p0)**(R/cp)

    // real(kind_phys),  intent(inout) :: phis(:)    ! Surface geopotential (m^2/s^2)

    // Needed for kessler_run
    view_2d<Pack>  rho;
    view_2d<Pack>  dz; // Helper, doesn't need a fortran view
    view_2d<Pack>  pk;
    view_2d<Pack>  z_mid;
    view_2d<Pack>  z_int; // Helper, doesn't need a fortran view and is an interface variable. 
    // Needed for the kessler_update
    view_1d<Scalar>  phis;

    // kessler_run Fortran holders/in Fortran format
    fview_2dl<Real>  f_rho;
    fview_2dl<Real>  f_pk;
    fview_2dl<Real>  f_z_mid;
    // kessler_update Fortran holders/in Fortran format
    fview_1d<Real>   f_phis;

    #if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
      // Host mirror views for passing to Fortran code on CPU

      view_2dh<Real>   h_rho;
      view_2dh<Real>   h_pk;
      view_2dh<Real>   h_z_mid;
      view_1dh<Real>   h_phis;
    #endif

    // Set number of variables for ATMBufferManager
    static constexpr int num_1d_intgr   = 0;  // number of 1D integer views
    static constexpr int num_1d_scalr   = 1;  // number of 1D scalar views (phis or f_phis)
    static constexpr int num_2d_c       = 4;  // number of 2D fields (dz, rho, pk, z_mid)
    static constexpr int num_2d_f       = 3;  // number of 2D fields (f_rho, f_pk, f_z_mid)
    static constexpr int num_2d_intlv_c = 1;  // for z_int, which is an interface variable.

    // Fortran-layout copies of the inputs. They are inputs only, so nothing
    // is copied back (f2c is a no-op), and they are not zero-filled first:
    // c2f overwrites every element.
    template <ekat::TransposeDirection::Enum D>
    void transpose(int ncol_in, int pver_in) {
      if (D == ekat::TransposeDirection::c2f) {
        Kokkos::parallel_for(
          "transpose c2f helpers", C2FPolicy({0, 0}, {ncol_in, pver_in}),
          KOKKOS_CLASS_LAMBDA(const int icol, const int klev) {
            f_rho(icol, klev)   = rho(icol, klev / Pack::n)[klev % Pack::n];
            f_pk(icol, klev)    = pk(icol, klev / Pack::n)[klev % Pack::n];
            f_z_mid(icol, klev) = z_mid(icol, klev / Pack::n)[klev % Pack::n];
            if (klev == 0) f_phis(icol) = phis(icol);
          }
        );
      }
    }; // End transpose
}; // End Struct params_helpers

struct params_computed {
    params_computed() = default;
    // Parameters that are computed/updated by kessler 

    // real(kind_phys),  intent(inout) :: theta(:,:) ! Potential temperature (K)
    // real(kind_phys),  intent(inout) :: qv(:,:)    ! Water vapor mixing ratio wrt dry air (kg/kg)
    // real(kind_phys),  intent(inout) :: qc(:,:)    ! Cloud water mixing ratio wrt dry air (kg/kg)
    // real(kind_phys),  intent(inout) :: qr(:,:)    ! Rain water mixing ratio wrt dry air (kg/kg)
    // real(kind_phys),  intent(out)   :: precl(:)   ! Precipitation rate (m_water / s)
    // real(kind_phys),  intent(out)   :: relhum(:,:)! Relative humidity in percent

    // real(kind_phys),  intent(inout) :: temp_prev(:,:)  ! Previous temperature (K)
    // real(kind_phys),  intent(inout) :: temp(:,:)       ! Current temperature (K)
    // real(kind_phys),  intent(inout) :: temp_tend(:,:)  ! Temperature tendency (K/s)
    // real(kind_phys),  intent(inout) :: st_energy(:,:)  ! Surface energy (J/m^2)

    // kessler_run C++ view
    view_2d<Pack>  theta;
    view_2d<Pack>  qv;
    view_2d<Pack>  qc;
    view_2d<Pack>  qr;
    view_1d<Scalar> precl;
    view_2d<Pack>  relhum;
    // kessler_update C++ view
    view_2d<Pack>   temp_prev;
    view_2d<Pack>   temp;
    view_2d<Pack>   temp_tend;
    view_2d<Pack>   st_energy;
  
    // kessler_run fortran (left layout) versions
    fview_2dl<Real>  f_theta;
    fview_2dl<Real>  f_qv;
    fview_2dl<Real>  f_qc;
    fview_2dl<Real>  f_qr;
    fview_1d<Real>   f_precl;  
    fview_2dl<Real>  f_relhum;
    // kessler_update fortran (left layout) versions
    fview_2dl<Real>  f_temp_prev;
    fview_2dl<Real>  f_temp;
    fview_2dl<Real>  f_temp_tend;
    fview_2dl<Real>  f_st_energy;

    #if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
      // Host mirror views for passing to Fortran code on CPU
      
      view_2dh<Real>   h_theta;
      view_2dh<Real>   h_qv;
      view_2dh<Real>   h_qc;
      view_2dh<Real>   h_qr;
      view_1dh<Real>   h_precl;
      view_2dh<Real>   h_relhum;
      view_2dh<Real>   h_temp_prev;
      view_2dh<Real>   h_temp;
      view_2dh<Real>   h_temp_tend;
      view_2dh<Real>   h_st_energy;
    #endif

    // Set number of variables for ATMBufferManager
    static constexpr int num_1d_intgr   = 0; // number of 1D integer views
    static constexpr int num_1d_scalr   = 1; // number of 1D scalar views (precl or f_precl)
    static constexpr int num_2d_c       = 9; // number of 2D fields (theta, qv, qc, qr, relhum, temp_prev, temp, temp_tend, st_energy)
    static constexpr int num_2d_f       = 9; // number of 2D fields (f_theta, f_qv, f_qc, f_qr, f_relhum, f_temp_prev, f_temp, f_temp_tend, f_st_energy)
    static constexpr int num_2d_intlv_c = 0; // no interface variables here.

    // Copies between EAMxx's arrays and the Fortran-layout ones, only of what
    // the Fortran needs: c2f the arrays it reads (theta, qv, qc, qr: inout in
    // kessler_run; temp: inout in the update), f2c the ones it writes (those,
    // plus relhum, precl, temp_prev, temp_tend, st_energy, which it only
    // writes, every element, so they are not copied in). Nothing is
    // zero-filled first: c2f overwrites every element it copies, and the
    // Fortran overwrites the rest.
    template <ekat::TransposeDirection::Enum D>
    void transpose(int ncol_in, int pver_in) {
      if (D == ekat::TransposeDirection::c2f) {
        Kokkos::parallel_for(
          "transpose c2f computed", C2FPolicy({0, 0}, {ncol_in, pver_in}),
          KOKKOS_CLASS_LAMBDA(const int icol, const int klev) {
            const int kp = klev / Pack::n, kl = klev % Pack::n;
            f_theta(icol, klev) = theta(icol, kp)[kl];
            f_qv(icol, klev)    = qv(icol, kp)[kl];
            f_qc(icol, klev)    = qc(icol, kp)[kl];
            f_qr(icol, klev)    = qr(icol, kp)[kl];
            f_temp(icol, klev)  = temp(icol, kp)[kl];
          }
        );
      }
      if (D == ekat::TransposeDirection::f2c) {
        Kokkos::parallel_for(
          "transpose f2c computed", F2CPolicy({0, 0}, {ncol_in, pver_in}),
          KOKKOS_CLASS_LAMBDA(const int icol, const int klev) {
            const int kp = klev / Pack::n, kl = klev % Pack::n;
            theta(icol, kp)[kl]     = f_theta(icol, klev);
            qv(icol, kp)[kl]        = f_qv(icol, klev);
            qc(icol, kp)[kl]        = f_qc(icol, klev);
            qr(icol, kp)[kl]        = f_qr(icol, klev);
            relhum(icol, kp)[kl]    = f_relhum(icol, klev);
            temp_prev(icol, kp)[kl] = f_temp_prev(icol, klev);
            temp(icol, kp)[kl]      = f_temp(icol, klev);
            temp_tend(icol, kp)[kl] = f_temp_tend(icol, klev);
            st_energy(icol, kp)[kl] = f_st_energy(icol, klev);
            if (klev == 0) precl(icol) = f_precl(icol);
          }
        );
      }
    }; // End transpose
}; // End params_computed 

}; // struct KesslerMicrophysicsFunctions
  
} // namespace kessler
} // namespace scream

#endif // KESSLER_FUNCTIONS_HPP