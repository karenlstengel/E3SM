#include "kessler_eamxx_bridge.hpp"
#include "kessler_functions.hpp"
#include "share/core/eamxx_types.hpp"
#include "share/util/eamxx_timing.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_workspace.hpp>

using scream::Real;
using scream::Int;

// A C++ interface to kessler fortran calls and vice versa

extern "C" {

    void kessler_eamxx_bridge_init_c(Int pcols, Int pver, Real lv_in, Real pref_in, Real rhoqr_in, Real cpair_in, Real rair_in);

    void kessler_eamxx_bridge_update_init_c(Int pcols, Int pver, Real gravit_in, Real cpair_in);

    void kessler_eamxx_bridge_run_c(Int pcols, Int nz, double dt, Int lyr_surf, Int lyr_toa, Real* rho, Real* z_mid, Real* pk, Real* theta, Real* qv, Real* qc, Real* qr, Real* precl, Real* relhum);

    void kessler_eamxx_bridge_update_c(Int pcols, Int nz, double dt, Real* pk, Real* theta, Real* temp_prev, Real* temp, Real* temp_tend, Real* z_mid, Real* phis, Real* st_energy);

    void kessler_eamxx_bridge_finalize_c();
} // extern "C" : end _c decls

namespace scream {
    namespace kessler {

    #if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
    // n consecutive Reals between device and host memory, as one copy
    // (Fortran on the host, e.g. the gpucpu build). The Fortran-layout arrays
    // are carved back to back in a fixed order on both sides
    // (KesslerMicrophysics::init_buffers / initialize_impl), so each group
    // the Fortran reads or writes is one range.
    template <typename DevView>
    void copy_range (Real* dst, const Real* src, const size_t n, const bool to_host) {
      using DevMem = typename DevView::memory_space;
      using DevRange  = Kokkos::View<Real*, DevMem, Kokkos::MemoryUnmanaged>;
      using HostRange = Kokkos::View<Real*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
      if (to_host) Kokkos::deep_copy(HostRange(dst, n), DevRange(const_cast<Real*>(src), n));
      else         Kokkos::deep_copy(DevRange(dst, n), HostRange(const_cast<Real*>(src), n));
    }
    #endif

    void kessler_eamxx_bridge_init( Int pcols, Int pver, Real lv_in, Real pref_in, Real rhoqr_in, Real cpair_in, Real rair_in, Real gravit_in){
        kessler_eamxx_bridge_init_c( pcols, pver, lv_in, pref_in, rhoqr_in, cpair_in, rair_in);
        kessler_eamxx_bridge_update_init_c( pcols, pver, gravit_in, cpair_in);
    }

    void kessler_eamxx_bridge_run( Int pcols, Int pver, double dt, Int lyr_surf, Int lyr_toa, KMF::params_helpers &params_helpers, KMF::params_computed &params_computed ){ 

        //----------------------------------------------------------------------------
        // Need to transpose to match how Fortran handles things
        start_timer("EAMxx::kessler::run::F90_run::transpose_c2f");
        params_helpers.transpose<ekat::TransposeDirection::c2f>(pcols,pver);
        params_computed.transpose<ekat::TransposeDirection::c2f>(pcols,pver);
        #if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
          // To the host the arrays the Fortran reads: f_rho .. f_temp (one
          // range, 8 arrays) and f_phis
          start_timer("EAMxx::kessler::run::sync_to_host");
          {
            using DV = decltype(params_helpers.f_rho);
            const size_t n2 = size_t(pcols) * pver;
            copy_range<DV>(params_helpers.h_rho.data(), params_helpers.f_rho.data(), 8 * n2, true);
            copy_range<DV>(params_helpers.h_phis.data(), params_helpers.f_phis.data(), pcols, true);
          }
          stop_timer("EAMxx::kessler::run::sync_to_host");
        #endif
        Kokkos::fence();
        stop_timer("EAMxx::kessler::run::F90_run::transpose_c2f");

        #if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
            
            kessler_eamxx_bridge_run_c(pcols, pver, dt, lyr_surf, lyr_toa, params_helpers.h_rho.data(),
                                                                       params_helpers.h_z_mid.data(),
                                                                       params_helpers.h_pk.data(),
                                                                       params_computed.h_theta.data(),
                                                                       params_computed.h_qv.data(),
                                                                       params_computed.h_qc.data(),
                                                                       params_computed.h_qr.data(),
                                                                       params_computed.h_precl.data(),
                                                                       params_computed.h_relhum.data());

            kessler_eamxx_bridge_update_c(pcols, pver, dt,
                                                       params_helpers.h_pk.data(),
                                                       params_computed.h_theta.data(), 
                                                       params_computed.h_temp_prev.data(),  
                                                       params_computed.h_temp.data(), 
                                                       params_computed.h_temp_tend.data(),
                                                       params_helpers.h_z_mid.data(), 
                                                       params_helpers.h_phis.data(), 
                                                       params_computed.h_st_energy.data());
        #else
            kessler_eamxx_bridge_run_c(pcols, pver, dt, lyr_surf, lyr_toa, params_helpers.f_rho.data(),
                                                                       params_helpers.f_z_mid.data(),
                                                                       params_helpers.f_pk.data(),
                                                                       params_computed.f_theta.data(),
                                                                       params_computed.f_qv.data(),
                                                                       params_computed.f_qc.data(),
                                                                       params_computed.f_qr.data(),
                                                                       params_computed.f_precl.data(),
                                                                       params_computed.f_relhum.data());

            kessler_eamxx_bridge_update_c(pcols, pver, dt,
                                                       params_helpers.f_pk.data(),
                                                       params_computed.f_theta.data(), 
                                                       params_computed.f_temp_prev.data(),  
                                                       params_computed.f_temp.data(), 
                                                       params_computed.f_temp_tend.data(),
                                                       params_helpers.f_z_mid.data(), 
                                                       params_helpers.f_phis.data(), 
                                                       params_computed.f_st_energy.data());
        #endif
        // Transpose back to C++ convention (params_helpers: inputs only, nothing to copy back)
        start_timer("EAMxx::kessler::run::F90_run::transpose_f2c");
        #if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
          // To the device the arrays the Fortran writes: f_theta .. f_st_energy
          // (one range, 9 arrays) and f_precl
          start_timer("EAMxx::kessler::run::sync_to_dev");
          {
            using DV = decltype(params_helpers.f_rho);
            const size_t n2 = size_t(pcols) * pver;
            copy_range<DV>(params_computed.f_theta.data(), params_computed.h_theta.data(), 9 * n2, false);
            copy_range<DV>(params_computed.f_precl.data(), params_computed.h_precl.data(), pcols, false);
          }
          stop_timer("EAMxx::kessler::run::sync_to_dev");
        #endif
        params_computed.transpose<ekat::TransposeDirection::f2c>(pcols,pver);
        Kokkos::fence();
        stop_timer("EAMxx::kessler::run::F90_run::transpose_f2c");

        //----------------------------------------------------------------------------
    } // end run

    void kessler_eamxx_bridge_finalize(){
        kessler_eamxx_bridge_finalize_c();
    }

    // end _c impls

    } // namespace kessler
} // namespace scream