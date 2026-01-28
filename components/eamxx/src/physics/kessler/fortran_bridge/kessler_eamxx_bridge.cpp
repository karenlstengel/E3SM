#include "kessler_eamxx_bridge.hpp"
#include "kessler_functions.hpp"
#include "share/core/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_workspace.hpp>

using scream::Real;
using scream::Int;

// A C++ interface to kessler fortran calls and vice versa

extern "C" {
    void set_log_file_name_f90_c(const char** fname);

    void kessler_eamxx_bridge_init_c(Int pcols, Int pver, Real lv_in, Real pref_in, Real rhoqr_in);

    void kessler_eamxx_bridge_run_c(Int pcols, Int nz, double dt, Int lyr_surf, Int lyr_toa, Real* cpair, Real* rair, Real* rho, Real* z, Real* pk, Real* theta, Real* qv, Real* qc, Real* qr, Real* precl, Real* relhum); 
} // extern "C" : end _c decls

namespace scream {
    namespace kessler {

    void set_log_file_name_f90(const char** logname){
        set_log_file_name_f90_c(logname);
    }

    void kessler_eamxx_bridge_init( Int pcols, Int pver, Real lv_in, Real pref_in, Real rhoqr_in){
        kessler_eamxx_bridge_init_c( pcols, pver, lv_in, pref_in, rhoqr_in);
    }

    void kessler_eamxx_bridge_run( Int pcols, Int pver, double dt, Int lyr_surf, Int lyr_toa, KMF::params_in &params_in, KMF::params_out &params_out ){ 

        //----------------------------------------------------------------------------
        // Need to transpose to match how Fortran handles things
        params_in.transpose<ekat::TransposeDirection::c2f>(pcols,pver);
        params_out.transpose<ekat::TransposeDirection::c2f>(pcols,pver); // needed for updated values

        for(int i=0; i<pcols; ++i) {
            for(int j=0; j<pver; ++j) {
                // Just to check values before calling fortran
                if (params_in.f_cpair(i,j) != params_in.f_cpair(0,0))
                {
                    printf("cpair(%d,%d) before f90 call: %f\n", i, j, params_in.f_cpair(i,j));
                }
            }
        }
        
        kessler_eamxx_bridge_run_c(pcols, pver, dt, lyr_surf, lyr_toa, params_in.f_cpair.data(), 
                                                                        params_in.f_rair.data(), 
                                                                        params_in.f_rho.data(), 
                                                                        params_in.f_dz.data(),
                                                                        params_in.f_pk.data(), 
                                                                        params_out.f_theta.data(), 
                                                                        params_out.f_qv.data(), 
                                                                        params_out.f_qc.data(), 
                                                                        params_out.f_qr.data(), 
                                                                        params_out.f_precl.data(), 
                                                                        params_out.f_relhum.data());

        // Transpose back to C++ convention
        params_out.transpose<ekat::TransposeDirection::f2c>(pcols,pver);

        //----------------------------------------------------------------------------
        }

    // end _c impls

    } // namespace kessler
} // namespace scream