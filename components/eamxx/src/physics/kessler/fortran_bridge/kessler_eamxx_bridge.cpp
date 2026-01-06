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

    void kessler_eamxx_bridge_init_c(const Real lv_in, const Real pref_in, const Real rhoqr_in, char* errmsg, Int& errflg);

    void kessler_eamxx_bridge_run_c(Int pcols, Int nz, Real* dt, Int lyr_surf, Int lyr_toa, Real* cpair, Real* rair, Real* rho, Real* z, Real* pk, Real* theta, Real* qv, Real* qc, Real* qr, Real* precl, Real* relhum, char* scheme_name, char* errmsg, Int& errflg); 
} // extern "C" : end _c decls

namespace scream {
    namespace kessler {

    void set_log_file_name_f90(const char** logname){
        set_log_file_name_f90_c(logname);
    }

    void kessler_eamxx_bridge_init( const Real lv_in, const Real pref_in, const Real rhoqr_in){
        char* errmsg;
        Int errflg = 0;
        kessler_eamxx_bridge_init_c( lv_in, pref_in, rhoqr_in, errmsg, &errflg );
    }

    void kessler_eamxx_bridge_run( Int pcols, Int pver, Real* dt, Int lyr_surf, Int lyr_toa, KMF:params_in &params_in, KMF:params_out &params_out ){ 

        // Error catching for Fortran Kessler
        char* scheme_name; 
        char* errmsg;
        Int errflg = 0;
        //----------------------------------------------------------------------------
        // Need to transpose to match how Fortran handles things
        params_in.transpose<ekat::TransposeDirection::c2f>(pcols,pver);
        params_out.transpose<ekat::TransposeDirection::c2f>(pcols,pver); // needed for updated values

        kessler_eamxx_bridge_run_c(pcols, pver, dt, lyr_surf, lyr_toa, params_in.f_cpair.data(), 
                                                                        params_in.f_rair.data(), 
                                                                        params_in.f_rho.data(), 
                                                                        params_in.f_z.data(),
                                                                        params_in.f_pk.data(), 
                                                                        params_out.f_theta.data(), 
                                                                        params_out.f_qv.data(), 
                                                                        params_out.f_qc.data(), 
                                                                        params_out.f_qr.data(), 
                                                                        params_out.f_precl.data(), 
                                                                        params_out.f_relhum.data(), 
                                                                        scheme_name, errmsg, &errflg);

        // Transpose back to C++ convention
        params_out.transpose<ekat::TransposeDirection::f2c>(pcols,pver);

        //----------------------------------------------------------------------------
        }

    // end _c impls

    } // namespace kessler
} // namespace scream