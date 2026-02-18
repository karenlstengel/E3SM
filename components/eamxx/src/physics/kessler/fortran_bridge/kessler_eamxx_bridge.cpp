#include "kessler_eamxx_bridge.hpp"
#include "kessler_functions.hpp"
#include "share/core/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_workspace.hpp>

using scream::Real;
using scream::Int;

// A C++ interface to kessler fortran calls and vice versa

extern "C" {

    void kessler_eamxx_bridge_init_c(Int pcols, Int pver, Real lv_in, Real pref_in, Real rhoqr_in);

    void kessler_eamxx_bridge_update_init_c(Int pcols, Int pver, Real gravit_in);

    void kessler_eamxx_bridge_run_c(Int pcols, Int nz, double dt, Int lyr_surf, Int lyr_toa, Real* cpair, Real* rair, Real* rho, Real* z, Real* pk, Real* theta, Real* qv, Real* qc, Real* qr, Real* precl, Real* relhum); 

    void kessler_eamxx_bridge_update_c(Int pcols, Int nz, double dt, Real* cpair, Real* pk, Real* theta, Real* temp_prev, Real* temp, Real* temp_tend, Real* z_mid, Real* phis, Real* st_energy);
} // extern "C" : end _c decls

namespace scream {
    namespace kessler {

    void kessler_eamxx_bridge_init( Int pcols, Int pver, Real lv_in, Real pref_in, Real rhoqr_in, Real gravit_in){
        kessler_eamxx_bridge_init_c( pcols, pver, lv_in, pref_in, rhoqr_in);
        kessler_eamxx_bridge_update_init_c( pcols, pver, gravit_in);
    }

    void kessler_eamxx_bridge_run( Int pcols, Int pver, double dt, Int lyr_surf, Int lyr_toa, KMF::params_in &params_in, KMF::params_out &params_out ){ 

        //----------------------------------------------------------------------------
        // Need to transpose to match how Fortran handles things
        params_in.transpose<ekat::TransposeDirection::c2f>(pcols,pver);
        params_out.transpose<ekat::TransposeDirection::c2f>(pcols,pver); // needed for updated values
        
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
        } // end run

    void kessler_eamxx_bridge_update(Int pcols, Int pver, double dt, KMF::params_in &params_in, KMF::params_out &params_out, KMF::params_update &params_update){ 
        //----------------------------------------------------------------------------
        // Need to transpose to match how Fortran handles things
        params_in.transpose<ekat::TransposeDirection::c2f>(pcols,pver);
        params_out.transpose<ekat::TransposeDirection::c2f>(pcols,pver);
        params_update.transpose<ekat::TransposeDirection::c2f>(pcols,pver); // needed for updated values

        // Call fortran update function here when available
        kessler_eamxx_bridge_update_c(pcols, pver, dt, params_in.f_cpair.data(),
                                                       params_in.f_pk.data(),
                                                       params_out.f_theta.data(), 
                                                       params_update.f_temp_prev.data(),  
                                                       params_update.f_temp.data(), 
                                                       params_update.f_temp_tend.data(),
                                                       params_update.f_z_mid.data(), 
                                                       params_update.f_phis.data(), 
                                                       params_update.f_st_energy.data());


        // Transpose back to C++ convention
        params_update.transpose<ekat::TransposeDirection::f2c>(pcols,pver);

        //----------------------------------------------------------------------------
        }


    // end _c impls

    } // namespace kessler
} // namespace scream