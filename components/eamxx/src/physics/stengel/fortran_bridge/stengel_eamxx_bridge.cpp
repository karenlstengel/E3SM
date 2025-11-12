#include "stengel_eamxx_bridge.hpp"
#include "share/core/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_workspace.hpp>

using scream::Real;
using scream::Int;

// A C++ interface to stengel fortran calls and vice versa

extern "C" {
    void stengel_eamxx_bridge_init_c(Int pcol, Int pver );

    void stengel_eamxx_bridge_run_c(Int ncol, Real p_mid, Real T_mid);
} // extern "C" : end _c decls

namespace scream {
    namespace stengel {

    void stengel_eamxx_bridge_init( Int pcol, Int pver ){
        stengel_eamxx_bridge_init_c( pcol, pver );
    }

    void stengel_eamxx_bridge_run( Int ncol, Int pver,
                            const view_2d<Spack>& p_mid,
                            const view_2d<Spack>& T_mid
                            ){ // TODO - required_arguments = p_mid, T_mid
        //----------------------------------------------------------------------------
        // Need to transpose to match how Fortran handles things
        // see zm_input.transpose<ekat::TransposeDirection::c2f>(ncol,pver);
        // see zm_output.transpose<ekat::TransposeDirection::c2f>(ncol,pver);

        view_2d<Real>& p_mid_f; 
        view_2d<Real>& T_mid_f;

        for (Int i=0; i<ncol; ++i) {
            for (Int j=0; j<pver; ++j) {
                p_mid_f(i,j) = p_mid   (i,j/Spack::n)[j%Spack::n];
                T_mid_f(i,j) = T_mid   (i,j/Spack::n)[j%Spack::n];
            }
        }

        stengel_eamxx_bridge_run_c(ncol, p_mid_f, T_mid_f); // TODO - required_arguments = p_mid, T_mid

        // TODO - print out max values in p_mid_f, T_mid_f

        // Transpose back to C++ convention
        // see zm_input.transpose<ekat::TransposeDirection::f2c>(ncol,pver);
        // see zm_output.transpose<ekat::TransposeDirection::f2c>(ncol,pver);

        //----------------------------------------------------------------------------
        }

        // end _c impls

    } // namespace stengel
} // namespace scream