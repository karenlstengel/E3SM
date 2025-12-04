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

    void kessler_eamxx_bridge_init_c(Int pcols, Int pver );

    void kessler_eamxx_bridge_run_c(Int pcols, TODO); // example: Real* p_mid
} // extern "C" : end _c decls

namespace scream {
    namespace kessler {

    void set_log_file_name_f90(const char** logname){
        set_log_file_name_f90_c(logname);
    }

    void kessler_eamxx_bridge_init( Int pcols, Int pver ){
        kessler_eamxx_bridge_init_c( pcols, pver );
    }

    void kessler_eamxx_bridge_run( Int pcols, Int pver, KMF:params &params){ 
        //----------------------------------------------------------------------------
        // Need to transpose to match how Fortran handles things
        params.transpose<ekat::TransposeDirection::c2f>(pcols,pver);

        kessler_eamxx_bridge_run_c(pcols, TODO); // example: params.f_p_mid.data()

        // Transpose back to C++ convention
        params.transpose<ekat::TransposeDirection::f2c>(pcols,pver);

        //----------------------------------------------------------------------------
        }

    // end _c impls

    } // namespace kessler
} // namespace scream