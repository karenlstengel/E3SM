#include "stengelF_eamxx_bridge.hpp"
#include "stengelF_functions.hpp"
#include "share/core/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_workspace.hpp>

using scream::Real;
using scream::Int;

// A C++ interface to stengelF fortran calls and vice versa

extern "C" {
    void set_log_file_name_f90_c(const char** fname);

    void stengelF_eamxx_bridge_init_c(Int pcols, Int pver );

    void stengelF_eamxx_bridge_run_c(Int pcols, Real* p_mid, Real* T_mid);
} // extern "C" : end _c decls

namespace scream {
    namespace stengelF {

    void set_log_file_name_f90(const char** logname){
        set_log_file_name_f90_c(logname);
    }

    void stengelF_eamxx_bridge_init( Int pcols, Int pver ){
        stengelF_eamxx_bridge_init_c( pcols, pver );
    }

    void stengelF_eamxx_bridge_run( Int pcols, Int pver, StengelFFunc::params &params){ 
        //----------------------------------------------------------------------------
        // Need to transpose to match how Fortran handles things
        params.transpose<ekat::TransposeDirection::c2f>(pcols,pver);

        stengelF_eamxx_bridge_run_c(pcols, params.f_p_mid.data(), params.f_T_mid.data()); 

        // Transpose back to C++ convention
        params.transpose<ekat::TransposeDirection::f2c>(pcols,pver);

        //----------------------------------------------------------------------------
        }

    // end _c impls

    } // namespace stengelF
} // namespace scream