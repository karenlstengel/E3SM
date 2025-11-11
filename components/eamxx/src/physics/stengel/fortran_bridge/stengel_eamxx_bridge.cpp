#include "stengel_eamxx_bridge.hpp"

using scream::Real;
using scream::Int;

// A C++ interface to stengel fortran calls and vice versa

extern "C" {stengel_eamxx_bridge_run_c( required_args);
    } // extern "C" : end _c decls

    namespace scream {
    namespace stengel {

    void stengel_eamxx_bridge_init( Int pcol, Int pver ){
    stengel_eamxx_bridge_init_c( pcol, pver );
    }

    void stengel_eamxx_bridge_run( Int ncol, Int pver,
                            PackageNameF::stengel_input_state& stengel_input,
                            PackageNameF::stengel_output_tend& stengel_output,
                            PackageNameF::stengel_runtime_opt& stengel_opts
                            ){ // TODO - required_arguments = p_mid, T_mid
    //----------------------------------------------------------------------------
    // Need to transpose to match how Fortran handles things
    // stengel_input.transpose<ekat::TransposeDirection::c2f>(ncol,pver);
    // stengel_output.transpose<ekat::TransposeDirection::c2f>(ncol,pver);

    stengel_eamxx_bridge_run_c( required_arguments); // TODO - required_arguments = p_mid, T_mid

    // Transpose back to C++ convention
    // stengel_input.transpose<ekat::TransposeDirection::f2c>(ncol,pver);
    // stengel_output.transpose<ekat::TransposeDirection::f2c>(ncol,pver);

    //----------------------------------------------------------------------------
    }

    // end _c impls

    } // namespace stengel
} // namespace scream