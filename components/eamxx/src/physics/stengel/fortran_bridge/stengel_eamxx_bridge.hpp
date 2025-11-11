#include "share/core/eamxx_types.hpp"

#include <array>
#include <utility>
#include <memory>   // for shared_ptr

#include "stengel_functions.hpp"

// Bridge functions to call fortran version of Stengel functions from C++

namespace scream {
namespace stengel {

    // using ZMF = stengel::Functions<Real, DefaultDevice>;

    // Glue functions to call fortran from from C++ with the Data struct
    void stengel_eamxx_bridge_init( Int pcols, Int pver );
    void stengel_eamxx_bridge_run( Int ncol, Int pver, required_arguments); // TODO - required_arguments = p_mid, T_mid

    extern "C" { // _f function decls
    }

    }  // namespace stengel
}  // namespace scream