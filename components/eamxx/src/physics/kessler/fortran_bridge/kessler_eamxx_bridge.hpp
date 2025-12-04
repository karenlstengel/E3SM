#include "share/core/eamxx_types.hpp"

#include <array>
#include <utility>
#include <memory>   // for shared_ptr

#include "kessler_functions.hpp"

// Bridge functions to call fortran version of Kessler functions from C++

namespace scream {
namespace kessler {

    using KesslerFunc = kessler::KesslerMicrophysicsFunctions<Real, DefaultDevice>;

    // Glue functions to call fortran from from C++ with the Data struct
    void kessler_eamxx_bridge_init( Int pcols, Int pver );
    void kessler_eamxx_bridge_run( Int pcols, Int pver, TODO); //example: KMF:params &params
    void set_log_file_name_f90(const char** logname);

    extern "C" { // _f function decls
    }

    }  // namespace kessler
}  // namespace scream