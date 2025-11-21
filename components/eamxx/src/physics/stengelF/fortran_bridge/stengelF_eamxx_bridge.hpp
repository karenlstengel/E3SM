#include "share/core/eamxx_types.hpp"

#include <array>
#include <utility>
#include <memory>   // for shared_ptr

#include "stengelF_functions.hpp"

// Bridge functions to call fortran version of StengelF functions from C++

namespace scream {
namespace stengelF {

    using StengelFFunc = stengelF::StengelFFunctions<Real, DefaultDevice>;

    // Glue functions to call fortran from from C++ with the Data struct
    void stengelF_eamxx_bridge_init( Int pcols, Int pver );
    void stengelF_eamxx_bridge_run( Int pcols, Int pver, StengelFFunc::params &params); 
    void set_log_file_name_f90(const char** logname);

    extern "C" { // _f function decls
    }

    }  // namespace stengelF
}  // namespace scream