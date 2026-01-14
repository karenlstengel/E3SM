#include "share/core/eamxx_types.hpp"

#include <array>
#include <utility>
#include <memory>   // for shared_ptr

#include "kessler_functions.hpp"

// Bridge functions to call fortran version of Kessler functions from C++

namespace scream {
namespace kessler {

    using KMF = kessler::KesslerMicrophysicsFunctions<Real, DefaultDevice>;

    // Glue functions to call fortran from from C++ with the Data struct
    void kessler_eamxx_bridge_init(const Real lv_in, const Real pref_in, const Real rhoqr_in);
    void kessler_eamxx_bridge_run(Int pcols, Int pver, const double dt, Int lyr_surf, Int lyr_toa, KMF::params_in &params_in, KMF::params_out &params_out); 
    void set_log_file_name_f90(const char** logname);

    extern "C" { // _f function decls
    }

    }  // namespace kessler
}  // namespace scream