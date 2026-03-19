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
    void kessler_eamxx_bridge_init(Int pcols, Int pver, Real lv_in, Real pref_in, Real rhoqr_in, Real gravit_in);
    void kessler_eamxx_bridge_run(Int pcols, Int pver, double dt, Int lyr_surf, Int lyr_toa, KMF::params_helpers &params_helpers, KMF::params_computed &params_computed); 
    // void kessler_eamxx_bridge_update(Int pcols, Int pver, double dt, KMF::params_helpers &params_helpers, KMF::params_computed &params_computed);

    extern "C" { // _f function decls
    }

    }  // namespace kessler
}  // namespace scream