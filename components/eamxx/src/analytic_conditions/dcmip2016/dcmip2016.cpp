// Explicit template instantiation for BaroclinicWaveFunctions<Real, DefaultDevice>.
//
// Including dcmip2016_functions_impl.hpp here provides the compiler with all
// template function definitions.  The explicit instantiation below creates a
// single strong definition of every member function for the <Real, DefaultDevice>
// specialization.  This avoids redundant weak instantiations across translation
// units and reduces both compile time and binary size.

#include "analytic_conditions/dcmip2016/dcmip2016_functions_impl.hpp"
#include "share/core/eamxx_types.hpp"

namespace scream {
namespace dcmip2016 {

template struct BaroclinicWaveFunctions<Real, DefaultDevice>;

} // namespace dcmip2016
} // namespace scream
