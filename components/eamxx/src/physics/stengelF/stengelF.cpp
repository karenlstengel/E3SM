#include "share/core/eamxx_types.hpp"

namespace scream {
namespace stengelF {

/*
 * Explicit instantiation for doing main on Reals using the
 * default device.
 */
// TODO - if using this then we need to do something like waht cld_fraction does
template struct StengelFFunctions<Real,DefaultDevice>;

} // namespace stengelF
} // namespace scream
