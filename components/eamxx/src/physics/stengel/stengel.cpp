#include "share/core/eamxx_types.hpp"

namespace scream {
namespace stengel {

/*
 * Explicit instantiation for doing main on Reals using the
 * default device.
 */

template struct StengelFunctions<Real,DefaultDevice>;

} // namespace stengel
} // namespace scream
