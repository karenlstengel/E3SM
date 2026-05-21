#include "kessler_functions_impl.hpp"

namespace scream {
namespace kessler {

/*
 * Explicit instantiation for the default device.
 * On GPU builds without relocatable device code the full implementation
 * is included via kessler_functions.hpp; on CPU and RDC builds this
 * translation unit provides the single instantiation.
 */
template struct KesslerFunctions<Real, DefaultDevice>;

} // namespace kessler
} // namespace scream
