# CMakeLists.txt for the DCMIP2016 analytic initial conditions process.
#
# Builds the dcmip2016_ic library, which contains:
#   - dcmip2016_functions.hpp       : BaroclinicWaveFunctions<S,D> struct declaration
#   - dcmip2016_functions_impl.hpp  : template function definitions
#   - dcmip2016.cpp                 : explicit instantiation for <Real, DefaultDevice>
#   - eamxx_dcmip2016_ic_process_interface.hpp/.cpp : AtmosphereProcess subclass
#
# The EAMXX_HAS_DCMIP2016_IC compile definition is used in
# register_analytic_conditions.hpp to conditionally include and register this
# process with the AtmosphereProcessFactory.
# It is propagated to all targets that link against dcmip2016_ic.

set(DCMIP2016_IC_SRCS
  dcmip2016.cpp
  eamxx_dcmip2016_ic_process_interface.cpp
)

add_library(dcmip2016_ic ${DCMIP2016_IC_SRCS})

# Propagate the feature flag to all targets that link this library
target_compile_definitions(dcmip2016_ic PUBLIC EAMXX_HAS_DCMIP2016_IC)

target_include_directories(dcmip2016_ic PUBLIC
  # Allows #include "analytic_conditions/dcmip2016/..." from any EAMxx source directory
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# EAMxx shared infrastructure: field manager, grid, KokkosTypes, units, etc.
target_link_libraries(dcmip2016_ic eamxx_physics_share scream_share)

# Wire into the top-level eamxx_physics interface target so this library is
# automatically pulled in by any EAMxx executable that uses physics.
if (TARGET eamxx_physics)
  target_link_libraries(eamxx_physics INTERFACE dcmip2016_ic)
endif()
