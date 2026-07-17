# List of all cpp source files 
set(KESSLER_SRCS
  eamxx_kessler_process_interface.cpp
)

# List of all hpp header files
set(KESSLER_HEADERS
  eamxx_kessler_process_interface.hpp
  kessler_functions.hpp
)

# this points to Iris's JAX translation
set(KESSLER_JAX_DIR
  /glade/derecho/scratch/kstengel/llm-fortran-modernization/fortran2jax-kessler/_officialJAX/kessler_jax
)

add_library(kessler ${KESSLER_SRCS})

message(STATUS "-- ${EAMXX_ENABLE_PYTHON} --")
# Add the library to eamxx_physics
message(STATUS "Using Python JAX version of Kessler")
find_package(Python COMPONENTS Interpreter Development REQUIRED)
target_include_directories(kessler SYSTEM PUBLIC ${PYTHON_INCLUDE_DIRS} ${KESSLER_JAX_DIR} )


target_compile_definitions(kessler PUBLIC EAMXX_HAS_KESSLER)

target_link_libraries(kessler Kokkos::kokkos)

target_link_libraries(kessler eamxx_physics_share scream_share Python::Python) # pybind11::pybind11 

target_compile_options(kessler PUBLIC)

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE kessler)
endif()

# Test functions go down here if we follow the format iirc