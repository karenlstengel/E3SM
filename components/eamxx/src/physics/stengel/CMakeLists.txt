# List of all cpp source files 
set(STENGEL_SRCS
  eamxx_stengel_process_interface.cpp
  stengel.cpp
)

# List of all hpp header files
set(STENGEL_HEADERS
  eamxx_stengel_process_interface.hpp
  stengel_functions.hpp
)

# Adds the library to eamxx_physics 
add_library(stengel ${STENGEL_SRCS})
target_compile_definitions(stengel PUBLIC EAMXX_HAS_STENGEL)
target_link_libraries(stengel eamxx_physics_share scream_share)
target_compile_options(stengel PUBLIC)

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE stengel)
endif()

# Test functions go down here if we follow the format iirc
