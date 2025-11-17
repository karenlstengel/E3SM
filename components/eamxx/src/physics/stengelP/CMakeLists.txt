# List of all cpp source files 
set(STENGELP_SRCS
  eamxx_stengelP_process_interface.cpp
  # stengelP.cpp
)

# List of all hpp header files
set(STENGELP_HEADERS
  eamxx_stengelP_process_interface.hpp
  stengelP_functions.hpp
)

# Adds the library to eamxx_physics 
add_library(stengelP ${STENGELP_SRCS})
target_compile_definitions(stengelP PUBLIC EAMXX_HAS_STENGELP)

target_include_directories(stengelP PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(stengelP eamxx_physics_share scream_share)
target_compile_options(stengelP PUBLIC)

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE stengelP)
endif()

# Test functions go down here if we follow the format iirc
