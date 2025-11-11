set(STENGEL_F90_SRCS
  # ----------------------------------------------------------------------------
  # EAMxx side C++ bridge methods
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge/stengel_eamxx_bridge.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge/stengel_eamxx_bridge_main.F90
  # list others as needed

  # ----------------------------------------------------------------------------
  # Fortran STENGEL code - legacy (DNE)
  # ${PATH_TO_LEGACY_STENGEL}/stengel_main.F90
  # List all needed 

  # ----------------------------------------------------------------------------
  # shared fotran modules
  ${SCREAM_BASE_DIR}/../../share/util/shr_sys_mod.F90
  ${SCREAM_BASE_DIR}/../../share/util/shr_kind_mod.F90
  ${SCREAM_BASE_DIR}/../../share/util/shr_assert_mod.F90.in

  # ----------------------------------------------------------------------------
  # misc other fortran modules if needed
  
  # ----------------------------------------------------------------------------
)


# List of all cpp source files 
set(STENGEL_SRCS
  eamxx_stengel_process_interface.cpp
  # stengel.cpp
)

# List of all hpp header files
set(STENGEL_HEADERS
  eamxx_stengel_process_interface.hpp
  stengel_functions.hpp
)

# Adds the library to eamxx_physics 
add_library(stengel ${STENGEL_F90_SRCS} ${STENGEL_SRCS})
target_compile_definitions(stengel PUBLIC EAMXX_HAS_STENGEL)

target_include_directories(stengel PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  # ${PATH_TO_LEGACY_STENGEL}
)

target_link_libraries(stengel eamxx_physics_share scream_share)
target_compile_options(stengel PUBLIC)

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE stengel)
endif()

# Test functions go down here if we follow the format iirc
