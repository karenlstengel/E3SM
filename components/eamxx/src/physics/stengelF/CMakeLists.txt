set(STENGELF_F90_SRCS
  # ----------------------------------------------------------------------------
  # EAMxx side C++ bridge methods
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge/stengelF_eamxx_bridge.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge/stengelF_eamxx_bridge_main.F90
  # list others as needed

  # ----------------------------------------------------------------------------
  # Fortran STENGELF code - legacy (DNE)
  # ${PATH_TO_LEGACY_STENGELF}/stengelF_main.F90
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
set(STENGELF_SRCS
  eamxx_stengelF_process_interface.cpp
  # stengelF.cpp
)

# List of all hpp header files
set(STENGELF_HEADERS
  eamxx_stengelF_process_interface.hpp
  stengelF_functions.hpp
)

# Adds the library to eamxx_physics 
add_library(stengelF ${STENGELF_F90_SRCS} ${STENGELF_SRCS})
target_compile_definitions(stengelF PUBLIC EAMXX_HAS_STENGELF)

target_include_directories(stengelF PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  # ${PATH_TO_LEGACY_STENGELF}
)

target_link_libraries(stengelF eamxx_physics_share scream_share)
target_compile_options(stengelF PUBLIC)

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE stengelF)
endif()

# Test functions go down here if we follow the format iirc
