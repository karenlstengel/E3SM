set(PATH_TO_LEGACY_EAM ${SCREAM_BASE_DIR}/../eam/src)
set(PATH_TO_LEGACY_CAM ${SCREAM_BASE_DIR}/../eam/src/pyhsics/cam)

set(KESSLER_F90_SRCS
  # ----------------------------------------------------------------------------
  # EAMxx side C++ bridge methods
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge/kessler_eamxx_bridge.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge/kessler_eamxx_bridge_main.F90
  # list others as needed

  # ----------------------------------------------------------------------------
  # Fortran KESSLER code - legacy (DNE)
  # ${PATH_TO_LEGACY_KESSLER}/kessler_main.F90
  # List all as needed 

  # ----------------------------------------------------------------------------
  # shared fortran modules
  ${SCREAM_BASE_DIR}/../../share/util/shr_sys_mod.F90
  ${SCREAM_BASE_DIR}/../../share/util/shr_kind_mod.F90
  ${SCREAM_BASE_DIR}/../../share/util/shr_assert_mod.F90.in
  ${PATH_TO_LEGACY_EAM}/physics/crm/openacc_utils.F90
  ${PATH_TO_LEGACY_EAM}/utils/spmd_utils.F90

  # ----------------------------------------------------------------------------
  # misc other fortran modules if needed
  ${SCREAM_BASE_DIR}/../eam/src/utils/cam_abortutils.F90
  ${SCREAM_BASE_DIR}/../eam/src/control/cam_logfile.F90
  # ----------------------------------------------------------------------------
)

# List of all cpp source files 
set(KESSLER_SRCS
  eamxx_kessler_process_interface.cpp
)

# List of all hpp header files
set(KESSLER_HEADERS
  eamxx_kessler_process_interface.hpp
  kessler_functions.hpp
)

# Adds the library to eamxx_physics 
add_library(kessler ${KESSLER_F90_SRCS} ${KESSLER_SRCS})

SET(USE_OPENACC TRUE) # TODO - check that this is getting set correctly
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OPENACC_Fortran_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OPENACC_Linker_FLAGS}")

set_target_properties(kessler PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_compile_definitions(kessler PUBLIC EAMXX_HAS_KESSLER)

target_include_directories(kessler PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  ${PATH_TO_LEGACY_EAM}
)

set_target_properties(kessler PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(kessler Kokkos::kokkos)

target_link_libraries(kessler eamxx_physics_share scream_share)
target_compile_options(kessler PUBLIC)

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE kessler)
endif()

# Test functions go down here if we follow the format iirc