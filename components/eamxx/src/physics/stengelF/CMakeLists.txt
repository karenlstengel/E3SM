set(PATH_TO_LEGACY_EAM ${SCREAM_BASE_DIR}/../eam/src)

set(STENGELF_F90_SRCS
  # ----------------------------------------------------------------------------
  # EAMxx side C++ bridge methods
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge/stengelF_eamxx_bridge.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge/stengelF_eamxx_bridge_main.F90
  # list others as needed

  # ----------------------------------------------------------------------------
  # Fortran STENGELF code - legacy (DNE)
  # ${PATH_TO_LEGACY_STENGELF}/stengelF_main.F90
  # List all as needed 

  # ----------------------------------------------------------------------------
  # shared fortran modules
  ${SCREAM_BASE_DIR}/../../share/util/shr_sys_mod.F90
  ${SCREAM_BASE_DIR}/../../share/util/shr_kind_mod.F90
  ${SCREAM_BASE_DIR}/../../share/util/shr_assert_mod.F90.in
  ${PATH_TO_LEGACY_EAM}/physics/crm/openacc_utils.F90
  # ${PATH_TO_LEGACY_EAM}/utils/spmd_utils.F90

  # ----------------------------------------------------------------------------
  # misc other fortran modules if needed
  ${SCREAM_BASE_DIR}/../eam/src/utils/cam_abortutils.F90
  ${SCREAM_BASE_DIR}/../eam/src/control/cam_logfile.F90
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

SET(USE_OPENACC TRUE)
SET(OPENACC_Fortran_FLAGS "-acc -gpu=cc80 -Minfo=accel")
SET(OPENACC_Linker_FLAGS "-acc -gpu=cc80 -Minfo=accel")
# -acc -gpu=cc80 -Minfo=accel 

# From HOMME build
# -DOPENACC_Fortran_FLAGS="-acc -ta=tesla,pin,cc35,cuda7.5 -Minfo=accel"
# -DOPENACC_Linker_FLAGS="-acc -ta=tesla,pin,cc35,cuda7.5 -Mcuda"

SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OPENACC_Fortran_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OPENACC_Linker_FLAGS}")

message(STATUS "USE_OPENACC = ${USE_OPENACC}")
message(STATUS "OPENACC_Fortran_FLAGS = ${OPENACC_Fortran_FLAGS}")
message(STATUS "OPENACC_Linker_FLAGS = ${OPENACC_Linker_FLAGS}")

set_target_properties(stengelF PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_compile_definitions(stengelF PUBLIC EAMXX_HAS_STENGELF)

target_include_directories(stengelF PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  ${PATH_TO_LEGACY_EAM}
)

set_target_properties(stengelF PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(stengelF Kokkos::kokkos)

target_link_libraries(stengelF eamxx_physics_share scream_share)
target_compile_options(stengelF PUBLIC)

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE stengelF)
endif()

# Test functions go down here if we follow the format iirc
