set(PATH_TO_LEGACY_EAM ${SCREAM_BASE_DIR}/../eam/src)
set(PATH_TO_LEGACY_CAM ${SCREAM_BASE_DIR}/../eam/src/pyhsics/cam)

# Decide on if we should use the OpenACC version of Kessler or not. 
if (EAMXX_ENABLE_OPENACC)
  message(STATUS "Using OpenACC version of Kessler")
  set(PATH_TO_LEGACY_CAM_SIMA ${SCREAM_BASE_DIR}/../../../GPU_ports/atmospheric_physics)
else()
  message(STATUS "Using non-OpenACC version of Kessler")
  set(PATH_TO_LEGACY_CAM_SIMA ${SCREAM_BASE_DIR}/../../../atmospheric_physics)
endif()

set(PATH_TO_LEGACY_KESSLER ${PATH_TO_LEGACY_CAM_SIMA}/schemes/kessler)

# From Dockerfile in CAM-SIMA, to get CCPP-framework for building Kessler
# Make CCPP-framework available before building test
# RUN cd atmospheric_physics/test \
#     && mkdir lib \
#     && cd lib \
#     && git clone -b develop --depth 1 https://github.com/NCAR/ccpp-framework.git
set(CCPP_SRC_PATH ${PATH_TO_LEGACY_CAM_SIMA}/test/lib/ccpp-framework/src)
set(CCPP_FORTRAN_TOOLS_PATH ${PATH_TO_LEGACY_CAM_SIMA}/test/lib/ccpp-framework/scripts/fortran_tools)

set(KESSLER_F90_SRCS
  # ----------------------------------------------------------------------------
  # EAMxx side C++ bridge methods
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge/kessler_eamxx_bridge.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge/kessler_eamxx_bridge_main.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge/kessler_eamxx_bridge_update.F90
  # ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge/ccpp_kinds.F90
  # list others as needed

  # ----------------------------------------------------------------------------
  # Fortran KESSLER code - legacy from CAM-SIMA
  ${PATH_TO_LEGACY_KESSLER}/kessler.F90
  ${PATH_TO_LEGACY_KESSLER}/kessler_update.F90
  ${PATH_TO_LEGACY_CAM_SIMA}/schemes/utilities/geopotential_temp.F90
  ${PATH_TO_LEGACY_CAM_SIMA}/test/include/ccpp_kinds.F90
  ${CCPP_SRC_PATH}/ccpp_constituent_prop_mod.F90
  ${CCPP_SRC_PATH}/ccpp_hash_table.F90
  ${CCPP_SRC_PATH}/ccpp_hashable.F90
  # ${CCPP_SRC_PATH}/ccpp_types.F90
  # ----------------------------------------------------------------------------
  # shared fortran modules
  ${SCREAM_BASE_DIR}/../../share/util/shr_sys_mod.F90
  ${SCREAM_BASE_DIR}/../../share/util/shr_kind_mod.F90
  ${SCREAM_BASE_DIR}/../../share/util/shr_assert_mod.F90.in
  ${PATH_TO_LEGACY_EAM}/physics/crm/openacc_utils.F90
  #${PATH_TO_LEGACY_EAM}/utils/spmd_utils.F90

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

set_target_properties(kessler PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_compile_definitions(kessler PUBLIC EAMXX_HAS_KESSLER)

target_include_directories(kessler PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_bridge
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  ${PATH_TO_LEGACY_EAM}
  ${CCPP_SRC_PATH}
  ${CCPP_FORTRAN_TOOLS_PATH}
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