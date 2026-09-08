include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

# Turn on the debugging logger
set(CMAKE_BUILD_TYPE TRUE)

set(EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)

if (USE_CUDA)
  include (${EKAT_MACH_FILES_PATH}/kokkos/nvidia-a100.cmake)
  include (${EKAT_MACH_FILES_PATH}/kokkos/cuda.cmake)
else()
  include (${EKAT_MACH_FILES_PATH}/kokkos/amd-zen3.cmake)
  include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
endif()
include (${EKAT_MACH_FILES_PATH}/mpi/other.cmake)

# NOTE: this must NOT be self-referential (e.g. "${EKAT_MPI_EXTRA_ARGS} ...").
# This machine file gets re-loaded (via -C) by every nested configure that
# FetchContent/ExternalProject spawns against the same CMake cache, so a
# self-appending "CACHE ... FORCE" here grows unboundedly across those
# reloads instead of just setting the value once.
if (USE_CUDA)
  set(EKAT_MPI_EXTRA_ARGS "--gpus-per-task=1" CACHE STRING "" FORCE)
endif()

#option(Kokkos_ARCH_AMPERE80 "" ON)
set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)

#message(STATUS "pm-cpu CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID} CMAKE_Fortran_COMPILER_VERSION=${CMAKE_Fortran_COMPILER_VERSION}")
# -fallow-argument-mismatch only works with gnu Fortran v10 and above; nvfortran
# and ifx don't recognize it at all, so only add it for an actual GNU Fortran build
# (independent of PROJECT_NAME/build path -- CIME vs standalone).
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" AND
    CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
  set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE)
endif()

# Set Python info
# check that this correctly finds conda install
# need to have pybind11 and mpi4py installed and Python >= 3.9.2

OPTION(EAMXX_ENABLE_PYTHON "" OFF)
if (EAMXX_ENABLE_PYTHON)
  # Sets Python_EXECUTABLE.
  if ("${CMAKE_VERSION}" VERSION_LESS "3.12.0")
    find_package(PythonInterp)
  else()
    find_package(Python COMPONENTS Interpreter Development REQUIRED)
    set(Python_EXECUTABLE ${Python_EXECUTABLE})
  endif()
endif()
message(STATUS "-- ${EAMXX_ENABLE_PYTHON} --")