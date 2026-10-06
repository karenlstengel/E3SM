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

set(EKAT_MPI_EXTRA_ARGS "${EKAT_MPI_EXTRA_ARGS} --gpus-per-task=1" CACHE STRING "" FORCE)

#option(Kokkos_ARCH_AMPERE80 "" ON)
set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)

# Numerics flags shared with the StormSPEED (CESM nvhpc) builds, so the two models
# are compiled the same way: -O2 -Mnofma -Mflushz -Kieee. They have to be set here:
# CMAKE_BUILD_TYPE=TRUE above means no CMAKE_<LANG>_FLAGS_RELEASE reach EAMxx targets
# (NVHPC then defaults to -O1), and the line above replaces all CIME C++ macro flags.
# On GPU, C++ goes through nvcc_wrapper, which passes -O2 to nvcc and the -M/-K flags
# to the host compiler (nvc++); --fmad=false turns off FMA in the CUDA kernels too
# (FMA is nvcc's default for device code).
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "NVHPC")
  set(EAMXX_NUMERICS_FLAGS "-O2 -Mnofma -Mflushz -Kieee")
  if (USE_CUDA)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EAMXX_NUMERICS_FLAGS} --fmad=false" CACHE STRING "" FORCE)
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EAMXX_NUMERICS_FLAGS}" CACHE STRING "" FORCE)
  endif()
  string(APPEND CMAKE_C_FLAGS " ${EAMXX_NUMERICS_FLAGS}")
  string(APPEND CMAKE_Fortran_FLAGS " ${EAMXX_NUMERICS_FLAGS}")
endif()

#message(STATUS "pm-cpu CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID} CMAKE_Fortran_COMPILER_VERSION=${CMAKE_Fortran_COMPILER_VERSION}")
if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
      set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
    endif()
  endif()
else()
  set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
endif()