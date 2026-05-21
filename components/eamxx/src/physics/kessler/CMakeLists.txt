set(KESSLER_SRCS
  eamxx_kessler_process_interface.cpp
)

# Add ETI source file if not on GPU (or if GPU with relocatable device code)
if (NOT EAMXX_ENABLE_GPU OR Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE
                         OR Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
  list(APPEND KESSLER_SRCS
    eti/kessler.cpp
  )
endif()

# Create the kessler library (must come before the driver executable)
add_library(kessler ${KESSLER_SRCS})
target_compile_definitions(kessler PUBLIC EAMXX_HAS_KESSLER)
target_link_libraries(kessler eamxx_physics_share scream_share)
target_include_directories(kessler PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/eti
)

if (TARGET eamxx_physics)
  target_link_libraries(eamxx_physics INTERFACE kessler)
endif()

# Standalone driver executable
if (NOT SCREAM_LIB_ONLY)
  add_executable(driver_kessler driver_kessler.cpp)
  target_link_libraries(driver_kessler kessler scream_share)
  target_include_directories(driver_kessler PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
  )
endif()
