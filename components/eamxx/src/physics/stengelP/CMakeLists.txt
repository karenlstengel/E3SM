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

# Python setup - from the ML correction package
include(ScreamUtils)
    if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.11.0")
    message(STATUS "Downloading Pybind11")
    include(FetchContent)

    FetchContent_Declare(pybind11 GIT_REPOSITORY https://github.com/pybind/pybind11.git GIT_TAG v3.0.1)
    FetchContent_MakeAvailable(pybind11)
else()
    message(FATAL_ERROR "pybind11 is missing. Use CMake >= 3.11 or download it")
endif()
find_package(Python REQUIRED COMPONENTS Interpreter Development)

# Adds the library to eamxx_physics 
add_library(stengelP ${STENGELP_SRCS})
target_compile_definitions(stengelP PUBLIC EAMXX_HAS_STENGELP)

target_include_directories(stengelP PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(stengelP SYSTEM PUBLIC ${PYTHON_INCLUDE_DIRS})

target_link_libraries(stengelP eamxx_physics_share scream_share pybind11::pybind11 Python::Python)
target_compile_options(stengelP PUBLIC)

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE stengelP)
endif()

# Test functions go down here if we follow the format iirc
