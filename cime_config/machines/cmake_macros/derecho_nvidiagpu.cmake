string(APPEND CONFIG_ARGS " --host=cray")
set(USE_CUDA "TRUE")
string(APPEND CPPDEFS " -DGPU -DMPAS_OPENACC")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_GETTIMEOFDAY")
endif()
string(APPEND CPPDEFS " -DTHRUST_IGNORE_CUB_VERSION_CHECK")
# string(APPEND CMAKE_C_FLAGS " -noacc")
string(APPEND CMAKE_CUDA_FLAGS " -ccbin CC -O2 -arch sm_80 --use_fast_math")
# Match the Fortran/OpenACC Kessler bridge (see FFLAGS in
# llm-fortran-modernization/laft-kessler/data/exp2_jd/src/Makefile, and the
# -Mnofma below in the OPENACC_GPU_OFFLOAD block), which disables FMA
# specifically to keep results reproducible/comparable across the Fortran,
# JAX, and C++ ports of Kessler. This applies to CMAKE_CXX_FLAGS (used for
# the Kokkos-based C++ Kessler port, compiled via nvcc_wrapper as C++), not
# CMAKE_CUDA_FLAGS above -- confirmed via the build log that CMAKE_CUDA_FLAGS
# does not actually reach the Kessler .cpp compile in this configuration.
string(APPEND CMAKE_CXX_FLAGS " -Mnofma")
string(APPEND KOKKOS_OPTIONS " -DKokkos_ARCH_AMPERE80=On -DKokkos_ENABLE_CUDA=On -DKokkos_ENABLE_CUDA_LAMBDA=On -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=Off -DKokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=Off")
set(CMAKE_CUDA_ARCHITECTURES "80")
if (OPENACC_GPU_OFFLOAD)
  # string(APPEND CMAKE_EXE_LINKER_FLAGS=" -noacc")
  # string(APPEND CMAKE_Fortran_FLAGS " -noacc")
  # string(APPEND CMAKE_EXE_LINKER_FLAGS " -noacc")
  # string(APPEND CMAKE_EXE_LINKER_FLAGS="-acc -gpu=cc80")
  set(EAMXX_ENABLE_OPENACC TRUE)
  string(APPEND CMAKE_Fortran_FLAGS " -acc -gpu=cc80 -Minfo=accel -Mnofma")
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -acc -gpu=cc80 -Minfo=accel -Mnofma")
endif()
set(HOMME_QUAD_PREC FALSE CACHE BOOL "" FORCE) # nvidia does not seem to support QUAD
set(SCC "cc")
set(SCXX "CC")
set(SFC "ftn")