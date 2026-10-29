#!/bin/bash
#
# Build driver_kessler_py.cpp against an existing JAX case build (the case's
# Kokkos/EKAT install, compile flags, and the Python/pybind11 it was built
# with), without touching CMakeLists.txt. Same approach as the C++ branch's
# build_driver_kessler.sh.
#
# Usage:
#   ./build_driver_kessler_py.sh [--case-build DIR] [--out EXE] [--cxxflags FLAGS]
#
#   --case-build DIR  EXEROOT of a JAX build (EAMXX_ENABLE_PYTHON ON)
#                     (default: builds/jax_gpu/build, the newest GPU build)
#   --out EXE         output executable (default: ./driver_kessler_py)
#   --cxxflags FLAGS  extra compiler flags, added after the case's flags
#
# Build and run on a compute node, with the case's modules and the JAX conda
# env active (the driver embeds that env's Python):
#   source <case>/.env_mach_specific.sh
#   module load conda && conda activate /glade/work/kstengel/conda-envs/jax-kessler
#   ./build_driver_kessler_py.sh --out ./driver_kessler_py_gpu
#   ./driver_kessler_py_gpu --grid ne16-gpu --backend host,cai --repeat 10
#   ./driver_kessler_py_gpu --help
#
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CASE_BUILD=/glade/derecho/scratch/kstengel/E3SM/e3sm_test/builds/jax_gpu/build
OUT=$HERE/driver_kessler_py
EXTRA_CXXFLAGS=""

while [ $# -gt 0 ]; do
  case $1 in
    --case-build) CASE_BUILD=$2; shift 2 ;;
    --out)        OUT=$2; shift 2 ;;
    --cxxflags)   EXTRA_CXXFLAGS=$2; shift 2 ;;
    -h|--help)    sed -n '2,22p' "$0"; exit 0 ;;
    *) echo "Unknown option: $1"; exit 2 ;;
  esac
done

CMB=$CASE_BUILD/cmake-bld
FLAGS_MAKE=$CMB/eamxx/src/physics/kessler/CMakeFiles/kessler.dir/flags.make
[ -f "$FLAGS_MAKE" ] || { echo "No kessler flags.make under $CMB"; exit 1; }
grep -q -- '-DEAMXX_HAS_PYTHON' "$FLAGS_MAKE" \
  || { echo "$CASE_BUILD was not built with EAMXX_ENABLE_PYTHON"; exit 1; }

# Kokkos/EKAT install prefix inside the case build
KOKKOS_LIB=$(find "$CASE_BUILD" -maxdepth 6 -name libkokkoscore.a | head -1)
[ -n "$KOKKOS_LIB" ] || { echo "libkokkoscore.a not found under $CASE_BUILD"; exit 1; }
INST=$(dirname "$(dirname "$KOKKOS_LIB")")

# --- Compiler and flags, taken from the kessler target -----------------------
get_var () { grep "^$1 = " "$FLAGS_MAKE" | sed "s/^$1 = //"; }
DEFINES=$(get_var CXX_DEFINES)
INCLUDES=$(get_var CXX_INCLUDES)
CXXFLAGS=$(get_var CXX_FLAGS)
HOST_CXX=$(grep '^CMAKE_CXX_COMPILER:' "$CMB/CMakeCache.txt" | cut -d= -f2)

if grep -q "^#define KOKKOS_ENABLE_CUDA\b" "$INST/include/kokkos/KokkosCore_config.h" 2>/dev/null \
   && [ -x "$INST/bin/nvcc_wrapper" ]; then
  export NVCC_WRAPPER_DEFAULT_COMPILER=$HOST_CXX
  CXX=$INST/bin/nvcc_wrapper
else
  CXX=$HOST_CXX
fi

# --- Python the case was built with (its include dir is in the flags) --------
PY_INC=$(echo "$INCLUDES" | tr ' ' '\n' | sed 's/^-isystem//; s/^-I//' | grep -E '/include/python3\.[0-9]+$' | head -1)
[ -n "$PY_INC" ] || { echo "No python include dir in $FLAGS_MAKE"; exit 1; }
PY_PREFIX=$(dirname "$(dirname "$PY_INC")")
PY_VER=$(basename "$PY_INC")          # python3.10
PY_LIBS="-L$PY_PREFIX/lib -l$PY_VER -Wl,-rpath,$PY_PREFIX/lib -lpthread -lutil -lm"

INCLUDES="-I$HERE $INCLUDES"
LIBS="$INST/lib64/libekat_kokkosutils.a $INST/lib64/libekat_core.a $INST/lib64/libspdlog.a
      $INST/lib64/libkokkoscontainers.a $INST/lib64/libkokkossimd.a $INST/lib64/libkokkoscore.a -ldl"

echo "Compiler: $CXX"
echo "Python:   $PY_PREFIX ($PY_VER)"
set -x
$CXX -O2 $DEFINES "-DKESSLER_PY_DIR=\"$HERE\"" $INCLUDES $CXXFLAGS $EXTRA_CXXFLAGS \
     -o "$OUT" "$HERE/driver_kessler_py.cpp" $LIBS $PY_LIBS
set +x
echo "Built $OUT"
