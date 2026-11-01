#!/bin/bash
#
# Build driver_bridge_transpose.cpp against an f90 case build (its Kokkos/EKAT
# install and the kessler target's C++ flags), and optionally compile-check
# the bridge sources with the same flags.
#
# Usage (on a compute node, with the case's modules loaded):
#   ./build_driver_bridge_transpose.sh --case-build <EXEROOT> --out EXE [--check-sources]
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CASE_BUILD=""; OUT=$HERE/driver_bridge_transpose; CHECK=0
while [ $# -gt 0 ]; do
  case $1 in
    --case-build)    CASE_BUILD=$2; shift 2 ;;
    --out)           OUT=$2; shift 2 ;;
    --check-sources) CHECK=1; shift ;;
    -h|--help)       sed -n '2,9p' "$0"; exit 0 ;;
    *) echo "Unknown option: $1"; exit 2 ;;
  esac
done
CMB=$CASE_BUILD/cmake-bld
FLAGS_MAKE=$CMB/eamxx/src/physics/kessler/CMakeFiles/kessler.dir/flags.make
[ -f "$FLAGS_MAKE" ] || { echo "No kessler flags.make under $CMB"; exit 1; }
KOKKOS_LIB=$(find "$CASE_BUILD" -maxdepth 6 -name libkokkoscore.a | head -1)
INST=$(dirname "$(dirname "$KOKKOS_LIB")")
get_var () { grep "^$1 = " "$FLAGS_MAKE" | sed "s/^$1 = //"; }
DEFINES=$(get_var CXX_DEFINES); INCLUDES=$(get_var CXX_INCLUDES); CXXFLAGS=$(get_var CXX_FLAGS)
HOST_CXX=$(grep '^CMAKE_CXX_COMPILER:' "$CMB/CMakeCache.txt" | cut -d= -f2)
if grep -q "^#define KOKKOS_ENABLE_CUDA\b" "$INST/include/kokkos/KokkosCore_config.h" 2>/dev/null; then
  export NVCC_WRAPPER_DEFAULT_COMPILER=$HOST_CXX; CXX=$INST/bin/nvcc_wrapper
else
  CXX=$HOST_CXX
fi
INCLUDES="-I$HERE $INCLUDES"
if [ $CHECK = 1 ]; then
  for f in eamxx_kessler_process_interface.cpp fortran_bridge/kessler_eamxx_bridge.cpp; do
    echo "compile check: $f"
    $CXX $DEFINES $INCLUDES $CXXFLAGS -c "$HERE/$f" -o /dev/null
  done
fi
LIBS="$INST/lib64/libekat_kokkosutils.a $INST/lib64/libekat_core.a $INST/lib64/libspdlog.a
      $INST/lib64/libkokkoscontainers.a $INST/lib64/libkokkossimd.a $INST/lib64/libkokkoscore.a -ldl"
$CXX $DEFINES $INCLUDES $CXXFLAGS -o "$OUT" "$HERE/driver_bridge_transpose.cpp" $LIBS
echo "Built $OUT"
