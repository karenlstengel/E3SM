#!/bin/bash
#
# Build driver_kessler.cpp against an existing EAMxx case build (the
# case's Kokkos/EKAT install and compile flags), without touching
# CMakeLists.txt.
#
# Usage:
#   ./build_driver_kessler.sh [--case-build DIR] [--src-from-case] [--out EXE]
#
#   --case-build DIR  case EXEROOT to borrow Kokkos/EKAT/flags from
#                     (default: CPP_ne16np4_1day_gpu, the failing case)
#   --src-from-case   compile against the exact Kessler source that case was
#                     built from (build/GIT_LOG commit + build/GIT_DIFF),
#                     reconstructed into ./driver_src_<case>/, instead of the
#                     working tree
#   --out EXE         output executable (default: ./driver_kessler)
#
# Run on a GPU node, e.g.:
#   qsub -I -A NTDD0004 -q main -l select=1:ncpus=4:ngpus=1 -l walltime=00:30:00
#   ./driver_kessler                   # device managed + carved + host reference
#   ./driver_kessler --poison-val 0    # does the answer depend on scratch contents?
#   ./driver_kessler --help
#
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CASE_BUILD=/glade/derecho/scratch/kstengel/E3SM/e3sm_test/JAX_v_Fortran_perf/CPP_ne16np4_1day_gpu/build
SRC_FROM_CASE=0
OUT=$HERE/driver_kessler

while [ $# -gt 0 ]; do
  case $1 in
    --case-build)    CASE_BUILD=$2; shift 2 ;;
    --src-from-case) SRC_FROM_CASE=1; shift ;;
    --out)           OUT=$2; shift 2 ;;
    -h|--help)       sed -n '2,22p' "$0"; exit 0 ;;
    *) echo "Unknown option: $1"; exit 2 ;;
  esac
done

CMB=$CASE_BUILD/cmake-bld
FLAGS_MAKE=$CMB/eamxx/src/physics/kessler/CMakeFiles/kessler.dir/flags.make
[ -f "$FLAGS_MAKE" ] || { echo "No kessler flags.make under $CMB"; exit 1; }

# Kokkos/EKAT install prefix inside the case build (compiler/mpi-specific path)
KOKKOS_LIB=$(find "$CASE_BUILD" -maxdepth 6 -name libkokkoscore.a | head -1)
[ -n "$KOKKOS_LIB" ] || { echo "libkokkoscore.a not found under $CASE_BUILD"; exit 1; }
INST=$(dirname "$(dirname "$KOKKOS_LIB")")

# --- Kessler source to compile against --------------------------------------
KSRC=$HERE
if [ $SRC_FROM_CASE -eq 1 ]; then
  REPO=$(git -C "$HERE" rev-parse --show-toplevel)
  COMMIT=$(head -1 "$CASE_BUILD/GIT_LOG" | awk '{print $1}')
  CASE_NAME=$(basename "$(dirname "$CASE_BUILD")")
  RECON=$HERE/driver_src_$CASE_NAME
  KREL=components/eamxx/src/physics/kessler
  echo "Reconstructing $KREL at $COMMIT + $CASE_BUILD/GIT_DIFF into $RECON"
  rm -rf "$RECON"
  mkdir -p "$RECON"
  for f in $(git -C "$REPO" ls-tree -r --name-only "$COMMIT" "$KREL"); do
    mkdir -p "$RECON/$(dirname "$f")"
    git -C "$REPO" show "$COMMIT:$f" > "$RECON/$f"
  done
  # Apply only the kessler part of the case's uncommitted diff
  python3 - "$CASE_BUILD/GIT_DIFF" "$RECON/kessler.patch" "$KREL" <<'EOF'
import re, sys
txt = open(sys.argv[1]).read()
parts = re.split(r'(?m)^(?=diff --git )', txt)
keep = [p for p in parts if p.startswith('diff --git a/' + sys.argv[3] + '/')]
open(sys.argv[2], 'w').write(''.join(keep))
print(len(keep), 'kessler file diffs')
EOF
  (cd "$RECON" && patch -p1 --quiet < kessler.patch)
  KSRC=$RECON/$KREL
  # Quoted includes search the including file's directory first, so the
  # driver must sit next to the reconstructed headers to pick them up.
  cp "$HERE/driver_kessler.cpp" "$KSRC/"
fi
DRIVER_CPP=$KSRC/driver_kessler.cpp
echo "Kessler source: $KSRC"

# --- Compiler and flags, taken from the kessler target -----------------------
get_var () { grep "^$1 = " "$FLAGS_MAKE" | sed "s/^$1 = //"; }
DEFINES=$(get_var CXX_DEFINES)
INCLUDES=$(get_var CXX_INCLUDES)
CXXFLAGS=$(get_var CXX_FLAGS)
HOST_CXX=$(grep '^CMAKE_CXX_COMPILER:' "$CMB/CMakeCache.txt" | cut -d= -f2)

# Kokkos installs nvcc_wrapper even in CPU-only builds, so key off whether
# this Kokkos install actually has the CUDA backend enabled.
if grep -q "^#define KOKKOS_ENABLE_CUDA\b" "$INST/include/kokkos/KokkosCore_config.h" 2>/dev/null \
   && [ -x "$INST/bin/nvcc_wrapper" ]; then
  export NVCC_WRAPPER_DEFAULT_COMPILER=$HOST_CXX
  CXX=$INST/bin/nvcc_wrapper
else
  CXX=$HOST_CXX
fi

# Put the chosen Kessler source ahead of the working-tree include paths
INCLUDES="-I$KSRC -I$KSRC/eti $INCLUDES"

LIBS="$INST/lib64/libekat_kokkosutils.a $INST/lib64/libekat_core.a $INST/lib64/libspdlog.a
      $INST/lib64/libkokkoscontainers.a $INST/lib64/libkokkossimd.a $INST/lib64/libkokkoscore.a -ldl"

echo "Compiler: $CXX"
set -x
$CXX -O2 $DEFINES $INCLUDES $CXXFLAGS -o "$OUT" "$DRIVER_CPP" $LIBS
set +x
echo "Built $OUT"
