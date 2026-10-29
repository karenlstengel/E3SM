#!/bin/bash
#
# Bit-for-bit check of the kessler::py_run_variant variants, through
# kessler.py, with driver_kessler_py: for each grid and rain fraction, one
# driver run over every backend x variant, each a chain of model-like steps
# repeated twice. A case passes when the driver exits 0, i.e.
#   - every backend/variant matches the first (host/current) bit for bit at
#     every step of the chain,
#   - the second chain repeats the first bit for bit (determinism),
#   - every output is finite and qv/qc/qr/precl are non-negative,
# and no JAX recompiles happen after kessler.warmup() (reported, not fatal).
#
# Run on a compute node with the case's modules and the JAX conda env active
# (as for build_driver_kessler_py.sh), after building the driver:
#   ./check_kessler_py_bfb.sh --hw gpu [--exe ./driver_kessler_py_gpu] [--out DIR]
#   ./check_kessler_py_bfb.sh --hw cpu [--exe ./driver_kessler_py_cpu] [--out DIR]
#   ./check_kessler_py_bfb.sh --hw gpu --quick     # ne16 only, 8% rain
# Exit code 0 only if every case passes.
#
set -uo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
HW=""
EXE=""
OUT=./kessler_py_bfb_$(date +%Y%m%d_%H%M%S)
QUICK=0
NSTEPS=3

while [ $# -gt 0 ]; do
  case $1 in
    --hw)     HW=$2; shift 2 ;;
    --exe)    EXE=$2; shift 2 ;;
    --out)    OUT=$2; shift 2 ;;
    --nsteps) NSTEPS=$2; shift 2 ;;
    --quick)  QUICK=1; shift ;;
    -h|--help) sed -n '2,22p' "$0"; exit 0 ;;
    *) echo "Unknown option: $1"; exit 2 ;;
  esac
done

case $HW in
  gpu)
    EXE=${EXE:-$HERE/driver_kessler_py_gpu}
    GRIDS="ne16-gpu ne60-gpu"
    BACKENDS="host,device,cai"
    # map/split_map run one column at a time: correct on a GPU but ~100x
    # slower there, so they are checked on the CPU only
    VARIANTS="current,vmap,split,pallas,auto"
    PIN="" ;;
  cpu)
    EXE=${EXE:-$HERE/driver_kessler_py_cpu}
    GRIDS="ne16-cpu ne60-cpu"
    BACKENDS="host,device"
    VARIANTS="current,vmap,split,map,split_map,auto"
    PIN="taskset -c 0" ;;   # one core, as one model rank
  *) echo "--hw gpu or --hw cpu"; exit 2 ;;
esac
RAINS="0 0.03 0.08 0.12"
if [ $QUICK = 1 ]; then GRIDS=${GRIDS%% *}; RAINS="0.08"; fi
[ -x "$EXE" ] || { echo "driver not found: $EXE (build it with build_driver_kessler_py.sh)"; exit 2; }
mkdir -p "$OUT"

echo "check_kessler_py_bfb: $HW, $EXE"
echo "  backends $BACKENDS; variants $VARIANTS; $NSTEPS chained steps, twice; logs in $OUT"
printf "\n  %-9s %-5s %-6s %-22s %-14s %-8s %-7s %s\n" grid rain result "bit for bit" repeats sanity "misses" log
fail=0
for g in $GRIDS; do
  for f in $RAINS; do
    log=$OUT/bfb_${g}_rain${f}.log
    # The driver writes its perf-log CSV to the working directory
    (cd "$OUT" && $PIN "$EXE" --grid "$g" --rain-frac "$f" --backend "$BACKENDS" \
        --run-variant "$VARIANTS" --py-warmup --warmup 0 --nsteps "$NSTEPS" \
        --repeat $((2 * NSTEPS))) > "$log" 2>&1
    rc=$?
    bfb=$(grep -o "[0-9]* of [0-9]* run/step comparisons match" "$log" | sed 's/ run\/step comparisons match//')
    rep_bad=$(grep -c "REPEAT DIFFERS" "$log")
    san_bad=$(grep -c "SANITY" "$log")
    misses=$(grep "JAX traces/lowerings during timed calls" "$log" | awk '{s += $NF} END {print s + 0}')
    result=PASS
    [ $rc -ne 0 ] && { result=FAIL; fail=1; }
    printf "  %-9s %-5s %-6s %-22s %-14s %-8s %-7s %s\n" "$g" "$f" "$result" "${bfb:-?}" \
      "$([ "$rep_bad" = 0 ] && echo ok || echo "$rep_bad differ")" \
      "$([ "$san_bad" = 0 ] && echo ok || echo "$san_bad bad")" "$misses" "$(basename "$log")"
  done
done
echo
if [ $fail = 0 ]; then echo "check_kessler_py_bfb: ALL PASS"; else echo "check_kessler_py_bfb: FAILURES (see logs)"; fi
exit $fail
