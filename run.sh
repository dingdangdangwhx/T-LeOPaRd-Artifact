#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
runs="${RUNS:-3}"
threads="${OMP_NUM_THREADS:-8}"
jobs="${JOBS:-4}"
for value in "$runs" "$threads" "$jobs"; do
  if [[ ! "$value" =~ ^[1-9][0-9]*$ ]] || (( value > 1000 )); then
    printf 'RUNS, OMP_NUM_THREADS and JOBS must be integers in [1,1000].\n' >&2
    exit 1
  fi
done
for tool in gcc g++ make cmake python3 ar flock lscpu; do
  command -v "$tool" >/dev/null || { printf 'Missing dependency: %s\n' "$tool" >&2; exit 1; }
done
[[ "$(uname -m)" == x86_64 ]] || { printf 'Linux x86-64 is required.\n' >&2; exit 1; }
flags=" $(awk '/^flags[[:space:]]*:/ {sub(/^[^:]*:/, ""); print; exit}' /proc/cpuinfo) "
for flag in avx2 aes avx512f avx512bw avx512dq avx512vl; do
  [[ "$flags" == *" $flag "* ]] || { printf 'Required CPU feature missing: %s\n' "$flag" >&2; exit 1; }
done
exec 9>"$root/.artifact.lock"
flock -n 9 || { printf 'Another artifact run is in progress.\n' >&2; exit 1; }
out="${1:-$root/results/$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$out"
out="$(cd "$out" && pwd)"
if [[ -e "$out/system.txt" ]]; then
  printf 'Output directory already contains a run: %s\n' "$out" >&2
  exit 1
fi
export OMP_NUM_THREADS="$threads"
export LEOPARD_PUBLIC_THREADS="${LEOPARD_PUBLIC_THREADS:-$threads}"
export LEOPARD_PARALLEL_BLINDEVAL=0
export LEOPARD_PARALLEL_VERIFY_REPLIES=1
export NSERVERS=3 THRESHOLD=2 Q_SPEC=0,2 SLOTS=10
{
  uname -a
  lscpu
  free -h
  gcc --version
  g++ --version
  cmake --version
  printf 'RUNS=%s\nOMP_NUM_THREADS=%s\nLEOPARD_PUBLIC_THREADS=%s\nLEOPARD_PARALLEL_BLINDEVAL=0\n' \
    "$runs" "$threads" "$LEOPARD_PUBLIC_THREADS"
} > "$out/system.txt"
trap 'printf "Run failed; inspect logs in %s\n" "$out" >&2' ERR

printf 'Building proof backend...\n'
make -C "$root" -j"$jobs" JOBS="$jobs" lib-static > "$out/build_backend.log" 2>&1
for profile in partial72 plain72; do
  printf 'Running malicious %s (one warmup + %s measurements)...\n' "$profile" "$runs"
  make -C "$root/src/protocol" PROFILE="$profile" all > "$out/build_$profile.log" 2>&1
  bash "$root/src/protocol/benchmark.sh" > "$out/${profile}_warmup.log" 2>&1
  for ((i=1; i<=runs; i++)); do
    bash "$root/src/protocol/benchmark.sh" > "$out/${profile}_run$i.log" 2>&1
  done
done
for profile in semi_partial56 semi_plain64; do
  printf 'Running semi-honest %s (one warmup + %s measurements)...\n' "$profile" "$runs"
  make -C "$root/src/protocol" PROFILE="$profile" all > "$out/build_$profile.log" 2>&1
  "$root/src/protocol/semi-honest-bench" "$runs" > "$out/$profile.jsonl" 2> "$out/$profile.stderr"
done
python3 "$root/scripts/summarize.py" "$out" "$runs"
printf '\nResults: %s\n' "$out"
