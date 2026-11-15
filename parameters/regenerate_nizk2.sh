#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root/parameters"
for profile in partial72 plain72; do
  target="$root/third_party/nizk2/params_double_${profile}_agg.h"
  trap 'rm -f "$target.tmp"' EXIT
  sage nizk2-params-codegen.sage "${target%.h}.py" > "$target.tmp"
  mv "$target.tmp" "$target"
done
