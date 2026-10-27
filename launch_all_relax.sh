#!/bin/bash
# Submit GS relax for all 12 models (or a subset via args).
# Usage:
#   ./launch_all_relax.sh                  # all LABELS
#   ./launch_all_relax.sh pris-Knear-noOH pris-Kfar-noOH
set -eo pipefail
ROOT=/global/home/users/jiaao/workspace/STO/32-Pt-DEDL-CDFT-extended
cd "$ROOT"
# shellcheck source=/dev/null
source "$ROOT/env.sh"
if [[ $# -gt 0 ]]; then
  LABELS=("$@")
else
  mapfile -t LABELS < <(python -c "from common import LABELS; print('\n'.join(LABELS))")
fi
mkdir -p logs
for lab in "${LABELS[@]}"; do
  if [[ ! -f "work/$lab/POSCAR.init" ]]; then
    echo "SKIP missing POSCAR.init: $lab" >&2
    continue
  fi
  jid=$(OPT_LABEL="$lab" sbatch --job-name="r32-$lab" \
        --export=ALL,OPT_LABEL="$lab" \
        --output="$ROOT/logs/relax-${lab}-%j.out" \
        --error="$ROOT/logs/relax-${lab}-%j.err" \
        "$ROOT/submit.sh" | awk '{print $4}')
  echo "submitted relax $lab → $jid"
done
