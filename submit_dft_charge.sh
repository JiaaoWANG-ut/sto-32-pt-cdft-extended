#!/bin/bash
#SBATCH -p lr8
#SBATCH -A lr_fstheory
#SBATCH -q fstheory_lr8_normal
#SBATCH -N 1
#SBATCH -n 128
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=128
#SBATCH --exclusive
#SBATCH -t 7-00:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err
# Unconstrained DFT charging (needs gs.gpw from pipeline).
# Usage: OPT_LABEL=pris-Knear-noOH sbatch submit_dft_charge.sh
set -eo pipefail
ROOT=/global/home/users/jiaao/workspace/STO/32-Pt-DEDL-CDFT-extended
: "${OPT_LABEL:?OPT_LABEL must be set}"
cd "$ROOT/work/$OPT_LABEL"
source "$ROOT/env.sh"
echo "=== dft_charge $OPT_LABEL $(hostname) $(date -Iseconds) MPI=${SLURM_NTASKS} ==="
if command -v srun >/dev/null 2>&1; then
  srun --mpi=pmix --cpu-bind=cores -n "${SLURM_NTASKS}" python "$ROOT/dft_charge.py"
else
  mpirun -np "${SLURM_NTASKS}" --bind-to none python "$ROOT/dft_charge.py"
fi
echo "=== done $OPT_LABEL $(date -Iseconds) ==="
