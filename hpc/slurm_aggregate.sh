#!/bin/bash
# chain:  AID=$(sbatch --parsable hpc/slurm_array.sh); sbatch --dependency=afterok:$AID hpc/slurm_aggregate.sh
#SBATCH --job-name=mhw-agg
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/aggregate_%j.out
set -euo pipefail
module load R/4.3.1 2>/dev/null || true
export MHW_ROOT="${SLURM_SUBMIT_DIR:-$PWD}" OMP_NUM_THREADS=1
cd "$MHW_ROOT"; Rscript scripts/aggregate.R; Rscript scripts/make_figures.R
