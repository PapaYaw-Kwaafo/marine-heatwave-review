#!/bin/bash
# slurm_array.sh -- submit: sbatch hpc/slurm_array.sh ; resume: resubmit (done panels skipped)
#SBATCH --job-name=mhw-audit
#SBATCH --array=1-40
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/audit_%A_%a.out
#SBATCH --error=logs/audit_%A_%a.err
set -euo pipefail
module load R/4.3.1 2>/dev/null || true
module load gdal geos proj udunits 2>/dev/null || true
export MHW_ROOT="${SLURM_SUBMIT_DIR:-$PWD}" MHW_NREP="${MHW_NREP:-10}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$MHW_ROOT"; mkdir -p logs results/panels results/tables results/figures
Rscript scripts/run_chunk.R --task "$SLURM_ARRAY_TASK_ID" --ntasks "$SLURM_ARRAY_TASK_COUNT" \
        --reps "$MHW_NREP" --models naive_glm,glm_lag1,gam,dlnm,sdmtmb,rf
