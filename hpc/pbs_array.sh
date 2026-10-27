#!/bin/bash
# pbs_array.sh -- PBS Pro / Torque equivalent. Submit: qsub hpc/pbs_array.sh
#PBS -N mhw-audit
#PBS -J 1-40
#PBS -l select=1:ncpus=1:mem=8GB
#PBS -l walltime=04:00:00
#PBS -j oe
#PBS -o logs/
set -euo pipefail
module load R/4.3.1 2>/dev/null || true
export MHW_ROOT="${PBS_O_WORKDIR:-$PWD}" MHW_NREP="${MHW_NREP:-10}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$MHW_ROOT"; mkdir -p logs results/panels results/tables
Rscript scripts/run_chunk.R --task "$PBS_ARRAY_INDEX" --ntasks 40 --reps "$MHW_NREP" \
        --models naive_glm,glm_lag1,gam,dlnm,sdmtmb,rf
