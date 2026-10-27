#!/bin/bash
# run_local_parallel.sh -- workstation: bash hpc/run_local_parallel.sh 8
# Budget ~2 GB RAM per worker; do not exceed (RAM in GB)/2 workers.
set -euo pipefail
NWORKERS="${1:-4}"
export MHW_ROOT="${MHW_ROOT:-$PWD}" MHW_NREP="${MHW_NREP:-10}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p logs results/panels results/tables results/figures
for i in $(seq 1 "$NWORKERS"); do
  Rscript scripts/run_chunk.R --task "$i" --ntasks "$NWORKERS" --reps "$MHW_NREP" \
      --models naive_glm,glm_lag1,gam,dlnm,sdmtmb,rf > "logs/local_task_${i}.log" 2>&1 &
done
wait; Rscript scripts/aggregate.R; Rscript scripts/make_figures.R
