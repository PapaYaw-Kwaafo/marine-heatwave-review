#!/usr/bin/env Rscript
# run_chunk.R -- one array task.
#   Rscript scripts/run_chunk.R --task 1 --ntasks 40 [--reps 10] [--models a,b] [--minutes 50] [--overwrite]
# SLURM_ARRAY_TASK_ID / _COUNT are read automatically when flags are absent.
suppressPackageStartupMessages(library(dplyr))
source(file.path(Sys.getenv("MHW_ROOT", unset = "."), "R", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default) { i <- which(args == flag); if (length(i) && length(args) > i) args[i + 1] else default }
task   <- as.integer(getarg("--task",   Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))
ntasks <- as.integer(getarg("--ntasks", Sys.getenv("SLURM_ARRAY_TASK_COUNT", "1")))
reps   <- as.integer(getarg("--reps", as.character(CFG$n_replicates)))
models <- strsplit(getarg("--models", paste(CFG$default_models, collapse = ",")), ",")[[1]]
mins   <- as.numeric(getarg("--minutes", "Inf"))
mhw_log("run_chunk.R | task %d/%d | reps %d | %s | %s", task, ntasks, reps, R.version.string, Sys.info()[["nodename"]])
print(available_models(models))
run_task(task = task, n_tasks = ntasks, models = models, n_replicates = reps,
         overwrite = "--overwrite" %in% args, max_minutes = mins)
