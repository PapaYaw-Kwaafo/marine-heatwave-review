#!/usr/bin/env Rscript
# aggregate.R -- run after all tasks complete.
suppressPackageStartupMessages(library(dplyr))
source(file.path(Sys.getenv("MHW_ROOT", unset = "."), "R", "load_all.R"))
missing <- run_status()
if (nrow(missing) > 0) cat("\nWARNING: design incomplete; summaries cover available panels only.\n\n")
raw <- collect_panels(); run_report(raw); scen <- write_all_summaries(raw)
cat("\n--- MODEL SUMMARY ---\n"); print(as.data.frame(summarise_models(scen)), row.names = FALSE, digits = 4)
cat("\n--- BIAS BY TRUE LAG ---\n"); print(as.data.frame(summarise_by_lag(scen)), row.names = FALSE, digits = 4)
writeLines(capture.output(sessionInfo()), file.path(CFG$dir_tables, "sessionInfo.txt"))
