#!/usr/bin/env Rscript
# diagnose_dgp.R -- RUN SECOND. Confirms the response is non-degenerate in every design cell.
suppressPackageStartupMessages(library(dplyr))
source(file.path(Sys.getenv("MHW_ROOT", unset = "."), "R", "load_all.R"))
cat("DGP:", paste(names(CFG$dgp), unlist(CFG$dgp), sep = "=", collapse = " "), "\n\n")
bad <- FALSE; all_sw <- list()
for (region in CFG$regions) for (ns in CFG$n_sites) {
  ps <- panel_seed(region, ns, 50, 1L)
  ex <- panel_exposure(build_spatial_panel(load_real_sst(region), ns, 50, seed = ps))
  sw <- sweep_dgp(ex, seed0 = ps) %>% mutate(region = region, n_sites = ns, .before = 1)
  all_sw[[paste(region, ns)]] <- sw
  cat(sprintf("--- %s n=%d ---\n", region, ns)); exposure_diagnostics(ex)
  print(as.data.frame(sw), row.names = FALSE, digits = 3)
  if (any(!sw$usable)) bad <- TRUE
}
utils::write.csv(bind_rows(all_sw), file.path(CFG$dir_tables, "dgp_calibration.csv"), row.names = FALSE)
cat(if (bad) "\nRESULT: DEGENERATE cells exist. Revise CFG$dgp before a production run.\n"
    else "\nRESULT: all cells usable. Safe to proceed.\n")
