#!/usr/bin/env Rscript
# verify_detection.R -- RUN FIRST. Regenerates detection reference values on
# this machine with this heatwaveR version and writes
# results/tables/detection_reference.csv, which the test suite checks against.
suppressPackageStartupMessages(library(dplyr))
source(file.path(Sys.getenv("MHW_ROOT", unset = "."), "R", "load_all.R"))
cat("heatwaveR", as.character(packageVersion("heatwaveR")), "|", R.version.string, "\n\n")
ref <- list(); sens <- list()
for (region in CFG$regions) {
  sst <- load_real_sst(region); res <- detect_mhw(sst)
  cat(sprintf("--- %s: %d days %s..%s, %d events ---\n", region, nrow(sst),
              format(min(sst$t)), format(max(sst$t)), nrow(res$events)))
  cat("top 3 by cumulative intensity:\n"); print(as.data.frame(top_events(res$events, "intensity_cumulative", 3)), row.names = FALSE)
  cat("top 3 by maximum intensity:\n");    print(as.data.frame(top_events(res$events, "intensity_max", 3)), row.names = FALSE)
  ex <- annual_exposure(res$daily); exposure_diagnostics(ex); cat("\n")
  te <- top_events(res$events, "intensity_max", 1)
  ref[[region]] <- data.frame(region = region, n_events = nrow(res$events), top_by = "intensity_max",
    date_start = as.character(te$date_start), date_end = as.character(te$date_end),
    duration = te$duration, intensity_max = round(te$intensity_max, 4),
    intensity_mean = round(te$intensity_mean, 4), intensity_cumulative = round(te$intensity_cumulative, 4),
    mean_i_cum = round(mean(ex$i_cum), 4), heatwaveR_version = as.character(packageVersion("heatwaveR")),
    R_version = paste(R.version$major, R.version$minor, sep = "."), generated = as.character(Sys.Date()),
    stringsAsFactors = FALSE)
  if (region == "WA") {
    sens[["WA"]] <- detect_mhw_sensitivity(sst) %>%
      filter(format(date_start, "%Y") == "2011") %>% group_by(baseline, pctile) %>%
      arrange(desc(intensity_max)) %>% slice(1) %>% ungroup() %>%
      dplyr::select(baseline, pctile, date_start, date_end, duration, intensity_max, intensity_mean, intensity_cumulative)
  }
}
out <- do.call(rbind, ref)
utils::write.csv(out, file.path(CFG$dir_tables, "detection_reference.csv"), row.names = FALSE)
utils::write.csv(sens[["WA"]], file.path(CFG$dir_tables, "detection_sensitivity_WA2011.csv"), row.names = FALSE)
cat("=== REFERENCE (written to results/tables/detection_reference.csv) ===\n")
print(out[, c("region", "date_start", "date_end", "duration", "intensity_max", "intensity_mean", "intensity_cumulative")], row.names = FALSE)
cat("\n=== WA 2011 event under alternative baselines/percentiles ===\n")
print(as.data.frame(sens[["WA"]]), row.names = FALSE, digits = 4)
