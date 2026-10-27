# =============================================================================
# 07_aggregate.R -- aggregation of panel outputs into reportable summaries
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tibble); library(purrr); library(tidyr) })

# Reads every results/panels* directory, so a model added in a later pass
# (same seeds, so identical datasets) merges cleanly with the main run.
collect_panels <- function(dirs = list.dirs(CFG$dir_results, recursive = FALSE)) {
  dirs <- dirs[grepl("panels", basename(dirs))]
  files <- unlist(lapply(dirs, list.files, pattern = "\\.csv$", full.names = TRUE))
  if (length(files) == 0L) stop("no panel files found under ", CFG$dir_results)
  raw <- map_dfr(files, function(f) utils::read.csv(f, stringsAsFactors = FALSE))
  raw <- raw %>% distinct(panel_key, true_beta1, true_lag, model_id, .keep_all = TRUE)
  mhw_log("collected %d rows from %d panel files in %d directories", nrow(raw), length(files), length(dirs))
  as_tibble(raw)
}

# Coverage is computed only over replicates that returned a finite interval,
# and n_ci records how many that was.
summarise_scenarios <- function(raw) {
  raw %>% filter(!is.na(model_id)) %>%
    group_by(scenario_id, region, model_id, true_beta1, true_lag, spatial_range_km, n_sites) %>%
    summarise(n_rep = dplyr::n(), n_converged = sum(converged, na.rm = TRUE),
              n_estimated = sum(is.finite(beta1_hat)),
              bias_mean = mean(bias, na.rm = TRUE), bias_median = stats::median(bias, na.rm = TRUE),
              bias_sd = stats::sd(bias, na.rm = TRUE), rmse = sqrt(mean(sq_err, na.rm = TRUE)),
              n_ci = sum(is.finite(beta1_lo) & is.finite(beta1_hi)),
              ci_coverage = if (sum(is.finite(beta1_lo)) > 0) mean(ci_covers[is.finite(beta1_lo)], na.rm = TRUE) else NA_real_,
              log_score_mean = mean(log_score, na.rm = TRUE),
              morans_I_mean = mean(morans_I, na.rm = TRUE),
              runtime_median = stats::median(runtime_sec, na.rm = TRUE), .groups = "drop")
}

summarise_models <- function(scen) {
  scen %>% group_by(model_id) %>%
    summarise(n_scenarios = dplyr::n(), mean_bias = mean(bias_mean, na.rm = TRUE),
              median_bias = stats::median(bias_median, na.rm = TRUE),
              median_rmse = stats::median(rmse, na.rm = TRUE),
              mean_coverage = mean(ci_coverage, na.rm = TRUE),
              mean_logscore = mean(log_score_mean, na.rm = TRUE),
              mean_morans_I = mean(morans_I_mean, na.rm = TRUE),
              median_runtime = stats::median(runtime_median, na.rm = TRUE),
              pct_converged = 100 * sum(n_converged, na.rm = TRUE) / sum(n_rep, na.rm = TRUE),
              .groups = "drop") %>% arrange(median_rmse)
}

summarise_by_lag <- function(scen) {
  scen %>% group_by(model_id, true_lag) %>%
    summarise(bias_mean = mean(bias_mean, na.rm = TRUE), bias_median = stats::median(bias_median, na.rm = TRUE),
              rmse = stats::median(rmse, na.rm = TRUE), ci_coverage = mean(ci_coverage, na.rm = TRUE),
              .groups = "drop") %>% arrange(model_id, true_lag)
}

summarise_by_region <- function(scen) {
  scen %>% group_by(model_id, region) %>%
    summarise(bias_mean = mean(bias_mean, na.rm = TRUE), rmse = stats::median(rmse, na.rm = TRUE),
              ci_coverage = mean(ci_coverage, na.rm = TRUE), .groups = "drop") %>% arrange(model_id, region)
}

summarise_by_factor <- function(scen, factor_col) {
  scen %>% group_by(model_id, .data[[factor_col]]) %>%
    summarise(bias_mean = mean(bias_mean, na.rm = TRUE), rmse = stats::median(rmse, na.rm = TRUE),
              ci_coverage = mean(ci_coverage, na.rm = TRUE), .groups = "drop")
}

# Percentage bias is unstable for small true_beta1, so it is reported per level.
summarise_relative_bias <- function(scen) {
  scen %>% mutate(rel_bias_pct = 100 * bias_mean / true_beta1) %>%
    group_by(model_id, true_beta1, true_lag) %>%
    summarise(rel_bias_pct = mean(rel_bias_pct, na.rm = TRUE), .groups = "drop") %>%
    arrange(model_id, true_beta1, true_lag)
}

write_all_summaries <- function(raw, dir = CFG$dir_tables) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  scen <- summarise_scenarios(raw)
  utils::write.csv(raw,  file.path(dir, "R_raw_results.csv"), row.names = FALSE)
  utils::write.csv(scen, file.path(dir, "R_scenario_summary.csv"), row.names = FALSE)
  utils::write.csv(summarise_models(scen),        file.path(dir, "R_model_summary.csv"), row.names = FALSE)
  utils::write.csv(summarise_by_lag(scen),        file.path(dir, "R_bias_by_lag.csv"), row.names = FALSE)
  utils::write.csv(summarise_by_region(scen),     file.path(dir, "R_bias_by_region.csv"), row.names = FALSE)
  utils::write.csv(summarise_relative_bias(scen), file.path(dir, "R_relative_bias.csv"), row.names = FALSE)
  utils::write.csv(summarise_by_factor(scen, "spatial_range_km"), file.path(dir, "R_bias_by_range.csv"), row.names = FALSE)
  utils::write.csv(summarise_by_factor(scen, "n_sites"),          file.path(dir, "R_bias_by_nsites.csv"), row.names = FALSE)
  utils::write.csv(summarise_by_factor(scen, "true_beta1"),       file.path(dir, "R_bias_by_beta.csv"), row.names = FALSE)
  mhw_log("wrote 9 summary tables to %s", dir)
  invisible(scen)
}

# Read this BEFORE quoting any number: a model that failed on half the design
# still produces a plausible-looking mean bias.
run_report <- function(raw) {
  cat("================ RUN REPORT ================\n")
  cat(sprintf("rows %d | panels %d | scenarios %d\n", nrow(raw),
              dplyr::n_distinct(raw$panel_key), dplyr::n_distinct(raw$scenario_id)))
  deg <- sum(is.na(raw$model_id))
  cat(sprintf("degenerate datasets (no models fitted): %d\n", deg))
  fits <- raw %>% filter(!is.na(model_id))
  tab <- fits %>% group_by(model_id) %>%
    summarise(n_fits = dplyr::n(), pct_estimate = 100 * mean(is.finite(beta1_hat)),
              pct_interval = 100 * mean(is.finite(beta1_lo)),
              pct_converge = 100 * mean(converged, na.rm = TRUE),
              median_sec = stats::median(runtime_sec, na.rm = TRUE), .groups = "drop")
  print(as.data.frame(tab), row.names = FALSE, digits = 4)
  notes <- fits$note[!is.na(fits$note) & !fits$note %in% c("linearised smooth", "cumulative log-RR/10", "finite-difference slope")]
  if (length(notes)) { cat("notes:\n"); print(utils::head(sort(table(notes), decreasing = TRUE), 10)) }
  cat("============================================\n")
  invisible(tab)
}
