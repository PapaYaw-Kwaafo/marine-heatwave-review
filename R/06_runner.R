# =============================================================================
# 06_runner.R -- resumable, checkpointed Monte Carlo runner
# =============================================================================
# THE PANEL IS THE UNIT OF WORK. One panel = one (region, n_sites, range,
# replicate). Its exposure is computed once (n_sites heatwaveR detections)
# and reused across all 9 (beta1 x lag) responses. Each finished panel is
# written atomically to results/panels/<key>.csv; existing files are skipped,
# so array tasks never collide and a killed job loses at most one panel.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tibble); library(purrr) })

mhw_log <- function(fmt, ...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), sprintf(fmt, ...), "\n", sep = "")
  utils::flush.console()
}

# Stride assignment spreads slow and fast panels evenly across tasks.
assign_panels <- function(task, n_tasks, n_replicates = CFG$n_replicates) {
  p <- build_panel_grid(n_replicates)
  if (n_tasks <= 1L) return(p)
  p[seq.int(from = task, to = nrow(p), by = n_tasks), , drop = FALSE]
}

run_one_panel <- function(panel_row, models = CFG$default_models,
                          scenarios = build_scenario_grid(), skip_degenerate = TRUE) {
  pr <- panel_row
  pseed <- panel_seed(pr$region, pr$n_sites, pr$spatial_range_km, pr$replicate)
  t_panel <- Sys.time()
  panel <- build_spatial_panel(load_real_sst(pr$region), pr$n_sites, pr$spatial_range_km, seed = pseed)
  ex    <- panel_exposure(panel)
  t_exp <- as.numeric(difftime(Sys.time(), t_panel, units = "secs"))
  out <- list(); k <- 0L
  for (b in CFG$true_beta1) for (lg in CFG$true_lag) {
    sid <- scenarios$scenario_id[scenarios$region == pr$region & scenarios$n_sites == pr$n_sites &
                                 scenarios$spatial_range_km == pr$spatial_range_km &
                                 scenarios$true_beta1 == b & scenarios$true_lag == lg]
    if (length(sid) != 1L) sid <- NA_character_
    bio <- simulate_response(ex, beta1 = b, lag = lg, seed = response_seed(pseed, b, lg))
    bio <- build_lag_panel(bio, lags = 0:CFG$dlnm_max_lag, cols = "i_cum")
    ck  <- check_response(bio)
    meta <- function(df) df %>% mutate(
      scenario_id = sid, panel_key = pr$panel_key, region = pr$region, replicate = pr$replicate,
      true_beta1 = b, true_lag = lg, spatial_range_km = pr$spatial_range_km, n_sites = pr$n_sites,
      zero_frac = unname(ck$stats["zero_frac"]), mean_y = unname(ck$stats["mean_y"]))
    k <- k + 1L
    if (skip_degenerate && !ck$ok) {
      out[[k]] <- meta(tibble(model_id = NA_character_, beta1_hat = NA_real_, beta1_lo = NA_real_,
                              beta1_hi = NA_real_, log_score = NA_real_, morans_I = NA_real_,
                              converged = FALSE, runtime_sec = NA_real_,
                              note = paste("degenerate response:", ck$reason))) %>%
        mutate(bias = NA_real_, sq_err = NA_real_, ci_covers = NA)
      next
    }
    out[[k]] <- meta(fit_all(bio, models = models)) %>%
      mutate(bias = beta1_hat - b, sq_err = (beta1_hat - b)^2,
             ci_covers = !is.na(beta1_lo) & !is.na(beta1_hi) & beta1_lo <= b & beta1_hi >= b)
  }
  res <- bind_rows(out)
  res$exposure_sec <- t_exp
  res$panel_sec    <- as.numeric(difftime(Sys.time(), t_panel, units = "secs"))
  res
}

run_task <- function(task = 1L, n_tasks = 1L, models = CFG$default_models,
                     n_replicates = CFG$n_replicates, overwrite = FALSE,
                     max_minutes = Inf) {
  mine <- assign_panels(task, n_tasks, n_replicates); scen <- build_scenario_grid()
  mhw_log("task %d/%d owns %d panels; models: %s", task, n_tasks, nrow(mine), paste(models, collapse = ", "))
  av <- available_models(models)
  if (any(!av$available)) mhw_log("WARNING: unavailable -> NA rows: %s", paste(av$model[!av$available], collapse = ", "))
  n_done <- 0L; t0 <- Sys.time()
  for (i in seq_len(nrow(mine))) {
    if (as.numeric(difftime(Sys.time(), t0, units = "mins")) > max_minutes) {
      mhw_log("time budget reached; stopping after %d panels", n_done); break
    }
    pr <- mine[i, ]; pfile <- file.path(CFG$dir_panels, paste0(pr$panel_key, ".csv"))
    if (file.exists(pfile) && !overwrite) next
    res <- tryCatch(run_one_panel(pr, models = models, scenarios = scen),
                    error = function(e) { mhw_log("  ERROR on %s: %s", pr$panel_key, conditionMessage(e)); NULL })
    if (is.null(res)) next
    tmp <- paste0(pfile, ".tmp"); utils::write.csv(res, tmp, row.names = FALSE); file.rename(tmp, pfile)
    n_done <- n_done + 1L
    mhw_log("  [%d/%d] %s (%d rows, %.1fs) elapsed %.1f min", i, nrow(mine), pr$panel_key,
            nrow(res), res$panel_sec[1], as.numeric(difftime(Sys.time(), t0, units = "mins")))
  }
  mhw_log("task %d: %d panels computed this call", task, n_done)
  invisible(n_done)
}

run_status <- function(n_replicates = CFG$n_replicates) {
  p <- build_panel_grid(n_replicates)
  done <- file.exists(file.path(CFG$dir_panels, paste0(p$panel_key, ".csv")))
  cat(sprintf("panels complete: %d / %d (%.1f%%)\n", sum(done), nrow(p), 100 * mean(done)))
  invisible(p[!done, , drop = FALSE])
}
