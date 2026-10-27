# =============================================================================
# 02_exposure.R -- annual exposure metrics (Step 3) and lag panels (Step 4)
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tibble); library(purrr); library(rlang) })

# i_max, i_mean: anomalies relative to the seasonal climatology (seas)
# i_cum: integrated exceedance above the threshold (thresh), floored at 0
# This mirrors heatwaveR's distinction between intensity and exceedance.
annual_exposure <- function(daily) {
  stopifnot(all(c("t", "temp", "seas", "thresh") %in% names(daily)))
  daily %>%
    mutate(year = as.integer(format(t, "%Y")), doy = as.integer(format(t, "%j")),
           anom = temp - seas, excess = pmax(temp - thresh, 0), over = temp > thresh) %>%
    group_by(year) %>%
    summarise(i_max = max(anom, na.rm = TRUE),
              i_mean = mean(pmax(anom, 0), na.rm = TRUE),
              i_cum = sum(excess, na.rm = TRUE),
              mhw_days = sum(over, na.rm = TRUE),
              n_events = sum(diff(c(0L, as.integer(over))) == 1L, na.rm = TRUE),
              onset_doy = if (any(over, na.rm = TRUE)) doy[which.max(over)] else NA_integer_,
              .groups = "drop")
}

# Lags filled with 0 (no prior exposure recorded) rather than NA so the first
# years remain usable and naive vs lagged comparisons stay like-for-like.
build_lag_panel <- function(d, lags = 0:3, cols = "i_cum") {
  stopifnot(all(c("site", "year") %in% names(d)), all(cols %in% names(d)))
  d <- d %>% arrange(site, year)
  for (col in cols) for (L in lags) {
    nm <- sprintf("%s_lag%d", col, L)
    d <- d %>% group_by(site) %>%
      mutate(!!nm := dplyr::lag(.data[[col]], n = L, default = 0)) %>% ungroup()
  }
  d
}

panel_exposure <- function(panel, ...) {
  stopifnot(all(c("t", "site", "X", "Y", "temp") %in% names(panel)))
  panel %>% group_split(site) %>%
    map_dfr(function(sdf) {
      res <- detect_mhw(sdf[, c("t", "temp")], ...)
      annual_exposure(res$daily) %>%
        mutate(site = sdf$site[1], X = sdf$X[1], Y = sdf$Y[1])
    })
}

exposure_diagnostics <- function(ex) {
  q <- quantile(ex$i_cum, c(0, .25, .5, .75, .95, 1), na.rm = TRUE)
  cat(sprintf("i_cum: n=%d  min=%.2f p25=%.2f med=%.2f p75=%.2f p95=%.2f max=%.2f mean=%.2f zero=%.1f%%\n",
              nrow(ex), q[1], q[2], q[3], q[4], q[5], q[6], mean(ex$i_cum),
              100 * mean(ex$i_cum == 0)))
  invisible(q)
}
