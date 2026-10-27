# =============================================================================
# 01_detection.R -- marine heatwave detection via heatwaveR (canonical)
# =============================================================================
# Reference values are NEVER hard-coded here. scripts/verify_detection.R
# regenerates them on the user's machine; the test suite diffs against that
# file. An earlier Python implementation circulated figures for the 2011
# Ningaloo Nino (peak 5.27 C, 62 d, 120.17 C-days) that canonical heatwaveR
# does not reproduce (it gives peak 6.51 C, 60 d, 192.71 C-days with the same
# Hobday parameters, heatwaveR 0.5.5.9000). That history is why the only
# supported reference is one generated and inspected by the analyst.
# =============================================================================

suppressPackageStartupMessages({ library(heatwaveR); library(dplyr); library(tibble) })

detect_mhw <- function(sst_df,
                       clim_start = CFG$clim_start, clim_end = CFG$clim_end,
                       pctile = CFG$pctile,
                       window_half_width = CFG$window_half_width,
                       smooth_pct_width = CFG$smooth_pct_width,
                       min_duration = CFG$min_duration, max_gap = CFG$max_gap) {
  stopifnot(all(c("t", "temp") %in% names(sst_df)))
  if (pctile <= 1) stop("pctile must be on the 0-100 scale, got ", pctile)
  sst_df <- sst_df %>% mutate(t = as.Date(t)) %>% arrange(t)
  clim <- heatwaveR::ts2clm(data = sst_df,
                            climatologyPeriod = c(clim_start, clim_end),
                            pctile = pctile, windowHalfWidth = window_half_width,
                            smoothPercentile = TRUE,
                            smoothPercentileWidth = smooth_pct_width)
  ev <- heatwaveR::detect_event(data = clim, minDuration = min_duration,
                                joinAcrossGaps = TRUE, maxGap = max_gap)
  list(daily = as_tibble(ev$climatology), events = as_tibble(ev$event))
}

detect_mhw_sensitivity <- function(sst_df,
                                   baselines = list(c("1983-01-01", "2012-12-31"),
                                                    c("1982-01-01", "2011-12-31"),
                                                    c("1982-01-01", "2022-12-31")),
                                   pctiles = c(90, 95)) {
  out <- list(); k <- 0L
  for (bw in baselines) for (pc in pctiles) {
    k <- k + 1L
    out[[k]] <- detect_mhw(sst_df, clim_start = bw[1], clim_end = bw[2],
                           pctile = pc)$events %>%
      mutate(baseline = paste(bw, collapse = "_"), pctile = pc)
  }
  bind_rows(out)
}

top_events <- function(events, by = "intensity_cumulative", n = 5) {
  by <- match.arg(by, c("intensity_cumulative", "intensity_max",
                        "duration", "intensity_mean"))
  events %>% arrange(desc(.data[[by]])) %>%
    dplyr::select(any_of(c("date_start", "date_peak", "date_end", "duration",
                           "intensity_max", "intensity_mean",
                           "intensity_cumulative"))) %>% head(n)
}

event_in_year <- function(events, year, by = "intensity_max") {
  e <- events %>% filter(as.integer(format(date_start, "%Y")) == as.integer(year))
  if (nrow(e) == 0L) return(NULL)
  e %>% arrange(desc(.data[[by]])) %>% slice(1)
}
