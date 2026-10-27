# =============================================================================
# 00_config.R -- design constants, paths, seeds. Single source of truth.
# =============================================================================
# Companion to: Afrifa-Yamoah, "Marine heatwave impact assessment on marine
# fauna: a methodological synthesis and simulation based roadmap".
# Pure R. Nothing here is re-implemented from another language.
# =============================================================================

CFG <- new.env(parent = emptyenv())

CFG$root        <- Sys.getenv("MHW_ROOT", unset = normalizePath(".", mustWork = FALSE))
CFG$dir_results <- file.path(CFG$root, "results")
CFG$dir_panels  <- file.path(CFG$dir_results, Sys.getenv("MHW_PANELS_SUBDIR", "panels"))  # per-model-set passes
CFG$dir_tables  <- file.path(CFG$dir_results, "tables")
CFG$dir_figures <- file.path(CFG$dir_results, "figures")
CFG$dir_logs    <- file.path(CFG$root, "logs")
for (d in c(CFG$dir_results, CFG$dir_panels, CFG$dir_tables,
            CFG$dir_figures, CFG$dir_logs))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# -- Detection: Hobday et al. (2016) defaults, passed verbatim to heatwaveR ---
CFG$clim_start        <- "1983-01-01"
CFG$clim_end          <- "2012-12-31"
CFG$pctile            <- 90          # heatwaveR convention: 0-100
CFG$window_half_width <- 5
CFG$smooth_pct_width  <- 31
CFG$min_duration      <- 5
CFG$max_gap           <- 2

# -- Factorial design ---------------------------------------------------------
CFG$regions          <- c("WA", "NW_Atl", "Med")   # heatwaveR bundled OISST
CFG$true_beta1       <- c(-0.04, -0.08, -0.12)
CFG$true_lag         <- c(0L, 1L, 2L)
CFG$spatial_range_km <- c(0, 50)
CFG$n_sites          <- c(9L, 25L)
CFG$n_replicates     <- as.integer(Sys.getenv("MHW_NREP", unset = "10"))

# -- Data generating process --------------------------------------------------
# eta = alpha + beta1*z + beta2*z^2 + e,  e ~ AR(1)(rho)*sigma per site
# y ~ NegBin(mu = exp(eta), size = nb_size); z = lagged i_cum (C-days, raw)
CFG$dgp <- list(alpha = 2.5, beta2 = -5e-4, nb_size = 4,
                ar1_rho = 0.4, ar1_sigma = 0.3)

# -- Spatial panel ------------------------------------------------------------
CFG$grid_spacing_km <- 25
CFG$site_noise_sd   <- 0.3

# -- Model fitting ------------------------------------------------------------
CFG$holdout_years <- 3L
CFG$dlnm_max_lag  <- 3L
CFG$dlnm_var_df   <- 3L
CFG$dlnm_lag_df   <- 2L
CFG$gam_k         <- 6L
CFG$mesh_cutoff   <- 25
CFG$rf_num_trees  <- 300L

CFG$all_models     <- c("naive_glm", "glm_lag1", "gam", "dlnm", "sdmtmb", "rf")
CFG$default_models <- CFG$all_models

# -- Reproducibility ----------------------------------------------------------
CFG$master_seed <- 2026L

panel_seed <- function(region, n_sites, spatial_range_km, replicate) {
  ri <- match(region, CFG$regions)
  if (is.na(ri)) stop("unknown region: ", region)
  as.integer((CFG$master_seed + ri * 1000000L + as.integer(n_sites) * 10000L +
              as.integer(spatial_range_km) * 100L + as.integer(replicate)) %%
             .Machine$integer.max)
}

response_seed <- function(pseed, beta1, lag) {
  as.integer((pseed + 7919L + as.integer(round(beta1 * 1000)) * 13L +
              as.integer(lag) * 101L) %% .Machine$integer.max)
}

build_scenario_grid <- function() {
  g <- expand.grid(region = CFG$regions, true_beta1 = CFG$true_beta1,
                   true_lag = CFG$true_lag, spatial_range_km = CFG$spatial_range_km,
                   n_sites = CFG$n_sites, KEEP.OUT.ATTRS = FALSE,
                   stringsAsFactors = FALSE)
  g$scenario_id <- sprintf("S%03d", seq_len(nrow(g)))
  g[, c("scenario_id", "region", "true_beta1", "true_lag",
        "spatial_range_km", "n_sites")]
}

# Exposure depends only on (region, n_sites, range, replicate): 120 panels
# feed 9 (beta1 x lag) response simulations each.
build_panel_grid <- function(n_replicates = CFG$n_replicates) {
  p <- expand.grid(region = CFG$regions, n_sites = CFG$n_sites,
                   spatial_range_km = CFG$spatial_range_km,
                   replicate = seq_len(n_replicates),
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  p$panel_key <- sprintf("%s_n%02d_sp%03d_r%03d", p$region, p$n_sites,
                         p$spatial_range_km, p$replicate)
  # Replicate-major order: a time-limited run completes whole replicates
  # across every region and geometry before starting the next replicate, so
  # a partial run is a balanced subset rather than a regional one.
  p[order(p$replicate, p$region, p$n_sites, p$spatial_range_km), ]
}

describe_design <- function(n_replicates = CFG$n_replicates,
                            models = CFG$default_models) {
  g <- build_scenario_grid(); p <- build_panel_grid(n_replicates)
  nd <- nrow(p) * length(CFG$true_beta1) * length(CFG$true_lag)
  cat(sprintf("scenarios %d | panels %d | datasets %d | models %d | fits %d\n",
              nrow(g), nrow(p), nd, length(models), nd * length(models)))
  invisible(list(scenarios = g, panels = p))
}
