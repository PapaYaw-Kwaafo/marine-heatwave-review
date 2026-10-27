# =============================================================================
# 03_panel.R -- real OISST loading and spatial panel construction
# =============================================================================
suppressPackageStartupMessages({ library(heatwaveR); library(dplyr); library(tibble); library(purrr) })

# The three series are heatwaveR's bundled example data (Schlegel & Smit 2018),
# traceable to NOAA OISST v2.1 (Banzon et al. 2016). Anyone with heatwaveR
# installed reproduces the panel exactly.
load_real_sst <- function(region = c("WA", "NW_Atl", "Med")) {
  region <- match.arg(region)
  nm <- switch(region, WA = "sst_WA", NW_Atl = "sst_NW_Atl", Med = "sst_Med")
  e <- new.env(parent = emptyenv())
  utils::data(list = nm, package = "heatwaveR", envir = e)
  d <- get(nm, envir = e)
  tibble(t = as.Date(d$t), temp = as.numeric(d$temp))
}

make_site_grid <- function(n_sites, spacing_km = CFG$grid_spacing_km) {
  side <- as.integer(round(sqrt(n_sites)))
  if (side * side != n_sites) stop("n_sites must be a perfect square; got ", n_sites)
  expand.grid(ix = seq_len(side), iy = seq_len(side)) %>% as_tibble() %>%
    mutate(site = sprintf("S%02d", row_number()),
           X = (ix - 1) * spacing_km, Y = (iy - 1) * spacing_km) %>%
    dplyr::select(site, X, Y)
}

# Matern nu = 1/2 (exponential). range_km <= 0 gives independent sites.
matern_corr <- function(coords, range_km, nugget = 1e-6) {
  n <- nrow(coords)
  if (range_km <= 0) return(diag(n))
  D <- as.matrix(dist(coords[, c("X", "Y")]))
  exp(-D / range_km) + nugget * diag(n)
}

# Every site carries the same real OISST trajectory (true seasonality,
# autocorrelation, trend and pulse statistics preserved) plus a spatially
# correlated daily noise field. Noise is independent across days; temporal
# dependence in the response is introduced by the AR(1) term in 04_dgp.R.
build_spatial_panel <- function(sst_base, n_sites, spatial_range_km,
                                site_noise_sd = CFG$site_noise_sd, seed = 1L) {
  set.seed(seed)
  co <- make_site_grid(n_sites); nd <- nrow(sst_base)
  L  <- chol(matern_corr(co, spatial_range_km))
  Z  <- matrix(rnorm(nd * n_sites), nrow = nd, ncol = n_sites)
  noise <- (Z %*% L) * site_noise_sd
  map_dfr(seq_len(n_sites), function(s)
    tibble(t = sst_base$t, site = co$site[s], X = co$X[s], Y = co$Y[s],
           temp = sst_base$temp + noise[, s]))
}

check_panel_correlation <- function(panel, spatial_range_km) {
  wide <- panel %>% dplyr::select(t, site, temp) %>%
    tidyr::pivot_wider(names_from = site, values_from = temp)
  emp <- cor(as.matrix(wide[, -1]))
  co  <- panel %>% distinct(site, X, Y) %>% arrange(site)
  tgt <- matern_corr(co, spatial_range_km, nugget = 0)
  ut  <- upper.tri(emp)
  cat(sprintf("mean off-diagonal corr: empirical %.4f target %.4f max|diff| %.4f\n",
              mean(emp[ut]), mean(tgt[ut]), max(abs(emp[ut] - tgt[ut]))))
  invisible(list(empirical = emp, target = tgt))
}
