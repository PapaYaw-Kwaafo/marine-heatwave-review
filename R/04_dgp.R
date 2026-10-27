# =============================================================================
# 04_dgp.R -- data generating process for the biological response
# =============================================================================
# eta_{s,y} = alpha + beta1 * z_{s,y-L} + beta2 * z_{s,y-L}^2 + e_{s,y}
# e_{s,.}  ~ AR(1)(rho) * sigma, independent across sites
# y_{s,y}  ~ NegBin(mu = exp(eta), size = k)
# z = lagged cumulative intensity i_cum on its raw C-day scale.
#
# beta2 is small but non-zero on purpose: it produces the curvature that
# linear-in-exposure models cannot absorb. Run scripts/diagnose_dgp.R before
# any production run; a degenerate (all-zero) response is reported, not fitted.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tibble) })

simulate_response <- function(ex, beta1, lag = 1L,
                              beta2 = CFG$dgp$beta2, alpha = CFG$dgp$alpha,
                              nb_size = CFG$dgp$nb_size, rho = CFG$dgp$ar1_rho,
                              sigma = CFG$dgp$ar1_sigma, seed = 1L) {
  set.seed(seed)
  d <- ex %>% arrange(site, year) %>% group_by(site) %>%
    mutate(i_cum_lag = dplyr::lag(i_cum, n = lag, default = 0)) %>% ungroup()
  d$z <- d$i_cum_lag
  d <- d %>% group_by(site) %>%
    mutate(e = as.numeric(stats::arima.sim(model = list(ar = rho), n = dplyr::n(),
                                           n.start = 50, sd = sqrt(1 - rho^2))) * sigma) %>%
    ungroup()
  d %>% mutate(lin = alpha + beta1 * z + beta2 * z^2 + e,
               mu  = exp(pmin(lin, 20)),
               y   = stats::rnbinom(dplyr::n(), mu = mu, size = nb_size))
}

check_response <- function(bio, max_zero_frac = 0.90, min_unique = 5L, max_mean = 1e5) {
  zf <- mean(bio$y == 0, na.rm = TRUE); nu <- length(unique(bio$y)); mn <- mean(bio$y, na.rm = TRUE)
  ok <- TRUE; reason <- "ok"
  if (zf > max_zero_frac)   { ok <- FALSE; reason <- sprintf("%.0f%% zeros", 100 * zf) }
  else if (nu < min_unique) { ok <- FALSE; reason <- sprintf("only %d unique values", nu) }
  else if (mn > max_mean)   { ok <- FALSE; reason <- sprintf("mean = %.3g", mn) }
  list(ok = ok, reason = reason,
       stats = c(zero_frac = zf, n_unique = nu, mean_y = mn,
                 max_y = max(bio$y, na.rm = TRUE), mean_mu = mean(bio$mu, na.rm = TRUE)))
}

sweep_dgp <- function(ex, seed0 = 1L) {
  out <- list(); k <- 0L
  for (b in CFG$true_beta1) for (lg in CFG$true_lag) {
    k <- k + 1L
    bio <- simulate_response(ex, beta1 = b, lag = lg, seed = response_seed(seed0, b, lg))
    ck <- check_response(bio)
    out[[k]] <- tibble(true_beta1 = b, true_lag = lg,
                       zero_frac = unname(ck$stats["zero_frac"]),
                       mean_y = unname(ck$stats["mean_y"]), max_y = unname(ck$stats["max_y"]),
                       n_unique = unname(ck$stats["n_unique"]), usable = ck$ok, reason = ck$reason)
  }
  bind_rows(out)
}
