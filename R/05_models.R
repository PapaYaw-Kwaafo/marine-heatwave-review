# =============================================================================
# 05_models.R -- model fitters with a unified return interface
# =============================================================================
# Every fit_*() returns list(model_id, beta1_hat, beta1_lo, beta1_hi,
#   log_score, morans_I, converged, runtime_sec, note).
# Required input columns: y, site, year, X, Y, i_cum, i_cum_lag0..3, mhw_days.
#
# ESTIMAND NOTE. beta1_hat is the coefficient on lagged exposure for the GLMs;
# an OLS linearisation of the fitted smooth over the central 90% of exposure
# for the GAM; the cumulative log-RR at exposure 10 divided by 10 for the DLNM;
# the sdmTMB fixed effect; and a finite-difference marginal effect for the
# random forest. These coincide only under a linear response. State this in
# any table that compares bias across model classes.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tibble); library(MASS) })

split_train_test <- function(d, holdout_years = CFG$holdout_years) {
  ys <- sort(unique(d$year))
  if (length(ys) <= holdout_years + 5L) stop("too few years for holdout")
  te <- utils::tail(ys, holdout_years)
  list(train = dplyr::filter(d, !year %in% te), test = dplyr::filter(d, year %in% te))
}

log_score_nb <- function(y, mu, size = CFG$dgp$nb_size) {
  mu <- pmax(as.numeric(mu), 1e-8)
  ll <- stats::dnbinom(y, mu = mu, size = size, log = TRUE)
  keep <- is.finite(ll)
  if (any(keep)) mean(ll[keep]) else NA_real_
}

# Residuals averaged to site level: the question is whether neighbouring SITES
# share residual sign, not whether observations within a site do.
residual_morans_I <- function(train_df, resid_vec, k = 4L) {
  if (!requireNamespace("spdep", quietly = TRUE)) return(NA_real_)
  if (length(resid_vec) != nrow(train_df)) return(NA_real_)
  bs <- train_df %>% mutate(.r = as.numeric(resid_vec)) %>%
    group_by(site, X, Y) %>% summarise(r = mean(.r, na.rm = TRUE), .groups = "drop")
  if (nrow(bs) < 4L) return(NA_real_)
  tryCatch({
    nb <- spdep::knn2nb(spdep::knearneigh(as.matrix(bs[, c("X", "Y")]),
                                          k = min(k, nrow(bs) - 1L)))
    lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
    unname(spdep::moran.test(bs$r, lw, zero.policy = TRUE)$estimate[1])
  }, error = function(e) NA_real_)
}

elapsed  <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))
null_fit <- function(model_id, t0, note = NA_character_)
  list(model_id = model_id, beta1_hat = NA_real_, beta1_lo = NA_real_, beta1_hi = NA_real_,
       log_score = NA_real_, morans_I = NA_real_, converged = FALSE,
       runtime_sec = elapsed(t0), note = note)

# -- 1. Naive NB GLM: contemporaneous exposure --------------------------------
fit_naive_glm <- function(d) {
  t0 <- Sys.time(); s <- tryCatch(split_train_test(d), error = function(e) NULL)
  if (is.null(s)) return(null_fit("naive_glm", t0, "split failed"))
  # Warnings (typically "iteration limit reached" for theta) are muffled and
  # recorded; the fit is kept and m$converged carries the verdict.
  wm <- NA_character_
  m <- tryCatch(withCallingHandlers(
         MASS::glm.nb(y ~ i_cum + factor(site), data = s$train, control = glm.control(maxit = 100)),
         warning = function(w) { wm <<- conditionMessage(w); invokeRestart("muffleWarning") }),
       error = function(e) NULL)
  if (is.null(m)) return(null_fit("naive_glm", t0, "glm.nb failed"))
  co <- summary(m)$coefficients
  if (!"i_cum" %in% rownames(co)) return(null_fit("naive_glm", t0, "term absent"))
  b <- unname(co["i_cum", 1]); se <- unname(co["i_cum", 2])
  mu_te <- tryCatch(predict(m, newdata = s$test, type = "response"),
                    error = function(e) rep(NA_real_, nrow(s$test)))
  list(model_id = "naive_glm", beta1_hat = b, beta1_lo = b - 1.96 * se, beta1_hi = b + 1.96 * se,
       log_score = log_score_nb(s$test$y, mu_te, size = m$theta),
       morans_I = residual_morans_I(s$train, residuals(m, type = "pearson")),
       converged = isTRUE(m$converged) && is.na(wm), runtime_sec = elapsed(t0),
       note = if (is.na(wm)) NA_character_ else substr(wm, 1, 60))
}

# -- 2. Lag-1 NB GLM -----------------------------------------------------------
fit_glm_lag1 <- function(d) {
  t0 <- Sys.time(); s <- tryCatch(split_train_test(d), error = function(e) NULL)
  if (is.null(s)) return(null_fit("glm_lag1", t0, "split failed"))
  wm <- NA_character_
  m <- tryCatch(withCallingHandlers(
         MASS::glm.nb(y ~ i_cum_lag1 + factor(site), data = s$train, control = glm.control(maxit = 100)),
         warning = function(w) { wm <<- conditionMessage(w); invokeRestart("muffleWarning") }),
       error = function(e) NULL)
  if (is.null(m)) return(null_fit("glm_lag1", t0, "glm.nb failed"))
  co <- summary(m)$coefficients
  if (!"i_cum_lag1" %in% rownames(co)) return(null_fit("glm_lag1", t0, "term absent"))
  b <- unname(co["i_cum_lag1", 1]); se <- unname(co["i_cum_lag1", 2])
  mu_te <- tryCatch(predict(m, newdata = s$test, type = "response"),
                    error = function(e) rep(NA_real_, nrow(s$test)))
  list(model_id = "glm_lag1", beta1_hat = b, beta1_lo = b - 1.96 * se, beta1_hi = b + 1.96 * se,
       log_score = log_score_nb(s$test$y, mu_te, size = m$theta),
       morans_I = residual_morans_I(s$train, residuals(m, type = "pearson")),
       converged = isTRUE(m$converged) && is.na(wm), runtime_sec = elapsed(t0),
       note = if (is.na(wm)) NA_character_ else substr(wm, 1, 60))
}

# -- 3. mgcv GAM, negative binomial family ------------------------------------
fit_gam <- function(d) {
  t0 <- Sys.time()
  if (!requireNamespace("mgcv", quietly = TRUE)) return(null_fit("gam", t0, "mgcv absent"))
  s <- tryCatch(split_train_test(d), error = function(e) NULL)
  if (is.null(s)) return(null_fit("gam", t0, "split failed"))
  s$train$site <- factor(s$train$site); s$test$site <- factor(s$test$site, levels = levels(s$train$site))
  m <- tryCatch(mgcv::gam(y ~ s(i_cum_lag1, k = CFG$gam_k) + s(site, bs = "re"),
                          family = mgcv::nb(), data = s$train, method = "REML"),
                error = function(e) NULL)
  if (is.null(m)) return(null_fit("gam", t0, "gam failed"))
  rg <- stats::quantile(d$i_cum_lag1, c(0.05, 0.95), na.rm = TRUE)
  if (!all(is.finite(rg)) || diff(rg) <= 0) return(null_fit("gam", t0, "degenerate range"))
  xs <- seq(rg[1], rg[2], length.out = 40)
  nd <- data.frame(i_cum_lag1 = xs, site = levels(s$train$site)[1])
  pr <- tryCatch(mgcv::predict.gam(m, nd, type = "link", se.fit = TRUE), error = function(e) NULL)
  if (is.null(pr)) return(null_fit("gam", t0, "predict failed"))
  sl <- unname(stats::coef(stats::lm(pr$fit ~ xs))[2]); se <- stats::sd(pr$se.fit) / sqrt(length(xs))
  theta <- tryCatch(m$family$getTheta(TRUE), error = function(e) CFG$dgp$nb_size)
  mu_te <- tryCatch(as.numeric(mgcv::predict.gam(m, s$test, type = "response")),
                    error = function(e) rep(NA_real_, nrow(s$test)))
  list(model_id = "gam", beta1_hat = sl, beta1_lo = sl - 1.96 * se, beta1_hi = sl + 1.96 * se,
       log_score = log_score_nb(s$test$y, mu_te, size = theta),
       morans_I = residual_morans_I(s$train, residuals(m, type = "pearson")),
       converged = isTRUE(m$converged), runtime_sec = elapsed(t0), note = "linearised smooth")
}

# -- 4. DLNM cross-basis, lags 0..3 -------------------------------------------
fit_dlnm <- function(d, max_lag = CFG$dlnm_max_lag) {
  t0 <- Sys.time()
  if (!requireNamespace("dlnm", quietly = TRUE)) return(null_fit("dlnm", t0, "dlnm absent"))
  s <- tryCatch(split_train_test(d), error = function(e) NULL)
  if (is.null(s)) return(null_fit("dlnm", t0, "split failed"))
  lag_cols <- sprintf("i_cum_lag%d", 0:max_lag)
  if (!all(lag_cols %in% names(s$train))) return(null_fit("dlnm", t0, "lag columns missing"))
  Q <- as.matrix(s$train[, lag_cols])
  if (any(!is.finite(Q))) return(null_fit("dlnm", t0, "non-finite exposure"))
  cb <- tryCatch(dlnm::crossbasis(Q, lag = c(0, max_lag),
                                  argvar = list(fun = "ns", df = CFG$dlnm_var_df),
                                  arglag = list(fun = "ns", df = CFG$dlnm_lag_df)),
                 error = function(e) NULL)
  if (is.null(cb)) return(null_fit("dlnm", t0, "crossbasis failed"))
  m <- tryCatch(stats::glm(y ~ cb + factor(site), family = stats::quasipoisson(), data = s$train),
                error = function(e) NULL)
  if (is.null(m)) return(null_fit("dlnm", t0, "glm failed"))
  pr <- tryCatch(dlnm::crosspred(cb, m, at = c(0, 10), cen = 0, cumul = TRUE), error = function(e) NULL)
  if (is.null(pr)) return(null_fit("dlnm", t0, "crosspred failed"))
  fit10 <- tryCatch(as.numeric(pr$allRRfit["10"]), error = function(e) NA_real_)
  hi10  <- tryCatch(as.numeric(pr$allRRhigh["10"]), error = function(e) NA_real_)
  if (!is.finite(fit10) || fit10 <= 0) return(null_fit("dlnm", t0, "invalid RR"))
  sl <- log(fit10) / 10
  se <- if (is.finite(hi10) && hi10 > 0) (log(hi10) - log(fit10)) / 1.96 / 10 else NA_real_
  list(model_id = "dlnm", beta1_hat = sl,
       beta1_lo = if (is.finite(se)) sl - 1.96 * se else NA_real_,
       beta1_hi = if (is.finite(se)) sl + 1.96 * se else NA_real_,
       log_score = NA_real_,
       morans_I = residual_morans_I(s$train, residuals(m, type = "pearson")),
       converged = isTRUE(m$converged), runtime_sec = elapsed(t0), note = "cumulative log-RR/10")
}

# -- 5. sdmTMB: Matern SPDE spatial field + AR1 spatiotemporal field -----------
fit_sdmtmb <- function(d) {
  t0 <- Sys.time()
  if (!requireNamespace("sdmTMB", quietly = TRUE)) return(null_fit("sdmtmb", t0, "sdmTMB not installed"))
  s <- tryCatch(split_train_test(d), error = function(e) NULL)
  if (is.null(s)) return(null_fit("sdmtmb", t0, "split failed"))
  tr <- as.data.frame(s$train); te <- as.data.frame(s$test)
  tr$year <- as.integer(tr$year); te$year <- as.integer(te$year)
  mesh <- tryCatch(sdmTMB::make_mesh(tr, xy_cols = c("X", "Y"), cutoff = CFG$mesh_cutoff),
                   error = function(e) NULL)
  if (is.null(mesh)) return(null_fit("sdmtmb", t0, "make_mesh failed"))
  # Warnings (e.g. non-positive-definite Hessian) are muffled, not fatal: the
  # fit is retained and sanity() decides the convergence flag that the run
  # report tallies. Hiding a flagged fit as a missing row would bias the
  # convergence statistics upward.
  warn_msg <- NA_character_
  m <- tryCatch(withCallingHandlers(
         sdmTMB::sdmTMB(y ~ i_cum_lag1, data = tr, mesh = mesh, time = "year",
                        family = sdmTMB::nbinom2(link = "log"),
                        spatial = "on", spatiotemporal = "ar1", silent = TRUE),
         warning = function(w) { warn_msg <<- conditionMessage(w); invokeRestart("muffleWarning") }),
       error = function(e) NULL)
  if (is.null(m)) return(null_fit("sdmtmb", t0, "sdmTMB fit failed"))
  san  <- tryCatch(sdmTMB::sanity(m, silent = TRUE), error = function(e) NULL)
  conv <- if (is.null(san)) FALSE else isTRUE(san$all_ok)
  ti <- tryCatch(sdmTMB::tidy(m, effects = "fixed", conf.int = TRUE), error = function(e) NULL)
  if (is.null(ti) || !"i_cum_lag1" %in% ti$term) return(null_fit("sdmtmb", t0, "tidy failed"))
  rw <- ti[ti$term == "i_cum_lag1", ]
  mu_te <- tryCatch({ p <- predict(m, newdata = te, type = "response")
                      if (is.data.frame(p)) as.numeric(p$est) else as.numeric(p) },
                    error = function(e) rep(NA_real_, nrow(te)))
  phi <- tryCatch(exp(m$sd_report$par.fixed[["ln_phi"]]), error = function(e) CFG$dgp$nb_size)
  rs  <- tryCatch(as.numeric(residuals(m)), error = function(e) rep(NA_real_, nrow(tr)))
  list(model_id = "sdmtmb", beta1_hat = as.numeric(rw$estimate),
       beta1_lo = as.numeric(rw$conf.low), beta1_hi = as.numeric(rw$conf.high),
       log_score = log_score_nb(te$y, mu_te, size = phi),
       morans_I = residual_morans_I(s$train, rs), converged = conv,
       runtime_sec = elapsed(t0),
       note = if (conv) NA_character_ else paste0("sanity() flagged",
                                                  if (!is.na(warn_msg)) paste0(": ", substr(warn_msg, 1, 60)) else ""))
}

# -- 6. Random forest (ranger): forecast oriented baseline ---------------------
fit_rf <- function(d) {
  t0 <- Sys.time()
  if (!requireNamespace("ranger", quietly = TRUE)) return(null_fit("rf", t0, "ranger absent"))
  s <- tryCatch(split_train_test(d), error = function(e) NULL)
  if (is.null(s)) return(null_fit("rf", t0, "split failed"))
  m <- tryCatch(ranger::ranger(y ~ i_cum + i_cum_lag1 + mhw_days + X + Y, data = s$train,
                               num.trees = CFG$rf_num_trees, min.node.size = 5, seed = 1L,
                               num.threads = 1L),
                error = function(e) NULL)
  if (is.null(m)) return(null_fit("rf", t0, "ranger failed"))
  pr_tr <- tryCatch(predict(m, s$train)$predictions, error = function(e) NULL)
  pr_te <- tryCatch(predict(m, s$test)$predictions,  error = function(e) NULL)
  if (is.null(pr_tr) || is.null(pr_te)) return(null_fit("rf", t0, "predict failed"))
  bumped <- s$train; bumped$i_cum_lag1 <- bumped$i_cum_lag1 + 1
  pr_b <- tryCatch(predict(m, bumped)$predictions, error = function(e) NULL)
  slope <- if (is.null(pr_b)) NA_real_ else log(max(mean(pr_b), 1e-8) / max(mean(pr_tr), 1e-8))
  list(model_id = "rf", beta1_hat = slope, beta1_lo = NA_real_, beta1_hi = NA_real_,
       log_score = log_score_nb(s$test$y, pmax(pr_te, 1e-8)),
       morans_I = residual_morans_I(s$train, (s$train$y - pr_tr) / sqrt(pmax(pr_tr, 1e-8))),
       converged = TRUE, runtime_sec = elapsed(t0), note = "finite-difference slope")
}

FITTERS <- list(naive_glm = fit_naive_glm, glm_lag1 = fit_glm_lag1, gam = fit_gam,
                dlnm = fit_dlnm, sdmtmb = fit_sdmtmb, rf = fit_rf)

available_models <- function(models = CFG$all_models) {
  need <- c(naive_glm = "MASS", glm_lag1 = "MASS", gam = "mgcv", dlnm = "dlnm",
            sdmtmb = "sdmTMB", rf = "ranger")
  ok <- vapply(models, function(m) requireNamespace(need[[m]], quietly = TRUE), logical(1))
  data.frame(model = models, package = unname(need[models]), available = unname(ok), row.names = NULL)
}

fit_all <- function(d, models = CFG$default_models) {
  bind_rows(lapply(models, function(mn)
    as_tibble(tryCatch(FITTERS[[mn]](d),
                       error = function(e) null_fit(mn, Sys.time(), conditionMessage(e))))))
}
