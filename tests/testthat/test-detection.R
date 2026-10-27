# Checks detection against results/tables/detection_reference.csv, generated
# on THIS machine by scripts/verify_detection.R. Nothing is hard-coded: an
# earlier circular test asserted values the same code had produced.
library(testthat)
root <- Sys.getenv("MHW_ROOT", unset = "."); source(file.path(root, "R", "load_all.R"))
ref_file <- file.path(CFG$dir_tables, "detection_reference.csv")

test_that("reference file exists (run scripts/verify_detection.R first)", {
  expect_true(file.exists(ref_file))
})
test_that("bundled OISST series load with expected structure", {
  for (rg in CFG$regions) { s <- load_real_sst(rg)
    expect_true(all(c("t", "temp") %in% names(s))); expect_gt(nrow(s), 10000)
    expect_true(all(is.finite(s$temp))); expect_false(is.unsorted(s$t)) }
})
test_that("detection reproduces the stored reference", {
  skip_if_not(file.exists(ref_file))
  ref <- read.csv(ref_file, stringsAsFactors = FALSE)
  if (!identical(ref$heatwaveR_version[1], as.character(packageVersion("heatwaveR"))))
    skip("heatwaveR version changed since the reference was generated; re-run verify_detection.R")
  for (i in seq_len(nrow(ref))) { r <- ref[i, ]; res <- detect_mhw(load_real_sst(r$region))
    expect_equal(nrow(res$events), r$n_events)
    te <- top_events(res$events, "intensity_max", 1)
    expect_equal(as.character(te$date_start), r$date_start)
    expect_equal(as.integer(te$duration), as.integer(r$duration))
    expect_equal(te$intensity_max, r$intensity_max, tolerance = 1e-3)
    expect_equal(te$intensity_cumulative, r$intensity_cumulative, tolerance = 1e-2) }
})
test_that("the 2011 WA event is top by MAX, not CUMULATIVE, intensity", {
  res <- detect_mhw(load_real_sst("WA"))
  by_max <- top_events(res$events, "intensity_max", 1); by_cum <- top_events(res$events, "intensity_cumulative", 1)
  expect_equal(format(as.Date(by_max$date_start), "%Y"), "2011")
  expect_false(identical(as.character(by_max$date_start), as.character(by_cum$date_start)))
})
test_that("exposure metrics are non-negative and finite", {
  ex <- annual_exposure(detect_mhw(load_real_sst("WA"))$daily)
  expect_true(all(ex$i_cum >= 0)); expect_true(all(ex$mhw_days >= 0)); expect_true(all(is.finite(ex$i_cum)))
})
test_that("panel construction is reproducible for a fixed seed", {
  sst <- load_real_sst("WA")
  a <- build_spatial_panel(sst, 9L, 50, seed = 42L); b <- build_spatial_panel(sst, 9L, 50, seed = 42L)
  expect_equal(a$temp, b$temp); expect_equal(dplyr::n_distinct(a$site), 9L)
})
test_that("simulated response is non-degenerate at design defaults", {
  ex <- panel_exposure(build_spatial_panel(load_real_sst("WA"), 9L, 50, seed = 1L))
  expect_true(check_response(simulate_response(ex, beta1 = -0.08, lag = 1L, seed = 1L))$ok)
})
test_that("every fitter returns the unified interface", {
  ex <- panel_exposure(build_spatial_panel(load_real_sst("WA"), 9L, 0, seed = 3L))
  bio <- build_lag_panel(simulate_response(ex, beta1 = -0.08, lag = 1L, seed = 3L))
  r <- fit_all(bio, models = c("naive_glm", "glm_lag1", "gam", "dlnm", "rf"))
  expect_equal(nrow(r), 5L)
  expect_true(all(c("model_id","beta1_hat","beta1_lo","beta1_hi","log_score","morans_I","converged","runtime_sec") %in% names(r)))
})
