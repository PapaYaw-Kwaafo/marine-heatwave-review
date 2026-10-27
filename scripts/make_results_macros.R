#!/usr/bin/env Rscript
# =============================================================================
# make_results_macros.R -- write manuscript/results.tex from results/tables/
# =============================================================================
# Every number quoted in the manuscript is a macro defined here, computed
# from the deposited CSVs. Re-running the simulation and this script updates
# the text; nothing is typed by hand.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
source(file.path(Sys.getenv("MHW_ROOT", unset = "."), "R", "load_all.R"))
tabs <- CFG$dir_tables
out  <- file.path(CFG$root, "manuscript", "results.tex")

scen <- read.csv(file.path(tabs, "R_scenario_summary.csv"), stringsAsFactors = FALSE)
raw  <- read.csv(file.path(tabs, "R_raw_results.csv"), stringsAsFactors = FALSE)
ref  <- read.csv(file.path(tabs, "detection_reference.csv"), stringsAsFactors = FALSE)
mods <- read.csv(file.path(tabs, "R_model_summary.csv"), stringsAsFactors = FALSE)
lag  <- read.csv(file.path(tabs, "R_bias_by_lag.csv"), stringsAsFactors = FALSE)
rel  <- read.csv(file.path(tabs, "R_relative_bias.csv"), stringsAsFactors = FALSE)

fits <- raw %>% filter(!is.na(model_id))
main_models <- c("naive_glm", "glm_lag1", "gam", "dlnm", "rf")
m <- list()
m$nScenarios  <- dplyr::n_distinct(scen$scenario_id)
m$nPanels     <- nrow(build_panel_grid())
m$nReplicates <- CFG$n_replicates
m$nDatasets   <- dplyr::n_distinct(paste(raw$panel_key, raw$true_beta1, raw$true_lag))
m$nFits       <- nrow(fits)
m$nFitsMain   <- sum(fits$model_id %in% main_models)
m$waMax       <- sprintf("%.2f", ref$intensity_max[ref$region == "WA"])
m$waCum       <- sprintf("%.1f", ref$intensity_cumulative[ref$region == "WA"])
m$waDur       <- ref$duration[ref$region == "WA"]

g <- function(model, lg) lag %>% filter(model_id == model, true_lag == lg)
pct <- function(x, d = 0) formatC(100 * x, format = "f", digits = d)
num <- function(x, d = 3) formatC(x, format = "f", digits = d)

# naive GLM bias by lag (absolute, signed)
m$naiveBiasLagZero <- num(g("naive_glm", 0)$bias_mean)
m$naiveBiasLagOne  <- num(g("naive_glm", 1)$bias_mean)
m$naiveBiasLagTwo  <- num(g("naive_glm", 2)$bias_mean)
m$lagoneBiasLagZero <- num(g("glm_lag1", 0)$bias_mean)
m$lagoneBiasLagOne  <- num(g("glm_lag1", 1)$bias_mean)
m$lagoneBiasLagTwo  <- num(g("glm_lag1", 2)$bias_mean)
m$dlnmBiasLagZero <- num(g("dlnm", 0)$bias_mean)
m$dlnmBiasLagOne  <- num(g("dlnm", 1)$bias_mean)
m$dlnmBiasLagTwo  <- num(g("dlnm", 2)$bias_mean)
m$gamBiasLagOne   <- num(g("gam", 1)$bias_mean)

# relative bias of naive GLM at lag >= 1, averaged over beta levels (negative = under estimation)
rb <- rel %>% filter(model_id == "naive_glm")
m$naiveRelBiasLagOne <- formatC(abs(mean(rb$rel_bias_pct[rb$true_lag == 1])), format = "f", digits = 0)
m$naiveRelBiasLagTwo <- formatC(abs(mean(rb$rel_bias_pct[rb$true_lag == 2])), format = "f", digits = 0)
m$naiveRelBiasLagZero <- formatC(mean(rb$rel_bias_pct[rb$true_lag == 0]), format = "f", digits = 0)
rbl <- rel %>% filter(model_id == "glm_lag1")
m$lagoneRelBiasLagOne <- formatC(mean(rbl$rel_bias_pct[rbl$true_lag == 1]), format = "f", digits = 0)
m$lagoneRelBiasLagTwo <- formatC(abs(mean(rbl$rel_bias_pct[rbl$true_lag == 2])), format = "f", digits = 0)
rbd <- rel %>% filter(model_id == "dlnm")
m$dlnmRelBiasRange <- paste0(formatC(min(rbd$rel_bias_pct), format = "f", digits = 0), "\\% to ",
                             formatC(max(rbd$rel_bias_pct), format = "f", digits = 0), "\\%")

# coverage
cov <- scen %>% filter(is.finite(ci_coverage)) %>% group_by(model_id) %>%
  summarise(cov = mean(ci_coverage), .groups = "drop")
getcov <- function(mm) cov$cov[cov$model_id == mm]
lin <- c(getcov("naive_glm"), getcov("glm_lag1"))
m$linearCoverageRange <- paste0(pct(min(lin)), "\\% to ", pct(max(lin)), "\\%")
m$naiveCoverage  <- paste0(pct(getcov("naive_glm")), "\\%")
m$lagoneCoverage <- paste0(pct(getcov("glm_lag1")), "\\%")
m$gamCoverage    <- paste0(pct(getcov("gam")), "\\%")
m$dlnmCoverage   <- paste0(pct(getcov("dlnm")), "\\%")
# lag-1 GLM coverage when its assumed lag is correct
m$lagoneCoverageCorrectLag <- paste0(pct(g("glm_lag1", 1)$ci_coverage), "\\%")
m$lagoneCoverageWrongLag <- paste0(pct(mean(c(g("glm_lag1", 0)$ci_coverage, g("glm_lag1", 2)$ci_coverage))), "\\%")
m$dlnmCoverageMin <- paste0(pct(min(lag$ci_coverage[lag$model_id == "dlnm"])), "\\%")
m$dlnmCoverageMax <- paste0(pct(max(lag$ci_coverage[lag$model_id == "dlnm"])), "\\%")
m$dlnmCoverageLagOne <- paste0(pct(g("dlnm", 1)$ci_coverage), "\\%")
m$gamCoverageLagOne  <- paste0(pct(g("gam", 1)$ci_coverage), "\\%")
m$gamCoverageMin     <- paste0(pct(min(lag$ci_coverage[lag$model_id == "gam"])), "\\%")
m$naiveCoverageLagZero <- paste0(pct(g("naive_glm", 0)$ci_coverage), "\\%")
m$naiveCoverageLagOne  <- paste0(pct(g("naive_glm", 1)$ci_coverage), "\\%")
# precision: median across-replicate SD of bias at lag 1
sdv <- scen %>% filter(true_lag == 1) %>% group_by(model_id) %>% summarise(s = median(bias_sd, na.rm = TRUE), .groups = "drop")
m$dlnmSdLagOne   <- num(sdv$s[sdv$model_id == "dlnm"])
m$lagoneSdLagOne <- num(sdv$s[sdv$model_id == "glm_lag1"])
m$naiveSdLagOne  <- num(sdv$s[sdv$model_id == "naive_glm"])
m$rfSdLagOne     <- num(sdv$s[sdv$model_id == "rf"])
# RMSE by lag for dlnm vs glm_lag1
m$dlnmRmseLagOne   <- num(g("dlnm", 1)$rmse);   m$dlnmRmseLagTwo <- num(g("dlnm", 2)$rmse)
m$lagoneRmseLagOne <- num(g("glm_lag1", 1)$rmse); m$lagoneRmseLagTwo <- num(g("glm_lag1", 2)$rmse)
m$naiveRmseLagOne  <- num(g("naive_glm", 1)$rmse)
# log score by lag
lsl <- scen %>% filter(is.finite(log_score_mean)) %>% group_by(model_id, true_lag) %>% summarise(ls = mean(log_score_mean), .groups = "drop")
gl <- function(mm, lg) num(lsl$ls[lsl$model_id == mm & lsl$true_lag == lg], 2)
m$naiveLsLagZero <- gl("naive_glm", 0); m$naiveLsLagOne <- gl("naive_glm", 1); m$naiveLsLagTwo <- gl("naive_glm", 2)
m$lagoneLsLagZero <- gl("glm_lag1", 0); m$lagoneLsLagOne <- gl("glm_lag1", 1); m$lagoneLsLagTwo <- gl("glm_lag1", 2)
m$gamLsLagOne <- gl("gam", 1); m$gamLsLagTwo <- gl("gam", 2)
m$rfLsLagZero <- gl("rf", 0); m$rfLsLagOne <- gl("rf", 1); m$rfLsLagTwo <- gl("rf", 2)
m$rfBiasLagOne <- num(g("rf", 1)$bias_mean); m$rfBiasLagZero <- num(g("rf", 0)$bias_mean)
m$rfRelBiasLagOne <- formatC(mean(rel$rel_bias_pct[rel$model_id == "rf" & rel$true_lag == 1]), format = "f", digits = 0)
# relative bias of glm_lag1 when the lag is wrong (lag 0 and 2)
m$lagoneRelBiasWrongLag <- paste0(formatC(abs(mean(rbl$rel_bias_pct[rbl$true_lag == 0])), format = "f", digits = 0), "\\% and ", formatC(abs(mean(rbl$rel_bias_pct[rbl$true_lag == 2])), format = "f", digits = 0), "\\%")
m$dlnmRelBiasLagOne <- formatC(mean(rbd$rel_bias_pct[rbd$true_lag == 1]), format = "f", digits = 0)
m$dlnmRelBiasLagZero <- formatC(mean(rbd$rel_bias_pct[rbd$true_lag == 0]), format = "f", digits = 0)
m$dlnmRelBiasLagTwo  <- formatC(mean(rbd$rel_bias_pct[rbd$true_lag == 2]), format = "f", digits = 0)
m$glmPctEstimate <- formatC(100 * mean(is.finite(fits$beta1_hat[fits$model_id %in% c("naive_glm","glm_lag1")])), format = "f", digits = 1)

# RMSE
rm_ <- mods %>% filter(model_id %in% main_models)
m$dlnmRmse   <- num(rm_$median_rmse[rm_$model_id == "dlnm"])
m$naiveRmse  <- num(rm_$median_rmse[rm_$model_id == "naive_glm"])
m$lagoneRmse <- num(rm_$median_rmse[rm_$model_id == "glm_lag1"])
m$gamRmse    <- num(rm_$median_rmse[rm_$model_id == "gam"])
m$rfRmse     <- num(rm_$median_rmse[rm_$model_id == "rf"])

# prediction: log score
ls_ <- mods %>% filter(is.finite(mean_logscore))
m$bestPredModel  <- ls_$model_id[which.max(ls_$mean_logscore)]
m$bestPredScore  <- num(max(ls_$mean_logscore), 2)
m$naiveLogScore  <- num(ls_$mean_logscore[ls_$model_id == "naive_glm"], 2)
m$lagoneLogScore <- num(ls_$mean_logscore[ls_$model_id == "glm_lag1"], 2)
m$gamLogScore    <- num(ls_$mean_logscore[ls_$model_id == "gam"], 2)
m$rfLogScore     <- num(ls_$mean_logscore[ls_$model_id == "rf"], 2)

# convergence
cv <- fits %>% group_by(model_id) %>% summarise(pc = 100 * mean(converged), .groups = "drop")
m$mainConvergedMin <- formatC(min(cv$pc[cv$model_id %in% main_models]), format = "f", digits = 1)

# region invariance: rank order of |bias| at lag 1 by region
reg <- scen %>% filter(true_lag == 1, model_id %in% main_models) %>%
  group_by(region, model_id) %>% summarise(b = mean(bias_mean), .groups = "drop") %>%
  group_by(region) %>% mutate(rk = rank(abs(b))) %>% ungroup()
rk_tab <- reg %>% dplyr::select(region, model_id, rk) %>% pivot_wider(names_from = region, values_from = rk)
m$regionRankInvariant <- if (all(apply(rk_tab[, -1], 1, function(r) length(unique(r)) == 1))) "identical" else "nearly identical"
regbias <- scen %>% filter(model_id == "naive_glm", true_lag == 1) %>% group_by(region) %>%
  summarise(b = mean(bias_mean), .groups = "drop")
m$naiveBiasLagOneRegionRange <- paste0(num(min(regbias$b)), " to ", num(max(regbias$b)))

# Moran's I
mor <- scen %>% filter(is.finite(morans_I_mean), model_id %in% c("naive_glm", "glm_lag1", "gam", "dlnm")) %>%
  group_by(spatial_range_km) %>% summarise(mi = mean(morans_I_mean), .groups = "drop")
m$moransIndep <- num(mor$mi[mor$spatial_range_km == 0], 3)
m$moransCorr  <- num(mor$mi[mor$spatial_range_km == 50], 3)

# sdmTMB coverage of design
sd_ <- fits %>% filter(model_id == "sdmtmb")
if (nrow(sd_) > 0) {
  m$sdmtmbPanelsN      <- dplyr::n_distinct(sd_$panel_key)
  m$sdmtmbScenarioFrac <- paste0(pct(dplyr::n_distinct(sd_$scenario_id) / m$nScenarios), "\\%")
  m$sdmtmbFits         <- nrow(sd_)
  m$sdmtmbNonConvPct   <- formatC(100 * mean(!sd_$converged), format = "f", digits = 0)
  m$sdmtmbNonConvIndep <- formatC(100 * mean(!sd_$converged[sd_$spatial_range_km == 0]), format = "f", digits = 0)
  m$sdmtmbNonConvCorr  <- formatC(100 * mean(!sd_$converged[sd_$spatial_range_km == 50]), format = "f", digits = 0)
  sdc <- scen %>% filter(model_id == "sdmtmb")
  m$sdmtmbBiasLagOne <- num(mean(sdc$bias_mean[sdc$true_lag == 1], na.rm = TRUE))
  m$sdmtmbBiasLagZero <- num(mean(sdc$bias_mean[sdc$true_lag == 0], na.rm = TRUE))
  m$sdmtmbBiasLagTwo <- num(mean(sdc$bias_mean[sdc$true_lag == 2], na.rm = TRUE))
  m$sdmtmbCoverage <- paste0(pct(mean(sdc$ci_coverage, na.rm = TRUE)), "\\%")
  m$sdmtmbCoverageCorr <- paste0(pct(mean(sdc$ci_coverage[sdc$spatial_range_km == 50], na.rm = TRUE)), "\\%")
  m$sdmtmbCoverageConvOnly <- paste0(pct(mean(sd_$ci_covers[sd_$converged], na.rm = TRUE)), "\\%")
  m$sdmtmbRmse <- num(median(sdc$rmse, na.rm = TRUE))
  m$sdmtmbMedianSec <- formatC(median(sd_$runtime_sec), format = "f", digits = 0)
  m$sdmtmbLogScore <- num(mean(sdc$log_score_mean, na.rm = TRUE), 2)
  m$sdmtmbReplicates <- dplyr::n_distinct(sd_$replicate)
} else {
  for (k in c("sdmtmbPanelsN","sdmtmbScenarioFrac","sdmtmbFits","sdmtmbNonConvPct","sdmtmbNonConvIndep",
              "sdmtmbNonConvCorr","sdmtmbBiasLagOne","sdmtmbBiasLagZero","sdmtmbBiasLagTwo","sdmtmbCoverage",
              "sdmtmbCoverageCorr","sdmtmbCoverageConvOnly","sdmtmbRmse","sdmtmbMedianSec","sdmtmbLogScore","sdmtmbReplicates")) m[[k]] <- "NA"
}

# runtimes
rt <- fits %>% group_by(model_id) %>% summarise(s = median(runtime_sec), .groups = "drop")
m$naiveSec <- num(rt$s[rt$model_id == "naive_glm"], 2)
m$dlnmSec  <- num(rt$s[rt$model_id == "dlnm"], 2)
m$gamSec   <- num(rt$s[rt$model_id == "gam"], 2)

# session info
m$heatwaveRVersion <- as.character(packageVersion("heatwaveR"))
m$dlnmVersion      <- as.character(packageVersion("dlnm"))
m$mgcvVersion      <- as.character(packageVersion("mgcv"))
m$sdmtmbVersion    <- if (requireNamespace("sdmTMB", quietly = TRUE)) as.character(packageVersion("sdmTMB")) else "not installed"
m$rangerVersion    <- as.character(packageVersion("ranger"))
m$rVersion         <- paste(R.version$major, R.version$minor, sep = ".")

lines <- c("% results.tex -- GENERATED by scripts/make_results_macros.R. Do not edit.",
           sprintf("%% generated %s", Sys.time()))
for (k in names(m)) lines <- c(lines, sprintf("\\newcommand{\\%s}{%s}", k, m[[k]]))
writeLines(lines, out)
mhw_log("wrote %d macros to %s", length(m), out)
cat(paste(lines[-(1:2)], collapse = "\n"), "\n")
