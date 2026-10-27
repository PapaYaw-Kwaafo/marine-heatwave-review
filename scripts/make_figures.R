#!/usr/bin/env Rscript
# =============================================================================
# make_figures.R -- manuscript figures from results/tables/*.csv
# =============================================================================
# Conventions (do not change without changing every figure):
#   * Okabe-Ito palette; one colour per MODEL CLASS, identical in every figure.
#   * No titles inside figures; captions live in the LaTeX source.
#   * One shared legend per figure, never repeated per panel.
#   * Output: PDF (vector, for LaTeX) and PNG (preview) in results/figures/.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
source(file.path(Sys.getenv("MHW_ROOT", unset = "."), "R", "load_all.R"))

tabs <- CFG$dir_tables; figs <- CFG$dir_figures
scen <- read.csv(file.path(tabs, "R_scenario_summary.csv"), stringsAsFactors = FALSE)
raw  <- read.csv(file.path(tabs, "R_raw_results.csv"), stringsAsFactors = FALSE) %>% filter(!is.na(model_id))

# -- shared aesthetics --------------------------------------------------------
model_levels <- c("naive_glm", "glm_lag1", "gam", "dlnm", "sdmtmb", "rf")
model_labels <- c(naive_glm = "Naive NB GLM", glm_lag1 = "Lag-1 NB GLM", gam = "NB GAM",
                  dlnm = "DLNM", sdmtmb = "sdmTMB (SPDE + AR1)", rf = "Random forest")
# Okabe-Ito, fixed assignment per model class across the whole document
model_cols <- c(naive_glm = "#D55E00", glm_lag1 = "#E69F00", gam = "#009E73",
                dlnm = "#0072B2", sdmtmb = "#CC79A7", rf = "#56B4E9")
present <- model_levels[model_levels %in% unique(scen$model_id)]
scen <- scen %>% mutate(model = factor(model_id, levels = present, labels = model_labels[present]),
                        region = factor(region, levels = c("WA", "NW_Atl", "Med"),
                                        labels = c("Western Australia", "Northwest Atlantic", "Mediterranean")))
raw  <- raw  %>% mutate(model = factor(model_id, levels = present, labels = model_labels[present]))
pal  <- setNames(model_cols[present], model_labels[present])

theme_ms <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.title = element_blank(), strip.background = element_rect(fill = "grey95"))

save_fig <- function(p, name, w, h) {
  ggsave(file.path(figs, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(figs, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  mhw_log("wrote %s", name)
}

# -- Figure: bias by true lag (headline) --------------------------------------
# Distribution of per-scenario mean bias, split by true lag; the lag omission
# cost is the collapse of the naive GLM toward zero bias of the WRONG sign.
d1 <- scen %>% filter(is.finite(bias_mean)) %>%
  mutate(lag_lab = factor(true_lag, labels = c("True lag 0", "True lag 1", "True lag 2")))
p1 <- ggplot(d1, aes(x = model, y = bias_mean, fill = model)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_boxplot(outlier.size = 0.6, width = 0.65, alpha = 0.85) +
  facet_wrap(~ lag_lab, nrow = 1) +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = expression(hat(beta)[1] - beta[1])) +
  theme_ms + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_fig(p1, "fig_bias_by_lag", 8, 3.6)

# -- Figure: relative bias vs true effect size, by lag ------------------------
d2 <- scen %>% filter(is.finite(bias_mean)) %>%
  group_by(model, true_beta1, true_lag) %>%
  summarise(rel = 100 * mean(bias_mean) / first(true_beta1), .groups = "drop") %>%
  mutate(lag_lab = factor(true_lag, labels = c("True lag 0", "True lag 1", "True lag 2")))
p2 <- ggplot(d2, aes(x = factor(true_beta1), y = rel, colour = model, group = model)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  facet_wrap(~ lag_lab, nrow = 1) + scale_colour_manual(values = pal) +
  labs(x = expression(paste("True ", beta[1])), y = "Relative bias (%)") + theme_ms
save_fig(p2, "fig_relative_bias", 8, 3.6)

# -- Figure: CI coverage and RMSE ----------------------------------------------
d3 <- scen %>% filter(is.finite(ci_coverage)) %>%
  mutate(lag_lab = factor(true_lag, labels = c("Lag 0", "Lag 1", "Lag 2")))
p3a <- ggplot(d3, aes(x = model, y = ci_coverage, fill = model)) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey40") +
  geom_boxplot(outlier.size = 0.6, width = 0.65, alpha = 0.85) +
  facet_wrap(~ lag_lab, nrow = 1) + scale_fill_manual(values = pal) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Empirical 95% CI coverage") +
  theme_ms + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_fig(p3a, "fig_coverage", 8, 3.4)

d3b <- scen %>% filter(is.finite(rmse)) %>%
  mutate(lag_lab = factor(true_lag, labels = c("Lag 0", "Lag 1", "Lag 2")))
p3b <- ggplot(d3b, aes(x = model, y = rmse, fill = model)) +
  geom_boxplot(outlier.size = 0.6, width = 0.65, alpha = 0.85) +
  facet_wrap(~ lag_lab, nrow = 1) + scale_fill_manual(values = pal) +
  labs(x = NULL, y = expression(paste("RMSE of ", hat(beta)[1]))) +
  theme_ms + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_fig(p3b, "fig_rmse", 8, 3.4)

# -- Figure: cross-regional robustness ----------------------------------------
# Region x lag grid with model on the x axis: no x tick text, so nothing overlaps.
d4 <- scen %>% filter(is.finite(bias_mean)) %>%
  mutate(lag_lab = factor(true_lag, labels = c("True lag 0", "True lag 1", "True lag 2")))
p4 <- ggplot(d4, aes(x = model, y = bias_mean, fill = model)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_boxplot(outlier.size = 0.5, width = 0.7, alpha = 0.85) +
  facet_grid(lag_lab ~ region) + scale_fill_manual(values = pal) +
  labs(x = NULL, y = expression(hat(beta)[1] - beta[1])) +
  theme_ms + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_fig(p4, "fig_region", 8.5, 7.2)

# -- Figure: inference vs prediction -------------------------------------------
d5 <- scen %>% filter(is.finite(rmse), is.finite(log_score_mean)) %>%
  group_by(model) %>% summarise(rmse = median(rmse), ls = mean(log_score_mean), .groups = "drop")
if (nrow(d5) > 1) {
  p5 <- ggplot(d5, aes(x = rmse, y = ls, colour = model)) +
    geom_point(size = 3.5) + scale_colour_manual(values = pal) +
    labs(x = expression(paste("Median RMSE of ", hat(beta)[1], " (inference)")),
         y = "Mean out-of-sample log score (prediction)") + theme_ms
  save_fig(p5, "fig_inference_prediction", 5.2, 4.2)
}

# -- Figure: residual Moran's I by spatial range -------------------------------
d6 <- scen %>% filter(is.finite(morans_I_mean)) %>%
  mutate(range_lab = factor(spatial_range_km, labels = c("Independent sites", "Range 50 km")))
p6 <- ggplot(d6, aes(x = model, y = morans_I_mean, fill = model)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_boxplot(outlier.size = 0.6, width = 0.65, alpha = 0.85) +
  facet_wrap(~ range_lab, nrow = 1) + scale_fill_manual(values = pal) +
  labs(x = NULL, y = "Residual Moran's I") +
  theme_ms + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_fig(p6, "fig_morans", 7, 3.4)

# -- Figure: event detection illustration (heatwaveR style) -------------------
# Exceedance fill uses a colour that carries no other meaning in the document
# (the model palette is reserved for model classes).
EXC <- "#CC3311"
ev_panel <- function(region, label) {
  sst <- load_real_sst(region); res <- detect_mhw(sst)
  ev  <- top_events(res$events, "intensity_max", 1)
  lo  <- as.Date(ev$date_start) - 120; hi <- as.Date(ev$date_end) + 120
  win <- res$daily %>% filter(t >= lo, t <= hi) %>% mutate(region = label)
  list(win = win, ev = ev %>% mutate(region = label))
}
parts <- list(ev_panel("WA", "Western Australia"), ev_panel("NW_Atl", "Northwest Atlantic"), ev_panel("Med", "Mediterranean"))
win <- bind_rows(lapply(parts, `[[`, "win")); evs <- bind_rows(lapply(parts, `[[`, "ev"))
win$region <- factor(win$region, levels = c("Western Australia", "Northwest Atlantic", "Mediterranean"))
evs$region <- factor(evs$region, levels = levels(win$region))
evs <- evs %>% mutate(label = sprintf("%s to %s\n%d d, peak %.2f, cumulative %.1f",
                                      format(as.Date(date_start), "%d %b %Y"), format(as.Date(date_end), "%d %b %Y"),
                                      duration, intensity_max, intensity_cumulative),
                      x = as.Date(date_start) - 110)
ymax <- win %>% group_by(region) %>% summarise(y = max(temp), .groups = "drop")
evs <- evs %>% left_join(ymax, by = "region")
p7 <- ggplot(win, aes(t)) +
  geom_rect(data = evs, aes(xmin = as.Date(date_start), xmax = as.Date(date_end), ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey85", alpha = 0.6) +
  geom_ribbon(aes(ymin = thresh, ymax = pmax(temp, thresh)), fill = EXC, alpha = 0.55) +
  geom_line(aes(y = seas), colour = "grey35", linewidth = 0.45) +
  geom_line(aes(y = thresh), colour = "grey35", linetype = 2, linewidth = 0.45) +
  geom_line(aes(y = temp), colour = "black", linewidth = 0.4) +
  geom_text(data = evs, aes(x = x, y = y, label = label), hjust = 0, vjust = 1, size = 2.8, inherit.aes = FALSE) +
  facet_wrap(~ region, ncol = 1, scales = "free") +
  labs(x = NULL, y = expression(paste("Sea surface temperature (", degree, "C)"))) + theme_ms
save_fig(p7, "fig_detection_events", 7, 7.5)

# Full 41-year WA series with every detected event shaded, for context.
sst <- load_real_sst("WA"); res <- detect_mhw(sst)
p8 <- ggplot(res$daily, aes(t)) +
  geom_ribbon(aes(ymin = thresh, ymax = pmax(temp, thresh)), fill = EXC, alpha = 0.7) +
  geom_line(aes(y = thresh), colour = "grey35", linetype = 2, linewidth = 0.3) +
  geom_line(aes(y = temp), colour = "black", linewidth = 0.2) +
  labs(x = NULL, y = expression(paste("SST (", degree, "C)"))) + theme_ms
save_fig(p8, "fig_detection_wa_series", 8, 2.6)

mhw_log("all figures written to %s", figs)
