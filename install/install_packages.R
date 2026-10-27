#!/usr/bin/env Rscript
# install_packages.R -- standard CRAN install.
#   Rscript install/install_packages.R            # core (five fast models)
#   Rscript install/install_packages.R --spatial  # + sdmTMB
#   Rscript install/install_packages.R --all
args <- commandArgs(trailingOnly = TRUE); want <- function(f) f %in% args || "--all" %in% args
local({ r <- getOption("repos"); if (is.null(r) || is.na(r["CRAN"]) || identical(unname(r["CRAN"]), "@CRAN@")) {
  r["CRAN"] <- "https://cloud.r-project.org"; options(repos = r) } })
need <- function(p) { m <- p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]
  if (length(m)) { message("installing: ", paste(m, collapse = ", ")); install.packages(m) } }
need(c("heatwaveR", "dplyr", "tidyr", "purrr", "tibble", "rlang", "MASS", "mgcv", "dlnm",
       "spdep", "ranger", "ggplot2", "testthat"))
if (want("--spatial")) need(c("sdmTMB"))
for (p in c("heatwaveR","dlnm","mgcv","MASS","spdep","ranger","sdmTMB","ggplot2"))
  message(sprintf("  %-10s %s", p, if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "NOT INSTALLED"))
