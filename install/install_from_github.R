#!/usr/bin/env Rscript
# install_from_github.R -- fallback when CRAN is blocked but GitHub is reachable.
# VERIFIED on Ubuntu 24.04, R 4.3.3 (this is the recipe that produced the
# deposited results). The `cran` GitHub organisation mirrors every CRAN
# package including release tags, so specific versions can be fetched.
#
# Two version constraints that this recipe already handles:
#  * dlnm master requires R >= 4.4; on R 4.3.x install release tag 2.4.7
#    (a genuine CRAN release requiring only R >= 3.2). Never edit DESCRIPTION.
#  * sdmTMB main and v1.2.0 require RTMB, which in turn needs a newer TMB
#    than Ubuntu ships. sdmTMB tag v1.1.0 (DESCRIPTION 1.0.0.9020) has no
#    RTMB dependency and is what the deposited results used.
tmp <- file.path(tempdir(), "rsrc"); dir.create(tmp, showWarnings = FALSE)
gh <- function(pkg, repo = paste0("cran/", pkg), ref = "master", subdir = NULL) {
  if (requireNamespace(pkg, quietly = TRUE)) { message("present: ", pkg); return(invisible()) }
  d <- file.path(tmp, pkg); unlink(d, recursive = TRUE)
  st <- system2("git", c("clone", "-q", "--depth", "1", "--branch", ref, paste0("https://github.com/", repo), d))
  if (st != 0) stop("clone failed: ", repo, "@", ref)
  install.packages(if (is.null(subdir)) d else file.path(d, subdir), repos = NULL, type = "source")
}
gh("fasttime", "s-u/fasttime");           gh("heatwaveR", "robwschlegel/heatwaveR")
gh("tsModel");                             gh("dlnm", ref = if (getRversion() >= "4.4.0") "master" else "2.4.7")
gh("assertthat"); gh("fishMod"); gh("reformulas"); gh("splancs"); gh("fmesher")
# sdmTMB compile needs ~4-6 GB RAM. If OOM-killed the failure is SILENT (log
# truncates mid-warning). Lower optimisation reduces peak memory:
#   mkdir -p ~/.R && printf 'CXXFLAGS=-O1 -g0\n' > ~/.R/Makevars   (remove afterwards)
gh("sdmTMB", "pbs-assess/sdmTMB", ref = "v1.1.0")
for (p in c("heatwaveR","dlnm","mgcv","MASS","spdep","ranger","sdmTMB"))
  message(sprintf("  %-10s %s", p, if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "NOT INSTALLED"))
