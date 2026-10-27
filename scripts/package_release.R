#!/usr/bin/env Rscript
# package_release.R -- build the release tarball (R helper; no git identity written).
root <- Sys.getenv("MHW_ROOT", unset = "."); out <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out)) out <- file.path(dirname(normalizePath(root)), "mhw-R-audit-v2.0.0.tar.gz")
excl <- c("--exclude=logs/*.log", "--exclude=manuscript/*.aux", "--exclude=manuscript/*.bbl",
          "--exclude=manuscript/*.blg", "--exclude=manuscript/*.fls", "--exclude=manuscript/*.fdb_latexmk",
          "--exclude=manuscript/*.out", "--exclude=manuscript/*.toc", "--exclude=manuscript/build_*.log",
          "--exclude=results/panels*/*.tmp")
st <- system2("tar", c("-czf", shQuote(out), excl, "-C", shQuote(dirname(normalizePath(root))),
                       shQuote(basename(normalizePath(root)))))
cat(if (st == 0) sprintf("wrote %s (%.1f MB)\n", out, file.size(out) / 1e6) else "tar failed\n")
