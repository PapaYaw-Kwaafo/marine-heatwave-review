#!/usr/bin/env Rscript
# =============================================================================
# prepare_literature_figures.R -- stage the author's literature synthesis
# figure panels for the manuscript.
#
# Input : a folder of the author's exported panels (PNG and PDF), as supplied
#         in Lit_rev_figs.zip. Usage:
#           Rscript scripts/prepare_literature_figures.R <input_folder>
# Output: manuscript/figures/lit_<panel>.png, each trimmed of white margins and
#         downsampled to at most 2400 px on the long side (about 8 in at 300 dpi),
#         so that the compiled PDFs stay a sensible size. Panels are combined
#         into multi panel figures in LaTeX (minipages), not here, so that panel
#         letters live in the document rather than inside the image.
#
# The 12800 x 7200 px map exceeds the default ImageMagick disk cache policy on
# some systems; if image_read() fails with "cache resources exhausted", point
# MAGICK_CONFIGURE_PATH at a folder holding a policy.xml with larger resource
# limits before running this script.
# =============================================================================
suppressPackageStartupMessages(library(magick))

args <- commandArgs(trailingOnly = TRUE)
src  <- if (length(args)) args[1] else stop("give the folder holding the exported panels")
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
dst  <- file.path(root, "manuscript", "figures")

# Map from the author's file names to the panel names used in the .tex files.
panels <- c(
  "PRISMA.png"                           = "lit_prisma",                   # Fig. S1
  "Fig3 Co-occurrence.png"               = "lit_cooccurrence_submetric",   # Fig. 1a
  "Fig3 Family Co-occurrence.png"        = "lit_cooccurrence_family",      # Fig. 1b
  "word cloud.pdf"                       = "lit_wordcloud",                # Fig. 1c
  "Fig3 p_submetric_facets.png"          = "lit_submetric_counts",         # Fig. 1d
  "Fig4_ChoroplethEvents_v23.png"        = "lit_map",                      # Fig. 2
  "Fig5 Horizontal_Bars_clean.png"       = "lit_approach_by_year",         # Fig. 3a
  "Fig5 p_methods.png"                   = "lit_approach_by_family",       # Fig. 3b
  "_Fig4 Temporal Trends.png"            = "lit_publications_citations",   # Fig. S2a
  "Fig4 Response_Heatmap_wide.png"       = "lit_family_by_year",           # Fig. S2b
  "figure 4a_new.png"                    = "lit_network_institutions",     # Fig. S3a
  "Fig6 Chord_International_Collabs.png" = "lit_network_countries"         # Fig. S3b
)

read_panel <- function(path) {
  if (grepl("\\.pdf$", path, ignore.case = TRUE)) {
    # rasterise the PDF at 300 dpi through pdftoppm (poppler), then read the PNG
    tmp <- tempfile()
    system2("pdftoppm", c("-r", "300", "-png", "-singlefile", shQuote(path), shQuote(tmp)))
    image_read(paste0(tmp, ".png"))
  } else image_read(path)
}

for (f in names(panels)) {
  p <- file.path(src, f)
  if (!file.exists(p)) { message("missing: ", f); next }
  im <- read_panel(p)
  im <- image_background(im, "white")          # flatten any transparency
  im <- image_trim(im, fuzz = 2)               # drop white margins
  im <- image_border(im, "white", "24x24")     # restore a small uniform margin
  w  <- image_info(im)$width; h <- image_info(im)$height
  if (max(w, h) > 2400) im <- image_resize(im, if (w >= h) "2400x" else "x2400")
  out <- file.path(dst, paste0(panels[[f]], ".png"))
  image_write(im, out, format = "png", density = "300")
  cat(sprintf("%-40s -> %s (%d x %d)\n", f, basename(out), image_info(im)$width, image_info(im)$height))
}
