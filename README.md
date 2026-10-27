# mhw-R-audit

Factorial Monte Carlo evaluation of modelling frameworks for marine heatwave (MHW) impact assessment, implemented entirely in R. Companion to:

> Afrifa-Yamoah, E. *Marine heatwave impact assessment on marine fauna and fisheries: a methodological synthesis and a simulation based evaluation of modelling frameworks.* Submitted to Reviews in Fish Biology and Fisheries.

Detection runs through `heatwaveR`; models through `MASS`, `mgcv`, `dlnm`, `sdmTMB` and `ranger`. Nothing is re-implemented and there are no cross language substitutions. Every number in the manuscript is a LaTeX macro generated from the deposited CSVs by `scripts/make_results_macros.R`.

## Quick start

```bash
export MHW_ROOT=$PWD
Rscript install/install_packages.R --all        # or install/install_from_github.R if CRAN is blocked
Rscript scripts/verify_detection.R               # REQUIRED: writes detection_reference.csv
Rscript scripts/diagnose_dgp.R                   # REQUIRED: confirms non-degenerate responses
Rscript -e 'library(testthat); test_dir("tests/testthat")'
bash hpc/run_local_parallel.sh 8                 # or: sbatch hpc/slurm_array.sh
Rscript scripts/aggregate.R
Rscript scripts/make_figures.R
Rscript scripts/make_results_macros.R
Rscript scripts/make_main_tables.R
Rscript scripts/make_supp_tables.R
cd manuscript && latexmk -pdf main.tex && latexmk -pdf supplement.tex
```

## Design

| Factor | Levels |
|---|---|
| Region (real OISST, heatwaveR example data) | WA, NW Atlantic, Mediterranean |
| True β₁ (per °C·day of lagged cumulative intensity) | −0.04, −0.08, −0.12 |
| True lag (years) | 0, 1, 2 |
| Spatial correlation range (km) | 0, 50 |
| Grid size (sites) | 9, 25 |
| Replicates | 10 |
| Models | naive NB GLM, lag-1 NB GLM, NB GAM, DLNM, sdmTMB, random forest |

108 scenarios; 120 unique exposure panels (exposure depends only on region × grid × range × replicate); 1,080 datasets.

DGP: `log μ = 2.5 + β₁ z + (−5e−4) z² + AR1(ρ = 0.4, σ = 0.3)`, `y ~ NegBin(μ, size = 4)`, with `z` the lagged annual cumulative intensity above threshold.

## Detection reference values (regenerate, never copy)

`scripts/verify_detection.R` writes `results/tables/detection_reference.csv` on your machine with your heatwaveR version, and the test suite checks against that file. With heatwaveR 1.0.0 and 0.5.5.9000 alike, the canonical algorithm gives, for the 2011 Western Australian event: start 2011-02-06, end 2011-04-06, 60 days, peak 6.51 °C, cumulative 192.71 °C·days. It is the largest event in that series by peak anomaly; a 1999 event is larger on the cumulative metric.

## Layout

```
R/00_config.R      design, seeds, paths           scripts/verify_detection.R
R/01_detection.R   heatwaveR wrappers             scripts/diagnose_dgp.R
R/02_exposure.R    metrics, lag panels            scripts/run_chunk.R
R/03_panel.R       OISST + Matern panel           scripts/aggregate.R
R/04_dgp.R         response simulator             scripts/make_figures.R
R/05_models.R      six fitters                    scripts/make_results_macros.R
R/06_runner.R      checkpointed runner            scripts/make_supp_tables.R
R/07_aggregate.R   summaries, run report          hpc/slurm_array.sh, pbs_array.sh, run_local_parallel.sh
tests/testthat/test-detection.R                   manuscript/main.tex, supplement.tex, refs.bib
```

## Running a slow model as a separate pass

Panels are deterministic in their seeds, so a model can be run later and merged:

```bash
MHW_PANELS_SUBDIR=panels_sdmtmb Rscript scripts/run_chunk.R --task 1 --ntasks 8 --models sdmtmb
```

`scripts/aggregate.R` reads every `results/panels*` directory and de-duplicates on (panel, β₁, lag, model). This is how the deposited sdmTMB results were produced.

## Before quoting a number

Read the run report printed by `scripts/aggregate.R` (also Table S9 in the supplement). A model that failed on part of the design still produces a plausible looking mean. Check `pct_estimate`, `pct_converge`, and `n_ci` in `R_scenario_summary.csv`.

## Licence

MIT. The SST series are heatwaveR's example data (Schlegel & Smit 2018), traceable to NOAA OISST v2.1 (Banzon et al. 2016).
