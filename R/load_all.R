# load_all.R -- source every module in dependency order. source("R/load_all.R")
.mhw_here <- function(...) file.path(Sys.getenv("MHW_ROOT", unset = "."), ...)
for (f in c("00_config.R", "01_detection.R", "02_exposure.R", "03_panel.R",
            "04_dgp.R", "05_models.R", "06_runner.R", "07_aggregate.R"))
  source(.mhw_here("R", f))
invisible(TRUE)
