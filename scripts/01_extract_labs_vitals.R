# =============================================================================
# 01_extract_labs_vitals.R — Per-admission extraction windows, then long-format labevents and
#   ICU chartevents before T0. -> output/extraction/derived_data/
# Usage: Rscript scripts/01_extract_labs_vitals.R
# =============================================================================
suppressPackageStartupMessages({library(data.table)})
PKG_AVAILABLE <- function(p) isTRUE(requireNamespace(p, quietly = TRUE))
source("R/01_paths_and_io.R"); source("R/53_extract.R")
cfg <- load_config("config/config.yml"); tbl_paths <- detect_mimic_tables(cfg$paths$mimic_root)
meds <- readRDS("output/extraction/derived_data/medications_all_v3.rds")
win <- build_extraction_windows_v31(meds)
saveRDS(win, "output/extraction/derived_data/extraction_windows_v31.rds")
lab <- extract_labs_long_v31(tbl_paths, win)
ce  <- extract_chartevents_long_v31(tbl_paths, win)
cat("DONE\n")
