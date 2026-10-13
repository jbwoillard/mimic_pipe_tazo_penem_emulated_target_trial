# =============================================================================
# 00_extract_base.R — Base MIMIC-IV extractions: all valid creatinine values and all
#   medications (eMAR > pharmacy > prescriptions).
#   -> output/extraction/derived_data/{creatinine_all, medications_all_v3}.rds
# Usage: Rscript scripts/00_extract_base.R
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(dplyr) })
PKG_AVAILABLE <- function(p) isTRUE(requireNamespace(p, quietly = TRUE))
for (f in c("R/01_paths_and_io.R", "R/02_codebooks.R", "R/03_extract_medications.R",
            "R/03b_extract_medications_extended.R", "R/04_extract_labs.R")) source(f)
cfg <- load_config("config/config.yml"); tbl_paths <- detect_mimic_tables(cfg$paths$mimic_root)
dir.create("output/extraction/derived_data", recursive = TRUE, showWarnings = FALSE)
creat <- extract_creatinine(tbl_paths, cfg)          # écrit output/extraction/derived_data/creatinine_all.rds
meds  <- extract_medications_v3(tbl_paths, cfg)      # écrit output/extraction/derived_data/medications_all_v3.rds
cat("DONE : créatinine", nrow(creat), "lignes ; médicaments", nrow(meds), "lignes\n")
