# =============================================================================
# 02_build_initiator_cohort.R — Documented antibiotic administrations and unified cohort of
#   initiators (PTZ / cefepime / meropenem) with covariates and outcomes.
#   -> output/extraction/derived_data/antipseudo_v31.rds (+ flow, administrations)
# Usage: Rscript scripts/02_build_initiator_cohort.R
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(dplyr)})
PKG_AVAILABLE <- function(p) isTRUE(requireNamespace(p, quietly = TRUE))
for (f in c("R/01_paths_and_io.R","R/02_codebooks.R","R/04_extract_labs.R","R/05_extract_icu_covariates.R",
            "R/20_bootstrap.R","R/31_baseline_creat.R","R/32_covariates.R","R/40_outcomes.R",
            "R/53_extract.R","R/55_outcomes.R","R/57_antipseudo_cohort.R")) source(f)
options(v3.tables_dir = "output/extraction/tables")
cfg <- load_config("config/config.yml"); tbl_paths <- detect_mimic_tables(cfg$paths$mimic_root)
creat <- readRDS("output/extraction/derived_data/creatinine_all.rds")
meds <- readRDS("output/extraction/derived_data/medications_all_v3.rds")
ab <- extract_administered_antibiotics_v31(tbl_paths)
print(ab[, .(n = .N, iv = round(mean(iv), 3), route_na = round(mean(is.na(route)), 3)), by = abx_class][order(-n)])
c0 <- build_antipseudo_initiators_v31(ab, tbl_paths, creat)
print(attr(c0, "flow"))
ce <- readRDS("output/extraction/derived_data/chartevents_long_v31.rds"); lab <- readRDS("output/extraction/derived_data/labs_long_v31.rds")
win <- readRDS("output/extraction/derived_data/extraction_windows_v31.rds")
cw <- merge(c0[, .(hadm_id, T0)], win, by = "hadm_id", all.x = TRUE)
cat("T0 hors fenêtre d'extraction (constantes) :", sum(is.na(cw$t_hi) | cw$T0 > cw$t_hi - 0), "/", nrow(cw), "\n")
d <- add_covariates_outcomes_antipseudo_v31(c0, ab, meds, ce, lab, creat, tbl_paths)
saveRDS(d, "output/extraction/derived_data/antipseudo_v31.rds")
saveRDS(attr(c0, "flow"), "output/extraction/derived_data/antipseudo_flow_v31.rds")
print(d[, .N, by = treatment])
