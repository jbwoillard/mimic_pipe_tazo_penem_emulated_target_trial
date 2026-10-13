# =============================================================================
# 08_exposure_flow_cumulative_incidence.R — Overlap population, exposure (doses) and microbiology, STROBE flow
#   by arm, KDIGO trigger criterion, cumulative incidence with CI bands and
#   numbers at risk. Descriptive only. -> output/final/
# Usage: Rscript scripts/08_exposure_flow_cumulative_incidence.R
# =============================================================================
options(warn = 1, scipen = 10, v31.cores = 8L)
suppressPackageStartupMessages({ library(data.table); library(glmnet); library(ggplot2) })
PKG_AVAILABLE <- function(p) isTRUE(requireNamespace(p, quietly = TRUE))
invisible(lapply(c("mice", "cobalt", "jsonlite", "survival", "patchwork", "pROC", "duckdb", "DBI"), loadNamespace))
IN <- "output/analysis"; OUT <- "output/final"
TAB <- function(f) file.path(OUT, "tables", f); FIG <- function(f) file.path(OUT, "figures", f)
FAST <- Sys.getenv("PIPELINE_FAST") == "1"; B_CIF <- if (FAST) 5L else 200L
T_START <- Sys.time()
for (f in c("R/01_paths_and_io.R", "R/02_codebooks.R", "R/03_extract_medications.R", "R/04_extract_labs.R", "R/05_extract_icu_covariates.R",
            "R/10_estimators.R", "R/20_bootstrap.R", "R/30_results_master.R", "R/31_baseline_creat.R",
            "R/50_master_cohort.R", "R/51_analysis_cohort.R", "R/52_estimation.R", "R/53_extract.R",
            "R/55_outcomes.R", "R/56_results_tables.R", "R/57_antipseudo_cohort.R", "R/58_designs.R",
            "R/71_registry.R", "R/75_positivity_flow.R", "R/76_cif.R", "R/77_tables_figs.R",
            "R/80_robustness.R", "R/81_exposure_flow_cif.R")) source(f)
cfg <- load_config("config/config.yml"); tbl_paths <- detect_mimic_tables(cfg$paths$mimic_root)
creat_dt <- readRDS("output/extraction/derived_data/creatinine_all.rds")
ab <- readRDS("output/extraction/derived_data/abx_administered_v31.rds")
spec <- readRDS(file.path(IN, "derived_data", "analysis_spec_v32.rds")); PS <- spec$PS; METHOD <- spec$FINAL_METHODS
master <- setup_registry_v32(readRDS(file.path(IN, "derived_data", "master_v32.rds")))$master; setkey(master, hadm_id)
get_analysis_cohort(list(id = "early_eligible", parent = "eligible_2arm", filter = EARLY_FILTER_V32), master)
IMP <- readRDS(file.path(IN, "derived_data", "mi_C_stratified_m20.rds"))$imp[seq_len(if (FAST) 2 else 20)]
D <- function(id) get_analysis_data(id, IMP)
RAW <- function(id) master[J(.COHORT_REGISTRY_V31$ids[[id]])]
W <- list()
for (p in c("early", "broad")) { e <- ato_mi_v31(D(p), PS, "outcome_primary", METHOD[[p]], B = 1, label = paste(p, "weights"))
  W[[p]] <- Reduce(`+`, e$weights_by_imputation) / length(IMP); W[[paste0(p, "_e")]] <- e }
nicelab <- function(x) { extra <- c(creat_last_pre_t0 = "Last pre-T0 creatinine, mg/dL", creat_delta_last_min7d = "Last minus 7-day minimum creatinine, mg/dL",
  creat_ratio_last_min7d = "Last / 7-day minimum creatinine", creat_slope_48h = "Creatinine slope over 48 h, mg/dL per day",
  creat_slope_48h_available = "Creatinine slope available", creat_single_pre_value = "Single pre-T0 creatinine value")
  ifelse(x %in% names(extra), extra[x], x) }

cat("=== B1 Population d'overlap (précoce et large) ===\n")
for (p in c("early", "broad")) {
  op <- build_overlap_population_table_v31(D(p), W[[paste0(p, "_e")]]$weights_by_imputation, PS)
  op[, characteristic := nicelab(characteristic)][, difference := overlap_population - full_cohort]
  fwrite(op, TAB(sprintf("overlap_population_%s_v33.csv", p)))
}

cat("=== B2 Exposition (doses emar_detail) et microbiologie, population précoce ===\n")
re <- RAW("early")
doses_raw <- extract_index_doses_v33(re$subject_id, cfg$paths$mimic_root, file.path(OUT, "derived_data", "doses_index_v33.rds"))
fwrite(doses_raw[, .N, by = .(drug = fifelse(grepl("piperacillin", tolower(medication)), "ptz", "meropenem"), event_txt, dose_given_unit)][order(-N)], TAB("dose_units_v33.csv"))
doses <- dose_grams_v33(doses_raw)
ex <- exposure_patient_v33(re, ab, doses)
ex <- ex[match(re$hadm_id, hadm_id)]; stopifnot(identical(ex$hadm_id, re$hadm_id))
saveRDS(ex, file.path(OUT, "derived_data", "exposure_patient_v33.rds"))
exs <- summ_exposure_v33(ex, W$early); fwrite(exs, TAB("exposure_description_v33.csv")); print(exs[, .(characteristic, arm, n_available, unweighted, overlap_weighted)])
mb <- microbiology_v33(re, tbl_paths, W$early)
fwrite(mb$summary, TAB("microbiology_v33.csv")); fwrite(mb$top, TAB("microbiology_top_pathogens_v33.csv")); print(mb$summary); print(mb$top)

cat("=== B3 Flux STROBE par bras (constructeur rejoué) ===\n")
c0 <- build_antipseudo_initiators_v31(ab, tbl_paths, creat_dt)
fla <- attr(c0, "flow_by_arm"); fl31 <- readRDS("output/extraction/derived_data/antipseudo_flow_v31.rds")
chk_tot <- merge(fla[, .(n_v33 = sum(n)), by = step], fl31, by = "step")
fwrite(fla, file.path(OUT, "derived_data", "flow_by_arm_raw_v33.csv"))
fm <- flow_main_steps_v33(fla, master, .COHORT_REGISTRY_V31$ids)
fwrite(fm, TAB("flowchart_detailed_v33.csv")); print(fm[, .(step, n, PTZ, meropenem, cefepime, carbapenem_other, excl_n)])
reg <- cohort_registry_table()
branches <- data.table(from = c("no_aki_t0", "ptz_or_mero", "evaluable"), label = c(
  sprintf("ACORN emulation (benchmark): PTZ vs cefepime,\nED or medical ICU, first dose <=12 h, culture drawn\nN = %s (PTZ %s / FEP %s)", format(reg[cohort_id == "acorn", n], big.mark = ","),
          format(reg[cohort_id == "acorn", n_ptz], big.mark = ","), format(reg[cohort_id == "acorn", n_comparator], big.mark = ",")),
  sprintf("BROAD population (secondary): all timings,\nAKI status by day 7 evaluable\nN = %s (PTZ %s / MEM %s)", format(reg[cohort_id == "broad", n], big.mark = ","),
          format(reg[cohort_id == "broad", n_ptz], big.mark = ","), format(reg[cohort_id == "broad", n_comparator], big.mark = ",")),
  sprintf("THREE-ARM early trial: same criteria,\ncefepime initiators added\nN = %s (PTZ %s / MEM + FEP %s)", format(reg[cohort_id == "three_arm_early", n], big.mark = ","),
          format(reg[cohort_id == "three_arm_early", n_ptz], big.mark = ","), format(reg[cohort_id == "three_arm_early", n_comparator], big.mark = ","))))
for (ext in c(".png", ".pdf")) file.copy(file.path(IN, "figures", paste0("fig1_flow", ext)), FIG(paste0("efig_flow_populations", ext)), overwrite = TRUE)   # flux par populations (v3.2)
save_fig_v32(fig_flow_detailed_v33(copy(fm), branches), FIG("fig1_flow"), 22, 30)
# exclus sans créatinine de baseline vs cohorte unifiée
xb <- attr(c0, "excluded_no_baseline")
desc_grp <- function(d, grp) data.table(group = grp, n = nrow(d), age_median = stats::median(d$age), female_pct = 100 * mean(toupper(d$gender) == "F"),
  ed_arrival_pct = 100 * mean(!is.na(d$edregtime)), emergency_admission_pct = 100 * mean(grepl("EMER|URGENT|EW EMER|DIRECT EMER", d$admission_type)),
  ptz_pct = 100 * mean(d$treatment == "PTZ"), meropenem_pct = 100 * mean(d$treatment == "meropenem"), cefepime_pct = 100 * mean(d$treatment == "cefepime"),
  los_days_median = stats::median(as.numeric(difftime(d$dischtime, d$admittime, units = "days")), na.rm = TRUE),
  t0_from_admission_h_median = stats::median(as.numeric(difftime(d$T0, d$admittime, units = "hours")), na.rm = TRUE),
  in_hospital_death_pct = 100 * mean(!is.na(d$deathtime)))
excl_cmp <- rbind(desc_grp(xb, "Excluded: no baseline creatinine"), desc_grp(c0, "Included: unified cohort"))
fwrite(excl_cmp, TAB("excluded_no_baseline_vs_included_v33.csv")); print(excl_cmp)

cat("=== B4 Critère KDIGO déclencheur ===\n")
trig <- rbind(aki_trigger_table_v33(RAW("early"), W$early)[, population := "Early (primary)"],
              aki_trigger_table_v33(RAW("broad"), W$broad)[, population := "Broad (secondary)"])
fwrite(trig, TAB("aki_trigger_v33.csv")); print(trig[population == "Early (primary)"])

cat("=== B5 Incidence cumulée : bandes d'IC (B =", B_CIF, ") et effectifs à risque ===\n")
cb_cache <- file.path(OUT, "derived_data", sprintf("cif_bands_B%d_v33.rds", B_CIF))   # le bootstrap n'est refait que si le cache manque
CB <- if (file.exists(cb_cache)) readRDS(cb_cache) else list(`Early (primary)` = cif_bands_v33(D("early")[[1]], PS, METHOD[["early"]], B = B_CIF),
                                                             `Broad (secondary)` = cif_bands_v33(D("broad")[[1]], PS, METHOD[["broad"]], B = B_CIF))
saveRDS(CB, cb_cache)
fwrite(rbindlist(lapply(names(CB), function(p) CB[[p]]$curves[, population := p])), TAB("cumulative_incidence_bands_v33.csv"))
fwrite(rbindlist(lapply(names(CB), function(p) CB[[p]]$at_risk[, population := p])), TAB("cumulative_incidence_at_risk_v33.csv"))
save_fig_v32(fig_cif_bands_v33(CB), FIG("fig4_cumulative_incidence"), 18, 11)

cat("=== QA phase B ===\n")
qa <- rbind(
  data.table(id = "QA-V33-5", check = "Overlap population described for the early (primary) population",
             ok = file.exists(TAB("overlap_population_early_v33.csv")) && nrow(fread(TAB("overlap_population_early_v33.csv"))) > 20, details = ""),
  data.table(id = "QA-V33-6", check = "Flow decomposed step by step with N by arm; totals identical to v3.1/v3.2; arms sum to totals",
             ok = all(chk_tot$n_v33 == chk_tot$n) && fm[step == "evaluable", n] == n_of("early") && fm[step == "no_esrd", n] == n_of("all") &&
               fm[step == "no_iv_48h", n] == n_of("early_eligible") && fm[step == "ptz_or_mero", n] == n_of("eligible_2arm"),
             details = sprintf("%d steps; unified %s; early %s", nrow(fm), format(fm[step == "no_esrd", n], big.mark = ","), format(fm[step == "evaluable", n], big.mark = ","))),
  data.table(id = "QA-V33-B2", check = "Index-drug doses available for most early patients (emar_detail)",
             ok = mean(!is.na(ex$dose_g_first24h)) > 0.8, details = sprintf("%.1f%% with a first-24 h dose", 100 * mean(!is.na(ex$dose_g_first24h)))),
  data.table(id = "QA-V33-B4", check = "Trigger categories sum to the AKI events of each arm (early)",
             ok = all(trig[population == "Early (primary)" & !grepl("^AKI by|^Confirmed", criterion), sum(n), by = treatment][order(treatment)]$V1 ==
               RAW("early")[outcome_primary == 1L, .N, by = treatment][order(treatment)]$N), details = ""),
  data.table(id = "QA-V33-B5", check = "Cumulative incidence with pointwise bootstrap bands and numbers at risk (both populations)",
             ok = all(vapply(CB, function(z) nrow(z$curves) > 0 && !anyNA(z$curves$lo) && nrow(z$at_risk) == 10, logical(1))),
             details = paste(sprintf("%s: B=%d, %d failed", names(CB), B_CIF, vapply(CB, `[[`, numeric(1), "n_failed")), collapse = "; ")))
qa[, status := fifelse(ok, "PASS", "FAIL")][, ok := NULL]
fwrite(qa, TAB("qa_phaseB_v33.csv")); print(qa[, .(id, status, details)])
cat(sprintf("PHASE B v3.3 TERMINÉE (%.1f min)\n", as.numeric(difftime(Sys.time(), T_START, units = "mins"))))
