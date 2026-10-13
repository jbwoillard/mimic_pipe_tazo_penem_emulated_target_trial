# =============================================================================
# 03_master_cohort_imputation_primary.R — Master cohort (single definition for all populations), creatinine
#   trajectory, ED vital signs, multiple imputation stratified by care location
#   (m = 20), propensity-score rule, primary and three-arm estimates,
#   two-arm / three-arm reconciliation. -> output/analysis/
# Usage: Rscript scripts/03_master_cohort_imputation_primary.R          (PIPELINE_FAST=1: few imputations/replicates, debugging)
# =============================================================================
options(warn = 1, scipen = 10, v31.cores = 8L)
suppressPackageStartupMessages({ library(data.table); library(glmnet); library(dplyr) })
PKG_AVAILABLE <- function(p) isTRUE(requireNamespace(p, quietly = TRUE))
invisible(lapply(c("mice", "cobalt", "grf", "jsonlite", "nnet", "survival"), loadNamespace))
FAST <- Sys.getenv("PIPELINE_FAST") == "1"
OUT <- if (FAST) "output/analysis_fast" else "output/analysis"
M_IMP <- if (FAST) 3L else 20L; B_MAIN <- if (FAST) 10L else 100L; B_SEC <- if (FAST) 10L else 50L
options(v3.tables_dir = file.path(OUT, "tables"))
for (d in c("results", "tables", "figures", "derived_data", "diagnostics")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
TAB <- function(f) file.path(OUT, "tables", f)
T_START <- Sys.time()

for (f in c("R/01_paths_and_io.R", "R/02_codebooks.R", "R/03_extract_medications.R", "R/04_extract_labs.R",
            "R/05_extract_icu_covariates.R", "R/10_estimators.R", "R/15_positivity_diagnostics.R",
            "R/20_bootstrap.R", "R/30_results_master.R", "R/31_baseline_creat.R", "R/32_covariates.R",
            "R/40_outcomes.R", "R/50_master_cohort.R", "R/51_analysis_cohort.R", "R/52_estimation.R",
            "R/53_extract.R", "R/54_mi.R", "R/55_outcomes.R", "R/56_results_tables.R",
            "R/57_antipseudo_cohort.R", "R/58_designs.R", "R/61_hte.R",
            "R/70_master.R", "R/71_registry.R", "R/72_ed_vitals.R", "R/73_mi_stratified.R")) source(f)
cfg <- load_config("config/config.yml"); tbl_paths <- detect_mimic_tables(cfg$paths$mimic_root)
creat_dt <- readRDS("output/extraction/derived_data/creatinine_all.rds")

cat("\n=== P1.1 Cohorte maîtresse v3.2 ===\n")
mpath <- "output/analysis/derived_data/master_v32.rds"
dir.create(dirname(mpath), recursive = TRUE, showWarnings = FALSE)
dir.create("output/analysis/tables", recursive = TRUE, showWarnings = FALSE)
if (file.exists(mpath)) master <- readRDS(mpath) else {
  ap <- readRDS("output/extraction/derived_data/antipseudo_v31.rds")
  ce <- readRDS("output/extraction/derived_data/chartevents_long_v31.rds")
  ed_long <- extract_ed_vitals_long_v32(unique(ap$hadm_id))
  master <- build_master_v32(ap, ce, ed_long, creat_dt, tbl_paths, out_tables = "output/analysis/tables")
  dir.create("output/analysis/tables", recursive = TRUE, showWarnings = FALSE)
  saveRDS(master, mpath); rm(ce); gc()
}
if (FAST) file.copy("output/analysis/tables/vitals_missingness_ed_v32.csv", TAB("vitals_missingness_ed_v32.csv"), overwrite = TRUE)

cat("\n=== P1.2 Registre ===\n")
master <- setup_registry_v32(master)$master
reg <- cohort_registry_table(); print(reg[, .(cohort_id, parent_id, n, n_ptz, n_comparator)])
QA_EARLY <- qa_early_ids_v32(master)
cat("QA early 2-arm IDs == three-arm IDs excluding cefepime :", QA_EARLY, "\n")

cat("\n=== P1.3 Imputations (A : constantes réa ; B : + ED ; C : + stratification par lieu = FINALE) ===\n")
EXTRA <- c(MI_EXTRA_PREDICTORS_V31, "acorn_ordinal_d14", "make14", "death_14d", "early_creat_rise_6h")
root <- master[J(.COHORT_REGISTRY_V31$ids[["all"]])]
mi_A <- impute_cohort_v31(root, icu_vitals_covs(PS_COVARIATES_V32), m = M_IMP, maxit = if (FAST) 5 else 20, extra_predictors = EXTRA,
                          cache = file.path(OUT, "derived_data", sprintf("mi_A_icu_vitals_m%d.rds", M_IMP)))
mi_B <- impute_cohort_v31(root, PS_COVARIATES_V32, m = M_IMP, maxit = if (FAST) 5 else 20, extra_predictors = EXTRA,
                          cache = file.path(OUT, "derived_data", sprintf("mi_B_ed_vitals_m%d.rds", M_IMP)))
mi_C <- impute_stratified_v32(root, PS_COVARIATES_V32, "care_location_t0", m = M_IMP, maxit = if (FAST) 5 else 20,
                              extra_predictors = EXTRA, cache = file.path(OUT, "derived_data", sprintf("mi_C_stratified_m%d.rds", M_IMP)))
IMP <- list(A = lapply(mi_A$imp, function(d) setkey(d, hadm_id)), B = lapply(mi_B$imp, function(d) setkey(d, hadm_id)), C = mi_C$imp)
D <- function(cohort_id, which = "C") get_analysis_data(cohort_id, IMP[[which]])
fwrite(rbind(observed_vs_imputed_v32(root, IMP$B, c("hr_mean", "map_mean", "lactate_max"), "care_location_t0")[, model := "non-stratified"],
             observed_vs_imputed_v32(root, IMP$C, c("hr_mean", "map_mean", "lactate_max"), "care_location_t0")[, model := "stratified by care location"]),
       TAB("observed_vs_imputed_by_location_v32.csv"))

cat("\n=== P1.4 Estimation : règle A2 par population, waterfall, contrôles négatifs ===\n")
rm_tbl <- init_results_master()
rm_tbl[, `:=`(cohort_id = character(), ps_method = character(), m_imputations = integer(), fmi_rd = numeric(), ci_method = character())]
ADD <- function(id, est, cohort_id, role, outcome, description, ci_method, estimand = "ATO", comparator = "meropenem") {
  rm_tbl <<- add_row_v31(rm_tbl, id, est, cohort_id, phase = "V32P1", role = role, estimand = estimand, outcome = outcome,
                         description = description, comparator = comparator)
  set(rm_tbl, nrow(rm_tbl), "ci_method", ci_method)
}
CI_ATO <- function(B) sprintf("Bootstrap SE within each of m=%d imputations (B=%d replicates, PS re-fitted in every replicate; ridge lambda fixed at the CV value), combined by Rubin's rules with Barnard-Rubin df", M_IMP, B)
# règle préspécifiée A2 appliquée dans chaque population et à chaque étape
est_rule <- function(cohort_id, which, covs, outcome, B, label) {
  dd <- D(cohort_id, which)
  er <- ato_mi_v31(dd, covs, outcome, "ridge", B = B, label = paste(label, "ridge"))
  br <- compute_balance_v31(dd, covs, er$weights_by_imputation)
  if (max(abs(br$Diff.Adj), na.rm = TRUE) > 0.10) {
    eg <- ato_mi_v31(dd, covs, outcome, "glm", B = B, label = paste(label, "glm"))
    list(sel = eg, alt = er, method = "glm", max_smd_ridge = max(abs(br$Diff.Adj), na.rm = TRUE))
  } else list(sel = er, alt = NULL, method = "ridge", max_smd_ridge = max(abs(br$Diff.Adj), na.rm = TRUE))
}
STEPS <- list(
  w1_unified_cohort = list(which = "A", covs = icu_vitals_covs(PS_COVARIATES_V31), desc = "Unified cohort of documented first administrations (B1); v3.1 PS covariates; ICU-only vital signs; non-stratified MI"),
  w2_creat_trajectory = list(which = "A", covs = icu_vitals_covs(PS_COVARIATES_V32), desc = "+ pre-T0 creatinine trajectory in the PS (B2)"),
  w3_ed_vitals = list(which = "B", covs = PS_COVARIATES_V32, desc = "+ emergency-department vital signs (MIMIC-IV-ED) (B3)"),
  w4_stratified_mi = list(which = "C", covs = PS_COVARIATES_V32, desc = "+ multiple imputation stratified by care location (B3) = FINAL v3.2"))
v31 <- readRDS("output/preliminary/results_master.rds")
wf <- list(); nc <- list(); FINAL <- list()
for (pop in c("broad", "early")) {
  v31_id <- if (pop == "broad") "primary_ato_v31" else "coprimary_early_empirical"
  r <- v31[analysis_id == v31_id]
  wf[[paste(pop, "w0")]] <- data.table(population = pop, step = "w0_v31_published", description = "v3.1 published estimate (cohort inherited from v2)",
                                       n = r$n, n_comp = r$n_carb, rd = r$rd, rd_lo = r$rd_ci_low, rd_hi = r$rd_ci_high, ps_method = r$ps_method)
  for (s in names(STEPS)) {
    st <- STEPS[[s]]
    e <- est_rule(pop, st$which, st$covs, "outcome_primary", if (s == "w4_stratified_mi") B_MAIN else B_SEC, paste(pop, s))
    wf[[paste(pop, s)]] <- data.table(population = pop, step = s, description = st$desc, n = e$sel$n, n_comp = e$sel$n_carb,
                                      rd = e$sel$rd, rd_lo = e$sel$rd_ci_low, rd_hi = e$sel$rd_ci_high, ps_method = e$method,
                                      max_smd_ridge = e$max_smd_ridge)
    if (s %in% c("w1_unified_cohort", "w2_creat_trajectory", "w4_stratified_mi")) for (ncv in list(c("nc_a12", "early_creat_rise_12h"), c("nc_a6", "early_creat_rise_6h"), c("nc_b", "nc_hypoglycemia_d7"))) {
      en <- ato_mi_v31(D(paste0(pop, "_", ncv[1]), st$which), st$covs, ncv[2], e$method, B = B_SEC, label = paste(pop, s, ncv[1]))
      nc[[paste(pop, s, ncv[1])]] <- data.table(population = pop, step = s, control = ncv[1], n = en$n, rd = en$rd, rd_lo = en$rd_ci_low, rd_hi = en$rd_ci_high)
      if (s == "w4_stratified_mi") ADD(paste0(pop, "_negative_control_", ncv[1]), en, paste0(pop, "_", ncv[1]), "NEGATIVE_CONTROL", ncv[2],
                                       paste0(pop, " population: negative control ", ncv[1], " (final v3.2 PS and imputations)."), CI_ATO(B_SEC))
    }
    if (s == "w4_stratified_mi") FINAL[[pop]] <- e
  }
}
fwrite(rbindlist(wf, fill = TRUE), TAB("waterfall_v32.csv")); fwrite(rbindlist(nc), TAB("negative_controls_by_step_v32.csv"))
print(rbindlist(wf, fill = TRUE)[, .(population, step, n, rd = round(100 * rd, 2), lo = round(100 * rd_lo, 2), hi = round(100 * rd_hi, 2), ps_method)])
print(rbindlist(nc)[, .(population, step, control, rd = round(100 * rd, 2), lo = round(100 * rd_lo, 2), hi = round(100 * rd_hi, 2))])

ADD("primary_early_v32", FINAL$early$sel, "early", "PRIMARY", "outcome_primary",
    sprintf("PRIMARY (v3.2): early empirical initiation, PTZ vs meropenem, documented first administrations, ATO, %s PS with natural splines including pre-T0 creatinine trajectory, MI (m=%d) stratified by care location.", FINAL$early$method, M_IMP),
    CI_ATO(B_MAIN))
ADD("secondary_broad_v32", FINAL$broad$sel, "broad", "SECONDARY", "outcome_primary",
    sprintf("Broad population (all PTZ or meropenem new users, any timing), same specification, %s PS.", FINAL$broad$method), CI_ATO(B_MAIN))
for (pop in c("early", "broad")) if (!is.null(FINAL[[pop]]$alt))
  ADD(paste0(pop, "_alt_ps_ridge_v32"), FINAL[[pop]]$alt, pop, "SENSITIVITY", "outcome_primary", paste0(pop, ": ridge PS (not selected by the A2 rule)."), CI_ATO(B_MAIN))
fwrite(data.table(population = names(FINAL), ps_method = sapply(FINAL, `[[`, "method"), max_abs_smd_ridge = sapply(FINAL, `[[`, "max_smd_ridge")),
       TAB("ps_method_decision_v32.csv"))

cat("\n=== P1.5 Essai à 3 bras (cohorte précoce unique) ===\n")
g3 <- gow_mi_v31(D("three_arm_early"), PS_COVARIATES_V32, "outcome_primary", B = B_SEC)
for (cn in names(g3$contrasts))
  ADD(paste0("three_arm_", cn, "_v32"), g3$contrasts[[cn]], "three_arm_early", "THREE_ARM", "outcome_primary",
      paste0("Three-arm target trial, early population (same definition as the primary): ", cn, "; multinomial PS, generalized overlap weights."),
      sprintf("Bootstrap SE within each of m=%d imputations (B=%d, stratified by arm, multinomial PS re-fitted), Rubin's rules", M_IMP, B_SEC),
      estimand = "ATO_generalized", comparator = c(ptz_vs_mero = "meropenem", ptz_vs_cef = "cefepime", cef_vs_mero = "meropenem")[[cn]])

cat("\n=== P1.6 Réconciliation 2 bras / 3 bras (B1-iii) ===\n")
rec <- list()
add_rec <- function(step, desc, rd, lo, hi, n, n_mero) rec[[length(rec) + 1]] <<- data.table(step = step, description = desc, n = n, n_carbapenem = n_mero, rd = rd, rd_lo = lo, rd_hi = hi, ci_half_width = (hi - lo) / 2)
r <- v31[analysis_id == "coprimary_early_empirical"]
add_rec("r1", "v3.1 published early estimate (cohort inherited from v2: first eMAR event, any carbapenem, cefepime-exposed patients possible)", r$rd, r$rd_ci_low, r$rd_ci_high, r$n, r$n_carb)
# r1b : cohorte héritée restreinte aux séjours aussi présents dans la cohorte précoce unifiée (mêmes covariables/imputations v3.1)
mv31 <- setkey(readRDS("output/preliminary/master_v31.rds"), hadm_id)
imp31 <- lapply(readRDS("output/preliminary/mi_eligible_m20.rds")$imp[seq_len(M_IMP)], function(d) setkey(d, hadm_id))
leg_ids <- intersect(.COHORT_REGISTRY_V31$ids[["early"]], mv31[!is.na(outcome_primary_v31), hadm_id])
e1b <- ato_mi_v31(lapply(imp31, function(d) d[J(leg_ids), nomatch = NULL][!is.na(outcome_primary_v31)]), PS_COVARIATES_V31, "outcome_primary_v31",
                  "glm", B = B_SEC, label = "r1b legacy∩unified")
add_rec("r1b", "v3.1 cohort and imputations, restricted to stays also in the unified early cohort (isolates cohort composition)", e1b$rd, e1b$rd_ci_low, e1b$rd_ci_high, e1b$n, e1b$n_carb)
w1e <- rbindlist(wf, fill = TRUE)[population == "early" & step == "w1_unified_cohort"]
add_rec("r2", "Unified early cohort, PTZ vs meropenem (cefepime excluded), binary overlap weights, v3.1 PS covariates (B1-i/ii)", w1e$rd, w1e$rd_lo, w1e$rd_hi, w1e$n, w1e$n_comp)
g3A <- gow_mi_v31(D("three_arm_early", "A"), icu_vitals_covs(PS_COVARIATES_V31), "outcome_primary", B = B_SEC)
x <- g3A$contrasts$ptz_vs_mero
add_rec("r3", "Same unified early cohort, three-arm generalized overlap weights (target population includes cefepime initiators), v3.1 PS covariates", x$rd, x$rd_ci_low, x$rd_ci_high, x$n, x$n_carb)
r4 <- v31[analysis_id == "three_arm_ptz_vs_mero"]
add_rec("r4", "v3.1 published three-arm PTZ minus meropenem", r4$rd, r4$rd_ci_low, r4$rd_ci_high, r4$n, r4$n_carb)
fe <- FINAL$early$sel; f3 <- g3$contrasts$ptz_vs_mero
add_rec("r5", "FINAL v3.2 primary (binary overlap weights, creatinine trajectory, ED vitals, stratified MI)", fe$rd, fe$rd_ci_low, fe$rd_ci_high, fe$n, fe$n_carb)
add_rec("r6", "FINAL v3.2 three-arm PTZ minus meropenem (same patients plus cefepime, generalized overlap weights)", f3$rd, f3$rd_ci_low, f3$rd_ci_high, f3$n, f3$n_carb)
rec <- rbindlist(rec); fwrite(rec, TAB("reconciliation_2arm_3arm_v32.csv"))
print(rec[, .(step, n, n_carbapenem, rd = round(100 * rd, 2), lo = round(100 * rd_lo, 2), hi = round(100 * rd_hi, 2), hw = round(100 * ci_half_width, 2))])

cat("\n=== P1.7 Sauvegardes, table des méthodes de variance, QA ===\n")
fwrite(cohort_registry_table(), TAB("n_by_analysis_v32.csv"))
saveRDS(rm_tbl, file.path(OUT, "results", "results_master.rds")); fwrite(rm_tbl, file.path(OUT, "results", "results_master.csv"))
fwrite(rm_tbl[, .(analysis_id, cohort_id, estimand, ps_method, m_imputations, ci_method)], TAB("variance_methods_v32.csv"))
saveRDS(list(PS = PS_COVARIATES_V32, FINAL_METHODS = sapply(FINAL, `[[`, "method"), M_IMP = M_IMP), file.path(OUT, "derived_data", "analysis_spec_v32.rds"))
cc <- qa_cohort_consistency_v31()
qa <- rbind(
  data.table(id = "QA-P1-1", check = "Every cohort is a subset of its parent; strata partition their parent exactly", ok = all(cc$ok), details = sprintf("%d/%d checks passed", sum(cc$ok), nrow(cc))),
  data.table(id = "QA-P1-2", check = "Single early-cohort definition: early 2-arm IDs equal three-arm IDs excluding cefepime", ok = QA_EARLY, details = ""),
  data.table(id = "QA-P1-3", check = "Every results_master row has a CI method", ok = all(nzchar(rm_tbl$ci_method)), details = sprintf("%d rows", nrow(rm_tbl))),
  data.table(id = "QA-P1-4", check = "Every results_master row has the N of its registered cohort", ok = all(rm_tbl$n == cohort_registry_table()[match(rm_tbl$cohort_id, cohort_id), n]), details = ""))
qa[, status := fifelse(ok, "PASS", "FAIL")][, ok := NULL]
fwrite(qa, TAB("qa_phase1_v32.csv")); print(qa)
cat(sprintf("PHASE 1 v3.2 TERMINÉE (%.1f min)\n", as.numeric(difftime(Sys.time(), T_START, units = "mins"))))
