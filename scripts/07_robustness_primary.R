# =============================================================================
# 07_robustness_primary.R — Robustness of the primary estimate: PS complexity, parsimonious PS,
#   augmented overlap weighting, cross-fitted AIPW, B = 200 bootstrap per
#   imputation, rule of 20 events for the RR. The master cohort and the
#   imputations of step 3 are read-only. -> output/final/
# Usage: Rscript scripts/07_robustness_primary.R
# =============================================================================
options(warn = 1, scipen = 10, v31.cores = 8L)
suppressPackageStartupMessages({ library(data.table); library(glmnet); library(ggplot2) })
PKG_AVAILABLE <- function(p) isTRUE(requireNamespace(p, quietly = TRUE))
invisible(lapply(c("mice", "cobalt", "jsonlite", "nnet", "survival", "dagitty", "ggdag", "patchwork", "pROC", "EValue", "MASS", "PSweight", "diptest"), loadNamespace))
IN <- "output/analysis"; OUT <- "output/final"
TAB <- function(f) file.path(OUT, "tables", f); FIG <- function(f) file.path(OUT, "figures", f)
FAST <- Sys.getenv("PIPELINE_FAST") == "1"
B_EARLY <- if (FAST) 5L else 200L; B_BROAD <- if (FAST) 5L else 100L; B_SINGLE <- if (FAST) 10L else 500L
T_START <- Sys.time()
for (f in c("R/01_paths_and_io.R", "R/02_codebooks.R", "R/03_extract_medications.R", "R/04_extract_labs.R", "R/05_extract_icu_covariates.R",
            "R/10_estimators.R", "R/15_positivity_diagnostics.R", "R/20_bootstrap.R", "R/30_results_master.R", "R/43_sensitivity.R",
            "R/50_master_cohort.R", "R/51_analysis_cohort.R", "R/52_estimation.R", "R/53_extract.R",
            "R/54_mi.R", "R/55_outcomes.R", "R/56_results_tables.R", "R/57_antipseudo_cohort.R", "R/58_designs.R",
            "R/71_registry.R", "R/75_positivity_flow.R", "R/77_tables_figs.R", "R/78_reporting.R",
            "R/80_robustness.R")) source(f)

cat("=== A.0 Sorties v3.3 initialisées depuis v3.2 (copie), registre, données gelées ===\n")
for (d in c("tables", "figures", "results", "diagnostics", "hte", "derived_data")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file.path(OUT, "results", "results_master_v32_frozen.rds"))) {
  for (d in c("tables", "figures", "results")) file.copy(list.files(file.path(IN, d), full.names = TRUE), file.path(OUT, d), overwrite = TRUE)
  file.copy(list.files(file.path(IN, "hte"), full.names = TRUE, pattern = "\\.(csv|png|pdf|txt)$"), file.path(OUT, "hte"), overwrite = TRUE)
  file.copy(file.path(IN, "results", "results_master.rds"), file.path(OUT, "results", "results_master_v32_frozen.rds"))
}
rm32 <- readRDS(file.path(OUT, "results", "results_master_v32_frozen.rds"))
cfg <- load_config("config/config.yml"); tbl_paths <- detect_mimic_tables(cfg$paths$mimic_root)
creat_dt <- readRDS("output/extraction/derived_data/creatinine_all.rds")
ab <- readRDS("output/extraction/derived_data/abx_administered_v31.rds")
spec <- readRDS(file.path(IN, "derived_data", "analysis_spec_v32.rds")); PS <- spec$PS; METHOD <- spec$FINAL_METHODS
master <- setup_registry_v32(readRDS(file.path(IN, "derived_data", "master_v32.rds")))$master
setkey(master, hadm_id)
reg32 <- fread(file.path(IN, "tables", "n_by_analysis_v32.csv"))
get_analysis_cohort(list(id = "early_eligible", parent = "eligible_2arm", filter = EARLY_FILTER_V32,
                         label = reg32[cohort_id == "early_eligible", label]), master)
for (cid in setdiff(reg32$cohort_id, cohort_registry_table()$cohort_id)) .COHORT_REGISTRY_V31$rows[[cid]] <- reg32[cohort_id == cid]
stopifnot(isTRUE(all.equal(cohort_registry_table()[match(reg32$cohort_id, cohort_id), n], reg32$n)))
IMP <- lapply(readRDS(file.path(IN, "derived_data", "mi_C_stratified_m20.rds"))$imp[seq_len(if (FAST) 2 else 20)], add_parsimonious_vars_v33)
D <- function(id) get_analysis_data(id, IMP)
RAW <- function(id) master[J(.COHORT_REGISTRY_V31$ids[[id]])]
M <- length(IMP)
CI_MI <- function(B, extra = "") sprintf("Bootstrap SE within each of m=%d imputations (B=%d, PS re-fitted per replicate%s), Rubin's rules (Barnard-Rubin df)", M, B, extra)

# results_master v3.3 = lignes v3.2 ; les lignes recalculées remplacent les anciennes à la même position
rmt <- copy(rm32); NEW <- list()
ADD <- function(id, est, cohort_id, role, outcome, description, ci_method, population, estimand = "ATO", comparator = "meropenem", arm_stats = NULL) {
  tmp <- add_row_v31(rmt[0], id, est, cohort_id, phase = "V33", role = role, estimand = estimand, outcome = outcome,
                     description = description, comparator = comparator)
  set(tmp, 1L, "ci_method", ci_method); set(tmp, 1L, "population", population)
  if (!is.null(arm_stats)) for (k in names(arm_stats)) set(tmp, 1L, k, arm_stats[[k]])
  if (!is.na(est$rd)) { pc <- precision_cols_v32(est$rd_ci_low, est$rd_ci_high); for (k in names(pc)) set(tmp, 1L, k, pc[[k]]) }
  if (id %in% rmt$analysis_id) { i <- which(rmt$analysis_id == id); rmt <<- rbind(rmt[seq_len(i - 1)], tmp, rmt[-seq_len(i)], fill = TRUE) }
  else rmt <<- rbind(rmt, tmp, fill = TRUE)
  NEW[[id]] <<- est
  invisible(est)
}
ESTS <- list()

cat("\n=== A.3a Population précoce : toutes les analyses bootstrap sous IM avec B =", B_EARLY, "===\n")
early_mi <- rm32[population == "Early (primary)" & grepl("^Bootstrap SE within each of m=20 imputations \\(B=[0-9]+, PS re-fitted per replicate\\)", ci_method)]
for (j in seq_len(nrow(early_mi))) {
  r <- early_mi[j]; covs <- if (grepl("_vanco_", r$analysis_id)) setdiff(PS, "vanco_6h") else PS
  dl <- D(r$cohort_id)
  e <- tryCatch(ato_mi_v31(dl, covs, r$outcome, r$ps_method, B = B_EARLY, label = r$analysis_id), error = function(err) { message("[", r$analysis_id, "] FAILED: ", conditionMessage(err)); NULL })
  if (is.null(e)) next
  ESTS[[r$analysis_id]] <- e
  ADD(r$analysis_id, e, r$cohort_id, r$role, r$outcome, r$description, CI_MI(B_EARLY), r$population, estimand = r$estimand, comparator = r$comparator,
      arm_stats = enrich_arm_stats_v32(e, dl, r$outcome))
}
P <- "primary_early_ato"; e_prim <- ESTS[[P]]

cat("\n=== A.3b Médiane + indicateurs (jeu unique) : IC normal fondé sur le SE bootstrap, B =", B_SINGLE, "===\n")
r <- rm32[analysis_id == "early_sens_median_indicators"]
raw <- RAW("early"); mi <- paste0(VARS_WITH_MISSING_V32, "_missing"); mi <- mi[vapply(mi, function(v) stats::var(raw[[v]]) > 0, logical(1))]
mi <- setdiff(mi, paste0(unname(ICU_VITALS_V32), "_missing"))
em <- ato_boot_v31(raw, c(PS, mi), "outcome_primary", r$ps_method, B = B_SINGLE)
pct <- c(em$rd_ci_low, em$rd_ci_high)
em$rd_ci_low <- em$rd - 1.96 * em$se_rd; em$rd_ci_high <- em$rd + 1.96 * em$se_rd
em$rr_ci_low <- exp(log(em$rr) - 1.96 * em$se_logrr); em$rr_ci_high <- exp(log(em$rr) + 1.96 * em$se_logrr)
ADD(r$analysis_id, em, "early", "SENSITIVITY", "outcome_primary", r$description,
    sprintf("Normal CI from the bootstrap SE (B=%d, PS re-fitted per replicate), single dataset; percentile CI %+.1f to %+.1f pp shown for reference", B_SINGLE, 100 * pct[1], 100 * pct[2]),
    r$population, arm_stats = list(events_ptz = sum(raw$outcome_primary[raw$treatment == "PTZ"]), events_comp = sum(raw$outcome_primary[raw$treatment != "PTZ"])))

cat("\n=== A.3c Per-protocol IPCW, B =", B_EARLY, "===\n")
ppan <- build_pp_panel_v31(raw, ab, creat_dt, tbl_paths)
ppf <- function(d, ps_covs, outcome_col, seed, ...) pp_ipcw_boot_v31(d, ppan, ps_covs, METHOD[["early"]], B = B_EARLY, seed = seed)
epp <- pool_estimator_mi_v31(D("early"), ppf, PS, "outcome_primary", "early per-protocol"); epp$ps_method <- METHOD[["early"]]; epp$m <- M
r <- rm32[analysis_id == "early_per_protocol_ipcw"]
ADD(r$analysis_id, epp, "early", "SENSITIVITY", "outcome_primary", r$description,
    sprintf("Bootstrap SE within each of m=%d imputations (B=%d; PS and censoring models re-fitted), Rubin's rules", M, B_EARLY), r$population, estimand = r$estimand,
    arm_stats = list(events_ptz = sum(raw$outcome_primary[raw$treatment == "PTZ"]), events_comp = sum(raw$outcome_primary[raw$treatment != "PTZ"])))

cat("\n=== A.1 Complexité du PS et trois spécifications alternatives (sensibilités) ===\n")
PARS <- unlist(PS_PARSIMONIOUS_V33, use.names = FALSE)
d1 <- D("early")[[1]]; b1 <- D("broad")[[1]]
cx <- rbind(ps_complexity_v33(d1, PS, "glm", 3, "Early (primary)", "Full PS, unpenalized (PRIMARY)"),
            ps_complexity_v33(d1, PS, "ridge", 3, "Early (primary)", "Full PS, ridge"),
            ps_complexity_v33(d1, PARS, "glm", PS_PARSIMONIOUS_DF_V33, "Early (primary)", "Parsimonious PS (DAG minimal set 1, df<=2), unpenalized"),
            ps_complexity_v33(b1, PS, "ridge", 3, "Broad (secondary)", "Full PS, ridge (selected)"),
            ps_complexity_v33(b1, PS, "glm", 3, "Broad (secondary)", "Full PS, unpenalized"))
fwrite(cx, TAB("ps_complexity_v33.csv")); print(cx[, .(population, specification, n_parameters, n_collinear, comparator_per_parameter, auc_in_sample, auc_out_of_fold)])
fwrite(data.table(block = rep(names(PS_PARSIMONIOUS_V33), lengths(PS_PARSIMONIOUS_V33)), variable = PARS), TAB("ps_parsimonious_variables_v33.csv"))

dl <- D("early"); arm <- enrich_arm_stats_v32(list(weights_by_imputation = NULL), dl, "outcome_primary")[c("events_ptz", "events_comp")]
# (a) PS parcimonieux
ea <- ato_mi_v31(dl, PARS, "outcome_primary", "glm", B = B_EARLY, label = "early parsimonious PS", ps_df = PS_PARSIMONIOUS_DF_V33)
ESTS$early_sens_ps_parsimonious <- ea
ADD("early_sens_ps_parsimonious", ea, "early", "SENSITIVITY", "outcome_primary",
    sprintf("Early (primary): parsimonious PS (DAG minimal adjustment set 1, %d pre-specified variables, natural splines df<=%d), unpenalized, overlap weights.", length(PARS), PS_PARSIMONIOUS_DF_V33),
    CI_MI(B_EARLY), "Early (primary)", arm_stats = enrich_arm_stats_v32(ea, dl, "outcome_primary"))
# (b) ATO augmenté — contrôle de l'implémentation contre PSweight sur l'imputation 1
chk <- tryCatch(check_augmented_vs_psweight_v33(dl[[1]], PS, "outcome_primary"), error = function(err) { message("PSweight check failed: ", conditionMessage(err)); NULL })
if (!is.null(chk)) { fwrite(chk, TAB("augmented_ato_psweight_check_v33.csv")); print(chk) }
pool_custom <- function(fun, label, ...) {
  res <- mc_retry_v31(M, function(i) fun(dl[[i]], PS, "outcome_primary", seed = SEED_V3 + i, ...))
  n <- res[[1]]$n
  prd <- rubin_pool_v31(sapply(res, `[[`, "rd"), sapply(res, `[[`, "se_rd")^2, n - length(PS))
  prr <- rubin_pool_v31(log(sapply(res, `[[`, "rr")), sapply(res, `[[`, "se_logrr")^2, n - length(PS))
  message(sprintf("[A1] %-28s RD=%+.2f%% [%.2f;%.2f]", label, 100 * prd$est, 100 * prd$lo, 100 * prd$hi))
  list(n = n, n_ptz = sum(dl[[1]]$treatment == "PTZ"), n_carb = sum(dl[[1]]$treatment != "PTZ"), n_events = sum(dl[[1]]$outcome_primary),
       weighted_risk_ptz = mean(sapply(res, `[[`, "mu1")), weighted_risk_carb = mean(sapply(res, `[[`, "mu0")),
       crude_risk_ptz = NA_real_, crude_risk_carb = NA_real_, rd = prd$est, rd_ci_low = prd$lo, rd_ci_high = prd$hi,
       rr = exp(prr$est), rr_ci_low = exp(prr$lo), rr_ci_high = exp(prr$hi), fmi_rd = prd$fmi, ess_pct = NA_real_, max_smd = NA_real_,
       n_boot = NA_integer_, seed = SEED_V3, alpha = NA_real_, lambda = NA_real_, weight_scheme = NA_character_, m = M, res = res)
}
eb <- pool_custom(ato_augmented_v33, "augmented ATO", B = B_EARLY); eb$ps_method <- "ridge"; eb$n_boot <- B_EARLY * M; eb$weight_scheme <- "overlap (augmented)"
ESTS$early_sens_ato_augmented <- eb
ADD("early_sens_ato_augmented", eb, "early", "SENSITIVITY", "outcome_primary",
    "Early (primary): augmented (doubly robust) overlap-weighting estimator (Mao, Li & Greene 2019), ridge PS and ridge outcome models by arm.",
    sprintf("Bootstrap SE within each of m=%d imputations (B=%d; PS and outcome models re-fitted, ridge lambdas fixed at their CV values), Rubin's rules (Barnard-Rubin df)", M, B_EARLY),
    "Early (primary)", arm_stats = c(arm, list(ess_ptz = mean(sapply(eb$res, `[[`, "ess_ptz")), ess_comp = mean(sapply(eb$res, `[[`, "ess_comp")))))
# (c) AIPW hors pli (ATE)
ec <- pool_custom(aipw_crossfit_v33, "cross-fitted AIPW"); ec$ps_method <- "ridge (cross-fitted)"
ESTS$early_sens_aipw_crossfit <- ec
ADD("early_sens_aipw_crossfit", ec, "early", "SENSITIVITY", "outcome_primary",
    "Early (primary): cross-fitted AIPW (5 folds, ridge PS and ridge outcome models by arm, PS truncated at 0.01/0.99); estimand ATE, not ATO.",
    sprintf("Influence-function SE within each of m=%d imputations (5-fold cross-fitting), Rubin's rules (Barnard-Rubin df)", M),
    "Early (primary)", estimand = "ATE", arm_stats = arm)

alt <- function(id, label, auc) { r <- rmt[analysis_id == id]
  data.table(analysis_id = id, specification = label, estimand = r$estimand, n_parameters = NA_integer_, rd = r$rd, rd_lo = r$rd_ci_low, rd_hi = r$rd_ci_high,
             ci_half_width = r$ci_half_width, ess_ptz = r$ess_ptz, ess_comp = r$ess_comp, auc_out_of_fold = auc) }
cxe <- cx[population == "Early (primary)"]
tab_alt <- rbind(alt(P, "Unpenalized logistic PS, full covariate set (PRIMARY)", cxe[ps_method == "glm" & grepl("^Full", specification), auc_out_of_fold])[, n_parameters := cxe[ps_method == "glm" & grepl("^Full", specification), n_parameters]],
                 alt("early_ato_ridge_ps", "Ridge PS, full covariate set", cxe[ps_method == "ridge", auc_out_of_fold])[, n_parameters := cxe[ps_method == "ridge", n_parameters]],
                 alt("early_sens_ps_parsimonious", "(a) Parsimonious unpenalized PS (DAG minimal set, df<=2)", cxe[grepl("^Parsimonious", specification), auc_out_of_fold])[, n_parameters := cxe[grepl("^Parsimonious", specification), n_parameters]],
                 alt("early_sens_ato_augmented", "(b) Augmented overlap weighting (doubly robust)", cxe[ps_method == "ridge", auc_out_of_fold])[, n_parameters := cxe[ps_method == "ridge", n_parameters]],
                 alt("early_sens_aipw_crossfit", "(c) Cross-fitted AIPW (ATE)", mean(sapply(ec$res, `[[`, "auc_oof")))[, n_parameters := cxe[ps_method == "ridge", n_parameters]])
fwrite(tab_alt, TAB("primary_alt_ps_v33.csv")); print(tab_alt[, .(specification, rd = round(100 * rd, 2), lo = round(100 * rd_lo, 2), hi = round(100 * rd_hi, 2), auc_out_of_fold)])

cat("\n=== A.2 Règle du PS ré-appliquée (valeurs non arrondies) + sensibilité PS non pénalisé, population large ===\n")
rule_row <- function(pop, er) { b <- compute_balance_v31(D(pop), PS, er$weights_by_imputation); i <- which.max(abs(b$Diff.Adj))
  data.table(population = pop, max_abs_smd_ridge = abs(b$Diff.Adj[i]), term_at_max = b$term[i], worst_imputation_max_abs_smd = max(b$worst_abs_diff_adj, na.rm = TRUE),
             threshold = 0.10, rule = "unpenalized PS if max |SMD| (averaged over imputations, all terms) > 0.10",
             decision_recomputed = if (abs(b$Diff.Adj[i]) > 0.10) "glm" else "ridge", decision_v32 = METHOD[[pop]]) }
er_broad <- ato_mi_v31(D("broad"), PS, "outcome_primary", "ridge", B = 1, label = "broad ridge (balance only)")
psr <- rbind(rule_row("early", ESTS$early_ato_ridge_ps), rule_row("broad", er_broad))
psr[, consistent := decision_recomputed == decision_v32]
fwrite(psr, TAB("ps_rule_v33.csv")); print(psr)
dlb <- D("broad"); eg <- ato_mi_v31(dlb, PS, "outcome_primary", "glm", B = B_BROAD, label = "broad glm PS")
ESTS$broad_ato_glm_ps <- eg
ADD("broad_ato_glm_ps", eg, "broad", "SENSITIVITY", "outcome_primary", "Broad (secondary): unpenalized PS (not selected by the A2 rule; decision at the threshold).",
    CI_MI(B_BROAD), "Broad (secondary)", arm_stats = enrich_arm_stats_v32(eg, dlb, "outcome_primary"))

cat("\n=== A.3d Distribution bootstrap du primaire ===\n")
bl <- boot_long_v33(e_prim)
saveRDS(bl, file.path(OUT, "diagnostics", "bootstrap_primary_v33.rds"))
dip <- rbind(data.table(statistic = "risk difference", n_replicates = nrow(bl), dip = diptest::dip.test(bl$rd)$statistic, p = diptest::dip.test(bl$rd)$p.value),
             data.table(statistic = "log risk ratio", n_replicates = nrow(bl), dip = diptest::dip.test(bl$log_rr)$statistic, p = diptest::dip.test(bl$log_rr)$p.value))
dip[, `:=`(skewness = c(mean(scale(bl$rd)^3), mean(scale(bl$log_rr)^3)), unimodal = p > 0.05)]
dip[, n_failed := sum(sapply(e_prim$boot_by_imputation, function(b) sum(is.na(b$boot_rd))))]
fwrite(dip, TAB("bootstrap_primary_unimodality_v33.csv")); print(dip)
save_fig_v32(fig_boot_dist_v33(bl), FIG("efig_bootstrap_primary"), 18, 14)

cat("\n=== A.3e Règle des 20 événements, précision, sauvegardes ===\n")
# per-protocol de la population large : événements de la cohorte ITT (borne supérieure), pour appliquer la règle
rb <- RAW("broad"); rmt[analysis_id == "broad_per_protocol_ipcw", `:=`(events_ptz = sum(rb$outcome_primary[rb$treatment == "PTZ"]), events_comp = sum(rb$outcome_primary[rb$treatment != "PTZ"]))]
rmt <- apply_rr_rule_v33(rmt)
for (i in which(!is.na(rmt$rd))) { pc <- precision_cols_v32(rmt$rd_ci_low[i], rmt$rd_ci_high[i]); for (k in names(pc)) set(rmt, i, k, pc[[k]]) }
fwrite(cohort_registry_table(), TAB("n_by_analysis_v32.csv"))
saveRDS(rmt, file.path(OUT, "results", "results_master.rds")); fwrite(rmt, file.path(OUT, "results", "results_master.csv"))
fwrite(rmt[, .(analysis_id, population, cohort_id, estimand, outcome, ps_method, m_imputations, n_boot, ci_method)], TAB("variance_methods_v32.csv"))
fwrite(rmt[, .(analysis_id, population, n, n_carb, rd, rd_ci_low, rd_ci_high, ci_half_width, mde80)][!is.na(rd)], TAB("precision_v32.csv"))
fwrite(rmt[, .(analysis_id, population, events_comp, rr_raw, rr_raw_ci_low, rr_raw_ci_high, rr_status)], TAB("rr_reporting_v33.csv"))
build_manuscript_numbers_json_v31(rmt, file.path(OUT, "results"))
cmp <- merge(rm32[population == "Early (primary)", .(analysis_id, rd32 = rd, lo32 = rd_ci_low, hi32 = rd_ci_high, B32 = n_boot)],
             rmt[, .(analysis_id, rd33 = rd, lo33 = rd_ci_low, hi33 = rd_ci_high, B33 = n_boot, rr_status)], by = "analysis_id")
fwrite(cmp, TAB("early_v32_vs_v33_v33.csv"))

dev <- rbind(protocol_deviations_v32()[, where_reported := where_reported],
  data.table(id = sprintf("D%02d", 11:15), version = "v3.3",
    deviation = c("Bootstrap replicates raised to 200 per imputation for every early-population analysis (CI of the primary estimate recomputed; point estimate unchanged)",
                  "Risk ratios not reported when the comparator arm has fewer than 20 events",
                  "Three sensitivity analyses of the primary propensity score: parsimonious PS, augmented overlap weighting, cross-fitted AIPW (ATE)",
                  "Unpenalized PS in the broad population added as a sensitivity analysis",
                  "Median-plus-indicators sensitivity: normal CI from the bootstrap SE instead of the percentile interval"),
    reason = c("Monte Carlo precision of the intervals of rare outcomes", "Log risk-ratio variance unstable with few events",
               "Overfitting of the unpenalized PS (low comparator-per-parameter ratio, gap between in-sample and out-of-fold discrimination)",
               "Ridge decision at the threshold of the pre-specified rule", "Consistency with the other analyses"),
    where_reported = c("Primary result; variance table", "All result tables", "Primary estimate under alternative PS specifications", "Broad population", "Sensitivity analyses")))
fwrite(dev, TAB("protocol_deviations_v33.csv"))

cat("\n=== A.4 QA phase A ===\n")
p32 <- rm32[analysis_id == P]; p33 <- rmt[analysis_id == P]
early_boot <- rmt[population == "Early (primary)" & grepl("^Bootstrap|bootstrap SE", ci_method)]
qa <- rbind(
  data.table(id = "QA-V33-0", check = "Primary analysis unchanged: same analysis_id, cohort, N, estimand, PS method and point estimate as v3.2 (1e-12)",
             ok = p33$cohort_id == p32$cohort_id && p33$n == p32$n && p33$estimand == p32$estimand && p33$ps_method == p32$ps_method &&
               abs(p33$rd - p32$rd) < 1e-12 && abs(p33$weighted_risk_ptz - p32$weighted_risk_ptz) < 1e-12 && abs(p33$weighted_risk_carb - p32$weighted_risk_carb) < 1e-12,
             details = sprintf("RD v3.2 %+.4f [%.4f; %.4f] (B=%d) -> v3.3 %+.4f [%.4f; %.4f] (B=%d)", 100 * p32$rd, 100 * p32$rd_ci_low, 100 * p32$rd_ci_high, p32$n_boot,
                               100 * p33$rd, 100 * p33$rd_ci_low, 100 * p33$rd_ci_high, p33$n_boot)),
  data.table(id = "QA-V33-1", check = "New v3.3 analyses have role SENSITIVITY and a protocol-deviation entry",
             ok = all(rmt[analysis_id %in% c("early_sens_ps_parsimonious", "early_sens_ato_augmented", "early_sens_aipw_crossfit", "broad_ato_glm_ps"), role] == "SENSITIVITY"),
             details = "D11-D15 in protocol_deviations_v33.csv"),
  data.table(id = "QA-V33-2", check = "No risk ratio reported with fewer than 20 events in the comparator arm",
             ok = rmt[!is.na(rr) & !grepl("^ordinal_OR", estimand), all(!is.na(events_comp) & events_comp >= 20)],
             details = sprintf("%d rows with RR not estimated", sum(grepl("^not estim", rmt$rr_status)))),
  data.table(id = "QA-V33-3", check = "Every bootstrap-based early analysis uses B >= 200 per imputation (single-dataset analyses: B >= 200)",
             ok = all(grepl("B=(2[0-9]{2}|[3-9][0-9]{2})", early_boot$ci_method)), details = sprintf("%d early bootstrap rows", nrow(early_boot))),
  data.table(id = "QA-V33-7", check = "PS rule re-applied (> 0.10 strict) with unrounded values; decisions identical to v3.2",
             ok = all(psr$consistent), details = paste(sprintf("%s: %.6f -> %s", psr$population, psr$max_abs_smd_ridge, psr$decision_recomputed), collapse = "; ")),
  data.table(id = "QA-V33-A3", check = "Bootstrap distribution of the primary RD and log RR unimodal (Hartigan dip test)",
             ok = all(dip$unimodal), details = paste(sprintf("%s p=%.2f", dip$statistic, dip$p), collapse = "; ")),
  data.table(id = "QA-V33-A1b", check = "Augmented ATO implementation equals PSweight (same PS and outcome predictions, |diff| < 1e-6)",
             ok = !is.null(chk) && abs(diff(chk$rd)) < 1e-6, details = if (is.null(chk)) "check not run" else sprintf("|diff| = %.2e", abs(diff(chk$rd)))))
qa[, status := fifelse(ok, "PASS", "FAIL")][, ok := NULL]
fwrite(qa, TAB("qa_phaseA_v33.csv")); print(qa[, .(id, status, details)])
print(cmp[, .(analysis_id, rd32 = round(100 * rd32, 2), lo32 = round(100 * lo32, 2), hi32 = round(100 * hi32, 2), rd33 = round(100 * rd33, 2), lo33 = round(100 * lo33, 2), hi33 = round(100 * hi33, 2), rr_status)], nrows = 100)
cat(sprintf("PHASE A v3.3 TERMINÉE (%.1f min)\n", as.numeric(difftime(Sys.time(), T_START, units = "mins"))))
