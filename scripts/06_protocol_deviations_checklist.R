# =============================================================================
# 06_protocol_deviations_checklist.R — Protocol deviations, TARGET checklist, HTE rows added to
#   results_master. No new effect estimate. -> output/analysis/
# Usage: Rscript scripts/06_protocol_deviations_checklist.R
# =============================================================================
options(warn = 1, scipen = 10)
suppressPackageStartupMessages({ library(data.table) })
OUT <- "output/analysis"; TAB <- function(f) file.path(OUT, "tables", f); FIG <- function(f) file.path(OUT, "figures", f)
MAN <- file.path(OUT, "manuscript")
for (d in c(MAN, file.path(MAN, c("tables", "figures", "etables", "efigures")))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
for (f in c("R/01_paths_and_io.R", "R/02_codebooks.R", "R/20_bootstrap.R", "R/30_results_master.R", "R/50_master_cohort.R", "R/51_analysis_cohort.R",
            "R/56_results_tables.R", "R/52_estimation.R", "R/71_registry.R", "R/77_tables_figs.R", "R/78_reporting.R")) source(f)

cat("=== P3.0 Registre des cohortes (reconstruit à l'identique) ===\n")
master <- setup_registry_v32(readRDS(file.path(OUT, "derived_data", "master_v32.rds")))$master
setkey(master, hadm_id)
reg_p2 <- fread(TAB("n_by_analysis_v32.csv"))[!grepl("^hte_tertile_", cohort_id)]   # tertiles ré-enregistrés en P3.1
get_analysis_cohort(list(id = "early_eligible", parent = "eligible_2arm", filter = EARLY_FILTER_V32,
                         label = reg_p2[cohort_id == "early_eligible", label]), master)
# acorn_ambbase dépend d'une issue recalculée en phase 2 : ligne reprise telle quelle du registre phase 2
.COHORT_REGISTRY_V31$rows[["acorn_ambbase"]] <- reg_p2[cohort_id == "acorn_ambbase"]
stopifnot(all.equal(cohort_registry_table()[match(reg_p2$cohort_id, cohort_id), n], reg_p2$n))

cat("=== P3.1 HTE v3.2 : lignes par tertile dans results_master ===\n")
rmt <- readRDS(file.path(OUT, "results", "results_master.rds"))
rmt <- rmt[!grepl("^hte_tertile_", analysis_id)]
pp <- readRDS(file.path(OUT, "hte", "derived", "hte_per_patient_v32.rds"))$pp
stopifnot(setequal(pp$hadm_id, .COHORT_REGISTRY_V31$ids[["broad"]]))
mini <- setkey(merge(master[J(pp$hadm_id), .(hadm_id, treatment)], pp[, .(hadm_id, tertile = as.character(tertile))], by = "hadm_id"), hadm_id)
for (tz in c("lower", "middle", "upper"))
  get_analysis_cohort(list(id = paste0("hte_tertile_", tz), parent = "broad", filter = bquote(tertile == .(tz)),
                           label = paste0("Broad, ", tz, " third of out-of-fold predicted effect (causal forests, v3.2 predictors without delay)")), mini)
register_strata_family("hte_tertiles_v32", "broad", paste0("hte_tertile_", c("lower", "middle", "upper")))
tt <- fread(file.path(OUT, "hte", "hte_tertile_effects_v32.csv"))
for (i in seq_len(nrow(tt))) {
  r <- tt[i]; id <- paste0("hte_tertile_", r$tertile)
  est <- list(n = r$n, n_ptz = r$n_ptz, n_carb = r$n_carb, n_events = r$events_ptz + r$events_carb,
              crude_risk_ptz = r$risk_obs_ptz, crude_risk_carb = r$risk_obs_carb, weighted_risk_ptz = r$dr_risk_ptz, weighted_risk_carb = r$dr_risk_carb,
              rd = r$rd, rd_ci_low = r$rd_lo, rd_ci_high = r$rd_hi, rr = r$rr, rr_ci_low = r$rr_lo, rr_ci_high = r$rr_hi,
              ess = NA_real_, ess_pct = NA_real_, max_smd = NA_real_, n_boot = NA_integer_, seed = NA_integer_, alpha = NA_real_, lambda = NA_real_,
              weight_scheme = "AIPW", ps_method = "grf", m = 1L)
  rmt <- add_row_v31(rmt, id, est, id, phase = "V32", role = "HTE", estimand = "Group ATE (doubly robust)", outcome = "outcome_primary",
                     description = sprintf("Broad (secondary): %s third of predicted effect, out-of-fold AIPW (causal forests, 13 predictors without delay).", r$tertile),
                     comparator = "meropenem")
  j <- nrow(rmt)
  set(rmt, j, "ci_method", "Out-of-fold AIPW scores (5-fold cross-fitting, nuisances fitted on other folds), influence-function SE, normal CI; imputed dataset 1")
  set(rmt, j, "population", "Broad (secondary)")
  set(rmt, j, "events_ptz", r$events_ptz); set(rmt, j, "events_comp", r$events_carb)
  pc <- precision_cols_v32(r$rd_lo, r$rd_hi); for (k in names(pc)) set(rmt, j, k, pc[[k]])
}

cat("=== P3.2 Sauvegardes results_master, registre, JSON ===\n")
fwrite(cohort_registry_table(), TAB("n_by_analysis_v32.csv"))
saveRDS(rmt, file.path(OUT, "results", "results_master.rds")); fwrite(rmt, file.path(OUT, "results", "results_master.csv"))
fwrite(rmt[, .(analysis_id, population, cohort_id, estimand, outcome, ps_method, m_imputations, n_boot, ci_method)], TAB("variance_methods_v32.csv"))
fwrite(rmt[, .(analysis_id, population, n, n_carb, rd, rd_ci_low, rd_ci_high, ci_half_width, mde80)][!is.na(rd)], TAB("precision_v32.csv"))
build_manuscript_numbers_json_v31(rmt, file.path(OUT, "results"))

cat("=== P3.3 C1 déviations, C2 TARGET, effet publié ACORN ===\n")
fwrite(protocol_deviations_v32(), TAB("protocol_deviations_v32.csv"))
fwrite(target_checklist_v32(), TAB("target_checklist_v32.csv"))
fwrite(ACORN_PUBLISHED_EFFECT_V32, TAB("acorn_published_effect_v32.csv"))
# clés uniques pour les nombres inline du rapport (ci_txt lit une seule colonne-clé)
cifk <- fread(TAB("cumulative_incidence_v32.csv"))[, key := paste0(fifelse(grepl("^Early", population), "early", "broad"), "_d", day)]
fwrite(cifk, TAB("cumulative_incidence_keyed_v32.csv"))
nck <- fread(TAB("negative_controls_by_step_v32.csv"))[, key := paste(population, step, control, sep = "_")]
fwrite(nck, TAB("negative_controls_by_step_keyed_v32.csv"))
# sessionInfo de l'environnement d'analyse (paquets du pipeline chargés)
invisible(lapply(c("glmnet", "mice", "cobalt", "nnet", "survival", "dagitty", "ggdag", "ggraph", "patchwork", "pROC", "EValue", "MASS", "grf",
                   "xgboost", "ggplot2", "ggtern", "PSweight", "flextable", "officer", "splines", "jsonlite"),
                 function(p) if (requireNamespace(p, quietly = TRUE)) loadNamespace(p)))
save_session_info_v31(file.path(OUT, "results"), SEED_V3)
file.rename(file.path(OUT, "results", c("session_info_v31.rds", "session_info_v31.txt")), file.path(OUT, "results", c("session_info_v32.rds", "session_info_v32.txt")))

cat("=== P3.4 C3 paquet manuscrit ===\n")
t1 <- fread(TAB("table1_early_v32.csv"))
tab1 <- t1[, .(Characteristic = characteristic, Statistic = type, `PTZ (unweighted)` = raw_ptz, `Meropenem (unweighted)` = raw_carb,
               `SMD unweighted` = sprintf("%.2f", smd_raw), `PTZ (overlap-weighted)` = w_ptz, `Meropenem (overlap-weighted)` = w_carb,
               `SMD weighted` = sprintf("%.2f", smd_weighted), `Missing, % PTZ / meropenem` = pct_missing_ptz_carb)]
tab2 <- manuscript_results_table_v32(rmt, names(TABLE2_IDS_V32), unname(TABLE2_IDS_V32))
tab3 <- rbind(data.table(Analysis = "ACORN trial (published): highest AKI stage or death by day 14, OR cefepime vs PTZ",
                         N = format(ACORN_PUBLISHED_EFFECT_V32$n, big.mark = ","),
                         `RR or OR (95% CI)` = with(ACORN_PUBLISHED_EFFECT_V32, sprintf("OR %.2f (%.2f to %.2f)", or_cefepime_vs_ptz, or_ci_low, or_ci_high))),
              manuscript_results_table_v32(rmt, names(TABLE3_IDS_V32), unname(TABLE3_IDS_V32)), fill = TRUE)
for (k in names(tab3)) set(tab3, which(is.na(tab3[[k]])), k, "")
cap <- c(
  t1 = sprintf("Table 1. Baseline characteristics of the early empirical population (primary), N = %s (PTZ %s, meropenem %s). Unweighted values are observed before imputation; overlap-weighted values and SMDs are averaged over %d imputations.",
               format(n_of("early"), big.mark = ","), format(cohort_registry_table()[cohort_id == "early", n_ptz], big.mark = ","),
               format(cohort_registry_table()[cohort_id == "early", n_comparator], big.mark = ","), rmt[analysis_id == "primary_early_ato", m_imputations]),
  t2 = "Table 2. AKI by day 7 and secondary outcomes, PTZ vs meropenem. First arm = PTZ; second arm = meropenem. Risks are overlap-weighted; RD = risk difference (percentage points); MDE80 = minimal detectable difference with 80% power (2.8 x SE).",
  t3 = "Table 3. Three-arm target trial (early population, generalized overlap weights) and ACORN benchmark (PTZ vs cefepime). First arm = first-named drug of the contrast. ACORN odds ratios are for cefepime vs PTZ (trial direction).")
TABS <- list(Table1_baseline_early = list(tab1, cap["t1"]), Table2_main_results = list(tab2, cap["t2"]), Table3_three_arm_and_benchmark = list(tab3, cap["t3"]))
for (nm in names(TABS)) { fwrite(TABS[[nm]][[1]], file.path(MAN, "tables", paste0(nm, ".csv")))
  save_docx_table_v32(TABS[[nm]][[1]], file.path(MAN, "tables", paste0(nm, ".docx")), TABS[[nm]][[2]]) }

HFIG <- function(f) file.path(OUT, "hte", f)
figs <- list(Figure1_flow = FIG("fig1_flow"), Figure2_dag = FIG("fig2_dag"), Figure3A_ps_distribution = FIG("fig3a_ps_distribution"),
             Figure3B_love_plot = FIG("fig3b_love_plot"), Figure4_cumulative_incidence = FIG("fig4_cumulative_incidence"),
             Figure5_forest_main = FIG("fig5_forest_main"), Figure6_trajectories = FIG("fig6_trajectories"))
efigs <- list(eFigure01_dag_detailed = FIG("efig_dag_detailed"), eFigure02_weights_distribution = FIG("efig_weights_distribution"),
              eFigure03_three_arm_gps = FIG("efig_positivity_three_arm_gps"), eFigure04_three_arm_ternary = FIG("efig_positivity_three_arm_ternary"),
              eFigure05_forest_secondary = FIG("efig_forest_secondary"), eFigure06_forest_sensitivity = FIG("efig_forest_sensitivity"),
              eFigure07_forest_negative_controls = FIG("efig_forest_negative_controls"), eFigure08_hte_tertiles = HFIG("fig_hte_tertiles_v32"),
              eFigure09_hte_period_validation = HFIG("efig_hte_tertiles_period_v32"), eFigure10_hte_calibration = HFIG("efig_hte_calibration_v32"),
              eFigure11_hte_importance = HFIG("efig_hte_importance_v32"), eFigure12_hte_dependence = HFIG("efig_hte_dependence_v32"))
cp <- function(lst, dir) for (nm in names(lst)) for (ext in c(".png", ".pdf")) {
  src <- paste0(lst[[nm]], ext); if (!file.exists(src)) stop("[P3] figure absente : ", src)
  file.copy(src, file.path(MAN, dir, paste0(nm, ext)), overwrite = TRUE) }
cp(figs, "figures"); cp(efigs, "efigures")

etabs <- c(eTable01_protocol_deviations = "protocol_deviations_v32.csv", eTable02_target_checklist = "target_checklist_v32.csv",
           eTable03_cohort_registry = "n_by_analysis_v32.csv", eTable04_dag_nodes = "dag_node_table_v32.csv",
           eTable05_dag_adjustment_sets = "dag_adjustment_sets_v32.csv", eTable06_positivity = "positivity_v32.csv",
           eTable07_positivity_three_arm = "positivity_three_arm_v32.csv", eTable08_table1_broad = "table1_broad_v32.csv",
           eTable09_table1_three_arm = "table1_three_arm_v32.csv", eTable10_table1_acorn = "table1_acorn_v32.csv",
           eTable11_acorn_trial_vs_emulation = "acorn_trial_vs_emulation_v32.csv", eTable12_acorn_ordinal_distribution = "acorn_ordinal_distribution_v32.csv",
           eTable13_overlap_population_early = "overlap_population_early_v32.csv", eTable14_overlap_population_broad = "overlap_population_broad_v32.csv",
           eTable15_cumulative_incidence = "cumulative_incidence_v32.csv", eTable16_reconciliation_2arm_3arm = "reconciliation_2arm_3arm_v32.csv",
           eTable17_waterfall = "waterfall_v32.csv", eTable18_negative_controls_by_step = "negative_controls_by_step_v32.csv",
           eTable19_per_protocol_censoring_early = "per_protocol_censoring_early.csv", eTable20_per_protocol_censoring_broad = "per_protocol_censoring_broad.csv",
           eTable21_evalues = "evalues_v32.csv", eTable22_precision = "precision_v32.csv", eTable23_variance_methods = "variance_methods_v32.csv",
           eTable24_observed_vs_imputed = "observed_vs_imputed_by_location_v32.csv", eTable25_vitals_missingness = "vitals_missingness_ed_v32.csv",
           eTable26_ps_method_decision = "ps_method_decision_v32.csv", eTable27_hte_tertiles = "hte_tertile_effects_v32.csv",
           eTable28_hte_metrics = "hte_metrics_v32.csv", eTable29_hte_period = "hte_period_tertile_effects_v32.csv",
           eTable30_hte_early_feasibility = "hte_early_feasibility_v32.csv", eTable31_trajectories = "trajectories_v32.csv")
for (nm in names(etabs)) file.copy(TAB(etabs[[nm]]), file.path(MAN, "etables", paste0(nm, ".csv")), overwrite = TRUE)
save_docx_table_v32(protocol_deviations_v32()[, .(ID = id, Version = version, Deviation = deviation, Reason = reason)],
                    file.path(MAN, "etables", "eTable01_protocol_deviations.docx"), "eTable 1. Deviations from the initial protocol.")
save_docx_table_v32(target_checklist_v32()[, .(Item = item, Section = section, Topic = topic, Location = location, Status = status)],
                    file.path(MAN, "etables", "eTable02_target_checklist.docx"), "eTable 2. TARGET reporting checklist (Cashin et al., JAMA 2025).")
fwrite(rbind(data.table(type = "table", file = names(TABS)), data.table(type = "figure", file = names(figs)),
             data.table(type = "etable", file = names(etabs), source = unname(etabs)),
             data.table(type = "efigure", file = names(efigs), source = basename(unlist(efigs))), fill = TRUE), file.path(MAN, "manifest.csv"))

# legends.md : chaque nombre lu dans results_master / le registre / les tables
R1 <- function(id) { r <- rmt[analysis_id == id]; stopifnot(nrow(r) == 1); r }
RDs <- function(id) with(R1(id), sprintf("%+.1f percentage points (95%% CI %.1f to %.1f)", 100 * rd, 100 * rd_ci_low, 100 * rd_ci_high))
NN <- function(id, f = "n") format(cohort_registry_table()[cohort_id == id][[f]], big.mark = ",")
pos <- fread(TAB("positivity_v32.csv")); cif <- fread(TAB("cumulative_incidence_v32.csv"))
leg <- c(
  "# Figure legends (generated from results_master; do not edit by hand)", "",
  sprintf("**Figure 1. Cohort flow.** Unified cohort of documented first administrations of PTZ, cefepime or meropenem (N = %s). Early population (primary): N = %s (PTZ %s, meropenem %s). Broad population (secondary): N = %s. Three-arm early trial: N = %s. ACORN emulation: N = %s.",
          NN("all"), NN("early"), NN("early", "n_ptz"), NN("early", "n_comparator"), NN("broad"), NN("three_arm_early"), NN("acorn")), "",
  "**Figure 2. Directed acyclic graph.** Measured nodes enter the propensity score; dashed nodes are unmeasured. The minimal adjustment sets and their coverage by the propensity score are listed in eTable 5.", "",
  sprintf("**Figure 3. Positivity and balance.** A, propensity-score distributions by arm (out-of-fold AUC %.2f in the early population, %.2f in the broad population). B, absolute standardized mean differences before and after overlap weighting (dashed line, 0.10).",
          pos[population == "Early (primary)", auc_out_of_fold], pos[population == "Broad (secondary)", auc_out_of_fold]), "",
  sprintf("**Figure 4. Cumulative incidence of AKI.** Overlap-weighted Aalen-Johansen estimator with death as a competing risk, follow-up censored at discharge. Difference at day 7 (PTZ minus meropenem): early %+.1f percentage points (95%% CI %.1f to %.1f); broad %+.1f (%.1f to %.1f). Bootstrap CIs.",
          100 * cif[population == "Early (primary)" & day == 7, rd], 100 * cif[population == "Early (primary)" & day == 7, rd_lo], 100 * cif[population == "Early (primary)" & day == 7, rd_hi],
          100 * cif[population == "Broad (secondary)" & day == 7, rd], 100 * cif[population == "Broad (secondary)" & day == 7, rd_lo], 100 * cif[population == "Broad (secondary)" & day == 7, rd_hi]), "",
  sprintf("**Figure 5. Main estimates.** Risk difference in AKI by day 7. Early (primary): %s. Broad (secondary): %s. Three-arm PTZ vs meropenem: %s. ACORN emulation, PTZ vs cefepime: %s. Multiple imputation (m = %d), bootstrap within imputations, Rubin's rules.",
          RDs("primary_early_ato"), RDs("secondary_broad_ato"), RDs("three_arm_ptz_vs_mero"), RDs("acorn_aki_d7"), R1("primary_early_ato")$m_imputations), "",
  "**Figure 6. Creatinine and urea trajectories.** Overlap-weighted daily means by arm (top) and difference with 95% CI (bottom), day 0 to 7; numbers of patients with a measurement shown per day.", "",
  sprintf("**eFigure 8. Heterogeneity (broad population).** Out-of-fold predicted effects (causal forests, 13 predictors, delay to T0 removed) and doubly robust effect by tertile: lower %s, middle %s, upper %s.",
          RDs("hte_tertile_lower"), RDs("hte_tertile_middle"), RDs("hte_tertile_upper")))
writeLines(leg, file.path(MAN, "legends.md"))
cat("Paquet manuscrit :", MAN, "\n"); print(tab2)
cat("PHASE 3 v3.2 (données) TERMINÉE\n")
