# =============================================================================
# 04_secondary_sensitivity_figures.R — DAG, positivity, flow chart, cumulative incidence, secondary
#   outcomes, sensitivity analyses, negative controls, ACORN benchmark,
#   Table 1, per-arm risks, forest plots. -> output/analysis/
# Usage: Rscript scripts/04_secondary_sensitivity_figures.R
# =============================================================================
options(warn = 1, scipen = 10, v31.cores = 8L)
suppressPackageStartupMessages({ library(data.table); library(glmnet); library(ggplot2) })
PKG_AVAILABLE <- function(p) isTRUE(requireNamespace(p, quietly = TRUE))
invisible(lapply(c("mice", "cobalt", "jsonlite", "nnet", "survival", "dagitty", "ggdag", "patchwork", "pROC", "EValue", "MASS"), loadNamespace))
OUT <- "output/analysis"; TAB <- function(f) file.path(OUT, "tables", f); FIG <- function(f) file.path(OUT, "figures", f)
options(v3.tables_dir = file.path(OUT, "tables"))
FAST <- Sys.getenv("PIPELINE_FAST") == "1"   # débogage : 2 imputations, B minimal (sorties écrasées par le calcul complet)
B_MAIN <- if (FAST) 5L else 100L; B_SEC <- if (FAST) 5L else 50L; B_SINGLE <- if (FAST) 5L else 300L; B_CIF <- if (FAST) 5L else 200L
T_START <- Sys.time()
for (f in c("R/01_paths_and_io.R", "R/02_codebooks.R", "R/03_extract_medications.R", "R/04_extract_labs.R", "R/05_extract_icu_covariates.R",
            "R/10_estimators.R", "R/15_positivity_diagnostics.R", "R/20_bootstrap.R", "R/30_results_master.R", "R/43_sensitivity.R",
            "R/50_master_cohort.R", "R/51_analysis_cohort.R", "R/52_estimation.R", "R/53_extract.R",
            "R/54_mi.R", "R/55_outcomes.R", "R/56_results_tables.R", "R/57_antipseudo_cohort.R", "R/58_designs.R",
            "R/61_hte.R", "R/71_registry.R", "R/74_dag.R", "R/75_positivity_flow.R", "R/76_cif.R",
            "R/77_tables_figs.R")) source(f)
cfg <- load_config("config/config.yml"); tbl_paths <- detect_mimic_tables(cfg$paths$mimic_root)
creat_dt <- readRDS("output/extraction/derived_data/creatinine_all.rds")
ab <- readRDS("output/extraction/derived_data/abx_administered_v31.rds")
spec <- readRDS(file.path(OUT, "derived_data", "analysis_spec_v32.rds")); PS <- spec$PS; METHOD <- spec$FINAL_METHODS   # early = glm, broad = ridge (règle A2, option 1 validée)

cat("=== P2.0 Données, registre, colonnes complémentaires ===\n")
master <- setup_registry_v32(readRDS(file.path(OUT, "derived_data", "master_v32.rds")))$master
# B6 : baseline « comme dans l'essai » (valeur ambulatoire <= 365 j prioritaire, sinon min 7 j) -> outcome ordinal ACORN
m_alt <- copy(master)[, baseline_creat_min7d := fcoalesce(baseline_creat_ambulatory_365d, baseline_creat_min7d)]
k14 <- compute_kdigo_v31(creat_dt, m_alt, 14, "rrt_14d", suffix = "_d14_ambbase")
master <- merge(master, k14[, .(subject_id, hadm_id, aki_stage_max_d14_ambbase)], by = c("subject_id", "hadm_id"), all.x = TRUE)
master[, acorn_ordinal_d14_ambbase := fifelse(death_14d == 1L, 4L, aki_stage_max_d14_ambbase)]
setkey(master, hadm_id)
C <- function(id, parent, filter = NULL, outcome = NULL, label = id)
  invisible(get_analysis_cohort(list(id = id, parent = parent, filter = filter, outcome = outcome, label = label), master))
C("early_eligible", "eligible_2arm", EARLY_FILTER_V32, label = "Early population before requiring an evaluable outcome (missing-outcome bounds)")
C("acorn_ambbase", "acorn", outcome = "acorn_ordinal_d14_ambbase", label = "ACORN, outcome with outpatient-priority baseline evaluable")
IMP <- lapply(readRDS(file.path(OUT, "derived_data", "mi_C_stratified_m20.rds"))$imp[seq_len(if (FAST) 2 else 20)], function(d) {
  d <- merge(d, master[, .(hadm_id, acorn_ordinal_d14_ambbase)], by = "hadm_id", all.x = TRUE); setkey(d, hadm_id) })
D <- function(id) get_analysis_data(id, IMP)
RAW <- function(id) master[J(.COHORT_REGISTRY_V31$ids[[id]])]

rm_tbl <- init_results_master()
rm_tbl[, `:=`(cohort_id = character(), ps_method = character(), m_imputations = integer(), fmi_rd = numeric(), ci_method = character(),
              events_ptz = integer(), events_comp = integer(), ess_ptz = numeric(), ess_comp = numeric(),
              ci_half_width = numeric(), mde80 = numeric(), population = character())]
CI_MI <- function(B, extra = "") sprintf("Bootstrap SE within each of m=%d imputations (B=%d, PS re-fitted per replicate%s), Rubin's rules (Barnard-Rubin df)", length(IMP), B, extra)
ADD <- function(id, est, cohort_id, role, outcome, description, ci_method, population, estimand = "ATO", comparator = "meropenem", arm_stats = NULL) {
  rm_tbl <<- add_row_v31(rm_tbl, id, est, cohort_id, phase = "V32", role = role, estimand = estimand, outcome = outcome,
                         description = description, comparator = comparator)
  i <- nrow(rm_tbl); set(rm_tbl, i, "ci_method", ci_method); set(rm_tbl, i, "population", population)
  if (!is.null(arm_stats)) for (k in names(arm_stats)) set(rm_tbl, i, k, arm_stats[[k]])
  if (!is.na(est$rd)) { pc <- precision_cols_v32(est$rd_ci_low, est$rd_ci_high); for (k in names(pc)) set(rm_tbl, i, k, pc[[k]]) }
}
MI2 <- function(id, cohort_id, outcome, role, description, population, method, covs = PS, B = B_SEC, estimand = "ATO") {
  dl <- D(cohort_id)
  e <- tryCatch(ato_mi_v31(dl, covs, outcome, method, B = B, label = id), error = function(err) { message("[", id, "] FAILED: ", conditionMessage(err)); NULL })
  if (is.null(e)) return(invisible(NULL))
  ADD(id, e, cohort_id, role, outcome, description, CI_MI(B), population, estimand = estimand, arm_stats = enrich_arm_stats_v32(e, dl, outcome))
  invisible(e)
}

POPS <- list(early = list(label = "Early (primary)", method = METHOD[["early"]], eligible = "early_eligible"),
             broad = list(label = "Broad (secondary)", method = METHOD[["broad"]], eligible = "eligible_2arm"))
SEC <- c(aki_stage2_d7_v31 = "AKI stage >=2 by day 7", aki_stage3_d7_v31 = "AKI stage 3 by day 7", aki_d14_v31 = "AKI by day 14",
         rrt_7d = "RRT within 7 days", rrt_30d = "RRT within 30 days", death_7d = "Death within 7 days", death_30d = "Death within 30 days",
         make30_v3 = "MAKE-30")
MAIN <- list(); POS <- list(); BAL <- list(); T1 <- list(); TRAJ <- list(); CIF <- list(); EV <- list()
creat_vals <- creat_dt[, .(subject_id, hadm_id, charttime, value = creat_mgdl)]
labs <- readRDS("output/extraction/derived_data/labs_long_v31.rds")
bun_vals <- labs[itemid == LAB_ITEMS_V31$bun & valuenum > 0 & valuenum < 300, .(subject_id, hadm_id, charttime, value = valuenum)]

for (p in names(POPS)) {
  cat("\n=== P2.1 Population", p, "===\n"); pp <- POPS[[p]]; mth <- pp$method
  e <- MI2(paste0(if (p == "early") "primary_" else "secondary_", p, "_ato"), p, "outcome_primary", if (p == "early") "PRIMARY" else "SECONDARY",
           sprintf("%s: PTZ vs meropenem, AKI by day 7, ATO, %s PS (A2 rule), MI stratified by care location.", pp$label, mth), pp$label, mth, B = B_MAIN)
  MAIN[[p]] <- e
  if (mth == "glm") MI2(paste0(p, "_ato_ridge_ps"), p, "outcome_primary", "SENSITIVITY", paste0(pp$label, ": ridge PS (not selected by the A2 rule)."), pp$label, "ridge", B = B_MAIN)
  for (oc in names(SEC)) MI2(paste0(p, "_sec_", oc), p, oc, "SECONDARY", paste0(pp$label, ": ", SEC[[oc]], "."), pp$label, mth)
  MI2(paste0(p, "_sec_persistent48"), paste0(p, "_persist48"), "aki_persistent_48h_v31", "SECONDARY", paste0(pp$label, ": persistent AKI >=48 h."), pp$label, mth)
  MI2(paste0(p, "_sec_persistent72"), paste0(p, "_persist72"), "aki_persistent_72h_v31", "SECONDARY", paste0(pp$label, ": persistent AKI >=72 h."), pp$label, mth)
  MI2(paste0(p, "_vanco_with"), paste0(p, "_vanco_with"), "outcome_primary", "STRATIFIED", paste0(pp$label, ": with vancomycin at T0."), pp$label, mth, covs = setdiff(PS, "vanco_6h"))
  MI2(paste0(p, "_vanco_without"), paste0(p, "_vanco_without"), "outcome_primary", "STRATIFIED", paste0(pp$label, ": without vancomycin at T0."), pp$label, mth, covs = setdiff(PS, "vanco_6h"))
  for (ncv in list(c("nc_a12", "early_creat_rise_12h", "creatinine rise within 12 h"), c("nc_a6", "early_creat_rise_6h", "creatinine rise within 6 h"),
                   c("nc_b", "nc_hypoglycemia_d7", "incident hypoglycemia days 1-7")))
    MI2(paste0(p, "_negative_control_", ncv[1]), paste0(p, "_", ncv[1]), ncv[2], "NEGATIVE_CONTROL", paste0(pp$label, ": negative control, ", ncv[3], "."), pp$label, mth)
  MI2(paste0(p, "_sens_ge2_creat"), paste0(p, "_ge2_creat"), "outcome_primary", "SENSITIVITY", paste0(pp$label, ": >=2 pre-T0 creatinine values."), pp$label, mth)
  MI2(paste0(p, "_sens_complete_vitals"), paste0(p, "_complete_vitals"), "outcome_primary", "SENSITIVITY", paste0(pp$label, ": vital signs measured (complete case on vitals)."), pp$label, mth)
  MI2(paste0(p, "_sens_daily_max"), paste0(p, "_dailymax"), "aki_d7_dailymax_v31", "SENSITIVITY", paste0(pp$label, ": highest creatinine per 24 h only."), pp$label, mth)
  MI2(paste0(p, "_moss_against_ptz"), pp$eligible, "outcome_moss_s3b", "SENSITIVITY", paste0(pp$label, ": missing-outcome bound against PTZ excess."), pp$label, mth)
  MI2(paste0(p, "_moss_for_ptz"), pp$eligible, "outcome_moss_s4b", "SENSITIVITY", paste0(pp$label, ": missing-outcome bound for PTZ excess."), pp$label, mth)
  # médiane + indicateurs de manquants (un seul jeu, non imputé)
  raw <- RAW(p); mi <- paste0(VARS_WITH_MISSING_V32, "_missing"); mi <- mi[vapply(mi, function(v) stats::var(raw[[v]]) > 0, logical(1))]
  mi <- setdiff(mi, paste0(unname(ICU_VITALS_V32), "_missing"))
  em <- ato_boot_v31(raw, c(PS, mi), "outcome_primary", mth, B = B_SINGLE)
  ADD(paste0(p, "_sens_median_indicators"), em, p, "SENSITIVITY", "outcome_primary", paste0(pp$label, ": median imputation + missing indicators (no MI)."),
      sprintf("Percentile bootstrap (B=%d, PS re-fitted per replicate), single dataset", B_SINGLE), pp$label,
      arm_stats = list(events_ptz = sum(raw$outcome_primary[raw$treatment == "PTZ"]), events_comp = sum(raw$outcome_primary[raw$treatment != "PTZ"])))
  # per-protocol IPCW (censure artificielle avant 48 h)
  ppan <- build_pp_panel_v31(raw, ab, creat_dt, tbl_paths)
  fwrite(merge(ppan$person, raw[, .(hadm_id, treatment)], by = "hadm_id")[, .(n = .N, pct_censored = round(100 * mean(censored), 1),
         pct_stop = round(100 * mean(censored == 1L & stop_h <= switch_h), 1), pct_switch = round(100 * mean(censored == 1L & switch_h < stop_h), 1)),
         by = treatment], TAB(paste0("per_protocol_censoring_", p, ".csv")))
  ppf <- function(d, ps_covs, outcome_col, seed, ...) pp_ipcw_boot_v31(d, ppan, ps_covs, mth, B = B_SEC, seed = seed)
  epp <- pool_estimator_mi_v31(D(p), ppf, PS, "outcome_primary", paste(p, "per-protocol"))
  epp$ps_method <- mth; epp$m <- length(IMP)
  ADD(paste0(p, "_per_protocol_ipcw"), epp, p, "SENSITIVITY", "outcome_primary", paste0(pp$label, ": per-protocol (continue >=48 h), artificial censoring + stabilized IPCW (P99), x baseline overlap weights."),
      sprintf("Bootstrap SE within each of m=%d imputations (B=%d; PS and censoring models re-fitted), Rubin's rules", length(IMP), B_SEC), pp$label, estimand = "per_protocol_ATO")
  # E-value du RR principal
  EV[[p]] <- data.table(population = pp$label, rr = e$rr, rr_lo = e$rr_ci_low, rr_hi = e$rr_ci_high,
                        evalue_point = EValue::evalues.RR(e$rr, e$rr_ci_low, e$rr_ci_high)[2, 1],
                        evalue_ci = EValue::evalues.RR(e$rr, e$rr_ci_low, e$rr_ci_high)[2, 2])
  # balance, Table 1, positivité, trajectoires, incidence cumulée
  BAL[[p]] <- compute_balance_v31(D(p), PS, e$weights_by_imputation)
  T1[[p]] <- build_table1_v31(raw, D(p), e$weights_by_imputation, PS, BAL[[p]])
  POS[[p]] <- positivity_binary_v32(D(p)[[1]], PS, mth, pp$label, comparator_label = "meropenem")
  wmean <- Reduce(`+`, e$weights_by_imputation) / length(IMP)
  TRAJ[[p]] <- rbind(trajectory_v32(raw, wmean, creat_vals, "value", "Creatinine (mg/dL)"),
                     trajectory_v32(raw, wmean, bun_vals, "value", "Urea nitrogen (mg/dL)"))[, population := pp$label]
  CIF[[pp$label]] <- cif_weighted_v32(D(p)[[1]], PS, mth, B = B_CIF)
  fwrite(BAL[[p]], TAB(paste0("balance_", p, "_v32.csv"))); fwrite(T1[[p]], TAB(paste0("table1_", p, "_v32.csv")))
  fwrite(build_overlap_population_table_v31(D(p), e$weights_by_imputation, PS), TAB(paste0("overlap_population_", p, "_v32.csv")))
}

cat("\n=== P2.2 Essai à 3 bras ===\n")
g3 <- gow_mi_v31(D("three_arm_early"), PS, "outcome_primary", B = B_SEC)
d3 <- D("three_arm_early")[[1]]; raw3 <- RAW("three_arm_early")
ess_arm <- function(k) mean(vapply(seq_along(g3$weights_by_imputation), function(i) { s <- D("three_arm_early")[[i]]$treatment == k
  w <- g3$weights_by_imputation[[i]]; sum(w[s])^2 / sum(w[s]^2) }, numeric(1)))
for (cn in names(g3$contrasts)) {
  a1 <- c(ptz_vs_mero = "PTZ", ptz_vs_cef = "PTZ", cef_vs_mero = "cefepime")[[cn]]; a0 <- c(ptz_vs_mero = "meropenem", ptz_vs_cef = "cefepime", cef_vs_mero = "meropenem")[[cn]]
  ADD(paste0("three_arm_", cn), g3$contrasts[[cn]], "three_arm_early", "THREE_ARM", "outcome_primary",
      paste0("Three-arm target trial, early population: ", a1, " minus ", a0, "; multinomial PS, generalized overlap weights."),
      sprintf("Bootstrap SE within each of m=%d imputations (B=%d, stratified by arm, multinomial PS re-fitted), Rubin's rules", length(IMP), B_SEC),
      "Three-arm (early)", estimand = "ATO_generalized", comparator = a0,
      arm_stats = list(events_ptz = sum(d3$outcome_primary[d3$treatment == a1]), events_comp = sum(d3$outcome_primary[d3$treatment == a0]),
                       ess_ptz = ess_arm(a1), ess_comp = ess_arm(a0)))
}
bal3 <- compute_balance_3arm_v31(D("three_arm_early"), PS, g3$weights_by_imputation)
fwrite(bal3$by_pair, TAB("balance_three_arm_by_pair_v32.csv")); fwrite(bal3$max, TAB("balance_three_arm_v32.csv"))
fwrite(table1_3arm_v32(raw3, D("three_arm_early"), g3$weights_by_imputation, PS, bal3$by_pair), TAB("table1_three_arm_v32.csv"))
pos3 <- positivity_3arm_v32(d3, PS); fwrite(pos3$table, TAB("positivity_three_arm_v32.csv"))
save_fig_v32(fig_gps_v32(pos3$long), FIG("efig_positivity_three_arm_gps"), 18, 8)
tern <- fig_ternary_v32(pos3$e)
if (!is.null(tern)) tryCatch({ ggplot2::ggsave(FIG("efig_positivity_three_arm_ternary.png"), tern, width = 12, height = 11, units = "cm", dpi = 300, bg = "white")
  ggplot2::ggsave(FIG("efig_positivity_three_arm_ternary.pdf"), tern, width = 12, height = 11, units = "cm", bg = "white") },
  error = function(err) message("[P2.2] ternary plot failed: ", conditionMessage(err)))

cat("\n=== P2.3 ACORN ===\n")
ab_ac <- D("acorn_aki7")
er <- ato_mi_v31(ab_ac, PS, "aki_d7_v31", "ridge", B = 1, label = "acorn balance probe")
br <- compute_balance_v31(ab_ac, PS, er$weights_by_imputation)
M_AC <- if (max(abs(br$Diff.Adj), na.rm = TRUE) > 0.10) "glm" else "ridge"
ordrow <- function(o, coh) list(weight_scheme = "overlap", n = o$n, n_ptz = n_of(coh) - cohort_registry_table()[cohort_id == coh, n_comparator],
  n_carb = cohort_registry_table()[cohort_id == coh, n_comparator], n_events = NA_integer_, weighted_risk_ptz = NA_real_, weighted_risk_carb = NA_real_,
  crude_risk_ptz = NA_real_, crude_risk_carb = NA_real_, rd = NA_real_, rd_ci_low = NA_real_, rd_ci_high = NA_real_,
  rr = o$or, rr_ci_low = o$or_lo, rr_ci_high = o$or_hi, fmi_rd = o$fmi, ess_pct = NA_real_, max_smd = NA_real_, n_boot = o$B * o$m,
  seed = SEED_V3, alpha = NA_real_, lambda = NA_real_, m = o$m, ps_method = M_AC)
CI_OR <- function(B) sprintf("Bootstrap SE of log OR within each of m=%d imputations (B=%d, PS re-fitted), Rubin's rules", length(IMP), B)
for (x in list(c("acorn_ordinal_or", "acorn", "acorn_ordinal_d14", "pre-specified (AKI at T0 excluded)"),
               c("acorn_ordinal_or_with_aki_t0", "acorn_with_aki_t0", "acorn_ordinal_d14", "including AKI at T0, as in the trial"),
               c("acorn_ordinal_or_outpatient_baseline", "acorn_ambbase", "acorn_ordinal_d14_ambbase", "baseline creatinine prioritizing an outpatient value within 365 days, as in the trial"))) {
  B <- if (x[1] == "acorn_ordinal_or") B_MAIN else B_SEC
  o <- ordinal_ow_mi_v31(D(x[2]), PS, x[3], M_AC, B = B)
  ADD(x[1], ordrow(o, x[2]), x[2], "BENCHMARK", x[3], paste0("ACORN emulation, ", x[4], ": proportional-odds OR for cefepime vs PTZ (trial direction), stored in rr columns."),
      CI_OR(B), "ACORN emulation", estimand = "ordinal_OR_cefepime_vs_PTZ", comparator = "cefepime")
  if (x[1] == "acorn_ordinal_or") fwrite(data.table(level = c("no_AKI", "stage1", "stage2", "stage3", "death"), t(o$dist)), TAB("acorn_ordinal_distribution_v32.csv"))
}
for (x in list(c("acorn_aki_d7", "acorn_aki7", "aki_d7_v31", "AKI by day 7"), c("acorn_make14", "acorn_make14", "make14", "MAKE14"))) {
  dl <- D(x[2]); e <- ato_mi_v31(dl, PS, x[3], M_AC, B = B_SEC, label = x[1])
  ADD(x[1], e, x[2], "BENCHMARK", x[3], paste0("ACORN emulation: ", x[4], ", PTZ minus cefepime."), CI_MI(B_SEC), "ACORN emulation", comparator = "cefepime",
      arm_stats = enrich_arm_stats_v32(e, dl, x[3]))
  if (x[1] == "acorn_aki_d7") { BAL$acorn <- compute_balance_v31(dl, PS, e$weights_by_imputation)
    T1$acorn <- build_table1_v31(RAW("acorn_aki7"), dl, e$weights_by_imputation, PS, BAL$acorn)
    POS$acorn <- positivity_binary_v32(dl[[1]], PS, M_AC, "ACORN emulation", comparator_label = "cefepime") }
}
fwrite(BAL$acorn, TAB("balance_acorn_v32.csv")); fwrite(T1$acorn, TAB("table1_acorn_v32.csv"))
fwrite(acorn_trial_vs_emulation_v32(RAW("acorn")), TAB("acorn_trial_vs_emulation_v32.csv"))

cat("\n=== P2.4 Figures et tableaux transversaux ===\n")
# A2 positivité
pos_tab <- rbindlist(lapply(POS, `[[`, "table")); fwrite(pos_tab, TAB("positivity_v32.csv"))
psd <- rbindlist(lapply(POS, `[[`, "data")); psd[, population := factor(population, levels = c("Early (primary)", "Broad (secondary)", "ACORN emulation"))]
save_fig_v32(fig_ps_mirror_v32(psd), FIG("fig3a_ps_distribution"), 18, 7)
save_fig_v32(fig_weights_v32(psd), FIG("efig_weights_distribution"), 18, 9)
# love plots, toutes populations
p_love <- fig_love_all_v32(BAL, bal3$max)
save_fig_v32(p_love, FIG("fig3b_love_plot"), 26, 22)
# A3 flux
reg <- cohort_registry_table(); fl <- readRDS("output/extraction/derived_data/antipseudo_flow_v31.rds")
nby <- function(id) { m <- RAW(id); c(n = nrow(m), ptz = sum(m$treatment == "PTZ"), mero = sum(m$treatment == "meropenem"), cef = sum(m$treatment == "cefepime")) }
mk <- function(id, label, x, y) { v <- nby(id); data.table(id = id, label = label, n = v[["n"]], n_ptz = v[["ptz"]], n_mero = v[["mero"]], n_cef = v[["cef"]], x = x, y = y) }
steps <- rbind(data.table(id = "stays", label = "Hospital stays with a documented first administration\nof PTZ, cefepime or a carbapenem", n = fl[step == "all", n],
                          n_ptz = NA, n_mero = NA, n_cef = NA, x = 0, y = 12),
               mk("all", "Unified cohort: first episode per patient, adult,\nbaseline creatinine, no prior RRT or ESRD", 0, 9.6),
               mk("eligible", "No AKI at T0", 0, 7.2),
               mk("acorn", "ACORN emulation: PTZ vs cefepime\n(ED or medical ICU, <=12 h, culture drawn)", 4.3, 7.2),
               mk("eligible_2arm", "PTZ or meropenem initiators", -3.2, 4.8),
               mk("three_arm_early", "THREE-ARM (early):\nPTZ, cefepime, meropenem", 3.4, 4.8),
               mk("broad", "BROAD (secondary):\nAKI by day 7 evaluable", -3.2, 2.4),
               mk("early", "EARLY (primary): T0 <=24 h after arrival,\nnew user, no IV antibiotic in prior 48 h", -3.2, 0), fill = TRUE)
edges <- data.table(from = c("stays", "all", "all", "eligible", "eligible", "eligible_2arm", "broad"),
                    to = c("all", "eligible", "acorn", "eligible_2arm", "three_arm_early", "broad", "early"),
                    excl_label = c(sprintf("Excluded (n = %s): co-initiation <=24 h, non-meropenem\ncarbapenem, later episode, age <18, no baseline\ncreatinine, prior RRT or ESRD", format(fl[step == "all", n] - n_of("all"), big.mark = ",")),
                                   sprintf("Excluded: AKI at T0 (n = %s)", format(n_of("all") - n_of("eligible"), big.mark = ",")), "", "", "",
                                   sprintf("Excluded: no evaluable\noutcome (n = %s)", format(n_of("eligible_2arm") - n_of("broad"), big.mark = ",")),
                                   sprintf("Excluded (n = %s): later initiation,\nprior anti-pseudomonal or IV antibiotic", format(n_of("broad") - n_of("early"), big.mark = ","))),
                    ex_x = c(0.5, 0.5, NA, NA, NA, -2.7, -2.7), ex_y = c(10.8, 8.4, NA, NA, NA, 3.6, 1.2))
fwrite(steps, TAB("flowchart_v32.csv"))
save_fig_v32(fig_flow_v32(steps, edges), FIG("fig1_flow"), 24, 18)
# A1 DAG
g <- build_dag_v32(); adj <- dag_adjustment_v32(g, PS)
fwrite(adj$sets, TAB("dag_adjustment_sets_v32.csv")); fwrite(dag_node_table_v32(PS), TAB("dag_node_table_v32.csv"))
save_fig_v32(plot_dag_v32(g) + ggplot2::coord_cartesian(xlim = c(-2.8, 7.6), clip = "off"), FIG("fig2_dag"), 18, 11)
DAG_NODES_V32[measured == TRUE & !node %in% c("Treatment", "AKI"), label := paste0(gsub("\n", " ", label), ":\n",
  vapply(variables, function(s) paste(strwrap(gsub(";", ",", s), 34), collapse = "\n"), character(1)))]
save_fig_v32(plot_dag_v32(build_dag_v32()) + ggplot2::coord_cartesian(xlim = c(-2.8, 7.6), clip = "off"), FIG("efig_dag_detailed"), 30, 22)
# A4 incidence cumulée
fwrite(rbindlist(lapply(names(CIF), function(n) CIF[[n]]$table[, population := n])), TAB("cumulative_incidence_v32.csv"))
fwrite(rbindlist(lapply(names(CIF), function(n) CIF[[n]]$curves[, population := n])), TAB("cumulative_incidence_curves_v32.csv"))
save_fig_v32(fig_cif_v32(CIF), FIG("fig4_cumulative_incidence"), 18, 8)
# B4 trajectoires
tr <- rbindlist(TRAJ); fwrite(tr, TAB("trajectories_v32.csv"))
save_fig_v32(fig_trajectory_v32(tr), FIG("fig6_trajectories"), 18, 20)
# E-values
fwrite(rbindlist(EV), TAB("evalues_v32.csv"))

cat("\n=== P2.5 Forest plots (B4) ===\n")
rmt <- rm_tbl
row <- function(id, lab, prim = FALSE) { r <- rmt[analysis_id == id]; data.table(label = lab, rd = r$rd, rd_lo = r$rd_ci_low, rd_hi = r$rd_ci_high, is_primary = prim) }
main_rows <- rbind(row("primary_early_ato", "Early (primary) | PTZ vs meropenem | ATO", TRUE),
                   row("secondary_broad_ato", "Broad | PTZ vs meropenem | ATO"),
                   row("three_arm_ptz_vs_mero", "Three-arm (early) | PTZ vs meropenem | generalized ATO"),
                   row("three_arm_ptz_vs_cef", "Three-arm (early) | PTZ vs cefepime | generalized ATO"),
                   row("three_arm_cef_vs_mero", "Three-arm (early) | cefepime vs meropenem | generalized ATO"),
                   row("acorn_aki_d7", "ACORN emulation | PTZ vs cefepime | ATO"))
fwrite(main_rows, TAB("forest_main_rows_v32.csv"))
save_fig_v32(fig_forest_v32(main_rows), FIG("fig5_forest_main"), 18, 7)
sens_ids <- rmt[population %in% c("Early (primary)", "Broad (secondary)") & role %in% c("SENSITIVITY", "STRATIFIED") & outcome %in% c("outcome_primary", "aki_d7_dailymax_v31", "outcome_moss_s3b", "outcome_moss_s4b")]
sens_rows <- sens_ids[, .(label = paste(population, "|", sub("^[^:]*: ", "", description)), rd, rd_lo = rd_ci_low, rd_hi = rd_ci_high, is_primary = FALSE)]
sens_rows[, label := substr(label, 1, 95)]
save_fig_v32(fig_forest_v32(sens_rows), FIG("efig_forest_sensitivity"), 22, 16)
nc_rows <- rmt[role == "NEGATIVE_CONTROL", .(label = paste(population, "|", sub("^[^:]*: negative control, ", "", description)), rd, rd_lo = rd_ci_low, rd_hi = rd_ci_high, is_primary = FALSE)]
save_fig_v32(fig_forest_v32(nc_rows, "Risk difference for the negative-control outcome, percentage points (95% CI)"), FIG("efig_forest_negative_controls"), 18, 7)
sec_rows <- rmt[role == "SECONDARY" & population %in% c("Early (primary)", "Broad (secondary)") & !grepl("_ato$", analysis_id),
                .(label = paste(population, "|", sub("^[^:]*: ", "", description)), rd, rd_lo = rd_ci_low, rd_hi = rd_ci_high, is_primary = FALSE)]
save_fig_v32(fig_forest_v32(sec_rows, "Risk difference for secondary outcomes, percentage points (95% CI)"), FIG("efig_forest_secondary"), 20, 14)

cat("\n=== P2.6 Sauvegardes, QA ===\n")
fwrite(cohort_registry_table(), TAB("n_by_analysis_v32.csv"))
saveRDS(rm_tbl, file.path(OUT, "results", "results_master.rds")); fwrite(rm_tbl, file.path(OUT, "results", "results_master.csv"))
fwrite(rm_tbl[, .(analysis_id, population, cohort_id, estimand, outcome, ps_method, m_imputations, n_boot, ci_method)], TAB("variance_methods_v32.csv"))
fwrite(rm_tbl[, .(analysis_id, population, n, n_carb, rd, rd_ci_low, rd_ci_high, ci_half_width, mde80)][!is.na(rd)], TAB("precision_v32.csv"))
build_manuscript_numbers_json_v31(rm_tbl, file.path(OUT, "results"))
pngs <- list.files(file.path(OUT, "figures"), pattern = "\\.png$", full.names = TRUE)
hashes <- tools::md5sum(pngs)
need <- list(early = c("table1_early_v32.csv", "balance_early_v32.csv"), broad = c("table1_broad_v32.csv", "balance_broad_v32.csv"),
             acorn = c("table1_acorn_v32.csv", "balance_acorn_v32.csv"), three_arm = c("table1_three_arm_v32.csv", "balance_three_arm_v32.csv", "positivity_three_arm_v32.csv"))
qa <- rbind(
  data.table(id = "QA-P2-1", check = "Each population has Table 1, balance (love plot) and positivity outputs",
             ok = all(file.exists(TAB(unlist(need)))) && all(c("Early (primary)", "Broad (secondary)", "ACORN emulation") %in% pos_tab$population),
             details = paste(names(need), collapse = ", ")),
  data.table(id = "QA-P2-2", check = "Each population's main result reports events and weighted risks by arm",
             ok = all(!is.na(rm_tbl[analysis_id %in% c("primary_early_ato", "secondary_broad_ato", "acorn_aki_d7", paste0("three_arm_", names(g3$contrasts))), events_ptz])), details = ""),
  data.table(id = "QA-P2-3", check = "DAG present; minimal adjustment sets computed; gaps between DAG sets and the PS listed",
             ok = file.exists(FIG("fig2_dag.png")) && nrow(adj$sets) > 0, details = sprintf("%d minimal sets; all covered by the PS: %s", nrow(adj$sets), all(adj$sets$all_in_ps))),
  data.table(id = "QA-P2-4", check = "No duplicated figure (MD5 of PNG files)", ok = !anyDuplicated(hashes), details = sprintf("%d PNG files", length(pngs))),
  data.table(id = "QA-P2-5", check = "Every results_master row has a CI method", ok = all(nzchar(rm_tbl$ci_method)), details = sprintf("%d rows", nrow(rm_tbl))),
  data.table(id = "QA-P2-6", check = "Every results_master row has the N of its registered cohort",
             ok = all(rm_tbl$n == cohort_registry_table()[match(rm_tbl$cohort_id, cohort_id), n]), details = ""),
  data.table(id = "QA-P2-7", check = "Single early-cohort definition for 2- and 3-arm analyses", ok = qa_early_ids_v32(master), details = ""))
qa[, status := fifelse(ok, "PASS", "FAIL")][, ok := NULL]
fwrite(qa, TAB("qa_phase2_v32.csv")); print(qa[, .(id, status, details)])
print(rm_tbl[, .(analysis_id, n, rd = round(100 * rd, 2), lo = round(100 * rd_ci_low, 2), hi = round(100 * rd_ci_high, 2), rr = round(rr, 3), mde = round(100 * mde80, 1))], nrows = 200)
cat(sprintf("PHASE 2 v3.2 TERMINÉE (%.1f min)\n", as.numeric(difftime(Sys.time(), T_START, units = "mins"))))
