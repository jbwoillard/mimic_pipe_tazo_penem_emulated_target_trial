# =============================================================================
# 83_text_formatter.R — Formatage UNIQUE des phrases de résultat (C4) : utilisé par
#   le rapport (helpers inline tracés), legends.md et key_numbers.md, pour que
#   le résumé, la conclusion et les légendes reprennent exactement les mêmes
#   chaînes (QA-COH).
# =============================================================================

.f1_v33 <- function(x, d) formatC(x, format = "f", digits = d)

#' "+5.7 percentage points [95% CI -3.0 to 14.4]"
fmt_rd_v33 <- function(r) sprintf("%+.1f percentage points [95%% CI %s to %s]", 100 * r$rd, .f1_v33(100 * r$rd_ci_low, 1), .f1_v33(100 * r$rd_ci_high, 1))
#' "1.22 [95% CI 0.89 to 1.66]" ou "not estimated (...)"
fmt_rr_v33 <- function(r) if (is.na(r$rr)) { if (!is.null(r$rr_status) && !is.na(r$rr_status)) r$rr_status else "not estimated" } else
  sprintf("%s [95%% CI %s to %s]", .f1_v33(r$rr, 2), .f1_v33(r$rr_ci_low, 2), .f1_v33(r$rr_ci_high, 2))
fmt_risk_v33 <- function(x) sprintf("%s%%", .f1_v33(100 * x, 1))

KEY_RESULTS_V33 <- c(primary_early_ato = "Early population (primary), PTZ vs meropenem",
                     secondary_broad_ato = "Broad population (secondary), PTZ vs meropenem",
                     three_arm_ptz_vs_mero = "Three-arm early trial, PTZ vs meropenem",
                     three_arm_ptz_vs_cef = "Three-arm early trial, PTZ vs cefepime",
                     acorn_aki_d7 = "ACORN emulation, AKI by day 7, PTZ vs cefepime")

#' Légendes des figures et tableaux principaux (les nombres viennent de results_master / registre / tables)
legends_v33 <- function(rmt, reg, tab) {
  R1 <- function(id) { r <- rmt[analysis_id == id]; stopifnot(nrow(r) == 1); r }
  NN <- function(id, f = "n") format(reg[cohort_id == id][[f]], big.mark = ",")
  pos <- tab("positivity_v32.csv"); cif <- tab("cumulative_incidence_v32.csv")
  cd <- function(p, k) 100 * cif[population == p & day == 7][[k]]
  c("# Figure and table legends (generated from results_master; do not edit by hand)", "",
    sprintf("**Figure 1. Cohort flow (STROBE).** From hospital stays with a first IV administration of a study drug to the early population (primary), with the number excluded at each step by arm. Unified cohort N = %s; early population N = %s (PTZ %s, meropenem %s). Dashed arrows: broad population (secondary, N = %s), three-arm early trial (N = %s) and ACORN emulation (benchmark, N = %s).",
            NN("all"), NN("early"), NN("early", "n_ptz"), NN("early", "n_comparator"), NN("broad"), NN("three_arm_early"), NN("acorn")), "",
    "**Figure 2. Directed acyclic graph.** Grey boxes: exposure and outcome. Dashed boxes and arrows: unmeasured causes (prescriber preference, beta-lactam allergy, undocumented colonization). Every arrow is an edge of the dagitty model used to derive the minimal adjustment sets.", "",
    sprintf("**Figure 3. Positivity and balance.** A, propensity-score distributions by arm (out-of-fold AUC %.2f in the early population, %.2f in the broad population). B, absolute standardized mean differences before and after overlap weighting (dashed line, 0.10).",
            pos[population == "Early (primary)", auc_out_of_fold], pos[population == "Broad (secondary)", auc_out_of_fold]), "",
    sprintf("**Figure 4. Cumulative incidence of AKI.** Overlap-weighted Aalen-Johansen estimator with death as a competing risk, censoring at discharge; shaded bands: pointwise 95%% bootstrap CI; numbers at risk (unweighted) below the axis. Difference at day 7 (PTZ minus meropenem): early %+.1f percentage points [95%% CI %s to %s]; broad %+.1f percentage points [95%% CI %s to %s].",
            cd("Early (primary)", "rd"), .f1_v33(cd("Early (primary)", "rd_lo"), 1), .f1_v33(cd("Early (primary)", "rd_hi"), 1),
            cd("Broad (secondary)", "rd"), .f1_v33(cd("Broad (secondary)", "rd_lo"), 1), .f1_v33(cd("Broad (secondary)", "rd_hi"), 1)), "",
    sprintf("**Figure 5. Main estimates.** Risk difference in AKI by day 7. %s. Multiple imputation (m = %d), bootstrap within imputations (B = %d per imputation for the early population), Rubin's rules.",
            paste(sprintf("%s: %s", unname(KEY_RESULTS_V33), vapply(names(KEY_RESULTS_V33), function(i) fmt_rd_v33(R1(i)), "")), collapse = ". "),
            R1("primary_early_ato")$m_imputations, R1("primary_early_ato")$n_boot / R1("primary_early_ato")$m_imputations), "",
    "**Figure 6. Creatinine, urea nitrogen and urea-to-creatinine ratio, days 0 to 7.** Overlap-weighted daily means by arm (top) and PTZ minus meropenem difference with 95% CI (bottom); numbers of patients with a value shown per day. Later days are selected (patients still in hospital and tested); a sensitivity analysis restricted to patients with a creatinine value on days 5-7 is in the supplement.", "",
    sprintf("**Table 2. Primary estimate under alternative propensity-score specifications (sensitivity).** Early population. Primary: %s. Augmented overlap weighting: %s. Cross-fitted AIPW (estimand ATE): %s.",
            fmt_rd_v33(R1("primary_early_ato")), fmt_rd_v33(R1("early_sens_ato_augmented")), fmt_rd_v33(R1("early_sens_aipw_crossfit"))), "",
    sprintf("**eFigure (exploratory). Heterogeneity of treatment effect, broad population.** Out-of-fold predicted effects (causal forests, thirteen predictors, delay to T0 removed) and doubly robust effect by tertile: lower %s, middle %s, upper %s. Exploratory analysis in the secondary population; not performed in the early population (insufficient meropenem initiators).",
            fmt_rd_v33(R1("hte_tertile_lower")), fmt_rd_v33(R1("hte_tertile_middle")), fmt_rd_v33(R1("hte_tertile_upper"))), "",
    sprintf("**Table 1. Baseline characteristics of the early population (primary).** N = %s (PTZ %s, meropenem %s). Unweighted: observed values before imputation; overlap-weighted values and SMDs averaged over %d imputations.",
            NN("early"), NN("early", "n_ptz"), NN("early", "n_comparator"), R1("primary_early_ato")$m_imputations), "",
    "**Table 3. Secondary population (broad), three-arm early trial and ACORN benchmark.** First arm = first-named drug of each contrast; ACORN odds ratios are for cefepime vs PTZ (trial direction); the published ACORN estimate is shown for comparison.", "",
    "**eFigures.** eFigure 1, flow by analytic population; eFigure 2, weights; eFigure 3, bootstrap distribution of the primary estimate; eFigures 4-6, forest plots of secondary outcomes, sensitivity analyses and negative controls; eFigure 7, trajectories restricted to patients with creatinine on days 5-7; eFigures 8-9, three-arm positivity; eFigures 10-14, exploratory heterogeneity (broad population).")
}
