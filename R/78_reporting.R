# =============================================================================
# 78_reporting.R — Restructuration pour publication (docs/v3.2_plan.md §4)
#   C1 déviations au protocole, C2 checklist TARGET, C3 paquet manuscrit,
#   QA « rapport intégralement en anglais »
# =============================================================================

#' C1 — Déviations au protocole (transparence obligatoire). Sources : plans et
#' changelogs v3, v3.1, v3.2 et HTE (docs/).
protocol_deviations_v32 <- function() data.table(
  id = sprintf("D%02d", 1:10),
  version = c("v3.1", "v3.1", "v3.1", "v3.1", "v3.2", "v3.2", "v3.2", "v3.2", "v3.2", "v3.2"),
  deviation = c(
    "Early empirical population added (co-primary) after the v3 results",
    "Strict KDIGO outcome (rolling 48 h window including pre-T0 values) replaced the v3 definition",
    "Per-protocol restriction to >=48 h of exposure replaced by artificial censoring with inverse probability of censoring weights",
    "ACORN benchmark redefined to emulate the trial eligibility, timing and ordinal day-14 outcome",
    "Cohort rebuilt from documented first administrations; comparator restricted to meropenem; patients exposed to cefepime before the index drug no longer eligible for the two-arm comparisons",
    "Early empirical population promoted to the single primary analysis; broad population becomes secondary",
    "Pre-T0 creatinine trajectory added to the propensity score",
    "Emergency-department vital signs added and multiple imputation stratified by care location",
    "Propensity-score rule (unpenalized model if the ridge model leaves a weighted SMD above 0.10) kept for the primary analysis despite instability; ridge model reported as sensitivity",
    "Delay from arrival to T0 removed from the effect modifiers of the heterogeneity analysis"),
  reason = c(
    "The broad cohort mixed empirical initiations with late escalations; the target trial is empirical initiation",
    "The v3 absolute criterion was referenced to the 7-day baseline when no value existed in the prior 48 h",
    "The v3 restriction conditioned on a post-baseline variable with differential selection",
    "The v3 benchmark applied the carbapenem pipeline to cefepime without the trial criteria",
    "The inherited cohort used any eMAR event as T0 (including doses not given) and included patients escalating from cefepime (not new users)",
    "Consistency between the stated target trial and the analysis hierarchy",
    "The early creatinine-rise negative control was non-null in the broad population",
    "Vital signs were missing by design outside the ICU and imputed from ICU patients",
    "Pre-specified rule; changing it after seeing the results would be a further deviation",
    "The delay is the main confounding axis identified by the design analyses"),
  where_reported = c("Early empirical population section; reconciliation table", "Outcome definition; v3.1 changelog",
                     "Sensitivity analyses", "ACORN benchmark section", "Cohort construction; flow diagram; reconciliation table",
                     "Target trial protocol; summary of findings", "Methods; waterfall table", "Missing data section",
                     "Primary result; positivity section", "Heterogeneity section"))

#' C2 — Checklist TARGET (Cashin et al., JAMA 2025; structure vérifiée sur jamanetwork.com/journals/jama/fullarticle/2837724).
#' Localisation = titres de section du rapport v3.2.
target_checklist_v32 <- function() data.table(
  item = c("1a", "1b", "1c", "2", "3", "4", "5", "6a", "6b", "6c", "6d", "6e", "6f", "6g", "6h",
           "7a", "7b", "7c", "7d", "7e", "7f", "7g.i", "7g.ii", "7h.i", "7h.ii",
           "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21"),
  section = c(rep("Abstract", 3), rep("Introduction", 3), "Methods", rep("Methods: target trial specification", 8),
              rep("Methods: target trial emulation", 10), rep("Results", 7), rep("Discussion", 2), rep("Other information", 5)),
  topic = c("Identification as a target trial emulation", "Data sources", "Assumptions, methods, findings, conclusions",
            "Background and knowledge gap", "Causal question", "Rationale for emulation", "Data sources and setting",
            "Eligibility criteria", "Treatment strategies", "Assignment procedures", "Start and end of follow-up", "Outcomes",
            "Causal contrasts", "Identifying assumptions", "Analysis plan",
            "Eligibility (emulation)", "Treatment strategies (emulation)", "Assignment (emulation)", "Follow-up (emulation)", "Outcomes (emulation)",
            "Causal contrasts (emulation)", "Identifying assumptions (emulation)", "Variables for the assumptions", "Analysis and missing data",
            "Sensitivity analyses", "Participant flow", "Baseline characteristics", "Follow-up", "Missing data", "Outcome frequency",
            "Effect estimates with precision", "Sensitivity and additional analyses", "Interpretation", "Limitations",
            "Ethics approval", "Protocol registration", "Data and code availability", "Funding", "Conflicts of interest"),
  location = c("Title; Summary of findings (manuscript abstract to be written from it)", "Summary of findings; Methods", "Summary of findings; Conclusion",
               "Manuscript introduction (to be written by the authors)", "Target trial protocol", "Target trial protocol; Protocol deviations",
               "Methods: data sources and cohort construction", rep("Target trial protocol (table)", 8),
               "Methods: cohort construction; Cohort flow", "Target trial protocol; Methods", "Target trial protocol (no randomization: confounding adjustment)",
               "Target trial protocol", "Outcome definition", "Target trial protocol; Primary result", "Directed acyclic graph; Positivity",
               "Directed acyclic graph (node table)", "Methods: estimation, missing data, variance", "Sensitivity analyses; Protocol deviations",
               "Cohort flow (Figure 1)", "Table 1 (each population)", "Cumulative incidence; Methods", "Missing data", "Primary result; Secondary outcomes",
               "Summary of findings; Primary result (absolute and relative)", "Sensitivity analyses; Benchmark; Heterogeneity", "Conclusion", "Limitations",
               "To be completed by the authors (MIMIC-IV credentialed access)", "To be completed by the authors", "Software; Methods (code in the project repository; data via PhysioNet)",
               "To be completed by the authors", "To be completed by the authors"),
  status = c("partial", "reported", "reported", "author", "reported", "reported", "reported", rep("reported", 8), rep("reported", 10),
             rep("reported", 7), "reported", "reported", "author", "author", "reported", "author", "author"))

#' Effet publié d'ACORN (Qian et al., JAMA 2023;330:1557-1567, résultat principal :
#' plus haut stade d'AKI ou décès à J14, OR céfépime vs PTZ) — repris du tableau
#' de référence v3.1 (acorn_benchmark_v31.csv)
ACORN_PUBLISHED_EFFECT_V32 <- data.table(analysis = "ACORN trial (published)", n = 2511L, or_cefepime_vs_ptz = 0.95, or_ci_low = 0.80, or_ci_high = 1.13)

#' Tableaux manuscrit (C3), lus dans results_master
fmt_res_v32 <- function(r) data.table(
  N = format(r$n, big.mark = ","),
  `Events / N, first arm` = ifelse(is.na(r$events_ptz), "", sprintf("%d / %d", r$events_ptz, r$n_ptz)),
  `Events / N, second arm` = ifelse(is.na(r$events_comp), "", sprintf("%d / %d", r$events_comp, r$n_carb)),
  `Weighted risk, first arm` = ifelse(is.na(r$weighted_risk_ptz), "", sprintf("%.1f%%", 100 * r$weighted_risk_ptz)),
  `Weighted risk, second arm` = ifelse(is.na(r$weighted_risk_carb), "", sprintf("%.1f%%", 100 * r$weighted_risk_carb)),
  `RD, pp (95% CI)` = ifelse(is.na(r$rd), "", sprintf("%+.1f (%.1f to %.1f)", 100 * r$rd, 100 * r$rd_ci_low, 100 * r$rd_ci_high)),
  `RR or OR (95% CI)` = ifelse(is.na(r$rr) | !is.finite(r$rr_ci_low) | !is.finite(r$rr_ci_high), "not estimable",
                               sprintf("%s %.2f (%.2f to %.2f)", ifelse(grepl("ordinal", r$analysis_id), "OR", "RR"), r$rr, r$rr_ci_low, r$rr_ci_high)),
  `ESS, first / second arm` = ifelse(is.na(r$ess_ptz), "", sprintf("%.0f / %.0f", r$ess_ptz, r$ess_comp)),
  `MDE80, pp` = ifelse(is.na(r$mde80), "", sprintf("%.1f", 100 * r$mde80)))

manuscript_results_table_v32 <- function(rmt, ids, labels) {
  r <- rmt[match(ids, analysis_id)]
  if (anyNA(r$analysis_id)) stop("[78] analysis_id absent : ", paste(ids[is.na(r$analysis_id)], collapse = ", "))
  cbind(data.table(Analysis = labels), fmt_res_v32(r))
}

TABLE2_IDS_V32 <- c(primary_early_ato = "Early (PRIMARY): PTZ vs meropenem, overlap weights",
                    early_ato_ridge_ps = "Early: ridge propensity score (sensitivity)",
                    early_per_protocol_ipcw = "Early: per-protocol, artificial censoring + IPCW",
                    early_sec_aki_stage2_d7_v31 = "Early: AKI stage >=2 by day 7",
                    early_sec_rrt_30d = "Early: RRT by day 30", early_sec_death_30d = "Early: death by day 30",
                    secondary_broad_ato = "Broad (secondary): PTZ vs meropenem, overlap weights",
                    broad_per_protocol_ipcw = "Broad: per-protocol, artificial censoring + IPCW")
TABLE3_IDS_V32 <- c(three_arm_ptz_vs_mero = "Three-arm (early): PTZ vs meropenem", three_arm_ptz_vs_cef = "Three-arm (early): PTZ vs cefepime",
                    three_arm_cef_vs_mero = "Three-arm (early): cefepime vs meropenem",
                    acorn_ordinal_or = "ACORN emulation: ordinal day-14 outcome, OR cefepime vs PTZ (pre-specified, no AKI at T0)",
                    acorn_ordinal_or_with_aki_t0 = "ACORN emulation: ordinal outcome including AKI at T0",
                    acorn_ordinal_or_outpatient_baseline = "ACORN emulation: ordinal outcome, outpatient-priority baseline creatinine",
                    acorn_aki_d7 = "ACORN emulation: AKI by day 7, PTZ vs cefepime", acorn_make14 = "ACORN emulation: MAKE14, PTZ vs cefepime")

#' Export .docx d'un data.table (flextable + officer)
save_docx_table_v32 <- function(dt, path, caption) {
  ft <- flextable::flextable(as.data.frame(dt))
  ft <- flextable::fontsize(ft, size = 8, part = "all"); ft <- flextable::autofit(ft)
  doc <- officer::read_docx(); doc <- officer::body_add_par(doc, caption, style = "Normal"); doc <- flextable::body_add_flextable(doc, ft)
  print(doc, target = path)
}

#' QA « rapport intégralement en anglais » : mots français fréquents et
#' lettres accentuées dans le texte visible du HTML (hors code et sorties R)
qa_english_v32 <- function(html_path) {
  s <- paste(readLines(html_path, warn = FALSE), collapse = "\n")
  s <- gsub("(?s)<script.*?</script>|<style.*?</style>|<pre.*?</pre>|<code.*?</code>", " ", s, perl = TRUE)
  s <- sub("(?s)<div id=\"refs\".*$", " ", s, perl = TRUE)   # bibliographie : noms d'auteurs accentués légitimes
  t <- gsub("<[^>]+>", " ", s); t <- gsub("&[a-z#0-9]+;", " ", t)
  fr <- c("\\bles\\b", "\\bdes\\b", "\\bune\\b", "\\bpour\\b", "\\bavec\\b", "\\bdans\\b", "\\best\\b", "\\bsont\\b", "\\bcohorte\\b",
          "\\bdonn[ée]es\\b", "\\bpatients? de\\b", "\\bmanquants?\\b", "\\bn'est\\b", "\\bpas de\\b")
  hits <- unlist(lapply(fr, function(p) regmatches(t, gregexpr(p, t, ignore.case = TRUE, perl = TRUE))[[1]]))
  acc <- regmatches(t, gregexpr("[A-Za-z]*[éèêàùçôîâ][A-Za-z]*", t, perl = TRUE))[[1]]
  acc <- setdiff(unique(acc), c("Aalen", "°C"))
  list(ok = length(hits) == 0 && length(acc) == 0, french_words = unique(hits), accented = acc)
}
