# =============================================================================
# 05_heterogeneity.R — Heterogeneity of treatment effect: cross-fitted causal forests,
#   AIPW scores, tertiles of predicted effect, temporal validation
#   (derivation <= 2019, test 2020-2022). -> output/analysis/hte/
# Usage: Rscript scripts/05_heterogeneity.R
# =============================================================================
options(warn = 1, scipen = 10, v31.cores = 8L)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(glmnet) })
PKG_AVAILABLE <- function(p) isTRUE(requireNamespace(p, quietly = TRUE))
invisible(lapply(c("grf", "xgboost", "patchwork", "sandwich", "lmtest", "survey", "jsonlite"), loadNamespace))
FAST <- Sys.getenv("HTE_FAST") == "1"
OUT <- "output/analysis"; HOUT <- file.path(OUT, "hte"); dir.create(file.path(HOUT, "derived"), recursive = TRUE, showWarnings = FALSE)
H <- function(f) file.path(HOUT, f); TAB <- function(f) file.path(OUT, "tables", f)
T0 <- Sys.time()
for (f in c("R/01_paths_and_io.R", "R/02_codebooks.R", "R/20_bootstrap.R", "R/30_results_master.R", "R/50_master_cohort.R",
            "R/51_analysis_cohort.R", "R/52_estimation.R", "R/53_extract.R", "R/56_results_tables.R",
            "R/61_hte.R", "R/62_hte_figures.R", "R/71_registry.R", "R/75_positivity_flow.R", "R/77_tables_figs.R")) source(f)
# (ii) modificateurs : 14 prédicteurs v3.1 moins le délai arrivée -> T0 (le service n'en faisait pas partie)
HTE_PREDICTORS_V31 <- setdiff(HTE_PREDICTORS_V31, "time_admission_to_T0_h")
cfg <- hte_cfg_v31(FAST)

master <- setup_registry_v32(readRDS(file.path(OUT, "derived_data", "master_v32.rds")))$master
IMP <- readRDS(file.path(OUT, "derived_data", "mi_C_stratified_m20.rds"))$imp
add_vars <- function(d, ids) { d <- copy(setkey(d, hadm_id)[J(ids)])
  d[, other_nephrotoxin_any := as.integer(direct_nephrotoxin_pre_t0 == 1L | calcineurin_inhibitor_pre_t0 == 1L | aminoglycoside_pre_t0 == 1L | nsaid_pre_t0 == 1L)]; d }
broad_ids <- .COHORT_REGISTRY_V31$ids[["broad"]]
d <- add_vars(IMP[[1]], broad_ids); PS <- PS_COVARIATES_V32
cat("HTE v3.2, population large : N =", nrow(d), "| prédicteurs :", length(HTE_PREDICTORS_V31), "\n")

cat("=== H.1 Cross-fitting 5 plis ===\n")
cf <- run_hte_crossfit_v31(d, PS, cfg); pp <- cf$per_patient
tz <- assign_tertiles_v31(pp$tau); pp[, tertile := tz$tertile]
for (k in c("tau_rlearner", "tau_tlearner")) pp[, (paste0("tertile_", sub("tau_", "", k))) := assign_tertiles_v31(get(k))$tertile]
saveRDS(list(pp = pp, folds = cf$folds), file.path(HOUT, "derived", "hte_per_patient_v32.rds"))
tt <- dr_by_group_v31(pp, "tertile")[order(tertile)]; wt <- wald_equality_v31(tt); sv <- svy_interaction_v31(pp)
tt[, `:=`(p_interaction_wald2df = wt$p, p_interaction_svy_wald = sv$p_wald)]
overall <- dr_by_group_v31(copy(pp)[, all := "overall"], "all")
met <- hte_metrics_v31(cbind(pp, hte_X14(d[match(pp$hadm_id, d$hadm_id)][, baseline_aki_risk := pp$baseline_aki_risk])[, setdiff(HTE_PREDICTORS_V31, "baseline_aki_risk")]), cfg$seed + 5)
met <- rbind(met, data.table(metric = c("interaction_wald_2df", "interaction_svyglm_wald"), estimate = c(wt$chi2, NA), p = c(wt$p, sv$p_wald)), fill = TRUE)
fwrite(tt, H("hte_tertile_effects_v32.csv")); fwrite(met, H("hte_metrics_v32.csv")); fwrite(overall, H("hte_overall_v32.csv"))
fwrite(data.table(variable = names(cf$vi), importance = cf$vi)[order(-importance)], H("hte_variable_importance_v32.csv"))
print(tt[, .(tertile, n, rd = round(100 * rd, 2), lo = round(100 * rd_lo, 2), hi = round(100 * rd_hi, 2))]); cat("P interaction =", signif(wt$p, 3), "\n")

cat("=== H.2 Validation par période rééquilibrée (<= 2019 -> 2020-2022) ===\n")
d[, approx_admit_year := master[J(d$hadm_id)]$approx_admit_year]
per <- run_hte_period_v31(d, PS, cfg, train_max_year = 2019, test_min_year = 2020); ppp <- per$per_patient
tp <- dr_by_group_v31(ppp, "tertile")[order(tertile)]; wtp <- wald_equality_v31(tp); tp[, p_interaction_wald2df := wtp$p]
fwrite(tp, H("hte_period_tertile_effects_v32.csv"))
fwrite(data.table(rule = "derivation: approximate admission year <= 2019; test: 2020-2022", n_derivation = per$n_derivation, n_test = per$n_test,
                  test_n_ptz = sum(ppp$treatment == "PTZ"), test_n_meropenem = sum(ppp$treatment != "PTZ"),
                  test_events_ptz = sum(ppp$y[ppp$treatment == "PTZ"]), test_events_meropenem = sum(ppp$y[ppp$treatment != "PTZ"]),
                  exploratory = sum(ppp$treatment != "PTZ") < 150), H("hte_period_summary_v32.csv"))

cat("=== H.3 (i) Population précoce : faisabilité ===\n")
rm_tbl <- readRDS(file.path(OUT, "results", "results_master.rds"))
pe <- rm_tbl[analysis_id == "primary_early_ato"]
se_all <- (pe$rd_ci_high - pe$rd_ci_low) / (2 * 1.96)
fwrite(data.table(population = "Early (primary)", n = pe$n, n_meropenem = pe$n_carb, expected_meropenem_per_tertile = round(pe$n_carb / 3),
                  se_overall = se_all, expected_se_per_tertile = se_all * sqrt(3),
                  expected_ci_half_width_per_tertile = 1.96 * se_all * sqrt(3), mde80_per_tertile = 2.8 * se_all * sqrt(3),
                  decision = "Not performed: per-tertile estimates would be uninformative (see expected half-width and MDE)"),
       H("hte_early_feasibility_v32.csv"))

cat("=== H.4 Figures (PDP ordonnés, importance unique) ===\n")
ato <- rm_tbl[analysis_id == "secondary_broad_ato", rd]
save_hte_fig_v31(fig_hte_tertiles_v31(pp, tt, ato, wt$p), H("fig_hte_tertiles_v32"))
save_hte_fig_v31(fig_hte_tertiles_v31(ppp, tp, ato, wtp$p), H("efig_hte_tertiles_period_v32"))
writeLines(hte_legend_text_v31(tt, nrow(pp), wt$p, sv$p_wald, "crossfit", thresholds = tz$thresholds), H("fig_hte_legend_v32.txt"))
fmt_num <- function(x) format(signif(x, 3), scientific = FALSE, trim = TRUE, big.mark = ",")
vi <- cf$vi; top5 <- names(sort(vi, decreasing = TRUE))[1:5]
xx <- cbind(pp[, .(tau)], d[match(pp$hadm_id, d$hadm_id), setdiff(top5, "baseline_aki_risk"), with = FALSE], baseline_aki_risk = pp$baseline_aki_risk)
dep <- rbindlist(lapply(top5, function(v) {
  z <- if (v == "gender") as.numeric(toupper(xx[[v]]) == "F") else as.numeric(xx[[v]])
  if (length(unique(z)) <= 2) { lev_id <- as.integer(z) + 1L; labs <- if (v == "gender") c("Male", "Female") else c("No", "Yes") }
  else { q <- unique(stats::quantile(z, 0:5 / 5)); lev_id <- findInterval(z, q, rightmost.closed = TRUE, all.inside = TRUE)
         labs <- sprintf("%s to %s", fmt_num(q[-length(q)]), fmt_num(q[-1])) }
  data.table(level_id = lev_id, tau = xx$tau)[, .(mean_tau = mean(tau), lo = mean(tau) - 1.96 * stats::sd(tau) / sqrt(.N),
                                                  hi = mean(tau) + 1.96 * stats::sd(tau) / sqrt(.N), n = .N), by = level_id][order(level_id)][
    , `:=`(variable = if (v %in% names(VAR_LABELS_V31)) VAR_LABELS_V31[[v]] else "Predicted baseline AKI risk", level = labs[level_id])]
}))
dep[, level := factor(level, levels = unique(level))]   # ordre = ordre numérique des classes (level_id croissant)
fwrite(dep, H("hte_dependence_top5_v32.csv"))
save_fig_v32(fig_hte_dependence_v31(dep) + ggplot2::scale_x_discrete(drop = TRUE), H("efig_hte_dependence_v32"), 30, 9)
save_fig_v32(fig_hte_importance_v31(vi), H("efig_hte_importance_v32"), 14, 9)
cal <- rbind(tt[, .(grouping = "tertiles", group = as.character(tertile), mean_pred_tau, rd, rd_lo, rd_hi)],
             { pp[, quintile := factor(cut(rank(tau, ties.method = "first"), breaks = quantile(seq_len(.N), 0:5 / 5), include.lowest = TRUE, labels = FALSE))]
               dr_by_group_v31(pp, "quintile")[, .(grouping = "quintiles", group = as.character(quintile), mean_pred_tau, rd, rd_lo, rd_hi)] })
save_fig_v32(fig_hte_calibration_v31(cal), H("efig_hte_calibration_v32"), 12, 9)

cat("=== H.5 QA HTE v3.2 ===\n")
ord_ok <- dep[, { lb <- level_id; all(diff(lb) > 0) }, by = variable]$V1
pngs <- c(list.files(file.path(OUT, "figures"), "\\.png$", full.names = TRUE), list.files(HOUT, "\\.png$", full.names = TRUE))
qa <- data.table(id = c("QA-H1", "QA-H2", "QA-H3", "QA-H4", "QA-H5"),
                 check = c("No leakage: each patient predicted once by models trained without them",
                           "N per tertile = N/3 +/- 1", "Dependence-plot levels ordered numerically (no alphabetical ordering)",
                           "No duplicated figure (MD5 of all PNG files, v3.2 figures and HTE)",
                           "Delay from arrival to T0 excluded from the effect modifiers"),
                 ok = c(all(cf$leak$overlap_ids == 0) && per$overlap_ids == 0 && !anyDuplicated(pp$hadm_id),
                        max(abs(tt$n - nrow(pp) / 3)) <= 1, all(ord_ok), !anyDuplicated(tools::md5sum(pngs)),
                        !"time_admission_to_T0_h" %in% HTE_PREDICTORS_V31),
                 details = c("", paste(tt$n, collapse = " / "), "", sprintf("%d PNG files", length(pngs)), paste(HTE_PREDICTORS_V31, collapse = ", ")))
qa[, status := fifelse(ok, "PASS", "FAIL")][, ok := NULL]
fwrite(qa, H("qa_hte_v32.csv")); print(qa[, .(id, status)])
for (f in c("hte_tertile_effects_v32.csv", "hte_metrics_v32.csv", "hte_period_tertile_effects_v32.csv", "hte_period_summary_v32.csv",
            "hte_early_feasibility_v32.csv", "hte_variable_importance_v32.csv", "qa_hte_v32.csv")) file.copy(H(f), TAB(f), overwrite = TRUE)
cat(sprintf("HTE v3.2 TERMINÉ (%.1f min)\n", as.numeric(difftime(Sys.time(), T0, units = "mins"))))
