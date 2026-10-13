# =============================================================================
# 56_results_tables.R — results_master v3.1, Table 1 publication,
#   population d'overlap, flowchart, sessionInfo complet (A1, A3, A5, A10)
# =============================================================================

#' Ajoute une ligne à results_master en traçant la cohorte d'analyse (A1)
add_row_v31 <- function(master, analysis_id, est, cohort_id, phase = "A31", role, estimand = "ATO",
                        outcome, description, comparator = "carbapenem", notes = "") {
  if (!is.null(cohort_id) && est$n != n_of(cohort_id))
    stop("[56_v31] N de l'estimation (", est$n, ") != N de la cohorte ", cohort_id, " (", n_of(cohort_id), ")")
  master <- add_result_row(master, analysis_id, est, comparator = comparator, phase = phase, role = role,
                           estimand = estimand, outcome = outcome, description = description, notes = notes)
  i <- nrow(master)
  set(master, i, "cohort_id", cohort_id)
  set(master, i, "ps_method", if (is.null(est$ps_method)) NA_character_ else est$ps_method)
  set(master, i, "m_imputations", if (is.null(est$m)) NA_integer_ else as.integer(est$m))
  set(master, i, "fmi_rd", if (is.null(est$fmi_rd)) NA_real_ else est$fmi_rd)
  master
}

# -----------------------------------------------------------------------
# Libellés lisibles (A10)
# -----------------------------------------------------------------------
VAR_LABELS_V31 <- c(
  age = "Age, years", gender = "Female sex", baseline_creat_min7d = "Baseline creatinine (7-day minimum), mg/dL",
  bun_baseline = "Blood urea nitrogen, mg/dL", ckd_priorhosp = "Chronic kidney disease (prior stays)",
  diabetes_priorhosp = "Diabetes (prior stays)", heart_failure_priorhosp = "Heart failure (prior stays)",
  liver_disease_priorhosp = "Liver disease (prior stays)", cancer_priorhosp = "Cancer (prior stays)",
  charlson_score_priorhosp = "Charlson score (prior stays)", hematologic_malignancy = "Hematologic malignancy (prior stays)",
  transplant_history = "Transplant history (prior stays)", vasopressor_pre_t0 = "Vasopressor before T0",
  vent_pre_t0 = "Mechanical ventilation before T0", map_mean = "Mean arterial pressure, mmHg (24 h mean)",
  hr_mean = "Heart rate, /min (24 h mean)", temp_max = "Temperature, °C (24 h max)",
  spo2_min = "SpO2, % (24 h min)", gcs_min = "Glasgow Coma Scale (24 h min)",
  sepsis_composite_v3 = "Sepsis organ-dysfunction composite", icu_before_t0 = "ICU stay before T0",
  lactate_max = "Lactate, mmol/L (24 h max)", neutropenia_pre_t0 = "Neutropenia (ANC < 0.5 G/L)",
  platelets_min = "Platelets, G/L (24 h min)", bilirubin_max = "Total bilirubin, mg/dL (24 h max)",
  weight_kg_v31 = "Weight, kg", direct_nephrotoxin_pre_t0 = "Direct nephrotoxin (48 h)",
  calcineurin_inhibitor_pre_t0 = "Calcineurin inhibitor (48 h)", loop_diuretic_pre_t0 = "Loop diuretic (48 h)",
  acei_arb_pre_t0 = "ACE inhibitor / ARB (48 h)", nsaid_pre_t0 = "NSAID (48 h)",
  aminoglycoside_pre_t0 = "Aminoglycoside (48 h)", vanco_6h = "Vancomycin at T0 (±6 h)",
  positive_culture_pre_t0_storetime = "Positive culture result available at T0",
  prior_antibiotic = "Other study antibiotic in prior 48 h", esbl_resistant_history_12mo = "ESBL / 3GC resistance (12 months)",
  site_blood = "Culture drawn before T0: blood", site_urine = "Culture drawn before T0: urine",
  site_respiratory = "Culture drawn before T0: respiratory", site_abdominal = "Culture drawn before T0: abdominal",
  site_skin_soft_tissue_bone = "Culture drawn before T0: skin / soft tissue / bone",
  site_other = "Culture drawn before T0: other site", time_admission_to_T0_h = "Arrival (ED or admission) to T0, hours",
  anchor_year_group = "Anchor year group", service_group = "Care unit at T0"
)
lab_v31 <- function(v) ifelse(v %in% names(VAR_LABELS_V31), VAR_LABELS_V31[v], v)

.wq <- function(x, w, p) { o <- order(x); x <- x[o]; w <- w[o]; cw <- cumsum(w) / sum(w); x[which(cw >= p)[1]] }

# -----------------------------------------------------------------------
# Table 1 publication (A3, A5, A10) : brut sur valeurs OBSERVÉES (avant
# imputation, % manquants réels), pondéré sur les jeux imputés (moyenne des m
# imputations), SMD brut et pondéré lus dans l'objet de balance cobalt UNIQUE.
# -----------------------------------------------------------------------
build_table1_v31 <- function(raw, imp_list, weights_list, ps_covs, bal) {
  trt <- raw$treatment == "PTZ"
  rows <- list()
  for (v in ps_covs) {
    x <- raw[[v]]
    miss <- sprintf("%.1f / %.1f", 100 * mean(is.na(x[trt])), 100 * mean(is.na(x[!trt])))
    is_cat <- v != "gender" && (is.character(x) || is.factor(x))
    is_bin <- !is_cat && (v == "gender" || all(stats::na.omit(unique(x)) %in% c(0, 1)))
    if (v == "gender") x <- as.integer(toupper(x) == "F")
    if (is_cat) {
      for (l in sort(unique(stats::na.omit(as.character(x))))) {
        tm <- paste0(v, "_", l)
        b <- bal[bal$term == tm]
        wv <- mean(vapply(seq_along(imp_list), function(i) {
          d <- imp_list[[i]]; w <- weights_list[[i]]; z <- d[[v]] == l; t1 <- d$treatment == "PTZ"
          c(sum(w[t1] * z[t1]) / sum(w[t1]))}, numeric(1)))
        wc <- mean(vapply(seq_along(imp_list), function(i) {
          d <- imp_list[[i]]; w <- weights_list[[i]]; z <- d[[v]] == l; t0 <- d$treatment != "PTZ"
          sum(w[t0] * z[t0]) / sum(w[t0])}, numeric(1)))
        rows[[length(rows) + 1]] <- data.table(
          variable = v, term = tm, characteristic = paste0(lab_v31(v), ": ", l), type = "n (%)",
          raw_ptz = sprintf("%d (%.1f)", sum(x[trt] == l, na.rm = TRUE), 100 * mean(x[trt] == l, na.rm = TRUE)),
          raw_carb = sprintf("%d (%.1f)", sum(x[!trt] == l, na.rm = TRUE), 100 * mean(x[!trt] == l, na.rm = TRUE)),
          smd_raw = b$Diff.Un, w_ptz = sprintf("%.1f%%", 100 * wv), w_carb = sprintf("%.1f%%", 100 * wc),
          smd_weighted = b$Diff.Adj, ks_weighted = b$KS.Adj, pct_missing_ptz_carb = miss)
      }
      next
    }
    b <- bal[bal$term == v]
    wstat <- function(arm) {
      vals <- vapply(seq_along(imp_list), function(i) {
        d <- imp_list[[i]]; w <- weights_list[[i]]; s <- if (arm) d$treatment == "PTZ" else d$treatment != "PTZ"
        z <- if (v == "gender") as.numeric(toupper(d[[v]]) == "F") else as.numeric(d[[v]])
        if (is_bin) sum(w[s] * z[s]) / sum(w[s]) else c(.wq(z[s], w[s], .5))
      }, numeric(1))
      if (is_bin) sprintf("%.1f%%", 100 * mean(vals)) else {
        q <- vapply(seq_along(imp_list), function(i) {
          d <- imp_list[[i]]; w <- weights_list[[i]]; s <- if (arm) d$treatment == "PTZ" else d$treatment != "PTZ"
          z <- as.numeric(d[[v]]); c(.wq(z[s], w[s], .25), .wq(z[s], w[s], .75))
        }, numeric(2))
        sprintf("%.1f [%.1f; %.1f]", mean(vals), mean(q[1, ]), mean(q[2, ]))
      }
    }
    rawstat <- function(z) if (is_bin) sprintf("%d (%.1f)", sum(z, na.rm = TRUE), 100 * mean(z, na.rm = TRUE)) else {
      q <- stats::quantile(z, c(.5, .25, .75), na.rm = TRUE); sprintf("%.1f [%.1f; %.1f]", q[1], q[2], q[3]) }
    rows[[length(rows) + 1]] <- data.table(
      variable = v, term = v, characteristic = lab_v31(v),
      type = if (is_bin) "n (%)" else "median [IQR]",
      raw_ptz = rawstat(as.numeric(x[trt])), raw_carb = rawstat(as.numeric(x[!trt])),
      smd_raw = b$Diff.Un, w_ptz = wstat(TRUE), w_carb = wstat(FALSE),
      smd_weighted = b$Diff.Adj, ks_weighted = b$KS.Adj, pct_missing_ptz_carb = miss)
  }
  rbindlist(rows, fill = TRUE)
}

#' Population d'overlap (A10, §5.1) : caractéristiques de la cohorte totale
#' (non pondérée) vs de la population d'overlap (pondération OW, deux bras
#' réunis), moyennées sur les imputations
build_overlap_population_table_v31 <- function(imp_list, weights_list, ps_covs) {
  rows <- lapply(ps_covs, function(v) {
    d1 <- imp_list[[1]]
    if (is.character(d1[[v]]) && v != "gender") {
      lv <- sort(unique(d1[[v]]))
      return(rbindlist(lapply(lv, function(l) {
        u <- mean(vapply(imp_list, function(d) mean(d[[v]] == l), numeric(1)))
        o <- mean(vapply(seq_along(imp_list), function(i) stats::weighted.mean(imp_list[[i]][[v]] == l, weights_list[[i]]), numeric(1)))
        data.table(characteristic = paste0(lab_v31(v), ": ", l), stat = "%", full_cohort = 100 * u, overlap_population = 100 * o)
      })))
    }
    z <- function(d) if (v == "gender") as.numeric(toupper(d[[v]]) == "F") else as.numeric(d[[v]])
    bin <- v == "gender" || all(unique(z(d1)) %in% c(0, 1))
    u <- mean(vapply(imp_list, function(d) mean(z(d)), numeric(1)))
    o <- mean(vapply(seq_along(imp_list), function(i) stats::weighted.mean(z(imp_list[[i]]), weights_list[[i]]), numeric(1)))
    data.table(characteristic = lab_v31(v), stat = if (bin) "%" else "mean",
               full_cohort = if (bin) 100 * u else u, overlap_population = if (bin) 100 * o else o)
  })
  rbindlist(rows)
}

#' Flowchart v3.1 : flowchart v2 traduit + étapes v3/v3.1
build_flowchart_v31 <- function(flowchart_v2, aki_t0_tbl, master) {
  tr <- c("Épisodes PTZ/carb (après règle assignation)" = "PTZ/carbapenem episodes (after treatment-assignment rule)",
          "Exclus âge < 18" = "Excluded: age < 18", "Exclus créatinine baseline manquante" = "Excluded: missing baseline creatinine",
          "Exclus RRT avant T0" = "Excluded: RRT before T0", "Exclus ESRD ICD" = "Excluded: ESRD (ICD code)",
          "Exclus épisodes ultérieurs (même patient)" = "Excluded: later episodes (same patient)",
          "COHORTE FINALE" = "Cohort before AKI-at-T0 exclusion")
  fc <- copy(flowchart_v2)
  fc[, step := ifelse(step %in% names(tr), tr[step], step)]
  n_pre <- fc[step == "Cohort before AKI-at-T0 exclusion", n]
  add <- data.table(
    step = c("Excluded: AKI already present at T0 (KDIGO vs. 7-day pre-T0 minimum)",
             "Eligible cohort",
             "Excluded: no creatinine and no RRT in (T0; T0+7 d] (primary outcome not evaluable)",
             "PRIMARY ANALYTIC COHORT"),
    n = c(n_pre - n_of("eligible"), n_of("eligible"), n_of("eligible") - n_of("primary"), n_of("primary")),
    n_ptz = c(aki_t0_tbl[treatment == "PTZ", n_aki_at_t0], cohort_registry_table()[cohort_id == "eligible", n_ptz],
              NA, cohort_registry_table()[cohort_id == "primary", n_ptz]),
    n_carb = c(aki_t0_tbl[treatment == "carbapenem", n_aki_at_t0], cohort_registry_table()[cohort_id == "eligible", n_comparator],
               NA, cohort_registry_table()[cohort_id == "primary", n_comparator]))
  rbind(fc, add, fill = TRUE)
}

#' sessionInfo complet (A10 §12.5) : version de R, plateforme, et TOUS les
#' packages chargés (attachés + chargés via ::), avec versions
save_session_info_v31 <- function(dir, seed) {
  si <- utils::sessionInfo()
  pk <- c(si$otherPkgs, si$loadedOnly)
  pkg_tbl <- data.table(package = names(pk), version = vapply(pk, function(p) p$Version, character(1)),
                        attached = names(pk) %in% names(si$otherPkgs))[order(package)]
  saveRDS(list(session_info = si, r_version = si$R.version$version.string, platform = si$platform,
               packages = pkg_tbl, seed = seed, timestamp = Sys.time()), file.path(dir, "session_info_v31.rds"))
  writeLines(utils::capture.output(print(si)), file.path(dir, "session_info_v31.txt"))
  pkg_tbl
}

#' manuscript_numbers.json v3.1 : results_master + N du registre de cohortes
build_manuscript_numbers_json_v31 <- function(master, dir) {
  build_manuscript_numbers_json(master, dir = dir)
  p <- file.path(dir, "manuscript_numbers.json")
  nums <- jsonlite::read_json(p)
  reg <- cohort_registry_table()
  for (i in seq_len(nrow(reg))) {
    nums[[paste0("cohort_", reg$cohort_id[i], "_n")]] <- reg$n[i]
    nums[[paste0("cohort_", reg$cohort_id[i], "_n_ptz")]] <- reg$n_ptz[i]
    nums[[paste0("cohort_", reg$cohort_id[i], "_n_comparator")]] <- reg$n_comparator[i]
  }
  jsonlite::write_json(nums, p, auto_unbox = TRUE, pretty = TRUE, na = "null", digits = NA)
  invisible(nums)
}
