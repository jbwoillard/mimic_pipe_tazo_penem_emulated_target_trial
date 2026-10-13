# =============================================================================
# 31_baseline_creat.R — Créatinine de baseline et AKI prévalente à T0
#                            (docs/v3_plan.md §2.1 — point le plus critique)
# =============================================================================
# v2 : baseline = "closest to T0" (dernière valeur avant T0, 82% des cas).
# Chez des patients septiques, cette valeur est souvent déjà élevée, ce qui
# biaise le critère KDIGO ×1.5 vers le haut (régression vers la moyenne) —
# observé en v2 : gradient RD 7.2 -> 4.4 -> 2.0 pts par tertile de créatinine
# baseline (p=0.031 pour l'interaction).
#
# v3 (primaire) : baseline = créatinine MINIMALE sur [T0-7j ; T0].
# aki_at_t0 est désormais réellement calculée (en v2 : cohort[, aki_at_t0 := 0L]
# pour tous les patients, jamais évaluée) et les patients en AKI à T0 sont
# exclus de l'analyse principale (rapportés par bras avant exclusion).
#
# Dépend de : compute_egfr_ckd_epi() (R/04_extract_labs.R).

#' @param creat_dt  data.table (subject_id, hadm_id, charttime, creat_mgdl) —
#'                  output/extraction/derived_data/creatinine_all.rds (cache v2, lu seul)
#' @param T0_dt     data.table (subject_id, hadm_id, T0, admit_time, age, gender)
#' @return data.table avec toutes les définitions de baseline + aki_at_t0
compute_baseline_creatinine_v3 <- function(creat_dt, T0_dt) {

  base <- T0_dt[, .(subject_id, hadm_id, T0, admit_time,
                     age = if ("age" %in% names(T0_dt)) age else NA_real_,
                     gender = if ("gender" %in% names(T0_dt)) gender else NA_character_)]

  # -----------------------------------------------------------------
  # 1. PRIMAIRE : minimum de créatinine sur [T0-7j ; T0], même hospitalisation
  # -----------------------------------------------------------------
  creat_hadm <- merge(creat_dt, base[, .(subject_id, hadm_id, T0)],
                       by = c("subject_id", "hadm_id"), all = FALSE,
                       allow.cartesian = TRUE)
  creat_hadm[, time_before_t0_h := as.numeric(difftime(T0, charttime, units = "hours"))]

  win7d <- creat_hadm[time_before_t0_h >= 0 & time_before_t0_h <= 7 * 24]
  min7d <- win7d[, .(
    baseline_creat_min7d      = min(creat_mgdl),
    baseline_creat_min7d_time = charttime[which.min(creat_mgdl)],
    n_creat_in_7d_window      = .N
  ), by = .(subject_id, hadm_id)]

  # Dernière valeur avant T0 dans la fenêtre 7j (= "creatinine at T0", sert au
  # calcul de aki_at_t0 et correspond à la définition "closest to T0" v2)
  closest <- win7d[, .(
    creat_at_t0             = creat_mgdl[which.min(time_before_t0_h)],
    creat_at_t0_time        = charttime[which.min(time_before_t0_h)],
    baseline_creat_closest  = creat_mgdl[which.min(time_before_t0_h)]  # alias, pour sensibilité v2
  ), by = .(subject_id, hadm_id)]

  out <- merge(base[, .(subject_id, hadm_id, T0, admit_time, age, gender)],
               min7d, by = c("subject_id", "hadm_id"), all.x = TRUE)
  out <- merge(out, closest, by = c("subject_id", "hadm_id"), all.x = TRUE)

  # -----------------------------------------------------------------
  # 2. Référence 48h pour le critère absolu KDIGO (+0.3 mg/dL) : minimum des
  #    valeurs STRICTEMENT avant creat_at_t0 dans les 48h précédentes
  # -----------------------------------------------------------------
  win48h_prior <- creat_hadm[time_before_t0_h > 0 & time_before_t0_h <= 48]
  # exclure la valeur elle-même : on ne garde que les mesures AVANT le point
  # utilisé comme creat_at_t0 (qui a time_before_t0_h minimal)
  win48h_prior <- merge(win48h_prior, closest[, .(subject_id, hadm_id, creat_at_t0_time)],
                         by = c("subject_id", "hadm_id"), all.x = TRUE)
  win48h_prior <- win48h_prior[charttime < creat_at_t0_time]
  ref48 <- win48h_prior[, .(creat_48h_ref = min(creat_mgdl)), by = .(subject_id, hadm_id)]
  out <- merge(out, ref48, by = c("subject_id", "hadm_id"), all.x = TRUE)

  # -----------------------------------------------------------------
  # 3. AKI à T0 (KDIGO créatinine, calculée pour de vrai — v2 mettait 0
  #    partout sans jamais l'évaluer)
  # -----------------------------------------------------------------
  out[, aki_at_t0 := as.integer(
    (!is.na(creat_at_t0) & !is.na(baseline_creat_min7d) &
       creat_at_t0 >= 1.5 * baseline_creat_min7d) |
      (!is.na(creat_at_t0) & !is.na(creat_48h_ref) &
         (creat_at_t0 - creat_48h_ref) >= 0.3)
  )]
  out[is.na(aki_at_t0), aki_at_t0 := 0L]  # pas de mesure interprétable -> non évaluable, traité comme 0 mais signalé
  out[, aki_at_t0_evaluable := as.integer(!is.na(creat_at_t0) & !is.na(baseline_creat_min7d))]

  # -----------------------------------------------------------------
  # 4. Sensibilité : valeur ambulatoire/pré-admission la plus proche sur 365j
  #    (labevents avec ou sans hadm_id — ~40% des lignes créatinine ont
  #    hadm_id manquant = mesures hors de l'hospitalisation index, rattachées
  #    au patient uniquement par date)
  # -----------------------------------------------------------------------
  creat_subj <- merge(creat_dt[, .(subject_id, charttime, creat_mgdl)],
                       base[, .(subject_id, T0)],
                       by = "subject_id", all = FALSE, allow.cartesian = TRUE)
  creat_subj[, days_before_t0 := as.numeric(difftime(T0, charttime, units = "days"))]
  # Fenêtre stricte "ambulatoire/pré-admission" : > 7j avant T0 (pour ne pas
  # doublonner la fenêtre primaire) et <= 365j avant T0
  amb <- creat_subj[days_before_t0 > 7 & days_before_t0 <= 365]
  amb_best <- amb[, .(
    baseline_creat_ambulatory_365d = creat_mgdl[which.min(days_before_t0)],
    baseline_creat_ambulatory_365d_days_before = round(min(days_before_t0), 1)
  ), by = subject_id]
  out <- merge(out, amb_best, by = "subject_id", all.x = TRUE)

  # -----------------------------------------------------------------
  # 5. Sensibilité : imputation par eGFR supposée = 75 mL/min/1.73m² (CKD-EPI
  #    2021 inversée) quand aucune créatinine n'est disponible sur 7j
  # -----------------------------------------------------------------------
  invert_ckd_epi_for_egfr75 <- function(age, sex, egfr_target = 75) {
    kappa <- ifelse(toupper(sex) == "F", 0.7, 0.9)
    alpha <- ifelse(toupper(sex) == "F", -0.241, -0.302)
    sex_m <- ifelse(toupper(sex) == "F", 1.012, 1.0)
    # Résolution numérique simple (la formule CKD-EPI 2021 n'est pas
    # inversible en forme close à cause du min/max piecewise) : recherche par
    # dichotomie sur creat in [0.3, 5]
    solve_one <- function(a, s) {
      f <- function(cr) compute_egfr_ckd_epi(cr, a, s) - egfr_target
      tryCatch(stats::uniroot(f, interval = c(0.2, 8))$root, error = function(e) NA_real_)
    }
    mapply(solve_one, age, sex)
  }
  out[, baseline_creat_mdrd_impute := invert_ckd_epi_for_egfr75(
    age, ifelse(toupper(gender) == "F", "F", "M"))]

  message("[31_baseline_v3] Baseline min-7j disponible pour ",
          sum(!is.na(out$baseline_creat_min7d)), "/", nrow(out), " épisodes")
  message("[31_baseline_v3] AKI à T0 (KDIGO, calculée) : ",
          sum(out$aki_at_t0 == 1), "/", sum(out$aki_at_t0_evaluable == 1),
          " évaluables (", nrow(out) - sum(out$aki_at_t0_evaluable), " non évaluables -> aki_at_t0 mis à 0 par défaut, signalé séparément)")

  out[, .(subject_id, hadm_id, baseline_creat_min7d, baseline_creat_min7d_time,
          n_creat_in_7d_window, creat_at_t0, creat_at_t0_time, creat_48h_ref,
          aki_at_t0, aki_at_t0_evaluable,
          baseline_creat_closest,
          baseline_creat_ambulatory_365d, baseline_creat_ambulatory_365d_days_before,
          baseline_creat_mdrd_impute)]
}

#' Rapport de distribution de aki_at_t0 par bras (AVANT exclusion), requis
#' explicitement par le cahier des charges §2.1
report_aki_at_t0_by_arm <- function(cohort_v3) {
  tbl <- cohort_v3[, .(
    n = .N,
    n_evaluable   = sum(aki_at_t0_evaluable, na.rm = TRUE),
    n_aki_at_t0   = sum(aki_at_t0, na.rm = TRUE),
    pct_aki_at_t0 = round(mean(aki_at_t0, na.rm = TRUE) * 100, 1),
    pct_aki_at_t0_among_evaluable = round(
      sum(aki_at_t0[aki_at_t0_evaluable == 1]) / sum(aki_at_t0_evaluable) * 100, 1)
  ), by = treatment]
  data.table::fwrite(tbl, v3_tables_path("aki_at_t0_by_arm.csv"))
  message("[31_baseline_v3] aki_at_t0_by_arm.csv sauvegardé")
  print(tbl)
  tbl
}
