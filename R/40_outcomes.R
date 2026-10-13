# =============================================================================
# 40_outcomes.R — Batterie d'outcomes v3 (§3.1) + correction de cohérence
#   avec la baseline créatinine min-7j (§2.1)
# =============================================================================
# IMPORTANT : la Phase B a changé la baseline créatinine utilisée comme
# CONFONDANT et comme critère d'exclusion (aki_at_t0), mais outcome_primary
# (AKI à J7, hérité du cache v2) était encore calculé contre l'ANCIENNE
# baseline ("closest to T0", R/07_define_outcomes.R v2). Incohérence
# méthodologique : le critère KDIGO x1.5 de l'OUTCOME doit être évalué
# contre la MÊME baseline (min-7j) que celle utilisée pour définir l'AKI à
# T0. Corrigé ici : tous les horizons (D2/D3/D7/D14/D30), les stades, la
# persistance (48h ET 72h, la 72h est nouvelle) et la récupération à J14
# sont recalculés contre baseline_creat_min7d.

#' Recalcule toute la batterie d'outcomes créatinine contre baseline_creat_min7d
#' @param cohort   cohort_phaseB (doit contenir baseline_creat_min7d, T0)
#' @param creat_dt data.table créatinine (cache v2, lecture seule)
rebuild_aki_outcomes_v3 <- function(cohort, creat_dt) {
  abs_thr <- 0.3; rel_thr <- 1.5

  creat_post <- merge(creat_dt, cohort[, .(subject_id, hadm_id, T0, baseline_creat_min7d)],
                       by = c("subject_id", "hadm_id"), all = FALSE)
  creat_post <- creat_post[!is.na(T0) & !is.na(baseline_creat_min7d)]
  creat_post[, time_h := as.numeric(difftime(charttime, T0, units = "hours"))]
  creat_post[, time_d := time_h / 24]
  creat_post <- creat_post[time_h > 0]

  setorder(creat_post, subject_id, hadm_id, time_h)

  creat_post[, creat_48h_ref := {
    res <- numeric(.N)
    for (i in seq_len(.N)) {
      prev_vals <- creat_mgdl[time_h >= time_h[i] - 48 & time_h < time_h[i]]
      res[i] <- if (length(prev_vals) > 0) min(prev_vals) else baseline_creat_min7d[i]
    }
    res
  }, by = .(subject_id, hadm_id)]

  creat_post[, aki_absolute := as.integer(creat_mgdl - creat_48h_ref >= abs_thr)]
  creat_post[, aki_relative := as.integer(creat_mgdl >= rel_thr * baseline_creat_min7d)]
  creat_post[, aki_any      := as.integer(aki_absolute == 1 | aki_relative == 1)]

  creat_post[, aki_stage := dplyr::case_when(
    creat_mgdl >= 3.0 * baseline_creat_min7d                      ~ 3L,
    creat_mgdl >= 4.0 & (creat_mgdl - creat_48h_ref) >= 0.3        ~ 3L,
    creat_mgdl >= 2.0 * baseline_creat_min7d                      ~ 2L,
    creat_mgdl >= 1.5 * baseline_creat_min7d                      ~ 1L,
    (creat_mgdl - creat_48h_ref) >= 0.3                            ~ 1L,
    TRUE                                                           ~ 0L
  )]
  creat_post[, aki_stage2 := as.integer(aki_stage >= 2)]
  creat_post[, aki_stage3 := as.integer(aki_stage >= 3)]

  result <- cohort[, .(subject_id, hadm_id)]
  for (day in c(2, 3, 7, 14, 30)) {
    col_aki <- paste0("aki_d", day, "_v3")
    col_s2  <- paste0("aki_stage2_d", day, "_v3")
    col_s3  <- paste0("aki_stage3_d", day, "_v3")
    col_mcr <- paste0("max_creat_d", day, "_v3")
    subset_d <- creat_post[time_d > 0 & time_d <= day]
    if (nrow(subset_d) == 0) {
      result[, (col_aki) := NA_integer_]; result[, (col_s2) := NA_integer_]
      result[, (col_s3) := NA_integer_];  result[, (col_mcr) := NA_real_]
      next
    }
    agg <- subset_d[, .(
      aki_val = as.integer(any(aki_any == 1, na.rm = TRUE)),
      s2_val  = as.integer(any(aki_stage2 == 1, na.rm = TRUE)),
      s3_val  = as.integer(any(aki_stage3 == 1, na.rm = TRUE)),
      mcr     = max(creat_mgdl, na.rm = TRUE)
    ), by = .(subject_id, hadm_id)]
    setnames(agg, c("aki_val", "s2_val", "s3_val", "mcr"), c(col_aki, col_s2, col_s3, col_mcr))
    result <- merge(result, agg, by = c("subject_id", "hadm_id"), all.x = TRUE)
  }

  # Persistance >=48h ET >=72h (nouveau)
  aki_persist <- creat_post[aki_any == 1,
                             .(n_aki_meas = .N, first_aki_h = min(time_h), last_aki_h = max(time_h)),
                             by = .(subject_id, hadm_id)]
  aki_persist[, aki_persistent_48h_v3 := as.integer(n_aki_meas >= 2 & (last_aki_h - first_aki_h) >= 48)]
  aki_persist[, aki_persistent_72h_v3 := as.integer(n_aki_meas >= 2 & (last_aki_h - first_aki_h) >= 72)]
  result <- merge(result, aki_persist[, .(subject_id, hadm_id, aki_persistent_48h_v3, aki_persistent_72h_v3)],
                   by = c("subject_id", "hadm_id"), all.x = TRUE)
  result[is.na(aki_persistent_48h_v3), aki_persistent_48h_v3 := 0L]
  result[is.na(aki_persistent_72h_v3), aki_persistent_72h_v3 := 0L]

  result[, outcome_primary_v3 := aki_d7_v3]

  n_changed <- sum(result$outcome_primary_v3 != cohort$outcome_primary, na.rm = TRUE)
  message("[40_outcomes_v3] Recalcul AKI (baseline min-7j) : ",
          sum(!is.na(result$outcome_primary_v3)), "/", nrow(result),
          " évaluables. Reclassifiés vs v2/Phase-B outcome_primary (ancienne baseline) : ",
          n_changed, " (", round(n_changed / sum(!is.na(result$outcome_primary_v3)) * 100, 1), "%)")

  result
}

#' RRT sur une fenêtre paramétrable (v2 n'avait que rrt_7d)
extract_rrt_window_v3 <- function(tbl_paths, cohort, window_days = 30) {
  paths <- tbl_paths$paths
  if (is.na(paths[["procedureevents"]])) {
    message("[40_outcomes_v3] procedureevents indisponible — rrt_", window_days, "d ignoré")
    return(cohort[, .(subject_id, hadm_id, rrt = 0L)][, setnames(.SD, "rrt", paste0("rrt_", window_days, "d"))])
  }
  pe <- read_mimic_table(paths[["procedureevents"]],
                          col_select = c("subject_id", "hadm_id", "itemid", "starttime"))
  rrt_item_ids <- c(225441L, 225802L, 225803L, 225805L, 225809L, 225955L, 226118L, 227290L)
  rrt_pe <- pe[itemid %in% rrt_item_ids]
  rrt_pe[, starttime := lubridate::ymd_hms(starttime, quiet = TRUE)]
  rrt_pe[, subject_id := as.integer(subject_id)]
  rrt_pe[, hadm_id := as.integer(hadm_id)]

  rrt_t0 <- merge(rrt_pe, cohort[, .(subject_id, hadm_id, T0)],
                   by = c("subject_id", "hadm_id"), all = FALSE, allow.cartesian = TRUE)
  col_name <- paste0("rrt_", window_days, "d")
  rrt_win <- rrt_t0[!is.na(starttime) & !is.na(T0) & starttime >= T0 &
                       starttime <= T0 + window_days * 24 * 3600,
                     .(flag = 1L), by = .(subject_id, hadm_id)]
  setnames(rrt_win, "flag", col_name)
  result <- merge(cohort[, .(subject_id, hadm_id)], rrt_win, by = c("subject_id", "hadm_id"), all.x = TRUE)
  result[is.na(get(col_name)), (col_name) := 0L]
  message("[40_outcomes_v3] ", col_name, " : ", sum(result[[col_name]]), " événements")
  result
}

#' Décès sur une fenêtre paramétrable (deathtime précis + dod en secours)
compute_death_window_v3 <- function(cohort, window_days) {
  col_name <- paste0("death_", window_days, "d_v3")
  dt <- cohort[, .(subject_id, hadm_id, T0, deathtime, dod)]
  dt[, (col_name) := as.integer(
    (!is.na(deathtime) & deathtime <= T0 + as.difftime(window_days, units = "days")) |
      (!is.na(dod) & as.Date(dod) <= as.Date(T0) + window_days)
  )]
  dt[is.na(get(col_name)), (col_name) := 0L]
  dt[, .SD, .SDcols = c("subject_id", "hadm_id", col_name)]
}

#' Récupération rénale à J14 : parmi les patients en AKI à J7 (v3), créatinine
#' revenue <= 1.2x baseline_min7d autour de J14 (fenêtre J12-J16)
compute_renal_recovery_d14_v3 <- function(cohort_with_outcomes, creat_dt) {
  aki_d7_cases <- cohort_with_outcomes[aki_d7_v3 == 1, .(subject_id, hadm_id, T0, baseline_creat_min7d)]
  if (nrow(aki_d7_cases) == 0) {
    return(cohort_with_outcomes[, .(subject_id, hadm_id, renal_recovery_d14_v3 = NA_integer_)])
  }
  cr14 <- merge(creat_dt, aki_d7_cases, by = c("subject_id", "hadm_id"), all = FALSE)
  cr14[, time_d := as.numeric(difftime(charttime, T0, units = "days"))]
  cr14_win <- cr14[time_d >= 12 & time_d <= 16]
  recov <- cr14_win[, .(closest_to_d14 = creat_mgdl[which.min(abs(time_d - 14))],
                         baseline_creat_min7d = baseline_creat_min7d[1]),
                     by = .(subject_id, hadm_id)]
  recov[, renal_recovery_d14_v3 := as.integer(closest_to_d14 <= 1.2 * baseline_creat_min7d)]

  result <- merge(cohort_with_outcomes[, .(subject_id, hadm_id, aki_d7_v3)],
                   recov[, .(subject_id, hadm_id, renal_recovery_d14_v3)],
                   by = c("subject_id", "hadm_id"), all.x = TRUE)
  n_eval <- sum(!is.na(result$renal_recovery_d14_v3) & result$aki_d7_v3 == 1)
  message("[40_outcomes_v3] Récupération rénale J14 évaluable pour ", n_eval,
          "/", sum(result$aki_d7_v3 == 1, na.rm = TRUE), " cas d'AKI-J7 (fenêtre J12-J16, ",
          "manquant = pas de créatinine dans cette fenêtre, ex. sortie précoce)")
  result[, .(subject_id, hadm_id, renal_recovery_d14_v3)]
}

#' MAKE-30 v3 (mêmes composantes que v2, mais death_30d/rrt_30d/renal_dysfunction
#' peuvent différer légèrement selon la baseline utilisée pour renal_dysfunction)
compute_make30_v3 <- function(cohort, creat_dt, rrt_30d_v3) {
  creat_30 <- merge(creat_dt, cohort[, .(subject_id, hadm_id, T0, baseline_creat_min7d)],
                     by = c("subject_id", "hadm_id"), all = FALSE)
  creat_30[, time_d := as.numeric(difftime(charttime, T0, units = "days"))]
  creat_30 <- creat_30[time_d > 0 & time_d <= 30]
  cr_persist <- creat_30[creat_mgdl >= 2 * baseline_creat_min7d,
                          .(renal_dysfunction_30d_v3 = 1L), by = .(subject_id, hadm_id)]
  out <- merge(cohort[, .(subject_id, hadm_id, death_30d)], cr_persist,
               by = c("subject_id", "hadm_id"), all.x = TRUE)
  out <- merge(out, rrt_30d_v3, by = c("subject_id", "hadm_id"), all.x = TRUE)
  out[is.na(renal_dysfunction_30d_v3), renal_dysfunction_30d_v3 := 0L]
  out[, make30_v3 := as.integer(death_30d == 1 | rrt_30d == 1 | renal_dysfunction_30d_v3 == 1)]
  message("[40_outcomes_v3] MAKE-30 (v3) incidence : ", round(mean(out$make30_v3, na.rm = TRUE) * 100, 1), "%")
  out[, .(subject_id, hadm_id, make30_v3)]
}
