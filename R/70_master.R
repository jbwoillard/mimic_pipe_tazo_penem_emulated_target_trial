# =============================================================================
# 70_master.R — Cohorte maîtresse v3.2 (docs/v3.2_plan.md §1.1, D1)
# =============================================================================
# UNE seule cohorte pour toutes les populations : premières ADMINISTRATIONS
# documentées (eMAR, événements d'administration, voie IV) de PTZ, céfépime
# ou méropénème (R/57_), un épisode par patient. La cohorte v3.1 (héritée de
# la v2) mélangeait initiations et escalades : 42 % des « initiateurs de
# carbapénème » avaient reçu du céfépime avant T0 dans le séjour.
# Ajoute aux covariables / outcomes de R/57_ :
#   - trajectoire de créatinine pré-T0 (B2) ;
#   - constantes ED fusionnées avec chartevents (B3) et lieu de soins à T0 ;
#   - outcomes manquants (persistance, MAKE30, max quotidien, temps
#     d'événements pour l'incidence cumulée).

#' Trajectoire de créatinine pré-T0 (B2). Toutes les valeurs <= T0.
creat_trajectory_v32 <- function(creat_dt, key) {
  cr <- merge(creat_dt[, .(subject_id, hadm_id, charttime, creat_mgdl)], key[, .(subject_id, hadm_id, T0)],
              by = c("subject_id", "hadm_id"))
  cr[, h := as.numeric(difftime(charttime, T0, units = "hours"))]
  w7 <- cr[h <= 0 & h >= -7 * 24][order(subject_id, hadm_id, h)]
  tr <- w7[, .(creat_last_pre_t0 = creat_mgdl[.N], min7 = min(creat_mgdl), n7 = .N), by = .(subject_id, hadm_id)]
  w48 <- cr[h <= 0 & h >= -48]
  sl <- w48[, .(creat_slope_48h = if (.N >= 2 && stats::var(h) > 0) unname(stats::coef(stats::lm(creat_mgdl ~ h))[2]) * 24 else NA_real_),
            by = .(subject_id, hadm_id)]
  out <- merge(key[, .(subject_id, hadm_id)], tr, by = c("subject_id", "hadm_id"), all.x = TRUE)
  out <- merge(out, sl, by = c("subject_id", "hadm_id"), all.x = TRUE)
  out[, `:=`(creat_delta_last_min7d = creat_last_pre_t0 - min7, creat_ratio_last_min7d = creat_last_pre_t0 / min7,
             creat_slope_48h_available = as.integer(!is.na(creat_slope_48h)),
             creat_single_pre_value = as.integer(is.na(n7) | n7 <= 1L))]
  # pente non calculable (< 2 valeurs en 48 h) : 0 + indicateur (pas d'imputation)
  out[is.na(creat_slope_48h), creat_slope_48h := 0]
  out[, c("min7", "n7") := NULL]
  out
}

#' Temps d'événements (heures depuis T0) pour l'incidence cumulée à J14 (A4)
event_times_v32 <- function(m, tbl_paths) {
  pe <- read_mimic_table(tbl_paths$paths[["procedureevents"]], col_select = c("subject_id", "hadm_id", "itemid", "starttime"))
  pe <- pe[itemid %in% c(225441L, 225802L, 225803L, 225805L, 225809L, 225955L, 226118L, 227290L)]
  pe[, `:=`(subject_id = as.integer(subject_id), hadm_id = as.integer(hadm_id), starttime = lubridate::ymd_hms(starttime, quiet = TRUE))]
  r <- merge(pe, m[, .(subject_id, hadm_id, T0)], by = c("subject_id", "hadm_id"))[starttime > T0]
  r <- r[, .(rrt_first_h = as.numeric(difftime(min(starttime), T0[1], units = "hours"))), by = .(subject_id, hadm_id)]
  x <- merge(m[, .(subject_id, hadm_id, T0, deathtime, dod, dischtime, t_first_aki_h_d14_v31)], r,
             by = c("subject_id", "hadm_id"), all.x = TRUE)
  x[, death_h := fifelse(!is.na(deathtime), as.numeric(difftime(deathtime, T0, units = "hours")),
                         fifelse(!is.na(dod), as.numeric(difftime(as.POSIXct(dod, tz = "UTC") + 86399, T0, units = "hours")), NA_real_))]
  x[, disch_h := as.numeric(difftime(dischtime, T0, units = "hours"))]
  x[, aki_h := pmin(t_first_aki_h_d14_v31, rrt_first_h, na.rm = TRUE)]
  x[is.infinite(aki_h), aki_h := NA_real_]
  x[, .(subject_id, hadm_id, aki_h, rrt_first_h, death_h, disch_h)]
}

build_master_v32 <- function(ap, ce, ed_long, creat_dt, tbl_paths, out_tables) {
  m <- copy(ap)
  key <- m[, .(subject_id, hadm_id, T0, treatment)]
  # --- B2 : trajectoire de créatinine
  m <- merge(m, creat_trajectory_v32(creat_dt, key), by = c("subject_id", "hadm_id"), all.x = TRUE)
  # --- B3 : constantes ED + chartevents ; versions chartevents seules conservées (waterfall)
  setnames(m, c("hr_mean", "map_mean", "temp_max", "spo2_min"), c("hr_mean_icu", "map_mean_icu", "temp_max_icu", "spo2_min_icu"))
  m <- merge(m, summarise_vitals_combined_v32(ce, ed_long, key), by = c("subject_id", "hadm_id"), all.x = TRUE)
  m[, care_location_t0 := fcase(grepl("^ICU", service_group), "ICU", service_group == "ED_observation", "ED", default = "ward")]
  # --- outcomes complémentaires
  m[, egfr_baseline_v31 := compute_egfr_ckd_epi(baseline_creat_min7d, age, ifelse(toupper(gender) == "F", "F", "M"))]
  for (h in c(48, 72)) {
    p <- compute_aki_persistence_v31(creat_dt, m, "aki_d7_v31", "t_first_aki_h_d7_v31", h)
    setnames(p, "out", paste0("aki_persistent_", h, "h_v31")); m <- merge(m, p, by = c("subject_id", "hadm_id"), all.x = TRUE)
  }
  cr30 <- merge(creat_dt, m[, .(subject_id, hadm_id, T0, baseline_creat_min7d)], by = c("subject_id", "hadm_id"))
  rd30 <- cr30[charttime > T0 & charttime <= T0 + 30 * 86400 & creat_mgdl >= 2 * baseline_creat_min7d,
               .(renal_dysfunction_30d = 1L), by = .(subject_id, hadm_id)]
  m <- merge(m, rd30, by = c("subject_id", "hadm_id"), all.x = TRUE)
  m[is.na(renal_dysfunction_30d), renal_dysfunction_30d := 0L]
  m[, make30_v3 := as.integer(death_30d == 1L | rrt_30d == 1L | renal_dysfunction_30d == 1L)]
  dm <- copy(cr30)[, day := floor(as.numeric(difftime(charttime, T0, units = "hours")) / 24)][order(-creat_mgdl)][, .SD[1], by = .(subject_id, hadm_id, day)]
  kdm <- compute_kdigo_v31(dm[, .(subject_id, hadm_id, charttime, creat_mgdl)], m, 7, "rrt_7d", suffix = "_d7_dailymax_v31")
  m <- merge(m, kdm[, .(subject_id, hadm_id, aki_d7_dailymax_v31)], by = c("subject_id", "hadm_id"), all.x = TRUE)
  nd <- cr30[charttime > T0 & charttime <= T0 + 7 * 86400, .(n_creat_draws_7d = .N), by = .(subject_id, hadm_id)]
  m <- merge(m, nd, by = c("subject_id", "hadm_id"), all.x = TRUE); m[is.na(n_creat_draws_7d), n_creat_draws_7d := 0L]
  # contrôle négatif (a) variante 0-6 h (B7)
  e6 <- compute_early_creat_rise_v31(creat_dt, m, window_h = 6); setnames(e6, "early_creat_rise_12h", "early_creat_rise_6h")
  m <- merge(m, e6, by = c("subject_id", "hadm_id"), all.x = TRUE)
  # temps d'événements (A4)
  m <- merge(m, event_times_v32(m, tbl_paths), by = c("subject_id", "hadm_id"), all.x = TRUE)
  # année d'admission approximative (validation par période HTE)
  m[, approx_admit_year := approx_admit_year_v31(m)]
  # --- manquants avant / après ajout des constantes ED, par lieu (rapport B3)
  miss <- rbindlist(lapply(c("hr_mean", "map_mean", "temp_max", "spo2_min"), function(v)
    m[, .(variable = v, pct_missing_chartevents_only = round(100 * mean(is.na(get(paste0(v, "_icu")))), 1),
          pct_missing_with_ed = round(100 * mean(is.na(get(v))), 1), n = .N), by = care_location_t0]))
  fwrite(miss, file.path(out_tables, "vitals_missingness_ed_v32.csv"))
  m[, outcome_primary := outcome_primary_v31]
  setkey(m, hadm_id)
  message("[70_master_v32] master v3.2 : n=", nrow(m), " (", paste(names(table(m$treatment)), table(m$treatment), sep = "=", collapse = ", "), ")")
  m
}
