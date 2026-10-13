# =============================================================================
# 72_ed_vitals.R — Constantes des urgences (MIMIC-IV-ED 2.2) fusionnées
#   avec chartevents (docs/v3.2_plan.md §1.3, B3)
# =============================================================================
# En v3.1, FC / PAM / T° / SpO2 venaient uniquement de chartevents (réanimation) :
# 63-69 % de manquants, structurellement (patients de salle ou des urgences),
# puis imputés à partir des patients de réanimation. v3.2 ajoute :
#   - ed.triage (horodaté à edstays.intime), ed.vitalsign (charttime) ;
#   - fenêtre [T0 - 24 h ; T0], même hospitalisation (edstays.hadm_id) ;
#   - température °F -> °C ; PAM = (PAS + 2 x PAD) / 3 ;
#   - mêmes bornes physiologiques que la v3.1 (PHYSIO_RANGES_V31).
# Le GCS n'existe pas dans le module ED (reste chartevents seul).

ED_ROOT_V32 <- "../mimic-iv-ed-2.2/ed"

#' Mesures ED en format long (variable, valeur, heure) pour les hadm_id donnés
extract_ed_vitals_long_v32 <- function(hadm_ids, root = ED_ROOT_V32,
                                       cache = "output/analysis/derived_data/ed_vitals_long_v32.rds") {
  if (file.exists(cache)) return(readRDS(cache))
  st <- fread(file.path(root, "edstays.csv.gz"), select = c("subject_id", "hadm_id", "stay_id", "intime"))
  st <- st[hadm_id %in% hadm_ids][, intime := as.POSIXct(intime, tz = "UTC")]
  tr <- fread(file.path(root, "triage.csv.gz"), select = c("stay_id", "temperature", "heartrate", "o2sat", "sbp", "dbp"))
  tr <- merge(tr, st, by = "stay_id")[, charttime := intime][, source := "ed_triage"]
  vs <- fread(file.path(root, "vitalsign.csv.gz"), select = c("stay_id", "charttime", "temperature", "heartrate", "o2sat", "sbp", "dbp"))
  vs <- merge(vs, st, by = "stay_id")[, charttime := as.POSIXct(charttime, tz = "UTC")][, source := "ed_vitalsign"]
  x <- rbind(tr[, .(subject_id, hadm_id, charttime, source, temperature, heartrate, o2sat, sbp, dbp)],
             vs[, .(subject_id, hadm_id, charttime, source, temperature, heartrate, o2sat, sbp, dbp)])
  num <- function(v) suppressWarnings(as.numeric(v))
  x[, `:=`(temperature = num(temperature), heartrate = num(heartrate), o2sat = num(o2sat), sbp = num(sbp), dbp = num(dbp))]
  long <- rbind(
    x[!is.na(heartrate), .(subject_id, hadm_id, charttime, source, var = "hr", value = heartrate)],
    x[!is.na(sbp) & !is.na(dbp), .(subject_id, hadm_id, charttime, source, var = "map", value = (sbp + 2 * dbp) / 3)],
    # °F dans le module ED ; valeurs < 50 supposées déjà en °C (rares saisies)
    x[!is.na(temperature), .(subject_id, hadm_id, charttime, source, var = "temp",
                             value = fifelse(temperature > 50, (temperature - 32) * 5 / 9, temperature))],
    x[!is.na(o2sat), .(subject_id, hadm_id, charttime, source, var = "spo2", value = o2sat)])
  dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
  saveRDS(long, cache)
  message("[72_ed_v32] ED : ", nrow(long), " mesures pour ", uniqueN(long$hadm_id), " séjours")
  long
}

#' Constantes chartevents en format long (mêmes conversions que R/53_ v3.1)
.chartevents_vitals_long_v32 <- function(ce) {
  ids <- CHART_ITEMS_V31
  rbind(
    ce[itemid %in% ids$hr & !is.na(valuenum), .(subject_id, hadm_id, charttime, source = "chartevents", var = "hr", value = valuenum)],
    ce[itemid %in% ids$map & !is.na(valuenum), .(subject_id, hadm_id, charttime, source = "chartevents", var = "map", value = valuenum)],
    ce[itemid == ids$temp_f & !is.na(valuenum), .(subject_id, hadm_id, charttime, source = "chartevents", var = "temp", value = (valuenum - 32) * 5 / 9)],
    ce[itemid == ids$temp_c & !is.na(valuenum), .(subject_id, hadm_id, charttime, source = "chartevents", var = "temp", value = valuenum)],
    ce[itemid %in% ids$spo2 & !is.na(valuenum), .(subject_id, hadm_id, charttime, source = "chartevents", var = "spo2", value = valuenum)])
}

#' Résumé pré-T0 (24 h) des constantes, toutes sources confondues, après
#' bornes physiologiques. Renvoie aussi la source (ED seul / chartevents seul / les deux).
summarise_vitals_combined_v32 <- function(ce, ed_long, key, window_h = 24) {
  all <- rbind(.chartevents_vitals_long_v32(ce), ed_long[, .(subject_id, hadm_id, charttime, source, var, value)])
  all <- merge(all, key[, .(subject_id, hadm_id, T0)], by = c("subject_id", "hadm_id"))
  all <- all[charttime <= T0 & charttime >= T0 - window_h * 3600]
  keep <- rep(FALSE, nrow(all))
  for (v in unique(all$var)) { i <- all$var == v; keep[i] <- .in_range(all$value[i], v) }
  all <- all[keep]
  s <- all[, .(hr_mean = mean(value[var == "hr"]), map_mean = mean(value[var == "map"]),
               temp_max = suppressWarnings(max(value[var == "temp"])), spo2_min = suppressWarnings(min(value[var == "spo2"])),
               vitals_source = paste(sort(unique(ifelse(source == "chartevents", "ICU", "ED"))), collapse = "+")),
           by = .(subject_id, hadm_id)]
  for (v in c("hr_mean", "map_mean", "temp_max", "spo2_min")) s[!is.finite(get(v)), (v) := NA_real_]
  out <- merge(key[, .(subject_id, hadm_id)], s, by = c("subject_id", "hadm_id"), all.x = TRUE)
  out[is.na(vitals_source), vitals_source := "none"]
  out
}
