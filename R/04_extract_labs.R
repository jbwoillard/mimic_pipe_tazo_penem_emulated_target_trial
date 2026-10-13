# =============================================================================
# 04_extract_labs.R — Extraction créatinine et autres biomarqueurs depuis MIMIC-IV
# =============================================================================

#' Extraire les valeurs de créatinine depuis labevents
#' @return data.table (subject_id, hadm_id, charttime, creat_mgdl)
extract_creatinine <- function(tbl_paths, cfg) {

  paths <- tbl_paths$paths

  if (is.na(paths[["labevents"]]) || is.na(paths[["d_labitems"]])) {
    stop("[04_labs] labevents ou d_labitems indisponible — créatinine impossible à extraire")
  }

  # Identifier les item_id créatinine
  d_lab <- read_mimic_table(paths[["d_labitems"]],
                              col_select = c("itemid", "label", "fluid",
                                             "category", "valueuom"))

  if (is.null(d_lab)) stop("[04_labs] d_labitems vide")

  creat_pat <- paste(cfg$labs$creatinine_patterns, collapse = "|")
  creat_ids <- d_lab[grepl(creat_pat, label, ignore.case = TRUE) &
                       grepl("blood|serum", fluid, ignore.case = TRUE) &
                       !grepl("urine|csf|pleural|dialysate", fluid,
                              ignore.case = TRUE), itemid]

  if (length(creat_ids) == 0) {
    # Fallback : créatinine sanguine sans filtre fluide
    creat_ids <- d_lab[grepl("^creatinine$", label, ignore.case = TRUE), itemid]
  }

  if (length(creat_ids) == 0) {
    # Dernier recours
    creat_ids <- d_lab[grepl(creat_pat, label, ignore.case = TRUE), itemid]
  }

  message("[04_labs] item_id créatinine détectés : ",
          paste(creat_ids, collapse = ", "),
          " (labels : ", paste(d_lab[itemid %in% creat_ids, label],
                                collapse = ", "), ")")

  # Lire labevents (grand fichier — utiliser duckdb si disponible)
  lab <- read_mimic_table(paths[["labevents"]],
                           col_select = c("subject_id", "hadm_id", "itemid",
                                          "charttime", "valuenum", "valueuom",
                                          "flag"),
                           use_duckdb = TRUE)

  if (is.null(lab)) stop("[04_labs] labevents vide")

  creat <- lab[itemid %in% creat_ids & !is.na(valuenum) & valuenum > 0]

  # Conversion µmol/L → mg/dL si nécessaire
  # Seuil heuristique : si médiane > 20, probablement µmol/L
  median_val <- median(creat$valuenum, na.rm = TRUE)
  if (median_val > 20) {
    message("[04_labs] Conversion µmol/L -> mg/dL détectée (médiane = ",
            round(median_val, 1), ")")
    creat[, valuenum := valuenum / 88.42]
  }

  # Filtrer valeurs physiologiquement impossibles
  creat <- creat[valuenum >= 0.1 & valuenum <= 30]

  creat[, charttime := lubridate::ymd_hms(charttime, quiet = TRUE)]
  creat[, creat_mgdl := valuenum]
  creat[, subject_id := as.integer(subject_id)]
  creat[, hadm_id    := as.integer(hadm_id)]

  setorder(creat, subject_id, charttime)

  message("[04_labs] Mesures créatinine valides : ", nrow(creat),
          " pour ", uniqueN(creat$subject_id), " patients")

  save_derived(creat[, .(subject_id, hadm_id, charttime, creat_mgdl)],
               "creatinine_all", format = "rds")

  creat[, .(subject_id, hadm_id, charttime, creat_mgdl)]
}

#' Calculer créatinine baseline par patient-épisode
#' Priorité 1 : dernière valeur avant T0 (admission → T0)
#' Priorité 2 : valeur la plus proche dans les 24h avant T0
#' Priorité 3 : minimum dans les 7j avant T0
#' @param creat_dt  data.table créatinine (subject_id, hadm_id, charttime, creat_mgdl)
#' @param T0_dt     data.table avec subject_id, hadm_id, T0, admit_time
#' @return data.table avec baseline_creat, baseline_creat_strategy, baseline_creat_time
compute_baseline_creatinine <- function(creat_dt, T0_dt) {

  result <- T0_dt[, .(subject_id, hadm_id, T0,
                       admit_time = if ("admit_time" %in% names(T0_dt))
                         admit_time else T0 - 24 * 3600)]

  # Merge creatinine avec T0/admit_time (merge explicite, plus safe que data.table join i.)
  creat_with_t0 <- merge(creat_dt, result, by = c("subject_id", "hadm_id"),
                          all = FALSE, allow.cartesian = TRUE)

  creat_with_t0[, time_before_t0_h := as.numeric(
    difftime(T0, charttime, units = "hours")
  )]

  # Stratégie 1 : après admission et avant T0
  strat1 <- creat_with_t0[
    time_before_t0_h >= 0 & charttime >= admit_time,
    .(baseline_creat      = creat_mgdl[which.min(time_before_t0_h)],
      baseline_creat_time = charttime[which.min(time_before_t0_h)],
      baseline_creat_strategy = "last_before_T0"),
    by = .(subject_id, hadm_id)
  ]

  # Stratégie 2 : dans les 24h avant T0 (même hors admission actuelle)
  strat2 <- creat_with_t0[
    time_before_t0_h >= 0 & time_before_t0_h <= 24,
    .(baseline_creat      = creat_mgdl[which.min(time_before_t0_h)],
      baseline_creat_time = charttime[which.min(time_before_t0_h)],
      baseline_creat_strategy = "last_24h_before_T0"),
    by = .(subject_id, hadm_id)
  ]

  # Stratégie 3 : minimum dans les 7j avant T0
  strat3 <- creat_with_t0[
    time_before_t0_h >= 0 & time_before_t0_h <= 7 * 24,
    .(baseline_creat      = min(creat_mgdl),
      baseline_creat_time = charttime[which.min(creat_mgdl)],
      baseline_creat_strategy = "min_7d_before_T0"),
    by = .(subject_id, hadm_id)
  ]

  # Priorité 1 → 2 → 3
  baseline <- copy(result[, .(subject_id, hadm_id)])
  baseline <- merge(baseline, strat1, by = c("subject_id", "hadm_id"), all.x = TRUE)

  # Remplir avec strat2 si manquant
  missing_idx <- is.na(baseline$baseline_creat)
  if (any(missing_idx)) {
    fill2 <- strat2[baseline[missing_idx, .(subject_id, hadm_id)],
                     on = .(subject_id, hadm_id)]
    baseline[missing_idx, c("baseline_creat", "baseline_creat_time",
                              "baseline_creat_strategy") :=
               fill2[, .(baseline_creat, baseline_creat_time,
                          baseline_creat_strategy)]]
  }

  # Remplir avec strat3 si encore manquant
  missing_idx2 <- is.na(baseline$baseline_creat)
  if (any(missing_idx2)) {
    fill3 <- strat3[baseline[missing_idx2, .(subject_id, hadm_id)],
                     on = .(subject_id, hadm_id)]
    baseline[missing_idx2, c("baseline_creat", "baseline_creat_time",
                               "baseline_creat_strategy") :=
               fill3[, .(baseline_creat, baseline_creat_time,
                          baseline_creat_strategy)]]
  }

  message("[04_labs] Créatinine baseline disponible pour ",
          sum(!is.na(baseline$baseline_creat)), "/", nrow(baseline),
          " épisodes")
  message("[04_labs] Stratégies utilisées :")
  print(baseline[, .N, by = baseline_creat_strategy])

  baseline
}

#' Calculer eGFR CKD-EPI 2021 (sans race)
#' @param creat_mgdl vecteur créatinine mg/dL
#' @param age vecteur âge
#' @param sex vecteur sexe ("M" ou "F")
compute_egfr_ckd_epi <- function(creat_mgdl, age, sex) {
  # CKD-EPI 2021 race-free
  kappa  <- ifelse(toupper(sex) == "F", 0.7, 0.9)
  alpha  <- ifelse(toupper(sex) == "F", -0.241, -0.302)
  sex_m  <- ifelse(toupper(sex) == "F", 1.012, 1.0)

  scr_k  <- creat_mgdl / kappa
  term1  <- pmin(scr_k, 1)^alpha
  term2  <- pmax(scr_k, 1)^(-1.200)

  142 * term1 * term2 * (0.9938^age) * sex_m
}

#' Extraire BUN (urée sanguine)
extract_bun <- function(tbl_paths, cfg) {

  paths <- tbl_paths$paths
  if (is.na(paths[["labevents"]]) || is.na(paths[["d_labitems"]])) return(NULL)

  d_lab <- read_mimic_table(paths[["d_labitems"]],
                              col_select = c("itemid", "label", "fluid"))

  bun_pat <- paste(cfg$labs$bun_patterns, collapse = "|")
  bun_ids <- d_lab[grepl(bun_pat, label, ignore.case = TRUE) &
                     !grepl("urine", fluid, ignore.case = TRUE), itemid]

  if (length(bun_ids) == 0) return(NULL)

  lab <- read_mimic_table(paths[["labevents"]],
                           col_select = c("subject_id", "hadm_id", "itemid",
                                          "charttime", "valuenum"),
                           use_duckdb = TRUE)

  bun_dt <- lab[itemid %in% bun_ids & !is.na(valuenum) & valuenum > 0 &
                   valuenum < 300]
  bun_dt[, charttime := lubridate::ymd_hms(charttime, quiet = TRUE)]
  bun_dt[, subject_id := as.integer(subject_id)]
  bun_dt[, hadm_id    := as.integer(hadm_id)]

  message("[04_labs] Mesures BUN valides : ", nrow(bun_dt))
  bun_dt[, .(subject_id, hadm_id, charttime, bun_mgdl = valuenum)]
}
