# =============================================================================
# 32_covariates.R — Corrections de fuite d'information + nouvelles
#   covariables (docs/v3_plan.md §2.2, §2.3, §1.13, §2.5)
# =============================================================================
# Dépend de : R/02_codebooks.R (get_icd_codes, get_charlson_codes,
#             CHARLSON_WEIGHTS), R/01_paths_and_io.R (read_mimic_table).

# -----------------------------------------------------------------------
# §1.13 — calendar_year (bruit, dates décalées 2110-2212) -> anchor_year_group
# -----------------------------------------------------------------------
#' Lit patients.csv.gz pour récupérer anchor_year_group (non lu par
#' 06_define_cohort.R en v2 — seule anchor_year l'était, d'où le bug §1.13).
read_anchor_year_group <- function(tbl_paths) {
  paths <- tbl_paths$paths
  pat <- read_mimic_table(paths[["patients"]],
                           col_select = c("subject_id", "anchor_year_group"))
  pat[, subject_id := as.integer(subject_id)]
  message("[32_cov_v3] anchor_year_group lu pour ", nrow(pat), " patients : ",
          paste(sort(unique(pat$anchor_year_group)), collapse = " | "))
  pat
}

# -----------------------------------------------------------------------
# §2.2 — Comorbidités / Charlson / sepsis_icd à partir des hospitalisations
#         ANTÉRIEURES uniquement (jamais l'hadm_id index, codé à la sortie)
# -----------------------------------------------------------------------
#' @param cohort     cohorte v3 en cours de construction (subject_id, hadm_id,
#'                   T0, admittime = date d'admission de l'épisode index)
#' @param tbl_paths  sortie de detect_mimic_tables
#' @param icd_codes  sortie de get_icd_codes()
build_icd_covariates_prior_v3 <- function(cohort, tbl_paths, icd_codes) {
  paths <- tbl_paths$paths
  result <- cohort[, .(subject_id, hadm_id)]

  # Toutes les admissions du patient (nécessaire pour dater les hadm_id
  # antérieurs à l'épisode index)
  adm_all <- read_mimic_table(paths[["admissions"]],
                               col_select = c("subject_id", "hadm_id", "admittime"))
  adm_all[, subject_id := as.integer(subject_id)]
  adm_all[, hadm_id    := as.integer(hadm_id)]
  adm_all[, admittime  := lubridate::ymd_hms(admittime, quiet = TRUE)]

  index_admit <- cohort[, .(subject_id, hadm_id_index = hadm_id,
                             admittime_index = admittime)]
  adm_prior <- merge(adm_all, index_admit, by = "subject_id", all = FALSE,
                      allow.cartesian = TRUE)
  # Hospitalisations STRICTEMENT antérieures à l'admission index (jamais l'index
  # lui-même, et jamais une hospitalisation qui commence après ou pendant)
  adm_prior <- adm_prior[hadm_id != hadm_id_index & admittime < admittime_index]
  prior_hadm_ids <- unique(adm_prior[, .(subject_id, hadm_id)])
  n_with_prior <- uniqueN(prior_hadm_ids$subject_id)
  message("[32_cov_v3] Patients avec >=1 hospitalisation MIMIC antérieure : ",
          n_with_prior, "/", nrow(cohort), " (",
          round(n_with_prior / nrow(cohort) * 100, 1), "%) — les autres ont ",
          "des comorbidités ICD non évaluables (mises à 0, sous-estimation ",
          "documentée en limite : MIMIC ne voit pas les hospitalisations hors ",
          "de ce centre).")

  if (is.na(paths[["diagnoses_icd"]]) || nrow(prior_hadm_ids) == 0) {
    for (col in c("ckd", "diabetes", "heart_failure", "liver_disease",
                  "cancer", "immunosuppression", "sepsis_icd", "esrd_icd_prior",
                  "charlson_score")) {
      result[, (col) := 0L]
    }
    result[, has_prior_hospitalization := as.integer(subject_id %in% prior_hadm_ids$subject_id)]
    return(result)
  }

  diag <- read_mimic_table(paths[["diagnoses_icd"]],
                            col_select = c("subject_id", "hadm_id", "icd_code", "icd_version"))
  diag[, subject_id := as.integer(subject_id)]
  diag[, hadm_id    := as.integer(hadm_id)]
  diag[, icd_code   := trimws(icd_code)]
  # Restreindre AUX hadm_id antérieurs uniquement — c'est la correction
  diag_prior <- diag[prior_hadm_ids, on = .(subject_id, hadm_id), nomatch = NULL]
  message("[32_cov_v3] Lignes ICD retenues (hospitalisations antérieures) : ",
          nrow(diag_prior), " / ", nrow(diag), " au total dans diagnoses_icd")

  flag_icd <- function(codes_icd9, codes_icd10) {
    diag_prior[
      (icd_version == 9 & icd_code %in% codes_icd9) |
        (icd_version == 10 & Reduce(`|`, lapply(codes_icd10, function(p) startsWith(icd_code, p)))),
      .(flag = 1L), by = subject_id
    ]
  }

  for (cond in c("ckd", "diabetes", "heart_failure", "liver_disease",
                 "cancer", "immunosuppression", "sepsis")) {
    col_name <- if (cond == "sepsis") "sepsis_icd" else cond
    codes <- icd_codes[[cond]]
    flags <- flag_icd(codes$icd9, codes$icd10)
    setnames(flags, "flag", col_name)
    result <- merge(result, flags[, c("subject_id", col_name), with = FALSE],
                     by = "subject_id", all.x = TRUE)
    result[is.na(get(col_name)), (col_name) := 0L]
  }

  esrd_codes <- icd_codes$esrd
  esrd_flags <- flag_icd(esrd_codes$icd9, esrd_codes$icd10)
  setnames(esrd_flags, "flag", "esrd_icd_prior")
  result <- merge(result, esrd_flags, by = "subject_id", all.x = TRUE)
  result[is.na(esrd_icd_prior), esrd_icd_prior := 0L]

  # Charlson (codes antérieurs uniquement)
  charlson_codes <- get_charlson_codes()
  charlson_dt <- unique(diag_prior[, .(subject_id)])
  for (cond in names(charlson_codes)) {
    codes <- charlson_codes[[cond]]
    flags <- flag_icd(codes$icd9, codes$icd10)
    setnames(flags, "flag", cond)
    charlson_dt <- merge(charlson_dt, flags, by = "subject_id", all.x = TRUE)
    charlson_dt[is.na(get(cond)), (cond) := 0L]
  }
  charlson_dt[, charlson_score := 0L]
  for (cond in names(CHARLSON_WEIGHTS)) {
    if (cond %in% names(charlson_dt)) {
      charlson_dt[, charlson_score := charlson_score + CHARLSON_WEIGHTS[[cond]] * get(cond)]
    }
  }
  result <- merge(result, charlson_dt[, .(subject_id, charlson_score)],
                   by = "subject_id", all.x = TRUE)
  result[is.na(charlson_score), charlson_score := 0L]

  result[, has_prior_hospitalization := as.integer(subject_id %in% prior_hadm_ids$subject_id)]
  message("[32_cov_v3] Comorbidités (hospitalisations antérieures) — prévalences : ",
          "ckd=", round(mean(result$ckd) * 100, 1), "%, diabetes=",
          round(mean(result$diabetes) * 100, 1), "%, cancer=",
          round(mean(result$cancer) * 100, 1), "%, charlson médian=",
          median(result$charlson_score))
  result
}

# -----------------------------------------------------------------------
# §2.3 — Microbiologie : storetime (résultat disponible) au lieu de
#         charttime (prélèvement) + antécédent BLSE/résistance C3G (12 mois,
#         toutes hospitalisations antérieures incluses) + site infectieux
#         présumé (§2.5, dérivé de spec_type_desc)
# -----------------------------------------------------------------------
build_microbiology_v3 <- function(cohort, tbl_paths) {
  paths <- tbl_paths$paths
  result <- cohort[, .(subject_id, hadm_id)]
  if (is.na(paths[["microbiologyevents"]])) {
    result[, `:=`(positive_culture_pre_t0_storetime = 0L,
                   esbl_resistant_history_12mo = 0L,
                   infection_site_proxy = "unknown")]
    return(result)
  }

  micro <- read_mimic_table(paths[["microbiologyevents"]],
                             col_select = c("subject_id", "hadm_id", "charttime",
                                            "storetime", "spec_type_desc",
                                            "org_name", "ab_name", "interpretation"))
  micro[, subject_id := as.integer(subject_id)]
  micro[, hadm_id    := as.integer(hadm_id)]
  micro[, charttime  := lubridate::ymd_hms(charttime, quiet = TRUE)]
  micro[, storetime  := lubridate::ymd_hms(storetime, quiet = TRUE)]

  micro_t0 <- merge(micro, cohort[, .(subject_id, hadm_id, T0)],
                     by = c("subject_id", "hadm_id"), all = FALSE, allow.cartesian = TRUE)

  # Positive culture avant T0, en n'utilisant QUE l'information disponible
  # (storetime) à T0 — corrige la fuite charttime de la v2
  pos_cult <- micro_t0[
    !is.na(T0) & !is.na(storetime) &
      storetime >= T0 - 48 * 3600 & storetime <= T0 &
      !is.na(org_name) & org_name != "",
    .(positive_culture_pre_t0_storetime = 1L), by = .(subject_id, hadm_id)
  ]
  result <- merge(result, pos_cult, by = c("subject_id", "hadm_id"), all.x = TRUE)
  result[is.na(positive_culture_pre_t0_storetime), positive_culture_pre_t0_storetime := 0L]

  # Site infectieux présumé : dernier spec_type_desc avec culture positive
  # dans les 48h avant T0 (toutes hospitalisations pertinentes, storetime <= T0)
  site_map <- micro_t0[
    !is.na(T0) & !is.na(storetime) & storetime >= T0 - 48 * 3600 & storetime <= T0 &
      !is.na(org_name) & org_name != "",
    .(infection_site_proxy = dplyr::case_when(
      grepl("blood", spec_type_desc, ignore.case = TRUE) ~ "bacteremia",
      grepl("urine", spec_type_desc, ignore.case = TRUE) ~ "urinary",
      grepl("sputum|bronch|lung|respiratory", spec_type_desc, ignore.case = TRUE) ~ "respiratory",
      grepl("abdo|perито|bile|peritoneal", spec_type_desc, ignore.case = TRUE) ~ "intra_abdominal",
      grepl("skin|wound|tissue|abscess", spec_type_desc, ignore.case = TRUE) ~ "skin_soft_tissue",
      TRUE ~ "other_unknown"
    )), by = .(subject_id, hadm_id, storetime)
  ]
  site_last <- site_map[order(subject_id, hadm_id, -storetime)][
    , .SD[1], by = .(subject_id, hadm_id)][, .(subject_id, hadm_id, infection_site_proxy)]
  result <- merge(result, site_last, by = c("subject_id", "hadm_id"), all.x = TRUE)
  result[is.na(infection_site_proxy), infection_site_proxy := "no_positive_culture_pre_t0"]

  # Antécédent BLSE / résistance C3G dans les 12 mois précédents (toutes
  # hospitalisations antérieures) : interprétation "R" à ceftriaxone,
  # ceftazidime, cefotaxime, cefepime, OU org_name mentionnant ESBL/BLSE
  index_t0 <- cohort[, .(subject_id, hadm_id_index = hadm_id, T0)]
  micro_hist <- merge(micro, index_t0, by = "subject_id", all = FALSE, allow.cartesian = TRUE)
  micro_hist <- micro_hist[hadm_id != hadm_id_index | is.na(hadm_id_index)]
  micro_hist[, days_before_t0 := as.numeric(difftime(T0, storetime, units = "days"))]
  esbl_flag <- micro_hist[
    !is.na(days_before_t0) & days_before_t0 > 0 & days_before_t0 <= 365 &
      ((toupper(interpretation) == "R" &
          grepl("ceftriaxone|ceftazidime|cefotaxime|cefepime", ab_name, ignore.case = TRUE)) |
         grepl("ESBL|extended.spectrum", org_name, ignore.case = TRUE)),
    .(esbl_resistant_history_12mo = 1L), by = subject_id
  ]
  result <- merge(result, esbl_flag, by = "subject_id", all.x = TRUE)
  result[is.na(esbl_resistant_history_12mo), esbl_resistant_history_12mo := 0L]

  message("[32_cov_v3] positive_culture_pre_t0 (storetime) : ",
          round(mean(result$positive_culture_pre_t0_storetime) * 100, 1),
          "% (vs charttime v2, fuite d'information corrigée)")
  message("[32_cov_v3] esbl_resistant_history_12mo : ",
          round(mean(result$esbl_resistant_history_12mo) * 100, 1), "%")
  print(result[, .N, by = infection_site_proxy][order(-N)])

  result
}

# -----------------------------------------------------------------------
# §2.2 — Sepsis à T0 : composite pragmatique (recommandation retenue §7.4)
#         — pas de code ICD de l'hospitalisation index. Suspicion d'infection
#         = par construction (bêta-lactame large spectre initié) ; dysfonction
#         d'organe = proxy composite (vasopresseurs, ventilation mécanique,
#         lactate élevé si disponible). Sensibilité ICU-only : score SOFA-lite
#         avec les composantes réellement mesurables (voir build_sofa_lite_icu_v3).
# -----------------------------------------------------------------------
build_sepsis_composite_v3 <- function(cohort) {
  dt <- copy(cohort)
  has_lactate <- "lactate_max" %in% names(dt)
  dt[, sepsis_composite_v3 := as.integer(
    vasopressor_pre_t0 == 1 | vent_pre_t0 == 1 |
      (has_lactate && !is.na(lactate_max) && lactate_max >= 2) |
      FALSE
  )]
  dt[is.na(sepsis_composite_v3), sepsis_composite_v3 := 0L]
  dt[, .(subject_id, hadm_id, sepsis_composite_v3)]
}

#' SOFA-lite, composantes disponibles uniquement (sensibilité ICU, taux de
#' calculabilité rapporté par bras plutôt que de fabriquer les composantes
#' manquantes — plaquettes/bilirubine ne sont pas extraites dans ce pipeline)
build_sofa_lite_icu_v3 <- function(cohort) {
  dt <- copy(cohort)
  dt[, sofa_cv   := fifelse(vasopressor_pre_t0 == 1, 2L, fifelse(!is.na(map_mean) & map_mean < 70, 1L, 0L))]
  dt[, sofa_resp := fifelse(vent_pre_t0 == 1, 2L, 0L)]
  dt[, sofa_cns  := fifelse(!is.na(gcs_min) & icu_before_t0 == 1,
                             fifelse(gcs_min <= 2, 3L, fifelse(gcs_min <= 3, 1L, 0L)), NA_integer_)]
  dt[, sofa_renal := fifelse(!is.na(baseline_creat_min7d),
                              fifelse(baseline_creat_min7d >= 5, 4L,
                                      fifelse(baseline_creat_min7d >= 3.5, 3L,
                                              fifelse(baseline_creat_min7d >= 2, 2L,
                                                      fifelse(baseline_creat_min7d >= 1.2, 1L, 0L)))),
                              NA_integer_)]
  dt[, sofa_lite_calculable := as.integer(icu_before_t0 == 1 & !is.na(sofa_cns))]
  dt[, sofa_lite_score := sofa_cv + sofa_resp + fifelse(is.na(sofa_cns), 0L, sofa_cns) +
       fifelse(is.na(sofa_renal), 0L, sofa_renal)]
  message("[32_cov_v3] SOFA-lite calculable (ICU + GCS dispo) pour ",
          round(mean(dt$sofa_lite_calculable) * 100, 1), "% de la cohorte")
  dt[, .(subject_id, hadm_id, sofa_lite_score, sofa_lite_calculable)]
}

# -----------------------------------------------------------------------
# §2.5 — Nouvelles covariables (celles ne nécessitant pas de ré-extraction :
#         neutropénie, hémopathie/greffe [ICD antérieurs], lactate max,
#         poids/IMC [omr], service à T0 [transfers/services])
# -----------------------------------------------------------------------
build_new_covariates_v3 <- function(cohort, tbl_paths) {
  paths <- tbl_paths$paths
  result <- cohort[, .(subject_id, hadm_id, T0)]

  # --- Neutropénie (PNN < 500) : labevents itemid 51256 (Neutrophils, %) et
  #     52075 (Absolute Neutrophil Count) — on utilise l'ANC absolu si dispo
  if (!is.na(paths[["labevents"]]) && !is.na(paths[["d_labitems"]])) {
    lab <- read_mimic_table(paths[["labevents"]],
                             col_select = c("subject_id", "hadm_id", "itemid",
                                            "charttime", "valuenum"),
                             use_duckdb = TRUE)
    lab[, subject_id := as.integer(subject_id)]
    lab[, hadm_id    := as.integer(hadm_id)]
    lab[, charttime  := lubridate::ymd_hms(charttime, quiet = TRUE)]

    # itemid 52075 (Absolute Neutrophil Count) est rapporté en K/uL dans
    # labevents (médiane ~4.48 K/uL, vérifié via d_labitems/valueuom), pas en
    # cellules/uL : le seuil clinique de neutropénie sévère (<500 cellules/uL)
    # correspond donc à valuenum < 0.5, PAS < 500 (sans cette conversion,
    # valuenum < 500 est vrai pour la quasi-totalité des mesures -> la
    # variable capturait "ANC mesuré dans les 48h", pas la neutropénie).
    anc <- lab[itemid == 52075 & !is.na(valuenum)]
    anc_t0 <- merge(anc, cohort[, .(subject_id, hadm_id, T0)],
                     by = c("subject_id", "hadm_id"), all = FALSE, allow.cartesian = TRUE)
    neutropenia <- anc_t0[
      charttime <= T0 & charttime >= T0 - 48 * 3600 & valuenum < 0.5,
      .(neutropenia_pre_t0 = 1L), by = .(subject_id, hadm_id)]
    result <- merge(result, neutropenia, by = c("subject_id", "hadm_id"), all.x = TRUE)
    result[is.na(neutropenia_pre_t0), neutropenia_pre_t0 := 0L]

    # --- Lactate max (itemid 50813) dans les 24h avant T0
    lact <- lab[itemid == 50813 & !is.na(valuenum) & valuenum > 0 & valuenum < 30]
    lact_t0 <- merge(lact, cohort[, .(subject_id, hadm_id, T0)],
                      by = c("subject_id", "hadm_id"), all = FALSE, allow.cartesian = TRUE)
    lact_max <- lact_t0[
      charttime <= T0 & charttime >= T0 - 24 * 3600,
      .(lactate_max = max(valuenum), lactate_measured = 1L), by = .(subject_id, hadm_id)]
    result <- merge(result, lact_max, by = c("subject_id", "hadm_id"), all.x = TRUE)
    result[is.na(lactate_measured), lactate_measured := 0L]
    rm(lab); gc()
  } else {
    result[, `:=`(neutropenia_pre_t0 = 0L, lactate_max = NA_real_, lactate_measured = 0L)]
  }

  # --- Hémopathie maligne / greffe (ICD ANTÉRIEURS uniquement, cohérent §2.2)
  if (!is.na(paths[["diagnoses_icd"]])) {
    adm_all <- read_mimic_table(paths[["admissions"]],
                                 col_select = c("subject_id", "hadm_id", "admittime"))
    adm_all[, subject_id := as.integer(subject_id)]
    adm_all[, hadm_id    := as.integer(hadm_id)]
    adm_all[, admittime  := lubridate::ymd_hms(admittime, quiet = TRUE)]
    index_admit <- cohort[, .(subject_id, hadm_id_index = hadm_id, admittime_index = admittime)]
    adm_prior <- merge(adm_all, index_admit, by = "subject_id", all = FALSE, allow.cartesian = TRUE)
    adm_prior <- adm_prior[hadm_id != hadm_id_index & admittime < admittime_index]
    prior_ids <- unique(adm_prior[, .(subject_id, hadm_id)])

    diag <- read_mimic_table(paths[["diagnoses_icd"]],
                              col_select = c("subject_id", "hadm_id", "icd_code", "icd_version"))
    diag[, subject_id := as.integer(subject_id)]
    diag[, hadm_id    := as.integer(hadm_id)]
    diag_prior <- diag[prior_ids, on = .(subject_id, hadm_id), nomatch = NULL]

    heme_icd9  <- c("200","201","202","203","204","205","206","207","208")
    heme_icd10 <- c("C81","C82","C83","C84","C85","C88","C90","C91","C92","C93","C94","C95","C96")
    transplant_icd9  <- c("V420","V421","V427")
    transplant_icd10 <- c("Z940","Z941","Z947")

    heme <- diag_prior[
      (icd_version == 9 & icd_code %in% heme_icd9) |
        (icd_version == 10 & Reduce(`|`, lapply(heme_icd10, function(p) startsWith(icd_code, p)))),
      .(hematologic_malignancy = 1L), by = subject_id]
    transpl <- diag_prior[
      (icd_version == 9 & icd_code %in% transplant_icd9) |
        (icd_version == 10 & Reduce(`|`, lapply(transplant_icd10, function(p) startsWith(icd_code, p)))),
      .(transplant_history = 1L), by = subject_id]

    result <- merge(result, heme, by = "subject_id", all.x = TRUE)
    result <- merge(result, transpl, by = "subject_id", all.x = TRUE)
    result[is.na(hematologic_malignancy), hematologic_malignancy := 0L]
    result[is.na(transplant_history), transplant_history := 0L]
  } else {
    result[, `:=`(hematologic_malignancy = 0L, transplant_history = 0L)]
  }

  # --- Poids / IMC (omr) — mesure la plus proche avant T0, sur 365j
  if (!is.na(paths[["omr"]])) {
    omr <- read_mimic_table(paths[["omr"]],
                             col_select = c("subject_id", "chartdate", "result_name", "result_value"))
    omr[, subject_id := as.integer(subject_id)]
    omr[, chartdate := lubridate::ymd(chartdate, quiet = TRUE)]
    wt <- omr[result_name %in% c("Weight (Lbs)", "Weight")]
    wt[, value_kg := suppressWarnings(as.numeric(result_value)) *
         ifelse(result_name == "Weight (Lbs)", 0.453592, 1)]
    bmi <- omr[result_name == "BMI (kg/m2)"]
    bmi[, value_bmi := suppressWarnings(as.numeric(result_value))]

    T0_dates <- cohort[, .(subject_id, T0, t0_date = as.Date(T0))]
    wt_t0 <- merge(wt[!is.na(value_kg)], T0_dates, by = "subject_id", all = FALSE, allow.cartesian = TRUE)
    wt_t0[, days_before := as.numeric(difftime(t0_date, chartdate, units = "days"))]
    wt_best <- wt_t0[days_before >= 0 & days_before <= 365][
      order(subject_id, days_before)][, .SD[1], by = subject_id][, .(subject_id, weight_kg = value_kg)]

    bmi_t0 <- merge(bmi[!is.na(value_bmi)], T0_dates, by = "subject_id", all = FALSE, allow.cartesian = TRUE)
    bmi_t0[, days_before := as.numeric(difftime(t0_date, chartdate, units = "days"))]
    bmi_best <- bmi_t0[days_before >= 0 & days_before <= 365][
      order(subject_id, days_before)][, .SD[1], by = subject_id][, .(subject_id, bmi = value_bmi)]

    result <- merge(result, wt_best, by = "subject_id", all.x = TRUE)
    result <- merge(result, bmi_best, by = "subject_id", all.x = TRUE)
  } else {
    result[, `:=`(weight_kg = NA_real_, bmi = NA_real_)]
  }
  result[, weight_measured := as.integer(!is.na(weight_kg))]
  result[, bmi_measured    := as.integer(!is.na(bmi))]

  # --- Service / unité à T0 (transfers)
  if (!is.na(paths[["transfers"]])) {
    tr <- read_mimic_table(paths[["transfers"]],
                            col_select = c("subject_id", "hadm_id", "eventtype",
                                           "careunit", "intime", "outtime"))
    tr[, subject_id := as.integer(subject_id)]
    tr[, hadm_id    := as.integer(hadm_id)]
    tr[, intime  := lubridate::ymd_hms(intime,  quiet = TRUE)]
    tr[, outtime := lubridate::ymd_hms(outtime, quiet = TRUE)]
    tr_t0 <- merge(tr[eventtype != "discharge"], cohort[, .(subject_id, hadm_id, T0)],
                    by = c("subject_id", "hadm_id"), all = FALSE, allow.cartesian = TRUE)
    tr_at_t0 <- tr_t0[intime <= T0 & (is.na(outtime) | outtime >= T0)]
    svc_t0 <- tr_at_t0[order(subject_id, hadm_id, -intime)][
      , .SD[1], by = .(subject_id, hadm_id)][, .(subject_id, hadm_id,
                                                    service_at_t0 = careunit)]
    result <- merge(result, svc_t0, by = c("subject_id", "hadm_id"), all.x = TRUE)
    result[is.na(service_at_t0), service_at_t0 := "unknown_or_ward"]
  } else {
    result[, service_at_t0 := "unknown"]
  }

  result[, T0 := NULL]
  message("[32_cov_v3] Nouvelles covariables : neutropenia=",
          round(mean(result$neutropenia_pre_t0) * 100, 1),
          "%, hematologic_malignancy=", round(mean(result$hematologic_malignancy) * 100, 1),
          "%, transplant_history=", round(mean(result$transplant_history) * 100, 1),
          "%, lactate mesuré=", round(mean(result$lactate_measured) * 100, 1),
          "%, poids mesuré=", round(mean(result$weight_measured) * 100, 1), "%")
  result
}

# -----------------------------------------------------------------------
# §2.5 — Néphrotoxiques décomposés + ACEI/ARB + NSAID (à partir de
#         medications_all_v3.rds — nécessite la ré-extraction v3)
# -----------------------------------------------------------------------
build_nephrotoxin_covariates_v3 <- function(meds_v3, T0_dt, window_h_before = 48) {
  result <- copy(T0_dt[, .(subject_id, hadm_id, T0)])

  add_class_flag <- function(result, meds_v3, class_filter, col_name, window_h = window_h_before) {
    dt_cls <- meds_v3[class_filter, .(subject_id, hadm_id, starttime)]
    if (nrow(dt_cls) == 0) { result[, (col_name) := 0L]; return(result) }
    merged <- merge(dt_cls, T0_dt[, .(subject_id, hadm_id, T0)],
                     by = c("subject_id", "hadm_id"), all = FALSE, allow.cartesian = TRUE)
    flags <- merged[!is.na(T0) & !is.na(starttime) &
                       starttime >= T0 - window_h * 3600 & starttime < T0,
                     .(flag = 1L), by = .(subject_id, hadm_id)]
    setnames(flags, "flag", col_name)
    result <- merge(result, flags, by = c("subject_id", "hadm_id"), all.x = TRUE)
    result[is.na(get(col_name)), (col_name) := 0L]
    result
  }

  result <- add_class_flag(result, meds_v3,
                            meds_v3$drug_class == "other_nephrotoxin" &
                              meds_v3$nephrotoxin_subclass == "direct_nephrotoxin",
                            "direct_nephrotoxin_pre_t0")
  result <- add_class_flag(result, meds_v3,
                            meds_v3$drug_class == "other_nephrotoxin" &
                              meds_v3$nephrotoxin_subclass == "calcineurin_inhibitor",
                            "calcineurin_inhibitor_pre_t0")
  result <- add_class_flag(result, meds_v3,
                            meds_v3$drug_class == "other_nephrotoxin" &
                              meds_v3$nephrotoxin_subclass == "loop_diuretic",
                            "loop_diuretic_pre_t0")
  result <- add_class_flag(result, meds_v3, meds_v3$drug_class == "acei_arb",
                            "acei_arb_pre_t0")
  result <- add_class_flag(result, meds_v3, meds_v3$drug_class == "nsaid",
                            "nsaid_pre_t0")

  result[, T0 := NULL]
  message("[32_cov_v3] Néphrotoxiques décomposés : direct=",
          round(mean(result$direct_nephrotoxin_pre_t0) * 100, 1),
          "%, calcineurine=", round(mean(result$calcineurin_inhibitor_pre_t0) * 100, 1),
          "%, diurétique de l'anse=", round(mean(result$loop_diuretic_pre_t0) * 100, 1),
          "%, IEC/ARA2=", round(mean(result$acei_arb_pre_t0) * 100, 1),
          "%, AINS=", round(mean(result$nsaid_pre_t0) * 100, 1), "%")
  result
}
