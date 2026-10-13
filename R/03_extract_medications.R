# =============================================================================
# 03_extract_medications.R — Extraction des médicaments depuis MIMIC-IV
# Priorité : emar_detail > emar > pharmacy > prescriptions
# =============================================================================

#' Extraire toutes les administrations de médicaments d'intérêt
#' @param tbl_paths   liste des chemins (sortie de detect_mimic_tables)
#' @param cfg         configuration (sortie de load_config)
#' @return data.table avec colonnes :
#'   subject_id, hadm_id, starttime, drug_name, drug_class
#'   (ptz | carbapenem | vancomycin | aminoglycoside | other_nephrotoxin)
extract_medications <- function(tbl_paths, cfg) {

  paths      <- tbl_paths$paths
  med_source <- tbl_paths$med_source

  # Patterns regex (case-insensitive)
  ptz_pat   <- paste(cfg$medications$ptz_patterns,   collapse = "|")
  carb_pat  <- paste(cfg$medications$carbapenem_patterns, collapse = "|")
  vanco_pat <- paste(cfg$medications$vanco_patterns, collapse = "|")
  amino_pat <- paste(cfg$medications$aminoglycoside_patterns, collapse = "|")
  nephr_pat <- paste(cfg$medications$other_nephrotoxins, collapse = "|")

  classify_drug <- function(drug_name) {
    dn <- tolower(drug_name)
    dplyr::case_when(
      grepl(ptz_pat,   dn, ignore.case = TRUE) &
        !grepl("^tazobactam$", dn)             ~ "ptz",
      grepl(carb_pat,  dn, ignore.case = TRUE) ~ "carbapenem",
      grepl(vanco_pat, dn, ignore.case = TRUE) ~ "vancomycin",
      grepl(amino_pat, dn, ignore.case = TRUE) ~ "aminoglycoside",
      grepl(nephr_pat, dn, ignore.case = TRUE) ~ "other_nephrotoxin",
      TRUE                                     ~ NA_character_
    )
  }

  meds <- NULL

  # -----------------------------------------------------------
  # Stratégie 1 : emar_detail + emar (administrations réelles)
  # -----------------------------------------------------------
  if (med_source %in% c("emar", "emar_detail") &&
      !is.na(paths[["emar"]]) && !is.na(paths[["emar_detail"]])) {

    message("[03_meds] Utilisation de emar + emar_detail (administrations réelles)")

    emar <- read_mimic_table(paths[["emar"]],
                              col_select = c("subject_id", "hadm_id", "emar_id",
                                             "medication", "scheduletime",
                                             "charttime", "event_txt"))
    ed   <- read_mimic_table(paths[["emar_detail"]],
                              col_select = c("subject_id", "emar_id",
                                             "product_description"))

    if (!is.null(emar) && !is.null(ed)) {
      # ---- OPTIMISATION : pre-filtrer emar avant la jointure avec emar_detail ----
      # Construire le pattern global une seule fois
      all_pat <- paste(unique(c(
        cfg$medications$ptz_patterns, cfg$medications$carbapenem_patterns,
        cfg$medications$vanco_patterns, cfg$medications$aminoglycoside_patterns,
        cfg$medications$other_nephrotoxins
      )), collapse = "|")

      # Filtrer emar par le champ medication (rapide, avant jointure)
      emar_filt <- emar[grepl(all_pat, medication, ignore.case = TRUE)]
      message("[03_meds] emar après pre-filtrage par nom : ",
              nrow(emar_filt), " lignes (sur ", nrow(emar), ")")

      # Filtrer emar_detail pour ne garder que les emar_id d'intérêt
      # Seulement emar_id + product_description (évite le conflit subject_id.x)
      ed_filt <- ed[emar_id %in% emar_filt$emar_id, .(emar_id, product_description)]
      message("[03_meds] emar_detail filtré : ", nrow(ed_filt), " lignes")

      # Jointure sur les sous-ensembles filtrés (beaucoup plus rapide)
      emar_filt[, drug_name := medication]
      emar_joined <- merge(emar_filt, ed_filt, by = "emar_id", all.x = TRUE)
      emar_joined[!is.na(product_description) & product_description != "",
                  drug_name := product_description]

      # Libérer la mémoire des grands tableaux maintenant inutiles
      rm(emar, ed, ed_filt)
      gc()

      emar_joined[, starttime := dplyr::coalesce(
        lubridate::ymd_hms(charttime, quiet = TRUE),
        lubridate::ymd_hms(scheduletime, quiet = TRUE)
      )]

      emar_joined[, drug_class := classify_drug(drug_name)]
      meds <- emar_joined[!is.na(drug_class),
                           .(subject_id, hadm_id, starttime, drug_name, drug_class,
                             source = "emar")]

    } else if (!is.null(emar)) {
      # emar sans emar_detail
      emar[, starttime := lubridate::ymd_hms(charttime, quiet = TRUE)]
      emar[, drug_class := classify_drug(medication)]
      meds <- emar[!is.na(drug_class),
                    .(subject_id, hadm_id, starttime,
                      drug_name = medication, drug_class, source = "emar")]
    }
  }

  # -----------------------------------------------------------
  # Stratégie 2 : pharmacy
  # -----------------------------------------------------------
  if (is.null(meds) && !is.na(paths[["pharmacy"]])) {
    message("[03_meds] Utilisation de pharmacy")
    ph <- read_mimic_table(paths[["pharmacy"]],
                            col_select = c("subject_id", "hadm_id",
                                           "medication", "starttime", "stoptime"))
    if (!is.null(ph)) {
      ph[, starttime := lubridate::ymd_hms(starttime, quiet = TRUE)]
      ph[, drug_class := classify_drug(medication)]
      meds <- ph[!is.na(drug_class),
                  .(subject_id, hadm_id, starttime,
                    drug_name = medication, drug_class, source = "pharmacy")]
    }
  }

  # -----------------------------------------------------------
  # Stratégie 3 : prescriptions (dates ordonnances, moins précis)
  # -----------------------------------------------------------
  if (is.null(meds) && !is.na(paths[["prescriptions"]])) {
    message("[03_meds] Utilisation de prescriptions (fallback — dates d'ordonnance)")
    pr <- read_mimic_table(paths[["prescriptions"]],
                            col_select = c("subject_id", "hadm_id", "drug",
                                           "starttime", "stoptime"))
    if (!is.null(pr)) {
      pr[, starttime := lubridate::ymd_hms(starttime, quiet = TRUE)]
      pr[, drug_class := classify_drug(drug)]
      meds <- pr[!is.na(drug_class),
                  .(subject_id, hadm_id, starttime,
                    drug_name = drug, drug_class, source = "prescriptions")]
    }
  }

  if (is.null(meds)) stop("[03_meds] Aucune source de médicaments disponible.")

  # Nettoyer
  meds <- meds[!is.na(starttime)]
  meds[, subject_id := as.integer(subject_id)]
  meds[, hadm_id    := as.integer(hadm_id)]
  setorder(meds, subject_id, hadm_id, starttime)

  # Identifier la molécule exacte pour les carbapénèmes
  meds[drug_class == "carbapenem",
       carbapenem_molecule := dplyr::case_when(
         grepl("meropenem",   drug_name, ignore.case = TRUE) ~ "meropenem",
         grepl("imipenem",    drug_name, ignore.case = TRUE) ~ "imipenem",
         grepl("doripenem",   drug_name, ignore.case = TRUE) ~ "doripenem",
         grepl("ertapenem",   drug_name, ignore.case = TRUE) ~ "ertapenem",
         TRUE                                                ~ "carbapenem_other"
       )]

  message("[03_meds] Administrations extraites : ", nrow(meds))
  message("[03_meds] Répartition par classe :")
  print(meds[, .N, by = drug_class])

  save_derived(meds, "medications_all", format = "rds")
  meds
}

#' Identifier T0 = première administration PTZ ou carbapénème par épisode
#' @param meds  data.table sortie de extract_medications
#' @param exclude_ertapenem  exclure ertapénème de la comparaison principale
#' @return data.table (subject_id, hadm_id, T0, treatment, molecule,
#'                     erta_sensitivity)
identify_t0 <- function(meds, exclude_ertapenem = TRUE) {

  # PTZ : première dose
  ptz_first <- meds[drug_class == "ptz",
                     .(T0 = min(starttime, na.rm = TRUE), treatment = "PTZ",
                       molecule = drug_name[which.min(starttime)]),
                     by = .(subject_id, hadm_id)]

  # Carbapénème : première dose (éventuellement sans ertapénème)
  carb_dt <- meds[drug_class == "carbapenem"]

  carb_all_first <- carb_dt[,
                              .(T0 = min(starttime, na.rm = TRUE),
                                molecule = drug_name[which.min(starttime)]),
                              by = .(subject_id, hadm_id)]
  carb_all_first[, treatment := "carbapenem_all"]

  if (exclude_ertapenem) {
    carb_no_erta <- carb_dt[!grepl("ertapenem", drug_name, ignore.case = TRUE)]
    if (nrow(carb_no_erta) > 0) {
      carb_noerta_first <- carb_no_erta[,
                                         .(T0 = min(starttime, na.rm = TRUE),
                                           molecule = drug_name[which.min(starttime)]),
                                         by = .(subject_id, hadm_id)]
      carb_noerta_first[, treatment := "carbapenem_noerta"]
    } else {
      carb_noerta_first <- carb_all_first[, .(subject_id, hadm_id, T0, molecule,
                                               treatment = "carbapenem_noerta")]
    }
  }

  # Carbapenem meropenem-only
  carb_mero <- carb_dt[grepl("meropenem", drug_name, ignore.case = TRUE)]
  if (nrow(carb_mero) > 0) {
    carb_mero_first <- carb_mero[,
                                   .(T0 = min(starttime, na.rm = TRUE),
                                     molecule = drug_name[which.min(starttime)]),
                                   by = .(subject_id, hadm_id)]
    carb_mero_first[, treatment := "meropenem_only"]
  } else {
    carb_mero_first <- data.table(subject_id = integer(0), hadm_id = integer(0),
                                   T0 = as.POSIXct(character(0)),
                                   molecule = character(0),
                                   treatment = character(0))
  }

  # Pour l'analyse principale : PTZ vs carbapenem_noerta (si exclude_ertapenem)
  # Un patient peut être dans les deux groupes → confusion d'exposition
  # Règle : new-user design = premier épisode éligible
  # On retourne toutes les versions pour les analyses de vibration

  list(
    ptz             = ptz_first,
    carb_all        = carb_all_first,
    carb_noerta     = if (exclude_ertapenem) carb_noerta_first else carb_all_first,
    carb_mero_only  = carb_mero_first
  )
}

#' Identifier la co-administration de vancomycine par fenêtre temporelle
#' @param meds  data.table
#' @param T0_dt data.table avec subject_id, hadm_id, T0
#' @param windows_h vecteur de fenêtres en heures (ex: c(6,12,24,48))
#' @return data.table avec colonnes vanco_Nh pour chaque fenêtre
add_vancomycin_windows <- function(meds, T0_dt, windows_h = c(6, 12, 24, 48)) {

  vanco_dt <- meds[drug_class == "vancomycin",
                    .(subject_id, hadm_id, starttime)]

  result <- copy(T0_dt)

  for (w in windows_h) {
    col_name <- paste0("vanco_", w, "h")
    # Merge vancomycin avec T0 pour avoir T0 dans le même data.table
    vanco_merged <- merge(vanco_dt, T0_dt[, .(subject_id, hadm_id, T0)],
                           by = c("subject_id", "hadm_id"),
                           all = FALSE, allow.cartesian = TRUE)

    vanco_in_window <- vanco_merged[
      !is.na(T0) & !is.na(starttime) &
        starttime >= T0 - 24 * 3600 &
        starttime <= T0 + w  * 3600,
      .(vanco_flag = 1L),
      by = .(subject_id, hadm_id)
    ]

    setnames(vanco_in_window, "vanco_flag", col_name)
    result <- merge(result, vanco_in_window,
                    by = c("subject_id", "hadm_id"), all.x = TRUE)
    result[is.na(get(col_name)), (col_name) := 0L]
  }

  result
}

#' Identifier autres néphrotoxiques avant T0
add_other_nephrotoxins <- function(meds, T0_dt, window_h_before = 48) {

  result <- copy(T0_dt)

  for (cls in c("aminoglycoside", "other_nephrotoxin")) {
    col_name <- paste0(cls, "_pre_t0")
    dt_cls <- meds[drug_class == cls, .(subject_id, hadm_id, starttime)]

    if (nrow(dt_cls) == 0) {
      result[, (col_name) := 0L]
      next
    }

    dt_merged <- merge(dt_cls, T0_dt[, .(subject_id, hadm_id, T0)],
                        by = c("subject_id", "hadm_id"),
                        all = FALSE, allow.cartesian = TRUE)
    flags <- dt_merged[
      !is.na(T0) & !is.na(starttime) &
        starttime >= T0 - window_h_before * 3600 & starttime < T0,
      .(flag = 1L),
      by = .(subject_id, hadm_id)
    ]

    setnames(flags, "flag", col_name)
    result <- merge(result, flags, by = c("subject_id", "hadm_id"), all.x = TRUE)
    result[is.na(get(col_name)), (col_name) := 0L]
  }

  result
}
