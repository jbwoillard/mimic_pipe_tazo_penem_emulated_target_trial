# =============================================================================
# 57_antipseudo_cohort.R — Administrations réelles d'antibiotiques et
#   cohorte unifiée des initiateurs PTZ / céfépime / méropénème
#   (docs/v3.1_plan.md B1-B3)
# =============================================================================
# Défaut v2/v3 identifié en phase B : medications_all_v3.rds garde TOUS les
# événements eMAR (Not Given, Stopped, Hold Dose, ...). Le T0 de 537 patients
# de la cohorte primaire (6,8 %) ne correspond pas à une administration ; 104
# n'ont jamais reçu le médicament index (99 PTZ, 5 carbapénèmes). Les
# cohortes construites ici n'utilisent que les administrations effectives.

EMAR_ADMIN_EVENTS_V31 <- c("Administered", "Delayed Administered", "Started", "Delayed Started",
                           "Administered in Other Location", "Started in Other Location",
                           "Partial Administered", "Restarted", "Documented in O.R. Holding")

# Antibactériens systémiques (les antinéoplasiques en -mycine sont exclus par nom)
ABX_REGEX_V31 <- paste0(
  "cillin|bactam|cef|ceph|penem|aztreonam|vancomycin|gentamicin|tobramycin|amikacin|streptomycin|",
  "floxacin|cycline|sulfamethoxazole|trimethoprim|metronidazole|clindamycin|linezolid|tedizolid|",
  "daptomycin|colistin|polymyxin|azithromycin|erythromycin|clarithromycin|rifampin|nitrofurantoin|",
  "fosfomycin|dalbavancin|oritavancin|telavancin|chloramphenicol|fidaxomicin")
ABX_EXCLUDE_V31 <- "bleomycin|mitomycin|dactinomycin|daunorubicin|doxorubicin|ophth|otic|topical|ointment|cream|neomycin|bacitracin"

classify_abx_v31 <- function(med) {
  m <- tolower(med)
  fcase(grepl("piperacillin|zosyn", m), "ptz",
        grepl("meropenem", m), "meropenem",
        grepl("imipenem|ertapenem|doripenem", m), "carbapenem_other",
        grepl("cefepime", m), "cefepime",
        grepl("ceftazidime", m), "ceftazidime",
        grepl("aztreonam", m), "aztreonam",
        grepl("vancomycin", m), "vancomycin",
        grepl("gentamicin|tobramycin|amikacin|streptomycin", m), "aminoglycoside",
        default = "other_antibacterial")
}

#' Administrations eMAR effectives d'antibactériens, avec la voie (pharmacy)
extract_administered_antibiotics_v31 <- function(tbl_paths,
                                                 cache = "output/extraction/derived_data/abx_administered_v31.rds") {
  if (file.exists(cache)) return(readRDS(cache))
  t0 <- Sys.time()
  con <- duckdb::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(duckdb::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  sql <- sprintf(
    "SELECT e.subject_id, e.hadm_id, coalesce(e.charttime, e.scheduletime) AS t, e.medication, e.event_txt, p.route
     FROM read_csv_auto('%s') e LEFT JOIN read_csv_auto('%s') p ON e.pharmacy_id = p.pharmacy_id
     WHERE regexp_matches(lower(e.medication), '%s') AND NOT regexp_matches(lower(e.medication), '%s')
       AND e.event_txt IN (%s)",
    tbl_paths$paths[["emar"]], tbl_paths$paths[["pharmacy"]], ABX_REGEX_V31, ABX_EXCLUDE_V31,
    paste0("'", EMAR_ADMIN_EVENTS_V31, "'", collapse = ", "))
  ab <- as.data.table(DBI::dbGetQuery(con, sql))
  ab[, `:=`(subject_id = as.integer(subject_id), hadm_id = as.integer(hadm_id),
            t = as.POSIXct(t, tz = "UTC"))]
  ab <- ab[!is.na(t)]
  ab[, abx_class := classify_abx_v31(medication)]
  # Voie : IV si pharmacy.route = IV ; voie inconnue -> IV pour les classes
  # sans forme orale/locale utilisée en pratique (bêta-lactamines
  # anti-pyocyaniques, carbapénèmes, aminosides)
  ab[, iv := as.integer(toupper(route) %in% c("IV", "IV DRIP", "IVPB", "IV BOLUS") |
                          (is.na(route) & abx_class %in% c("ptz", "meropenem", "carbapenem_other", "cefepime",
                                                            "ceftazidime", "aztreonam", "aminoglycoside")))]
  saveRDS(ab, cache)
  message("[57_abx_v31] ", nrow(ab), " administrations d'antibactériens (",
          round(difftime(Sys.time(), t0, units = "mins"), 1), " min) ; IV : ", round(100 * mean(ab$iv), 1), "%")
  ab
}

#' Cohorte unifiée : initiation (première ADMINISTRATION effective) de PTZ,
#' céfépime ou méropénème. Règles :
#'  - T0 = première administration parmi les trois ; bras = médicament de T0 ;
#'  - exclusion si un AUTRE des trois est administré dans les 24 h suivant T0
#'    (co-initiation, bras ambigu) ou si le premier carbapénème est un autre
#'    que le méropénème avant T0+24 h ;
#'  - un épisode par patient (le premier T0 éligible) ;
#'  - adultes, créatinine de baseline (min 7 j) disponible, pas de RRT avant
#'    T0, pas d'ESRD (codes ICD du séjour index, comme v2/v3).
#'  L'AKI à T0 est CALCULÉE mais pas exclue ici (exclusion = spec de cohorte).
build_antipseudo_initiators_v31 <- function(ab, tbl_paths, creat_dt) {
  st <- ab[abx_class %in% c("ptz", "cefepime", "meropenem", "carbapenem_other") & iv == 1L & !is.na(hadm_id)]
  first <- st[, .(t_first = min(t)), by = .(subject_id, hadm_id, abx_class)]
  w <- dcast(first, subject_id + hadm_id ~ abx_class, value.var = "t_first")
  for (cl in c("ptz", "cefepime", "meropenem", "carbapenem_other")) if (!cl %in% names(w)) w[, (cl) := as.POSIXct(NA)]
  w[, T0 := pmin(ptz, cefepime, meropenem, carbapenem_other, na.rm = TRUE)]
  w[, treatment := fcase(!is.na(ptz) & ptz == T0, "PTZ", !is.na(cefepime) & cefepime == T0, "cefepime",
                         !is.na(meropenem) & meropenem == T0, "meropenem", default = "carbapenem_other")]
  co <- function(x) !is.na(x) & x > w$T0 & x <= w$T0 + 24 * 3600
  w[, coinit := (treatment != "PTZ" & co(ptz)) | (treatment != "cefepime" & co(cefepime)) |
      (treatment != "meropenem" & co(meropenem)) | co(carbapenem_other) |
      (treatment == "PTZ" & !is.na(cefepime) & cefepime == T0) | (treatment == "PTZ" & !is.na(meropenem) & meropenem == T0)]
  fl <- list(); fl$all <- nrow(w)
  # v3.3 (STROBE) : décomptes par bras à chaque étape, sans changer la cohorte
  fla <- list(); arm <- function(step, x) x[, .(step = step, n = .N), by = treatment]
  fla$all <- arm("all", w)
  w <- w[treatment != "carbapenem_other" & !coinit]; fl$single_arm <- nrow(w); fla$single_arm <- arm("single_arm", w)
  w <- w[order(subject_id, T0)][, .SD[1], by = subject_id]; fl$first_episode <- nrow(w); fla$first_episode <- arm("first_episode", w)

  pat <- read_mimic_table(tbl_paths$paths[["patients"]], col_select = c("subject_id", "gender", "anchor_age", "anchor_year", "dod"))
  adm <- read_mimic_table(tbl_paths$paths[["admissions"]],
                          col_select = c("subject_id", "hadm_id", "admittime", "dischtime", "deathtime", "edregtime", "admission_type"))
  pat[, `:=`(subject_id = as.integer(subject_id), dod = as.Date(substr(dod, 1, 10)))]
  adm[, `:=`(subject_id = as.integer(subject_id), hadm_id = as.integer(hadm_id),
             admittime = lubridate::ymd_hms(admittime, quiet = TRUE), dischtime = lubridate::ymd_hms(dischtime, quiet = TRUE),
             deathtime = lubridate::ymd_hms(deathtime, quiet = TRUE), edregtime = lubridate::ymd_hms(edregtime, quiet = TRUE))]
  c0 <- merge(w[, .(subject_id, hadm_id, T0, treatment)], adm, by = c("subject_id", "hadm_id"))
  c0 <- merge(c0, pat, by = "subject_id")
  c0[, age := anchor_age + as.numeric(difftime(admittime, lubridate::make_date(anchor_year, 7, 1), units = "days")) / 365.25]
  c0 <- c0[age >= 18]; fl$adult <- nrow(c0); fla$adult <- arm("adult", c0)
  b <- compute_baseline_creatinine_v3(creat_dt, c0[, .(subject_id, hadm_id, T0, admit_time = admittime, age, gender)])
  c0 <- merge(c0, b, by = c("subject_id", "hadm_id"))
  excl_baseline <- c0[is.na(baseline_creat_min7d)]
  c0 <- c0[!is.na(baseline_creat_min7d)]; fl$baseline <- nrow(c0); fla$baseline <- arm("baseline", c0)
  rrt <- extract_rrt(tbl_paths, c0[, .(subject_id, hadm_id, T0)], followup_days = 7)
  c0 <- merge(c0, rrt, by = c("subject_id", "hadm_id"), all.x = TRUE)
  c0[is.na(rrt_pre_t0), rrt_pre_t0 := 0L]; c0[is.na(rrt_7d), rrt_7d := 0L]
  c0 <- c0[rrt_pre_t0 == 0]; fl$no_rrt <- nrow(c0); fla$no_rrt <- arm("no_rrt", c0)
  icd <- get_icd_codes()
  diag <- read_mimic_table(tbl_paths$paths[["diagnoses_icd"]], col_select = c("subject_id", "hadm_id", "icd_code", "icd_version"))
  diag[, `:=`(subject_id = as.integer(subject_id), hadm_id = as.integer(hadm_id))]
  esrd <- diag[(icd_version == 9 & icd_code %in% icd$esrd$icd9) |
                 (icd_version == 10 & Reduce(`|`, lapply(icd$esrd$icd10, function(p) startsWith(icd_code, p)))),
               .(esrd = 1L), by = .(subject_id, hadm_id)]
  c0 <- merge(c0, esrd, by = c("subject_id", "hadm_id"), all.x = TRUE)
  c0 <- c0[is.na(esrd)][, esrd := NULL]; fl$no_esrd <- nrow(c0); fla$no_esrd <- arm("no_esrd", c0)
  c0[, aki_at_t0 := as.integer(aki_at_t0 == 1L)]
  attr(c0, "flow") <- data.table(step = names(fl), n = unlist(fl))
  attr(c0, "flow_by_arm") <- rbindlist(fla)
  attr(c0, "excluded_no_baseline") <- excl_baseline
  message("[57_cohort_v31] Initiateurs PTZ/céfépime/méropénème : ", nrow(c0), " (",
          paste(names(table(c0$treatment)), table(c0$treatment), sep = "=", collapse = ", "), ")")
  c0
}

#' Covariables du PS v3.1 + outcomes pour la cohorte unifiée. Reproduit les
#' définitions de la cohorte primaire (R/32_, R/50_, R/53_, R/55_), avec les
#' administrations effectives pour les expositions antibiotiques.
add_covariates_outcomes_antipseudo_v31 <- function(c0, ab, meds_v3, ce, lab, creat_dt, tbl_paths) {
  d <- copy(c0)
  key <- d[, .(subject_id, hadm_id, T0, treatment)]
  T0_dt <- d[, .(subject_id, hadm_id, T0)]
  d[, arrival_time := pmin(admittime, fcoalesce(edregtime, admittime))]
  d[, h_arrival_to_T0 := as.numeric(difftime(T0, arrival_time, units = "hours"))]
  d[, time_admission_to_T0_h := pmax(0, h_arrival_to_T0)]   # même définition que master (R/50_)

  # Comorbidités (hospitalisations antérieures), réanimation, sévérité
  icdp <- build_icd_covariates_prior_v3(d[, .(subject_id, hadm_id, admittime)], tbl_paths, get_icd_codes())
  setnames(icdp, c("ckd", "diabetes", "heart_failure", "liver_disease", "cancer", "immunosuppression", "sepsis_icd", "charlson_score"),
           paste0(c("ckd", "diabetes", "heart_failure", "liver_disease", "cancer", "immunosuppression", "sepsis_icd", "charlson_score"), "_priorhosp"))
  d <- merge(d, icdp, by = c("subject_id", "hadm_id"), all.x = TRUE)
  icu <- extract_icustays(tbl_paths)
  icu1 <- icu[, .(icu_intime = min(intime)), by = .(subject_id, hadm_id)]
  d <- merge(d, icu1, by = c("subject_id", "hadm_id"), all.x = TRUE)
  d[, icu_before_t0 := as.integer(!is.na(icu_intime) & icu_intime <= T0)]
  sev <- extract_severity_chartevents(tbl_paths, T0_dt)
  d <- merge(d, sev[, .(subject_id, hadm_id, vasopressor_pre_t0, vent_pre_t0)], by = c("subject_id", "hadm_id"), all.x = TRUE)
  ie <- extract_inputevents(tbl_paths, T0_dt)
  d <- merge(d, ie, by = c("subject_id", "hadm_id"), all.x = TRUE)
  d[, vasopressor_pre_t0 := as.integer(vasopressor_pre_t0 %in% 1L | vasopressor_inputevents_pre %in% 1L)]
  d[, vent_pre_t0 := as.integer(vent_pre_t0 %in% 1L)]

  # Antibiotiques (administrations effectives)
  abm <- merge(ab[iv == 1L], T0_dt, by = c("subject_id", "hadm_id"))
  flag <- function(cond, name) { f <- abm[eval(cond), .(v = 1L), by = .(subject_id, hadm_id)]; setnames(f, "v", name)
    d <<- merge(d, f, by = c("subject_id", "hadm_id"), all.x = TRUE); d[is.na(get(name)), (name) := 0L] }
  flag(quote(abx_class == "vancomycin" & t >= T0 - 24 * 3600 & t <= T0 + 6 * 3600), "vanco_6h")
  flag(quote(abx_class == "aminoglycoside" & t >= T0 - 48 * 3600 & t < T0), "aminoglycoside_pre_t0")
  flag(quote(abx_class %in% c("ptz", "meropenem", "carbapenem_other", "vancomycin", "aminoglycoside") &
               t >= T0 - 48 * 3600 & t < T0), "prior_antibiotic")
  flag(quote(t >= T0 - 48 * 3600 & t < T0), "any_iv_abx_48h_pre")
  # vrai new-user 90 j (toutes hospitalisations, administrations effectives)
  ap <- merge(ab[abx_class %in% c("ptz", "meropenem", "carbapenem_other", "cefepime", "ceftazidime", "aztreonam")],
              d[, .(subject_id, T0)], by = "subject_id", allow.cartesian = TRUE)
  pu <- ap[t < T0 & t >= T0 - 90 * 86400, .(prior_ap_90d = 1L), by = subject_id]
  p7 <- ap[t < T0 & t >= T0 - 7 * 86400, .(n_doses_7d = .N), by = subject_id]
  d <- merge(d, pu, by = "subject_id", all.x = TRUE); d <- merge(d, p7, by = "subject_id", all.x = TRUE)
  d[, true_new_user_90d := as.integer(is.na(prior_ap_90d))]
  d[, acorn_new_user_7d := as.integer(is.na(n_doses_7d) | n_doses_7d <= 1L)]   # ACORN : <= 1 dose dans les 7 j
  d[, c("prior_ap_90d", "n_doses_7d") := NULL]

  d <- merge(d, build_nephrotoxin_covariates_v3(meds_v3, T0_dt), by = c("subject_id", "hadm_id"), all.x = TRUE)
  mic <- build_microbiology_v3(T0_dt, tbl_paths)
  d <- merge(d, mic[, .(subject_id, hadm_id, positive_culture_pre_t0_storetime, esbl_resistant_history_12mo)],
             by = c("subject_id", "hadm_id"), all.x = TRUE)
  d <- merge(d, build_infection_site_v31(tbl_paths, T0_dt), by = c("subject_id", "hadm_id"), all.x = TRUE)
  # suspicion d'infection (ACORN) : culture prélevée dans [T0 - 24 h ; T0 + 24 h]
  sus <- build_infection_site_v31(tbl_paths, d[, .(subject_id, hadm_id, T0 = T0 + 24 * 3600)], window_h = 48)
  d <- merge(d, sus[, .(subject_id, hadm_id, culture_drawn_pm24h = culture_drawn_pre_t0)], by = c("subject_id", "hadm_id"), all.x = TRUE)
  nc <- build_new_covariates_v3(d[, .(subject_id, hadm_id, T0, admittime)], tbl_paths)
  d <- merge(d, nc, by = c("subject_id", "hadm_id"), all.x = TRUE)
  d[, service_group := group_service_v31(service_at_t0)]
  d <- merge(d, read_anchor_year_group(tbl_paths), by = "subject_id", all.x = TRUE)
  d[is.na(anchor_year_group), anchor_year_group := "unknown"]
  d[, lactate_max := NULL]
  pre <- summarise_pre_t0_v31(ce, lab, key)
  d <- merge(d, pre$cov, by = c("subject_id", "hadm_id"), all.x = TRUE)
  d[, weight_kg_v31 := fcoalesce(weight_icu_kg, fifelse(.in_range(weight_kg, "weight"), weight_kg, NA_real_))]
  d[, sepsis_composite_v3 := as.integer(vasopressor_pre_t0 == 1 | vent_pre_t0 == 1 | (!is.na(lactate_max) & lactate_max >= 2))]
  d <- merge(d, build_contrast_proxy_v31(tbl_paths, T0_dt), by = c("subject_id", "hadm_id"), all.x = TRUE)
  d <- merge(d, compute_creat_counts_pre_t0_v31(creat_dt, T0_dt), by = c("subject_id", "hadm_id"), all.x = TRUE)

  # Outcomes
  for (wd in c(14, 30)) d <- merge(d, extract_rrt_window_v3(tbl_paths, d, window_days = wd), by = c("subject_id", "hadm_id"), all.x = TRUE)
  d <- merge(d, compute_kdigo_v31(creat_dt, d, 7, "rrt_7d"), by = c("subject_id", "hadm_id"), all.x = TRUE)
  d <- merge(d, compute_kdigo_v31(creat_dt, d, 14, "rrt_14d"), by = c("subject_id", "hadm_id"), all.x = TRUE)
  d <- merge(d, compute_early_creat_rise_v31(creat_dt, d), by = c("subject_id", "hadm_id"), all.x = TRUE)
  d[, nc_hypoglycemia_d7 := fifelse(!is.na(n_glucose_d7) & n_glucose_d7 > 0 & !(hypoglycemia_pre_t0 %in% 1L), hypoglycemia_d7, NA_integer_)]
  for (days in c(7, 14, 30))
    d[, (paste0("death_", days, "d")) := as.integer((!is.na(deathtime) & deathtime <= T0 + days * 86400) |
                                                      (!is.na(dod) & dod <= as.Date(T0) + days))]
  d[, outcome_primary_v31 := aki_d7_v31]
  # Échelle ordinale ACORN à J14 : 0 pas d'AKI, 1-3 stade KDIGO max, 4 décès
  d[, acorn_ordinal_d14 := fifelse(death_14d == 1L, 4L, aki_stage_max_d14_v31)]
  # MAKE14 : décès, nouvelle RRT, ou dernière créatinine à J14 >= 2 x baseline
  cr <- merge(creat_dt, d[, .(subject_id, hadm_id, T0, baseline_creat_min7d)], by = c("subject_id", "hadm_id"))
  last14 <- cr[charttime > T0 & charttime <= T0 + 14 * 86400][order(charttime)][, .(cr_last14 = creat_mgdl[.N]), by = .(subject_id, hadm_id)]
  d <- merge(d, last14, by = c("subject_id", "hadm_id"), all.x = TRUE)
  d[, make14 := fifelse(death_14d == 1L | rrt_14d == 1L, 1L,
                        fifelse(!is.na(cr_last14), as.integer(cr_last14 >= 2 * baseline_creat_min7d), NA_integer_))]
  setkey(d, hadm_id)
  d
}
