# =============================================================================
# 53_extract.R — Ré-extraction des constantes et des biologies pré-T0
#   (docs/v3.1_plan.md A4, A5, A11, A7b)
# =============================================================================
# Corrige trois défauts d'extraction de R/05_extract_icu_covariates.R (v2,
# repris tel quel en v3) :
#   1. gcs_min = min() sur les itemids 220739/223900/223901, qui sont les
#      COMPOSANTES Yeux (1-4) / Verbal (1-5) / Moteur (1-6) du GCS et non un
#      total -> moyenne observée 1.43. v3.1 : total E+V+M par charttime.
#   2. Température : conversion °F->°C appliquée à TOUTES les valeurs dès que
#      le max dépasse 50, y compris aux mesures déjà en °C (itemid 223762) ->
#      P1 = 1.7 °C. v3.1 : conversion par itemid.
#   3. Aucune borne de plausibilité (SpO2 = 0, PAM = 270, poids = 830 kg).
#      v3.1 : bornes physiologiques appliquées mesure par mesure
#      (PHYSIO_RANGES_V31), valeurs hors bornes comptées puis mises à NA.
#
# Lecture unique (DuckDB, filtre SQL sur itemid + fenêtre temporelle par
# hadm_id) puis mise en cache. Les fenêtres couvrent toutes les initiations
# PTZ / carbapénème / céfépime, pour que les cohortes de la phase B (empirique
# précoce, ACORN, 3 bras) réutilisent le même cache sans relecture.

CHART_ITEMS_V31 <- list(
  hr = 220045L, map = c(220052L, 220181L, 225312L),
  temp_f = 223761L, temp_c = 223762L, spo2 = 220277L,
  gcs_eye = 220739L, gcs_verbal = 223900L, gcs_motor = 223901L,
  weight_admit_kg = 226512L, weight_daily_kg = 224639L, fio2 = 223835L
)

LAB_ITEMS_V31 <- list(
  platelets = 51265L, bilirubin = 50885L, glucose = c(50931L, 50809L),
  lactate = 50813L, bun = 51006L
)

# Bornes physiologiques (inclusives) appliquées AVANT résumé et imputation.
PHYSIO_RANGES_V31 <- list(
  hr = c(20, 250), map = c(20, 200), temp = c(30, 43), spo2 = c(50, 100),
  gcs = c(3, 15), weight = c(25, 300), bmi = c(10, 80), platelets = c(1, 2000),
  bilirubin = c(0.1, 50), glucose = c(10, 1500), lactate = c(0.1, 30),
  bun = c(1, 300), creatinine = c(0.1, 25)
)

.in_range <- function(x, var) {
  r <- PHYSIO_RANGES_V31[[var]]
  !is.na(x) & x >= r[1] & x <= r[2]
}

#' Fenêtres d'extraction par hadm_id, couvrant toutes les initiations
#' anti-pyocyaniques possibles (T0 de n'importe quelle cohorte v3.1).
#' @param pre_h  heures avant la première initiation
#' @param post_h heures après la DERNIÈRE « première initiation » d'une classe
build_extraction_windows_v31 <- function(meds_v3, pre_h = 72, post_h = 48) {
  x <- meds_v3[drug_class %in% c("ptz", "carbapenem", "cefepime") & !is.na(starttime),
               .(first_start = min(starttime)), by = .(hadm_id, drug_class)]
  w <- x[, .(t_lo = min(first_start) - pre_h * 3600,
             t_hi = max(first_start) + post_h * 3600), by = hadm_id]
  w[!is.na(hadm_id)]
}

.duckdb_window_query <- function(path, item_ids, windows, select_cols, time_col = "charttime") {
  con <- duckdb::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(duckdb::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  win_df <- as.data.frame(windows[, .(hadm_id = as.integer(hadm_id), t_lo, t_hi)])
  duckdb::duckdb_register(con, "win", win_df)
  sql <- sprintf(
    "SELECT %s FROM read_csv_auto('%s') c
     JOIN win w ON c.hadm_id = w.hadm_id
     WHERE c.itemid IN (%s) AND c.%s >= w.t_lo AND c.%s <= w.t_hi",
    paste0("c.", select_cols, collapse = ", "), path,
    paste(item_ids, collapse = ", "), time_col, time_col)
  data.table::as.data.table(DBI::dbGetQuery(con, sql))
}

#' Lecture chartevents (constantes, composantes GCS, poids) — cache
extract_chartevents_long_v31 <- function(tbl_paths, windows,
                                          cache = "output/extraction/derived_data/chartevents_long_v31.rds") {
  if (file.exists(cache)) {
    message("[53_extract_v31] chartevents : cache ", cache)
    return(readRDS(cache))
  }
  t <- Sys.time()
  message("[53_extract_v31] Lecture chartevents (DuckDB, ", nrow(windows), " séjours)...")
  ce <- .duckdb_window_query(tbl_paths$paths[["chartevents"]], unlist(CHART_ITEMS_V31), windows,
                             c("subject_id", "hadm_id", "itemid", "charttime", "value", "valuenum"))
  ce[, `:=`(subject_id = as.integer(subject_id), hadm_id = as.integer(hadm_id),
            itemid = as.integer(itemid))]
  saveRDS(ce, cache)
  message("[53_extract_v31] chartevents : ", nrow(ce), " lignes (",
          round(difftime(Sys.time(), t, units = "mins"), 1), " min)")
  ce
}

#' Lecture labevents (plaquettes, bilirubine, glucose, lactate, urée) — cache
extract_labs_long_v31 <- function(tbl_paths, windows,
                                   cache = "output/extraction/derived_data/labs_long_v31.rds") {
  if (file.exists(cache)) {
    message("[53_extract_v31] labevents : cache ", cache)
    return(readRDS(cache))
  }
  t <- Sys.time()
  # fenêtre post-T0 étendue à 8 jours pour l'hypoglycémie J1-J7 (A7b)
  win_lab <- copy(windows)[, t_hi := t_hi + 6 * 24 * 3600]
  message("[53_extract_v31] Lecture labevents (DuckDB, ", nrow(win_lab), " séjours)...")
  lab <- .duckdb_window_query(tbl_paths$paths[["labevents"]], unlist(LAB_ITEMS_V31), win_lab,
                              c("subject_id", "hadm_id", "itemid", "charttime", "valuenum"))
  lab[, `:=`(subject_id = as.integer(subject_id), hadm_id = as.integer(hadm_id),
             itemid = as.integer(itemid))]
  lab <- lab[!is.na(valuenum)]
  saveRDS(lab, cache)
  message("[53_extract_v31] labevents : ", nrow(lab), " lignes (",
          round(difftime(Sys.time(), t, units = "mins"), 1), " min)")
  lab
}

# -----------------------------------------------------------------------
# GCS total (A4) — convention MIMIC-Code (mimic-iv/concepts/measurement/gcs.sql)
# -----------------------------------------------------------------------
#' Les 3 composantes sont regroupées par charttime ; une composante absente
#' est reprise de la dernière valeur de la même composante dans les 6 h
#' précédentes. Verbal « No Response-ETT » (patient intubé, verbal non
#' évaluable) : GCS total fixé à 15, comme MIMIC-Code (convention SOFA
#' « non évaluable sous sédation = normal ») ; l'intubation est captée
#' séparément par vent_pre_t0. Total retenu seulement si E, V et M sont
#' tous disponibles (après report 6 h) ou si ETT.
compute_gcs_total_v31 <- function(ce) {
  g <- ce[itemid %in% c(CHART_ITEMS_V31$gcs_eye, CHART_ITEMS_V31$gcs_verbal,
                        CHART_ITEMS_V31$gcs_motor)]
  g[, comp := fcase(itemid == CHART_ITEMS_V31$gcs_eye, "E",
                    itemid == CHART_ITEMS_V31$gcs_verbal, "V",
                    itemid == CHART_ITEMS_V31$gcs_motor, "M")]
  g[, ett := as.integer(comp == "V" & grepl("ETT", value, ignore.case = TRUE))]
  g <- g[, .(val = suppressWarnings(max(valuenum, na.rm = TRUE)), ett = max(ett)),
         by = .(subject_id, hadm_id, charttime, comp)]
  g[is.infinite(val), val := NA_real_]
  w <- dcast(g, subject_id + hadm_id + charttime ~ comp, value.var = c("val", "ett"), fill = NA)
  for (cc in c("val_E", "val_V", "val_M", "ett_V")) if (!cc %in% names(w)) w[, (cc) := NA_real_]
  setorder(w, subject_id, hadm_id, charttime)

  # report 6 h de la dernière valeur observée, par composante
  carry6h <- function(v, tt) {
    out <- v; last_v <- NA_real_; last_t <- NA_real_
    tn <- as.numeric(tt)
    for (i in seq_along(v)) {
      if (!is.na(v[i])) { last_v <- v[i]; last_t <- tn[i] }
      else if (!is.na(last_t) && tn[i] - last_t <= 6 * 3600) out[i] <- last_v
    }
    out
  }
  w[, `:=`(E = carry6h(val_E, charttime), V = carry6h(val_V, charttime),
           M = carry6h(val_M, charttime), ETT = carry6h(ett_V, charttime)),
    by = .(subject_id, hadm_id)]
  w[, gcs_total := fifelse(!is.na(ETT) & ETT == 1, 15,
                           fifelse(!is.na(E) & !is.na(V) & !is.na(M), E + V + M, NA_real_))]
  w[, gcs_ett := as.integer(!is.na(ETT) & ETT == 1)]
  w[!is.na(gcs_total), .(subject_id, hadm_id, charttime, gcs_total, gcs_ett)]
}

# -----------------------------------------------------------------------
# Résumé pré-T0 des constantes et biologies, avec bornes physiologiques
# -----------------------------------------------------------------------
#' @param T0_dt data.table(subject_id, hadm_id, T0, treatment)
#' @return list(cov = data.table une ligne par épisode, range_report = data.table)
summarise_pre_t0_v31 <- function(ce, lab, T0_dt, window_h = 24) {
  key <- T0_dt[, .(subject_id, hadm_id, T0, treatment)]
  range_rows <- list()
  note_range <- function(dt, var, raw_col = "x") {
    # dt contient treatment et la valeur brute ; compte les mesures hors bornes
    rr <- dt[, .(n_meas = .N, n_out_of_range = sum(!.in_range(get(raw_col), var)),
                 min_raw = suppressWarnings(min(get(raw_col), na.rm = TRUE)),
                 max_raw = suppressWarnings(max(get(raw_col), na.rm = TRUE))), by = treatment]
    rr[, variable := var]
    range_rows[[length(range_rows) + 1]] <<- rr
  }

  # --- constantes (chartevents) dans [T0 - window_h ; T0]
  cv <- merge(ce[!is.na(valuenum)], key, by = c("subject_id", "hadm_id"))
  cv <- cv[charttime <= T0 & charttime >= T0 - window_h * 3600]
  vit <- function(ids, var, transform = identity) {
    d <- cv[itemid %in% ids, .(subject_id, hadm_id, treatment, x = transform(valuenum), itemid)]
    note_range(d, var)
    d[.in_range(x, var)]
  }
  hr   <- vit(CHART_ITEMS_V31$hr, "hr")
  mp   <- vit(CHART_ITEMS_V31$map, "map")
  tf   <- cv[itemid == CHART_ITEMS_V31$temp_f, .(subject_id, hadm_id, treatment, x = (valuenum - 32) * 5 / 9)]
  tc   <- cv[itemid == CHART_ITEMS_V31$temp_c, .(subject_id, hadm_id, treatment, x = valuenum)]
  tmp  <- rbind(tf, tc); note_range(tmp, "temp"); tmp <- tmp[.in_range(x, "temp")]
  sp   <- vit(CHART_ITEMS_V31$spo2, "spo2")

  gcs <- compute_gcs_total_v31(ce)
  gcs <- merge(gcs, key, by = c("subject_id", "hadm_id"))
  gcs <- gcs[charttime <= T0 & charttime >= T0 - window_h * 3600]
  gcs_d <- gcs[, .(subject_id, hadm_id, treatment, x = gcs_total, gcs_ett)]
  note_range(gcs_d, "gcs"); gcs_d <- gcs_d[.in_range(x, "gcs")]

  out <- key[, .(subject_id, hadm_id)]
  add <- function(out, d, fun, name) {
    if (nrow(d) == 0) { out[, (name) := NA_real_]; return(out) }
    a <- d[, .(v = fun(x)), by = .(subject_id, hadm_id)]
    setnames(a, "v", name)
    merge(out, a, by = c("subject_id", "hadm_id"), all.x = TRUE)
  }
  out <- add(out, hr, mean, "hr_mean")
  out <- add(out, mp, mean, "map_mean")
  out <- add(out, tmp, max, "temp_max")
  out <- add(out, sp, min, "spo2_min")
  out <- add(out, gcs_d, min, "gcs_min")
  ett <- gcs_d[, .(gcs_ett_any = as.integer(any(gcs_ett == 1))), by = .(subject_id, hadm_id)]
  out <- merge(out, ett, by = c("subject_id", "hadm_id"), all.x = TRUE)
  out[is.na(gcs_ett_any), gcs_ett_any := 0L]

  # --- poids d'admission en réa (chartevents, tout le séjour avant T0)
  wt <- merge(ce[itemid %in% c(CHART_ITEMS_V31$weight_admit_kg, CHART_ITEMS_V31$weight_daily_kg) &
                   !is.na(valuenum)], key, by = c("subject_id", "hadm_id"))
  wt <- wt[charttime <= T0, .(subject_id, hadm_id, treatment, x = valuenum, charttime)]
  note_range(wt, "weight"); wt <- wt[.in_range(x, "weight")]
  wt_last <- wt[order(subject_id, hadm_id, -charttime)][, .SD[1], by = .(subject_id, hadm_id)]
  out <- merge(out, wt_last[, .(subject_id, hadm_id, weight_icu_kg = x)],
               by = c("subject_id", "hadm_id"), all.x = TRUE)

  # --- biologies (labevents)
  lb <- merge(lab, key, by = c("subject_id", "hadm_id"))
  lb_pre  <- lb[charttime <= T0 & charttime >= T0 - window_h * 3600]
  labv <- function(ids, var, d = lb_pre) {
    x <- d[itemid %in% ids, .(subject_id, hadm_id, treatment, x = valuenum)]
    note_range(x, var); x[.in_range(x, var)]
  }
  out <- add(out, labv(LAB_ITEMS_V31$platelets, "platelets"), min, "platelets_min")
  out <- add(out, labv(LAB_ITEMS_V31$bilirubin, "bilirubin"), max, "bilirubin_max")
  out <- add(out, labv(LAB_ITEMS_V31$lactate, "lactate"), max, "lactate_max")
  out <- add(out, labv(LAB_ITEMS_V31$bun, "bun"), mean, "bun_baseline")

  # --- glucose : hypoglycémie pré-T0 (48 h) et incidente (T0 ; T0+7 j] (A7b)
  glu <- lb[itemid %in% LAB_ITEMS_V31$glucose, .(subject_id, hadm_id, treatment, T0, charttime, x = valuenum)]
  note_range(glu, "glucose"); glu <- glu[.in_range(x, "glucose")]
  glu_pre  <- glu[charttime <= T0 & charttime >= T0 - 48 * 3600,
                  .(hypoglycemia_pre_t0 = as.integer(any(x < 70))), by = .(subject_id, hadm_id)]
  glu_post <- glu[charttime > T0 & charttime <= T0 + 7 * 24 * 3600,
                  .(hypoglycemia_d7 = as.integer(any(x < 70)), n_glucose_d7 = .N),
                  by = .(subject_id, hadm_id)]
  out <- merge(out, glu_pre, by = c("subject_id", "hadm_id"), all.x = TRUE)
  out <- merge(out, glu_post, by = c("subject_id", "hadm_id"), all.x = TRUE)

  range_report <- rbindlist(range_rows, fill = TRUE)
  setcolorder(range_report, c("variable", "treatment"))
  range_report[is.infinite(min_raw), min_raw := NA_real_]
  range_report[is.infinite(max_raw), max_raw := NA_real_]
  list(cov = out, range_report = range_report)
}

# -----------------------------------------------------------------------
# Site infectieux présumé (A11) — cultures PRÉLEVÉES avant T0
# -----------------------------------------------------------------------
#' Le prélèvement (charttime / chartdate du spécimen) est une information
#' connue du prescripteur à T0 ; seul le RÉSULTAT (storetime) ne l'est pas.
#' v3 exigeait un résultat positif disponible avant T0 -> 98 % de « aucune ».
#' Indicateurs non exclusifs par site, fenêtre [T0 - 24 h ; T0].
build_infection_site_v31 <- function(tbl_paths, T0_dt, window_h = 24) {
  micro <- read_mimic_table(tbl_paths$paths[["microbiologyevents"]],
                            col_select = c("subject_id", "hadm_id", "chartdate",
                                           "charttime", "spec_type_desc"))
  micro[, `:=`(subject_id = as.integer(subject_id), hadm_id = as.integer(hadm_id))]
  micro[, t := lubridate::ymd_hms(charttime, quiet = TRUE)]
  # chartdate seule (heure inconnue) : on prend la fin de journée, i.e. on
  # ne la compte comme pré-T0 que si le jour entier précède T0
  micro[is.na(t), t := lubridate::ymd(substr(chartdate, 1, 10), quiet = TRUE) + 86399]
  m <- merge(micro, T0_dt[, .(subject_id, hadm_id, T0)], by = c("subject_id", "hadm_id"))
  m <- m[!is.na(t) & t <= T0 & t >= T0 - window_h * 3600]
  s <- tolower(m$spec_type_desc)
  m[, site := fcase(
    grepl("blood culture", s), "blood",
    grepl("urine", s), "urine",
    grepl("sputum|bronch|bal|tracheal|respiratory|pleural|lung", s), "respiratory",
    grepl("peritoneal|abdom|bile|biliary|ascites|stool|c\\. diff|fecal", s), "abdominal",
    grepl("swab|wound|tissue|abscess|skin|bone|joint|fluid received in blood culture bottles", s), "skin_soft_tissue_bone",
    grepl("mrsa screen|screen|serology|blood \\(", s), NA_character_,   # dépistages, sérologies : pas un site
    default = "other")]
  m <- m[!is.na(site)]
  res <- T0_dt[, .(subject_id, hadm_id)]
  for (st in c("blood", "urine", "respiratory", "abdominal", "skin_soft_tissue_bone", "other")) {
    col <- paste0("site_", st)
    f <- unique(m[site == st, .(subject_id, hadm_id)])[, (col) := 1L]
    res <- merge(res, f, by = c("subject_id", "hadm_id"), all.x = TRUE)
    res[is.na(get(col)), (col) := 0L]
  }
  res[, culture_drawn_pre_t0 := as.integer(rowSums(.SD) > 0), .SDcols = patterns("^site_")]
  message("[53_extract_v31] Culture prélevée dans les ", window_h, " h pré-T0 : ",
          round(mean(res$culture_drawn_pre_t0) * 100, 1), "%")
  res
}

# -----------------------------------------------------------------------
# Service à T0 regroupé (A3) — 40 careunit bruts -> 7 niveaux cliniques
# -----------------------------------------------------------------------
group_service_v31 <- function(careunit) {
  s <- tolower(ifelse(is.na(careunit), "", careunit))
  # règles ordonnées : la première qui correspond l'emporte
  fcase(
    grepl("emergency|^observation", s), "ED_observation",
    grepl("\\(micu\\)|\\(ccu\\)|micu/sicu|^intensive care unit", s), "ICU_medical",
    grepl("\\(sicu\\)|\\(cvicu\\)|\\(tsicu\\)|neuro sicu", s), "ICU_surgical_cardiac_neuro",
    grepl("hematology|oncology|transplant", s), "ward_heme_onc_transplant",
    grepl("^surg|cardiac surgery|cardiology surgery|trauma|pacu|vascular", s), "ward_surgical",
    grepl("medicine|^cardiology|neuro|^med/surg", s), "ward_medical",
    default = "other")
}

# -----------------------------------------------------------------------
# Proxy « produit de contraste iodé » (A11, SENSIBILITÉ UNIQUEMENT)
# -----------------------------------------------------------------------
#' Non mesurable de façon fiable dans MIMIC-IV : eMAR ne trace quasiment pas
#' le contraste radiologique (3 399 lignes dans toute la base) et les
#' procédures ICD sont datées au jour près, codées à la sortie et
#' sous-déclarées. Proxy : procédure avec contraste iodé datée dans les 3
#' jours civils précédant T0 (inclus), même hospitalisation.
#'   ICD-9  : 88.4x (artériographie), 88.5x (angiocardiographie), 37.21-37.23
#'            (cathétérisme cardiaque)
#'   ICD-10 : section B imagerie, type 0/1/2 (radio/fluoro/TDM), 5e caractère
#'            de contraste 0/1/Y ; 4A023N6-8 (cathétérisme) ; 0270-0273 (ATL
#'            coronaire)
build_contrast_proxy_v31 <- function(tbl_paths, T0_dt) {
  pr <- read_mimic_table(tbl_paths$paths[["procedures_icd"]],
                         col_select = c("subject_id", "hadm_id", "chartdate", "icd_code", "icd_version"))
  pr[, `:=`(subject_id = as.integer(subject_id), hadm_id = as.integer(hadm_id),
            chartdate = as.Date(substr(chartdate, 1, 10)))]
  code <- pr$icd_code
  is_contrast <- (pr$icd_version == 9 & (grepl("^884|^885", code) | code %in% c("3721", "3722", "3723"))) |
    (pr$icd_version == 10 & ((substr(code, 1, 1) == "B" & substr(code, 3, 3) %in% c("0", "1", "2") &
                                substr(code, 5, 5) %in% c("0", "1", "Y")) |
                               grepl("^4A023N[678]", code) | grepl("^027[0-3]", code)))
  pr <- pr[is_contrast]
  m <- merge(pr, T0_dt[, .(subject_id, hadm_id, T0)], by = c("subject_id", "hadm_id"))
  m <- m[chartdate <= as.Date(T0) & chartdate >= as.Date(T0) - 3]
  res <- T0_dt[, .(subject_id, hadm_id)]
  f <- unique(m[, .(subject_id, hadm_id)])[, contrast_proxy_pre_t0 := 1L]
  res <- merge(res, f, by = c("subject_id", "hadm_id"), all.x = TRUE)
  res[is.na(contrast_proxy_pre_t0), contrast_proxy_pre_t0 := 0L]
  message("[53_extract_v31] Proxy contraste (procédure ICD datée J-3..J0) : ",
          round(mean(res$contrast_proxy_pre_t0) * 100, 2), "%")
  res
}

# -----------------------------------------------------------------------
# Nombre de créatinines pré-T0 (A6)
# -----------------------------------------------------------------------
compute_creat_counts_pre_t0_v31 <- function(creat_dt, T0_dt) {
  m <- merge(creat_dt, T0_dt[, .(subject_id, hadm_id, T0)], by = c("subject_id", "hadm_id"))
  m[, h := as.numeric(difftime(T0, charttime, units = "hours"))]
  res <- m[, .(n_creat_pre_7d = sum(h >= 0 & h <= 7 * 24),
               n_creat_pre_48h = sum(h >= 0 & h <= 48)), by = .(subject_id, hadm_id)]
  res <- merge(T0_dt[, .(subject_id, hadm_id)], res, by = c("subject_id", "hadm_id"), all.x = TRUE)
  res[is.na(n_creat_pre_7d), n_creat_pre_7d := 0L]
  res[is.na(n_creat_pre_48h), n_creat_pre_48h := 0L]
  res
}
