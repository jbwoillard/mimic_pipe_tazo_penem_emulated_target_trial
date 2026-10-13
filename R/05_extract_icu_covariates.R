# =============================================================================
# 05_extract_icu_covariates.R — Covariables ICU (chartevents, inputevents, etc.)
# =============================================================================

#' Extraire les données de séjour ICU
extract_icustays <- function(tbl_paths) {
  paths <- tbl_paths$paths
  if (is.na(paths[["icustays"]])) {
    message("[05_icu] icustays non disponible")
    return(NULL)
  }

  icu <- read_mimic_table(paths[["icustays"]],
                           col_select = c("subject_id", "hadm_id", "stay_id",
                                          "intime", "outtime", "los",
                                          "first_careunit", "last_careunit"))
  if (is.null(icu)) return(NULL)

  icu[, intime  := lubridate::ymd_hms(intime,  quiet = TRUE)]
  icu[, outtime := lubridate::ymd_hms(outtime, quiet = TRUE)]
  icu[, subject_id := as.integer(subject_id)]
  icu[, hadm_id    := as.integer(hadm_id)]
  icu[, stay_id    := as.integer(stay_id)]

  message("[05_icu] Séjours ICU : ", nrow(icu))
  icu
}

#' Extraire indicateurs de sévérité depuis chartevents (SOFA-proxy)
#' Pour chaque patient-épisode, récupérer les valeurs dans la fenêtre [T0-24h; T0]
extract_severity_chartevents <- function(tbl_paths, T0_dt) {
  paths <- tbl_paths$paths
  if (is.na(paths[["chartevents"]])) {
    message("[05_icu] chartevents non disponible — sévérité depuis ICU limitée")
    return(NULL)
  }

  item_ids <- get_item_ids()
  all_ids  <- unique(c(
    item_ids$hr_itemids,
    item_ids$map_itemids,
    item_ids$temp_itemids,
    item_ids$spo2_itemids,
    item_ids$gcs_itemids,
    item_ids$fio2_itemids,
    item_ids$vasopressor_itemids,
    item_ids$vent_itemids
  ))

  # Filtrage DuckDB SQL pour éviter le chargement de 3.5 GB en mémoire
  # On filtre directement sur item_id et les patients de la cohorte
  subject_ids <- unique(T0_dt$subject_id)

  message("[05_icu] Lecture chartevents via DuckDB (filtrage SQL, ",
          length(subject_ids), " patients, ", length(all_ids), " item_ids) ...")

  ce <- tryCatch({
    if (PKG_AVAILABLE("duckdb")) {
      con <- duckdb::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
      on.exit(duckdb::dbDisconnect(con, shutdown = TRUE), add = TRUE)

      # Construire le filtre SQL
      ids_str    <- paste(all_ids,     collapse = ", ")
      subj_str   <- paste(subject_ids, collapse = ", ")

      sql <- sprintf(
        "SELECT subject_id, hadm_id, itemid, charttime, valuenum, warning
         FROM read_csv_auto('%s')
         WHERE itemid IN (%s)
           AND subject_id IN (%s)
           AND valuenum IS NOT NULL",
        paths[["chartevents"]], ids_str, subj_str
      )
      res <- DBI::dbGetQuery(con, sql)
      data.table::as.data.table(res)
    } else {
      message("[05_icu] DuckDB non disponible — lecture complète chartevents")
      read_mimic_table(paths[["chartevents"]],
                        col_select = c("subject_id", "hadm_id",
                                       "itemid", "charttime", "valuenum", "warning"),
                        use_duckdb = FALSE)
    }
  }, error = function(e) {
    message("[05_icu] Erreur lecture chartevents : ", e$message)
    NULL
  })

  if (is.null(ce) || nrow(ce) == 0) {
    message("[05_icu] chartevents non disponible — sévérité ignorée")
    return(NULL)
  }

  ce[, charttime := lubridate::ymd_hms(charttime, quiet = TRUE)]
  ce[, subject_id := as.integer(subject_id)]
  ce[, hadm_id    := as.integer(hadm_id)]

  # Joindre T0 via merge (évite les problèmes de référence i. hors jointure)
  ce_t0 <- merge(ce, T0_dt[, .(subject_id, hadm_id, T0)],
                  by = c("subject_id", "hadm_id"), all.x = FALSE)

  # Garder les valeurs dans les 24h avant T0
  ce_pre <- ce_t0[!is.na(T0) & charttime <= T0 & charttime >= T0 - 24 * 3600]

  if (nrow(ce_pre) == 0) return(NULL)

  # Calculer statistiques par patient-épisode et variable
  severity <- ce_pre[, .(
    hr_mean  = mean(valuenum[itemid %in% item_ids$hr_itemids], na.rm = TRUE),
    hr_max   = max(valuenum[itemid %in% item_ids$hr_itemids],  na.rm = TRUE),
    map_mean = mean(valuenum[itemid %in% item_ids$map_itemids], na.rm = TRUE),
    map_min  = min(valuenum[itemid %in% item_ids$map_itemids],  na.rm = TRUE),
    temp_max = max(valuenum[itemid %in% item_ids$temp_itemids], na.rm = TRUE),
    temp_min = min(valuenum[itemid %in% item_ids$temp_itemids], na.rm = TRUE),
    spo2_min = min(valuenum[itemid %in% item_ids$spo2_itemids], na.rm = TRUE),
    gcs_min  = min(valuenum[itemid %in% item_ids$gcs_itemids],  na.rm = TRUE),
    fio2_max = max(valuenum[itemid %in% item_ids$fio2_itemids], na.rm = TRUE),
    vasopressor_pre_t0 = as.integer(any(itemid %in%
                                          item_ids$vasopressor_itemids)),
    vent_pre_t0 = as.integer(any(itemid %in% item_ids$vent_itemids))
  ), by = .(subject_id, hadm_id)]

  # Remplacer Inf / -Inf par NA
  for (col in names(severity)) {
    if (is.numeric(severity[[col]])) {
      severity[is.infinite(get(col)), (col) := NA_real_]
    }
  }

  # Normaliser température (°F -> °C si max > 50)
  if (!is.na(max(severity$temp_max, na.rm = TRUE)) &&
      max(severity$temp_max, na.rm = TRUE) > 50) {
    severity[, temp_max := (temp_max - 32) * 5 / 9]
    severity[, temp_min := (temp_min - 32) * 5 / 9]
  }

  message("[05_icu] Variables de sévérité extraites pour ",
          nrow(severity), " épisodes")
  severity
}

#' Extraire vasopresseurs depuis inputevents
extract_inputevents <- function(tbl_paths, T0_dt) {
  paths <- tbl_paths$paths
  if (is.na(paths[["inputevents"]])) {
    message("[05_icu] inputevents non disponible")
    return(NULL)
  }

  item_ids <- get_item_ids()

  subject_ids <- unique(T0_dt$subject_id)
  vaso_ids    <- item_ids$vasopressor_input_itemids

  ie <- tryCatch({
    if (PKG_AVAILABLE("duckdb")) {
      con <- duckdb::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
      on.exit(duckdb::dbDisconnect(con, shutdown = TRUE), add = TRUE)
      sql <- sprintf(
        "SELECT subject_id, hadm_id, itemid, starttime, amount, rate
         FROM read_csv_auto('%s')
         WHERE itemid IN (%s) AND subject_id IN (%s)",
        paths[["inputevents"]],
        paste(vaso_ids, collapse = ", "),
        paste(subject_ids, collapse = ", ")
      )
      data.table::as.data.table(DBI::dbGetQuery(con, sql))
    } else {
      tmp <- read_mimic_table(paths[["inputevents"]],
                               col_select = c("subject_id", "hadm_id", "stay_id",
                                              "itemid", "starttime", "amount",
                                              "amountuom", "rate"),
                               use_duckdb = FALSE)
      if (!is.null(tmp)) tmp[itemid %in% vaso_ids] else NULL
    }
  }, error = function(e) { message("[05_icu] Erreur inputevents : ", e$message); NULL })

  if (is.null(ie)) return(NULL)
  ie <- ie[itemid %in% vaso_ids]
  ie[, starttime := lubridate::ymd_hms(starttime, quiet = TRUE)]
  ie[, subject_id := as.integer(subject_id)]
  ie[, hadm_id    := as.integer(hadm_id)]

  # Joindre T0 via merge
  ie_t0 <- merge(ie, T0_dt[, .(subject_id, hadm_id, T0)],
                  by = c("subject_id", "hadm_id"), all.x = FALSE)

  vaso_pre <- ie_t0[!is.na(starttime) & !is.na(T0) & starttime <= T0,
                     .(vasopressor_inputevents_pre = 1L),
                     by = .(subject_id, hadm_id)]

  message("[05_icu] Vasopresseurs inputevents : ",
          nrow(vaso_pre), " épisodes avec vasopresseurs avant T0")
  vaso_pre
}

#' Extraire diurèse depuis outputevents
extract_urine_output <- function(tbl_paths, T0_dt,
                                  window_h_before = 24,
                                  window_h_after  = 7 * 24) {
  paths <- tbl_paths$paths
  if (is.na(paths[["outputevents"]])) {
    message("[05_icu] outputevents non disponible")
    return(NULL)
  }

  item_ids <- get_item_ids()

  subject_ids <- unique(T0_dt$subject_id)
  uo_ids      <- item_ids$urine_output_itemids

  oe <- tryCatch({
    if (PKG_AVAILABLE("duckdb")) {
      con <- duckdb::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
      on.exit(duckdb::dbDisconnect(con, shutdown = TRUE), add = TRUE)
      sql <- sprintf(
        "SELECT subject_id, hadm_id, itemid, charttime, value
         FROM read_csv_auto('%s')
         WHERE itemid IN (%s) AND subject_id IN (%s) AND value IS NOT NULL",
        paths[["outputevents"]],
        paste(uo_ids, collapse = ", "),
        paste(subject_ids, collapse = ", ")
      )
      data.table::as.data.table(DBI::dbGetQuery(con, sql))
    } else {
      tmp <- read_mimic_table(paths[["outputevents"]],
                               col_select = c("subject_id", "hadm_id", "stay_id",
                                              "itemid", "charttime", "value"),
                               use_duckdb = FALSE)
      if (!is.null(tmp)) tmp[itemid %in% uo_ids & !is.na(value)] else NULL
    }
  }, error = function(e) { message("[05_icu] Erreur outputevents : ", e$message); NULL })

  if (is.null(oe)) return(NULL)
  oe <- oe[itemid %in% uo_ids & !is.na(value)]
  oe[, charttime := lubridate::ymd_hms(charttime, quiet = TRUE)]
  oe[, subject_id := as.integer(subject_id)]
  oe[, hadm_id    := as.integer(hadm_id)]
  oe[, value      := as.numeric(value)]

  # Joindre T0 via merge
  oe_t0 <- merge(oe, T0_dt[, .(subject_id, hadm_id, T0)],
                  by = c("subject_id", "hadm_id"), all.x = FALSE)

  # Diurèse dans les 24h avant T0
  uo_pre <- oe_t0[!is.na(charttime) & !is.na(T0) &
                    charttime >= T0 - window_h_before * 3600 &
                    charttime <= T0,
                   .(uo_pre_24h_ml = sum(value, na.rm = TRUE)),
                   by = .(subject_id, hadm_id)]

  # Diurèse dans les 7j après T0 (outcome secondaire)
  uo_post <- oe_t0[!is.na(charttime) & !is.na(T0) &
                     charttime > T0 &
                     charttime <= T0 + window_h_after * 3600,
                    .(uo_post_7d_ml = sum(value, na.rm = TRUE)),
                    by = .(subject_id, hadm_id)]

  result <- merge(uo_pre, uo_post, by = c("subject_id", "hadm_id"), all = TRUE)
  message("[05_icu] Diurèse extraite pour ", nrow(result), " épisodes")
  result
}

#' Extraire informations sur la RRT (dialyse) depuis procedureevents
extract_rrt <- function(tbl_paths, T0_dt, followup_days = 7) {
  paths <- tbl_paths$paths

  rrt_records <- NULL

  # procedureevents ICU
  if (!is.na(paths[["procedureevents"]])) {
    pe <- read_mimic_table(paths[["procedureevents"]],
                            col_select = c("subject_id", "hadm_id", "stay_id",
                                           "itemid", "starttime", "value",
                                           "location"))
    if (!is.null(pe)) {
      rrt_item_ids <- c(
        225441L, 225802L, 225803L, 225805L, 225809L,
        225955L, 226118L, 227290L
      )

      rrt_pe <- pe[itemid %in% rrt_item_ids]
      if (nrow(rrt_pe) > 0) {
        rrt_pe[, starttime := lubridate::ymd_hms(starttime, quiet = TRUE)]
        rrt_pe[, subject_id := as.integer(subject_id)]
        rrt_pe[, hadm_id    := as.integer(hadm_id)]
        rrt_records <- rrt_pe[, .(subject_id, hadm_id, starttime,
                                   rrt_type = "ICU_procedure")]
      }
    }
  }

  if (is.null(rrt_records)) {
    message("[05_icu] Aucune donnée RRT trouvée via procedureevents")
    return(T0_dt[, .(subject_id, hadm_id, rrt_pre_t0 = 0L, rrt_7d = 0L)])
  }

  # Joindre T0 via merge
  rrt_t0 <- merge(rrt_records, T0_dt[, .(subject_id, hadm_id, T0)],
                   by = c("subject_id", "hadm_id"), all.x = FALSE)

  rrt_pre <- rrt_t0[!is.na(starttime) & !is.na(T0) & starttime < T0,
                     .(rrt_pre_t0 = 1L),
                     by = .(subject_id, hadm_id)]

  rrt_post <- rrt_t0[!is.na(starttime) & !is.na(T0) &
                       starttime >= T0 &
                       starttime <= T0 + followup_days * 24 * 3600,
                      .(rrt_7d = 1L),
                      by = .(subject_id, hadm_id)]

  result <- T0_dt[, .(subject_id, hadm_id)]
  result <- merge(result, rrt_pre,  by = c("subject_id", "hadm_id"), all.x = TRUE)
  result <- merge(result, rrt_post, by = c("subject_id", "hadm_id"), all.x = TRUE)
  result[is.na(rrt_pre_t0), rrt_pre_t0 := 0L]
  result[is.na(rrt_7d),     rrt_7d     := 0L]

  message("[05_icu] RRT avant T0 : ", sum(result$rrt_pre_t0),
          " | RRT 7j post-T0 : ", sum(result$rrt_7d))
  result
}
