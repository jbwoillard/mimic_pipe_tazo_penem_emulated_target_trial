# =============================================================================
# 01_paths_and_io.R — Détection automatique des fichiers MIMIC-IV et I/O
# =============================================================================

# Charger la config YAML
load_config <- function(config_path = "config/config.yml") {
  if (!file.exists(config_path)) stop("config.yml introuvable : ", config_path)
  yaml::read_yaml(config_path)
}

# -----------------------------------------------------------------------------
# Détection récursive des tables MIMIC-IV
# -----------------------------------------------------------------------------
detect_mimic_tables <- function(mimic_root) {
  if (!dir.exists(mimic_root)) {
    stop("[01_paths] Répertoire MIMIC introuvable : ", mimic_root)
  }

  # Extensions supportées
  exts <- c("csv", "csv.gz", "parquet", "feather")
  pat  <- paste0("\\.(", paste(exts, collapse = "|"), ")$")

  all_files <- list.files(mimic_root, pattern = pat, recursive = TRUE,
                           full.names = TRUE)

  # Tables MIMIC connues
  hosp_tables <- c("patients", "admissions", "diagnoses_icd", "procedures_icd",
                   "d_icd_diagnoses", "d_icd_procedures", "prescriptions",
                   "pharmacy", "emar", "emar_detail", "microbiologyevents",
                   "labevents", "d_labitems", "transfers", "services",
                   "omr", "poe", "poe_detail", "drgcodes")

  icu_tables  <- c("icustays", "chartevents", "inputevents", "procedureevents",
                   "d_items", "outputevents", "datetimeevents", "ingredientevents",
                   "caregiver")

  all_known <- c(hosp_tables, icu_tables)

  tbl_map <- list()
  for (tname in all_known) {
    # Chercher le fichier dont le nom (sans ext) correspond à tname
    matches <- all_files[grepl(paste0("(^|[/\\\\])", tname,
                                       "\\.(csv\\.gz|csv|parquet|feather)$"),
                                all_files)]
    if (length(matches) > 0) {
      # Préférer csv.gz > csv > parquet > feather
      pref_order <- c("csv.gz", "csv", "parquet", "feather")
      chosen <- NA_character_
      for (ext in pref_order) {
        cand <- matches[grepl(paste0("\\.", ext, "$"), matches)]
        if (length(cand) > 0) { chosen <- cand[1]; break }
      }
      tbl_map[[tname]] <- chosen
    }
  }

  message("[01_paths] Tables MIMIC-IV détectées (", sum(!is.na(tbl_map)), "/",
          length(all_known), ") :")
  for (nm in names(tbl_map)) {
    message("  ", nm, " -> ", basename(tbl_map[[nm]]))
  }

  # Vérifier tables indispensables
  critical <- c("patients", "admissions", "labevents", "d_labitems")
  med_sources <- c("emar", "emar_detail", "prescriptions", "pharmacy")
  med_avail <- any(med_sources %in% names(tbl_map) &
                     !is.na(unlist(tbl_map[med_sources[med_sources %in%
                                                          names(tbl_map)]])))

  missing_critical <- critical[!critical %in% names(tbl_map) |
                                  is.na(unlist(tbl_map[critical]))]
  if (length(missing_critical) > 0) {
    stop("[01_paths] Tables CRITIQUES manquantes : ",
         paste(missing_critical, collapse = ", "))
  }
  if (!med_avail) {
    stop("[01_paths] Aucune source de médicaments disponible (emar/prescriptions/pharmacy)")
  }

  # Documenter la source médicaments
  med_source_used <- med_sources[med_sources %in% names(tbl_map) &
                                    !is.na(unlist(tbl_map[med_sources[
                                      med_sources %in% names(tbl_map)]]))][1]
  message("[01_paths] Source médicaments utilisée en priorité : ", med_source_used)

  list(
    paths          = tbl_map,
    med_source     = med_source_used,
    mimic_root     = mimic_root,
    hosp_available = hosp_tables[hosp_tables %in% names(tbl_map)],
    icu_available  = icu_tables[icu_tables %in% names(tbl_map)]
  )
}

# -----------------------------------------------------------------------------
# Lecture robuste selon l'extension
# -----------------------------------------------------------------------------
read_mimic_table <- function(path, col_select = NULL, nrows = Inf,
                              use_duckdb = FALSE) {
  if (is.na(path) || is.null(path)) return(NULL)
  if (!file.exists(path)) {
    warning("[01_paths] Fichier introuvable : ", path)
    return(NULL)
  }

  ext <- tolower(sub(".*\\.(.+)$", "\\1", path))
  # Gérer .csv.gz (double extension)
  if (grepl("\\.csv\\.gz$", path, ignore.case = TRUE)) ext <- "csv.gz"

  msg <- paste0("[01_paths] Lecture de ", basename(path),
                " (", round(file.size(path) / 1e6, 1), " MB)")
  message(msg)

  dt <- tryCatch({
    if (ext %in% c("parquet") && PKG_AVAILABLE("arrow")) {
      tbl <- arrow::read_parquet(path,
               col_select = if (!is.null(col_select))
                 dplyr::all_of(col_select) else dplyr::everything())
      data.table::as.data.table(tbl)

    } else if (ext %in% c("feather") && PKG_AVAILABLE("arrow")) {
      tbl <- arrow::read_feather(path)
      data.table::as.data.table(tbl)

    } else if (ext %in% c("csv", "csv.gz")) {
      if (use_duckdb && PKG_AVAILABLE("duckdb") && is.infinite(nrows)) {
        # DuckDB pour grands fichiers
        con <- duckdb::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
        on.exit(duckdb::dbDisconnect(con, shutdown = TRUE), add = TRUE)
        sql <- paste0("SELECT ", if (is.null(col_select)) "*"
                      else paste(col_select, collapse = ", "),
                      " FROM read_csv_auto('", path, "')")
        res <- DBI::dbGetQuery(con, sql)
        data.table::as.data.table(res)
      } else {
        data.table::fread(path,
                          select    = col_select,
                          nrows     = if (is.infinite(nrows)) Inf else nrows,
                          na.strings = c("NA", "", ".", "NULL"),
                          showProgress = FALSE)
      }
    } else {
      # Fallback readr
      readr::read_csv(path, col_select = col_select,
                      n_max = if (is.infinite(nrows)) Inf else nrows,
                      show_col_types = FALSE)
    }
  }, error = function(e) {
    warning("[01_paths] Erreur lecture ", path, " : ", e$message)
    NULL
  })

  if (!is.null(dt)) {
    message("[01_paths]   -> ", nrow(dt), " lignes, ", ncol(dt), " colonnes")
  }
  dt
}

# -----------------------------------------------------------------------------
# Sauvegarde / rechargement des objets intermédiaires
# -----------------------------------------------------------------------------
save_derived <- function(obj, name, dir = "output/extraction/derived_data",
                          format = "rds") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(name, ".", format))
  if (format == "rds") {
    saveRDS(obj, path)
  } else if (format == "csv") {
    data.table::fwrite(as.data.frame(obj), path)
  }
  message("[01_paths] Sauvegardé : ", path)
  invisible(path)
}

load_derived <- function(name, dir = "output/extraction/derived_data", format = "rds") {
  path <- file.path(dir, paste0(name, ".", format))
  if (!file.exists(path)) {
    message("[01_paths] Fichier dérivé introuvable : ", path)
    return(NULL)
  }
  if (format == "rds") readRDS(path)
  else data.table::fread(path)
}

# Logger simple
pipeline_log <- function(msg, level = "INFO", log_file = NULL) {
  ts  <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  out <- paste0("[", level, "] ", ts, " | ", msg)
  message(out)
  if (!is.null(log_file)) {
    cat(out, "\n", file = log_file, append = TRUE)
  }
}
