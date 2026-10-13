# =============================================================================
# 30_results_master.R — Source unique de vérité pour tous les chiffres v3
# =============================================================================
# Corrige la règle #3 du cahier des charges : chaque estimation est calculée
# UNE FOIS et stockée ici ; le Rmd v3 ne fera que lire ce fichier (R inline),
# sans aucun nombre codé en dur. Voir docs/v3_plan.md §4.
#
# Schéma fixe (une ligne = une analyse identifiable de façon unique par
# `analysis_id`) :
#   analysis_id, comparator, phase, role, estimand, weight_scheme, outcome,
#   description,
#   n, n_ptz, n_carb, n_events,
#   crude_risk_ptz, crude_risk_carb,
#   weighted_risk_ptz, weighted_risk_carb,
#   rd, rd_ci_low, rd_ci_high,
#   rr, rr_ci_low, rr_ci_high,
#   ess, ess_pct, max_smd,
#   n_boot, seed, alpha, lambda,
#   notes

RESULTS_MASTER_SCHEMA <- c(
  "analysis_id", "comparator", "phase", "role", "estimand", "weight_scheme",
  "outcome", "description",
  "n", "n_ptz", "n_carb", "n_events",
  "crude_risk_ptz", "crude_risk_carb",
  "weighted_risk_ptz", "weighted_risk_carb",
  "rd", "rd_ci_low", "rd_ci_high",
  "rr", "rr_ci_low", "rr_ci_high",
  "ess", "ess_pct", "max_smd",
  "n_boot", "seed", "alpha", "lambda",
  "notes"
)

#' Créer une table results_master vide, avec le schéma fixé
init_results_master <- function() {
  dt <- data.table::data.table(matrix(nrow = 0, ncol = length(RESULTS_MASTER_SCHEMA)))
  data.table::setnames(dt, RESULTS_MASTER_SCHEMA)
  dt
}

#' Ajouter une ligne à results_master à partir de la sortie de boot_estimate()
#'
#' @param master        data.table existant (init_results_master() ou accumulé)
#' @param analysis_id   identifiant unique et stable (p.ex. "primary_ato",
#'                       "vanco_window_12h", "design_cohort_B") — utilisé pour
#'                       vérifier qu'une même analyse n'apparaît jamais deux
#'                       fois avec des valeurs différentes (QA-1.3)
#' @param est           liste retournée par boot_estimate()
#' @param ...           autres champs du schéma (comparator, phase, role,
#'                       estimand, outcome, description, notes, ...)
add_result_row <- function(master, analysis_id, est, comparator = "carbapenem",
                            phase = "A", role = "OTHER", estimand = NA_character_,
                            outcome = "outcome_primary", description = "",
                            notes = "") {
  if (analysis_id %in% master$analysis_id) {
    stop("[30_results_master_v3] analysis_id déjà présent (violerait QA-1.3, ",
         "une estimation = un IC) : ", analysis_id)
  }
  row <- data.table::data.table(
    analysis_id = analysis_id, comparator = comparator, phase = phase,
    role = role, estimand = estimand,
    weight_scheme = if (!is.null(est$weight_scheme)) est$weight_scheme else NA_character_,
    outcome = outcome, description = description,
    n = est$n, n_ptz = est$n_ptz, n_carb = est$n_carb, n_events = est$n_events,
    crude_risk_ptz = est$crude_risk_ptz, crude_risk_carb = est$crude_risk_carb,
    weighted_risk_ptz = est$weighted_risk_ptz, weighted_risk_carb = est$weighted_risk_carb,
    rd = est$rd, rd_ci_low = est$rd_ci_low, rd_ci_high = est$rd_ci_high,
    rr = est$rr, rr_ci_low = est$rr_ci_low, rr_ci_high = est$rr_ci_high,
    ess = est$ess, ess_pct = est$ess_pct, max_smd = est$max_smd,
    n_boot = est$n_boot, seed = est$seed, alpha = est$alpha, lambda = est$lambda,
    notes = notes
  )
  data.table::rbindlist(list(master, row), fill = TRUE)
}

#' Sauvegarder results_master en .rds et .csv (une seule source de vérité,
#' relue par le Rmd v3 et par R/qa_checks_v3.R)
save_results_master <- function(master, dir = "output/extraction/results") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(master, file.path(dir, "results_master.rds"))
  data.table::fwrite(master, file.path(dir, "results_master.csv"))
  message("[30_results_master_v3] ", nrow(master), " lignes sauvegardées dans ",
          dir, "/results_master.{rds,csv}")
  invisible(master)
}

#' Lire results_master (v3)
load_results_master <- function(dir = "output/extraction/results") {
  path <- file.path(dir, "results_master.rds")
  if (!file.exists(path)) {
    message("[30_results_master_v3] introuvable : ", path)
    return(NULL)
  }
  readRDS(path)
}

#' Construire manuscript_numbers.json : tous les chiffres cités dans le texte
#' du manuscrit, sous forme de clés nommées (N, %, RD, IC, ESS, SMD max, ...),
#' pour que le Rmd v3 les lise par nom plutôt que par position dans un tableau
#' (contrainte du cahier des charges §"une seule source de vérité").
build_manuscript_numbers_json <- function(master, dir = "output/extraction/results") {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    message("[30_results_master_v3] jsonlite indisponible — manuscript_numbers.json non généré")
    return(invisible(NULL))
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  numbers <- list()
  for (i in seq_len(nrow(master))) {
    r <- master[i]
    key_prefix <- r$analysis_id
    numbers[[paste0(key_prefix, "_n")]]              <- r$n
    numbers[[paste0(key_prefix, "_n_ptz")]]           <- r$n_ptz
    numbers[[paste0(key_prefix, "_n_carb")]]          <- r$n_carb
    numbers[[paste0(key_prefix, "_n_events")]]        <- r$n_events
    numbers[[paste0(key_prefix, "_weighted_risk_ptz")]]  <- r$weighted_risk_ptz
    numbers[[paste0(key_prefix, "_weighted_risk_carb")]] <- r$weighted_risk_carb
    numbers[[paste0(key_prefix, "_rd_pct")]]          <- round(r$rd * 100, 2)
    numbers[[paste0(key_prefix, "_rd_ci_low_pct")]]   <- round(r$rd_ci_low * 100, 2)
    numbers[[paste0(key_prefix, "_rd_ci_high_pct")]]  <- round(r$rd_ci_high * 100, 2)
    numbers[[paste0(key_prefix, "_rr")]]              <- round(r$rr, 3)
    numbers[[paste0(key_prefix, "_rr_ci_low")]]       <- round(r$rr_ci_low, 3)
    numbers[[paste0(key_prefix, "_rr_ci_high")]]      <- round(r$rr_ci_high, 3)
    numbers[[paste0(key_prefix, "_ess_pct")]]         <- round(r$ess_pct, 1)
    numbers[[paste0(key_prefix, "_max_smd")]]         <- round(r$max_smd, 3)
  }

  path <- file.path(dir, "manuscript_numbers.json")
  jsonlite::write_json(numbers, path, auto_unbox = TRUE, pretty = TRUE,
                       na = "null", digits = NA)
  message("[30_results_master_v3] ", length(numbers), " clés sauvegardées dans ", path)
  invisible(numbers)
}

#' Enregistrer sessionInfo() + graine (critère d'acceptation QA-5.1)
save_v3_session_info <- function(dir = "output/extraction/results", seed = SEED_V3) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  si <- list(session_info = utils::sessionInfo(), seed_v3 = seed,
             timestamp = Sys.time())
  saveRDS(si, file.path(dir, "session_info.rds"))
  message("[30_results_master_v3] session_info.rds sauvegardé (seed=", seed, ")")
  invisible(si)
}
