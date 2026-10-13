# =============================================================================
# 51_analysis_cohort.R — Porte d'entrée unique vers les cohortes
#   d'analyse + registre « N par analyse » (docs/v3.1_plan.md A1)
# =============================================================================
# En v3, chaque analyse filtrait implicitement sur son propre outcome
# (boot_estimate : !is.na(outcome_col)) et deux outcomes coexistaient
# (outcome_primary phase B / outcome_primary_v3 phase C) -> N impossibles
# (méropénème seul 7 672 > primaire 7 670 ; strates vanco 4 205 + 3 482 =
# 7 687). En v3.1, toute cohorte est définie par une spec, dérive d'un
# parent déclaré, et est enregistrée ; QA-A1 vérifie les inclusions.

.COHORT_REGISTRY_V31 <- new.env()
.COHORT_REGISTRY_V31$ids  <- list()    # cohort_id -> vecteur hadm_id
.COHORT_REGISTRY_V31$rows <- list()    # cohort_id -> ligne de registre
.COHORT_REGISTRY_V31$families <- list() # famille de strates -> (parent, ids des strates)

#' @param spec list(id, parent, filter = expression (évaluée dans master),
#'                  outcome = colonne exigée non NA (ou NULL), label)
#'             parent = "master" pour la racine.
#' @return vecteur de hadm_id (enregistré)
get_analysis_cohort <- function(spec, master) {
  stopifnot(!is.null(spec$id), !is.null(spec$parent))
  if (!is.null(.COHORT_REGISTRY_V31$ids[[spec$id]])) return(.COHORT_REGISTRY_V31$ids[[spec$id]])
  parent_ids <- if (spec$parent == "master") master$hadm_id else {
    p <- .COHORT_REGISTRY_V31$ids[[spec$parent]]
    if (is.null(p)) stop("[51_cohort_v31] parent non enregistré : ", spec$parent)
    p
  }
  d <- master[J(parent_ids)]
  keep <- rep(TRUE, nrow(d))
  if (!is.null(spec$filter)) keep <- keep & (eval(spec$filter, d) %in% TRUE)
  if (!is.null(spec$outcome)) keep <- keep & !is.na(d[[spec$outcome]])
  ids <- d$hadm_id[keep]
  sub <- d[keep]
  .COHORT_REGISTRY_V31$ids[[spec$id]] <- ids
  .COHORT_REGISTRY_V31$rows[[spec$id]] <- data.table(
    cohort_id = spec$id, parent_id = spec$parent,
    label = if (is.null(spec$label)) spec$id else spec$label,
    filter = if (is.null(spec$filter)) "" else paste(deparse(spec$filter), collapse = " "),
    outcome_required = if (is.null(spec$outcome)) "" else spec$outcome,
    n = length(ids), n_ptz = sum(sub$treatment == "PTZ"),
    n_comparator = sum(sub$treatment != "PTZ"), n_parent = length(parent_ids))
  ids
}

#' Déclare une famille de strates qui doit partitionner exactement son parent
register_strata_family <- function(family_id, parent_id, strata_ids) {
  .COHORT_REGISTRY_V31$families[[family_id]] <- list(parent = parent_id, strata = strata_ids)
}

cohort_registry_table <- function() rbindlist(.COHORT_REGISTRY_V31$rows, fill = TRUE)
n_of <- function(cohort_id) {
  r <- .COHORT_REGISTRY_V31$rows[[cohort_id]]
  if (is.null(r)) stop("[51_cohort_v31] cohorte inconnue : ", cohort_id)
  r$n
}

#' Sous-ensemble des m jeux imputés correspondant à une cohorte enregistrée
get_analysis_data <- function(cohort_id, imp_list) {
  ids <- .COHORT_REGISTRY_V31$ids[[cohort_id]]
  if (is.null(ids)) stop("[51_cohort_v31] cohorte inconnue : ", cohort_id)
  lapply(imp_list, function(d) d[J(ids)])
}

#' QA-A1 : inclusion dans le parent et partition exacte des familles de strates
qa_cohort_consistency_v31 <- function() {
  reg <- cohort_registry_table()
  rows <- lapply(seq_len(nrow(reg)), function(i) {
    r <- reg[i]
    ok <- if (r$parent_id == "master") TRUE else
      all(.COHORT_REGISTRY_V31$ids[[r$cohort_id]] %in% .COHORT_REGISTRY_V31$ids[[r$parent_id]]) &&
      r$n <= r$n_parent
    data.table(check = paste0("subset:", r$cohort_id), ok = ok,
               details = sprintf("n=%d <= n(%s)=%d", r$n, r$parent_id, r$n_parent))
  })
  fams <- lapply(names(.COHORT_REGISTRY_V31$families), function(f) {
    fm <- .COHORT_REGISTRY_V31$families[[f]]
    all_ids <- unlist(lapply(fm$strata, function(s) .COHORT_REGISTRY_V31$ids[[s]]))
    par <- .COHORT_REGISTRY_V31$ids[[fm$parent]]
    ok <- length(all_ids) == length(par) && !anyDuplicated(all_ids) && setequal(all_ids, par)
    data.table(check = paste0("partition:", f), ok = ok,
               details = sprintf("sum(strata)=%d vs n(%s)=%d", length(all_ids), fm$parent, length(par)))
  })
  rbindlist(c(rows, fams))
}
