# =============================================================================
# 54_mi.R — Imputation multiple, analyse PRIMAIRE (docs/v3.1_plan.md A5)
# =============================================================================
# Une seule imputation (m = 20) sur une cohorte PARENTE (la cohorte éligible
# v3.1, avant restriction à l'outcome évaluable). Toutes les analyses qui en
# dérivent (primaire, sous-cohortes, outcomes secondaires, contrôles
# négatifs) sous-échantillonnent LES MÊMES jeux imputés via
# get_analysis_data(), ce qui garantit la cohérence des N (A1).
#
# Modèle d'imputation : pmm pour toutes les variables incomplètes ; prédicteurs
# = toutes les covariables du PS + traitement + outcomes (primaire,
# secondaires, contrôles négatifs). Les outcomes incomplets sont imputés
# uniquement pour servir de prédicteurs, puis les patients sans outcome
# observé sont retirés de chaque analyse (« multiple imputation, then
# deletion », von Hippel 2007).

MI_EXTRA_PREDICTORS_V31 <- c(
  "aki_d7_v31", "aki_stage2_d7_v31", "death_30d", "rrt_7d",
  "early_creat_rise_12h", "nc_hypoglycemia_d7"
)

impute_cohort_v31 <- function(parent, ps_covs, m = 20, maxit = 20, seed = SEED_V3,
                               cache = NULL, extra_predictors = MI_EXTRA_PREDICTORS_V31) {
  if (!is.null(cache) && file.exists(cache)) {
    message("[54_mi_v31] Imputations : cache ", cache)
    return(readRDS(cache))
  }
  vars <- unique(c(ps_covs, intersect(extra_predictors, names(parent)), "treatment"))
  df <- as.data.frame(parent[, vars, with = FALSE])
  for (v in names(df)) if (is.character(df[[v]])) df[[v]] <- factor(df[[v]])
  meth <- mice::make.method(df)
  meth[meth != ""] <- "pmm"
  pred <- mice::make.predictorMatrix(df)
  t0 <- Sys.time()
  mids <- mice::mice(df, m = m, maxit = maxit, method = meth, predictorMatrix = pred,
                     seed = seed, printFlag = FALSE)
  message("[54_mi_v31] mice : m=", m, ", maxit=", maxit, ", variables imputées : ",
          paste(names(meth)[meth != ""], collapse = ", "), " (",
          round(difftime(Sys.time(), t0, units = "mins"), 1), " min)")
  # Jeux complétés : covariables imputées, OUTCOMES OBSERVÉS restaurés (NA
  # conservés pour la suppression post-imputation), identifiants conservés
  imp_cols <- intersect(ps_covs, names(meth)[meth != ""])
  imp_list <- lapply(seq_len(m), function(i) {
    ci <- mice::complete(mids, i)
    d <- copy(parent)
    for (v in imp_cols) d[[v]] <- if (is.factor(ci[[v]])) as.character(ci[[v]]) else ci[[v]]
    d
  })
  out <- list(imp = imp_list, mids_loggedEvents = mids$loggedEvents, method = meth,
              m = m, maxit = maxit, parent_n = nrow(parent))
  if (!is.null(cache)) saveRDS(out, cache)
  out
}
