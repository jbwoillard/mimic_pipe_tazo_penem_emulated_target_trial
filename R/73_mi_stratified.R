# =============================================================================
# 73_mi_stratified.R — Imputation multiple stratifiée par lieu de soins
#   (docs/v3.2_plan.md §1.3, B3)
# =============================================================================
# En v3.1, un seul modèle d'imputation pour toute la cohorte : les constantes
# des patients de salle (manquantes car non surveillés en réanimation) étaient
# imputées à partir de patients de réanimation (FC médiane observée 91 ->
# 102 après imputation et pondération). v3.2 : un modèle mice par lieu de
# soins à T0 (réanimation / urgences / salle), chaque strate imputée
# uniquement à partir de ses propres patients ; assemblage par imputation.
# MAR reste une hypothèse (conditionnelle au lieu et aux covariables) : elle
# est discutée dans les méthodes et encadrée par les sensibilités « médiane +
# indicateurs » et « cas complets sur constantes ».

impute_stratified_v32 <- function(parent, ps_covs, strata_col, m = 20, maxit = 20, extra_predictors,
                                  cache = NULL) {
  if (!is.null(cache) && file.exists(cache)) { message("[73_mi_v32] cache ", cache); return(readRDS(cache)) }
  strata <- split(parent, by = strata_col, keep.by = TRUE)
  per <- lapply(names(strata), function(s) {
    message("[73_mi_v32] strate ", s, " : n = ", nrow(strata[[s]]))
    impute_cohort_v31(strata[[s]], ps_covs, m = m, maxit = maxit, seed = SEED_V3 + match(s, names(strata)),
                      extra_predictors = extra_predictors)
  })
  imp <- lapply(seq_len(m), function(i) setkey(rbindlist(lapply(per, function(p) p$imp[[i]]), use.names = TRUE), hadm_id))
  out <- list(imp = imp, strata = names(strata), n_by_stratum = vapply(strata, nrow, integer(1)), m = m, maxit = maxit,
              logged = lapply(per, `[[`, "mids_loggedEvents"))
  if (!is.null(cache)) saveRDS(out, cache)
  out
}

#' Distribution observée vs imputée, par lieu de soins (documente le
#' déplacement des distributions dû à l'imputation)
observed_vs_imputed_v32 <- function(raw, imp_list, vars, strata_col) {
  rbindlist(lapply(vars, function(v) {
    miss <- is.na(raw[[v]])
    rbindlist(lapply(unique(raw[[strata_col]]), function(s) {
      i <- raw[[strata_col]] == s
      imputed_vals <- unlist(lapply(imp_list, function(d) d[J(raw$hadm_id[i & miss])][[v]]))
      data.table(variable = v, stratum = s, n = sum(i), pct_missing = round(100 * mean(miss[i]), 1),
                 median_observed = stats::median(raw[[v]][i & !miss]), median_imputed = if (length(imputed_vals)) stats::median(imputed_vals) else NA_real_)
    }))
  }))
}
