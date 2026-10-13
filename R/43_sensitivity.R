# =============================================================================
# 43_sensitivity.R — E-value (§3.5) + analyse de biais quantitative (§3.5)
# =============================================================================

#' E-value pour le RR pondéré et sa borne d'IC la plus proche de 1
compute_evalue_v3 <- function(rr, rr_ci_low, rr_ci_high) {
  if (!requireNamespace("EValue", quietly = TRUE)) {
    message("[43_sensitivity_v3] Package EValue indisponible — E-value non calculée")
    return(NULL)
  }
  ev <- EValue::evalues.RR(est = rr, lo = rr_ci_low, hi = rr_ci_high)
  tbl <- as.data.table(ev, keep.rownames = "component")
  data.table::fwrite(tbl, v3_tables_path("evalue_primary_v3.csv"))
  message("[43_sensitivity_v3] E-value (point) = ", round(ev["E-values", "point"], 2),
          " | E-value (borne IC la plus proche de 1) = ",
          round(ev["E-values", ifelse(!is.na(ev["E-values","lower"]), "lower", "upper")], 2))
  print(ev)
  list(table = tbl, evalue_point = ev["E-values", "point"], raw = ev)
}

#' Analyse de biais quantitative simple : confondeur binaire U non mesuré,
#' avec prévalence p_u, RR(U->outcome) = rr_uy, RR(U->PTZ) = rr_ua (différence
#' de prévalence de U entre bras). Formule du "bias factor" de VanderWeele &
#' Ding (2017), la même que sous-tend l'E-value, appliquée à une grille de
#' scénarios concrets (docs/v3_plan.md §3.5 : BLSE, allergie pénicilline).
run_quantitative_bias_analysis_v3 <- function(rr_observed, scenarios) {
  bias_factor <- function(rr_ua, rr_uy) (rr_ua * rr_uy) / (rr_ua + rr_uy - 1)
  rows <- lapply(seq_len(nrow(scenarios)), function(i) {
    s <- scenarios[i]
    bf <- bias_factor(s$rr_confounder_exposure, s$rr_confounder_outcome)
    rr_adjusted <- rr_observed / bf
    data.table(scenario = s$scenario, rr_confounder_exposure = s$rr_confounder_exposure,
               rr_confounder_outcome = s$rr_confounder_outcome, bias_factor = round(bf, 3),
               rr_observed = round(rr_observed, 3), rr_bias_adjusted = round(rr_adjusted, 3),
               would_nullify = rr_adjusted <= 1)
  })
  tbl <- rbindlist(rows)
  data.table::fwrite(tbl, v3_tables_path("quantitative_bias_analysis_v3.csv"))
  message("[43_sensitivity_v3] Analyse de biais quantitative :")
  print(tbl)
  tbl
}
