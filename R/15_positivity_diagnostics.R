# =============================================================================
# 15_positivity_diagnostics.R — Audit de positivité et sélection de l'estimand
# =============================================================================

# PS principal : régression logistique pénalisée (moins flexible que ranger)
# → évite la séparation quasi-parfaite artificielle

# Variables PS primaires (vanco_6h au lieu de vancomycin_concomitant=vanco_48h)
PS_COVARIATES_PRIMARY <- c(
  "age", "gender",
  "baseline_creat", "egfr_baseline", "bun_baseline",
  "ckd", "diabetes", "heart_failure", "liver_disease", "cancer",
  "charlson_score",
  "vasopressor_pre_t0", "vent_pre_t0",
  "map_mean", "hr_mean", "temp_max", "spo2_min", "gcs_min",
  "sepsis_composite",
  "icu_before_t0",
  "vanco_6h",                        # grace period 6h seulement (vs 48h)
  "aminoglycoside_pre_t0", "other_nephrotoxin_pre_t0",
  "positive_culture_pre_t0", "prior_antibiotic",
  "time_admission_to_T0_h",          # délai admission → traitement
  "calendar_year"
)

# -----------------------------------------------------------------------
# Ajouter variables de timing
# -----------------------------------------------------------------------
add_timing_vars <- function(cohort) {
  dt <- copy(cohort)

  # Délai admission → T0 (en heures)
  if ("admittime" %in% names(dt) && "T0" %in% names(dt)) {
    dt[, time_admission_to_T0_h := as.numeric(
      difftime(T0, admittime, units = "hours")
    )]
    # Valeurs négatives ou extrêmes → NA
    dt[time_admission_to_T0_h < 0 | time_admission_to_T0_h > 30 * 24,
       time_admission_to_T0_h := NA_real_]
    dt[is.na(time_admission_to_T0_h),
       time_admission_to_T0_h := median(dt$time_admission_to_T0_h, na.rm = TRUE)]
  } else {
    dt[, time_admission_to_T0_h := NA_real_]
  }

  # Délai ICU → T0 (en heures)
  if ("icu_intime" %in% names(dt) && "T0" %in% names(dt)) {
    dt[, time_icu_to_T0_h := as.numeric(difftime(T0, icu_intime, units = "hours"))]
    dt[time_icu_to_T0_h < 0, time_icu_to_T0_h := NA_real_]
    dt[time_icu_to_T0_h > 30 * 24, time_icu_to_T0_h := NA_real_]
    med_icu <- median(dt$time_icu_to_T0_h, na.rm = TRUE)
    dt[is.na(time_icu_to_T0_h), time_icu_to_T0_h := med_icu]
  } else {
    dt[, time_icu_to_T0_h := NA_real_]
  }

  dt
}

# -----------------------------------------------------------------------
# PS logistique (principal — moins flexible)
# -----------------------------------------------------------------------
fit_ps_logistic <- function(dt, ps_covs = PS_COVARIATES_PRIMARY) {
  vars_avail <- ps_covs[ps_covs %in% names(dt)]
  trt        <- as.integer(dt$treatment == "PTZ")
  X          <- prepare_X(dt, vars_avail)
  df         <- data.frame(trt = trt, X)

  ps <- NULL

  # Essai 1 : glmnet elastic net (évite la séparation parfaite)
  if (PKG_AVAILABLE("glmnet")) {
    tryCatch({
      cv_fit <- glmnet::cv.glmnet(X, trt, family = "binomial",
                                   alpha = 0.1, nfolds = 5)
      ps <- as.vector(predict(cv_fit, X, s = "lambda.min",
                               type = "response"))
      message("[15_pos] PS logistique (glmnet ridge) ajusté")
    }, error = function(e) message("[15_pos] glmnet erreur : ", e$message))
  }

  # Fallback : GLM logistique standard
  if (is.null(ps)) {
    fmla <- as.formula(paste("trt ~", paste(vars_avail, collapse = " + ")))
    tryCatch({
      fit_glm <- glm(fmla, data = df, family = binomial())
      ps      <- predict(fit_glm, type = "response")
      message("[15_pos] PS logistique (GLM standard) ajusté")
    }, error = function(e) {
      message("[15_pos] GLM erreur : ", e$message)
      ps <- rep(mean(trt), nrow(dt))
    })
  }

  # Clip pour stabilité
  ps <- pmax(pmin(ps, 0.99), 0.01)
  ps
}

# -----------------------------------------------------------------------
# PS ranger (secondaire — plus flexible)
# -----------------------------------------------------------------------
fit_ps_ranger <- function(dt, ps_covs = PS_COVARIATES_PRIMARY) {
  vars_avail <- ps_covs[ps_covs %in% names(dt)]
  trt        <- as.integer(dt$treatment == "PTZ")
  X          <- prepare_X(dt, vars_avail)
  ps         <- NULL

  if (PKG_AVAILABLE("ranger")) {
    tryCatch({
      fit_rf <- ranger::ranger(y = factor(trt), x = data.frame(X),
                                num.trees = 500, probability = TRUE,
                                min.node.size = 10)
      ps <- predict(fit_rf, data.frame(X))$predictions[, "1"]
      message("[15_pos] PS ranger ajusté")
    }, error = function(e) message("[15_pos] ranger erreur : ", e$message))
  }

  if (is.null(ps)) ps <- fit_ps_logistic(dt, ps_covs)
  pmax(pmin(ps, 0.99), 0.01)
}

# -----------------------------------------------------------------------
# Fonction principale : diagnostic de positivité
# -----------------------------------------------------------------------
diagnose_positivity <- function(data, ps, covariates = PS_COVARIATES_PRIMARY,
                                 outcome_col = "outcome_primary",
                                 label = "primary") {
  dt  <- data[!is.na(data[[outcome_col]]), ]
  trt <- as.integer(dt$treatment == "PTZ")
  y   <- as.integer(dt[[outcome_col]])
  n   <- nrow(dt)
  n1  <- sum(trt); n0 <- n - n1

  ps1 <- ps[trt == 1]; ps0 <- ps[trt == 0]

  # Quantiles du PS par groupe
  qq <- c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1.0)
  q1 <- quantile(ps1, qq, na.rm = TRUE)
  q0 <- quantile(ps0, qq, na.rm = TRUE)

  # Support commun brut
  cs_lo  <- max(min(ps1), min(ps0))
  cs_hi  <- min(max(ps1), max(ps0))
  in_cs  <- sum(ps >= cs_lo & ps <= cs_hi)
  pct_cs <- round(in_cs / n * 100, 1)

  # Support commun stable (p5-p95)
  cs_lo5  <- max(quantile(ps1, 0.05), quantile(ps0, 0.05))
  cs_hi5  <- min(quantile(ps1, 0.95), quantile(ps0, 0.95))
  in_cs5  <- sum(ps >= cs_lo5 & ps <= cs_hi5)
  pct_cs5 <- round(in_cs5 / n * 100, 1)

  cs_lo10  <- max(quantile(ps1, 0.10), quantile(ps0, 0.10))
  cs_hi10  <- min(quantile(ps1, 0.90), quantile(ps0, 0.90))
  in_cs10  <- sum(ps >= cs_lo10 & ps <= cs_hi10)
  pct_cs10 <- round(in_cs10 / n * 100, 1)

  # Extremes
  p_lo05 <- round(mean(ps < 0.05) * 100, 1)
  p_hi95 <- round(mean(ps > 0.95) * 100, 1)
  p_lo10 <- round(mean(ps < 0.10) * 100, 1)
  p_hi90 <- round(mean(ps > 0.90) * 100, 1)

  # ESS IPTW
  p_trt <- mean(trt)
  w_iptw <- ifelse(trt == 1, p_trt / ps, (1 - p_trt) / (1 - ps))
  ess_iptw <- (sum(w_iptw))^2 / sum(w_iptw^2)
  ess_iptw_pct <- round(ess_iptw / n * 100, 1)

  # ESS overlap
  w_ow <- ifelse(trt == 1, 1 - ps, ps)
  ess_ow <- (sum(w_ow))^2 / sum(w_ow^2)
  ess_ow_pct <- round(ess_ow / n * 100, 1)

  # AUC du PS
  auc_ps <- NA_real_
  if (PKG_AVAILABLE("pROC")) {
    tryCatch({
      roc_obj <- pROC::roc(trt, ps, quiet = TRUE)
      auc_ps  <- as.numeric(roc_obj$auc)
    }, error = function(e) NULL)
  }

  # SMD pour variables clés
  smd_fn <- function(v) {
    if (!v %in% names(dt)) return(NA_real_)
    x <- as.numeric(dt[[v]])
    x1 <- x[trt == 1]; x0 <- x[trt == 0]
    sd_pool <- sqrt((var(x1, na.rm = TRUE) + var(x0, na.rm = TRUE)) / 2)
    if (is.na(sd_pool) || sd_pool == 0) return(NA_real_)
    abs(mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / sd_pool
  }

  key_vars <- c("vanco_6h", "vancomycin_concomitant", "vanco_48h",
                "time_admission_to_T0_h", "hosp_los_before_t0_h",
                "time_icu_to_T0_h", "prior_antibiotic",
                "sepsis_composite", "baseline_creat", "egfr_baseline",
                "vasopressor_pre_t0", "vent_pre_t0")
  smd_vals <- sapply(key_vars, smd_fn)

  # SMD global (max et médiane)
  cov_avail <- covariates[covariates %in% names(dt)]
  smd_all   <- sapply(cov_avail, smd_fn)
  smd_max   <- max(smd_all, na.rm = TRUE)
  smd_med   <- median(smd_all, na.rm = TRUE)

  # SMD après IPTW
  smd_iptw_all <- tryCatch({
    sapply(cov_avail, function(v) {
      if (!v %in% names(dt)) return(NA_real_)
      x <- as.numeric(dt[[v]])
      x1 <- x[trt == 1]; x0 <- x[trt == 0]
      w1 <- w_iptw[trt == 1]; w0 <- w_iptw[trt == 0]
      m1 <- weighted.mean(x1, w1, na.rm = TRUE)
      m0 <- weighted.mean(x0, w0, na.rm = TRUE)
      sd_pool <- sqrt((var(x1, na.rm = TRUE) + var(x0, na.rm = TRUE)) / 2)
      if (is.na(sd_pool) || sd_pool == 0) return(NA_real_)
      abs(m1 - m0) / sd_pool
    })
  }, error = function(e) rep(NA_real_, length(cov_avail)))
  smd_iptw_max <- max(smd_iptw_all, na.rm = TRUE)

  # SMD après overlap weighting
  smd_ow_all <- tryCatch({
    sapply(cov_avail, function(v) {
      if (!v %in% names(dt)) return(NA_real_)
      x <- as.numeric(dt[[v]])
      x1 <- x[trt == 1]; x0 <- x[trt == 0]
      w1 <- w_ow[trt == 1]; w0 <- w_ow[trt == 0]
      m1 <- weighted.mean(x1, w1, na.rm = TRUE)
      m0 <- weighted.mean(x0, w0, na.rm = TRUE)
      sd_pool <- sqrt((var(x1, na.rm = TRUE) + var(x0, na.rm = TRUE)) / 2)
      if (is.na(sd_pool) || sd_pool == 0) return(NA_real_)
      abs(m1 - m0) / sd_pool
    })
  }, error = function(e) rep(NA_real_, length(cov_avail)))
  smd_ow_max <- max(smd_ow_all, na.rm = TRUE)

  # Decider de l'estimand
  est_choice <- select_primary_estimand(
    pct_cs = pct_cs5, ess_iptw_pct = ess_iptw_pct,
    smd_iptw_max = smd_iptw_max, smd_ow_max = smd_ow_max,
    auc_ps = auc_ps
  )

  # Construire la table summary
  summary_row <- data.table(
    label               = label,
    n_total             = n,
    n_ptz               = n1,
    n_carb              = n0,
    ps_min              = round(min(ps), 4),
    ps_p01              = round(quantile(ps, 0.01), 4),
    ps_p05              = round(quantile(ps, 0.05), 4),
    ps_p10              = round(quantile(ps, 0.10), 4),
    ps_p25              = round(quantile(ps, 0.25), 4),
    ps_p50              = round(quantile(ps, 0.50), 4),
    ps_p75              = round(quantile(ps, 0.75), 4),
    ps_p90              = round(quantile(ps, 0.90), 4),
    ps_p95              = round(quantile(ps, 0.95), 4),
    ps_p99              = round(quantile(ps, 0.99), 4),
    ps_max              = round(max(ps), 4),
    ps_min_ptz          = round(min(ps1), 4),
    ps_max_ptz          = round(max(ps1), 4),
    ps_min_carb         = round(min(ps0), 4),
    ps_max_carb         = round(max(ps0), 4),
    pct_ps_lt05         = p_lo05,
    pct_ps_gt95         = p_hi95,
    pct_ps_lt10         = p_lo10,
    pct_ps_gt90         = p_hi90,
    n_common_support    = in_cs,
    pct_common_support  = pct_cs,
    n_trim_p5p95        = in_cs5,
    pct_trim_p5p95      = pct_cs5,
    n_trim_p10p90       = in_cs10,
    pct_trim_p10p90     = pct_cs10,
    ess_iptw            = round(ess_iptw),
    ess_iptw_pct        = ess_iptw_pct,
    ess_ow              = round(ess_ow),
    ess_ow_pct          = ess_ow_pct,
    iptw_w_max          = round(max(w_iptw), 2),
    iptw_w_p95          = round(quantile(w_iptw, 0.95), 2),
    iptw_w_p99          = round(quantile(w_iptw, 0.99), 2),
    pct_extreme_weights = round(mean(w_iptw > quantile(w_iptw, 0.99)) * 100, 1),
    auc_ps              = round(auc_ps, 3),
    smd_max_crude       = round(smd_max, 3),
    smd_med_crude       = round(smd_med, 3),
    smd_max_iptw        = round(smd_iptw_max, 3),
    smd_max_ow          = round(smd_ow_max, 3),
    smd_vanco6h         = round(smd_vals["vanco_6h"], 3),
    smd_vanco48h        = round(smd_vals["vancomycin_concomitant"], 3),
    smd_time_to_T0      = round(smd_vals["time_admission_to_T0_h"], 3),
    smd_hosp_los_pre    = round(smd_vals["hosp_los_before_t0_h"], 3),
    smd_prior_abx       = round(smd_vals["prior_antibiotic"], 3),
    smd_sepsis          = round(smd_vals["sepsis_composite"], 3),
    smd_creat_baseline  = round(smd_vals["baseline_creat"], 3),
    smd_vasopressor     = round(smd_vals["vasopressor_pre_t0"], 3),
    smd_vent            = round(smd_vals["vent_pre_t0"], 3),
    primary_estimand    = est_choice$estimand,
    positivity_verdict  = est_choice$verdict
  )

  # Warnings
  warnings <- character(0)
  if (pct_cs5 < 50 || p_lo05 + p_hi95 > 20)
    warnings <- c(warnings, "SEVERE: Severe practical positivity violation")
  else if (pct_cs5 < 75 || ess_iptw_pct < 30)
    warnings <- c(warnings, "MODERATE: Moderate positivity concern")

  if (ess_iptw_pct < 30 || max(w_iptw) > 20)
    warnings <- c(warnings, "ATE not recommended: extreme IPTW weights")

  if (smd_ow_max < smd_iptw_max * 0.8)
    warnings <- c(warnings, "ATO recommended: OW improves balance over IPTW")

  if (length(warnings) == 0)
    warnings <- "OK: No major positivity concern detected"

  for (w in warnings) message("[15_pos] ", w)

  list(
    summary_table = summary_row,
    warnings      = warnings,
    ps            = ps,
    w_iptw        = w_iptw,
    w_ow          = w_ow,
    cs_bounds     = list(p0p100 = c(cs_lo, cs_hi),
                         p5p95  = c(cs_lo5, cs_hi5),
                         p10p90 = c(cs_lo10, cs_hi10)),
    estimand_recommendation = est_choice
  )
}

# -----------------------------------------------------------------------
# Sélection automatique de l'estimand
# -----------------------------------------------------------------------
select_primary_estimand <- function(pct_cs, ess_iptw_pct,
                                     smd_iptw_max, smd_ow_max, auc_ps = NA) {
  if (pct_cs >= 85 && ess_iptw_pct >= 50 && smd_iptw_max <= 0.10) {
    list(estimand = "ATE", verdict = "OK: Good overlap, ATE identified",
         rationale = "Common support > 85%, ESS > 50%, max SMD after IPTW ≤ 0.10")
  } else if (pct_cs >= 60 && smd_ow_max <= 0.10) {
    list(estimand = "ATO",
         verdict = "MODERATE: Limited overlap — ATO (overlap-weighted) recommended",
         rationale = "Common support 60-85% or high IPTW weights; OW achieves SMD ≤ 0.10")
  } else if (pct_cs >= 40) {
    list(estimand = "ATO+COMMON_SUPPORT",
         verdict = "POOR: Poor overlap — ATO primary; common-support restriction also recommended",
         rationale = "Common support < 60%; ATE requires extrapolation")
  } else {
    list(estimand = "UNKNOWN",
         verdict = "SEVERE: Insufficient overlap for reliable causal inference",
         rationale = "Common support < 40%; results should be interpreted with extreme caution")
  }
}

# -----------------------------------------------------------------------
# Audit de timing des covariables
# -----------------------------------------------------------------------
audit_covariate_timing <- function() {
  timing_tbl <- data.table(
    variable = c(
      "age", "gender", "baseline_creat", "egfr_baseline", "bun_baseline",
      "ckd", "diabetes", "heart_failure", "liver_disease", "cancer",
      "charlson_score", "vasopressor_pre_t0", "vent_pre_t0",
      "map_mean", "hr_mean", "temp_max", "spo2_min", "gcs_min",
      "sepsis_composite", "icu_before_t0",
      "vanco_6h", "vanco_12h", "vanco_24h", "vanco_48h",
      "vancomycin_concomitant",
      "aminoglycoside_pre_t0", "other_nephrotoxin_pre_t0",
      "positive_culture_pre_t0", "prior_antibiotic",
      "hosp_los_before_t0_h", "calendar_year",
      "time_admission_to_T0_h", "time_icu_to_T0_h",
      "rrt_pre_t0", "esrd_icd",
      "aki_at_t0", "aki_d2", "aki_d7", "make30",
      "death_30d", "time_from_T0_to_disch"
    ),
    used_in_primary_ps = c(
      TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE,
      TRUE, FALSE, FALSE, FALSE,   # vanco_6h only in primary PS
      FALSE,                        # vancomycin_concomitant = vanco_48h exclu
      TRUE, TRUE,
      TRUE, TRUE,
      TRUE, TRUE,
      TRUE, TRUE,
      FALSE, FALSE,                 # exclusion criteria, not PS covariates
      FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE
    ),
    measurement_window = c(
      "at admission", "at admission", "pre-T0", "derived from creat", "pre-T0",
      "ICD codes (lifetime)", "ICD codes", "ICD codes", "ICD codes", "ICD codes",
      "ICD codes", "T0-24h to T0", "T0-24h to T0",
      "T0-24h to T0", "T0-24h to T0", "T0-24h to T0", "T0-24h to T0", "T0-24h to T0",
      "ICD codes + vitals T0-24h", "at T0",
      "T0-24h to T0+6h", "T0-24h to T0+12h", "T0-24h to T0+24h", "T0-24h to T0+48h",
      "T0-24h to T0+48h (=vanco_48h)",
      "T0-48h to T0", "T0-48h to T0",
      "pre-T0", "pre-T0",
      "admission to T0", "year of T0",
      "admission to T0", "ICU admission to T0",
      "pre-T0 (exclusion criterion)", "ICD codes (exclusion criterion)",
      "at T0", "T0 to T0+2d", "T0 to T0+7d", "T0 to T0+30d",
      "T0 to T0+30d", "T0 to discharge"
    ),
    allowed_pre_t0 = c(
      TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE,
      TRUE, TRUE, TRUE, FALSE,     # vanco_48h has post-T0 component
      FALSE,                        # vancomycin_concomitant = vanco_48h
      TRUE, TRUE,
      TRUE, TRUE,
      TRUE, TRUE,
      TRUE, TRUE,
      TRUE, FALSE,                  # rrt_pre_t0 TRUE, esrd_icd FALSE (ICD ok)
      FALSE, FALSE, FALSE, FALSE,   # aki_at_t0, aki_d2, aki_d7, make30
      FALSE, FALSE                  # death_30d, time_from_T0_to_disch
    ),
    possible_post_t0_leakage = c(
      FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE,
      FALSE, FALSE, FALSE, TRUE,   # vanco_48h: 42h post-T0
      TRUE,                         # vancomycin_concomitant = vanco_48h
      FALSE, FALSE,
      FALSE, FALSE,
      FALSE, FALSE,
      FALSE, FALSE,
      FALSE, FALSE,                 # rrt_pre_t0, esrd_icd
      FALSE, TRUE, TRUE, TRUE,     # aki_at_t0, aki_d2, aki_d7, make30
      TRUE, TRUE                   # death_30d, time_from_T0_to_disch
    ),
    decision = c(
      "keep", "keep", "keep", "keep", "keep",
      "keep", "keep", "keep", "keep", "keep",
      "keep", "keep", "keep",
      "keep", "keep", "keep", "keep", "keep",
      "keep", "keep",
      "keep in primary PS", "sensitivity only", "sensitivity only", "EXCLUDED from primary PS",
      "EXCLUDED from primary PS (=vanco_48h)",
      "keep", "keep",
      "keep", "keep",
      "keep", "keep",
      "ADD to primary PS", "add if available",
      "exclusion criterion", "exclusion criterion",
      "baseline check only", "OUTCOME - exclude", "OUTCOME - exclude", "OUTCOME - exclude",
      "OUTCOME - exclude", "OUTCOME - exclude"
    )
  )

  dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(timing_tbl, "output/tables/covariate_timing_audit.csv")
  message("[15_pos] Audit timing covariables sauvegardé")
  timing_tbl
}

# -----------------------------------------------------------------------
# Figure : densité du PS par groupe
# -----------------------------------------------------------------------
plot_ps_density_v2 <- function(dt, ps, label = "primary",
                                 cs_bounds = NULL) {
  if (!PKG_AVAILABLE("ggplot2")) return(invisible(NULL))

  trt <- as.integer(dt$treatment == "PTZ")
  df  <- data.frame(
    ps        = ps,
    treatment = ifelse(trt == 1, "PTZ", "Carbapenems")
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = ps, fill = treatment)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::scale_fill_manual(values = c("PTZ" = "#E41A1C",
                                           "Carbapenems" = "#377EB8")) +
    ggplot2::labs(
      title    = paste0("Propensity score distribution (", label, ")"),
      subtitle = paste0("n = ", nrow(dt), " | PTZ: ",
                        sum(trt), " | Carb: ", sum(1 - trt)),
      x = "Propensity score (P[PTZ|covariates])",
      y = "Density",
      fill = "Treatment"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")

  if (!is.null(cs_bounds)) {
    p <- p +
      ggplot2::geom_vline(xintercept = cs_bounds[1],
                           linetype = "dashed", color = "gray40", alpha = 0.7) +
      ggplot2::geom_vline(xintercept = cs_bounds[2],
                           linetype = "dashed", color = "gray40", alpha = 0.7) +
      ggplot2::annotate("text", x = cs_bounds[1], y = Inf,
                         label = paste0("CS lower: ", round(cs_bounds[1], 2)),
                         vjust = 2, hjust = 0, size = 3, color = "gray40")
  }

  fname <- paste0("output/figures/ps_density_", gsub(" ", "_", label), ".png")
  tryCatch(
    ggplot2::ggsave(fname, p, width = 8, height = 5, dpi = 150),
    error = function(e) message("[15_pos] Erreur sauvegarde ", fname)
  )
  message("[15_pos] Figure PS sauvegardée : ", fname)
  p
}

# -----------------------------------------------------------------------
# Figure : comparaison poids IPTW vs overlap weighting
# -----------------------------------------------------------------------
plot_weights_comparison <- function(dt, w_iptw, w_ow) {
  if (!PKG_AVAILABLE("ggplot2")) return(invisible(NULL))

  # dt is already filtered to outcome_primary non-NA
  trt <- as.integer(dt$treatment == "PTZ")
  df  <- data.frame(
    iptw    = w_iptw,
    ow      = w_ow,
    treatment = ifelse(trt == 1, "PTZ", "Carbapenems")
  )

  # IPTW
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = iptw, fill = treatment)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    ggplot2::scale_fill_manual(values = c("PTZ" = "#E41A1C",
                                           "Carbapenems" = "#377EB8")) +
    ggplot2::labs(title = "IPTW weights (ATE)", x = "Weight", y = "Count") +
    ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "none")

  # OW
  p2 <- ggplot2::ggplot(df, ggplot2::aes(x = ow, fill = treatment)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    ggplot2::scale_fill_manual(values = c("PTZ" = "#E41A1C",
                                           "Carbapenems" = "#377EB8")) +
    ggplot2::labs(title = "Overlap weights (ATO)", x = "Weight", y = "Count",
                  fill = "Treatment") +
    ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "bottom")

  if (PKG_AVAILABLE("patchwork")) {
    p_combined <- p1 / p2
  } else {
    p_combined <- p2
  }

  tryCatch(
    ggplot2::ggsave("output/figures/weights_overlap_vs_iptw.png",
                    p_combined, width = 8, height = 7, dpi = 150),
    error = function(e) message("[15_pos] Erreur weight plot")
  )
  invisible(p_combined)
}

# -----------------------------------------------------------------------
# Love plot avec OW
# -----------------------------------------------------------------------
create_love_plot_ow <- function(dt, ps, w_iptw, w_ow,
                                 covariates = PS_COVARIATES_PRIMARY) {
  if (!PKG_AVAILABLE("cobalt")) return(invisible(NULL))

  trt        <- as.integer(dt$treatment == "PTZ")
  cov_avail  <- covariates[covariates %in% names(dt)]
  df         <- as.data.frame(dt[, cov_avail, with = FALSE])
  for (col in names(df)) df[[col]] <- as.numeric(df[[col]])

  tryCatch({
    bal_obj <- cobalt::bal.tab(
      x = df, treat = trt,
      weights = list(IPTW = w_iptw, Overlap = w_ow),
      stats   = "mean.diffs",
      estimand = "ATO"
    )

    p <- cobalt::love.plot(
      bal_obj,
      thresholds  = c(m = 0.1),
      var.order   = "unadjusted",
      shapes      = c("circle filled", "square filled", "triangle filled"),
      colors      = c("black", "#E41A1C", "#377EB8"),
      line        = FALSE,
      limits      = c(-0.05, max(abs(bal_obj$Balance$Diff.Un), na.rm = TRUE) * 1.1),
      title       = "Love plot: Unadjusted vs IPTW vs Overlap weighting"
    )

    tryCatch(
      ggplot2::ggsave("output/figures/love_plot_overlap_weighting.png",
                      p, width = 9, height = 7, dpi = 150),
      error = function(e) message("[15_pos] Erreur love plot OW")
    )
    message("[15_pos] Love plot OW sauvegardé")
  }, error = function(e) message("[15_pos] cobalt erreur : ", e$message))
}

# -----------------------------------------------------------------------
# Exécuter tout l'audit de positivité
# -----------------------------------------------------------------------
run_positivity_audit <- function(cohort, ps_logistic, ps_ranger = NULL,
                                  ps_original = NULL) {
  dir.create("output/diagnostics", recursive = TRUE, showWarnings = FALSE)
  dir.create("output/figures",     recursive = TRUE, showWarnings = FALSE)

  # Aligner tous les vecteurs PS sur les lignes avec outcome non-NA
  ok_idx <- !is.na(cohort$outcome_primary)
  dt     <- cohort[ok_idx]

  ps_log_f <- ps_logistic[ok_idx]
  ps_rng_f <- if (!is.null(ps_ranger))   ps_ranger[ok_idx]   else NULL
  ps_ori_f <- if (!is.null(ps_original)) ps_original[ok_idx] else NULL

  audit_rows <- list()

  # 1. PS original (ranger flexible)
  if (!is.null(ps_ori_f)) {
    diag_orig <- diagnose_positivity(dt, ps_ori_f, label = "original_ranger")
    audit_rows[["original"]] <- diag_orig$summary_table
    plot_ps_density_v2(dt, ps_ori_f, "original (ranger)",
                       diag_orig$cs_bounds$p5p95)
  }

  # 2. PS logistique (principal)
  diag_log <- diagnose_positivity(dt, ps_log_f, label = "logistic_primary")
  audit_rows[["logistic"]] <- diag_log$summary_table
  plot_ps_density_v2(dt, ps_log_f, "primary logistic",
                     diag_log$cs_bounds$p5p95)

  # 3. PS ranger (secondaire)
  if (!is.null(ps_rng_f)) {
    diag_rf <- diagnose_positivity(dt, ps_rng_f, label = "ranger_secondary")
    audit_rows[["ranger"]] <- diag_rf$summary_table
    plot_ps_density_v2(dt, ps_rng_f, "secondary ranger",
                       diag_rf$cs_bounds$p5p95)
  }

  # Figures poids (dt est déjà filtré sur outcome_primary non-NA)
  plot_weights_comparison(dt, diag_log$w_iptw, diag_log$w_ow)
  create_love_plot_ow(dt, ps_log_f, diag_log$w_iptw, diag_log$w_ow)

  # Sauvegarder audit
  audit_tbl <- rbindlist(audit_rows, fill = TRUE)
  data.table::fwrite(audit_tbl, "output/diagnostics/positivity_audit.csv")
  message("[15_pos] Audit de positivité complet sauvegardé")

  # Audit timing
  audit_covariate_timing()

  list(
    audit_table      = audit_tbl,
    diag_logistic    = diag_log,
    diag_ranger      = if (!is.null(ps_rng_f)) diag_rf else NULL,
    diag_original    = if (!is.null(ps_ori_f)) diag_orig else NULL,
    ps_logistic      = ps_logistic,
    ps_ranger        = ps_ranger,
    primary_estimand = diag_log$estimand_recommendation
  )
}
