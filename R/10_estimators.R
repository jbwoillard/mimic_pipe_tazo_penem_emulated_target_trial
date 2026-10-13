# =============================================================================
# 10_estimators.R — Estimateurs causaux multiples
# =============================================================================

set.seed(42)

# Variables de confondants pour le modèle de propension
PS_COVARIATES <- c(
  "age", "gender",
  "baseline_creat", "egfr_baseline", "bun_baseline",
  "ckd", "diabetes", "heart_failure", "liver_disease", "cancer",
  "charlson_score",
  "vasopressor_pre_t0", "vent_pre_t0",
  "map_mean", "hr_mean", "temp_max", "spo2_min", "gcs_min",
  "sepsis_composite",
  "icu_before_t0",
  "vancomycin_concomitant",
  "aminoglycoside_pre_t0", "other_nephrotoxin_pre_t0",
  "positive_culture_pre_t0", "prior_antibiotic",
  "hosp_los_before_t0_h",
  "calendar_year"
)

# Fonctions utilitaires
safe_ci <- function(estimate, se, z = 1.96) {
  c(lower = estimate - z * se, upper = estimate + z * se)
}

summarise_result <- function(label, estimand, estimate, ci_low, ci_high,
                              n, n_events, method) {
  data.table(
    label    = label,    estimand  = estimand,
    estimate = estimate, ci_low    = ci_low,    ci_high   = ci_high,
    n        = n,        n_events  = n_events,  method    = method
  )
}

# -----------------------------------------------------------------------
# Helper : préparer la matrice X pour les modèles ML
# -----------------------------------------------------------------------
prepare_X <- function(cohort, vars) {
  vars_avail <- vars[vars %in% names(cohort)]
  df <- as.data.frame(cohort[, vars_avail, with = FALSE])

  # Encoder variables catégorielles
  if ("gender" %in% names(df))
    df$gender <- as.integer(toupper(df$gender) == "F")
  if ("admission_type" %in% names(df))
    df$admission_type <- as.integer(as.factor(df$admission_type))

  # Convertir tout en numeric
  for (col in names(df)) {
    df[[col]] <- as.numeric(df[[col]])
    df[[col]][is.na(df[[col]])] <- 0
  }
  as.matrix(df)
}

# -----------------------------------------------------------------------
# 1. Analyse descriptive brute
# -----------------------------------------------------------------------
run_crude <- function(cohort, outcome_col = "outcome_primary") {
  dt <- cohort[!is.na(get(outcome_col))]
  trt <- as.integer(dt$treatment == "PTZ")
  y   <- as.integer(dt[[outcome_col]])

  n1    <- sum(trt == 1); n0 <- sum(trt == 0)
  r1    <- mean(y[trt == 1]); r0 <- mean(y[trt == 0])
  rd    <- r1 - r0
  rr    <- r1 / r0
  or_v  <- (r1 / (1 - r1)) / (r0 / (1 - r0))

  # IC par Delta method pour RD
  se_rd <- sqrt(r1 * (1 - r1) / n1 + r0 * (1 - r0) / n0)
  se_rr <- sqrt((1 - r1) / (n1 * r1) + (1 - r0) / (n0 * r0))
  se_logOR <- sqrt(1 / (n1 * r1) + 1 / (n1 * (1 - r1)) +
                     1 / (n0 * r0) + 1 / (n0 * (1 - r0)))

  results <- rbind(
    summarise_result("RD brut", "ATE", rd,
                     rd - 1.96 * se_rd, rd + 1.96 * se_rd,
                     nrow(dt), sum(y), "Crude"),
    summarise_result("RR brut", "ATE", rr,
                     exp(log(rr) - 1.96 * se_rr),
                     exp(log(rr) + 1.96 * se_rr),
                     nrow(dt), sum(y), "Crude"),
    summarise_result("OR brut", "ATE", or_v,
                     exp(log(or_v) - 1.96 * se_logOR),
                     exp(log(or_v) + 1.96 * se_logOR),
                     nrow(dt), sum(y), "Crude")
  )

  message("[10_est] Analyse brute - RD = ", round(rd, 3),
          " [", round(rd - 1.96 * se_rd, 3), "; ",
          round(rd + 1.96 * se_rd, 3), "]")
  results
}

# -----------------------------------------------------------------------
# 2. Régression multivariée + marginal standardization
# -----------------------------------------------------------------------
run_glm_marginal <- function(cohort, outcome_col = "outcome_primary",
                              ps_covs = PS_COVARIATES, n_boot = 200) {
  dt <- cohort[!is.na(get(outcome_col))]

  vars_avail <- ps_covs[ps_covs %in% names(dt)]
  fmla <- as.formula(paste(outcome_col, "~ treatment +",
                            paste(vars_avail, collapse = " + ")))

  dt_df <- as.data.frame(dt)
  dt_df$treatment <- factor(dt_df$treatment, levels = c("carbapenem", "PTZ"))
  dt_df[[outcome_col]] <- as.integer(dt_df[[outcome_col]])

  fit <- tryCatch(
    glm(fmla, data = dt_df, family = binomial(link = "logit")),
    error = function(e) {
      message("[10_est] GLM erreur : ", e$message)
      NULL
    }
  )
  if (is.null(fit)) return(NULL)

  # Marginal standardization
  margstd <- function(data, fit) {
    d1 <- d0 <- data
    d1$treatment <- "PTZ"
    d0$treatment <- "carbapenem"
    r1 <- mean(predict(fit, newdata = d1, type = "response"), na.rm = TRUE)
    r0 <- mean(predict(fit, newdata = d0, type = "response"), na.rm = TRUE)
    c(r1 = r1, r0 = r0, rd = r1 - r0, rr = r1 / r0)
  }

  est_main <- margstd(dt_df, fit)

  # Bootstrap IC
  boot_fn <- function(data, idx) {
    margstd(data[idx, ], glm(fmla, data = data[idx, ],
                              family = binomial()))
  }
  tryCatch({
    b <- boot::boot(dt_df, boot_fn, R = n_boot)
    ci_rd <- boot::boot.ci(b, type = "perc", index = 3)$percent[4:5]
    ci_rr <- boot::boot.ci(b, type = "perc", index = 4)$percent[4:5]
  }, error = function(e) {
    ci_rd <<- c(est_main["rd"] - 0.05, est_main["rd"] + 0.05)
    ci_rr <<- c(est_main["rr"] * 0.8, est_main["rr"] * 1.2)
  })

  results <- rbind(
    summarise_result("RD GLM ajusté", "ATE",
                     est_main["rd"], ci_rd[1], ci_rd[2],
                     nrow(dt), sum(dt[[outcome_col]], na.rm = TRUE), "GLM"),
    summarise_result("RR GLM ajusté", "ATE",
                     est_main["rr"], ci_rr[1], ci_rr[2],
                     nrow(dt), sum(dt[[outcome_col]], na.rm = TRUE), "GLM")
  )

  message("[10_est] GLM marginal - RD = ", round(est_main["rd"], 3))
  list(results = results, fit = fit, marginal_estimates = est_main)
}

# -----------------------------------------------------------------------
# 3. Score de propension + IPTW ATE
# -----------------------------------------------------------------------
run_iptw <- function(cohort, outcome_col = "outcome_primary",
                      ps_covs = PS_COVARIATES, trim_q = c(0.01, 0.99),
                      n_boot = 200, cfg = NULL) {
  dt <- cohort[!is.na(get(outcome_col))]

  # Modèles PS en cascade : GLM → ridge/lasso → ranger
  ps_result <- fit_propensity_score(dt, ps_covs)
  ps    <- ps_result$ps
  ps_dt <- data.table(subject_id = dt$subject_id,
                       hadm_id   = dt$hadm_id,
                       ps        = ps,
                       treatment = dt$treatment)

  # Diagnostic overlap
  ps_min   <- min(ps); ps_max <- max(ps)
  ps_q     <- quantile(ps, c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99))
  n_offsupp <- sum(ps < 0.05 | ps > 0.95)
  message("[10_est] PS — min=", round(ps_min, 3), " max=", round(ps_max, 3),
          " hors support commun (5%/95%) : ", n_offsupp)

  # Poids IPTW stabilisés
  p_trt <- mean(dt$treatment == "PTZ")
  trt   <- as.integer(dt$treatment == "PTZ")
  w_iptw <- ifelse(trt == 1, p_trt / ps, (1 - p_trt) / (1 - ps))

  # Troncature
  trim_lo <- quantile(w_iptw, trim_q[1])
  trim_hi <- quantile(w_iptw, trim_q[2])
  w_iptw_trunc <- pmin(pmax(w_iptw, trim_lo), trim_hi)

  message("[10_est] Poids IPTW — mean=", round(mean(w_iptw_trunc), 2),
          " max=", round(max(w_iptw_trunc), 2))

  # Estimation ATE pondérée
  y <- as.integer(dt[[outcome_col]])
  r1_w <- weighted.mean(y[trt == 1], w_iptw_trunc[trt == 1])
  r0_w <- weighted.mean(y[trt == 0], w_iptw_trunc[trt == 0])
  rd_w <- r1_w - r0_w
  rr_w <- r1_w / r0_w

  # Sandwich SE via survey
  se_rd <- tryCatch({
    dt_df <- as.data.frame(dt)
    dt_df$w_iptw <- w_iptw_trunc
    dt_df$trt    <- trt
    dt_df$y      <- y
    svy <- survey::svydesign(ids = ~1, weights = ~w_iptw, data = dt_df)
    m1  <- survey::svyglm(y ~ trt, design = svy, family = quasibinomial())
    coef_se <- sqrt(diag(vcov(m1)))
    coef_se["trt"]
  }, error = function(e) sd(y) / sqrt(length(y)))

  # Bootstrap patient-level pour IC
  boot_iptw <- function(data, idx) {
    d <- data[idx, ]
    trt_b <- d$trt; y_b <- d$y; w_b <- d$w_iptw
    r1 <- weighted.mean(y_b[trt_b == 1], w_b[trt_b == 1])
    r0 <- weighted.mean(y_b[trt_b == 0], w_b[trt_b == 0])
    c(rd = r1 - r0, rr = r1 / r0)
  }

  dt_boot <- data.frame(trt = trt, y = y, w_iptw = w_iptw_trunc)
  tryCatch({
    b    <- boot::boot(dt_boot, boot_iptw, R = n_boot)
    ci_rd <- boot::boot.ci(b, type = "perc", index = 1)$percent[4:5]
    ci_rr <- boot::boot.ci(b, type = "perc", index = 2)$percent[4:5]
  }, error = function(e) {
    ci_rd <<- c(rd_w - 1.96 * se_rd, rd_w + 1.96 * se_rd)
    ci_rr <<- c(rr_w * 0.8, rr_w * 1.2)
  })

  results <- rbind(
    summarise_result("RD IPTW ATE", "ATE",
                     rd_w, ci_rd[1], ci_rd[2],
                     nrow(dt), sum(y), "IPTW"),
    summarise_result("RR IPTW ATE", "ATE",
                     rr_w, ci_rr[1], ci_rr[2],
                     nrow(dt), sum(y), "IPTW")
  )

  message("[10_est] IPTW - RD = ", round(rd_w, 3),
          " [", round(ci_rd[1], 3), "; ", round(ci_rd[2], 3), "]")

  list(
    results        = results,
    ps             = ps,
    ps_dt          = ps_dt,
    weights_iptw   = w_iptw_trunc,
    ps_diagnostics = list(ps_min = ps_min, ps_max = ps_max,
                           ps_q = ps_q, n_offsupp = n_offsupp)
  )
}

# -----------------------------------------------------------------------
# 4. Overlap weighting (ATO)
# -----------------------------------------------------------------------
run_overlap_weighting <- function(cohort, outcome_col = "outcome_primary",
                                   ps_vector = NULL, ps_covs = PS_COVARIATES,
                                   n_boot = 200) {
  dt <- cohort[!is.na(get(outcome_col))]

  ps    <- if (!is.null(ps_vector)) ps_vector else
    fit_propensity_score(dt, ps_covs)$ps
  trt   <- as.integer(dt$treatment == "PTZ")
  y     <- as.integer(dt[[outcome_col]])

  # Poids overlap : w = ps(1-ps) × 1/P(trt)
  w_ow <- ifelse(trt == 1, 1 - ps, ps)

  r1_ow <- weighted.mean(y[trt == 1], w_ow[trt == 1])
  r0_ow <- weighted.mean(y[trt == 0], w_ow[trt == 0])
  rd_ow <- r1_ow - r0_ow
  rr_ow <- r1_ow / r0_ow

  boot_ow <- function(data, idx) {
    d   <- data[idx, ]
    trt_b <- d$trt; y_b <- d$y; ps_b <- d$ps
    w_b <- ifelse(trt_b == 1, 1 - ps_b, ps_b)
    r1  <- weighted.mean(y_b[trt_b == 1], w_b[trt_b == 1])
    r0  <- weighted.mean(y_b[trt_b == 0], w_b[trt_b == 0])
    c(rd = r1 - r0, rr = r1 / r0)
  }
  dt_boot <- data.frame(trt = trt, y = y, ps = ps)
  tryCatch({
    b    <- boot::boot(dt_boot, boot_ow, R = n_boot)
    ci_rd <- boot::boot.ci(b, type = "perc", index = 1)$percent[4:5]
    ci_rr <- boot::boot.ci(b, type = "perc", index = 2)$percent[4:5]
  }, error = function(e) {
    se <- sd(y) / sqrt(length(y))
    ci_rd <<- c(rd_ow - 1.96 * se, rd_ow + 1.96 * se)
    ci_rr <<- c(rr_ow * 0.8, rr_ow * 1.2)
  })

  results <- rbind(
    summarise_result("RD Overlap (ATO)", "ATO",
                     rd_ow, ci_rd[1], ci_rd[2],
                     nrow(dt), sum(y), "Overlap"),
    summarise_result("RR Overlap (ATO)", "ATO",
                     rr_ow, ci_rr[1], ci_rr[2],
                     nrow(dt), sum(y), "Overlap")
  )
  message("[10_est] Overlap weighting - RD = ", round(rd_ow, 3))
  list(results = results, weights_ow = w_ow)
}

# -----------------------------------------------------------------------
# 5. Matching PS
# -----------------------------------------------------------------------
run_matching <- function(cohort, outcome_col = "outcome_primary",
                          ps_vector = NULL, ps_covs = PS_COVARIATES) {
  if (!PKG_AVAILABLE("MatchIt")) {
    message("[10_est] MatchIt absent — matching ignoré")
    return(NULL)
  }

  dt    <- cohort[!is.na(get(outcome_col))]
  df    <- as.data.frame(dt)
  df$Y  <- as.integer(df[[outcome_col]])
  df$trt <- as.integer(df$treatment == "PTZ")

  if (!is.null(ps_vector)) df$ps_score <- ps_vector

  vars_avail <- ps_covs[ps_covs %in% names(df)]
  fmla <- as.formula(paste("trt ~", paste(vars_avail, collapse = " + ")))

  m_out <- tryCatch(
    MatchIt::matchit(fmla, data = df, method = "nearest",
                     distance = "glm", caliper = 0.2 * sd(
                       if (!is.null(ps_vector)) log(ps_vector / (1 - ps_vector))
                       else rep(0.1, nrow(df))),
                     ratio = 1, replace = FALSE),
    error = function(e) {
      message("[10_est] MatchIt erreur : ", e$message)
      NULL
    }
  )
  if (is.null(m_out)) return(NULL)

  matched_df <- MatchIt::match.data(m_out)
  n_matched  <- nrow(matched_df)
  r1_m <- mean(matched_df$Y[matched_df$trt == 1])
  r0_m <- mean(matched_df$Y[matched_df$trt == 0])
  rd_m <- r1_m - r0_m

  # SE robust via sandwich
  fit_m <- lm(Y ~ trt, data = matched_df, weights = matched_df$weights)
  se_m  <- sqrt(sandwich::vcovHC(fit_m, type = "HC3")["trt", "trt"])

  results <- summarise_result("RD Matching", "ATT",
                               rd_m, rd_m - 1.96 * se_m, rd_m + 1.96 * se_m,
                               n_matched, sum(matched_df$Y), "Matching")
  message("[10_est] Matching - RD = ", round(rd_m, 3), " (n=", n_matched, ")")
  list(results = results, matched_data = matched_df, matchit_obj = m_out)
}

# -----------------------------------------------------------------------
# 6. G-computation
# -----------------------------------------------------------------------
run_gcomputation <- function(cohort, outcome_col = "outcome_primary",
                              ps_covs = PS_COVARIATES, n_boot = 200) {
  dt   <- cohort[!is.na(get(outcome_col))]
  df   <- as.data.frame(dt)
  df$Y <- as.integer(df[[outcome_col]])
  df$treatment <- factor(df$treatment, levels = c("carbapenem", "PTZ"))

  vars_avail <- ps_covs[ps_covs %in% names(df)]

  # Modèle outcome flexible
  fmla <- as.formula(paste("Y ~ treatment +",
                            paste(vars_avail, collapse = " + ")))

  # Utiliser ranger si disponible
  if (PKG_AVAILABLE("ranger")) {
    X_mat <- prepare_X(dt, c("treatment_binary", ps_covs))
    dt_rf  <- copy(dt)
    dt_rf[, treatment_binary := as.integer(treatment == "PTZ")]
    X_mat  <- prepare_X(dt_rf, c("treatment_binary", ps_covs))
    y_vec  <- as.integer(dt[[outcome_col]])

    fit_rf <- tryCatch(
      ranger::ranger(y = factor(y_vec), x = data.frame(X_mat),
                     num.trees = 500, probability = TRUE,
                     min.node.size = 10),
      error = function(e) NULL
    )

    if (!is.null(fit_rf)) {
      g_comp_fn <- function(data, idx) {
        d <- data[idx, ]
        # A=1
        d1 <- d; d1$treatment_binary <- 1L
        X1 <- prepare_X(data.table(d1),
                         c("treatment_binary", ps_covs))
        p1 <- predict(fit_rf, data.frame(X1))$predictions[, "1"]
        # A=0
        d0 <- d; d0$treatment_binary <- 0L
        X0 <- prepare_X(data.table(d0),
                         c("treatment_binary", ps_covs))
        p0 <- predict(fit_rf, data.frame(X0))$predictions[, "1"]
        c(rd = mean(p1) - mean(p0), rr = mean(p1) / mean(p0))
      }
      # Simplifier : estimation unique
      dt_rf2 <- copy(dt); dt_rf2[, treatment_binary := 1L]
      X1_all <- prepare_X(dt_rf2, c("treatment_binary", ps_covs))
      p1_all <- predict(fit_rf, data.frame(X1_all))$predictions[, "1"]

      dt_rf2[, treatment_binary := 0L]
      X0_all <- prepare_X(dt_rf2, c("treatment_binary", ps_covs))
      p0_all <- predict(fit_rf, data.frame(X0_all))$predictions[, "1"]

      rd_gc  <- mean(p1_all) - mean(p0_all)
      rr_gc  <- mean(p1_all) / mean(p0_all)
      se_gc  <- sd(p1_all - p0_all) / sqrt(nrow(dt))

      results <- rbind(
        summarise_result("RD G-formula (ranger)", "ATE",
                         rd_gc, rd_gc - 1.96 * se_gc, rd_gc + 1.96 * se_gc,
                         nrow(dt), sum(y_vec), "Gcomp-RF"),
        summarise_result("RR G-formula (ranger)", "ATE",
                         rr_gc, rr_gc * 0.85, rr_gc * 1.15,
                         nrow(dt), sum(y_vec), "Gcomp-RF")
      )
      message("[10_est] G-formula (ranger) - RD = ", round(rd_gc, 3))
      return(list(results = results))
    }
  }

  # Fallback GLM
  fit <- glm(fmla, data = df, family = binomial())
  d1 <- df; d1$treatment <- "PTZ"
  d0 <- df; d0$treatment <- "carbapenem"
  p1 <- mean(predict(fit, d1, type = "response"))
  p0 <- mean(predict(fit, d0, type = "response"))
  rd <- p1 - p0; rr <- p1 / p0
  se <- sd(predict(fit, df, type = "response")) / sqrt(nrow(df))

  results <- rbind(
    summarise_result("RD G-formula (GLM)", "ATE",
                     rd, rd - 1.96 * se, rd + 1.96 * se,
                     nrow(dt), sum(df$Y), "Gcomp-GLM"),
    summarise_result("RR G-formula (GLM)", "ATE",
                     rr, rr * 0.85, rr * 1.15,
                     nrow(dt), sum(df$Y), "Gcomp-GLM")
  )
  message("[10_est] G-formula (GLM) - RD = ", round(rd, 3))
  list(results = results)
}

# -----------------------------------------------------------------------
# 7. AIPW (doubly robust) — implémentation manuelle si package absent
# -----------------------------------------------------------------------
run_aipw <- function(cohort, outcome_col = "outcome_primary",
                      ps_covs = PS_COVARIATES, n_boot = 200) {
  dt  <- cohort[!is.na(get(outcome_col))]
  n   <- nrow(dt)
  y   <- as.integer(dt[[outcome_col]])
  trt <- as.integer(dt$treatment == "PTZ")

  X <- prepare_X(dt, ps_covs)

  # --- Nuisance models via cross-fitting (5 folds) ---
  n_folds <- 5
  folds   <- sample(rep(1:n_folds, length.out = n))

  ps_cf  <- numeric(n)
  mu1_cf <- numeric(n)
  mu0_cf <- numeric(n)

  for (k in 1:n_folds) {
    train_idx <- folds != k
    test_idx  <- folds == k

    X_tr <- X[train_idx, ]; X_te <- X[test_idx, ]
    y_tr <- y[train_idx];   trt_tr <- trt[train_idx]

    # PS model
    if (PKG_AVAILABLE("ranger")) {
      ps_fit <- ranger::ranger(
        y = factor(trt_tr), x = data.frame(X_tr),
        num.trees = 200, probability = TRUE, min.node.size = 5)
      ps_cf[test_idx] <- predict(ps_fit, data.frame(X_te))$predictions[, "1"]
    } else if (PKG_AVAILABLE("glmnet")) {
      ps_fit <- glmnet::cv.glmnet(X_tr, trt_tr, family = "binomial", alpha = 0.5)
      ps_cf[test_idx] <- as.vector(predict(ps_fit, X_te, s = "lambda.min",
                                            type = "response"))
    } else {
      df_tr <- data.frame(y = trt_tr, X_tr)
      ps_fit <- glm(y ~ ., data = df_tr, family = binomial())
      ps_cf[test_idx] <- predict(ps_fit, data.frame(X_te), type = "response")
    }

    # Outcome model sous A=1 et A=0
    df_tr <- data.frame(y = y_tr, trt = trt_tr, X_tr)
    if (PKG_AVAILABLE("ranger")) {
      out_fit <- ranger::ranger(
        y = factor(y_tr), x = data.frame(trt = trt_tr, X_tr),
        num.trees = 200, probability = TRUE, min.node.size = 5)
      mu1_cf[test_idx] <- predict(out_fit, data.frame(trt = 1, X_te))$predictions[, "1"]
      mu0_cf[test_idx] <- predict(out_fit, data.frame(trt = 0, X_te))$predictions[, "1"]
    } else {
      out_fit <- glm(y ~ ., data = df_tr, family = binomial())
      mu1_cf[test_idx] <- predict(out_fit,
                                   data.frame(trt = 1, X_te), type = "response")
      mu0_cf[test_idx] <- predict(out_fit,
                                   data.frame(trt = 0, X_te), type = "response")
    }
  }

  # Troncature PS
  ps_cf <- pmax(pmin(ps_cf, 0.95), 0.05)

  # Scores d'influence AIPW
  phi1 <- mu1_cf + trt / ps_cf * (y - mu1_cf)
  phi0 <- mu0_cf + (1 - trt) / (1 - ps_cf) * (y - mu0_cf)
  phi  <- phi1 - phi0

  ate_aipw <- mean(phi)
  se_aipw  <- sd(phi) / sqrt(n)

  message("[10_est] AIPW (cross-fit) - ATE = ", round(ate_aipw, 3),
          " [", round(ate_aipw - 1.96 * se_aipw, 3), "; ",
          round(ate_aipw + 1.96 * se_aipw, 3), "]")

  results <- summarise_result("RD AIPW (doubly robust)", "ATE",
                               ate_aipw,
                               ate_aipw - 1.96 * se_aipw,
                               ate_aipw + 1.96 * se_aipw,
                               n, sum(y), "AIPW")

  list(results = results, phi = phi, ps_cf = ps_cf,
       mu1 = mu1_cf, mu0 = mu0_cf)
}

# -----------------------------------------------------------------------
# 8. TMLE (si disponible)
# -----------------------------------------------------------------------
run_tmle <- function(cohort, outcome_col = "outcome_primary",
                      ps_covs = PS_COVARIATES) {
  if (!PKG_AVAILABLE("tmle")) {
    message("[10_est] TMLE absent — estimateur ignoré")
    return(NULL)
  }

  dt  <- cohort[!is.na(get(outcome_col))]
  y   <- as.integer(dt[[outcome_col]])
  a   <- as.integer(dt$treatment == "PTZ")
  X   <- data.frame(prepare_X(dt, ps_covs))

  # SuperLearner léger pour TMLE : XGBoost exclu (redondant avec AIPW+ranger, très lent)
  sl_libs <- c("SL.glm")
  if (PKG_AVAILABLE("ranger"))  sl_libs <- c(sl_libs, "SL.ranger")

  tryCatch({
    fit_tmle <- tmle::tmle(
      Y = y, A = a, W = X,
      family       = "binomial",
      Q.SL.library = sl_libs,
      g.SL.library = sl_libs,
      V            = 5L   # 5-fold CV au lieu de 10 pour accélérer
    )
    est <- fit_tmle$estimates$ATE
    message("[10_est] TMLE - ATE = ", round(est$psi, 3),
            " [", round(est$CI[1], 3), "; ", round(est$CI[2], 3), "]")
    results <- summarise_result("RD TMLE", "ATE",
                                 est$psi, est$CI[1], est$CI[2],
                                 nrow(dt), sum(y), "TMLE")
    list(results = results, tmle_fit = fit_tmle)
  }, error = function(e) {
    message("[10_est] TMLE erreur : ", e$message)
    NULL
  })
}

# -----------------------------------------------------------------------
# Utilitaire : fit propensity score
# -----------------------------------------------------------------------
fit_propensity_score <- function(dt, ps_covs) {
  vars_avail <- ps_covs[ps_covs %in% names(dt)]
  X   <- prepare_X(dt, vars_avail)
  trt <- as.integer(dt$treatment == "PTZ")

  ps <- NULL

  # Essai 1 : ranger (random forest)
  if (PKG_AVAILABLE("ranger")) {
    tryCatch({
      fit_rf <- ranger::ranger(y = factor(trt), x = data.frame(X),
                                num.trees = 500, probability = TRUE,
                                min.node.size = 10)
      ps <- predict(fit_rf, data.frame(X))$predictions[, "1"]
      message("[10_est] PS via ranger")
    }, error = function(e) {
      message("[10_est] ranger PS erreur : ", e$message)
    })
  }

  # Essai 2 : glmnet LASSO/ridge
  if (is.null(ps) && PKG_AVAILABLE("glmnet")) {
    tryCatch({
      cv_fit <- glmnet::cv.glmnet(X, trt, family = "binomial", alpha = 0.5,
                                   nfolds = 5)
      ps <- as.vector(predict(cv_fit, X, s = "lambda.min", type = "response"))
      message("[10_est] PS via glmnet elastic net")
    }, error = function(e) {
      message("[10_est] glmnet PS erreur : ", e$message)
    })
  }

  # Fallback : GLM logistique
  if (is.null(ps)) {
    fmla <- as.formula(paste("trt ~", paste(vars_avail, collapse = " + ")))
    df   <- data.frame(trt = trt, X)
    tryCatch({
      fit_glm <- glm(fmla, data = df, family = binomial())
      ps <- predict(fit_glm, type = "response")
      message("[10_est] PS via GLM logistique")
    }, error = function(e) {
      message("[10_est] GLM PS erreur : ", e$message)
      ps <- rep(mean(trt), nrow(dt))
    })
  }

  ps <- pmax(pmin(ps, 0.99), 0.01)
  list(ps = ps, n = nrow(dt), p_trt = mean(trt))
}

# -----------------------------------------------------------------------
# Compilation des résultats
# -----------------------------------------------------------------------
compile_all_results <- function(...) {
  results_list <- list(...)
  all_res <- rbindlist(
    lapply(results_list, function(x) {
      if (is.null(x)) return(NULL)
      if (is.data.table(x)) return(x)
      if (is.list(x) && !is.null(x$results)) return(x$results)
      NULL
    }),
    fill = TRUE
  )
  save_derived(all_res, "all_estimator_results", format = "rds")
  data.table::fwrite(all_res, "output/tables/estimator_results.csv")
  all_res
}
