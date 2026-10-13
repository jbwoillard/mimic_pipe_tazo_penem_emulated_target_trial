# =============================================================================
# 80_robustness.R — Phase A v3.3 (docs/v3.3_plan.md) : complexité du PS,
#   PS parcimonieux, overlap weighting augmenté, AIPW hors pli, règle des
#   20 événements pour le RR, distribution bootstrap du primaire
# =============================================================================

# -----------------------------------------------------------------------
# Variables dérivées utilisées par le PS parcimonieux (ajoutées aux jeux imputés)
# -----------------------------------------------------------------------
add_parsimonious_vars_v33 <- function(d) {
  d <- copy(d)
  d[, other_nephrotoxin_any := as.integer(direct_nephrotoxin_pre_t0 == 1L | calcineurin_inhibitor_pre_t0 == 1L |
                                            aminoglycoside_pre_t0 == 1L | nsaid_pre_t0 == 1L)]
  d[, period_2017plus := as.integer(anchor_year_group %in% c("2017 - 2019", "2020 - 2022"))]
  d[, care_location_3 := as.character(care_location_t0)]   # ICU / ED / ward (catégorielle, 2 paramètres)
  setkey(d, hadm_id)
  d
}

#' PS parcimonieux (A1-a) : les 8 blocs de l'ensemble d'ajustement minimal n° 1
#' du DAG, 1 à 4 variables pré-spécifiées par bloc (figé avant tout calcul).
PS_PARSIMONIOUS_V33 <- list(
  BaselineRenal   = c("baseline_creat_min7d", "bun_baseline"),
  CreatTrajectory = c("creat_ratio_last_min7d", "creat_slope_48h"),
  Severity        = c("sepsis_composite_v3", "vasopressor_pre_t0", "vent_pre_t0", "lactate_max"),
  CareUnit        = c("care_location_3"),
  Comorbidity     = c("age", "ckd_priorhosp", "charlson_score_priorhosp"),
  Nephrotoxins    = c("loop_diuretic_pre_t0", "other_nephrotoxin_any"),
  Vancomycin      = c("vanco_6h"),
  Period          = c("period_2017plus"))
PS_PARSIMONIOUS_DF_V33 <- 2L

# -----------------------------------------------------------------------
# A1 — complexité du PS et surajustement
# -----------------------------------------------------------------------
ps_complexity_v33 <- function(d, covs, method, df, population, spec_label, seed = SEED_V3) {
  X <- build_ps_design_v31(d, covs, df = df); a <- as.integer(d$treatment == "PTZ")
  rk <- qr(cbind(1, X))$rank
  f <- fit_ps_v31(X, a, method, seed = seed)
  data.table(population = population, specification = spec_label, ps_method = method, spline_df = df,
             n_covariates = length(covs), n_parameters = ncol(X), n_collinear = ncol(X) + 1L - rk,
             n = nrow(d), n_comparator = sum(a == 0), comparator_per_parameter = sum(a == 0) / ncol(X),
             auc_in_sample = as.numeric(pROC::auc(pROC::roc(a, f$ps, quiet = TRUE))),
             auc_out_of_fold = .oof_auc(X, a, method, seed = seed),
             pct_ps_outside_0.1_0.9 = 100 * mean(f$ps < 0.1 | f$ps > 0.9))
}

# -----------------------------------------------------------------------
# A1-b — overlap weighting augmenté (Mao, Li & Greene 2019)
# -----------------------------------------------------------------------
.ridge_fit_v33 <- function(X, y, lambda = NULL, seed = SEED_V3) {
  if (is.null(lambda)) { set.seed(seed)
    lambda <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = 0, nfolds = 5)$lambda.min }
  list(fit = glmnet::glmnet(X, y, family = "binomial", alpha = 0, lambda = lambda), lambda = lambda)
}
.ridge_pred_v33 <- function(f, X) as.vector(stats::predict(f$fit, X, type = "response", s = f$lambda))

#' Estimateur ATO augmenté sur un jeu (PS ridge, outcome ridge par bras)
#' @return mu1, mu0 (risques augmentés), h, poids
.ato_aug_point_v33 <- function(X, a, y, lam = NULL) {
  fe <- fit_ps_v31(X, a, "ridge", lambda = lam$ps)
  e <- fe$ps
  f1 <- .ridge_fit_v33(X[a == 1, , drop = FALSE], y[a == 1], lam$m1); f0 <- .ridge_fit_v33(X[a == 0, , drop = FALSE], y[a == 0], lam$m0)
  m1 <- .ridge_pred_v33(f1, X); m0 <- .ridge_pred_v33(f0, X)
  h <- e * (1 - e); w <- ifelse(a == 1, 1 - e, e)
  mu1 <- sum(h * m1) / sum(h) + sum((w * (y - m1))[a == 1]) / sum(w[a == 1])
  mu0 <- sum(h * m0) / sum(h) + sum((w * (y - m0))[a == 0]) / sum(w[a == 0])
  list(mu1 = mu1, mu0 = mu0, e = e, m1 = m1, m0 = m0, w = w,
       lam = list(ps = fe$lambda, m1 = f1$lambda, m0 = f0$lambda))
}

#' Interface pool_estimator_mi_v31 : bootstrap (B) avec lambdas fixés à leur
#' valeur de validation croisée sur l'échantillon complet
ato_augmented_v33 <- function(dt, ps_covs, outcome_col, seed = SEED_V3, B = 200) {
  X <- build_ps_design_v31(dt, ps_covs); a <- as.integer(dt$treatment == "PTZ"); y <- as.integer(dt[[outcome_col]])
  set.seed(seed); pt <- .ato_aug_point_v33(X, a, y)
  set.seed(seed)
  bt <- t(vapply(seq_len(B), function(b) tryCatch({ i <- sample.int(length(a), replace = TRUE)
    r <- .ato_aug_point_v33(build_ps_design_v31(dt[i], ps_covs), a[i], y[i], pt$lam); c(r$mu1, r$mu0) },
    error = function(e) c(NA_real_, NA_real_)), numeric(2)))
  ess <- function(s) sum(pt$w[s])^2 / sum(pt$w[s]^2)
  list(n = length(a), mu1 = pt$mu1, mu0 = pt$mu0, rd = pt$mu1 - pt$mu0, rr = pt$mu1 / pt$mu0,
       se_rd = stats::sd(bt[, 1] - bt[, 2], na.rm = TRUE), se_logrr = stats::sd(log(bt[, 1] / bt[, 2]), na.rm = TRUE),
       ess_ptz = ess(a == 1), ess_comp = ess(a == 0), n_boot_failed = sum(is.na(bt[, 1])), point = pt)
}

#' Contrôle de l'implémentation contre PSweight (même PS et mêmes prédictions d'outcome)
check_augmented_vs_psweight_v33 <- function(dt, ps_covs, outcome_col) {
  X <- build_ps_design_v31(dt, ps_covs); a <- as.integer(dt$treatment == "PTZ"); y <- as.integer(dt[[outcome_col]])
  set.seed(SEED_V3); pt <- .ato_aug_point_v33(X, a, y)
  df <- data.frame(a = a, y = y, x1 = X[, 1])
  ps_mat <- cbind(`0` = 1 - pt$e, `1` = pt$e)
  ps <- PSweight::PSweight(ps.formula = a ~ x1, ps.estimate = ps_mat, yname = "y", data = df, weight = "overlap",
                           augmentation = TRUE, out.estimate = cbind(`0` = pt$m0, `1` = pt$m1), family = "binomial")
  mu <- ps$muhat
  data.table(estimator = c("own implementation", "PSweight"), mu1 = c(pt$mu1, mu[["1"]]), mu0 = c(pt$mu0, mu[["0"]]),
             rd = c(pt$mu1 - pt$mu0, mu[["1"]] - mu[["0"]]))
}

# -----------------------------------------------------------------------
# A1-c — AIPW hors pli (ATE), nuisances ridge, 5 plis
# -----------------------------------------------------------------------
aipw_crossfit_v33 <- function(dt, ps_covs, outcome_col, seed = SEED_V3, K = 5) {
  X <- build_ps_design_v31(dt, ps_covs); a <- as.integer(dt$treatment == "PTZ"); y <- as.integer(dt[[outcome_col]]); n <- length(a)
  set.seed(seed); folds <- sample(rep(seq_len(K), length.out = n))
  e <- m1 <- m0 <- numeric(n)
  for (k in seq_len(K)) {
    tr <- folds != k; te <- folds == k
    fe <- .ridge_fit_v33(X[tr, , drop = FALSE], a[tr], seed = seed + k); e[te] <- .ridge_pred_v33(fe, X[te, , drop = FALSE])
    f1 <- .ridge_fit_v33(X[tr & a == 1, , drop = FALSE], y[tr & a == 1], seed = seed + k); m1[te] <- .ridge_pred_v33(f1, X[te, , drop = FALSE])
    f0 <- .ridge_fit_v33(X[tr & a == 0, , drop = FALSE], y[tr & a == 0], seed = seed + k); m0[te] <- .ridge_pred_v33(f0, X[te, , drop = FALSE])
  }
  auc_oof <- as.numeric(pROC::auc(pROC::roc(a, e, quiet = TRUE)))
  e <- pmin(pmax(e, 0.01), 0.99)
  phi1 <- m1 + a * (y - m1) / e; phi0 <- m0 + (1 - a) * (y - m0) / (1 - e)
  mu1 <- mean(phi1); mu0 <- mean(phi0)
  list(n = n, mu1 = mu1, mu0 = mu0, rd = mu1 - mu0, rr = mu1 / mu0,
       se_rd = stats::sd((phi1 - mu1) - (phi0 - mu0)) / sqrt(n),
       se_logrr = stats::sd((phi1 - mu1) / mu1 - (phi0 - mu0) / mu0) / sqrt(n),
       auc_oof = auc_oof, pct_truncated = 100 * mean(e %in% c(0.01, 0.99)))
}

# -----------------------------------------------------------------------
# A3 — règle des 20 événements et distribution bootstrap
# -----------------------------------------------------------------------
MIN_EVENTS_RR_V33 <- 20L

#' RR non rapporté si le bras comparateur compte moins de 20 événements bruts.
#' La valeur calculée est conservée dans rr_raw* (traçabilité) ; rr* = NA.
apply_rr_rule_v33 <- function(rmt, min_events = MIN_EVENTS_RR_V33) {
  rmt <- copy(rmt)
  rmt[, `:=`(rr_raw = rr, rr_raw_ci_low = rr_ci_low, rr_raw_ci_high = rr_ci_high)]
  is_or <- grepl("^ordinal_OR", rmt$estimand)
  rmt[, rr_status := fcase(is_or, "odds ratio (ordinal outcome), rule not applicable",
                           is.na(events_comp), "events not available",
                           events_comp < min_events, sprintf("not estimated (<%d events in the comparator arm)", min_events),
                           default = "reported")]
  rmt[grepl("^not estimated", rr_status), `:=`(rr = NA_real_, rr_ci_low = NA_real_, rr_ci_high = NA_real_)]
  rmt[rr_status == "reported" & (!is.finite(rr_ci_low) | !is.finite(rr_ci_high)),
      `:=`(rr = NA_real_, rr_ci_low = NA_real_, rr_ci_high = NA_real_, rr_status = "not estimable (non-finite bootstrap interval)")]
  rmt
}

#' Réplications bootstrap de toutes les imputations, centrées sur l'estimation
#' de leur imputation puis recentrées sur l'estimation poolée
boot_long_v33 <- function(est) {
  rbindlist(lapply(seq_along(est$boot_by_imputation), function(i) { b <- est$boot_by_imputation[[i]]
    data.table(imputation = i, rd = b$boot_rd - b$rd + est$rd, log_rr = log(b$boot_rr) - log(b$rr) + log(est$rr)) }))[is.finite(rd) & is.finite(log_rr)]
}

fig_boot_dist_v33 <- function(bl) {
  long <- rbind(bl[, .(stat = "Risk difference, percentage points", value = 100 * rd)], bl[, .(stat = "log risk ratio", value = log_rr)])
  long[, stat := factor(stat, levels = c("Risk difference, percentage points", "log risk ratio"))]
  p1 <- ggplot2::ggplot(long, ggplot2::aes(value)) + ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 60, fill = "grey75", colour = "white") +
    ggplot2::geom_density(linewidth = 0.6) + ggplot2::facet_wrap(~stat, scales = "free") + ggplot2::labs(x = NULL, y = "Density") + theme_fig_v32()
  p2 <- ggplot2::ggplot(long, ggplot2::aes(sample = value)) + ggplot2::stat_qq(size = 0.4, alpha = 0.4) + ggplot2::stat_qq_line(colour = "grey40") +
    ggplot2::facet_wrap(~stat, scales = "free") + ggplot2::labs(x = "Normal quantiles", y = "Bootstrap quantiles") + theme_fig_v32()
  patchwork::wrap_plots(p1, p2, ncol = 1)
}
