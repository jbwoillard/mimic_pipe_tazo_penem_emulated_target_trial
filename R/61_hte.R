# =============================================================================
# 61_hte.R — Hétérogénéité de l'effet (docs/hte_plan.md)
#   cross-fitting 5 plis, risque de baseline emboîté, 5 causal forests,
#   scores AIPW hors pli, comparateurs (R-learner xgboost, T-learner),
#   métriques (RATE, calibration, BLP, interaction), validation par période
# =============================================================================
# Règle de non-fuite : pour tout patient, τ̂, le risque de baseline et les
# scores Γ proviennent de modèles ajustés SANS ce patient. Les bases de
# splines des matrices de design sont construites sur l'ensemble train + test
# (covariables seules, ni traitement ni outcome) pour que train et test
# partagent la même base ; c'est la seule information commune.

HTE_PREDICTORS_V31 <- c("age", "gender", "baseline_creat_min7d", "bun_baseline", "ckd_priorhosp",
                        "sepsis_composite_v3", "vasopressor_pre_t0", "vent_pre_t0", "icu_before_t0",
                        "vanco_6h", "loop_diuretic_pre_t0", "other_nephrotoxin_any", "time_admission_to_T0_h",
                        "baseline_aki_risk")

hte_cfg_v31 <- function(fast = FALSE) list(
  K = 5L, n_forests = if (fast) 1L else 5L, num_trees = if (fast) 500L else 4000L,
  tune = if (fast) "none" else "all", B_boot = if (fast) 50L else 1000L, seed = SEED_V3)

#' Matrice numérique des prédicteurs du CATE (gender -> indicatrice F)
hte_X14 <- function(d) {
  X <- as.matrix(data.frame(lapply(HTE_PREDICTORS_V31, function(v)
    if (v == "gender") as.numeric(toupper(d[[v]]) == "F") else as.numeric(d[[v]]))))
  colnames(X) <- HTE_PREDICTORS_V31
  X
}

.make_folds <- function(W, Y, K, seed) {
  set.seed(seed)
  f <- integer(length(W))
  for (s in split(seq_along(W), paste(W, Y))) f[s] <- sample(rep(seq_len(K), length.out = length(s)))
  f
}

# Ridge logistique sur un design pré-construit : ajuste sur `tr`, prédit sur `te`
.ridge_fit_predict <- function(X, y, tr, te, seed) {
  set.seed(seed)
  cv <- glmnet::cv.glmnet(X[tr, , drop = FALSE], y[tr], family = "binomial", alpha = 0, nfolds = 5)
  as.vector(stats::predict(cv, X[te, , drop = FALSE], s = "lambda.min", type = "response"))
}

#' Prédiction hors pli (validation croisée interne) d'un ridge sur les indices `idx`
.ridge_crossfit <- function(X, y, idx, K, seed) {
  set.seed(seed)
  f <- sample(rep(seq_len(K), length.out = length(idx)))
  out <- numeric(length(idx))
  for (k in seq_len(K)) out[f == k] <- .ridge_fit_predict(X, y, idx[f != k], idx[f == k], seed + k)
  out
}

#' Entraîne la chaîne CATE sur `tr` et prédit sur `te`
#' @return list(tau, var, risk_te, risk_tr, forests_vi, tau_tr_oob, models...)
.train_predict_cate <- function(d, Xps, Xmu, tr, te, cfg, seed) {
  Y <- d$outcome_primary_v31; W <- as.integer(d$treatment == "PTZ")
  # 1. risque de baseline (sans traitement), emboîté
  risk_tr <- .ridge_crossfit(Xps, Y, tr, cfg$K, seed + 11)
  risk_te <- .ridge_fit_predict(Xps, Y, tr, te, seed + 12)
  d_tr <- copy(d[tr]); d_tr[, baseline_aki_risk := risk_tr]
  d_te <- copy(d[te]); d_te[, baseline_aki_risk := risk_te]
  X14_tr <- hte_X14(d_tr); X14_te <- hte_X14(d_te)
  # 2. nuisances de la forêt, hors pli à l'intérieur de l'entraînement
  W_hat <- pmin(pmax(.ridge_crossfit(Xps, W, tr, cfg$K, seed + 21), 0.01), 0.99)
  rf_y <- grf::regression_forest(Xmu[tr, , drop = FALSE], Y[tr], num.trees = 2000, seed = seed + 22)
  Y_hat <- stats::predict(rf_y)$predictions            # OOB
  # 3. forêts causales (moyenne de n_forests graines)
  fits <- lapply(seq_len(cfg$n_forests), function(j)
    grf::causal_forest(X14_tr, Y[tr], W[tr], W.hat = W_hat, Y.hat = Y_hat, num.trees = cfg$num_trees,
                       honesty = TRUE, tune.parameters = cfg$tune, seed = seed + 100 * j))
  pr <- lapply(fits, function(f) stats::predict(f, X14_te, estimate.variance = TRUE))
  oob <- sapply(fits, function(f) stats::predict(f)$predictions)
  list(tau = rowMeans(sapply(pr, `[[`, "predictions")), var = rowMeans(sapply(pr, `[[`, "variance.estimates")),
       tau_tr_oob = if (is.matrix(oob)) rowMeans(oob) else oob,
       risk_te = risk_te, risk_tr = risk_tr, W_hat_tr = W_hat, Y_hat_tr = Y_hat,
       vi = rowMeans(sapply(fits, grf::variable_importance)), X14_tr = X14_tr, X14_te = X14_te)
}

#' Scores AIPW sur `te` avec des nuisances ajustées sur `tr`
.aipw_scores <- function(d, Xps, Xmu, tr, te, seed) {
  Y <- d$outcome_primary_v31; W <- as.integer(d$treatment == "PTZ")
  e <- pmin(pmax(.ridge_fit_predict(Xps, W, tr, te, seed + 31), 0.01), 0.99)
  mu <- lapply(0:1, function(a) {
    s <- tr[W[tr] == a]
    f <- grf::regression_forest(Xmu[s, , drop = FALSE], Y[s], num.trees = 2000, seed = seed + 32 + a)
    stats::predict(f, Xmu[te, , drop = FALSE])$predictions
  })
  mu0 <- mu[[1]]; mu1 <- mu[[2]]; w <- W[te]; y <- Y[te]
  phi1 <- mu1 + w * (y - mu1) / e; phi0 <- mu0 + (1 - w) * (y - mu0) / (1 - e)
  data.table(row = te, e = e, mu1 = mu1, mu0 = mu0, phi1 = phi1, phi0 = phi0, gamma = phi1 - phi0)
}

#' R-learner (Robinson + xgboost) et T-learner (GLM + splines) sur `tr` -> `te`
.comparators <- function(d, ctp, tr, te, seed) {
  Y <- d$outcome_primary_v31; W <- as.integer(d$treatment == "PTZ")
  e <- ctp$W_hat_tr; m <- ctp$Y_hat_tr; wt <- W[tr]
  ps <- (Y[tr] - m) / (wt - e); ww <- (wt - e)^2
  dtr <- xgboost::xgb.DMatrix(ctp$X14_tr, label = ps, weight = ww)
  params <- list(objective = "reg:squarederror", eta = 0.05, max_depth = 3, subsample = 0.8, colsample_bytree = 0.8)
  set.seed(seed + 41)
  cv <- xgboost::xgb.cv(params = params, data = dtr, nrounds = 1000, nfold = 5, early_stopping_rounds = 20, verbose = 0)
  nr <- if (!is.null(cv$best_iteration)) cv$best_iteration else which.min(cv$evaluation_log$test_rmse_mean)
  bst <- xgboost::xgb.train(params = params, data = dtr, nrounds = max(nr, 1), verbose = 0)
  tau_r <- stats::predict(bst, xgboost::xgb.DMatrix(ctp$X14_te))
  # T-learner : base de splines commune train + test
  d14 <- rbind(data.table(ctp$X14_tr), data.table(ctp$X14_te))
  cont <- intersect(c("age", "baseline_creat_min7d", "bun_baseline", "time_admission_to_T0_h", "baseline_aki_risk"), names(d14))
  B <- do.call(cbind, c(lapply(cont, function(v) splines::ns(d14[[v]], df = 3)),
                        list(as.matrix(d14[, setdiff(HTE_PREDICTORS_V31, cont), with = FALSE]))))
  ntr <- length(tr); Btr <- B[seq_len(ntr), ]; Bte <- B[-seq_len(ntr), ]
  pr <- sapply(0:1, function(a) {
    s <- which(wt == a)
    f <- suppressWarnings(stats::glm.fit(cbind(1, Btr[s, ]), Y[tr][s], family = stats::binomial()))
    cf <- f$coefficients; cf[is.na(cf)] <- 0
    as.vector(stats::plogis(cbind(1, Bte) %*% cf))
  })
  list(tau_rlearner = tau_r, tau_tlearner = pr[, 2] - pr[, 1])
}

#' Schéma 2a : cross-fitting complet. Renvoie un data.table une ligne par patient.
run_hte_crossfit_v31 <- function(d, ps_covs, cfg) {
  n <- nrow(d); Y <- d$outcome_primary_v31; W <- as.integer(d$treatment == "PTZ")
  Xps <- build_ps_design_v31(d, ps_covs)                      # même spécification que le PS primaire
  Xmu <- build_ps_design_v31(d, ps_covs, df = 1)              # forêts : covariables brutes + indicatrices
  folds <- .make_folds(W, Y, cfg$K, cfg$seed)
  res <- list(); vi <- list(); leak <- list()
  for (k in seq_len(cfg$K)) {
    t0 <- Sys.time()
    tr <- which(folds != k); te <- which(folds == k)
    ctp <- .train_predict_cate(d, Xps, Xmu, tr, te, cfg, cfg$seed + 1000 * k)
    sc <- .aipw_scores(d, Xps, Xmu, tr, te, cfg$seed + 1000 * k)
    cmp <- .comparators(d, ctp, tr, te, cfg$seed + 1000 * k)
    res[[k]] <- data.table(row = te, fold = k, tau = ctp$tau, tau_var = ctp$var, baseline_aki_risk = ctp$risk_te,
                           tau_rlearner = cmp$tau_rlearner, tau_tlearner = cmp$tau_tlearner)[sc, on = "row"]
    vi[[k]] <- ctp$vi
    leak[[k]] <- data.table(fold = k, n_train = length(tr), n_test = length(te),
                            overlap_ids = length(intersect(d$hadm_id[tr], d$hadm_id[te])))
    message(sprintf("[61_hte_v31] pli %d/%d : n_test=%d, τ̂ moyen=%+.2f%% (%.1f min)", k, cfg$K, length(te),
                    100 * mean(ctp$tau), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  out <- rbindlist(res)[order(row)]
  out[, `:=`(hadm_id = d$hadm_id[row], treatment = d$treatment[row], y = Y[row])]
  list(per_patient = out, vi = setNames(Reduce(`+`, vi) / cfg$K, HTE_PREDICTORS_V31), leak = rbindlist(leak),
       folds = folds)
}

#' Tiers : seuils globaux sur τ̂ hors pli (D3)
assign_tertiles_v31 <- function(tau, thresholds = NULL) {
  if (is.null(thresholds)) thresholds <- stats::quantile(tau, c(1/3, 2/3), names = FALSE)
  t <- fifelse(tau <= thresholds[1], "lower", fifelse(tau <= thresholds[2], "middle", "upper"))
  list(tertile = factor(t, levels = c("lower", "middle", "upper")), thresholds = thresholds)
}

#' Effets DR par groupe : moyenne des Γ, IC normal, RR par delta-méthode,
#' AKI observée par bras, effectifs
dr_by_group_v31 <- function(pp, group_col) {
  pp[, {
    n <- .N; m1 <- mean(phi1); m0 <- mean(phi0); g <- phi1 - phi0
    se <- stats::sd(g) / sqrt(n)
    vl <- stats::var(phi1) / (n * m1^2) + stats::var(phi0) / (n * m0^2) - 2 * stats::cov(phi1, phi0) / (n * m1 * m0)
    .(n = n, n_ptz = sum(treatment == "PTZ"), n_carb = sum(treatment != "PTZ"),
      events_ptz = sum(y[treatment == "PTZ"]), events_carb = sum(y[treatment != "PTZ"]),
      risk_obs_ptz = mean(y[treatment == "PTZ"]), risk_obs_carb = mean(y[treatment != "PTZ"]),
      mean_pred_tau = mean(tau), dr_risk_ptz = m1, dr_risk_carb = m0,
      rd = mean(g), rd_se = se, rd_lo = mean(g) - 1.96 * se, rd_hi = mean(g) + 1.96 * se,
      rr = m1 / m0, rr_lo = exp(log(m1 / m0) - 1.96 * sqrt(vl)), rr_hi = exp(log(m1 / m0) + 1.96 * sqrt(vl)))
  }, by = group_col]
}

#' Bootstrap des patients (τ̂ et Γ fixés) pour l'IC des effets par tiers
dr_by_group_boot_v31 <- function(pp, group_col, B, seed) {
  set.seed(seed)
  bt <- replicate(B, { s <- pp[, .SD[sample.int(.N, replace = TRUE)], by = group_col]
    s[, .(rd = mean(gamma)), by = group_col]$rd })
  data.table(group = levels(pp[[group_col]]), boot_lo = apply(bt, 1, stats::quantile, 0.025),
             boot_hi = apply(bt, 1, stats::quantile, 0.975))
}

#' ATO dans chaque tiers (overlap weighting avec le PS hors pli), IC bootstrap
ato_by_group_v31 <- function(pp, group_col, B, seed) {
  f <- function(s) { w <- ifelse(s$treatment == "PTZ", 1 - s$e, s$e); a <- s$treatment == "PTZ"
    sum(w[a] * s$y[a]) / sum(w[a]) - sum(w[!a] * s$y[!a]) / sum(w[!a]) }
  set.seed(seed)
  pp[, { est <- f(.SD); bt <- replicate(B, f(.SD[sample.int(.N, replace = TRUE)]))
         .(ato_rd = est, ato_lo = stats::quantile(bt, 0.025), ato_hi = stats::quantile(bt, 0.975)) }, by = group_col]
}

#' Test d'interaction principal : Wald à 2 ddl d'égalité des 3 effets DR
wald_equality_v31 <- function(tt) {
  v <- tt$rd_se^2; th <- sum(tt$rd / v) / sum(1 / v)
  chi <- sum((tt$rd - th)^2 / v); df <- nrow(tt) - 1
  list(chi2 = chi, df = df, p = stats::pchisq(chi, df, lower.tail = FALSE))
}

#' Test d'interaction complémentaire : logistique pondérée OW, traitement x tiers
svy_interaction_v31 <- function(pp) {
  dd <- data.frame(y = pp$y, ptz = as.integer(pp$treatment == "PTZ"), tertile = pp$tertile,
                   w = ifelse(pp$treatment == "PTZ", 1 - pp$e, pp$e))
  des <- survey::svydesign(ids = ~1, weights = ~w, data = dd)
  fit <- survey::svyglm(y ~ ptz * tertile, design = des, family = stats::quasibinomial())
  wald <- survey::regTermTest(fit, ~ptz:tertile, method = "Wald")
  lrt <- tryCatch(survey::regTermTest(fit, ~ptz:tertile, method = "LRT"), error = function(e) NULL)
  list(p_wald = as.numeric(wald$p), p_lrt = if (is.null(lrt)) NA_real_ else as.numeric(lrt$p))
}

#' Métriques hors échantillon : RATE (AUTOC, Qini), calibration, BLP
hte_metrics_v31 <- function(pp, seed) {
  set.seed(seed)
  rate <- lapply(c("AUTOC", "QINI"), function(tg) {
    r <- grf::rank_average_treatment_effect.fit(pp$gamma, pp$tau, target = tg, R = 500)
    data.table(metric = paste0("RATE_", tg), estimate = r$estimate, se = r$std.err,
               lo = r$estimate - 1.96 * r$std.err, hi = r$estimate + 1.96 * r$std.err,
               p = 2 * stats::pnorm(-abs(r$estimate / r$std.err)))
  })
  cal <- stats::lm(gamma ~ 0 + mean_term + diff_term, data = data.frame(gamma = pp$gamma, mean_term = mean(pp$tau),
                                                                       diff_term = pp$tau - mean(pp$tau)))
  ct <- lmtest::coeftest(cal, vcov. = sandwich::vcovHC(cal, type = "HC3"))
  calib <- data.table(metric = c("calibration_mean_forest_prediction", "calibration_differential_forest_prediction"),
                      estimate = ct[, 1], se = ct[, 2], lo = ct[, 1] - 1.96 * ct[, 2], hi = ct[, 1] + 1.96 * ct[, 2],
                      p = c(2 * stats::pnorm(-abs((ct[1, 1] - 1) / ct[1, 2])),        # H0 : coefficient = 1 (moyenne)
                            stats::pnorm(-ct[2, 1] / ct[2, 2])))                      # H0 : <= 0 (unilatéral, comme grf)
  blp <- function(vars, tag) {
    f <- stats::lm(stats::as.formula(paste("gamma ~", paste(vars, collapse = " + "))), data = pp)
    ct <- lmtest::coeftest(f, vcov. = sandwich::vcovHC(f, type = "HC3"))
    data.table(metric = paste0(tag, ":", rownames(ct)), estimate = ct[, 1], se = ct[, 2],
               lo = ct[, 1] - 1.96 * ct[, 2], hi = ct[, 1] + 1.96 * ct[, 2], p = ct[, 4])
  }
  rbindlist(c(rate, list(calib, blp(c("vanco_6h", "baseline_aki_risk"), "BLP_prespecified"),
                         blp(setdiff(HTE_PREDICTORS_V31, "gender"), "BLP_all"))))
}

#' grf::test_calibration et BLP sur une forêt ajustée sur toute la cohorte (OOB,
#' étiqueté « in-sample OOB » : complément, pas l'évaluation principale)
hte_full_forest_checks_v31 <- function(d, pp, ps_covs, cfg) {
  Xps <- build_ps_design_v31(d, ps_covs); W <- as.integer(d$treatment == "PTZ"); Y <- d$outcome_primary_v31
  dd <- copy(d); dd[, baseline_aki_risk := pp$baseline_aki_risk[match(d$hadm_id, pp$hadm_id)]]
  W_hat <- pmin(pmax(.ridge_crossfit(Xps, W, seq_len(nrow(d)), cfg$K, cfg$seed + 501), 0.01), 0.99)
  cf <- grf::causal_forest(hte_X14(dd), Y, W, W.hat = W_hat, num.trees = cfg$num_trees, tune.parameters = cfg$tune,
                           seed = cfg$seed + 502)
  tc <- grf::test_calibration(cf)
  bl <- grf::best_linear_projection(cf, hte_X14(dd)[, c("vanco_6h", "baseline_aki_risk")])
  rbind(data.table(metric = paste0("insample_oob_", rownames(tc)), estimate = tc[, 1], se = tc[, 2],
                   lo = tc[, 1] - 1.96 * tc[, 2], hi = tc[, 1] + 1.96 * tc[, 2], p = tc[, 4]),
        data.table(metric = paste0("insample_oob_BLP:", rownames(bl)), estimate = bl[, 1], se = bl[, 2],
                   lo = bl[, 1] - 1.96 * bl[, 2], hi = bl[, 1] + 1.96 * bl[, 2], p = bl[, 4]))
}

#' Année d'admission réelle approximative (anchor_year_group réaligné)
approx_admit_year_v31 <- function(d) {
  mid <- as.numeric(substr(d$anchor_year_group, 1, 4)) + 1
  mid + (lubridate::year(d$admittime) - d$anchor_year)
}

#' Schéma 2b : dérivation (≤ 2016) -> test (≥ 2017), modèle figé
run_hte_period_v31 <- function(d, ps_covs, cfg, train_max_year = 2016, test_min_year = 2017) {
  yr <- if ("approx_admit_year" %in% names(d)) d$approx_admit_year else approx_admit_year_v31(d)
  tr <- which(yr <= train_max_year); te <- which(yr >= test_min_year)
  Xps <- build_ps_design_v31(d, ps_covs); Xmu <- build_ps_design_v31(d, ps_covs, df = 1)
  ctp <- .train_predict_cate(d, Xps, Xmu, tr, te, cfg, cfg$seed + 7000)
  thr <- stats::quantile(ctp$tau_tr_oob, c(1/3, 2/3), names = FALSE)   # seuils figés sur la dérivation
  # Γ du test : nuisances par validation croisée à l'intérieur du test
  W <- as.integer(d$treatment == "PTZ"); Y <- d$outcome_primary_v31
  fte <- .make_folds(W[te], Y[te], cfg$K, cfg$seed + 7100)
  sc <- rbindlist(lapply(seq_len(cfg$K), function(k)
    .aipw_scores(d, Xps, Xmu, te[fte != k], te[fte == k], cfg$seed + 7200 + k)))
  pp <- data.table(row = te, tau = ctp$tau, baseline_aki_risk = ctp$risk_te)[sc, on = "row"]
  pp[, `:=`(hadm_id = d$hadm_id[row], treatment = d$treatment[row], y = Y[row], approx_year = yr[row])]
  pp[, tertile := assign_tertiles_v31(tau, thr)$tertile]
  list(per_patient = pp, thresholds = thr, n_derivation = length(tr), n_test = length(te),
       overlap_ids = length(intersect(d$hadm_id[tr], d$hadm_id[te])))
}

#' Sensibilité MI : tiers figés (imputation 1), Γ recalculés sur plusieurs
#' imputations (nuisances ré-ajustées, mêmes plis), effets par tiers poolés (Rubin)
mi_tertile_effects_v31 <- function(imp_list, tertile_by_hadm, folds_by_hadm, ps_covs, cfg) {
  per <- lapply(seq_along(imp_list), function(i) {
    d <- imp_list[[i]]; Xps <- build_ps_design_v31(d, ps_covs); Xmu <- build_ps_design_v31(d, ps_covs, df = 1)
    f <- folds_by_hadm[match(d$hadm_id, names(folds_by_hadm))]
    sc <- rbindlist(lapply(seq_len(cfg$K), function(k) .aipw_scores(d, Xps, Xmu, which(f != k), which(f == k), cfg$seed + 9000 + 10 * i + k)))
    sc[, tertile := tertile_by_hadm[as.character(d$hadm_id[row])]]
    sc[, .(rd = mean(gamma), v = stats::var(gamma) / .N), by = tertile][, imp := i]
  })
  x <- rbindlist(per)
  x[, { p <- rubin_pool_v31(rd, v, 1e6); .(rd = p$est, rd_lo = p$lo, rd_hi = p$hi, m = .N) }, by = tertile]
}
