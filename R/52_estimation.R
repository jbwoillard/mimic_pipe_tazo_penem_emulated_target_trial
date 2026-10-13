# =============================================================================
# 52_estimation.R — PS, ATO, imputation multiple (Rubin), balance unique
#   (docs/v3.1_plan.md A2, A5)
# =============================================================================
# - Design du PS : splines naturelles (3 ddl) sur les continues, indicatrices
#   pour les catégorielles (référence = premier niveau), binaires telles quelles.
# - Deux méthodes de PS sur la MÊME matrice de design :
#     "glm"   : logistique non pénalisée (glm.fit). Avec l'overlap weighting,
#               elle donne un équilibre exact des moyennes de chaque colonne du
#               design (Li, Morgan & Zaslavsky 2018), donc des moyennes des
#               covariables (la base ns() contient la fonction linéaire).
#     "ridge" : glmnet alpha = 0, lambda par CV 5 plis sur l'échantillon
#               complet puis figé dans le bootstrap (choix v3 inchangé).
# - Bootstrap : ré-estimation complète (design + PS + poids + risques) à
#   chaque réplique ; RD et RR sur les mêmes répliques.
# - Imputation multiple : estimation dans chaque jeu imputé, variance intra
#   = variance bootstrap, combinaison par les règles de Rubin (RD ; RR sur
#   l'échelle log), ddl de Barnard & Rubin (1999).

CONTINUOUS_PS_VARS_V31 <- c(
  "age", "baseline_creat_min7d", "bun_baseline", "map_mean", "hr_mean", "temp_max",
  "spo2_min", "gcs_min", "time_admission_to_T0_h", "time_admission_to_T0_h_v3", "lactate_max",
  "platelets_min", "bilirubin_max", "weight_kg_v31"
)
CATEGORICAL_PS_VARS_V31 <- c("gender", "anchor_year_group", "service_group")

build_ps_design_v31 <- function(dt, vars, df = 3) {
  cont_vars <- intersect(vars, CONTINUOUS_PS_VARS_V31)
  cat_vars  <- union(intersect(vars, CATEGORICAL_PS_VARS_V31),
                     vars[vapply(vars, function(v) is.character(dt[[v]]) || is.factor(dt[[v]]), logical(1))])
  bin_vars  <- setdiff(vars, c(cont_vars, cat_vars))
  mats <- list()
  if (length(bin_vars) > 0) {
    Xb <- vapply(bin_vars, function(v) {
      x <- suppressWarnings(as.numeric(dt[[v]])); x[is.na(x)] <- 0; x
    }, numeric(nrow(dt)))
    mats$bin <- matrix(Xb, nrow = nrow(dt), dimnames = list(NULL, bin_vars))
  }
  for (v in cont_vars) {
    x <- suppressWarnings(as.numeric(dt[[v]]))
    if (anyNA(x)) x[is.na(x)] <- stats::median(x, na.rm = TRUE)  # seulement pour la sensibilité médiane
    b <- tryCatch(suppressWarnings(splines::ns(x, df = df)), error = function(e) matrix(x, ncol = 1))
    colnames(b) <- paste0(v, "_ns", seq_len(ncol(b)))
    mats[[v]] <- as.matrix(b)
  }
  for (v in cat_vars) {
    if (v == "gender") {
      mats[[v]] <- matrix(as.numeric(toupper(dt[[v]]) == "F"), ncol = 1, dimnames = list(NULL, "gender_F"))
      next
    }
    x <- as.character(dt[[v]]); x[is.na(x)] <- "unknown"
    lv <- sort(unique(x))[-1]
    if (length(lv) == 0) next
    mats[[v]] <- matrix(vapply(lv, function(l) as.numeric(x == l), numeric(length(x))),
                        nrow = length(x), dimnames = list(NULL, paste0(v, "_", lv)))
  }
  X <- do.call(cbind, mats)
  storage.mode(X) <- "double"
  # retirer les colonnes constantes (p.ex. un niveau absent d'une réplique)
  keep <- apply(X, 2, function(z) stats::var(z) > 0)
  X[, keep, drop = FALSE]
}

#' Ajuste le PS ; renvoie ps + l'information pour refaire le même ajustement
fit_ps_v31 <- function(X, trt, method = c("glm", "ridge"), lambda = NULL, seed = SEED_V3) {
  method <- match.arg(method)
  if (method == "glm") {
    fit <- suppressWarnings(stats::glm.fit(cbind(1, X), trt, family = stats::binomial()))
    ps <- fit$fitted.values
    lam <- NA_real_
  } else {
    if (is.null(lambda)) {
      set.seed(seed)
      foldid <- sample(rep(1:5, length.out = nrow(X)))
      cv <- glmnet::cv.glmnet(X, trt, family = "binomial", alpha = 0, foldid = foldid)
      lambda <- cv$lambda.min
    }
    fit <- glmnet::glmnet(X, trt, family = "binomial", alpha = 0, lambda = lambda)
    ps <- as.vector(stats::predict(fit, X, type = "response", s = lambda))
    lam <- lambda
  }
  list(ps = pmin(pmax(ps, 1e-6), 1 - 1e-6), lambda = lam, method = method)
}

.ow_risks <- function(ps, y, trt) {
  w <- ifelse(trt == 1, 1 - ps, ps)
  mu1 <- sum(w[trt == 1] * y[trt == 1]) / sum(w[trt == 1])
  mu0 <- sum(w[trt == 0] * y[trt == 0]) / sum(w[trt == 0])
  list(mu1 = mu1, mu0 = mu0, rd = mu1 - mu0, rr = mu1 / mu0, w = w)
}

.smd_vector <- function(X, trt, w) {
  vapply(seq_len(ncol(X)), function(j) {
    x <- X[, j]
    sp <- sqrt((stats::var(x[trt == 1]) + stats::var(x[trt == 0])) / 2)
    if (sp == 0) return(NA_real_)
    abs(sum(w[trt == 1] * x[trt == 1]) / sum(w[trt == 1]) -
          sum(w[trt == 0] * x[trt == 0]) / sum(w[trt == 0])) / sp
  }, numeric(1))
}

#' ATO sur UN jeu de données complet. Interface et sortie compatibles avec
#' boot_estimate() v3 (alimente directement add_result_row()).
#' @param dt doit déjà être restreint à la cohorte d'analyse (outcome non NA)
ato_boot_v31 <- function(dt, ps_covs, outcome_col, ps_method = "glm", B = 300,
                          seed = SEED_V3, treated_level = "PTZ", verbose = FALSE, ps_df = 3) {   # ps_df : v3.3 (PS parcimonieux)
  stopifnot(!anyNA(dt[[outcome_col]]))
  y <- as.integer(dt[[outcome_col]]); trt <- as.integer(dt$treatment == treated_level)
  vars <- intersect(ps_covs, names(dt))
  X <- build_ps_design_v31(dt, vars, df = ps_df)
  fit <- fit_ps_v31(X, trt, ps_method, seed = seed)
  pt <- .ow_risks(fit$ps, y, trt)
  # max |SMD| sur les covariables d'origine (continues : échelle brute,
  # catégorielles : max sur les niveaux), cohérent avec cobalt
  Xorig <- build_ps_design_v31(dt, vars, df = 1)
  max_smd <- max(.smd_vector(Xorig, trt, pt$w), na.rm = TRUE)
  set.seed(seed)
  n <- length(y); brd <- brr <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, replace = TRUE)
    if (sum(trt[idx]) < 5 || sum(1 - trt[idx]) < 5) next
    r <- tryCatch({
      Xb <- build_ps_design_v31(dt[idx], vars, df = ps_df)
      fb <- fit_ps_v31(Xb, trt[idx], ps_method, lambda = if (ps_method == "ridge") fit$lambda else NULL)
      .ow_risks(fb$ps, y[idx], trt[idx])
    }, error = function(e) NULL)
    if (!is.null(r)) { brd[b] <- r$rd; brr[b] <- r$rr }
  }
  ci_rd <- stats::quantile(brd, c(.025, .975), na.rm = TRUE)
  ci_rr <- stats::quantile(brr, c(.025, .975), na.rm = TRUE)
  ess <- sum(pt$w)^2 / sum(pt$w^2)
  if (verbose) message(sprintf("[52_est_v31] %s n=%d RD=%+.2f%% [%.2f;%.2f]", outcome_col, n,
                               pt$rd * 100, ci_rd[1] * 100, ci_rd[2] * 100))
  list(weight_scheme = "overlap", n = n, n_ptz = sum(trt), n_carb = sum(1 - trt), n_events = sum(y),
       weighted_risk_ptz = pt$mu1, weighted_risk_carb = pt$mu0,
       crude_risk_ptz = mean(y[trt == 1]), crude_risk_carb = mean(y[trt == 0]),
       rd = pt$rd, rd_ci_low = unname(ci_rd[1]), rd_ci_high = unname(ci_rd[2]),
       rr = pt$rr, rr_ci_low = unname(ci_rr[1]), rr_ci_high = unname(ci_rr[2]),
       se_rd = stats::sd(brd, na.rm = TRUE), se_logrr = stats::sd(log(brr), na.rm = TRUE),
       ess = ess, ess_pct = ess / n * 100, max_smd = max_smd, ps = fit$ps, weights = pt$w,
       lambda = fit$lambda, alpha = if (ps_method == "ridge") 0 else NA_real_,
       n_boot = B, n_boot_failed = sum(is.na(brd)), seed = seed, boot_rd = brd, boot_rr = brr,
       ps_covs = vars, ps_method = ps_method)
}

#' mclapply robuste : un worker tué (mémoire) renvoie NULL sans erreur ; les
#' imputations concernées sont relancées en série (résultat identique : la
#' graine est fixée DANS chaque tâche)
mc_retry_v31 <- function(n, fun, cores = getOption("v31.cores", 1L)) {
  res <- if (cores > 1) parallel::mclapply(seq_len(n), fun, mc.cores = cores) else lapply(seq_len(n), fun)
  bad <- which(vapply(res, function(r) is.null(r) || inherits(r, "try-error"), logical(1)))
  if (length(bad) > 0) {
    message("[52_est_v31] ", length(bad), " tâche(s) parallèle(s) sans résultat -> relance en série : ", paste(bad, collapse = ","))
    for (i in bad) res[[i]] <- fun(i)
  }
  res
}

#' Combinaison de Rubin pour un scalaire (estimations q, variances u)
rubin_pool_v31 <- function(q, u, n_com) {
  m <- length(q); qbar <- mean(q); ubar <- mean(u); b <- stats::var(q)
  tt <- ubar + (1 + 1 / m) * b
  lambda <- (1 + 1 / m) * b / tt
  df_old <- (m - 1) / max(lambda, 1e-12)^2
  df_obs <- (n_com + 1) / (n_com + 3) * n_com * (1 - lambda)
  df <- 1 / (1 / df_old + 1 / df_obs)          # Barnard & Rubin (1999)
  tq <- stats::qt(0.975, df)
  list(est = qbar, se = sqrt(tt), lo = qbar - tq * sqrt(tt), hi = qbar + tq * sqrt(tt),
       df = df, fmi = lambda, within = ubar, between = b)
}

#' ATO sous imputation multiple : ato_boot_v31 dans chacun des m jeux,
#' combinaison de Rubin. Parallélisé sur les imputations (fork).
#' @param imp_list liste de m data.tables complétées, déjà restreintes à la cohorte d'analyse
ato_mi_v31 <- function(imp_list, ps_covs, outcome_col, ps_method = "glm", B = 100,
                        seed = SEED_V3, cores = getOption("v31.cores", 1L), label = outcome_col, ps_df = 3) {
  m <- length(imp_list)
  run_one <- function(i) ato_boot_v31(imp_list[[i]], ps_covs, outcome_col, ps_method,
                                      B = B, seed = seed + i, ps_df = ps_df)
  res <- mc_retry_v31(m, run_one, cores)
  bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r$rd), logical(1))
  if (any(bad)) stop("[52_est_v31] ", sum(bad), " imputation(s) en échec pour ", label)
  n <- res[[1]]$n
  p_rd <- rubin_pool_v31(sapply(res, `[[`, "rd"), sapply(res, `[[`, "se_rd")^2, n_com = n - length(ps_covs))
  p_rr <- rubin_pool_v31(log(sapply(res, `[[`, "rr")), sapply(res, `[[`, "se_logrr")^2, n_com = n - length(ps_covs))
  r1 <- res[[1]]
  message(sprintf("[52_est_v31] MI %-34s n=%d RD=%+.2f%% [%.2f;%.2f] RR=%.3f (m=%d, B=%d, FMI=%.2f)",
                  label, n, p_rd$est * 100, p_rd$lo * 100, p_rd$hi * 100, exp(p_rr$est), m, B, p_rd$fmi))
  list(weight_scheme = "overlap", n = n, n_ptz = r1$n_ptz, n_carb = r1$n_carb, n_events = r1$n_events,
       weighted_risk_ptz = mean(sapply(res, `[[`, "weighted_risk_ptz")),
       weighted_risk_carb = mean(sapply(res, `[[`, "weighted_risk_carb")),
       crude_risk_ptz = r1$crude_risk_ptz, crude_risk_carb = r1$crude_risk_carb,
       rd = p_rd$est, rd_ci_low = p_rd$lo, rd_ci_high = p_rd$hi,
       rr = exp(p_rr$est), rr_ci_low = exp(p_rr$lo), rr_ci_high = exp(p_rr$hi),
       se_rd = p_rd$se, fmi_rd = p_rd$fmi, df_rd = p_rd$df,
       ess = mean(sapply(res, `[[`, "ess")), ess_pct = mean(sapply(res, `[[`, "ess_pct")),
       max_smd = mean(sapply(res, `[[`, "max_smd")),
       lambda = NA_real_, alpha = r1$alpha, n_boot = B * m, seed = seed,
       m = m, B_per_imputation = B, ps_method = ps_method,
       rd_by_imputation = sapply(res, `[[`, "rd"),
       weights_by_imputation = lapply(res, `[[`, "weights"),
       ps_by_imputation = lapply(res, `[[`, "ps"),
       # v3.3 : réplications bootstrap conservées (distribution, eFigure A3)
       boot_by_imputation = lapply(res, function(r) list(rd = r$rd, rr = r$rr, boot_rd = r$boot_rd, boot_rr = r$boot_rr)),
       lambda_by_imputation = sapply(res, `[[`, "lambda"))
}

# -----------------------------------------------------------------------
# Balance (A2) : UN SEUL appel cobalt::bal.tab sur les poids finaux, toutes
# imputations empilées (statistiques moyennées sur les imputations).
# Table 1 (SMD) et §12.2 lisent cet objet.
# -----------------------------------------------------------------------
compute_balance_v31 <- function(imp_list, ps_covs, weights_list, treated_level = "PTZ") {
  stk <- rbindlist(lapply(seq_along(imp_list), function(i)
    imp_list[[i]][, c("treatment", ps_covs), with = FALSE][, `:=`(.imp = i, .w = weights_list[[i]])]))
  df <- as.data.frame(stk[, ps_covs, with = FALSE])
  for (v in ps_covs) {
    if (v == "gender") df[[v]] <- as.numeric(toupper(df[[v]]) == "F")
    else if (is.character(df[[v]]) || is.factor(df[[v]])) df[[v]] <- factor(df[[v]])
    else df[[v]] <- as.numeric(df[[v]])
  }
  bal <- cobalt::bal.tab(x = df, treat = as.integer(stk$treatment == treated_level),
                         weights = stk$.w, imp = stk$.imp, method = "weighting",
                         s.d.denom = "pooled", binary = "std", continuous = "std",
                         stats = c("mean.diffs", "variance.ratios", "ks.statistics"),
                         un = TRUE, disp = c("means"), abs = FALSE)
  b <- as.data.table(bal$Balance.Across.Imputations, keep.rownames = "term")
  # statistiques moyennées sur les imputations sous les noms usuels ; la pire
  # imputation est conservée (Max.Diff.Adj / Min.Diff.Adj, Max.KS.Adj)
  for (s in c("Diff.Un", "V.Ratio.Un", "KS.Un", "Diff.Adj", "V.Ratio.Adj", "KS.Adj"))
    setnames(b, paste0("Mean.", s), s)
  b[, worst_abs_diff_adj := pmax(abs(Min.Diff.Adj), abs(Max.Diff.Adj))]
  b <- b[, .(term, Type, Diff.Un, Diff.Adj, worst_abs_diff_adj, V.Ratio.Un, V.Ratio.Adj, KS.Un, KS.Adj, Max.KS.Adj)]
  # nom de la covariable d'origine pour chaque terme (niveaux des facteurs)
  b[, variable := vapply(term, function(t) {
    hit <- ps_covs[t == ps_covs | startsWith(t, paste0(ps_covs, "_"))]
    if (length(hit) == 0) NA_character_ else hit[which.max(nchar(hit))]
  }, character(1))]
  b
}

# -----------------------------------------------------------------------
# Estimateurs de contrôle (A5 : sous imputation multiple, Rubin)
# -----------------------------------------------------------------------
#' G-computation (standardisation marginale) : GLM logistique de l'outcome
#' sur traitement + même design que le PS ; estimand = ATE de la cohorte.
gcomp_boot_v31 <- function(dt, ps_covs, outcome_col, B = 50, seed = SEED_V3) {
  vars <- intersect(ps_covs, names(dt))
  one <- function(d) {
    X <- build_ps_design_v31(d, vars); a <- as.integer(d$treatment == "PTZ")
    y <- as.integer(d[[outcome_col]])
    fit <- suppressWarnings(stats::glm.fit(cbind(1, a, X), y, family = stats::binomial()))
    cf <- fit$coefficients; cf[is.na(cf)] <- 0
    p1 <- stats::plogis(cbind(1, 1, X) %*% cf); p0 <- stats::plogis(cbind(1, 0, X) %*% cf)
    c(mu1 = mean(p1), mu0 = mean(p0))
  }
  pt <- one(dt); set.seed(seed)
  bt <- t(vapply(seq_len(B), function(b) tryCatch(one(dt[sample.int(nrow(dt), replace = TRUE)]),
                                                   error = function(e) c(mu1 = NA_real_, mu0 = NA_real_)), numeric(2)))
  list(n = nrow(dt), mu1 = pt[["mu1"]], mu0 = pt[["mu0"]], rd = pt[["mu1"]] - pt[["mu0"]],
       rr = pt[["mu1"]] / pt[["mu0"]], se_rd = stats::sd(bt[, 1] - bt[, 2], na.rm = TRUE),
       se_logrr = stats::sd(log(bt[, 1] / bt[, 2]), na.rm = TRUE))
}

#' AIPW (ATE) avec nuisances par forêts aléatoires (ranger) en cross-fitting
#' 5 plis ; variance par la fonction d'influence (pas de bootstrap).
#' PS borné à [0,01 ; 0,99]. Classe de modèles différente du PS primaire :
#' vrai contrôle de robustesse à la spécification.
aipw_if_v31 <- function(dt, ps_covs, outcome_col, seed = SEED_V3, num_trees = 300) {
  vars <- intersect(ps_covs, names(dt))
  X <- as.data.frame(build_ps_design_v31(dt, vars, df = 1))
  a <- as.integer(dt$treatment == "PTZ"); y <- as.integer(dt[[outcome_col]]); n <- nrow(dt)
  set.seed(seed); folds <- sample(rep(1:5, length.out = n))
  e <- m1 <- m0 <- numeric(n)
  for (k in 1:5) {
    tr <- folds != k; te <- folds == k
    fe <- ranger::ranger(x = X[tr, ], y = factor(a[tr]), probability = TRUE, num.trees = num_trees,
                         min.node.size = 10, seed = seed, num.threads = 1)
    e[te] <- predict(fe, X[te, ])$predictions[, "1"]
    for (arm in 0:1) {
      idx <- tr & a == arm
      fy <- ranger::ranger(x = X[idx, ], y = factor(y[idx], levels = 0:1), probability = TRUE,
                           num.trees = num_trees, min.node.size = 10, seed = seed, num.threads = 1)
      p <- predict(fy, X[te, ])$predictions[, "1"]
      if (arm == 1) m1[te] <- p else m0[te] <- p
    }
  }
  e <- pmin(pmax(e, 0.01), 0.99)
  phi1 <- m1 + a * (y - m1) / e; phi0 <- m0 + (1 - a) * (y - m0) / (1 - e)
  mu1 <- mean(phi1); mu0 <- mean(phi0)
  if_rd <- (phi1 - mu1) - (phi0 - mu0)
  if_lrr <- (phi1 - mu1) / mu1 - (phi0 - mu0) / mu0
  list(n = n, mu1 = mu1, mu0 = mu0, rd = mu1 - mu0, rr = mu1 / mu0,
       se_rd = stats::sd(if_rd) / sqrt(n), se_logrr = stats::sd(if_lrr) / sqrt(n))
}

#' Applique un estimateur (gcomp_boot_v31 / aipw_if_v31) à chaque imputation
#' et combine par Rubin ; sortie compatible add_result_row()
pool_estimator_mi_v31 <- function(imp_list, fun, ps_covs, outcome_col, label,
                                   cores = getOption("v31.cores", 1L), ...) {
  dots <- list(...)
  res <- mc_retry_v31(length(imp_list), function(i)
    do.call(fun, c(list(imp_list[[i]], ps_covs, outcome_col, seed = SEED_V3 + i), dots)), cores)
  if (any(vapply(res, inherits, logical(1), "try-error"))) stop("[52_est_v31] échec ", label)
  n <- res[[1]]$n
  prd <- rubin_pool_v31(sapply(res, `[[`, "rd"), sapply(res, `[[`, "se_rd")^2, n - length(ps_covs))
  prr <- rubin_pool_v31(log(sapply(res, `[[`, "rr")), sapply(res, `[[`, "se_logrr")^2, n - length(ps_covs))
  message(sprintf("[52_est_v31] MI %-34s n=%d RD=%+.2f%% [%.2f;%.2f]", label, n,
                  prd$est * 100, prd$lo * 100, prd$hi * 100))
  d0 <- imp_list[[1]]
  list(n = n, n_ptz = sum(d0$treatment == "PTZ"), n_carb = sum(d0$treatment != "PTZ"),
       n_events = sum(d0[[outcome_col]]),
       weighted_risk_ptz = mean(sapply(res, `[[`, "mu1")), weighted_risk_carb = mean(sapply(res, `[[`, "mu0")),
       crude_risk_ptz = NA_real_, crude_risk_carb = NA_real_,
       rd = prd$est, rd_ci_low = prd$lo, rd_ci_high = prd$hi,
       rr = exp(prr$est), rr_ci_low = exp(prr$lo), rr_ci_high = exp(prr$hi),
       ess_pct = NA_real_, max_smd = NA_real_, n_boot = NA_integer_, seed = SEED_V3,
       alpha = NA_real_, lambda = NA_real_, weight_scheme = NA_character_)
}
