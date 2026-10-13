# =============================================================================
# 20_bootstrap.R — Estimateur PS + bootstrap uniques (v3, Phase A)
# =============================================================================
# Corrige trois bugs transversaux identifiés dans la v2 (docs/v3_plan.md §0.3-0.5,
# §1.1-1.3) :
#   1. Le PS primaire n'est ajusté QU'UNE SEULE FOIS (fit_ps_primary_v3), et le
#      même objet est réutilisé pour l'estimation, l'audit et les figures — plus
#      de refit implicite avec une graine de CV différente.
#   2. Une seule fonction de bootstrap (boot_estimate) est réutilisée pour toutes
#      les analyses ATO/ATE ; RD et RR sont calculés sur les MÊMES répliques
#      bootstrap (même tirage d'indices), avec ré-estimation complète PS + poids
#      + risques à chaque réplique.
#   3. Compromis de calcul validé par l'utilisateur (docs/v3_plan.md §7.2) : le
#      paramètre de pénalisation lambda est sélectionné UNE FOIS par validation
#      croisée sur l'échantillon complet (fit_ps_primary_v3), puis chaque
#      réplique bootstrap ré-ajuste le PS à lambda FIXE (glmnet::glmnet, pas de
#      cv.glmnet) — un refit complet par réplique serait hors de portée en temps
#      de calcul une fois combiné à l'imputation multiple (m=20, voir §1.14).
#
# Dépend de : R/10_estimators.R (prepare_X), R/15_positivity_diagnostics.R
#             (add_timing_vars, PS_COVARIATES_PRIMARY) — chargés par run_all_v3.R.

# -----------------------------------------------------------------------
# Graine globale documentée (v3)
# -----------------------------------------------------------------------
# Toutes les analyses v3 partagent cette graine ; elle est aussi enregistrée
# dans output/extraction/results/session_info.rds (voir R/30_results_master.R).
SEED_V3 <- 20260927L

# Répertoire des tables de sortie. Défaut = output/extraction/tables (pipeline v3,
# inchangé) ; la v3.1 le redirige via options(v3.tables_dir = "output/extraction/tables").
v3_tables_path <- function(f) file.path(getOption("v3.tables_dir", "output/extraction/tables"), f)

# -----------------------------------------------------------------------
# A. Ajustement du PS primaire — UNE SEULE FOIS par analyse
# -----------------------------------------------------------------------
#' @param dt       data.table (déjà restreint aux lignes pertinentes, p.ex.
#'                 outcome non-NA, si l'appelant le souhaite)
#' @param ps_covs  covariables du PS
#' @param alpha    paramètre glmnet (0.1 en Phase A = valeur v2 inchangée pour
#'                 permettre un contrôle de non-régression ; passera à 0 en
#'                 Phase B, voir docs/v3_plan.md §1.7)
#' @param nfolds   nombre de folds pour la sélection de lambda par CV
#' @return liste : ps (vecteur), lambda (scalaire figé pour le bootstrap),
#'         alpha, foldid (pour audit/reproductibilité), vars, X (matrice de
#'         design), trt, cv_fit (objet cv.glmnet complet, pour diagnostic)
fit_ps_primary_v3 <- function(dt, ps_covs, alpha = 0.1, nfolds = 5,
                               seed = SEED_V3, treatment_col = "treatment",
                               treated_level = "PTZ", design_fn = NULL) {
  vars_avail <- ps_covs[ps_covs %in% names(dt)]
  trt <- as.integer(dt[[treatment_col]] == treated_level)
  build_X <- if (is.null(design_fn)) function(d, v) prepare_X(d, v) else design_fn
  X   <- build_X(dt, vars_avail)

  set.seed(seed)
  foldid <- sample(rep(seq_len(nfolds), length.out = nrow(X)))

  cv_fit <- glmnet::cv.glmnet(X, trt, family = "binomial", alpha = alpha,
                               foldid = foldid)
  lambda_hat <- cv_fit$lambda.min

  ps <- as.vector(predict(cv_fit, X, s = "lambda.min", type = "response"))
  ps <- pmax(pmin(ps, 0.99), 0.01)

  list(
    ps = ps, lambda = lambda_hat, alpha = alpha, foldid = foldid,
    vars = vars_avail, X = X, trt = trt, cv_fit = cv_fit,
    n = nrow(X), seed = seed, design_fn = build_X
  )
}

#' Ré-ajustement rapide du PS à lambda FIXE (utilisé dans chaque réplique
#' bootstrap — pas de validation croisée, donc ~100x plus rapide que
#' fit_ps_primary_v3 sur le même n).
refit_ps_fixed_lambda <- function(X_b, trt_b, alpha, lambda) {
  fit_b <- glmnet::glmnet(X_b, trt_b, family = "binomial", alpha = alpha,
                            lambda = lambda)
  ps_b <- as.vector(predict(fit_b, X_b, type = "response", s = lambda))
  pmax(pmin(ps_b, 0.99), 0.01)
}

# -----------------------------------------------------------------------
# B. Poids à partir d'un vecteur de PS, selon le schéma choisi
# -----------------------------------------------------------------------
.compute_weighted_risks <- function(ps_v, y_v, trt_v, weight_scheme, trunc_pct) {
  if (weight_scheme == "overlap") {
    w <- ifelse(trt_v == 1, 1 - ps_v, ps_v)
  } else if (weight_scheme == "iptw_trunc") {
    p_trt <- mean(trt_v)
    w_raw <- ifelse(trt_v == 1, p_trt / ps_v, (1 - p_trt) / (1 - ps_v))
    lo <- stats::quantile(w_raw, trunc_pct)
    hi <- stats::quantile(w_raw, 1 - trunc_pct)
    w  <- pmax(pmin(w_raw, hi), lo)
  } else {
    stop("[20_bootstrap_v3] weight_scheme inconnu : ", weight_scheme)
  }
  mu1 <- stats::weighted.mean(y_v[trt_v == 1], w[trt_v == 1])
  mu0 <- stats::weighted.mean(y_v[trt_v == 0], w[trt_v == 0])
  list(mu1 = mu1, mu0 = mu0, rd = mu1 - mu0, rr = mu1 / mu0, w = w)
}

# -----------------------------------------------------------------------
# C. LA fonction de bootstrap unique — réutilisée par toutes les analyses
#    ATO / ATE-tronqué de la v3. RD et RR proviennent des mêmes répliques.
# -----------------------------------------------------------------------
#' @param dt              data.table complet (sera restreint sur outcome non-NA)
#' @param ps_covs         covariables du PS
#' @param outcome_col     colonne outcome binaire
#' @param weight_scheme   "overlap" (ATO, primaire) ou "iptw_trunc" (ATE tronqué)
#' @param trunc_pct       percentile de troncature bilatérale pour iptw_trunc
#' @param alpha           paramètre glmnet (voir fit_ps_primary_v3)
#' @param B               nombre de répliques bootstrap (>= 1000, docs/v3_plan.md §1.3)
#' @param seed            graine (SEED_V3 par défaut, documentée)
#' @param ps_primary      objet retourné par fit_ps_primary_v3(), optionnel — si
#'                        fourni, évite un refit redondant du PS ponctuel (mais
#'                        chaque réplique bootstrap est TOUJOURS ré-estimée en
#'                        interne, à lambda fixe, conformément à la règle
#'                        "une estimation = un IC" du cahier des charges)
#' @return liste de résultats prête à alimenter results_master (voir
#'         R/30_results_master.R)
boot_estimate <- function(dt, ps_covs, outcome_col = "outcome_primary",
                           treatment_col = "treatment", treated_level = "PTZ",
                           weight_scheme = c("overlap", "iptw_trunc"),
                           trunc_pct = 0.01, alpha = 0.1, B = 1000,
                           seed = SEED_V3, ps_primary = NULL,
                           design_fn = NULL, verbose = TRUE) {
  weight_scheme <- match.arg(weight_scheme)

  ok  <- !is.na(dt[[outcome_col]])
  sub <- dt[ok]
  y   <- as.integer(sub[[outcome_col]])
  trt <- as.integer(sub[[treatment_col]] == treated_level)
  n   <- length(y)

  vars_avail <- ps_covs[ps_covs %in% names(sub)]
  # design_fn : fonction (data.table, vars) -> matrice numérique. Par défaut
  # prepare_X() (encodage simple). En Phase B (§1.7), la ridge + splines
  # naturelles passe une fonction basée sur splines::ns() — voir
  # R/33_ps_design_splines_v3.R:build_ps_design_splines_v3(). La matrice est
  # reconstruite depuis les données brutes à CHAQUE réplique bootstrap (pas
  # un simple X[idx,]) pour que la base de spline soit elle aussi
  # ré-estimée sur l'échantillon rééchantillonné, conformément à la règle
  # "ré-estimation complète ... à chaque réplique" (docs/v3_plan.md §0.4/§1.3).
  build_X <- if (is.null(design_fn)) function(d, v) prepare_X(d, v) else design_fn
  X <- build_X(sub, vars_avail)

  if (is.null(ps_primary)) {
    fit0 <- fit_ps_primary_v3(sub, vars_avail, alpha = alpha, seed = seed,
                               treatment_col = treatment_col,
                               treated_level = treated_level, design_fn = build_X)
    ps         <- fit0$ps
    lambda_hat <- fit0$lambda
  } else {
    ps <- if (length(ps_primary$ps) == nrow(dt)) ps_primary$ps[ok] else ps_primary$ps
    lambda_hat <- ps_primary$lambda
    if (!is.null(ps_primary$design_fn)) build_X <- ps_primary$design_fn
  }

  point <- .compute_weighted_risks(ps, y, trt, weight_scheme, trunc_pct)

  ess  <- (sum(point$w))^2 / sum(point$w^2)

  .smd_numeric <- function(x) {
    m1 <- stats::weighted.mean(x[trt == 1], point$w[trt == 1], na.rm = TRUE)
    m0 <- stats::weighted.mean(x[trt == 0], point$w[trt == 0], na.rm = TRUE)
    sp <- sqrt((stats::var(x[trt == 1], na.rm = TRUE) +
                  stats::var(x[trt == 0], na.rm = TRUE)) / 2)
    if (is.na(sp) || sp == 0) return(NA_real_)
    abs(m1 - m0) / sp
  }
  smd_fn <- function(v) {
    raw <- sub[[v]]
    if (v == "gender") {
      # caractère "F"/"M" : as.numeric() direct produirait NA silencieusement
      # (même défaut que 10_estimators_v2_positivity.R:smd_fn_ow en v2, qui
      # excluait gender de max_smd_ow sans le signaler).
      return(.smd_numeric(as.numeric(toupper(raw) == "F")))
    }
    if (is.character(raw) || is.factor(raw)) {
      # Variable catégorielle multi-niveaux (ex. anchor_year_group,
      # infection_site_proxy, service_at_t0, §2.5/§1.13) : as.numeric()
      # produirait aussi NA silencieusement et l'exclurait du diagnostic
      # (même défaut). Un SMD scalaire n'a pas de sens pour une variable
      # nominale à >2 niveaux : on encode chaque niveau en indicatrice et on
      # retient le SMD maximal parmi les niveaux, comme le fait cobalt.
      x <- as.character(raw); x[is.na(x)] <- "unknown_or_missing"
      lv <- unique(x)
      if (length(lv) <= 1) return(NA_real_)
      lvl_smds <- vapply(lv, function(l) .smd_numeric(as.numeric(x == l)), numeric(1))
      return(max(lvl_smds, na.rm = TRUE))
    }
    .smd_numeric(as.numeric(raw))
  }
  smd_vals <- vapply(vars_avail, smd_fn, numeric(1))
  max_smd  <- max(smd_vals, na.rm = TRUE)

  if (verbose) {
    message("[20_bootstrap_v3] ", weight_scheme, " point estimate: RD=",
            round(point$rd, 4), " RR=", round(point$rr, 3),
            " (n=", n, ", B=", B, ")")
  }

  set.seed(seed)
  boot_rd <- numeric(B)
  boot_rr <- numeric(B)
  is_default_design <- is.null(design_fn) && is.null(ps_primary$design_fn)
  for (b in seq_len(B)) {
    idx   <- sample.int(n, replace = TRUE)
    trt_b <- trt[idx]
    y_b   <- y[idx]
    if (sum(trt_b) < 5 || sum(1 - trt_b) < 5) {
      boot_rd[b] <- NA_real_; boot_rr[b] <- NA_real_; next
    }
    # Design par defaut (prepare_X, encodage lineaire) : indexer X suffit et
    # evite de reconstruire la matrice a chaque replique (plus rapide).
    # Design custom (ex. splines) : reconstruire depuis sub[idx] pour que la
    # base de spline soit elle aussi re-estimee sur l'echantillon retire.
    X_b <- if (is_default_design) X[idx, , drop = FALSE] else build_X(sub[idx], vars_avail)
    ps_b <- tryCatch(refit_ps_fixed_lambda(X_b, trt_b, alpha, lambda_hat),
                      error = function(e) ps[idx])
    res_b <- tryCatch(.compute_weighted_risks(ps_b, y_b, trt_b, weight_scheme, trunc_pct),
                       error = function(e) list(rd = NA_real_, rr = NA_real_))
    boot_rd[b] <- res_b$rd
    boot_rr[b] <- res_b$rr
  }

  ci_rd <- stats::quantile(boot_rd, c(0.025, 0.975), na.rm = TRUE)
  ci_rr <- stats::quantile(boot_rr, c(0.025, 0.975), na.rm = TRUE)
  n_boot_failed <- sum(is.na(boot_rd))
  if (verbose && n_boot_failed > 0) {
    message("[20_bootstrap_v3]   ", n_boot_failed, "/", B,
            " répliques bootstrap ignorées (bras trop petit ou refit PS en échec)")
  }

  list(
    weight_scheme   = weight_scheme,
    n = n, n_ptz = sum(trt), n_carb = sum(1 - trt), n_events = sum(y),
    weighted_risk_ptz  = unname(point$mu1),
    weighted_risk_carb = unname(point$mu0),
    crude_risk_ptz     = mean(y[trt == 1]),
    crude_risk_carb    = mean(y[trt == 0]),
    rd = unname(point$rd), rd_ci_low = unname(ci_rd[1]), rd_ci_high = unname(ci_rd[2]),
    rr = unname(point$rr), rr_ci_low = unname(ci_rr[1]), rr_ci_high = unname(ci_rr[2]),
    se_rd = stats::sd(boot_rd, na.rm = TRUE),
    ess = ess, ess_pct = ess / n * 100,
    max_smd = max_smd,
    ps = ps, weights = point$w, lambda = lambda_hat, alpha = alpha,
    n_boot = B, n_boot_failed = n_boot_failed, seed = seed,
    boot_rd = boot_rd, boot_rr = boot_rr,
    ps_covs = vars_avail
  )
}
