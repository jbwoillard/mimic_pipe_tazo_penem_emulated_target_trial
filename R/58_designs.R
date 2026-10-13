# =============================================================================
# 58_designs.R — Estimateurs de la phase B (docs/v3.1_plan.md B2-B4)
# =============================================================================

# -----------------------------------------------------------------------
# B3 — Essai à 3 bras : PS multinomial + generalized overlap weights
# (Li & Li 2019, Ann Appl Stat) : w_i = h(x_i) / e_{a_i}(x_i), avec
# h(x) = 1 / sum_k 1/e_k(x). Avec un PS logistique multinomial non pénalisé,
# ces poids donnent un équilibre exact des moyennes entre les trois bras.
# -----------------------------------------------------------------------
ARMS_3_V31 <- c("PTZ", "cefepime", "meropenem")

fit_ps_multinom_v31 <- function(X, arm) {
  df <- data.frame(arm = factor(arm, levels = ARMS_3_V31), X, check.names = TRUE)
  fit <- nnet::multinom(arm ~ ., data = df, trace = FALSE, maxit = 1000, MaxNWts = 100000)
  e <- stats::predict(fit, newdata = df, type = "probs")
  e <- pmin(pmax(e[, ARMS_3_V31, drop = FALSE], 1e-6), 1)
  e / rowSums(e)
}

gow_weights_v31 <- function(e, arm) {
  h <- 1 / rowSums(1 / e)
  h / e[cbind(seq_along(arm), match(arm, ARMS_3_V31))]
}

.gow_contrasts <- function(dt, vars, outcome_col) {
  X <- build_ps_design_v31(dt, vars)
  e <- fit_ps_multinom_v31(X, dt$treatment)
  w <- gow_weights_v31(e, dt$treatment)
  y <- dt[[outcome_col]]
  mu <- vapply(ARMS_3_V31, function(k) { s <- dt$treatment == k; sum(w[s] * y[s]) / sum(w[s]) }, numeric(1))
  list(mu = mu, w = w, e = e,
       rd = c(ptz_vs_mero = mu[["PTZ"]] - mu[["meropenem"]], ptz_vs_cef = mu[["PTZ"]] - mu[["cefepime"]],
              cef_vs_mero = mu[["cefepime"]] - mu[["meropenem"]]),
       lrr = c(ptz_vs_mero = log(mu[["PTZ"]] / mu[["meropenem"]]), ptz_vs_cef = log(mu[["PTZ"]] / mu[["cefepime"]]),
               cef_vs_mero = log(mu[["cefepime"]] / mu[["meropenem"]])))
}

gow_boot_v31 <- function(dt, ps_covs, outcome_col, B = 50, seed = SEED_V3) {
  vars <- intersect(ps_covs, names(dt))
  pt <- .gow_contrasts(dt, vars, outcome_col)
  set.seed(seed)
  bt <- lapply(seq_len(B), function(b) {
    idx <- unlist(lapply(ARMS_3_V31, function(k) { i <- which(dt$treatment == k); i[sample.int(length(i), replace = TRUE)] }))
    r <- tryCatch(.gow_contrasts(dt[idx], vars, outcome_col), error = function(e) NULL)
    if (is.null(r)) return(NULL); c(r$rd, setNames(r$lrr, paste0("lrr_", names(r$lrr))))
  })
  bt <- do.call(rbind, Filter(Negate(is.null), bt))
  ess <- vapply(ARMS_3_V31, function(k) { s <- dt$treatment == k; sum(pt$w[s])^2 / sum(pt$w[s]^2) }, numeric(1))
  list(n = nrow(dt), n_arm = table(factor(dt$treatment, levels = ARMS_3_V31)), mu = pt$mu, rd = pt$rd, lrr = pt$lrr,
       var_rd = apply(bt[, names(pt$rd)], 2, stats::var), var_lrr = setNames(apply(bt[, paste0("lrr_", names(pt$lrr))], 2, stats::var), names(pt$lrr)),
       weights = pt$w, ess = ess)
}

#' Trois contrastes sous MI : renvoie une liste d'objets compatibles
#' add_result_row() (un par contraste) + les risques pondérés par bras
gow_mi_v31 <- function(imp_list, ps_covs, outcome_col, B = 50, cores = getOption("v31.cores", 1L)) {
  res <- mc_retry_v31(length(imp_list), function(i) gow_boot_v31(imp_list[[i]], ps_covs, outcome_col, B, SEED_V3 + i), cores)
  if (any(vapply(res, inherits, logical(1), "try-error"))) stop("[58_v31] échec GOW")
  n <- res[[1]]$n; na <- res[[1]]$n_arm
  mu <- Reduce(`+`, lapply(res, `[[`, "mu")) / length(res)
  out <- list()
  for (cn in c("ptz_vs_mero", "ptz_vs_cef", "cef_vs_mero")) {
    prd <- rubin_pool_v31(sapply(res, function(r) r$rd[[cn]]), sapply(res, function(r) r$var_rd[[cn]]), n - length(ps_covs))
    prr <- rubin_pool_v31(sapply(res, function(r) r$lrr[[cn]]), sapply(res, function(r) r$var_lrr[[cn]]), n - length(ps_covs))
    a1 <- switch(cn, ptz_vs_mero = "PTZ", ptz_vs_cef = "PTZ", cef_vs_mero = "cefepime")
    a0 <- switch(cn, ptz_vs_mero = "meropenem", ptz_vs_cef = "cefepime", cef_vs_mero = "meropenem")
    out[[cn]] <- list(weight_scheme = "generalized_overlap", n = n, n_ptz = as.integer(na[[a1]]), n_carb = as.integer(na[[a0]]),
                      n_events = NA_integer_, weighted_risk_ptz = mu[[a1]], weighted_risk_carb = mu[[a0]],
                      crude_risk_ptz = NA_real_, crude_risk_carb = NA_real_,
                      rd = prd$est, rd_ci_low = prd$lo, rd_ci_high = prd$hi,
                      rr = exp(prr$est), rr_ci_low = exp(prr$lo), rr_ci_high = exp(prr$hi), fmi_rd = prd$fmi,
                      ess_pct = NA_real_, max_smd = NA_real_, n_boot = B * length(res), seed = SEED_V3,
                      alpha = NA_real_, lambda = NA_real_, m = length(res), ps_method = "multinomial_logistic")
  }
  message(sprintf("[58_v31] GOW 3 bras n=%d : PTZ-méro %+.2f%% | PTZ-céf %+.2f%% | céf-méro %+.2f%%", n,
                  100 * out$ptz_vs_mero$rd, 100 * out$ptz_vs_cef$rd, 100 * out$cef_vs_mero$rd))
  list(contrasts = out, mu = mu, weights_by_imputation = lapply(res, `[[`, "weights"),
       ess_by_arm = Reduce(`+`, lapply(res, `[[`, "ess")) / length(res))
}

#' Balance 3 bras : pour chaque paire de bras, un appel cobalt (via
#' compute_balance_v31) sur les mêmes poids GOW, imputations empilées ;
#' renvoie les SMD par paire et le max |SMD| / max KS sur les paires
compute_balance_3arm_v31 <- function(imp_list, ps_covs, weights_list) {
  pairs <- list(c("PTZ", "cefepime"), c("PTZ", "meropenem"), c("cefepime", "meropenem"))
  per <- rbindlist(lapply(pairs, function(pr) {
    sel <- lapply(imp_list, function(d) d$treatment %in% pr)
    b <- compute_balance_v31(lapply(seq_along(imp_list), function(i) imp_list[[i]][sel[[i]]]), ps_covs,
                             lapply(seq_along(imp_list), function(i) weights_list[[i]][sel[[i]]]), treated_level = pr[1])
    b[, .(term, pair = paste(pr, collapse = "_vs_"), Diff.Un, Diff.Adj, KS.Adj)]
  }))
  agg <- per[, .(Max.Diff.Un = max(abs(Diff.Un), na.rm = TRUE), Max.Diff.Adj = max(abs(Diff.Adj), na.rm = TRUE),
                 Max.KS.Adj = max(KS.Adj, na.rm = TRUE)), by = term]
  list(by_pair = per, max = agg)
}

# -----------------------------------------------------------------------
# B2 — Outcome ordinal ACORN : OR à cotes proportionnelles pondéré OW
# OR = céfépime vs PTZ (sens de l'essai : OR < 1 = moins de « pire issue »
# sous céfépime). Poids OW du PS binaire céfépime vs PTZ.
# -----------------------------------------------------------------------
.polr_or <- function(y, cef, w) {
  fit <- suppressWarnings(MASS::polr(factor(y, levels = 0:4, ordered = TRUE) ~ cef, weights = w / mean(w),
                                     method = "logistic", Hess = FALSE))
  unname(stats::coef(fit)[["cef"]])
}

ordinal_ow_boot_v31 <- function(dt, ps_covs, outcome_col, ps_method, B = 100, seed = SEED_V3) {
  vars <- intersect(ps_covs, names(dt))
  one <- function(d, lam = NULL) {
    X <- build_ps_design_v31(d, vars); a <- as.integer(d$treatment == "cefepime")
    f <- fit_ps_v31(X, a, ps_method, lambda = lam)
    w <- ifelse(a == 1, 1 - f$ps, f$ps)
    list(lor = .polr_or(d[[outcome_col]], a, w), w = w, lambda = f$lambda,
         dist = rbind(PTZ = tapply(w[a == 0], factor(d[[outcome_col]][a == 0], levels = 0:4), sum) / sum(w[a == 0]),
                      cefepime = tapply(w[a == 1], factor(d[[outcome_col]][a == 1], levels = 0:4), sum) / sum(w[a == 1])))
  }
  pt <- one(dt)
  set.seed(seed)
  bl <- vapply(seq_len(B), function(b) tryCatch(one(dt[sample.int(nrow(dt), replace = TRUE)],
                                                    if (ps_method == "ridge") pt$lambda else NULL)$lor,
                                                error = function(e) NA_real_), numeric(1))
  list(n = nrow(dt), lor = pt$lor, var_lor = stats::var(bl, na.rm = TRUE), dist = pt$dist, weights = pt$w)
}

ordinal_ow_mi_v31 <- function(imp_list, ps_covs, outcome_col, ps_method, B = 100, cores = getOption("v31.cores", 1L)) {
  res <- mc_retry_v31(length(imp_list), function(i)
    ordinal_ow_boot_v31(imp_list[[i]], ps_covs, outcome_col, ps_method, B, SEED_V3 + i), cores)
  if (any(vapply(res, inherits, logical(1), "try-error"))) stop("[58_v31] échec polr")
  p <- rubin_pool_v31(sapply(res, `[[`, "lor"), sapply(res, `[[`, "var_lor"), res[[1]]$n - length(ps_covs))
  dist <- Reduce(`+`, lapply(res, `[[`, "dist")) / length(res)
  message(sprintf("[58_v31] Ordinal OW (céfépime vs PTZ) n=%d : OR=%.3f [%.3f ; %.3f]", res[[1]]$n, exp(p$est), exp(p$lo), exp(p$hi)))
  list(n = res[[1]]$n, or = exp(p$est), or_lo = exp(p$lo), or_hi = exp(p$hi), fmi = p$fmi, dist = dist, m = length(res), B = B)
}

# -----------------------------------------------------------------------
# B4 — Per-protocol : censure artificielle + IPCW
# Stratégie : poursuivre le médicament index au moins 48 h (ou jusqu'à la
# sortie / au décès si < 48 h). Écart au protocole à t < 48 h :
#   - arrêt : pas de nouvelle administration du médicament index pendant 24 h
#     après une dose (écart daté à dose + 24 h, si < 48 h et patient encore
#     hospitalisé) ;
#   - switch : administration d'une autre bêta-lactamine anti-pyocyanique.
# Intervalles de 12 h (k = 0..3). Modèle de censure : logistique poolée par
# bras, covariables baseline (prédicteur linéaire du PS) + covariables
# dépendant du temps au début de l'intervalle (dernière créatinine et son
# rapport à la baseline, vasopresseur et ventilation invasive en cours,
# résultat de culture positive disponible, vancomycine dans les 12 h).
# Poids stabilisés cumulés (numérateur : intervalle seul), tronqués au P99.
# -----------------------------------------------------------------------
build_pp_panel_v31 <- function(cohort, ab, creat_dt, tbl_paths) {
  key <- cohort[, .(subject_id, hadm_id, T0, treatment, dischtime, deathtime, baseline_creat_min7d)]
  idx_cls <- function(tr) if (tr == "PTZ") "ptz" else c("meropenem", "carbapenem_other")
  abk <- merge(ab[iv == 1L], key[, .(subject_id, hadm_id, T0, treatment)], by = c("subject_id", "hadm_id"))
  abk[, h := as.numeric(difftime(t, T0, units = "hours"))]
  abk <- abk[h >= 0 & h <= 72]
  own <- abk[(treatment == "PTZ" & abx_class == "ptz") | (treatment != "PTZ" & abx_class %in% c("meropenem", "carbapenem_other"))]
  other <- abk[abx_class %in% c("ptz", "meropenem", "carbapenem_other", "cefepime", "ceftazidime", "aztreonam") &
                 !((treatment == "PTZ" & abx_class == "ptz") | (treatment != "PTZ" & abx_class %in% c("meropenem", "carbapenem_other")))]
  # arrêt : première dose suivie d'un vide > 24 h (avant 48 h) ; sans dose du tout -> écart à 0 h
  stop_t <- own[order(hadm_id, h)][, .(stop_h = {
    hh <- c(h, Inf); gap <- diff(hh); i <- which(gap > 24 & h < 48)[1]; if (is.na(i)) Inf else h[i] + 24 }),
    by = hadm_id]
  sw_t <- other[, .(switch_h = min(h)), by = hadm_id]
  p <- merge(key, stop_t, by = "hadm_id", all.x = TRUE)
  p <- merge(p, sw_t, by = "hadm_id", all.x = TRUE)
  p[is.na(stop_h), stop_h := fifelse(hadm_id %in% own$hadm_id, Inf, 0)]
  p[is.na(switch_h), switch_h := Inf]
  p[, end_h := pmax(0, pmin(48, as.numeric(difftime(fcoalesce(deathtime, dischtime, T0 + 48 * 3600), T0, units = "hours"))), na.rm = TRUE)]
  p[, dev_h := pmin(stop_h, switch_h)]
  p[, censored := as.integer(dev_h < end_h)]
  p[, cens_k := fifelse(censored == 1L, pmin(3L, as.integer(floor(pmin(dev_h, 1e6) / 12))), NA_integer_)]
  p[, last_k := fifelse(censored == 1L, cens_k, pmin(3L, as.integer(ceiling(pmax(end_h, 0.001) / 12)) - 1L))]
  # panel long : une ligne par intervalle à risque de censure
  pl <- p[, .(k = 0:last_k), by = .(hadm_id, subject_id, treatment, T0, baseline_creat_min7d, cens_k)]
  pl[, cens := as.integer(!is.na(cens_k) & k == cens_k)]
  pl[, t_start := T0 + k * 12 * 3600]
  # covariables dépendant du temps au début de l'intervalle
  cr <- merge(creat_dt[, .(subject_id, hadm_id, charttime, creat_mgdl)], p[, .(hadm_id, T0)], by = "hadm_id")
  cr <- cr[charttime <= T0 + 48 * 3600 & charttime >= T0 - 48 * 3600]
  last_cr <- cr[pl[, .(hadm_id, k, t_start)], on = .(hadm_id, charttime <= t_start), mult = "last",
                .(hadm_id, k = i.k, cr_last = x.creat_mgdl)]
  pl <- merge(pl, last_cr, by = c("hadm_id", "k"), all.x = TRUE)
  # rapport à la baseline ; NA (pas de créatinine ou baseline manquante) -> 1
  pl[, cr_ratio := fifelse(is.na(cr_last) | is.na(baseline_creat_min7d), 1, cr_last / baseline_creat_min7d)]
  vaso <- .read_vaso_vent_intervals_v31(tbl_paths, unique(p$hadm_id))
  ov <- function(type) {
    x <- vaso[kind == type]
    m <- x[pl[, .(hadm_id, k, t_start)], on = .(hadm_id, st <= t_start, en >= t_start), nomatch = NULL,
           .(hadm_id, k = i.k)]
    unique(m)[, flag := 1L]
  }
  for (type in c("vaso", "vent")) {
    f <- ov(type); setnames(f, "flag", paste0(type, "_now"))
    pl <- merge(pl, f, by = c("hadm_id", "k"), all.x = TRUE)
    pl[is.na(get(paste0(type, "_now"))), (paste0(type, "_now")) := 0L]
  }
  mic <- read_mimic_table(tbl_paths$paths[["microbiologyevents"]], col_select = c("hadm_id", "storetime", "org_name"))
  mic <- mic[!is.na(org_name) & org_name != ""][, `:=`(hadm_id = as.integer(hadm_id), storetime = lubridate::ymd_hms(storetime, quiet = TRUE))]
  pos <- mic[hadm_id %in% p$hadm_id & !is.na(storetime), .(first_pos = min(storetime)), by = hadm_id]
  pl <- merge(pl, pos, by = "hadm_id", all.x = TRUE)
  pl[, culture_pos_known := as.integer(!is.na(first_pos) & first_pos <= t_start)]
  vc <- merge(abk[abx_class == "vancomycin", .(hadm_id, h)], pl[, .(hadm_id, k)], by = "hadm_id", allow.cartesian = TRUE)
  vc <- unique(vc[h >= k * 12 - 12 & h < k * 12, .(hadm_id, k)])[, vanco_prev12h := 1L]
  pl <- merge(pl, vc, by = c("hadm_id", "k"), all.x = TRUE)
  pl[is.na(vanco_prev12h), vanco_prev12h := 0L]
  pl[, c("first_pos", "cens_k", "t_start") := NULL]
  list(person = p[, .(hadm_id, censored, dev_h, end_h, stop_h, switch_h)], panel = pl)
}

.read_vaso_vent_intervals_v31 <- function(tbl_paths, hadm_ids) {
  con <- duckdb::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(duckdb::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  duckdb::duckdb_register(con, "h", data.frame(hadm_id = hadm_ids))
  vaso_ids <- c(221662L, 221653L, 221289L, 221906L, 222315L, 229617L)
  v <- as.data.table(DBI::dbGetQuery(con, sprintf(
    "SELECT i.hadm_id, i.starttime st, i.endtime en FROM read_csv_auto('%s') i JOIN h ON i.hadm_id = h.hadm_id WHERE i.itemid IN (%s)",
    tbl_paths$paths[["inputevents"]], paste(vaso_ids, collapse = ","))))[, kind := "vaso"]
  m <- as.data.table(DBI::dbGetQuery(con, sprintf(
    "SELECT p.hadm_id, p.starttime st, p.endtime en FROM read_csv_auto('%s') p JOIN h ON p.hadm_id = h.hadm_id WHERE p.itemid = 225792",
    tbl_paths$paths[["procedureevents"]])))[, kind := "vent"]
  x <- rbind(v, m)
  x[, `:=`(hadm_id = as.integer(hadm_id), st = as.POSIXct(st, tz = "UTC"), en = as.POSIXct(en, tz = "UTC"))]
  x[!is.na(st) & !is.na(en)]
}

#' IPCW + OW pour un jeu imputé ; renvoie risques, RD, RR
.pp_point <- function(d, panel, person, vars, ps_method, lambda = NULL, trunc = 0.99) {
  X <- build_ps_design_v31(d, vars); a <- as.integer(d$treatment == "PTZ")
  f <- fit_ps_v31(X, a, ps_method, lambda = lambda)
  ow <- ifelse(a == 1, 1 - f$ps, f$ps)
  lp <- data.table(hadm_id = d$hadm_id, .row = seq_len(nrow(d)), ps_lp = stats::qlogis(f$ps))
  pl <- merge(panel, lp, by = "hadm_id", allow.cartesian = TRUE)   # .row : réplique bootstrap (doublons)
  pl[, `:=`(log_cr_ratio = log(pmax(cr_ratio, 0.1)), kf = factor(k))]
  sw <- pl[, {
    den <- stats::glm(cens ~ kf + ps_lp + log_cr_ratio + vaso_now + vent_now + culture_pos_known + vanco_prev12h,
                      family = stats::binomial(), data = .SD)
    num <- stats::glm(cens ~ kf, family = stats::binomial(), data = .SD)
    .(.row = .row, k = k, cens = cens, pd = stats::predict(den, type = "response"), pn = stats::predict(num, type = "response"))
  }, by = treatment]
  sw <- sw[order(.row, k)][, .(ipcw = prod((1 - pn) / (1 - pd)), cens_any = max(cens)), by = .row]
  sw[, ipcw := pmin(ipcw, stats::quantile(ipcw, trunc))]
  w <- rep(0, nrow(d)); w[sw$.row] <- ifelse(sw$cens_any == 1L, 0, sw$ipcw)
  w <- w * ow
  y <- d$outcome_primary_v31
  mu1 <- sum(w[a == 1] * y[a == 1]) / sum(w[a == 1]); mu0 <- sum(w[a == 0] * y[a == 0]) / sum(w[a == 0])
  list(mu1 = mu1, mu0 = mu0, rd = mu1 - mu0, rr = mu1 / mu0, lambda = f$lambda,
       n_uncensored = sum(sw$cens_any == 0L), max_ipcw = max(sw$ipcw))
}

pp_ipcw_boot_v31 <- function(d, pp, ps_covs, ps_method, B = 50, seed = SEED_V3) {
  vars <- intersect(ps_covs, names(d))
  panel <- pp$panel[hadm_id %in% d$hadm_id]
  pt <- .pp_point(d, panel, pp$person, vars, ps_method)
  set.seed(seed)
  bt <- t(vapply(seq_len(B), function(b) {
    idx <- sample.int(nrow(d), replace = TRUE); db <- d[idx]
    # chaque copie d'un patient tiré plusieurs fois reçoit son propre panel via
    # la clé .row construite dans .pp_point (jointure panel x lignes de db)
    r <- tryCatch(.pp_point(db, panel, pp$person, vars, ps_method, if (ps_method == "ridge") pt$lambda else NULL),
                  error = function(e) NULL)
    if (is.null(r)) c(NA_real_, NA_real_) else c(r$rd, log(r$rr))
  }, numeric(2)))
  list(n = nrow(d), n_ptz = sum(d$treatment == "PTZ"), n_carb = sum(d$treatment != "PTZ"),
       mu1 = pt$mu1, mu0 = pt$mu0, rd = pt$rd, rr = pt$rr,
       se_rd = stats::sd(bt[, 1], na.rm = TRUE), se_logrr = stats::sd(bt[, 2], na.rm = TRUE),
       n_uncensored = pt$n_uncensored, max_ipcw = pt$max_ipcw)
}
