# =============================================================================
# 76_cif.R — Incidence cumulée d'AKI J0-J14, Aalen-Johansen pondéré
#   (overlap), décès = risque concurrent (A4)
# =============================================================================
# Événement AKI : premier critère KDIGO strict rempli ou RRT (aki_h). Décès :
# deathtime / dod. Le suivi de l'AKI (créatinine) s'arrête à la sortie : un
# patient sorti vivant sans AKI est censuré à la sortie ; un décès hospitalier
# avant toute AKI est l'événement concurrent. Horizon 14 jours (336 h).
# IC : bootstrap des patients (jeu imputé n° 1) avec ré-estimation du PS
# (même méthode) dans chaque réplique.

.cif_status_v32 <- function(d, horizon_h = 336) {
  end <- pmin(fifelse(is.na(d$disch_h), horizon_h, d$disch_h), horizon_h)
  aki <- !is.na(d$aki_h) & d$aki_h > 0 & d$aki_h <= end
  dth <- !aki & !is.na(d$death_h) & d$death_h <= end
  time <- fifelse(aki, d$aki_h, fifelse(dth, d$death_h, end))
  data.table(time = pmax(time, 0.01) / 24, status = factor(fifelse(aki, "aki", fifelse(dth, "death", "censor")), levels = c("censor", "aki", "death")))
}

.cif_fit_v32 <- function(d, w, days) {
  st <- .cif_status_v32(d)
  arm <- ifelse(d$treatment == "PTZ", "PTZ", "comparator")
  fit <- survival::survfit(survival::Surv(st$time, st$status) ~ arm, weights = w / mean(w), id = seq_len(nrow(d)))
  sm <- summary(fit, times = days, extend = TRUE)
  strata <- sub("arm=", "", as.character(sm$strata))
  j <- which(colnames(sm$pstate) == "aki")
  out <- data.table(arm = strata, day = sm$time, cif_aki = sm$pstate[, j], cif_death = sm$pstate[, which(colnames(sm$pstate) == "death")])
  list(points = out, fit = fit)
}

#' @return list(curves = courbes pondérées par bras, table = CIF et RD à J2/3/7/14 avec IC)
cif_weighted_v32 <- function(d, covs, method, B = 200, days = c(2, 3, 7, 14), seed = SEED_V3) {
  # lambda ridge choisi UNE fois sur l'échantillon complet puis fixé dans les
  # répliques (comme les autres bootstraps v3.x) : fit_ps_v31() réinitialise la
  # graine quand il fait la validation croisée, ce qui figerait le tirage.
  X0 <- build_ps_design_v31(d, covs); a0 <- as.integer(d$treatment == "PTZ")
  f0 <- fit_ps_v31(X0, a0, method)
  lam <- if (method == "ridge") f0$lambda else NULL
  wfun <- function(dd) { X <- build_ps_design_v31(dd, covs); a <- as.integer(dd$treatment == "PTZ")
    ps <- fit_ps_v31(X, a, method, lambda = lam)$ps; ifelse(a == 1, 1 - ps, ps) }
  w <- ifelse(a0 == 1, 1 - f0$ps, f0$ps)
  pt <- .cif_fit_v32(d, w, days)$points
  grid <- .cif_fit_v32(d, w, seq(0, 14, by = 0.25))$points
  set.seed(seed)
  bt <- rbindlist(lapply(seq_len(B), function(b) {
    i <- sample.int(nrow(d), replace = TRUE); db <- d[i]
    r <- tryCatch(.cif_fit_v32(db, wfun(db), days)$points, error = function(e) NULL)
    if (is.null(r)) return(NULL); r[, rep := b]
  }))
  rd_b <- dcast(bt, rep + day ~ arm, value.var = "cif_aki")[, rd := PTZ - comparator]
  ci_arm <- bt[, .(lo = stats::quantile(cif_aki, 0.025), hi = stats::quantile(cif_aki, 0.975)), by = .(arm, day)]
  ci_rd <- rd_b[, .(rd_lo = stats::quantile(rd, 0.025), rd_hi = stats::quantile(rd, 0.975)), by = day]
  tabw <- dcast(pt, day ~ arm, value.var = "cif_aki")[, rd := PTZ - comparator]
  tab <- merge(tabw, ci_rd, by = "day")
  dth <- dcast(pt, day ~ arm, value.var = "cif_death"); setnames(dth, c("comparator", "PTZ"), c("cif_death_comparator", "cif_death_PTZ"))
  tab <- merge(tab, dth, by = "day")
  setnames(tab, c("comparator", "PTZ"), c("cif_aki_comparator", "cif_aki_PTZ"))
  list(table = tab, curves = grid, ci_arm = ci_arm, B = B, n = nrow(d))
}

fig_cif_v32 <- function(curves_list) {
  x <- rbindlist(lapply(names(curves_list), function(p) curves_list[[p]]$curves[, population := p]))
  x[, arm := ifelse(arm == "PTZ", "PTZ", "Meropenem")]
  x[, population := factor(population, levels = intersect(c("Early (primary)", "Broad (secondary)"), unique(population)))]
  ggplot2::ggplot(x, ggplot2::aes(day, 100 * cif_aki, linetype = arm)) +
    ggplot2::geom_step(linewidth = 0.55) +
    ggplot2::facet_wrap(~population, nrow = 1) +
    ggplot2::scale_x_continuous(breaks = c(0, 2, 3, 7, 14)) +
    ggplot2::labs(x = "Days after T0", y = "Cumulative incidence of AKI, %", linetype = NULL,
                  caption = "Overlap-weighted Aalen-Johansen estimator, death as a competing risk; follow-up censored at hospital discharge.") +
    theme_fig_v32()
}
