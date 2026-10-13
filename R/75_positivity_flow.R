# =============================================================================
# 75_positivity_flow.R — Positivité (A2) et diagramme de flux (A3)
# =============================================================================
# Le PS montré est celui EFFECTIVEMENT utilisé par l'estimation (même design
# à splines, même méthode ridge/logistique selon la règle A2), ajusté sur le
# jeu imputé n° 1 ; l'AUC est calculée hors pli (5 plis, même spécification).

theme_fig_v32 <- function(base = 9) ggplot2::theme_classic(base_size = base) +
  ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold", size = base),
                 legend.position = "bottom", axis.text = ggplot2::element_text(size = 8), legend.text = ggplot2::element_text(size = 8))
save_fig_v32 <- function(p, path_noext, w_cm, h_cm) {
  ggplot2::ggsave(paste0(path_noext, ".png"), p, width = w_cm, height = h_cm, units = "cm", dpi = 300, bg = "white")
  ggplot2::ggsave(paste0(path_noext, ".pdf"), p, width = w_cm, height = h_cm, units = "cm", device = grDevices::cairo_pdf, bg = "white")
  invisible(paste0(path_noext, ".png"))
}

.oof_auc <- function(X, a, method, K = 5, seed = SEED_V3) {
  set.seed(seed); f <- sample(rep(seq_len(K), length.out = length(a))); p <- numeric(length(a))
  for (k in seq_len(K)) {
    tr <- f != k
    fit <- fit_ps_v31(X[tr, , drop = FALSE], a[tr], method, seed = seed + k)
    p[!tr] <- if (method == "glm") {
      g <- suppressWarnings(stats::glm.fit(cbind(1, X[tr, , drop = FALSE]), a[tr], family = stats::binomial()))
      cf <- g$coefficients; cf[is.na(cf)] <- 0; as.vector(stats::plogis(cbind(1, X[!tr, , drop = FALSE]) %*% cf))
    } else {
      cv <- glmnet::glmnet(X[tr, , drop = FALSE], a[tr], family = "binomial", alpha = 0, lambda = fit$lambda)
      as.vector(stats::predict(cv, X[!tr, , drop = FALSE], type = "response", s = fit$lambda))
    }
  }
  as.numeric(pROC::auc(pROC::roc(a, p, quiet = TRUE)))
}

.ess <- function(w) sum(w)^2 / sum(w^2)

#' Positivité d'une comparaison binaire (traité = `treated`)
positivity_binary_v32 <- function(d, covs, method, population, treated = "PTZ", comparator_label) {
  X <- build_ps_design_v31(d, covs); a <- as.integer(d$treatment == treated)
  ps <- fit_ps_v31(X, a, method)$ps
  ow <- ifelse(a == 1, 1 - ps, ps); iptw <- ifelse(a == 1, 1 / ps, 1 / (1 - ps))
  arm <- ifelse(a == 1, treated, comparator_label)
  tab <- data.table(population = population, ps_method = method, n = length(a),
                    auc_out_of_fold = .oof_auc(X, a, method),
                    pct_ps_lt_0.05_or_gt_0.95 = 100 * mean(ps < 0.05 | ps > 0.95), pct_ps_lt_0.1_or_gt_0.9 = 100 * mean(ps < 0.1 | ps > 0.9),
                    ps_min_treated = min(ps[a == 1]), ps_max_treated = max(ps[a == 1]), ps_min_comparator = min(ps[a == 0]), ps_max_comparator = max(ps[a == 0]),
                    ess_ow_treated = .ess(ow[a == 1]), ess_ow_comparator = .ess(ow[a == 0]),
                    ess_iptw_treated = .ess(iptw[a == 1]), ess_iptw_comparator = .ess(iptw[a == 0]),
                    max_w_ow = max(ow), p99_w_ow = stats::quantile(ow, 0.99), max_w_iptw = max(iptw), p99_w_iptw = stats::quantile(iptw, 0.99))
  list(table = tab, data = data.table(population = population, arm = arm, treated = a, ps = ps, w_ow = ow, w_iptw = iptw))
}

fig_ps_mirror_v32 <- function(psd) {
  psd <- copy(psd); psd[, side := ifelse(treated == 1, "PTZ (above)", "Comparator (below)")]
  ggplot2::ggplot(psd, ggplot2::aes(ps)) +
    ggplot2::geom_histogram(data = psd[treated == 1], ggplot2::aes(y = ggplot2::after_stat(count)), binwidth = 0.025,
                            boundary = 0, fill = "grey35", colour = "white", linewidth = 0.15) +
    ggplot2::geom_histogram(data = psd[treated == 0], ggplot2::aes(y = -ggplot2::after_stat(count)), binwidth = 0.025,
                            boundary = 0, fill = "grey75", colour = "white", linewidth = 0.15) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = c(0.1, 0.9), linetype = "dashed", linewidth = 0.3) +
    ggplot2::facet_wrap(~population, scales = "free_y", ncol = 3) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_y_continuous(labels = function(x) abs(x)) +
    ggplot2::labs(x = "Propensity score for PTZ (model used for estimation)", y = "Patients (PTZ above, comparator below)",
                  caption = "Dashed lines: Crump bounds (0.1, 0.9). Dark bars: PTZ; light bars: comparator (meropenem or cefepime).") +
    theme_fig_v32()
}

fig_weights_v32 <- function(psd) {
  x <- rbind(psd[, .(population, arm, scheme = "Overlap weights", w = w_ow)], psd[, .(population, arm, scheme = "Inverse probability weights", w = w_iptw)])
  ggplot2::ggplot(x, ggplot2::aes(w, linetype = arm)) +
    ggplot2::geom_density(linewidth = 0.4) + ggplot2::scale_x_log10() +
    ggplot2::facet_grid(scheme ~ population, scales = "free") +
    ggplot2::labs(x = "Weight (log scale)", y = "Density", linetype = NULL) + theme_fig_v32()
}

#' 3 bras : PS généralisés (multinomial) sur le jeu imputé n° 1
positivity_3arm_v32 <- function(d, covs) {
  X <- build_ps_design_v31(d, covs); e <- fit_ps_multinom_v31(X, d$treatment); w <- gow_weights_v31(e, d$treatment)
  arms <- ARMS_3_V31
  long <- rbindlist(lapply(arms, function(k) data.table(observed_arm = d$treatment, ps_arm = paste0("P(", k, ")"), ps = e[, k])))
  tab <- rbindlist(lapply(arms, function(k) data.table(population = "Three-arm (early)", arm = k, n = sum(d$treatment == k),
                                                      min_own_ps = min(e[d$treatment == k, k]), pct_own_ps_lt_0.05 = 100 * mean(e[d$treatment == k, k] < 0.05),
                                                      ess_gow = .ess(w[d$treatment == k]), max_w_gow = max(w[d$treatment == k]),
                                                      ess_iptw = .ess(1 / e[cbind(which(d$treatment == k), match(k, arms))]))))
  list(table = tab, long = long, e = data.table(e, observed_arm = d$treatment))
}

fig_gps_v32 <- function(long) {
  ggplot2::ggplot(long, ggplot2::aes(ps, linetype = observed_arm)) + ggplot2::geom_density(linewidth = 0.45) +
    ggplot2::facet_wrap(~ps_arm, ncol = 3) + ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Generalized propensity score", y = "Density", linetype = "Observed arm") + theme_fig_v32()
}

fig_ternary_v32 <- function(e) {
  if (!requireNamespace("ggtern", quietly = TRUE)) return(NULL)
  tryCatch(ggtern::ggtern(as.data.frame(e), ggplot2::aes(x = PTZ, y = cefepime, z = meropenem, shape = observed_arm)) +
             ggplot2::geom_point(size = 0.5, alpha = 0.5) + ggplot2::labs(shape = NULL) + ggtern::theme_bw(base_size = 8),
           error = function(err) { message("[75_v32] ggtern indisponible : ", conditionMessage(err)); NULL })
}

# -----------------------------------------------------------------------
# A3 — Diagramme de flux (ggplot), N par bras à chaque étape
# -----------------------------------------------------------------------
#' @param steps data.table(id, label, n, n_ptz, n_mero, n_cef, x, y)
#' @param edges data.table(from, to, excl_label, ex_x, ex_y) : boîte d'exclusion
#'        latérale optionnelle en (ex_x, ex_y), reliée à l'arête par un trait horizontal
fig_flow_v32 <- function(steps, edges, xlim = c(-5.6, 6.4)) {
  b <- copy(steps)
  b[, txt := paste0(label, "\nN = ", format(n, big.mark = ",", trim = TRUE),
                    ifelse(is.na(n_ptz), "", paste0("\nPTZ ", format(n_ptz, big.mark = ",", trim = TRUE),
                                                    ifelse(is.na(n_mero) | n_mero == 0, "", paste0(" | MEM ", format(n_mero, big.mark = ",", trim = TRUE))),
                                                    ifelse(is.na(n_cef) | n_cef == 0, "", paste0(" | FEP ", format(n_cef, big.mark = ",", trim = TRUE))))))]
  e <- merge(edges, b[, .(from = id, x0 = x, y0 = y)], by = "from")
  e <- merge(e, b[, .(to = id, x1 = x, y1 = y)], by = "to")
  ex <- e[!is.na(excl_label) & nzchar(excl_label)]
  # point de départ du trait vers la boîte d'exclusion : sur l'arête, à la hauteur ex_y
  ex[, xm := x0 + (x1 - x0) * (y0 - ex_y) / (y0 - y1)]
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = e, ggplot2::aes(x = x0, y = y0 - 0.42, xend = x1, yend = y1 + 0.42), linewidth = 0.35,
                          arrow = grid::arrow(length = grid::unit(4, "pt"), type = "closed")) +
    ggplot2::geom_segment(data = ex, ggplot2::aes(x = xm, y = ex_y, xend = ex_x - 0.05, yend = ex_y), linewidth = 0.3,
                          arrow = grid::arrow(length = grid::unit(3, "pt"), type = "closed")) +
    ggplot2::geom_label(data = ex, ggplot2::aes(x = ex_x, y = ex_y, label = excl_label), size = 2.4, hjust = 0,
                        label.r = grid::unit(0, "pt"), linewidth = 0.2, fill = "grey95", lineheight = 0.9) +
    ggplot2::geom_label(data = b, ggplot2::aes(x = x, y = y, label = txt), size = 2.6, label.r = grid::unit(1, "pt"),
                        linewidth = 0.35, lineheight = 0.9, fill = "white") +
    ggplot2::coord_cartesian(xlim = xlim, clip = "off") + ggplot2::theme_void() +
    ggplot2::labs(caption = "MEM, meropenem; FEP, cefepime. Each box lists the remaining patients by arm.") +
    ggplot2::theme(plot.caption = ggplot2::element_text(size = 8), plot.margin = ggplot2::margin(8, 8, 8, 8))
}
