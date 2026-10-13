# =============================================================================
# 62_hte_figures.R — Figure HTE à 2 panneaux (structure de Buell et al.,
#   JAMA 2024, Figure 2), légende générée depuis les résultats, Tables HTE-1/2
# =============================================================================
# Couleurs des tiers : imposées par le cahier des charges (bleu / gris / rouge).

HTE_COL <- c(lower = "#1f77b4", middle = "#8c8c8c", upper = "#d62728")
HTE_LAB <- c(upper = "Upper third", middle = "Middle third", lower = "Lower third")
HTE_XLAB <- c(lower = "Lower third\n(smallest predicted PTZ excess)", middle = "Middle third",
              upper = "Upper third\n(largest predicted PTZ excess)")

.panel_tag <- function(tag) ggplot2::annotate("label", x = -Inf, y = Inf, label = tag, hjust = -0.2, vjust = 1.2,
                                              fontface = "bold", size = 4.2, label.r = grid::unit(0, "pt"),
                                              linewidth = 0.5, fill = "white")

.favors_arrows <- function(x_arrow, ylim, dx) {
  up <- 0.8 * ylim[2]; dn <- 0.8 * ylim[1]
  list(ggplot2::annotate("segment", x = x_arrow, xend = x_arrow, y = 0, yend = up, linewidth = 0.4,
                         arrow = grid::arrow(length = grid::unit(0.12, "cm"), type = "closed")),
       ggplot2::annotate("segment", x = x_arrow, xend = x_arrow, y = 0, yend = dn, linewidth = 0.4,
                         arrow = grid::arrow(length = grid::unit(0.12, "cm"), type = "closed")),
       # libellés à droite de la flèche, centrés sur chaque demi-flèche
       ggplot2::annotate("text", x = x_arrow + dx, y = up / 2, label = "Favors carbapenem", angle = 90, size = 2.7),
       ggplot2::annotate("text", x = x_arrow + dx, y = dn / 2, label = "Favors PTZ", angle = 90, size = 2.7))
}

.theme_hte <- function() ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = 8.5), legend.title = ggplot2::element_blank(),
                 legend.background = ggplot2::element_blank(), legend.key.height = grid::unit(0.35, "cm"))

#' @param pp     data.table une ligne par patient (tau, tertile)
#' @param tt     effets DR par tiers (dr_by_group_v31) + p_interaction
#' @param ato    ATO global (fraction) pour la ligne pointillée
fig_hte_tertiles_v31 <- function(pp, tt, ato, p_interaction, show_ci = FALSE) {
  d <- copy(pp)[order(tau)][, rank := .I]
  ylim <- range(c(100 * d$tau, 100 * tt$rd_lo, 100 * tt$rd_hi, 0)); ylim <- ylim + c(-0.1, 0.1) * diff(ylim)
  N <- nrow(d)
  pA <- ggplot2::ggplot(d, ggplot2::aes(rank, 100 * tau, colour = factor(tertile, levels = c("upper", "middle", "lower")))) +
    { if (show_ci) ggplot2::geom_ribbon(ggplot2::aes(ymin = 100 * (tau - 1.96 * sqrt(tau_var)), ymax = 100 * (tau + 1.96 * sqrt(tau_var))),
                                        fill = "grey85", colour = NA) } +
    ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = 100 * ato, linetype = "dashed", colour = "grey30", linewidth = 0.3) +
    ggplot2::annotate("text", x = 0.02 * N, y = 100 * ato, label = "Average effect (ATO)", hjust = 0, vjust = -0.5, size = 2.6, colour = "grey30") +
    ggplot2::geom_point(size = 0.6) +
    ggplot2::scale_colour_manual(values = HTE_COL[c("upper", "middle", "lower")], labels = HTE_LAB, breaks = c("upper", "middle", "lower")) +
    .favors_arrows(1.05 * N, ylim, 0.035 * N) +
    ggplot2::scale_x_continuous(limits = c(0, 1.11 * N), breaks = scales::pretty_breaks(4), expand = c(0.01, 0)) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::labs(x = "Patients ranked by predicted individualized treatment effect", y = "Predicted individualized treatment effect, %") +
    .theme_hte() + ggplot2::theme(legend.position = c(0.2, 0.86)) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2))) + .panel_tag("A")
  b <- copy(tt)[, tertile := factor(tertile, levels = c("lower", "middle", "upper"))]
  b[, lab := sprintf("n = %s; events PTZ/carb = %d/%d", format(n, big.mark = ","), events_ptz, events_carb)]
  pB <- ggplot2::ggplot(b, ggplot2::aes(tertile, 100 * rd, colour = tertile)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = 100 * rd_lo, ymax = 100 * rd_hi), width = 0.15, linewidth = 0.6) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_text(ggplot2::aes(y = 100 * rd_lo - 0.04 * diff(ylim), label = lab), vjust = 1, size = 2.4, colour = "grey20") +
    ggplot2::scale_colour_manual(values = HTE_COL, guide = "none") +
    ggplot2::scale_x_discrete(labels = HTE_XLAB, expand = ggplot2::expansion(add = c(0.6, 0.9))) +
    .favors_arrows(3.65, ylim, 0.14) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::labs(x = NULL, y = "Observed risk difference in AKI by day 7, %\n(PTZ minus carbapenem, doubly robust)") +
    .theme_hte() + .panel_tag("B")
  patchwork::wrap_plots(pA, pB, nrow = 1)
}

save_hte_fig_v31 <- function(p, path_noext, w_cm = 36, h_cm = 9) {
  ggplot2::ggsave(paste0(path_noext, ".png"), p, width = w_cm, height = h_cm, units = "cm", dpi = 300, bg = "white")
  ggplot2::ggsave(paste0(path_noext, ".pdf"), p, width = w_cm, height = h_cm, units = "cm", device = grDevices::cairo_pdf, bg = "white")
}

#' Légende prête pour le manuscrit, entièrement générée depuis les résultats
hte_legend_text_v31 <- function(tt, n, p_int, p_svy, scheme = c("crossfit", "period"), n_test = NULL, thresholds) {
  scheme <- match.arg(scheme)
  f <- function(r) sprintf("%+.1f%% (95%% CI, %.1f%% to %.1f%%)", 100 * r$rd, 100 * r$rd_lo, 100 * r$rd_hi)
  meth <- if (scheme == "crossfit")
    sprintf("Predictions are out-of-sample: the cohort (N = %s) was split into 5 folds and each patient's effect was predicted by models trained on the other 4 folds (average of 5 causal forests; baseline AKI risk and nuisance models also cross-fitted). Tertile cut points (%.1f%% and %.1f%%) were set on the pooled out-of-fold predictions.",
            format(n, big.mark = ","), 100 * thresholds[1], 100 * thresholds[2])
  else sprintf("The model was derived in patients admitted in or before 2016 (approximate year realigned from anchor_year_group) and applied unchanged to %s patients admitted from 2017 onward; tertile cut points (%.1f%% and %.1f%%) were fixed on the derivation predictions.",
               format(n_test, big.mark = ","), 100 * thresholds[1], 100 * thresholds[2])
  paste0(
    "A, Predicted individualized treatment effect (difference in the probability of acute kidney injury [AKI] by day 7 with piperacillin-tazobactam [PTZ] vs a carbapenem) for each patient, ranked from smallest to largest and colored by tertile. Positive values indicate more AKI with PTZ. ", meth,
    " B, Observed difference in AKI by day 7 between PTZ and carbapenem initiators within each tertile, estimated by augmented inverse probability weighting with nuisance models fitted on other folds (doubly robust; observational data, not randomized): lower third ", f(tt[tertile == "lower"]),
    "; middle third ", f(tt[tertile == "middle"]), "; upper third ", f(tt[tertile == "upper"]),
    sprintf(". P for interaction (2-df Wald test of equal tertile effects) = %s; overlap-weighted logistic treatment-by-tertile test P = %s. Error bars indicate 95%% CIs. Both panels share the same y-axis scale.",
            format.pval(p_int, digits = 2, eps = 0.001), format.pval(p_svy, digits = 2, eps = 0.001)))
}

#' Table HTE-1 : caractéristiques par tiers
table_hte1_v31 <- function(d, pp) {
  x <- copy(d)[, `:=`(tertile = pp$tertile[match(hadm_id, pp$hadm_id)], baseline_aki_risk = pp$baseline_aki_risk[match(hadm_id, pp$hadm_id)],
                      ps_oof = pp$e[match(hadm_id, pp$hadm_id)], ptz = as.integer(treatment == "PTZ"))]
  vars <- c(HTE_PREDICTORS_V31, "ptz", "ps_oof")
  rbindlist(lapply(vars, function(v) {
    z <- if (v == "gender") as.integer(toupper(x[[v]]) == "F") else as.numeric(x[[v]])
    bin <- all(stats::na.omit(unique(z)) %in% c(0, 1))
    st <- vapply(levels(x$tertile), function(t) { zz <- z[x$tertile == t]
      if (bin) sprintf("%d (%.1f)", sum(zz), 100 * mean(zz)) else {
        q <- stats::quantile(zz, c(.5, .25, .75)); sprintf("%.2f [%.2f; %.2f]", q[1], q[2], q[3]) } }, character(1))
    p <- if (bin) suppressWarnings(stats::chisq.test(table(z, x$tertile))$p.value) else stats::kruskal.test(z ~ x$tertile)$p.value
    data.table(variable = v, label = if (v %in% names(VAR_LABELS_V31)) VAR_LABELS_V31[[v]] else
      c(baseline_aki_risk = "Predicted baseline AKI risk", other_nephrotoxin_any = "Other nephrotoxin (48 h)", ptz = "PTZ initiators",
        ps_oof = "Propensity score (out-of-fold)")[[v]],
      stat = if (bin) "n (%)" else "median [IQR]", lower = st[["lower"]], middle = st[["middle"]], upper = st[["upper"]], p = p)
  }))
}

fig_hte_calibration_v31 <- function(cal) {
  ggplot2::ggplot(cal, ggplot2::aes(100 * mean_pred_tau, 100 * rd)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = 100 * rd_lo, ymax = 100 * rd_hi), width = 0, colour = "grey30") +
    ggplot2::geom_point(ggplot2::aes(shape = grouping), size = 2.4) +
    ggplot2::labs(x = "Mean predicted effect, %", y = "Observed doubly robust effect, %", shape = NULL) + .theme_hte()
}

fig_hte_importance_v31 <- function(vi) {
  d <- data.table(var = names(vi), imp = vi)[order(imp)]
  d[, label := ifelse(var %in% names(VAR_LABELS_V31), VAR_LABELS_V31[var],
                      c(baseline_aki_risk = "Predicted baseline AKI risk", other_nephrotoxin_any = "Other nephrotoxin (48 h)")[var])]
  d[, label := factor(label, levels = label)]
  ggplot2::ggplot(d, ggplot2::aes(imp, label)) + ggplot2::geom_col(width = 0.6, fill = "#1f77b4") +
    ggplot2::labs(x = "Variable importance (causal forest, mean over folds and seeds)", y = NULL) + .theme_hte()
}

fig_hte_dependence_v31 <- function(pp_long) {
  ggplot2::ggplot(pp_long, ggplot2::aes(level, 100 * mean_tau)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = 100 * lo, ymax = 100 * hi), width = 0.15, colour = "#1f77b4") +
    ggplot2::geom_point(colour = "#1f77b4", size = 2) +
    ggplot2::facet_wrap(~variable, scales = "free_x", nrow = 1) +
    ggplot2::labs(x = NULL, y = "Mean out-of-fold predicted effect, %") + .theme_hte() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}
