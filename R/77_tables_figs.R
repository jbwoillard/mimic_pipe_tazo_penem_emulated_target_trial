# =============================================================================
# 77_tables_figs.R — Risques et ESS par bras (A6), précision / MDE (A8),
#   trajectoires (B4), Table 1 3 bras (A5), forest plots (B4), ACORN essai vs
#   émulation (B6)
# =============================================================================

#' Événements bruts et ESS par bras pour une estimation MI (poids par imputation)
enrich_arm_stats_v32 <- function(est, d_list, outcome) {
  d1 <- d_list[[1]]; a <- d1$treatment == "PTZ"
  ess <- function(w, s) sum(w[s])^2 / sum(w[s]^2)
  wl <- est$weights_by_imputation
  list(events_ptz = sum(d1[[outcome]][a]), events_comp = sum(d1[[outcome]][!a]),
       ess_ptz = if (is.null(wl)) NA_real_ else mean(vapply(wl, function(w) ess(w, a), numeric(1))),
       ess_comp = if (is.null(wl)) NA_real_ else mean(vapply(wl, function(w) ess(w, !a), numeric(1))))
}

#' Demi-largeur d'IC, SE et effet minimal détectable à 80 % (alpha 5 %, bilatéral)
precision_cols_v32 <- function(rd_lo, rd_hi) {
  se <- (rd_hi - rd_lo) / (2 * 1.96)
  list(ci_half_width = (rd_hi - rd_lo) / 2, se_rd_from_ci = se, mde80 = (stats::qnorm(0.975) + stats::qnorm(0.80)) * se)
}

#' Trajectoires pondérées J0-J7 : moyennes par bras (valeurs absolues),
#' différence avec IC bootstrap (poids fixés), N contributif par jour et bras.
#' Fenêtre identique pour tous les marqueurs : jour = floor(t), t dans [0 ; 8) j.
trajectory_v32 <- function(cohort, w, values, value_col, marker, B = 200, seed = SEED_V3) {
  v <- merge(values, cohort[, .(subject_id, hadm_id, T0, treatment)], by = c("subject_id", "hadm_id"))
  v[, t := as.numeric(difftime(charttime, T0, units = "days"))]
  v <- v[t >= 0 & t < 8][, day := floor(t)]
  daily <- v[, .(val = mean(get(value_col))), by = .(subject_id, hadm_id, treatment, day)]
  daily <- merge(daily, data.table(hadm_id = cohort$hadm_id, w = w), by = "hadm_id")
  set.seed(seed)
  rbindlist(lapply(0:7, function(dd) {
    s <- daily[day == dd]; if (nrow(s) < 10) return(NULL)
    wm <- function(x) c(ptz = stats::weighted.mean(x$val[x$treatment == "PTZ"], x$w[x$treatment == "PTZ"]),
                        comp = stats::weighted.mean(x$val[x$treatment != "PTZ"], x$w[x$treatment != "PTZ"]))
    p <- wm(s); bt <- replicate(B, { r <- wm(s[sample.int(.N, replace = TRUE)]); r[["ptz"]] - r[["comp"]] })
    data.table(marker = marker, day = dd, n_ptz = sum(s$treatment == "PTZ"), n_comp = sum(s$treatment != "PTZ"),
               mean_ptz = p[["ptz"]], mean_comp = p[["comp"]], diff = p[["ptz"]] - p[["comp"]],
               diff_lo = stats::quantile(bt, 0.025, na.rm = TRUE), diff_hi = stats::quantile(bt, 0.975, na.rm = TRUE))
  }))
}

fig_trajectory_v32 <- function(tr) {
  ab <- rbind(tr[, .(marker, population, day, arm = "PTZ", value = mean_ptz, n = n_ptz)],
              tr[, .(marker, population, day, arm = "Meropenem", value = mean_comp, n = n_comp)])
  p1 <- ggplot2::ggplot(ab, ggplot2::aes(day, value, linetype = arm, shape = arm)) +
    ggplot2::geom_line(linewidth = 0.5) + ggplot2::geom_point(size = 1.6) +
    ggplot2::facet_grid(marker ~ population, scales = "free_y") + ggplot2::scale_x_continuous(breaks = 0:7) +
    ggplot2::labs(x = "Day after T0", y = "Overlap-weighted mean", linetype = NULL, shape = NULL) + theme_fig_v32()
  p2 <- ggplot2::ggplot(tr, ggplot2::aes(day, diff)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = diff_lo, ymax = diff_hi), width = 0.2, linewidth = 0.4) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::geom_text(ggplot2::aes(y = -Inf, label = paste0(n_ptz, "\n", n_comp)), vjust = -0.25, size = 1.9, lineheight = 0.85, colour = "grey30") +
    ggplot2::facet_grid(marker ~ population, scales = "free_y") + ggplot2::scale_x_continuous(breaks = 0:7) +
    ggplot2::labs(x = "Day after T0", y = "PTZ minus meropenem (95% CI)",
                  caption = "Numbers: patients contributing a value that day (PTZ above, meropenem below).") + theme_fig_v32()
  patchwork::wrap_plots(p1, p2, ncol = 1)
}

#' Table 1 du 3 bras : brut par bras, pondéré GOW par bras, SMD par paire
table1_3arm_v32 <- function(raw, imp_list, weights_list, covs, bal_by_pair) {
  arms <- ARMS_3_V31
  rows <- lapply(covs, function(v) {
    x <- raw[[v]]; is_cat <- v != "gender" && (is.character(x) || is.factor(x))
    if (is_cat) return(NULL)
    z <- if (v == "gender") as.numeric(toupper(x) == "F") else as.numeric(x)
    bin <- all(stats::na.omit(unique(z)) %in% c(0, 1))
    st <- vapply(arms, function(k) { zz <- z[raw$treatment == k]
      if (bin) sprintf("%d (%.1f)", sum(zz, na.rm = TRUE), 100 * mean(zz, na.rm = TRUE)) else {
        q <- stats::quantile(zz, c(.5, .25, .75), na.rm = TRUE); sprintf("%.1f [%.1f; %.1f]", q[1], q[2], q[3]) } }, character(1))
    wst <- vapply(arms, function(k) mean(vapply(seq_along(imp_list), function(i) { d <- imp_list[[i]]; s <- d$treatment == k
      zz <- if (v == "gender") as.numeric(toupper(d[[v]]) == "F") else as.numeric(d[[v]]); stats::weighted.mean(zz[s], weights_list[[i]][s]) }, numeric(1))), numeric(1))
    sm <- bal_by_pair[term == v]
    data.table(variable = v, characteristic = lab_v31(v), stat = if (bin) "n (%)" else "median [IQR]",
               PTZ = st[["PTZ"]], cefepime = st[["cefepime"]], meropenem = st[["meropenem"]],
               w_PTZ = wst[["PTZ"]], w_cefepime = wst[["cefepime"]], w_meropenem = wst[["meropenem"]],
               smd_ptz_cef = sm[pair == "PTZ_vs_cefepime", Diff.Adj], smd_ptz_mero = sm[pair == "PTZ_vs_meropenem", Diff.Adj],
               smd_cef_mero = sm[pair == "cefepime_vs_meropenem", Diff.Adj],
               pct_missing = sprintf("%.1f / %.1f / %.1f", 100 * mean(is.na(x[raw$treatment == "PTZ"])), 100 * mean(is.na(x[raw$treatment == "cefepime"])),
                                     100 * mean(is.na(x[raw$treatment == "meropenem"]))))
  })
  rbindlist(rows)
}

#' Forest plot : une ligne par population x estimand, même outcome (AKI J7)
fig_forest_v32 <- function(rows, xlab = "Risk difference in AKI by day 7, percentage points (95% CI)") {
  d <- copy(rows)[, yid := factor(sprintf("%02d", .I), levels = rev(sprintf("%02d", .I)))]
  ggplot2::ggplot(d, ggplot2::aes(100 * rd, yid)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = 100 * rd_lo, xmax = 100 * rd_hi), height = 0, linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(shape = is_primary), size = 2.4, fill = "white") +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 23, `FALSE` = 16), guide = "none") +
    ggplot2::scale_y_discrete(labels = setNames(d$label, d$yid)) +
    ggplot2::labs(x = xlab, y = NULL) + theme_fig_v32()
}

#' ACORN publié (Qian et al., JAMA 2023;330:1557-1567, Table 1 et 2, PMC10576861)
#' vs émulation. Valeurs publiées recopiées telles que publiées (stades d'AKI à
#' l'inclusion séparés, non recombinés à la main).
ACORN_PUBLISHED_V32 <- data.table(
  characteristic = c("Age, median (IQR), y", "Female sex, %", "Enrolled in the ED, %", "Enrolled in the ICU, %", "Sepsis, %",
                     "Invasive mechanical ventilation, %", "Vancomycin at enrollment, %", "Baseline creatinine, median (IQR), mg/dL",
                     "Chronic kidney disease, %", "AKI stage 1 at enrollment, %", "AKI stage 2 at enrollment, %", "AKI stage 3 at enrollment, %",
                     "In-hospital death by day 14, %"),
  trial_cefepime = c("57 (42-68)", "43.1", "93.5", "6.5", "54.2", "9.1", "77.6", "0.7 (0.6-0.8)", "20.4", "19.0", "11.0", "12.2", "7.6"),
  trial_ptz = c("59 (44-69)", "42.3", "95.8", "4.2", "54.3", "7.3", "76.9", "0.8 (0.6-0.9)", "20.3", "24.0", "9.5", "11.1", "6.0"))

acorn_trial_vs_emulation_v32 <- function(em) {
  f <- function(s, arm) { x <- em[treatment == arm]
    q <- function(v) { z <- stats::quantile(v, c(.5, .25, .75), na.rm = TRUE); sprintf("%.1f (%.1f-%.1f)", z[1], z[2], z[3]) }
    p <- function(v) sprintf("%.1f", 100 * mean(v, na.rm = TRUE))
    c(q(x$age), p(toupper(x$gender) == "F"), p(x$service_group == "ED_observation"), p(x$service_group == "ICU_medical"),
      p(x$sepsis_composite_v3 == 1), p(x$vent_pre_t0 == 1), p(x$vanco_6h == 1), q(x$baseline_creat_min7d), p(x$ckd_priorhosp == 1),
      NA, NA, NA, p(x$death_14d == 1)) }
  out <- copy(ACORN_PUBLISHED_V32)
  out[, `:=`(emulation_cefepime = f(em, "cefepime"), emulation_ptz = f(em, "PTZ"))]
  out[, note := ""]
  out[characteristic == "Sepsis, %", note := "Emulation: pragmatic organ-dysfunction composite (vasopressor, ventilation or lactate >= 2), not the trial definition"]
  out[grepl("^AKI stage", characteristic), note := "Emulation excludes AKI at T0 (pre-specified); see the sensitivity analysis including it"]
  out[characteristic == "Chronic kidney disease, %", note := "Emulation: ICD codes from prior MIMIC-IV hospitalizations only (lower bound)"]
  out
}

#' Libellé lisible d'un terme de balance (variables, niveaux de facteurs,
#' indicateur de manquant cobalt « var:<NA> »)
nice_term_v32 <- function(t) vapply(t, function(x) {
  extra <- c(creat_last_pre_t0 = "Last pre-T0 creatinine, mg/dL", creat_delta_last_min7d = "Last minus 7-day minimum creatinine",
             creat_ratio_last_min7d = "Last / 7-day minimum creatinine", creat_slope_48h = "Creatinine slope over 48 h, mg/dL/day",
             creat_slope_48h_available = "Creatinine slope computable (>=2 values)", creat_single_pre_value = "Single pre-T0 creatinine value")
  labs <- c(VAR_LABELS_V31, extra)
  if (grepl(":<NA>$", x)) { v <- sub(":<NA>$", "", x)
    return(paste0(if (v %in% names(labs)) labs[[v]] else v, ": not measured (structurally missing, ED stratum)")) }
  if (x %in% names(labs)) return(unname(labs[[x]]))
  v <- names(labs)[startsWith(x, paste0(names(labs), "_"))]
  if (length(v) == 0) return(x)
  v <- v[which.max(nchar(v))]; paste0(labs[[v]], ": ", gsub("_", " ", substring(x, nchar(v) + 2)))
}, character(1))

#' Love plot toutes populations (A2/C3), libellés lisibles, population primaire en premier
fig_love_all_v32 <- function(bal_list, bal3_max) {
  pops <- c(early = "Early (primary)", broad = "Broad (secondary)", acorn = "ACORN emulation")
  lp <- rbind(rbindlist(lapply(names(pops), function(p) bal_list[[p]][, .(term, smd_raw = abs(Diff.Un), smd_w = abs(Diff.Adj), population = pops[[p]])])),
              bal3_max[, .(term, smd_raw = Max.Diff.Un, smd_w = Max.Diff.Adj, population = "Three-arm (max pairwise)")])
  lp[, label := nice_term_v32(term)]
  lp[, label := factor(label, levels = unique(lp[population == "Early (primary)"][order(smd_raw), label]))]
  lp[, population := factor(population, levels = c(unname(pops), "Three-arm (max pairwise)"))]
  lpl <- melt(lp, id.vars = c("label", "population"), measure.vars = c("smd_raw", "smd_w"), variable.name = "sample", value.name = "smd")
  lpl[, sample := factor(ifelse(sample == "smd_raw", "Unweighted", "Weighted"), levels = c("Unweighted", "Weighted"))]
  ggplot2::ggplot(lpl, ggplot2::aes(smd, label, shape = sample)) + ggplot2::geom_vline(xintercept = 0.1, linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_point(size = 1.4) + ggplot2::scale_shape_manual(values = c(Unweighted = 1, Weighted = 16)) +
    ggplot2::facet_wrap(~population, nrow = 1) + ggplot2::labs(x = "Absolute standardized mean difference", y = NULL, shape = NULL) +
    theme_fig_v32() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6.5))
}
