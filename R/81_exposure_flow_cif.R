# =============================================================================
# 81_exposure_flow_cif.R — Phase B v3.3 (docs/v3.3_plan.md) : exposition et
#   microbiologie (B2), flux STROBE par bras (B3), critère KDIGO déclencheur
#   (B4), incidence cumulée avec bandes d'IC et effectifs à risque (B5)
# =============================================================================

# -----------------------------------------------------------------------
# B2 — doses (emar_detail), exposition, vancomycine, microbiologie
# -----------------------------------------------------------------------
#' Doses administrées de PTZ et méropénème (emar JOIN emar_detail), DuckDB,
#' restreint aux patients de la cohorte ; cache rds
extract_index_doses_v33 <- function(subject_ids, mimic_root, cache) {
  if (file.exists(cache)) return(readRDS(cache))
  con <- duckdb::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(duckdb::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(con, "ids", data.frame(subject_id = as.character(unique(subject_ids))))
  f <- function(x) normalizePath(file.path(mimic_root, "hosp", x))
  q <- sprintf("
    WITH e AS (SELECT subject_id, hadm_id, emar_id, charttime, medication, event_txt
               FROM read_csv('%s', all_varchar = true)
               WHERE subject_id IN (SELECT subject_id FROM ids)
                 AND (lower(medication) LIKE '%%piperacillin%%' OR lower(medication) LIKE '%%meropenem%%'))
    SELECT e.*, d.parent_field_ordinal, d.dose_given, d.dose_given_unit, d.product_description
    FROM e JOIN read_csv('%s', all_varchar = true) d ON e.emar_id = d.emar_id
    WHERE d.dose_given IS NOT NULL", f("emar.csv.gz"), f("emar_detail.csv.gz"))
  x <- data.table::as.data.table(DBI::dbGetQuery(con, q))
  x[, `:=`(subject_id = as.integer(subject_id), hadm_id = as.integer(hadm_id), charttime = lubridate::ymd_hms(charttime, quiet = TRUE),
           dose_num = suppressWarnings(as.numeric(gsub(",", "", dose_given))))]
  saveRDS(x, cache); x
}

#' Dose en grammes par administration (une valeur par emar_id : le maximum des
#' lignes de détail exprimées en g ou mg ; PTZ en g de pipéracilline-tazobactam)
dose_grams_v33 <- function(x) {
  x <- x[!is.na(dose_num) & tolower(event_txt) %in% c("administered", "delayed administered", "started", "delayed started", "restarted", "new bag", "applied", "confirmed")]
  x[, g := fcase(tolower(dose_given_unit) == "g", dose_num, tolower(dose_given_unit) == "mg", dose_num / 1000, default = NA_real_)]
  x[!is.na(g), .(g = max(g)), by = .(subject_id, hadm_id, emar_id, charttime, drug = fifelse(grepl("piperacillin", tolower(medication)), "ptz", "meropenem"))]
}

NARROWER_REGEX_V33 <- "cefazolin|ceftriaxone|cefotaxime|ampicillin|amoxicillin|oxacillin|nafcillin|cephalexin|cefpodoxime|cefuroxime|penicillin|sulfameth|(cipro|levo|moxi)floxacin"

#' Description de l'exposition, par patient
#' @param cohort data.table (subject_id, hadm_id, T0, treatment, dischtime, deathtime)
exposure_patient_v33 <- function(cohort, ab, doses) {
  key <- cohort[, .(subject_id, hadm_id, T0, treatment, end_t = fcoalesce(deathtime, dischtime))]
  a <- merge(ab, key, by = c("subject_id", "hadm_id"))[, h := as.numeric(difftime(t, T0, units = "hours"))][h >= 0]
  own <- a[iv == 1L & ((treatment == "PTZ" & abx_class == "ptz") | (treatment == "meropenem" & abx_class == "meropenem"))][order(hadm_id, h)]
  # cycle initial : se termine à la première dose suivie d'un intervalle > 24 h
  crs <- own[, { hh <- c(h, Inf); i <- which(diff(hh) > 24)[1]; last <- h[i]
                 .(n_doses_course = i, last_dose_h = last, gap_after = hh[i + 1] - last) }, by = hadm_id]
  p <- merge(key, crs, by = "hadm_id", all.x = TRUE)
  p[, follow_h := as.numeric(difftime(end_t, T0, units = "hours"))]
  p[, `:=`(duration_days = last_dose_h / 24,
           stop_lt48h = as.integer(last_dose_h < 48 & (is.na(follow_h) | follow_h > last_dose_h + 24)))]
  # désescalade : début d'un agent plus étroit entre -12 h et +24 h autour de la dernière dose du cycle
  nar <- a[abx_class == "other_antibacterial" & grepl(NARROWER_REGEX_V33, tolower(medication))]
  nar <- merge(nar[, .(hadm_id, h)], p[, .(hadm_id, last_dose_h)], by = "hadm_id")
  de <- nar[h >= last_dose_h - 12 & h <= last_dose_h + 24, .(deescalation = 1L), by = hadm_id]
  # passage vers l'autre bêta-lactamine (PTZ <-> méropénème) ou le céfépime après T0 + 24 h
  sw <- a[iv == 1L & h > 24 & ((treatment == "PTZ" & abx_class %in% c("meropenem", "cefepime")) |
                                 (treatment == "meropenem" & abx_class %in% c("ptz", "cefepime"))), .(switch_other_bl = 1L), by = hadm_id]
  # vancomycine IV J0-J3 (jour = floor(h / 24))
  vc <- a[abx_class == "vancomycin" & iv == 1L & h < 96, .(vanco_days_d0_d3 = uniqueN(floor(h / 24))), by = hadm_id]
  # dose des 24 premières heures et du jour 1 complet (24-48 h)
  dg <- merge(doses, key[, .(hadm_id, T0, treatment)], by = "hadm_id")[, h := as.numeric(difftime(charttime, T0, units = "hours"))]
  dg <- dg[(treatment == "PTZ" & drug == "ptz") | (treatment == "meropenem" & drug == "meropenem")]
  d24 <- dg[h >= 0 & h < 24, .(dose_g_first24h = sum(g)), by = hadm_id]
  d1 <- dg[h >= 24 & h < 48, .(dose_g_day1 = sum(g)), by = hadm_id]
  for (x in list(de, sw, vc, d24, d1)) p <- merge(p, x, by = "hadm_id", all.x = TRUE)
  p[is.na(deescalation), deescalation := 0L]; p[is.na(switch_other_bl), switch_other_bl := 0L]; p[is.na(vanco_days_d0_d3), vanco_days_d0_d3 := 0L]
  p[, vanco_d0_d3 := as.integer(vanco_days_d0_d3 > 0)]
  p[, high_dose := fifelse(is.na(dose_g_day1), NA_integer_, as.integer((treatment == "PTZ" & dose_g_day1 >= 13.5) | (treatment == "meropenem" & dose_g_day1 >= 3)))]
  p
}

.wq <- function(x, w) { o <- order(x); x <- x[o]; w <- w[o]; cw <- cumsum(w) / sum(w); vapply(c(.5, .25, .75), function(q) x[which(cw >= q)[1]], numeric(1)) }
summ_exposure_v33 <- function(p, w) {
  rows <- list(
    c("Duration of the initial course, days", "median [IQR]", "duration_days", "cont"),
    c("Dose, first 24 h, g", "median [IQR]", "dose_g_first24h", "cont"),
    c("Dose, day 1 (24-48 h), g", "median [IQR]", "dose_g_day1", "cont"),
    c("High day-1 dose (PTZ >=13.5 g; meropenem >=3 g)", "%", "high_dose", "bin"),
    c("Initial course ended within 48 h", "%", "stop_lt48h", "bin"),
    c("De-escalation to a narrower agent", "%", "deescalation", "bin"),
    c("Switch to the other study beta-lactam or cefepime after 24 h", "%", "switch_other_bl", "bin"),
    c("IV vancomycin, days 0-3", "%", "vanco_d0_d3", "bin"),
    c("Days with IV vancomycin, days 0-3 (exposed)", "median [IQR]", "vanco_days_d0_d3", "cont_exposed"))
  rbindlist(lapply(rows, function(r) rbindlist(lapply(c("PTZ", "meropenem"), function(arm) {
    s <- p$treatment == arm & !is.na(p[[r[3]]]); if (r[4] == "cont_exposed") s <- s & p$vanco_d0_d3 == 1L
    x <- p[[r[3]]][s]; ww <- w[s]
    f <- function(x, ww) if (r[4] == "bin") sprintf("%.1f", 100 * sum(ww * x) / sum(ww)) else { q <- .wq(x, ww); sprintf("%.1f [%.1f; %.1f]", q[1], q[2], q[3]) }
    data.table(characteristic = r[1], stat = r[2], variable = r[3], arm = arm, n_available = sum(s),
               unweighted = f(x, rep(1, length(x))), overlap_weighted = f(x, ww))
  }))))
}

#' Microbiologie : cultures prélevées dans [T0 - 24 h ; T0 + 72 h]
ENTEROBACTERALES_REGEX_V33 <- "ESCHERICHIA|KLEBSIELLA|ENTEROBACTER|PROTEUS|SERRATIA|CITROBACTER|MORGANELLA|PROVIDENCIA|HAFNIA|KLUYVERA|RAOULTELLA|SALMONELLA|PANTOEA"
microbiology_v33 <- function(cohort, tbl_paths, w) {
  mic <- read_mimic_table(tbl_paths$paths[["microbiologyevents"]], col_select = c("subject_id", "charttime", "chartdate", "spec_type_desc", "org_name"))
  mic[, `:=`(subject_id = as.integer(subject_id), t = fcoalesce(lubridate::ymd_hms(charttime, quiet = TRUE), lubridate::ymd_hms(chartdate, quiet = TRUE)))]
  m <- merge(mic, cohort[, .(subject_id, hadm_id, T0, treatment)], by = "subject_id")
  m[, h := as.numeric(difftime(t, T0, units = "hours"))]
  m <- m[h >= -24 & h <= 72]
  pos <- m[!is.na(org_name) & org_name != "" & !grepl("CANCELLED|NO GROWTH", org_name)]
  per <- cohort[, .(hadm_id, treatment)][, `:=`(w = w, cultured = as.integer(hadm_id %in% m$hadm_id), positive = as.integer(hadm_id %in% pos$hadm_id),
                                               pseudomonas = as.integer(hadm_id %in% pos[grepl("PSEUDOMONAS AERUGINOSA", org_name), hadm_id]),
                                               enterobacterales = as.integer(hadm_id %in% pos[grepl(ENTEROBACTERALES_REGEX_V33, org_name), hadm_id]))]
  summ <- melt(per, id.vars = c("hadm_id", "treatment", "w"), variable.name = "indicator")[
    , .(n = .N, n_yes = sum(value), pct_unweighted = 100 * mean(value), pct_overlap_weighted = 100 * sum(w * value) / sum(w)), by = .(indicator, treatment)]
  orgs <- unique(pos[, .(hadm_id, treatment, org_name)])[, .(n_patients = .N), by = .(treatment, org_name)]
  orgs <- merge(orgs, per[, .(n_arm = .N), by = treatment], by = "treatment")[, pct := 100 * n_patients / n_arm]
  top <- orgs[, .(tot = sum(n_patients)), by = org_name][order(-tot)][1:10, org_name]
  list(summary = summ, top = dcast(orgs[org_name %in% top], org_name ~ treatment, value.var = c("n_patients", "pct"), fill = 0)[order(-(n_patients_PTZ + n_patients_meropenem))])
}

# -----------------------------------------------------------------------
# B3 — flux STROBE par bras
# -----------------------------------------------------------------------
ARM_ORDER_V33 <- c("PTZ", "meropenem", "cefepime", "carbapenem_other")
flow_main_steps_v33 <- function(fla, master, reg_ids) {
  w <- dcast(fla, step ~ treatment, value.var = "n", fill = 0)
  build <- data.table(step = c("all", "single_arm", "first_episode", "adult", "baseline", "no_rrt", "no_esrd"),
                      label = c("Hospital stays with a first IV administration of PTZ, cefepime, meropenem or another carbapenem",
                                "Single study drug at initiation", "First episode per patient", "Adult (>=18 years)",
                                "Baseline creatinine available (7-day minimum)", "No RRT before T0", "No ESRD = unified cohort"),
                      exclusion = c("", "Co-initiation of another study drug within 24 h, or another carbapenem", "Later episodes of the same patient",
                                    "Age <18 years", "No creatinine in the 7 days before T0", "RRT before T0", "End-stage renal disease"))
  b <- merge(build, w, by = "step", sort = FALSE)
  m <- master
  seq_steps <- list(
    list("no_aki_t0", "No AKI at T0", "AKI at T0", quote(aki_at_t0 == 0L)),
    list("ptz_or_mero", "PTZ or meropenem initiators", "Cefepime initiators (three-arm and ACORN analyses only)", quote(treatment %in% c("PTZ", "meropenem"))),
    list("t0_24h", "T0 within 24 h of arrival", "T0 more than 24 h after arrival", quote(h_arrival_to_T0 <= 24)),
    list("new_user", "New user (no anti-pseudomonal beta-lactam in 90 days)", "Anti-pseudomonal beta-lactam in the prior 90 days", quote(true_new_user_90d == 1L)),
    list("no_iv_48h", "No IV antibiotic in the prior 48 h", "IV antibiotic in the prior 48 h", quote(any_iv_abx_48h_pre == 0L)),
    list("evaluable", "AKI status by day 7 evaluable = EARLY POPULATION (primary)", "AKI status by day 7 not evaluable", quote(!is.na(outcome_primary))))
  cur <- m[J(reg_ids[["all"]])]; out <- list()
  for (s in seq_steps) { cur <- cur[eval(s[[4]], cur) %in% TRUE]
    out[[s[[1]]]] <- data.table(step = s[[1]], label = s[[2]], exclusion = s[[3]], PTZ = sum(cur$treatment == "PTZ"), meropenem = sum(cur$treatment == "meropenem"),
                                cefepime = sum(cur$treatment == "cefepime"), carbapenem_other = 0L) }
  f <- rbind(b, rbindlist(out), fill = TRUE)
  for (k in ARM_ORDER_V33) { if (!k %in% names(f)) f[, (k) := 0L]; set(f, which(is.na(f[[k]])), k, 0L) }
  f[, n := PTZ + meropenem + cefepime + carbapenem_other]
  for (k in c("n", ARM_ORDER_V33)) f[, (paste0("excl_", k)) := c(NA, -diff(get(k)))]
  f[]
}

fmt_arms_v33 <- function(r, pre = "") {
  parts <- c(PTZ = "PTZ", meropenem = "MEM", cefepime = "FEP", carbapenem_other = "other CBP")
  v <- vapply(names(parts), function(k) r[[paste0(pre, k)]], numeric(1))
  paste(sprintf("%s %s", parts[v > 0], format(v[v > 0], big.mark = ",", trim = TRUE)), collapse = " / ")
}

fig_flow_detailed_v33 <- function(f, branches) {
  nn <- nrow(f); f[, y := rev(seq_len(nn)) * 2.2]
  f[, box := sprintf("%s\nN = %s (%s)", label, format(n, big.mark = ",", trim = TRUE), vapply(seq_len(.N), function(i) fmt_arms_v33(f[i]), ""))]
  ex <- f[!is.na(excl_n) & excl_n > 0][, ex_lab := sprintf("Excluded n = %s: %s\n(%s)", format(excl_n, big.mark = ",", trim = TRUE), exclusion,
                                                         vapply(seq_len(.N), function(i) fmt_arms_v33(.SD[i], "excl_"), ""))]
  ex[, y_ex := y + 1.1]
  br <- merge(branches, f[, .(step, y)], by.x = "from", by.y = "step")
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = f[-.N], ggplot2::aes(x = 0, xend = 0, y = y - 0.55, yend = y - 1.65), arrow = grid::arrow(length = grid::unit(1.6, "mm"), type = "closed"), linewidth = 0.3) +
    ggplot2::geom_segment(data = ex, ggplot2::aes(x = 0, xend = 1.55, y = y_ex, yend = y_ex), arrow = grid::arrow(length = grid::unit(1.4, "mm"), type = "closed"), linewidth = 0.3) +
    ggplot2::geom_segment(data = br, ggplot2::aes(x = 0, xend = -1.55, y = y, yend = y), linetype = "dashed", arrow = grid::arrow(length = grid::unit(1.4, "mm"), type = "closed"), linewidth = 0.3) +
    ggplot2::geom_label(data = f, ggplot2::aes(0, y, label = box), size = 2.1, label.size = 0.3, lineheight = 0.95, fill = ifelse(f$step == "evaluable", "grey88", "white")) +
    ggplot2::geom_label(data = ex, ggplot2::aes(1.6, y_ex, label = ex_lab), size = 1.85, hjust = 0, label.size = 0.2, lineheight = 0.95, colour = "grey20") +
    ggplot2::geom_label(data = br, ggplot2::aes(-1.6, y, label = label), size = 1.85, hjust = 1, label.size = 0.2, lineheight = 0.95, fill = "grey96") +
    ggplot2::scale_x_continuous(limits = c(-4.2, 5.6)) + ggplot2::theme_void() +
    ggplot2::labs(caption = "PTZ, piperacillin-tazobactam; MEM, meropenem; FEP, cefepime; CBP, carbapenem. Dashed arrows: secondary, design and benchmark populations.") +
    ggplot2::theme(plot.caption = ggplot2::element_text(size = 6, hjust = 0), plot.margin = ggplot2::margin(4, 4, 4, 4))
}

# -----------------------------------------------------------------------
# B4 — critère KDIGO déclencheur et AKI confirmée
# -----------------------------------------------------------------------
aki_trigger_table_v33 <- function(raw, w) {
  d <- raw[, .(treatment, trig = aki_trigger_d7_v31, aki = outcome_primary, pers = aki_persistent_48h_v31)][, w := w]
  lv <- c(abs_0.3_48h_only = "+0.3 mg/dL within 48 h only", rel_1.5x_only = ">=1.5x the 7-day minimum only", both = "Both creatinine criteria", rrt_only = "RRT only")
  a <- d[aki == 1L]
  t1 <- a[, .(n = .N, pct_unweighted = 100 * .N / nrow(a[treatment == .BY$treatment]), pct_overlap_weighted = 100 * sum(w) / sum(a[treatment == .BY$treatment, w])), by = .(treatment, trig)]
  t1[, criterion := lv[trig]][, trig := factor(trig, levels = names(lv))]
  conf <- a[, .(n = sum(pers == 1L, na.rm = TRUE), n_evaluable = sum(!is.na(pers)), pct_unweighted = 100 * mean(pers == 1L, na.rm = TRUE),
                pct_overlap_weighted = 100 * sum((w * (pers == 1L))[!is.na(pers)]) / sum(w[!is.na(pers)])), by = treatment][, criterion := "Confirmed AKI (persistent >=48 h) among AKI"]
  inc <- d[, .(n = sum(aki), n_evaluable = .N, pct_unweighted = 100 * mean(aki), pct_overlap_weighted = 100 * sum(w * aki) / sum(w)), by = treatment][, criterion := "AKI by day 7 (all criteria)"]
  rbind(inc, t1[order(treatment, trig)][, trig := NULL], conf, fill = TRUE)[, .(treatment, criterion, n, n_evaluable, pct_unweighted, pct_overlap_weighted)]
}

# -----------------------------------------------------------------------
# B5 — incidence cumulée avec bandes d'IC et effectifs à risque
# -----------------------------------------------------------------------
cif_bands_v33 <- function(d, covs, method, B = 200, seed = SEED_V3, grid = seq(0, 14, by = 0.25), at_risk_days = c(0, 2, 3, 7, 14)) {
  X0 <- build_ps_design_v31(d, covs); a0 <- as.integer(d$treatment == "PTZ")
  f0 <- fit_ps_v31(X0, a0, method); lam <- if (method == "ridge") f0$lambda else NULL
  wfun <- function(dd) { X <- build_ps_design_v31(dd, covs); a <- as.integer(dd$treatment == "PTZ")
    ps <- fit_ps_v31(X, a, method, lambda = lam)$ps; ifelse(a == 1, 1 - ps, ps) }
  w0 <- ifelse(a0 == 1, 1 - f0$ps, f0$ps)
  pt <- .cif_fit_v32(d, w0, grid)$points
  set.seed(seed)
  bt <- rbindlist(lapply(seq_len(B), function(b) { i <- sample.int(nrow(d), replace = TRUE); db <- d[i]
    r <- tryCatch(.cif_fit_v32(db, wfun(db), grid)$points, error = function(e) NULL); if (is.null(r)) NULL else r[, rep := b] }))
  band <- bt[, .(lo = stats::quantile(cif_aki, 0.025), hi = stats::quantile(cif_aki, 0.975)), by = .(arm, day)]
  curves <- merge(pt[, .(arm, day, cif_aki)], band, by = c("arm", "day"))
  st <- .cif_status_v32(d)[, arm := ifelse(d$treatment == "PTZ", "PTZ", "comparator")]
  ar <- rbindlist(lapply(at_risk_days, function(t) st[, .(day = t, n_at_risk = sum(time >= t | (t == 0))), by = arm]))
  list(curves = curves, at_risk = ar, B = B, n_failed = B - uniqueN(bt$rep))
}

fig_cif_bands_v33 <- function(cl) {
  x <- rbindlist(lapply(names(cl), function(p) cl[[p]]$curves[, population := p]))
  ar <- rbindlist(lapply(names(cl), function(p) cl[[p]]$at_risk[, population := p]))
  lv <- intersect(c("Early (primary)", "Broad (secondary)"), names(cl))
  for (z in list(x, ar)) { z[, arm := factor(ifelse(arm == "PTZ", "PTZ", "Meropenem"), levels = c("PTZ", "Meropenem"))]; z[, population := factor(population, levels = lv)] }
  xs <- ggplot2::scale_x_continuous(breaks = c(0, 2, 3, 7, 14), limits = c(-0.6, 14.6), expand = c(0, 0))
  p1 <- ggplot2::ggplot(x, ggplot2::aes(day, 100 * cif_aki)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 100 * lo, ymax = 100 * hi, fill = arm), alpha = 0.18, stat = "identity") +
    ggplot2::geom_step(ggplot2::aes(linetype = arm), linewidth = 0.55) + ggplot2::facet_wrap(~population, nrow = 1) +
    ggplot2::scale_fill_manual(values = c(PTZ = "grey20", Meropenem = "grey55")) + xs +
    ggplot2::labs(x = NULL, y = "Cumulative incidence of AKI, %", linetype = NULL, fill = NULL) + theme_fig_v32() +
    ggplot2::theme(legend.position = "top")
  # effectifs à risque : J2 et J3 décalés verticalement pour éviter le chevauchement
  ar[, yy := as.numeric(arm) + fifelse(day == 2, 0.22, fifelse(day == 3, -0.22, 0))]
  p2 <- ggplot2::ggplot(ar, ggplot2::aes(day, yy, label = format(n_at_risk, big.mark = ",", trim = TRUE))) + ggplot2::geom_text(size = 1.9) +
    ggplot2::facet_wrap(~population, nrow = 1) + xs +
    ggplot2::scale_y_continuous(breaks = 1:2, labels = levels(ar$arm), limits = c(0.5, 2.5)) +
    ggplot2::labs(x = "Days after T0", y = "No. at risk") +
    theme_fig_v32() + ggplot2::theme(strip.text = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), axis.line.y = ggplot2::element_blank(),
                                     axis.ticks.y = ggplot2::element_blank())
  patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(4, 1.1)) +
    patchwork::plot_annotation(caption = paste("Overlap-weighted Aalen-Johansen estimator, death as a competing risk, censoring at discharge;",
                                               sprintf("shaded bands: pointwise 95%% bootstrap CI (B = %d). Numbers at risk are unweighted.", cl[[1]]$B), sep = "\n"),
                               theme = ggplot2::theme(plot.caption = ggplot2::element_text(size = 7, hjust = 0)))
}

# -----------------------------------------------------------------------
# C3 — figure des trajectoires (3 marqueurs) : libellés sur deux lignes,
#   effectifs placés sous les barres d'erreur (espace réservé par panneau)
# -----------------------------------------------------------------------
MARKER_SHORT_V33 <- c("Creatinine (mg/dL)" = "Creatinine,\nmg/dL", "Urea nitrogen (mg/dL)" = "Urea nitrogen,\nmg/dL",
                      "Urea-to-creatinine ratio" = "Urea-to-\ncreatinine ratio")
fig_trajectory_v33 <- function(tr) {
  tr <- copy(tr)[, marker := factor(MARKER_SHORT_V33[as.character(marker)], levels = unname(MARKER_SHORT_V33))]
  ab <- rbind(tr[, .(marker, population, day, arm = "PTZ", value = mean_ptz)], tr[, .(marker, population, day, arm = "Meropenem", value = mean_comp)])
  ab[, arm := factor(arm, levels = c("PTZ", "Meropenem"))]
  p1 <- ggplot2::ggplot(ab, ggplot2::aes(day, value, linetype = arm, shape = arm)) +
    ggplot2::geom_line(linewidth = 0.5) + ggplot2::geom_point(size = 1.6) +
    ggplot2::facet_grid(marker ~ population, scales = "free_y") + ggplot2::scale_x_continuous(breaks = 0:7) +
    ggplot2::labs(x = NULL, y = "Overlap-weighted mean", linetype = NULL, shape = NULL) + theme_fig_v32() +
    ggplot2::theme(legend.position = "top", strip.text.y = ggplot2::element_text(angle = 0, hjust = 0))
  rng <- tr[, .(lo = min(diff_lo, 0, na.rm = TRUE), hi = max(diff_hi, 0, na.rm = TRUE)), by = marker][, r := hi - lo]
  tr <- merge(tr, rng, by = "marker")[, `:=`(lab_y = lo - 0.14 * r, blank_y = lo - 0.30 * r)]
  p2 <- ggplot2::ggplot(tr, ggplot2::aes(day, diff)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = diff_lo, ymax = diff_hi), width = 0.2, linewidth = 0.4) +
    ggplot2::geom_point(size = 1.6) + ggplot2::geom_blank(ggplot2::aes(y = blank_y)) +
    ggplot2::geom_text(ggplot2::aes(y = lab_y, label = paste0(format(n_ptz, big.mark = ",", trim = TRUE), "\n", format(n_comp, big.mark = ",", trim = TRUE))),
                       size = 1.75, lineheight = 0.85, colour = "grey30") +
    ggplot2::facet_grid(marker ~ population, scales = "free_y") + ggplot2::scale_x_continuous(breaks = 0:7) +
    ggplot2::labs(x = "Day after T0", y = "PTZ minus meropenem (95% CI)",
                  caption = "Numbers under each panel: patients contributing a value that day (PTZ above, meropenem below).") +
    theme_fig_v32() + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0, hjust = 0))
  patchwork::wrap_plots(p1, p2, ncol = 1)
}
