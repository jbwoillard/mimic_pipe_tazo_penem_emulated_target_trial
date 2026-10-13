# =============================================================================
# 74_dag.R — DAG (A1) : figure principale (12 nœuds), version détaillée,
#   table nœud -> variables mesurées -> fenêtre -> dans le PS, ensemble
#   d'ajustement minimal (dagitty::adjustmentSets) et écarts avec le PS
# =============================================================================

DAG_NODES_V32 <- data.table::data.table(
  node = c("Treatment", "AKI", "Comorbidity", "BaselineRenal", "CreatTrajectory", "Severity", "InitContext",
           "CareUnit", "Infection", "Nephrotoxins", "Vancomycin", "Period",
           "PrescriberPreference", "Allergy", "Colonization"),
  label = c("PTZ vs\nmeropenem", "AKI\nby day 7", "Comorbidity &\ndemographics", "Baseline\nrenal function", "Pre-T0 creatinine\ntrajectory",
            "Acute severity", "Initiation context\n(timing, prior abx)", "Care unit\nat T0", "Infection site\n& cultures", "Nephrotoxic\nco-exposures",
            "Vancomycin", "Calendar\nperiod", "Prescriber\npreference*", "Beta-lactam\nallergy*", "Undocumented\ncolonization*"),
  measured = c(TRUE, TRUE, rep(TRUE, 10), FALSE, FALSE, FALSE),
  variables = c("treatment", "outcome_primary (strict KDIGO, day 7)",
                "age; gender; weight_kg_v31; ckd_priorhosp; diabetes_priorhosp; heart_failure_priorhosp; liver_disease_priorhosp; cancer_priorhosp; charlson_score_priorhosp; hematologic_malignancy; transplant_history",
                "baseline_creat_min7d; bun_baseline",
                "creat_last_pre_t0; creat_delta_last_min7d; creat_ratio_last_min7d; creat_slope_48h; creat_slope_48h_available; creat_single_pre_value",
                "sepsis_composite_v3; vasopressor_pre_t0; vent_pre_t0; map_mean; hr_mean; temp_max; spo2_min; gcs_min; lactate_max; platelets_min; bilirubin_max; neutropenia_pre_t0",
                "time_admission_to_T0_h (from ED/hospital arrival); prior_antibiotic; any_iv_abx_48h_pre (population restriction); true_new_user_90d (population restriction)",
                "service_group; icu_before_t0",
                "site_blood; site_urine; site_respiratory; site_abdominal; site_skin_soft_tissue_bone; site_other; positive_culture_pre_t0_storetime; esbl_resistant_history_12mo",
                "direct_nephrotoxin_pre_t0; calcineurin_inhibitor_pre_t0; loop_diuretic_pre_t0; acei_arb_pre_t0; nsaid_pre_t0; aminoglycoside_pre_t0; contrast_proxy_pre_t0 (sensitivity only)",
                "vanco_6h", "anchor_year_group", "not recorded in MIMIC-IV", "no allergy table in MIMIC-IV", "only documented resistance history is measured"),
  window = c("T0", "(T0; T0+7 d]", "prior hospitalizations; up to T0", "[T0-7 d; T0] / [T0-24 h; T0]", "[T0-7 d; T0] and [T0-48 h; T0]",
             "[T0-24 h; T0]", "arrival to T0; [T0-48 h; T0]", "at T0", "[T0-24 h; T0] (specimen drawn); 12 months (resistance)",
             "[T0-48 h; T0)", "[T0-24 h; T0+6 h]", "patient-level 3-year window", "-", "-", "-"))

DAG_EDGES_V32 <- c(
  "Comorbidity -> {BaselineRenal Severity Treatment AKI CareUnit}",
  "BaselineRenal -> {CreatTrajectory Treatment AKI}",
  "CreatTrajectory -> {Treatment AKI}",
  "Severity -> {CreatTrajectory Treatment AKI Vancomycin Nephrotoxins CareUnit}",
  "InitContext -> {Treatment CreatTrajectory Severity}",
  "CareUnit -> {Treatment AKI}",
  "Infection -> {Treatment Vancomycin Severity}",
  "Nephrotoxins -> {AKI Treatment}",
  "Vancomycin -> AKI",
  "Period -> {Treatment AKI}",
  "PrescriberPreference -> Treatment",
  "Allergy -> Treatment",
  "Colonization -> {Treatment Infection}",
  "Treatment -> AKI")

build_dag_v32 <- function() {
  lat <- DAG_NODES_V32[measured == FALSE, node]
  txt <- paste0("dag {\n", paste(paste0(lat, " [latent]"), collapse = "\n"), "\nTreatment [exposure]\nAKI [outcome]\n",
                paste(DAG_EDGES_V32, collapse = "\n"), "\n}")
  dagitty::dagitty(txt)
}

#' Ensembles d'ajustement minimaux et écart avec les nœuds couverts par le PS
dag_adjustment_v32 <- function(g, ps_covs) {
  sets <- dagitty::adjustmentSets(g, exposure = "Treatment", outcome = "AKI", type = "minimal")
  sets_l <- lapply(sets, as.character)
  # nœuds dont au moins une variable mesurée est dans le PS
  in_ps <- vapply(DAG_NODES_V32$node, function(n) {
    v <- trimws(unlist(strsplit(DAG_NODES_V32[node == n, variables], ";")))
    v <- sub(" .*", "", v); any(v %in% ps_covs)
  }, logical(1))
  covered <- names(in_ps)[in_ps]
  rows <- lapply(seq_along(sets_l), function(i) data.table::data.table(
    set_id = i, nodes = paste(sort(sets_l[[i]]), collapse = ", "),
    nodes_not_in_ps = paste(setdiff(sets_l[[i]], covered), collapse = ", "),
    all_in_ps = length(setdiff(sets_l[[i]], covered)) == 0))
  list(sets = data.table::rbindlist(rows), nodes_in_ps = covered,
       ps_nodes_not_required = setdiff(covered, unique(unlist(sets_l))))
}

#' Table nœud -> variables -> fenêtre -> dans le PS
dag_node_table_v32 <- function(ps_covs) {
  x <- copy(DAG_NODES_V32)
  x[, in_ps := vapply(variables, function(s) {
    v <- sub(" .*", "", trimws(unlist(strsplit(s, ";"))))
    hit <- v[v %in% ps_covs]
    if (length(hit) == 0) "no" else if (length(hit) == length(v)) "yes (all)" else paste0("yes (", length(hit), "/", length(v), ")")
  }, character(1))]
  x[node %in% c("Treatment", "AKI"), in_ps := "-"]
  x[, .(node = gsub("\n", " ", label), status = fifelse(measured, "measured", "UNMEASURED"), variables, window, in_ps)]
}

#' Figure du DAG (ggdag) : nœuds non mesurés en pointillés et texte gris
plot_dag_v32 <- function(g, detailed = FALSE) {
  set.seed(SEED_V3)
  coords <- data.table::data.table(
    name = c("Treatment", "AKI", "Comorbidity", "BaselineRenal", "CreatTrajectory", "Severity", "InitContext", "CareUnit",
             "Infection", "Nephrotoxins", "Vancomycin", "Period", "PrescriberPreference", "Allergy", "Colonization"),
    x = c(0, 7, 2, 3.6, 5.2, 3.6, 1.2, 5.6, 1.8, 5.4, 3.8, 3.2, -2, -2, -2),
    y = c(0, 0, 3.4, 3.4, 3.4, 1.8, 1.8, 1.8, -2, -2, -2.8, -1.2, 1.4, 0, -1.6))
  dagitty::coordinates(g) <- list(x = setNames(coords$x, coords$name), y = setNames(coords$y, coords$name))
  td <- ggdag::tidy_dagitty(g)
  df <- as.data.frame(td$data)
  df$latent <- df$name %in% DAG_NODES_V32[measured == FALSE, node]
  lab <- setNames(DAG_NODES_V32$label, DAG_NODES_V32$node)
  df$label <- lab[df$name]
  df$edge_type <- ifelse(df$latent, "dashed", "solid")
  df$role <- ifelse(df$name == "Treatment", "Exposure", ifelse(df$name == "AKI", "Outcome", ifelse(df$latent, "Unmeasured", "Measured confounder")))
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    ggdag::geom_dag_edges(ggplot2::aes(edge_linetype = edge_type), edge_width = 0.3, edge_colour = "grey30",
                          arrow_directed = grid::arrow(length = grid::unit(4, "pt"), type = "closed"),
                          start_cap = ggraph::circle(9, "mm"), end_cap = ggraph::circle(9, "mm")) +
    ggraph::scale_edge_linetype_manual(values = c(dashed = "dashed", solid = "solid"), guide = "none") +
    ggplot2::geom_label(ggplot2::aes(label = label, fill = role, colour = role), size = 2.8, label.r = grid::unit(2, "pt"),
                        linewidth = 0.3, lineheight = 0.9) +
    ggplot2::scale_fill_manual(values = c(Exposure = "#d9d9d9", Outcome = "#d9d9d9", `Measured confounder` = "white", Unmeasured = "white")) +
    ggplot2::scale_colour_manual(values = c(Exposure = "black", Outcome = "black", `Measured confounder` = "black", Unmeasured = "grey45")) +
    ggplot2::labs(caption = "Grey boxes: exposure and outcome. * Unmeasured in MIMIC-IV (grey text, dashed arrows).") +
    ggdag::theme_dag_blank(base_size = 9) + ggplot2::theme(legend.position = "none", plot.caption = ggplot2::element_text(size = 8))
}
