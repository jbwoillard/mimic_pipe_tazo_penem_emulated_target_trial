# =============================================================================
# 55_outcomes.R — AKI KDIGO stricte (docs/v3.1_plan.md A6, A7a, A8)
# =============================================================================
# Défaut corrigé (R/40_outcomes.R:rebuild_aki_outcomes_v3) :
#   - la référence du critère absolu (+0,3 mg/dL en 48 h) n'utilisait que des
#     valeurs POST-T0, et
#   - en l'absence de valeur antérieure dans les 48 h, retombait sur
#     baseline_creat_min7d.
# La première créatinine post-T0 était donc toujours jugée « +0,3 vs minimum
# sur 7 jours », ce qui n'est pas le critère KDIGO et est différentiel (les
# carbapénèmes, initiés plus tard dans le séjour, ont plus de valeurs
# pré-T0, donc un minimum plus bas).
#
# Définition v3.1 (KDIGO stricte), pour chaque créatinine j à t_j ∈ (T0 ; T0+H] :
#   critère absolu : ∃ i, t_i ∈ [t_j - 48 h ; t_j), t_i ≥ T0 - 48 h (valeurs
#                    pré-T0 incluses), Cr_j - Cr_i ≥ 0,3 mg/dL. Sans paire
#                    disponible : critère non rempli (pas de repli baseline).
#   critère relatif : Cr_j ≥ 1,5 × baseline_creat_min7d (min sur [T0-7 j ; T0]).
#   stade : 3 si Cr_j ≥ 3 × baseline, ou Cr_j ≥ 4,0 avec un critère rempli,
#           ou RRT initiée dans (T0 ; T0+H] ; 2 si Cr_j ≥ 2 × baseline ;
#           1 si un critère rempli.
#   évaluable : ≥ 1 créatinine dans (T0 ; T0+H], ou RRT dans (T0 ; T0+H].

ABS_THR_V31 <- 0.3
REL_THR_V31 <- 1.5

#' Évaluation mesure par mesure (fonction interne, testée directement)
#' @param creat  data.table(subject_id, hadm_id, charttime, creat_mgdl)
#' @param cohort data.table(subject_id, hadm_id, T0, baseline_creat_min7d)
#' @param horizon_h horizon post-T0 en heures
#' @return data.table une ligne par créatinine post-T0 avec abs/rel/stade
kdigo_measurements_v31 <- function(creat, cohort, horizon_h) {
  key <- cohort[, .(subject_id, hadm_id, T0, baseline = baseline_creat_min7d)]
  cr <- merge(creat[, .(subject_id, hadm_id, charttime, creat_mgdl)], key,
              by = c("subject_id", "hadm_id"))
  cr[, t := as.numeric(difftime(charttime, T0, units = "hours"))]
  cr <- cr[t >= -48 & t <= horizon_h & !is.na(creat_mgdl)]
  post <- cr[t > 0]
  if (nrow(post) == 0) return(post[, `:=`(ref48 = numeric(), abs_crit = integer(),
                                          rel_crit = integer(), stage = integer())])
  post[, `:=`(row_j = .I, lo = t - 48)]
  # référence 48 h glissante : min des valeurs (pré- ou post-T0) dans [t_j - 48 h ; t_j)
  ref <- cr[post, on = .(subject_id, hadm_id, t >= lo, t < t),
            .(row_j = i.row_j, cr_i = x.creat_mgdl), allow.cartesian = TRUE, nomatch = NULL]
  ref <- ref[, .(ref48 = min(cr_i)), by = row_j]
  post <- merge(post, ref, by = "row_j", all.x = TRUE)
  post[, abs_crit := as.integer(!is.na(ref48) & (creat_mgdl - ref48) >= ABS_THR_V31 - 1e-9)]
  post[, rel_crit := as.integer(!is.na(baseline) & creat_mgdl >= REL_THR_V31 * baseline - 1e-9)]
  post[, stage := fcase(
    !is.na(baseline) & creat_mgdl >= 3 * baseline - 1e-9, 3L,
    creat_mgdl >= 4.0 & (abs_crit == 1L | rel_crit == 1L), 3L,
    !is.na(baseline) & creat_mgdl >= 2 * baseline - 1e-9, 2L,
    abs_crit == 1L | rel_crit == 1L, 1L,
    default = 0L)]
  post[, c("lo", "row_j") := NULL]
  post[]
}

#' Outcomes patient-niveau sur un horizon donné
#' @param rrt_flag_col nom d'une colonne 0/1 de `cohort` : RRT initiée dans
#'                     (T0 ; T0+horizon] (NULL = ignorée)
compute_kdigo_v31 <- function(creat, cohort, horizon_d = 7, rrt_flag_col = NULL,
                               suffix = paste0("_d", horizon_d, "_v31")) {
  m <- kdigo_measurements_v31(creat, cohort, horizon_d * 24)
  agg <- m[, .(n_post = .N, abs_any = max(abs_crit), rel_any = max(rel_crit),
               stage_max = max(stage),
               t_first_aki = if (any(abs_crit == 1L | rel_crit == 1L))
                 min(t[abs_crit == 1L | rel_crit == 1L]) else NA_real_),
           by = .(subject_id, hadm_id)]
  out <- merge(cohort[, .(subject_id, hadm_id)], agg, by = c("subject_id", "hadm_id"), all.x = TRUE)
  rrt <- if (!is.null(rrt_flag_col)) as.integer(cohort[out, on = .(subject_id, hadm_id)][[rrt_flag_col]] == 1L) else rep(0L, nrow(out))
  rrt[is.na(rrt)] <- 0L
  evaluable <- (!is.na(out$n_post) & out$n_post > 0) | rrt == 1L
  aki <- as.integer((!is.na(out$abs_any) & (out$abs_any == 1L | out$rel_any == 1L)) | rrt == 1L)
  stg <- pmax(fifelse(is.na(out$stage_max), 0L, out$stage_max), 3L * rrt)
  trig <- fcase(rrt == 1L & !(out$abs_any %in% 1L) & !(out$rel_any %in% 1L), "rrt_only",
                out$abs_any %in% 1L & out$rel_any %in% 1L, "both",
                out$abs_any %in% 1L, "abs_0.3_48h_only",
                out$rel_any %in% 1L, "rel_1.5x_only",
                default = "none")
  res <- data.table(subject_id = out$subject_id, hadm_id = out$hadm_id)
  res[, paste0("aki", suffix)          := fifelse(evaluable, aki, NA_integer_)]
  res[, paste0("aki_abs48h", suffix)   := fifelse(evaluable, as.integer(out$abs_any %in% 1L), NA_integer_)]
  res[, paste0("aki_rel15", suffix)    := fifelse(evaluable, as.integer(out$rel_any %in% 1L), NA_integer_)]
  res[, paste0("aki_stage_max", suffix) := fifelse(evaluable, stg, NA_integer_)]
  res[, paste0("aki_stage2", suffix)   := fifelse(evaluable, as.integer(stg >= 2L), NA_integer_)]
  res[, paste0("aki_stage3", suffix)   := fifelse(evaluable, as.integer(stg >= 3L), NA_integer_)]
  res[, paste0("aki_trigger", suffix)  := fifelse(evaluable, trig, NA_character_)]
  res[, paste0("t_first_aki_h", suffix) := out$t_first_aki]
  res[, paste0("n_creat_post", suffix) := fifelse(is.na(out$n_post), 0L, out$n_post)]
  res
}

#' AKI persistante (A8) : parmi les AKI de l'horizon primaire, la première
#' créatinine mesurée ≥ `persist_h` après le premier critère rempli satisfait
#' encore Cr ≥ 1,5 × baseline ou Cr ≥ baseline + 0,3. Codage explicite :
#'   0  = pas d'AKI (outcome primaire évaluable et négatif) ;
#'   1  = AKI persistante ;
#'   0  = AKI résolue (première valeur ≥ persist_h sous les deux seuils) ;
#'   NA = AKI sans aucune créatinine ≥ persist_h après le début (sortie,
#'        décès précoce, arrêt des dosages) -> exclu de cette analyse, N rapporté.
#' Recherche des valeurs de suivi jusqu'à T0 + 14 j.
compute_aki_persistence_v31 <- function(creat, cohort, aki_col, t_first_col, persist_h) {
  key <- cohort[, .(subject_id, hadm_id, T0, baseline = baseline_creat_min7d,
                    aki = get(aki_col), t1 = get(t_first_col))]
  cr <- merge(creat[, .(subject_id, hadm_id, charttime, creat_mgdl)], key[aki == 1L],
              by = c("subject_id", "hadm_id"))
  cr[, t := as.numeric(difftime(charttime, T0, units = "hours"))]
  nxt <- cr[t >= t1 + persist_h & t <= 14 * 24][order(subject_id, hadm_id, t)][
    , .SD[1], by = .(subject_id, hadm_id)]
  nxt[, persist := as.integer(creat_mgdl >= REL_THR_V31 * baseline - 1e-9 |
                                creat_mgdl - baseline >= ABS_THR_V31 - 1e-9)]
  res <- merge(key[, .(subject_id, hadm_id, aki)], nxt[, .(subject_id, hadm_id, persist)],
               by = c("subject_id", "hadm_id"), all.x = TRUE)
  res[, out := fcase(is.na(aki), NA_integer_, aki == 0L, 0L, !is.na(persist), persist,
                     default = NA_integer_)]
  res[, .(subject_id, hadm_id, out)]
}

#' Contrôle négatif (A7a) : hausse précoce de créatinine dans (T0 ; T0+12 h],
#' même règle KDIGO stricte restreinte à 12 h. Évaluable si ≥ 1 créatinine
#' dans (T0 ; T0+12 h]. Trop précoce pour une toxicité tubulaire : un effet
#' non nul signale une confusion par la trajectoire de créatinine.
compute_early_creat_rise_v31 <- function(creat, cohort, window_h = 12) {
  m <- kdigo_measurements_v31(creat, cohort, window_h)
  agg <- m[, .(early_creat_rise_12h = as.integer(any(abs_crit == 1L | rel_crit == 1L))),
           by = .(subject_id, hadm_id)]
  merge(cohort[, .(subject_id, hadm_id)], agg, by = c("subject_id", "hadm_id"), all.x = TRUE)
}
