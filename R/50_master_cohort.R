# =============================================================================
# 50_master_cohort.R — Objet cohorte unique v3.1 (docs/v3.1_plan.md A1)
# =============================================================================
# `master` = cohorte ÉLIGIBLE v3 (après exclusion de l'AKI à T0, avant
# restriction à l'outcome évaluable), une ligne par épisode, enrichie :
#   - des constantes / biologies ré-extraites, NON imputées (A4, A5) ;
#   - des nouvelles covariables (A11) et du service regroupé (A3) ;
#   - des outcomes KDIGO stricts (A6), de persistance (A8) et des contrôles
#     négatifs (A7) ;
#   - de l'heure d'arrivée (urgences ou admission) pour la phase B.
# Toute cohorte d'analyse en dérive via get_analysis_cohort() (R/51_).

#' Variables de covariables du PS primaire v3.1. Changements vs v3 :
#'  - infection_site_proxy (98 % « aucune culture positive connue à T0 »)
#'    remplacé par les indicateurs de site des cultures PRÉLEVÉES pré-T0 ;
#'  - service_at_t0 (40 careunit bruts) remplacé par service_group (7 niveaux) ;
#'  - ajout : platelets_min, bilirubin_max (composantes SOFA non mesurées
#'    ailleurs), weight_kg_v31.
PS_COVARIATES_V31 <- c(
  "age", "gender", "baseline_creat_min7d", "bun_baseline",
  "ckd_priorhosp", "diabetes_priorhosp", "heart_failure_priorhosp", "liver_disease_priorhosp",
  "cancer_priorhosp", "charlson_score_priorhosp", "hematologic_malignancy", "transplant_history",
  "vasopressor_pre_t0", "vent_pre_t0", "map_mean", "hr_mean", "temp_max", "spo2_min", "gcs_min",
  "sepsis_composite_v3", "icu_before_t0", "lactate_max", "neutropenia_pre_t0",
  "platelets_min", "bilirubin_max", "weight_kg_v31",
  "direct_nephrotoxin_pre_t0", "calcineurin_inhibitor_pre_t0", "loop_diuretic_pre_t0",
  "acei_arb_pre_t0", "nsaid_pre_t0", "aminoglycoside_pre_t0",
  "vanco_6h", "positive_culture_pre_t0_storetime", "prior_antibiotic", "esbl_resistant_history_12mo",
  "site_blood", "site_urine", "site_respiratory", "site_abdominal", "site_skin_soft_tissue_bone", "site_other",
  "time_admission_to_T0_h", "anchor_year_group", "service_group"
)

#' Variables continues avec manquants réels (imputées) — pour les
#' indicateurs *_missing de la sensibilité « médiane + indicateurs »
VARS_WITH_MISSING_V31 <- c("hr_mean", "map_mean", "temp_max", "spo2_min", "gcs_min", "bun_baseline",
                           "lactate_max", "platelets_min", "bilirubin_max", "weight_kg_v31",
                           "baseline_creat_min7d")

build_master_cohort_v31 <- function(cohort_v3, ce, lab, creat_dt, tbl_paths,
                                     out_tables = "output/extraction/tables") {
  m <- copy(cohort_v3)
  key <- m[, .(subject_id, hadm_id, T0, treatment)]

  # --- 0. conserver l'outcome v3 (sensibilité « ancienne définition ») et
  #        retirer les colonnes imputées/erronées héritées
  m[, aki_d7_v3def := outcome_primary_v3]
  drop <- c("hr_mean", "hr_max", "map_mean", "map_min", "temp_max", "temp_min", "spo2_min",
            "gcs_min", "fio2_max", "bun_baseline", "lactate_max", grep("_missing$", names(m), value = TRUE))
  # infection_site_proxy (v3) conservée uniquement pour l'étape 0-1 du waterfall
  m[, (intersect(drop, names(m))) := NULL]

  # --- 1. constantes et biologies pré-T0 ré-extraites (A4, A5)
  pre <- summarise_pre_t0_v31(ce, lab, key)
  m <- merge(m, pre$cov, by = c("subject_id", "hadm_id"), all.x = TRUE)
  fwrite(pre$range_report, file.path(out_tables, "physio_range_report_v31.csv"))

  # --- 2. poids : réa (chartevents) sinon omr <= 365 j ; bornes physiologiques
  m[, weight_omr := fifelse(.in_range(weight_kg, "weight"), weight_kg, NA_real_)]
  m[, weight_kg_v31 := fcoalesce(weight_icu_kg, weight_omr)]
  m[, bmi_v31 := fifelse(.in_range(bmi, "bmi"), bmi, NA_real_)]
  m[, c("weight_omr") := NULL]

  # --- 3. sepsis composite recalculé sur le lactate brut (la v3 l'évaluait
  #        sur un lactate imputé à la médiane, < 2 -> sans effet ; identique ici)
  m[, sepsis_composite_v3 := as.integer(vasopressor_pre_t0 == 1 | vent_pre_t0 == 1 |
                                           (!is.na(lactate_max) & lactate_max >= 2))]

  # --- 4. site infectieux (cultures prélevées), service regroupé, proxy contraste
  m <- merge(m, build_infection_site_v31(tbl_paths, key), by = c("subject_id", "hadm_id"), all.x = TRUE)
  m[, service_group := group_service_v31(service_at_t0)]
  m <- merge(m, build_contrast_proxy_v31(tbl_paths, key), by = c("subject_id", "hadm_id"), all.x = TRUE)
  m <- merge(m, compute_creat_counts_pre_t0_v31(creat_dt, key), by = c("subject_id", "hadm_id"), all.x = TRUE)

  # --- 5. eGFR baseline (CKD-EPI 2021) depuis la baseline min-7 j — modificateur
  #        de la BLP uniquement (A10), pas dans le PS (colinéaire avec créat/âge/sexe)
  m[, egfr_baseline_v31 := compute_egfr_ckd_epi(baseline_creat_min7d, age,
                                                 ifelse(toupper(gender) == "F", "F", "M"))]

  # --- 6. RRT 14 j (pour le stade KDIGO à J14 et l'outcome ordinal ACORN)
  rrt14 <- extract_rrt_window_v3(tbl_paths, m, window_days = 14)
  m <- merge(m, rrt14, by = c("subject_id", "hadm_id"), all.x = TRUE)

  # --- 7. outcomes KDIGO stricts (A6)
  kd7  <- compute_kdigo_v31(creat_dt, m, horizon_d = 7,  rrt_flag_col = "rrt_7d")
  kd14 <- compute_kdigo_v31(creat_dt, m, horizon_d = 14, rrt_flag_col = "rrt_14d")
  m <- merge(m, kd7, by = c("subject_id", "hadm_id"), all.x = TRUE)
  m <- merge(m, kd14, by = c("subject_id", "hadm_id"), all.x = TRUE)

  # persistance (A8) : NA explicite si non évaluable
  p48 <- compute_aki_persistence_v31(creat_dt, m, "aki_d7_v31", "t_first_aki_h_d7_v31", 48)
  p72 <- compute_aki_persistence_v31(creat_dt, m, "aki_d7_v31", "t_first_aki_h_d7_v31", 72)
  setnames(p48, "out", "aki_persistent_48h_v31"); setnames(p72, "out", "aki_persistent_72h_v31")
  m <- merge(m, p48, by = c("subject_id", "hadm_id"), all.x = TRUE)
  m <- merge(m, p72, by = c("subject_id", "hadm_id"), all.x = TRUE)

  # sensibilité « une valeur max par 24 h » sous la règle stricte (biais de surveillance)
  creat_dm <- merge(creat_dt, key[, .(subject_id, hadm_id, T0)], by = c("subject_id", "hadm_id"))
  creat_dm[, day := floor(as.numeric(difftime(charttime, T0, units = "hours")) / 24)]
  creat_dm <- creat_dm[order(-creat_mgdl)][, .SD[1], by = .(subject_id, hadm_id, day)]
  kdm <- compute_kdigo_v31(creat_dm[, .(subject_id, hadm_id, charttime, creat_mgdl)], m,
                           horizon_d = 7, rrt_flag_col = "rrt_7d", suffix = "_d7_dailymax_v31")
  m <- merge(m, kdm[, .(subject_id, hadm_id, aki_d7_dailymax_v31)], by = c("subject_id", "hadm_id"), all.x = TRUE)

  # --- 8. contrôles négatifs (A7)
  m <- merge(m, compute_early_creat_rise_v31(creat_dt, m), by = c("subject_id", "hadm_id"), all.x = TRUE)
  m[, nc_hypoglycemia_d7 := fifelse(!is.na(n_glucose_d7) & n_glucose_d7 > 0 &
                                      !(hypoglycemia_pre_t0 %in% 1L), hypoglycemia_d7, NA_integer_)]

  # --- 9. arrivée (urgences ou admission) pour la phase B
  adm <- read_mimic_table(tbl_paths$paths[["admissions"]], col_select = c("hadm_id", "edregtime"))
  adm[, `:=`(hadm_id = as.integer(hadm_id), edregtime = lubridate::ymd_hms(edregtime, quiet = TRUE))]
  m <- merge(m, adm, by = "hadm_id", all.x = TRUE)
  m[, arrival_time := pmin(admittime, fcoalesce(edregtime, admittime))]
  m[, h_arrival_to_T0 := as.numeric(difftime(T0, arrival_time, units = "hours"))]
  # Délai jusqu'à T0 : add_timing_vars() (v2/v3) mettait à NA puis à la MÉDIANE
  # (41,5 h) tout délai admission -> T0 négatif, c'est-à-dire l'antibiotique
  # donné aux urgences AVANT l'admission formelle (16,7 % des PTZ vs 6,7 % des
  # carbapénèmes) : les patients traités le plus tôt recevaient un délai de
  # 41,5 h. v3.1 : délai depuis l'ARRIVÉE (urgences ou admission), borné à 0,
  # sans imputation. L'ancienne valeur est gardée pour le waterfall.
  m[, time_admission_to_T0_h_v3 := time_admission_to_T0_h]
  m[, time_admission_to_T0_h := pmax(0, h_arrival_to_T0)]

  # --- 10. outcome primaire et cohortes
  m[, outcome_primary_v31 := aki_d7_v31]
  setcolorder(m, c("subject_id", "hadm_id", "treatment", "T0"))
  setkey(m, hadm_id)
  message("[50_master_v31] master : n=", nrow(m), " (PTZ=", sum(m$treatment == "PTZ"),
          ", carbapénème=", sum(m$treatment == "carbapenem"), ") ; outcome primaire évaluable : ",
          sum(!is.na(m$outcome_primary_v31)))
  m
}

#' Distribution de l'AKI à T0 par bras AVANT exclusion (recalculée sur le
#' cache v2 avec la définition v3/v3.1 : Cr à T0 >= 1,5 x min 7 j, ou +0,3
#' vs le min des 48 h précédentes) — alimente le flowchart et §8.2.
compute_aki_at_t0_table_v31 <- function(cohort_v2, creat_dt, out_tables = "output/extraction/tables") {
  b <- compute_baseline_creatinine_v3(creat_dt, cohort_v2[, .(subject_id, hadm_id, T0,
                                                               admit_time = admittime, age, gender)])
  x <- merge(cohort_v2[, .(subject_id, hadm_id, treatment)], b, by = c("subject_id", "hadm_id"))
  tbl <- x[, .(n = .N, n_evaluable = sum(aki_at_t0_evaluable), n_aki_at_t0 = sum(aki_at_t0),
               pct_aki_at_t0_among_evaluable = round(100 * sum(aki_at_t0) / sum(aki_at_t0_evaluable), 1)),
           by = treatment]
  fwrite(tbl, file.path(out_tables, "aki_at_t0_by_arm_v31.csv"))
  tbl
}
