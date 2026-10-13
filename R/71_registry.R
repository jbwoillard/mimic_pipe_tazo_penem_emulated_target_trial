# =============================================================================
# 71_registry.R — Covariables du PS v3.2 et registre des populations
#   (docs/v3.2_plan.md §1.1, §1.2)
# =============================================================================

TRAJ_COVS_V32 <- c("creat_last_pre_t0", "creat_delta_last_min7d", "creat_ratio_last_min7d",
                   "creat_slope_48h", "creat_slope_48h_available", "creat_single_pre_value")
PS_COVARIATES_V32 <- c(PS_COVARIATES_V31, TRAJ_COVS_V32)
# variante « constantes de réanimation seules » (étapes w1-w2 du waterfall)
ICU_VITALS_V32 <- c(hr_mean = "hr_mean_icu", map_mean = "map_mean_icu", temp_max = "temp_max_icu", spo2_min = "spo2_min_icu")
icu_vitals_covs <- function(covs) { i <- covs %in% names(ICU_VITALS_V32); covs[i] <- ICU_VITALS_V32[covs[i]]; covs }

# splines aussi pour les nouvelles continues
CONTINUOUS_PS_VARS_V31 <- unique(c(CONTINUOUS_PS_VARS_V31, unname(ICU_VITALS_V32),
                                   "creat_last_pre_t0", "creat_delta_last_min7d", "creat_ratio_last_min7d", "creat_slope_48h"))
VARS_WITH_MISSING_V32 <- c(VARS_WITH_MISSING_V31, unname(ICU_VITALS_V32))

# comparateur des analyses à 2 bras : méropénème (le seul carbapénème de la cohorte unifiée)
EARLY_FILTER_V32 <- quote(h_arrival_to_T0 <= 24 & true_new_user_90d == 1L & any_iv_abx_48h_pre == 0L)
ACORN_FILTER_V32 <- quote(treatment %in% c("PTZ", "cefepime") & h_arrival_to_T0 <= 12 &
                            service_group %in% c("ED_observation", "ICU_medical") & culture_drawn_pm24h == 1L &
                            acorn_new_user_7d == 1L)

#' Colonnes dérivées + registre ; renvoie list(master, SUBGROUPS)
setup_registry_v32 <- function(master) {
  m <- copy(master)
  m[, outcome_moss_s3b := fifelse(!is.na(outcome_primary), outcome_primary, fifelse(treatment == "PTZ", 0L, 1L))]
  m[, outcome_moss_s4b := fifelse(!is.na(outcome_primary), outcome_primary, fifelse(treatment == "PTZ", 1L, 0L))]
  for (v in VARS_WITH_MISSING_V32) m[, (paste0(v, "_missing")) := as.integer(is.na(get(v)))]
  m[, egfr_group := fcase(is.na(egfr_baseline_v31), "missing", egfr_baseline_v31 < 60, "<60", egfr_baseline_v31 < 90, "60-89", default = ">=90")]
  setkey(m, hadm_id)
  C <- function(id, parent, filter = NULL, outcome = NULL, label = id)
    invisible(get_analysis_cohort(list(id = id, parent = parent, filter = filter, outcome = outcome, label = label), m))
  C("all", "master", label = "First documented administration of PTZ, cefepime or meropenem (one episode per patient)")
  C("eligible", "all", quote(aki_at_t0 == 0L), label = "No AKI at T0")
  C("eligible_2arm", "eligible", quote(treatment %in% c("PTZ", "meropenem")), label = "PTZ or meropenem initiators")
  C("broad", "eligible_2arm", outcome = "outcome_primary", label = "BROAD: PTZ vs meropenem, all initiations, AKI-D7 evaluable")
  C("early", "broad", EARLY_FILTER_V32,
    label = "EARLY (primary): T0 <=24 h after arrival, no anti-pseudomonal beta-lactam in 90 days, no IV antibiotic in 48 h")
  C("three_arm_early", "eligible", EARLY_FILTER_V32, outcome = "outcome_primary",
    label = "Three-arm (PTZ / cefepime / meropenem), early population")
  C("acorn_with_aki_t0", "all", ACORN_FILTER_V32, outcome = "acorn_ordinal_d14", label = "ACORN emulation including AKI at T0")
  C("acorn", "acorn_with_aki_t0", quote(aki_at_t0 == 0L), label = "ACORN emulation (pre-specified: no AKI at T0)")
  C("acorn_aki7", "acorn", outcome = "aki_d7_v31", label = "ACORN, AKI-D7 evaluable")
  C("acorn_make14", "acorn", outcome = "make14", label = "ACORN, MAKE14 evaluable")
  for (pop in c("broad", "early")) {
    C(paste0(pop, "_nc_a12"), pop, outcome = "early_creat_rise_12h", label = paste0(pop, ": negative control (a), 12 h"))
    C(paste0(pop, "_nc_a6"), pop, outcome = "early_creat_rise_6h", label = paste0(pop, ": negative control (a), 6 h"))
    C(paste0(pop, "_nc_b"), pop, outcome = "nc_hypoglycemia_d7", label = paste0(pop, ": negative control (b)"))
    C(paste0(pop, "_persist48"), pop, outcome = "aki_persistent_48h_v31", label = paste0(pop, ": persistent AKI 48 h evaluable"))
    C(paste0(pop, "_persist72"), pop, outcome = "aki_persistent_72h_v31", label = paste0(pop, ": persistent AKI 72 h evaluable"))
    C(paste0(pop, "_dailymax"), pop, outcome = "aki_d7_dailymax_v31", label = paste0(pop, ": daily-max outcome evaluable"))
    C(paste0(pop, "_ge2_creat"), pop, quote(n_creat_pre_7d >= 2L), label = paste0(pop, ": >=2 pre-T0 creatinine values"))
    C(paste0(pop, "_complete_vitals"), pop, quote(!is.na(hr_mean) & !is.na(map_mean) & !is.na(temp_max) & !is.na(spo2_min) & !is.na(gcs_min)),
      label = paste0(pop, ": vital signs measured"))
    C(paste0(pop, "_vanco_with"), pop, quote(vanco_6h == 1L), label = paste0(pop, ": vancomycin at T0"))
    C(paste0(pop, "_vanco_without"), pop, quote(vanco_6h == 0L), label = paste0(pop, ": no vancomycin at T0"))
    register_strata_family(paste0(pop, "_vancomycin"), pop, paste0(pop, c("_vanco_with", "_vanco_without")))
  }
  list(master = m)
}

#' QA B1 : une seule définition de la cohorte précoce pour 2 et 3 bras
qa_early_ids_v32 <- function(master) {
  e2 <- .COHORT_REGISTRY_V31$ids[["early"]]
  e3 <- .COHORT_REGISTRY_V31$ids[["three_arm_early"]]
  e3_no_cef <- e3[master[J(e3)]$treatment != "cefepime"]
  setequal(e2, e3_no_cef)
}
