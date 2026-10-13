# =============================================================================
# 03b_extract_medications_extended.R — Ré-extraction v3 des médicaments (nouvelles
#   classes) — docs/v3_plan.md §1.9 (vrai new-user), §2.5 (nouvelles
#   covariables), §3.7 (comparateur céfépime)
# =============================================================================
# Copie adaptée de R/03_extract_medications.R:extract_medications() — logique
# identique (emar > pharmacy > prescriptions, pré-filtrage par regex avant
# jointure), mais :
#   - classes supplémentaires : cefepime, ceftazidime, aztreonam (§1.9/§3.7),
#     acei_arb, nsaid (§2.5) ;
#   - écrit dans output/extraction/derived_data/ (JAMAIS dans output/extraction/derived_data/,
#     pour ne pas écraser le cache v2 — règle #2 du cahier des charges).
#
# N'est PAS appelé par run_all.R (v2). Appelé uniquement par run_all_v3.R.

extract_medications_v3 <- function(tbl_paths, cfg,
                                    out_path = "output/extraction/derived_data/medications_all_v3.rds") {

  paths      <- tbl_paths$paths
  med_source <- tbl_paths$med_source

  ptz_pat      <- paste(cfg$medications$ptz_patterns,          collapse = "|")
  carb_pat     <- paste(cfg$medications$carbapenem_patterns,   collapse = "|")
  vanco_pat    <- paste(cfg$medications$vanco_patterns,        collapse = "|")
  amino_pat    <- paste(cfg$medications$aminoglycoside_patterns, collapse = "|")
  nephr_pat    <- paste(cfg$medications$other_nephrotoxins,    collapse = "|")
  cefepime_pat <- paste(cfg$medications$cefepime_patterns,     collapse = "|")
  ceftaz_pat   <- paste(cfg$medications$ceftazidime_patterns,  collapse = "|")
  aztreo_pat   <- paste(cfg$medications$aztreonam_patterns,    collapse = "|")
  acei_pat     <- paste(cfg$medications$acei_arb_patterns,     collapse = "|")
  nsaid_pat    <- paste(cfg$medications$nsaid_patterns,        collapse = "|")

  classify_drug_v3 <- function(drug_name) {
    dn <- tolower(drug_name)
    dplyr::case_when(
      grepl(ptz_pat,      dn, ignore.case = TRUE) & !grepl("^tazobactam$", dn) ~ "ptz",
      grepl(carb_pat,     dn, ignore.case = TRUE) ~ "carbapenem",
      grepl(cefepime_pat, dn, ignore.case = TRUE) ~ "cefepime",
      grepl(ceftaz_pat,   dn, ignore.case = TRUE) ~ "ceftazidime",
      grepl(aztreo_pat,   dn, ignore.case = TRUE) ~ "aztreonam",
      grepl(vanco_pat,    dn, ignore.case = TRUE) ~ "vancomycin",
      grepl(amino_pat,    dn, ignore.case = TRUE) ~ "aminoglycoside",
      grepl(nephr_pat,    dn, ignore.case = TRUE) ~ "other_nephrotoxin",
      grepl(acei_pat,     dn, ignore.case = TRUE) ~ "acei_arb",
      grepl(nsaid_pat,    dn, ignore.case = TRUE) ~ "nsaid",
      TRUE                                        ~ NA_character_
    )
  }

  all_pat <- paste(unique(c(
    cfg$medications$ptz_patterns, cfg$medications$carbapenem_patterns,
    cfg$medications$vanco_patterns, cfg$medications$aminoglycoside_patterns,
    cfg$medications$other_nephrotoxins, cfg$medications$cefepime_patterns,
    cfg$medications$ceftazidime_patterns, cfg$medications$aztreonam_patterns,
    cfg$medications$acei_arb_patterns, cfg$medications$nsaid_patterns
  )), collapse = "|")

  meds <- NULL

  if (med_source %in% c("emar", "emar_detail") &&
      !is.na(paths[["emar"]]) && !is.na(paths[["emar_detail"]])) {

    message("[03b_meds_v3] Utilisation de emar + emar_detail (administrations réelles)")

    emar <- read_mimic_table(paths[["emar"]],
                              col_select = c("subject_id", "hadm_id", "emar_id",
                                             "medication", "scheduletime",
                                             "charttime", "event_txt"))
    ed   <- read_mimic_table(paths[["emar_detail"]],
                              col_select = c("subject_id", "emar_id",
                                             "product_description"))

    if (!is.null(emar) && !is.null(ed)) {
      emar_filt <- emar[grepl(all_pat, medication, ignore.case = TRUE)]
      message("[03b_meds_v3] emar après pre-filtrage par nom : ",
              nrow(emar_filt), " lignes (sur ", nrow(emar), ")")

      ed_filt <- ed[emar_id %in% emar_filt$emar_id, .(emar_id, product_description)]

      emar_filt[, drug_name := medication]
      emar_joined <- merge(emar_filt, ed_filt, by = "emar_id", all.x = TRUE)
      emar_joined[!is.na(product_description) & product_description != "",
                  drug_name := product_description]

      rm(emar, ed, ed_filt); gc()

      emar_joined[, starttime := dplyr::coalesce(
        lubridate::ymd_hms(charttime, quiet = TRUE),
        lubridate::ymd_hms(scheduletime, quiet = TRUE)
      )]

      emar_joined[, drug_class := classify_drug_v3(drug_name)]
      meds <- emar_joined[!is.na(drug_class),
                           .(subject_id, hadm_id, starttime, drug_name, drug_class,
                             source = "emar")]
    }
  }

  if (is.null(meds)) stop("[03b_meds_v3] Extraction eMAR v3 a échoué (voir logs ci-dessus).")

  meds <- meds[!is.na(starttime)]
  meds[, subject_id := as.integer(subject_id)]
  meds[, hadm_id    := as.integer(hadm_id)]
  setorder(meds, subject_id, hadm_id, starttime)

  meds[drug_class == "carbapenem",
       carbapenem_molecule := dplyr::case_when(
         grepl("meropenem",   drug_name, ignore.case = TRUE) ~ "meropenem",
         grepl("imipenem",    drug_name, ignore.case = TRUE) ~ "imipenem",
         grepl("doripenem",   drug_name, ignore.case = TRUE) ~ "doripenem",
         grepl("ertapenem",   drug_name, ignore.case = TRUE) ~ "ertapenem",
         TRUE                                                ~ "carbapenem_other"
       )]

  # Décomposition de other_nephrotoxin en sous-classes cliniques (docs/v3_plan.md
  # §2.5) — n'exige pas de ré-extraction supplémentaire, seulement une
  # re-classification du texte drug_name déjà extrait.
  meds[drug_class == "other_nephrotoxin",
       nephrotoxin_subclass := dplyr::case_when(
         grepl("amphotericin|colistin|polymyxin|acyclovir", drug_name, ignore.case = TRUE) ~ "direct_nephrotoxin",
         grepl("tacrolimus|cyclosporine",                    drug_name, ignore.case = TRUE) ~ "calcineurin_inhibitor",
         grepl("furosemide|bumetanide|torsemide",             drug_name, ignore.case = TRUE) ~ "loop_diuretic",
         TRUE                                                                                 ~ "other"
       )]

  message("[03b_meds_v3] Administrations extraites (v3) : ", nrow(meds))
  message("[03b_meds_v3] Répartition par classe :")
  print(meds[, .N, by = drug_class])
  message("[03b_meds_v3] Sous-classes other_nephrotoxin :")
  print(meds[drug_class == "other_nephrotoxin", .N, by = nephrotoxin_subclass])

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(meds, out_path)
  message("[03b_meds_v3] Sauvegardé (lecture seule, n'écrase pas le cache v2) : ", out_path)

  meds
}
