# Piperacillin–tazobactam vs meropenem and acute kidney injury — target trial emulation in MIMIC-IV

R code for the analysis reported in the manuscript *"Piperacillin–Tazobactam vs Meropenem and Acute Kidney Injury in MIMIC-IV: a Target Trial Emulation"*.

> **No data in this repository.** MIMIC-IV is distributed by PhysioNet under a credentialed data use agreement that forbids redistribution. This repository contains code only. All outputs (patient-level derived data, imputations, tables, figures) are regenerated locally and are excluded by `.gitignore`.

---

## 1. Study in one paragraph

New-user, active-comparator target trial emulation of adult hospital admissions in MIMIC-IV 3.1 (with MIMIC-IV-ED 2.2 for emergency department vital signs). Patients start piperacillin–tazobactam (PTZ) or meropenem, based on documented administrations (eMAR). Time zero (T0) is the first administration. Patients with acute kidney injury (AKI) at T0 are excluded.

- **Primary population:** early empirical initiators (T0 ≤ 24 h after arrival).
- **Primary outcome:** strict KDIGO creatinine-based AKI by day 7.
- **Primary estimand:** average treatment effect in the overlap population (ATO), estimated with overlap weights from an unpenalized logistic propensity score with splines.
- **Missing data:** multiple imputation (m = 20, stratified by care location), with a bootstrap (B = 200 per imputation) combined with Rubin's rules.

Secondary and sensitivity analyses include a broad population, a three-arm trial that adds cefepime, the ACORN benchmark, negative controls, cumulative incidence with death as a competing risk, and heterogeneity of treatment effect with causal forests.

## 2. Requirements

- **R ≥ 4.4.** The published results were produced with R 4.4.2 on macOS (aarch64).
- **Data access:** credentialed PhysioNet access to
  - [MIMIC-IV 3.1](https://physionet.org/content/mimiciv/3.1/) (`hosp/`, `icu/`)
  - [MIMIC-IV-ED 2.2](https://physionet.org/content/mimic-iv-ed/2.2/) (`ed/`)
- **R packages:**

```r
install.packages(c(
  "data.table", "dplyr", "lubridate", "yaml", "readr", "arrow", "duckdb", "DBI",
  "glmnet", "mice", "cobalt", "PSweight", "survival", "splines", "nnet", "MASS",
  "boot", "sandwich", "lmtest", "survey", "EValue", "pROC", "diptest",
  "grf", "xgboost", "ranger", "tmle", "MatchIt",
  "dagitty", "ggdag", "ggraph", "ggplot2", "ggtern", "patchwork", "scales",
  "flextable", "officer", "jsonlite", "knitr", "rmarkdown", "testthat"
))
```

## 3. Expected directory layout

The code uses paths relative to the repository root. MIMIC-IV sits **next to** the repository:

```
parent_folder/
├── mimic-iv-3.1/          # hosp/, icu/   (config/config.yml: paths$mimic_root)
├── mimic-iv-ed-2.2/ed/    # ED tables     (R/72_ed_vitals.R: ED_ROOT_V32)
└── mimic_pipe_tazo_penem_emulated_target_trial/   # this repository (working directory)
```

If your data live elsewhere, edit `paths$mimic_root` in `config/config.yml` and `ED_ROOT_V32` in `R/72_ed_vitals.R`.

## 4. Running the pipeline

Run the numbered scripts in `scripts/` **from the repository root** and **in order**. Each step caches its results, so a later step reads the outputs of the earlier ones.

| Step | Script | What it does | Writes to |
|------|--------|--------------|-----------|
| 0 | `scripts/00_extract_base.R` | Extracts all valid creatinine values and all medications (eMAR > pharmacy > prescriptions) from MIMIC-IV | `output/extraction/` |
| 1 | `scripts/01_extract_labs_vitals.R` | Sets per-admission extraction windows, then extracts long-format labs and ICU chartevents before T0 | `output/extraction/` |
| 2 | `scripts/02_build_initiator_cohort.R` | Builds the unified cohort of initiators (PTZ, cefepime, meropenem) from documented administrations, with covariates and outcomes | `output/extraction/` |
| 3 | `scripts/03_master_cohort_imputation_primary.R` | Master cohort (one definition for all populations), creatinine trajectory, ED vital signs, stratified multiple imputation (m = 20), propensity-score rule, primary and three-arm estimates | `output/analysis/` |
| 4 | `scripts/04_secondary_sensitivity_figures.R` | DAG, positivity, flow chart, cumulative incidence, secondary outcomes, sensitivity analyses, negative controls, ACORN benchmark, Table 1, forest plots | `output/analysis/` |
| 5 | `scripts/05_heterogeneity.R` | Heterogeneity of treatment effect: cross-fitted causal forests, AIPW scores, tertiles, temporal validation | `output/analysis/hte/` |
| 6 | `scripts/06_protocol_deviations_checklist.R` | Protocol deviations, TARGET checklist, HTE rows added to `results_master` | `output/analysis/` |
| 7 | `scripts/07_robustness_primary.R` | Robustness of the primary estimate (PS complexity, parsimonious PS, augmented overlap weighting, cross-fitted AIPW), B = 200 bootstrap per imputation, rule of 20 events for the RR | `output/final/` |
| 8 | `scripts/08_exposure_flow_cumulative_incidence.R` | Overlap population, exposure (doses) and microbiology, STROBE flow by arm, KDIGO trigger criterion, cumulative incidence with CI bands | `output/final/` |
| 9 | `scripts/09_dag_negative_controls_trajectories.R` | Final DAG (edges checked against `dagitty`), 6 h and 12 h negative controls, creatinine and urea trajectories, figure legends | `output/final/` |
| 10 | `scripts/10_final_figures_legends.R` | Final forest plots (main results, sensitivity, secondary outcomes, negative controls) from `results_master`, summary tables, figure legends | `output/final/` |

**Final results** are in `output/final/`:
- `results/results_master.csv`: every estimate, one row per analysis;
- `tables/` and `figures/`;
- `figure_legends/legends.md`.

Steps 0–3 build the analysis data. From step 7 on, the master cohort and the 20 imputations are only read, never rebuilt, so every later analysis uses exactly the same data.

For quick debugging, set the environment variable `PIPELINE_FAST=1` to use few imputations and bootstrap replicates. Step 3 then writes to `output/analysis_fast/`.

Optional unit test of the KDIGO definitions: `Rscript -e 'testthat::test_file("tests/test_kdigo.R")'`.

**Run time.** Steps 3, 7 and 8 bootstrap within each of the 20 imputations and take several hours. The scripts use `parallel` with `options(v31.cores = 8L)`.

## 5. Architecture

### 5.1 Data flow

```
MIMIC-IV 3.1 + MIMIC-IV-ED 2.2
        │  steps 0–1  extraction (R/01–05, R/53)
        ▼
creatinine, medications, labs, chartevents                         output/extraction/
        │  step 2     administered antibiotics → initiator cohort (R/57, R/31, R/32, R/40, R/55)
        ▼
initiator cohort (PTZ / cefepime / meropenem, covariates, outcomes)
        │  step 3     master cohort + population registry + stratified MI (R/70–73, R/50–54)
        ▼
master cohort + 20 imputations                                     output/analysis/derived_data/
        │  steps 3–6  estimation, secondary and sensitivity analyses, HTE (R/52, R/58, R/61–62, R/74–78)
        │  steps 7–9  robustness, descriptive analyses, DAG (R/80–82)
        │  step 10    final forest plots and figure legends (R/77, R/83)
        ▼
results_master + tables + figures + legends                        output/final/
```

### 5.2 Design principles

- **One source of truth.** Every effect estimate is one row of `results_master` (`analysis_id`, cohort, estimand, role PRIMARY / SECONDARY / SENSITIVITY / …, estimate, CI, CI method). Figures and legends read the numbers from this table, the cohort registry and the tables; no result is typed by hand.
- **One entry point for cohorts.** Every analysis population is a registered subset of the master cohort (`R/51`, `R/71`). The registry (`n_by_analysis_v32.csv`) checks that each cohort is a subset of its parent and that strata partition their parent exactly.
- **Single PS/weighting engine.** `R/20_bootstrap.R` and `R/52_estimation.R` hold the propensity-score fit, overlap weights and bootstrap. The PS is refitted in every bootstrap replicate. Results are combined across imputations with Rubin's rules.
- **Single text formatter.** `R/83_text_formatter.R` writes every result sentence used in the figure legends, so the same number is always formatted the same way.
- **Built-in checks.** Each step writes QA tables to `tables/qa_*.csv`: cohort consistency, PS rule, DAG edges, etc.

### 5.3 Module map (`R/`)

The scripts in `scripts/` only orchestrate the work. All functions live in `R/`, and the numeric prefix gives the layer.

| Layer | Modules | Role |
|-------|---------|------|
| I/O & codebooks | `01_paths_and_io`, `02_codebooks` | Config loading, MIMIC table detection, ICD / itemid code lists |
| Extraction | `03_extract_medications`, `03b_extract_medications_extended`, `04_extract_labs`, `05_extract_icu_covariates`, `53_extract`, `72_ed_vitals` | Medications, creatinine and labs, ICU chartevents, ED vital signs |
| Cohort & covariates | `31_baseline_creat`, `32_covariates`, `57_antipseudo_cohort`, `50_master_cohort`, `70_master` | Baseline creatinine and AKI at T0, pre-T0 covariates (no information leakage), initiator cohort, master cohort |
| Populations | `51_analysis_cohort`, `71_registry` | Population registry, PS covariate set |
| Outcomes | `40_outcomes`, `55_outcomes` | Strict KDIGO AKI, stages, persistence, MAKE14, ACORN ordinal outcome |
| Missing data | `54_mi`, `73_mi_stratified` | Multiple imputation (m = 20), stratified by care location |
| Estimation | `10_estimators`, `15_positivity_diagnostics`, `20_bootstrap`, `52_estimation`, `58_designs`, `43_sensitivity`, `80_robustness` | PS, overlap / generalized overlap weights, bootstrap with Rubin's rules, three-arm design, E-value and bias analysis, robustness analyses |
| Heterogeneity | `61_hte`, `62_hte_figures` | Cross-fitted causal forests, AIPW scores, comparator learners, calibration |
| Results & figures | `30_results_master`, `56_results_tables`, `74_dag`, `82_dag_figure`, `75_positivity_flow`, `76_cif`, `77_tables_figs`, `81_exposure_flow_cif`, `78_reporting`, `83_text_formatter` | `results_master`, Table 1, DAG, positivity, flow chart, cumulative incidence, forest plots, trajectories, protocol deviations, text formatter |

Some function names, variable names and output file names carry an internal suffix (`_v31`, `_v32`, `_v33`, for example `outcome_primary_v31` or `n_by_analysis_v32.csv`). These are fixed identifiers kept so that the code produces exactly the published numbers. They do **not** refer to separate analyses.

### 5.4 Other files

- `config/config.yml`: MIMIC path, drug regex patterns, KDIGO thresholds.
- `tests/test_kdigo.R`: unit tests of the KDIGO AKI algorithm.

## 6. Known limitation of stand-alone reproduction

One supplementary table, the **reconciliation of two-arm and three-arm estimates** (rows `r1`, `r1b` and `r4`), compares the final estimate with a preliminary analysis of the same question that is not part of this repository. For this block (P1.6), step 3 expects these files:

- `output/preliminary/results_master.rds`
- `output/preliminary/master_v31.rds`
- `output/preliminary/mi_eligible_m20.rds`

Without them, block P1.6 fails. Comment it out to reproduce every other result. The primary, secondary and sensitivity estimates do not depend on these files.

## 7. Reproducibility notes

- Seed: `SEED_V3` in `R/20_bootstrap.R`, used by the bootstrap and imputation functions.
- The complete `sessionInfo()` of the published run is written to `output/final/results/session_info_v32.txt`.

## 8. Citation & data

If you use this code, please cite the article (reference to be added on publication) together with MIMIC-IV, MIMIC-IV-ED and PhysioNet:

- Johnson A, et al. MIMIC-IV, a freely accessible electronic health record dataset. *Sci Data* 2023;10:1. doi:10.1038/s41597-022-01899-x
- Goldberger A, et al. PhysioBank, PhysioToolkit, and PhysioNet. *Circulation* 2000;101:e215–e220.

## 9. License

The code is released under the [MIT License](LICENSE). The license covers the code only, not MIMIC-IV data, which remain subject to the PhysioNet Credentialed Health Data Use Agreement.
