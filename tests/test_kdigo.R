# Tests unitaires de la définition KDIGO stricte v3.1 (docs/v3.1_plan.md A6).
# Trajectoires de créatinine construites à la main, résultat attendu fixé a priori.
# Lancement : testthat::test_file("tests/test_kdigo.R") depuis la racine du projet.
library(testthat); library(data.table)
source(file.path(Sys.getenv("PROJECT_ROOT", "."), "R/55_outcomes.R"))

T0 <- as.POSIXct("2150-01-10 12:00:00", tz = "UTC")
mk <- function(id, hours, values, baseline) {
  list(creat = data.table(subject_id = id, hadm_id = id, charttime = T0 + hours * 3600, creat_mgdl = values),
       cohort = data.table(subject_id = id, hadm_id = id, T0 = T0, baseline_creat_min7d = baseline))
}
cases <- list(
  # 1. stable : aucune AKI
  stable       = mk(1L, c(-20, 10, 40, 100), c(0.9, 1.0, 0.9, 1.0), 0.9),
  # 2. +0,3 en 36 h avec une valeur PRÉ-T0 dans la fenêtre glissante : AKI (absolu)
  abs_pre_t0   = mk(2L, c(-12, 24), c(1.0, 1.3), 1.0),
  # 3. +0,3 étalé sur 60 h, aucune paire <= 48 h (paliers de +0,15) : pas d'AKI
  slow_60h     = mk(3L, c(-12, 18, 48), c(1.00, 1.15, 1.30), 1.0),
  # 4. montée lente sur 5 j jusqu'à 1,5 x la baseline min 7 j : AKI (relatif seul)
  rel_slow     = mk(4L, c(-2, 24, 48, 72, 96, 120), c(1.0, 1.1, 1.2, 1.3, 1.4, 1.5), 1.0),
  # 5. bug v3 : unique dosage post-T0 = baseline 7 j + 0,3, aucune valeur dans
  #    les 48 h précédentes -> v3 = AKI (repli baseline), v3.1 = pas d'AKI
  v3_fallback  = mk(5L, c(-100, 30), c(1.0, 1.3), 1.0),
  # 6. stade 3 par Cr >= 4,0 (critère absolu rempli, < 3 x baseline)
  stage3_cr4   = mk(6L, c(-6, 20), c(3.6, 4.0), 2.0),
  # 7. aucune créatinine post-T0 : non évaluable
  no_post      = mk(7L, c(-6), c(1.0), 1.0)
)
creat  <- rbindlist(lapply(cases, `[[`, "creat"))
cohort <- rbindlist(lapply(cases, `[[`, "cohort"))
res <- compute_kdigo_v31(creat, cohort, horizon_d = 7)
setkey(res, subject_id)

test_that("1. trajectoire stable : pas d'AKI", {
  expect_equal(res[J(1L)]$aki_d7_v31, 0L)
})
test_that("2. +0,3 en 36 h avec référence pré-T0 : AKI par critère absolu", {
  expect_equal(res[J(2L)]$aki_d7_v31, 1L)
  expect_equal(res[J(2L)]$aki_trigger_d7_v31, "abs_0.3_48h_only")
})
test_that("3. +0,3 en 60 h sans paire <= 48 h : pas d'AKI", {
  expect_equal(res[J(3L)]$aki_d7_v31, 0L)
})
test_that("4. montée lente à 1,5 x baseline : AKI par critère relatif seul", {
  expect_equal(res[J(4L)]$aki_d7_v31, 1L)
  expect_equal(res[J(4L)]$aki_trigger_d7_v31, "rel_1.5x_only")
  expect_equal(res[J(4L)]$aki_stage_max_d7_v31, 1L)
})
test_that("5. pas de repli sur la baseline 7 j pour le critère absolu (bug v3)", {
  expect_equal(res[J(5L)]$aki_d7_v31, 0L)
})
test_that("6. Cr >= 4,0 avec critère absolu : stade 3", {
  expect_equal(res[J(6L)]$aki_stage_max_d7_v31, 3L)
})
test_that("7. aucune créatinine post-T0 : non évaluable (NA)", {
  expect_true(is.na(res[J(7L)]$aki_d7_v31))
})
