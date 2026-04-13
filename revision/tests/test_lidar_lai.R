# ---
# title:  test_lidar_lai.R
# desc:   Unit tests for compute_pai() and compute_lad_profile() in
#         revision/R/lidar_lai.R. Tests verify numerical equivalence with
#         the original functions (myPAI / myPAD / myPAD_dtm) from
#         02_CODES/libraries/functions_lidar.R, and correct behaviour of
#         the new h_min / k arguments.
#
# Run from the project root:
#   source("revision/R/lidar_lai.R")
#   source("02_CODES/libraries/functions_lidar.R")   # for reference values
#   source("revision/tests/test_lidar_lai.R")
# ---

library("lidR")

# ── helpers ────────────────────────────────────────────────────────────────────

ok  <- function(label) cat("  [PASS]", label, "\n")
fail <- function(label, msg) stop(paste("[FAIL]", label, "->", msg))

assert_near <- function(a, b, tol = 1e-9, label) {
  if (is.na(a) || is.na(b)) fail(label, paste("NA: a =", a, "b =", b))
  if (abs(a - b) > tol) fail(label, paste("a =", a, "!= b =", b))
  ok(label)
}

assert_true <- function(expr, label) {
  if (!isTRUE(expr)) fail(label, deparse(substitute(expr)))
  ok(label)
}

# ── shared test data ───────────────────────────────────────────────────────────

set.seed(42)
z_mixed <- c(
  runif(80,  0,  1.5),   # sub-canopy / gap returns (below 2 m)
  runif(120, 2,  15),    # mid-canopy
  runif(100, 15, 35)     # upper canopy
)

z_dense <- runif(300, 5, 35)   # no points below 2 m → PAI = Inf

z_empty <- numeric(0)

# ── Tests: compute_pai ─────────────────────────────────────────────────────────

cat("\n=== compute_pai() ===\n")

# 1. Numerical equivalence with myPAI (h_min = 2, k = 0.5)
pai_new <- compute_pai(z_mixed, h_min = 2, k = 0.5)
pai_old <- myPAI(z_mixed, zmin = 2, k = 0.5)
assert_near(pai_new, pai_old, tol = 1e-12,
            "compute_pai(h_min=2, k=0.5) == myPAI(zmin=2, k=0.5)")

# 2. k sensitivity: result changes with k
pai_k04 <- compute_pai(z_mixed, h_min = 2, k = 0.4)
pai_k06 <- compute_pai(z_mixed, h_min = 2, k = 0.6)
assert_true(pai_k04 > pai_new, "PAI increases as k decreases (k=0.4 > k=0.5)")
assert_true(pai_k06 < pai_new, "PAI decreases as k increases (k=0.6 < k=0.5)")

# 3. h_min sensitivity: result changes with h_min
pai_h3 <- compute_pai(z_mixed, h_min = 3, k = 0.5)
pai_h5 <- compute_pai(z_mixed, h_min = 5, k = 0.5)
assert_true(pai_h3 < pai_new, "PAI decreases as h_min increases (more gap)")
assert_true(pai_h5 < pai_h3,  "PAI decreases further at h_min=5")

# 4. Edge cases
assert_true(is.na(compute_pai(z_empty)),           "empty z returns NA")
assert_true(is.infinite(compute_pai(z_dense)),      "no gap returns Inf")
assert_true(is.na(compute_pai(c(NA, NA, NA))),     "all-NA z returns NA")

# 5. Single point below threshold
assert_true(compute_pai(1.0, h_min = 2) == 0,      "single point below h_min: PAI = 0")

cat("All compute_pai() tests passed.\n")

# ── Tests: compute_lad_profile (dsm) ──────────────────────────────────────────

cat("\n=== compute_lad_profile(norm='dsm') ===\n")

prof_new <- compute_lad_profile(z_mixed, norm = "dsm", z_max = 40, k = 0.5)
prof_old <- myPAD(z_mixed)

# 1. Output structure: names match between old and new
lad_names_new <- grep("^LAD_Layer_", names(prof_new), value = TRUE)
lad_names_old <- grep("^LAD_Layer_", names(prof_old), value = TRUE)
assert_true(
  setequal(lad_names_new, lad_names_old),
  "LAD_Layer_* names identical to myPAD output"
)

pai_names_new <- grep("^PAI_", names(prof_new), value = TRUE)
pai_names_old <- grep("^PAI_", names(prof_old), value = TRUE)
assert_true(
  setequal(pai_names_new, pai_names_old),
  "PAI_* names identical to myPAD output"
)

# 2. LAD values: numerical equivalence for non-NA layers
common_lad <- intersect(lad_names_new, lad_names_old)
mismatches <- sapply(common_lad, function(nm) {
  v_new <- as.numeric(prof_new[[nm]])
  v_old <- as.numeric(prof_old[[nm]])
  if (is.na(v_new) && is.na(v_old)) return(FALSE)   # both NA: ok
  if (is.na(v_new) || is.na(v_old)) return(TRUE)     # one NA: mismatch
  abs(v_new - v_old) > 1e-9
})
n_mismatch <- sum(mismatches)
assert_true(
  n_mismatch == 0,
  paste0("All LAD_Layer values match myPAD (", length(common_lad),
         " bins checked)")
)

# 3. Grid always starts at 2.5 (z0 = 2 fixed)
first_lad <- names(prof_new)[grepl("^LAD_Layer_", names(prof_new))][1]
assert_true(first_lad == "LAD_Layer_2.5", "First LAD bin always at 2.5 m")

# 4. k parameterisation: different k changes LAD values
prof_k04 <- compute_lad_profile(z_mixed, norm = "dsm", z_max = 40, k = 0.4)
lad_mid   <- "LAD_Layer_10.5"
assert_true(
  !isTRUE(all.equal(prof_new[[lad_mid]], prof_k04[[lad_mid]])),
  "k=0.4 gives different LAD_Layer_10.5 than k=0.5"
)

cat("All compute_lad_profile(dsm) tests passed.\n")

# ── Tests: compute_lad_profile (dtm) ──────────────────────────────────────────

cat("\n=== compute_lad_profile(norm='dtm') ===\n")

prof_dtm_new <- compute_lad_profile(z_mixed, norm = "dtm", z_max = 40, k = 0.5)
prof_dtm_old <- myPAD_dtm(z_mixed)

lad_names_dtm_new <- grep("^LAD_Layer_", names(prof_dtm_new), value = TRUE)
lad_names_dtm_old <- grep("^LAD_Layer_", names(prof_dtm_old), value = TRUE)
assert_true(
  setequal(lad_names_dtm_new, lad_names_dtm_old),
  "DTM: LAD_Layer_* names identical to myPAD_dtm output"
)

common_dtm <- intersect(lad_names_dtm_new, lad_names_dtm_old)
mismatches_dtm <- sapply(common_dtm, function(nm) {
  v_new <- as.numeric(prof_dtm_new[[nm]])
  v_old <- as.numeric(prof_dtm_old[[nm]])
  if (is.na(v_new) && is.na(v_old)) return(FALSE)
  if (is.na(v_new) || is.na(v_old)) return(TRUE)
  abs(v_new - v_old) > 1e-9
})
n_dtm_mismatch <- sum(mismatches_dtm)
assert_true(
  n_dtm_mismatch == 0,
  paste0("All LAD_Layer values match myPAD_dtm (", length(common_dtm),
         " bins checked)")
)

cat("All compute_lad_profile(dtm) tests passed.\n")

# ── Tests: normalize_dsm vs end_dsm_norm ──────────────────────────────────────
# !! AVANT DE LANCER CE BLOC : remplace LAS_TILE_PATH par un chemin réel !!
# Exemple : "/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/
#             3-las_normalized_utm/XXXX_YYYY.las"
# Nathan fournira le chemin. Ne pas hardcoder une valeur fictive.

cat("\n=== normalize_dsm() vs end_dsm_norm() ===\n")

LAS_TILE_PATH <- "/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/3-las_normalized_utm/LAS_367440.46_5270188.33.las"

if (is.null(LAS_TILE_PATH)) {
  cat("  [SKIP] LAS_TILE_PATH non défini — section normalize_dsm non exécutée.\n")
} else {
  las_raw <- lidR::readLAS(LAS_TILE_PATH, select = "xyzc")

  # Appliquer les deux fonctions sur la même tuile
  nlas_old <- end_dsm_norm(las_raw)
  nlas_new <- normalize_dsm(las_raw)   # h_target=40, pitfree, veg_class=4

  # 1. Même nombre de points
  assert_true(
    lidR::npoints(nlas_old) == lidR::npoints(nlas_new),
    paste0("normalize_dsm: même nb points que end_dsm_norm (",
           lidR::npoints(nlas_old), ")")
  )

  # 2. Même Classification (ordre des points préservé)
  assert_true(
    identical(nlas_old$Classification, nlas_new$Classification),
    "normalize_dsm: Classification identique à end_dsm_norm"
  )

  # 3. Z identiques à 1e-6 près (floating-point DSM interpolation)
  z_diff <- abs(nlas_old$Z - nlas_new$Z)
  max_diff <- max(z_diff, na.rm = TRUE)
  assert_true(
    max_diff < 1e-6,
    paste0("normalize_dsm: Z identiques à 1e-6 m (max diff = ",
           format(max_diff, scientific = TRUE), ")")
  )

  cat("All normalize_dsm() tests passed.\n")
}

# ── Tests: normalize_dtm vs end_dtm_norm ──────────────────────────────────────

cat("\n=== normalize_dtm() vs end_dtm_norm() ===\n")

if (is.null(LAS_TILE_PATH)) {
  cat("  [SKIP] LAS_TILE_PATH non défini — section normalize_dtm non exécutée.\n")
} else {
  # Recharger la tuile brute (non normalisée) pour le test DTM
  # Si end_dtm_norm attend un LAS non normalisé, utiliser las_raw.
  # Si la tuile LAS disponible est déjà normalisée DTM, Nathan précisera.
  las_raw2 <- lidR::readLAS(LAS_TILE_PATH, select = "xyzc")

  nlas_dtm_old <- end_dtm_norm(las_raw2)
  nlas_dtm_new <- normalize_dtm(las_raw2)   # h_min=2, z_max=40

  # 1. Même nombre de points
  assert_true(
    lidR::npoints(nlas_dtm_old) == lidR::npoints(nlas_dtm_new),
    paste0("normalize_dtm: même nb points que end_dtm_norm (",
           lidR::npoints(nlas_dtm_old), ")")
  )

  # 2. Même Classification
  assert_true(
    identical(nlas_dtm_old$Classification, nlas_dtm_new$Classification),
    "normalize_dtm: Classification identique à end_dtm_norm"
  )

  # 3. Z identiques à 1e-6 près
  z_diff_dtm <- abs(nlas_dtm_old$Z - nlas_dtm_new$Z)
  max_diff_dtm <- max(z_diff_dtm, na.rm = TRUE)
  assert_true(
    max_diff_dtm < 1e-6,
    paste0("normalize_dtm: Z identiques à 1e-6 m (max diff = ",
           format(max_diff_dtm, scientific = TRUE), ")")
  )

  cat("All normalize_dtm() tests passed.\n")
}

cat("\n=== All tests passed ===\n\n")
