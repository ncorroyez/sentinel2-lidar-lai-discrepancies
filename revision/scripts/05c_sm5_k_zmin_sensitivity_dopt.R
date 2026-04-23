# ---
# title:  05c_sm5_k_zmin_sensitivity_dopt.R
# desc:   Sensitivity of Pareto d_opt to the extinction coefficient k
#         (0.3 to 0.8) and to the minimum height in the LAD equation z_min
#         (2, 3, 5 m).
#
#         k sensitivity: Beer-Lambert gives LAI_ALS proportional to 1/k,
#         so LAI_ALS_k = LAI_ALS_{k_ref} x (k_ref / k). Scaled from the
#         existing PAD sampling CSVs; S2 estimates unchanged.
#
#         z_min sensitivity: the PAD rasters integrate from canopy top down
#         to x.5 m above ground (PAD_{x.5}_40.tif, depth = 40 - floor(x)).
#         The existing z_min is ~2 m (deepest file: PAD_2.5_40.tif).
#         For z_min > 2 m, each pixel's PAD at depth d is capped at the
#         PAD value for depth min(d, floor(max_height - z_min)), where
#         max_height per pixel is read from the sampling GPKG produced by 03a.
#         This correctly excludes the bottom z_min metres of each canopy.
#
#         h_min_pixel is fixed at 10 m (reference analysis value). Both
#         sensitivities are computed for DSM_keepTrees and DTM_keepTrees,
#         for each of the three sites.
#
# Prerequisites:
#   03a — Sampling_stratified_uniform_hmin10_nbSamples_5000.GPKG
#   03b — PAD CSVs (all depths 1..38)
#   04a — ATBD LAI estimated CSVs (hmin10)
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/05c_sm5_k_zmin_sensitivity_dopt.R")
# ---

library("here")
library("data.table")
library("terra")
library("cli")

source(here::here("revision", "R", "sm5_metrics.R"))
source(here::here("revision", "R", "sm5_dopt.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites           <- c("Aigoual", "Blois", "Mormal")
norms_select    <- c("DSM_keepTrees", "DTM_keepTrees")
depths          <- 1:38
sampling_method <- "stratified_uniform"
name_strategy   <- "LIDFa_lai_LMA_BROWN"
nb_samples      <- 5000L
h_min_pixel     <- 10L          # pixel selection threshold (fixed)
k_ref           <- 0.5          # k used to compute the existing PAD CSVs
k_values        <- seq(0.3, 0.8, by = 0.1)
z_min_values    <- c(2, 3, 5)   # metres above ground

# ── Helper: build PAD wide matrix + max_height for one (site, norm) ───────────
# Returns a list:
#   $pad_mat  : n × 38 numeric matrix (row = sample, col = depth)
#   $max_h    : n-vector of max_height per sample (matched to row order)
#   $s2_vals  : n-vector of ATBD LAI_S2 estimates
#   $sample_ids: n-vector of sample IDs

load_site_norm_data <- function(site, norm) {

  gpkg_path <- here::here(
    "03_RESULTS", site, "PROSAIL_Optimization", "sampling",
    paste0("Sampling_", sampling_method,
           "_hmin", h_min_pixel,
           "_nbSamples_", nb_samples, ".GPKG")
  )
  s2_csv <- here::here(
    "revision", "output", "intermediate", "PROSAIL_Models",
    site, name_strategy, "atbd",
    paste0("LAI_estimated_atbd_", sampling_method,
           "_hmin", h_min_pixel,
           "_nbSamples_", nb_samples, ".csv")
  )
  if (!file.exists(gpkg_path)) {
    cli::cli_warn("GPKG not found, skipping {site}: {gpkg_path}")
    return(NULL)
  }
  if (!file.exists(s2_csv)) {
    cli::cli_warn("S2 CSV not found, skipping {site}: {s2_csv}")
    return(NULL)
  }

  # Read max_height per sample from GPKG
  samp_vect  <- terra::vect(gpkg_path)
  sample_ids <- samp_vect$sample_id
  max_h      <- samp_vect$max_height
  n          <- length(sample_ids)

  # S2 LAI estimates (row order = sample 1..n)
  s2_dt   <- data.table::fread(s2_csv, header = TRUE, sep = "\t")
  s2_vals <- s2_dt[["LAI_atbd"]]

  # Load all depth CSVs → pad_mat (n × 38)
  pad_mat <- matrix(NA_real_, nrow = n, ncol = length(depths))
  colnames(pad_mat) <- as.character(depths)

  for (depth in depths) {
    lidar_csv <- here::here(
      "03_RESULTS", site, "PROSAIL_Optimization", "sampling",
      paste0("PAD_", norm, "_Depth_", depth,
             "_Samples_", sampling_method,
             "_hmin", h_min_pixel,
             "_nbSamples_", nb_samples, ".csv")
    )
    if (!file.exists(lidar_csv)) next
    lidar_dt  <- data.table::fread(lidar_csv, header = TRUE, sep = "\t")
    # Align by sample_id — lidar CSV rows may not be in the same order as GPKG
    idx <- match(sample_ids, lidar_dt[["samples_id"]])
    pad_mat[, depth] <- lidar_dt[["lidar_values"]][idx]
  }

  list(pad_mat = pad_mat, max_h = max_h, s2_vals = s2_vals,
       sample_ids = sample_ids)
}

# ── Helper: compute metrics data.table from a PAD matrix ──────────────────────
# Returns a data.table with columns Site, Norm, Depth, R, R2, RMSE, Bias,
# Slope, ATBD, Column, plus any extra columns passed via `extra`.

compute_metrics_from_mat <- function(pad_mat, s2_vals, site, norm,
                                      extra = list()) {
  rows <- vector("list", length(depths))
  for (depth in depths) {
    pad_d <- pad_mat[, depth]
    valid <- !is.na(pad_d) & !is.na(s2_vals)
    if (sum(valid) < 5L) next
    m <- compute_metrics_s2_lidar(s2_vals[valid], pad_d[valid])
    row <- data.table::data.table(
      Site   = site, Norm = norm, Depth = as.integer(depth),
      R = round(m$R, 3L), R2 = round(m$R2, 3L),
      RMSE = round(m$RMSE, 3L), Bias = round(m$Bias, 3L),
      Slope = round(m$Slope, 3L), ATBD = TRUE, Column = "LAI_atbd"
    )
    for (nm in names(extra)) row[[nm]] <- extra[[nm]]
    rows[[depth]] <- row
  }
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

# ── Main computation loop ─────────────────────────────────────────────────────

cli::cli_h1("Loading data and computing metrics...")

k_rows    <- list()
zmin_rows <- list()

for (site in sites) {
  for (norm in norms_select) {

    cli::cli_alert_info("{site} / {norm}")
    dat <- load_site_norm_data(site, norm)
    if (is.null(dat)) next

    # ── k sensitivity (z_min fixed at 2 m — existing PAD baseline) ────────────
    for (k_val in k_values) {
      scaled_mat <- dat$pad_mat * (k_ref / k_val)
      rows <- compute_metrics_from_mat(
        scaled_mat, dat$s2_vals, site, norm,
        extra = list(k = k_val, z_min = 2)
      )
      k_rows <- c(k_rows, list(rows))
    }

    # ── z_min sensitivity (k fixed at k_ref = 0.5) ────────────────────────────
    for (z_min in z_min_values) {
      # Per-pixel cap: depth beyond which the pixel's canopy is below z_min
      cap_depth <- as.integer(pmax(1L, floor(dat$max_h - z_min)))

      # Build corrected PAD matrix: for each pixel and depth d,
      # use PAD at min(d, cap_depth_i) — caps at the last valid layer
      corrected_mat <- dat$pad_mat
      for (depth in depths) {
        eff_d <- pmin(depth, cap_depth)          # per pixel, integer 1..38
        # Replace column with per-pixel value at effective depth
        corrected_mat[, depth] <- dat$pad_mat[cbind(seq_len(nrow(dat$pad_mat)),
                                                      eff_d)]
      }

      rows <- compute_metrics_from_mat(
        corrected_mat, dat$s2_vals, site, norm,
        extra = list(k = k_ref, z_min = z_min)
      )
      zmin_rows <- c(zmin_rows, list(rows))
    }
  }
}

k_dt    <- data.table::rbindlist(k_rows,    use.names = TRUE, fill = TRUE)
zmin_dt <- data.table::rbindlist(zmin_rows, use.names = TRUE, fill = TRUE)

# ── d_opt (Pareto) for k sensitivity ──────────────────────────────────────────

cli::cli_h1("d_opt sensitivity to k (z_min = 2 m, h_min_pixel = {h_min_pixel} m)")

dopt_k_rows <- list()
for (k_val in k_values) {
  sub <- k_dt[k == k_val & Norm %in% norms_select]
  if (nrow(sub) == 0L) next
  d <- select_dopt(sub, methods = "pareto", max_depth = h_min_pixel,
                   prosail_filter = "ATBD")
  d[method_dopt == "pareto"][, k := k_val]
  dopt_k_rows <- c(dopt_k_rows, list(d[method_dopt == "pareto"]))
}
dopt_k <- data.table::rbindlist(dopt_k_rows, use.names = TRUE, fill = TRUE)

wide_k <- data.table::dcast(
  dopt_k, Site + Norm ~ paste0("k=", k),
  value.var = "d_opt"
)
data.table::setorder(wide_k, Norm, Site)
print(wide_k)

# ── d_opt (Pareto) for z_min sensitivity ──────────────────────────────────────

cli::cli_h1("d_opt sensitivity to z_min (k = {k_ref}, h_min_pixel = {h_min_pixel} m)")

dopt_zmin_rows <- list()
for (z_min_val in z_min_values) {
  sub <- zmin_dt[z_min == z_min_val & Norm %in% norms_select]
  if (nrow(sub) == 0L) next
  d <- select_dopt(sub, methods = "pareto", max_depth = h_min_pixel,
                   prosail_filter = "ATBD")
  d[method_dopt == "pareto", z_min := z_min_val]
  dopt_zmin_rows <- c(dopt_zmin_rows, list(d[method_dopt == "pareto"]))
}
dopt_zmin <- data.table::rbindlist(dopt_zmin_rows, use.names = TRUE, fill = TRUE)

wide_zmin <- data.table::dcast(
  dopt_zmin, Site + Norm ~ paste0("zmin=", z_min),
  value.var = "d_opt"
)
data.table::setorder(wide_zmin, Norm, Site)
print(wide_zmin)

cat("Done.\n")
