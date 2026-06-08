# ---
# title:  05c_sm5_k_zmin_sensitivity_dopt.R
# desc:   Sensitivity of Pareto d_opt to:
#           - all (k, theta) combinations: k ∈ {0.3..0.8}, theta ∈ {0..30°}
#             Combined rescaling: PAD(k, θ) = PAD_ref × (k_ref/k) × cos(θ)
#           - minimum height z_min (2, 3, 5 m): per-pixel canopy depth cap
#             (k = k_ref, theta = 0°)
#         All sensitivities are computed for each individual site and for
#         "Sites_combined" (all 3 sites pooled), for DSM_keepTrees and
#         DTM_keepTrees.
#
#         Combined rescaling: PAD(k, θ) = PAD_ref × (k_ref / k) × cos(θ)
#         z_min cap: for depth d, effective depth = min(d, floor(max_h - z_min))
#
#         h_min_pixel is fixed at 10 m (reference analysis value).
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

source(here::here("revision", "R", "paths.R"))
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
theta_values    <- c(0, 5, 10, 15, 20, 25, 30)  # scan angles in degrees

out_dir <- file.path(paths$output, "intermediate", "sm5")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Helper: build PAD wide matrix + max_height for one (site, norm) ───────────
# Returns a list:
#   $pad_mat  : n × 38 numeric matrix (row = sample, col = depth)
#   $max_h    : n-vector of max_height per sample (matched to row order)
#   $s2_vals  : n-vector of ATBD LAI_S2 estimates
#   $sample_ids: n-vector of sample IDs

load_site_norm_data <- function(site, norm) {

  gpkg_path <- file.path(
    paths$ext_results, site, "PROSAIL_Optimization", "sampling",
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
    lidar_csv <- file.path(
      paths$ext_results, site, "PROSAIL_Optimization", "sampling",
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

# ── Step 1: load all per-site data ────────────────────────────────────────────

cli::cli_h1("Loading data...")

all_dat <- list()   # keyed by "<site> <norm>"

for (site in sites) {
  for (norm in norms_select) {
    cli::cli_alert_info("{site} / {norm}")
    dat <- load_site_norm_data(site, norm)
    all_dat[[paste(site, norm)]] <- dat
  }
}

# ── Step 2: build Sites_combined (pixel pool, always computed) ────────────────

cli::cli_h2("Building Sites_combined (pixel pool)...")
for (norm in norms_select) {
  site_dats <- Filter(Negate(is.null),
                      lapply(sites, function(s) all_dat[[paste(s, norm)]]))
  if (length(site_dats) == 0L) next
  combined <- list(
    pad_mat = do.call(rbind, lapply(site_dats, `[[`, "pad_mat")),
    max_h   = do.call(c,    lapply(site_dats, `[[`, "max_h")),
    s2_vals = do.call(c,    lapply(site_dats, `[[`, "s2_vals"))
  )
  all_dat[[paste("Sites_combined", norm)]] <- combined
  cli::cli_alert_success(
    "Sites_combined / {norm}: {nrow(combined$pad_mat)} samples"
  )
}
all_site_labels <- c(sites, "Sites_combined")

# ── Step 3: compute sensitivity metrics ───────────────────────────────────────

cli::cli_h1("Computing sensitivity metrics...")

kt_rows   <- list()   # all (k, theta) combinations
zmin_rows <- list()

n_kt <- length(k_values) * length(theta_values)

for (site_label in all_site_labels) {
  for (norm in norms_select) {

    dat <- all_dat[[paste(site_label, norm)]]
    if (is.null(dat)) next

    cli::cli_h3("{site_label} / {norm}")

    # ── k × theta full factorial (z_min = 2 m) ────────────────────────────────
    cli::cli_alert_info("{n_kt} (k, theta) combinations...")
    for (k_val in k_values) {
      for (theta_val in theta_values) {
        scale_factor <- (k_ref / k_val) * cos(theta_val * pi / 180)
        scaled_mat   <- dat$pad_mat * scale_factor
        rows <- compute_metrics_from_mat(
          scaled_mat, dat$s2_vals, site_label, norm,
          extra = list(k = k_val, theta = theta_val, z_min = 2)
        )
        kt_rows <- c(kt_rows, list(rows))
      }
    }

    # ── z_min sensitivity (k = k_ref, theta = 0°) ────────────────────────────
    for (z_min in z_min_values) {
      cap_depth <- as.integer(pmax(1L, floor(dat$max_h - z_min)))

      corrected_mat <- dat$pad_mat
      for (depth in depths) {
        eff_d <- pmin(depth, cap_depth)
        corrected_mat[, depth] <- dat$pad_mat[cbind(seq_len(nrow(dat$pad_mat)),
                                                      eff_d)]
      }

      rows <- compute_metrics_from_mat(
        corrected_mat, dat$s2_vals, site_label, norm,
        extra = list(k = k_ref, theta = 0, z_min = z_min)
      )
      zmin_rows <- c(zmin_rows, list(rows))
    }
  }
}

kt_dt   <- data.table::rbindlist(kt_rows,   use.names = TRUE, fill = TRUE)
zmin_dt <- data.table::rbindlist(zmin_rows, use.names = TRUE, fill = TRUE)

# ── Sites_averaged: average per-site metrics across the 3 sites (always) ──────

cli::cli_h2("Computing Sites_averaged (mean metrics across sites)...")

avg_by <- function(dt, grp_vars) {
  dt[Site %in% sites,
    .(R     = mean(R,     na.rm = TRUE),
      R2    = mean(R2,    na.rm = TRUE),
      RMSE  = mean(RMSE,  na.rm = TRUE),
      Bias  = mean(Bias,  na.rm = TRUE),
      Slope = mean(Slope, na.rm = TRUE),
      ATBD  = TRUE, Column = "LAI_atbd"),
    by = grp_vars
  ][, Site := "Sites_averaged"]
}

kt_dt   <- data.table::rbindlist(
  list(kt_dt,   avg_by(kt_dt,   c("Norm", "Depth", "k", "theta", "z_min"))),
  use.names = TRUE, fill = TRUE
)
zmin_dt <- data.table::rbindlist(
  list(zmin_dt, avg_by(zmin_dt, c("Norm", "Depth", "k", "theta", "z_min"))),
  use.names = TRUE, fill = TRUE
)

cli::cli_alert_success(
  "Sites_averaged appended — kt_dt: {nrow(kt_dt)} rows, zmin_dt: {nrow(zmin_dt)} rows"
)

# ── Write kt_dt to CSV (read by 05b / 06 for custom k and theta) ──────────────

kt_csv <- file.path(out_dir, "kt_sensitivity_atbd_LIDFa_lai_LMA_BROWN.csv")
data.table::fwrite(kt_dt, kt_csv)
cli::cli_alert_success("Written: {kt_csv}  ({nrow(kt_dt)} rows)")

# ── d_opt (Pareto) for k × theta factorial ────────────────────────────────────

cli::cli_h1("d_opt sensitivity to (k, theta) — {n_kt} combinations, h_min_pixel = {h_min_pixel} m")

dopt_kt_rows <- list()
for (k_val in k_values) {
  for (theta_val in theta_values) {
    sub <- kt_dt[k == k_val & theta == theta_val & Norm %in% norms_select]
    if (nrow(sub) == 0L) next
    d <- select_dopt(sub, methods = "pareto", max_depth = h_min_pixel,
                     prosail_filter = "ATBD")
    d[method_dopt == "pareto", `:=`(k = k_val, theta = theta_val)]
    dopt_kt_rows <- c(dopt_kt_rows, list(d[method_dopt == "pareto"]))
  }
}
dopt_kt <- data.table::rbindlist(dopt_kt_rows, use.names = TRUE, fill = TRUE)

wide_kt <- data.table::dcast(
  dopt_kt, Site + Norm ~ paste0("k=", k, "_theta=", theta),
  value.var = "d_opt"
)
data.table::setorder(wide_kt, Norm, Site)
print(wide_kt)

dopt_kt_csv <- file.path(out_dir, "dopt_k_theta_sensitivity.csv")
data.table::fwrite(dopt_kt, dopt_kt_csv)
cli::cli_alert_success("Written: {dopt_kt_csv}  ({nrow(dopt_kt)} rows)")

# ── d_opt (Pareto) for z_min sensitivity ──────────────────────────────────────

cli::cli_h1("d_opt sensitivity to z_min (k = {k_ref}, theta = 0°, h_min_pixel = {h_min_pixel} m)")

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

dopt_zmin_csv <- file.path(out_dir, "dopt_zmin_sensitivity.csv")
data.table::fwrite(dopt_zmin, dopt_zmin_csv)
cli::cli_alert_success("Written: {dopt_zmin_csv}  ({nrow(dopt_zmin)} rows)")

cat("Done.\n")
