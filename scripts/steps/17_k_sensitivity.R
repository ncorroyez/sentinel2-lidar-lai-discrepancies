# ---
# title:  17_k_sensitivity.R
# desc:   Orchestration — full factorial sensitivity of LAI_ALS vs LAI_S2_ATBD
#         to k × theta (42 combinations = 6 k values × 7 scan angles).
#
#         Combined rescaling (exact, no approximation):
#           LAI_ALS(k, θ) = LAI_ALS(k_ref, 0°) × (k_ref / k) × cos(θ)
#         where k_ref = 0.5 (Bouvier et al. 2015) and LAI_ALS(k_ref, 0°) is
#         the integrated LAD sum from the ladstack_classic.tif raster.
#
#         Addresses reviewer comment R3.MAJOR.
#         Analysis limited to: LAI_S2_ATBD vs LAI_ALS (integrated LAI),
#         3 sites, no heterogeneity stratification, no LAI_ALS_dopt.
#
#         Output:
#           revision/output/intermediate/reviewers/k_theta_sensitivity_atbd.csv
#           126 rows (3 sites × 42 combinations), 17 columns:
#             site, k_value, theta_deg, n_pixels,
#             lai_als_{mean,median,p95,max},
#             lai_s2_atbd_{mean,median,p95,max},
#             R, R2, RMSE, Bias, Slope
#
#         Console summary: one k × theta matrix per metric per site (15 tables).
#
#         Input rasters (Deciduous_Only, produced by legacy pipeline):
#           03_RESULTS/{site}/Metrics/Deciduous_Only/ladstack_classic.tif
#           03_RESULTS/{site}/Metrics/Deciduous_Only/s2lai_summer_atbd_res_10_m.tif
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/17_k_sensitivity.R")
# ---

library(here)
library(terra)
library(data.table)
library(cli)

source(here::here("revision", "R", "paths.R"))
source(here::here("revision", "R", "k_sensitivity.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites        <- c("Aigoual", "Blois", "Mormal")
k_values     <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
k_ref        <- 0.5                   # reference k embedded in ladstack_classic.tif
theta_values <- c(0, 5, 10, 15, 20, 25, 30)  # degrees from nadir

ext_dir  <- paths$ext_results
out_dir  <- file.path(paths$output, "intermediate", "reviewers")

# ── Input path helpers ─────────────────────────────────────────────────────────

ladstack_path <- function(site) {
  file.path(ext_dir, site, "Metrics", "Deciduous_Only", "ladstack_classic.tif")
}

sm6a_dir <- file.path(paths$output, "intermediate", "sm6")

s2_atbd_path <- function(site) {
  file.path(sm6a_dir, site, "s2lai_summer_atbd_T_res_10_m.tif")
}

# ── Pre-flight: verify all required input files ─────────────────────────────────

required_files <- list()

for (site in sites) {
  required_files[[paste0(site, "_ladstack")]] <- list(
    path     = ladstack_path(site),
    producer = "produced by 2.calculate_25m_metrics.R"
  )
  required_files[[paste0(site, "_s2_atbd")]] <- list(
    path     = s2_atbd_path(site),
    producer = "produced by 3_train_predict_prosail.R"
  )
}

missing_msgs <- character(0)
for (key in names(required_files)) {
  entry <- required_files[[key]]
  if (!file.exists(entry$path)) {
    missing_msgs <- c(
      missing_msgs,
      paste0("  MISSING [", entry$producer, "]\n    ", entry$path)
    )
  }
}

if (length(missing_msgs) > 0L) {
  stop(
    "k_sensitivity pre-flight failed — ", length(missing_msgs),
    " file(s) missing:\n",
    paste(missing_msgs, collapse = "\n"),
    call. = FALSE
  )
}

cli::cli_alert_success(
  "Pre-flight passed — all {length(required_files)} required files found."
)

# ── Output directory ───────────────────────────────────────────────────────────

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# ── Main loop: site × k × theta (full factorial) ──────────────────────────────

t0       <- proc.time()
all_rows <- list()
n_combo  <- length(k_values) * length(theta_values)

cli::cli_h2(
  "k × theta factorial — {length(sites)} sites × {n_combo} combinations ({length(k_values)} k × {length(theta_values)} theta)"
)

for (site in sites) {
  cli::cli_h3("Site: {site}")
  lai_s2_atbd_rast <- terra::rast(s2_atbd_path(site))

  # Load reference raster once per site (k_ref, theta=0°), then scale
  lai_als_ref <- compute_lai_als_at_k(
    ladstack_path = ladstack_path(site),
    k_new         = k_ref,
    k_ref         = k_ref
  )

  for (k in k_values) {
    for (theta in theta_values) {
      scale_factor <- (k_ref / k) * cos(theta * pi / 180)
      lai_als_rast <- lai_als_ref * scale_factor

      row <- compute_k_sensitivity_metrics(
        lai_als_rast     = lai_als_rast,
        lai_s2_atbd_rast = lai_s2_atbd_rast,
        site             = site,
        k_value          = k,
        theta_deg        = theta
      )
      all_rows <- c(all_rows, list(row))
    }
  }
  cli::cli_alert_success("{site} done — {n_combo} combinations computed")
}

# ── Combine, reorder columns, write ───────────────────────────────────────────

final_dt <- data.table::rbindlist(all_rows)

data.table::setcolorder(final_dt, c(
  "site", "k_value", "theta_deg", "n_pixels",
  "lai_als_mean", "lai_als_median", "lai_als_p95", "lai_als_max",
  "lai_s2_atbd_mean", "lai_s2_atbd_median", "lai_s2_atbd_p95", "lai_s2_atbd_max",
  "R", "R2", "RMSE", "Bias", "Slope"
))

out_csv <- file.path(out_dir, "k_theta_sensitivity_atbd.csv")
data.table::fwrite(final_dt, out_csv)

elapsed <- round((proc.time() - t0)[["elapsed"]], 1)

# ── Console summary ────────────────────────────────────────────────────────────

cli::cli_h1("k × theta sensitivity — résumé")
cli::cli_bullets(c(
  "v" = "Durée : {elapsed} s",
  "v" = "CSV   : {out_csv}  ({nrow(final_dt)} lignes — {length(sites)} sites × {n_combo} combos)"
))

# One k × theta matrix per metric per site
metrics_display <- c("R", "R2", "RMSE", "Bias", "Slope")

for (s in sites) {
  cli::cli_h2("Site: {s}")
  site_dt <- final_dt[site == s]

  for (metric in metrics_display) {
    wide <- data.table::dcast(
      site_dt, k_value ~ paste0("theta=", theta_deg),
      value.var = metric,
      fun.aggregate = function(x) round(x, 3)
    )
    # Sort theta columns numerically
    theta_cols <- paste0("theta=", sort(theta_values))
    data.table::setcolorder(wide, c("k_value", theta_cols))

    cat("\n  ", metric, ":\n")
    print(wide, row.names = FALSE)
  }
}
