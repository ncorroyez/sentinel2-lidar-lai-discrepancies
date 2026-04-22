# ---
# title:  10_k_sensitivity.R
# desc:   Orchestration — LAI_ALS sensitivity to extinction coefficient k.
#         Computes how LAI_ALS descriptive statistics and LAI_S2_ATBD vs
#         LAI_ALS comparison metrics vary with k ∈ {0.4, 0.5, 0.6}, using
#         the analytical rescaling LAI(k_new) = LAI(k=0.5) × 0.5 / k_new
#         (Bouvier et al. 2015). No pipeline rerun; no PROSAIL LUT retraining.
#
#         Addresses reviewer comment R3.MAJOR.
#         Analysis limited to: LAI_S2_ATBD vs LAI_ALS (integrated LAI),
#         3 sites, no heterogeneity stratification, no LAI_ALS_dopt.
#
#         Output:
#           revision/output/intermediate/reviewers/k_sensitivity_atbd.csv
#           9 rows (3 sites × 3 k values), 15 columns:
#             site, k_value, n_pixels,
#             lai_als_{mean,median,p95,max},
#             lai_s2_atbd_{mean,median,p95,max},
#             R, R2, RMSE, Bias, Slope
#
#         Input rasters (Deciduous_Only, produced by legacy pipeline):
#           03_RESULTS/{site}/Metrics/Deciduous_Only/ladstack_classic.tif
#           03_RESULTS/{site}/Metrics/Deciduous_Only/s2lai_summer_atbd_res_10_m.tif
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/10_k_sensitivity.R")
# ---

library(here)
library(terra)
library(data.table)
library(cli)

source(here::here("revision", "R", "k_sensitivity.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites    <- c("Aigoual", "Blois", "Mormal")
k_values <- c(0.4, 0.5, 0.6)    # extinction coefficient values to test
k_ref    <- 0.5                  # reference k embedded in ladstack_classic.tif

ext_dir  <- here::here("03_RESULTS")
out_dir  <- here::here("revision", "output", "intermediate", "reviewers")

# ── Input path helpers ─────────────────────────────────────────────────────────

ladstack_path <- function(site) {
  file.path(ext_dir, site, "Metrics", "Deciduous_Only", "ladstack_classic.tif")
}

s2_atbd_path <- function(site) {
  file.path(ext_dir, site, "Metrics", "Deciduous_Only",
            "s2lai_summer_atbd_res_10_m.tif")
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

# ── Main loop: site × k ────────────────────────────────────────────────────────

t0       <- proc.time()
n_total  <- length(sites) * length(k_values)
all_rows <- vector("list", n_total)
row_idx  <- 0L

cli::cli_h2(
  "k sensitivity — {n_total} combinations ",
  "({length(sites)} sites × {length(k_values)} k values)"
)

for (site in sites) {
  cli::cli_h3("Site: {site}")

  # LAI_S2_ATBD is constant across k — load once per site
  lai_s2_atbd_rast <- terra::rast(s2_atbd_path(site))

  for (k in k_values) {
    row_idx <- row_idx + 1L
    cli::cli_alert_info("  k = {k}")

    lai_als_rast <- compute_lai_als_at_k(
      ladstack_path = ladstack_path(site),
      k_new         = k,
      k_ref         = k_ref
    )

    all_rows[[row_idx]] <- compute_k_sensitivity_metrics(
      lai_als_rast      = lai_als_rast,
      lai_s2_atbd_rast  = lai_s2_atbd_rast,
      site              = site,
      k_value           = k
    )

    cli::cli_alert_success(
      "    n={all_rows[[row_idx]]$n_pixels}  ",
      "LAI_ALS_mean={round(all_rows[[row_idx]]$lai_als_mean, 3)}  ",
      "R2={round(all_rows[[row_idx]]$R2, 3)}"
    )
  }
}

# ── Combine, reorder columns, write ───────────────────────────────────────────

final_dt <- data.table::rbindlist(all_rows)

data.table::setcolorder(final_dt, c(
  "site", "k_value", "n_pixels",
  "lai_als_mean", "lai_als_median", "lai_als_p95", "lai_als_max",
  "lai_s2_atbd_mean", "lai_s2_atbd_median", "lai_s2_atbd_p95", "lai_s2_atbd_max",
  "R", "R2", "RMSE", "Bias", "Slope"
))

out_csv <- file.path(out_dir, "k_sensitivity_atbd.csv")
data.table::fwrite(final_dt, out_csv)

elapsed <- round((proc.time() - t0)[["elapsed"]], 1)

# ── Console summary ────────────────────────────────────────────────────────────

cli::cli_h1("k sensitivity — résumé")
cli::cli_bullets(c(
  "v" = "Durée : {elapsed} s",
  "v" = "CSV   : {out_csv}  ({nrow(final_dt)} lignes)"
))
cli::cli_text("Résultats (pivot site × k) :")
print(
  final_dt[, .(site, k_value, n_pixels,
               lai_als_mean = round(lai_als_mean, 3),
               R2           = round(R2, 3),
               RMSE         = round(RMSE, 3),
               Bias         = round(Bias, 3),
               Slope        = round(Slope, 3))]
)
