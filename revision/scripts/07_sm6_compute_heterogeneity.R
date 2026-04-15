# ---
# title:  07_sm6_compute_heterogeneity.R
# desc:   Orchestration — computes DSM_std and CHM_std heterogeneity rasters
#         at 10 m resolution for SM6, masked to deciduous-only pixels, and
#         writes them to revision/output/intermediate/sm6/{site}/.
#
#         Calls compute_heterogeneity() and save_heterogeneity_raster() from
#         revision/R/sm6_heterogeneity.R.
#
#         Mask used: artifacts_deciduous_only_low_vegetation_majority_90_p_res_10_m.envi
#         (pre-computed per-site raster in 03_RESULTS/{site}/LiDAR/Heterogeneity_Masks/).
#         This combines the BDForêt deciduous-only filter, low-vegetation exclusion
#         (fCover majority 90 %), and artifact removal — equivalent to the
#         "Deciduous_Only" branch of apply_and_save_masks() in the legacy pipeline.
#
#         Semantic equivalent of the legacy output:
#           03_RESULTS/{site}/Metrics/Deciduous_Only/dsm_sd_res_10_m.tif
#         (reference run: 2025-07-15). The legacy produced this in two steps
#         (writeRaster to Raw/, then apply_and_save_masks to Deciduous_Only/);
#         here aggregate + project + mask are fused in a single call to
#         compute_heterogeneity(), using the Deciduous_Only mask directly.
#
#         CHM source: chm/res_1_m/chm.tif (standard, not pit-free).
#         A pit-free CHM exists for Aigoual but not for Blois/Mormal; the
#         standard CHM is used for all three sites to ensure inter-site
#         consistency.
#
#         Prerequisites (all three sites must have):
#           DSM 1 m: 03_RESULTS/{site}/LiDAR/dsm/res_1_m/rasterize_canopy.vrt
#           CHM 1 m: 03_RESULTS/{site}/LiDAR/chm/res_1_m/chm.tif
#           Mask 10m: 03_RESULTS/{site}/LiDAR/Heterogeneity_Masks/
#                     artifacts_deciduous_only_low_vegetation_majority_90_p_res_10_m.envi
#
#         Addresses reviewer comment R2.8: CHM_std vs DSM_std comparison.
#         Both metrics are computed with the same pipeline (fact = 10, fun = "sd").
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/07_sm6_compute_heterogeneity.R")
# ---

source(here::here("revision", "R", "sm6_heterogeneity.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites   <- c("Aigoual", "Blois", "Mormal")
metrics <- c("DSM", "CHM")
out_dir <- here::here("revision", "output", "intermediate", "sm6")

# ── Input path helpers ─────────────────────────────────────────────────────────

source_path <- function(site, metric) {
  if (metric == "DSM") {
    here::here("03_RESULTS", site, "LiDAR", "dsm", "res_1_m",
               "rasterize_canopy.vrt")
  } else {
    here::here("03_RESULTS", site, "LiDAR", "chm", "res_1_m", "chm.tif")
  }
}

mask_path <- function(site) {
  here::here("03_RESULTS", site, "LiDAR", "Heterogeneity_Masks",
             "artifacts_deciduous_only_low_vegetation_majority_90_p_res_10_m.envi")
}

# ── Pre-flight: verify all input files exist before any computation ────────────

missing_files <- character(0)

for (site in sites) {
  for (metric in metrics) {
    f <- source_path(site, metric)
    if (!file.exists(f)) missing_files <- c(missing_files, f)
  }
  m <- mask_path(site)
  if (!file.exists(m)) missing_files <- c(missing_files, m)
}

if (length(missing_files) > 0L) {
  stop(
    "SM6: missing input file(s) — cannot proceed:\n",
    paste0("  ", missing_files, collapse = "\n")
  )
}

# ── Output directory ───────────────────────────────────────────────────────────

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# ── Compute and save ───────────────────────────────────────────────────────────

t0         <- proc.time()
out_paths  <- character(0)

for (site in sites) {
  cli::cli_alert_info("Site: {site}")

  for (metric in metrics) {
    cli::cli_alert_info("  Metric: {metric}")

    het <- compute_heterogeneity(
      source_raster_path = source_path(site, metric),
      mask_raster_path   = mask_path(site)
    )

    out_path <- save_heterogeneity_raster(
      het_raster = het,
      site       = site,
      metric     = metric,
      out_dir    = out_dir
    )

    out_paths <- c(out_paths, out_path)
    cli::cli_alert_success("    Written: {out_path}")
  }
}

elapsed <- round((proc.time() - t0)[["elapsed"]], 1)

# ── Console summary ────────────────────────────────────────────────────────────

file_sizes_kb <- round(file.info(out_paths)[["size"]] / 1024, 1)

cli::cli_h1("SM6 — résumé")
cli::cli_bullets(c(
  "v" = "Fichiers écrits : {length(out_paths)}",
  "i" = "Durée : {elapsed} s"
))
for (i in seq_along(out_paths)) {
  cli::cli_bullets(c("i" = "{basename(out_paths[i])} — {file_sizes_kb[i]} KB"))
}
