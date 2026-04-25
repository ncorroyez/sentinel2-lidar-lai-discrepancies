# ---
# title:  07_sm6_compute_heterogeneity.R
# desc:   Orchestration — SM6a passe 2. Computes DSM_std and CHM_std
#         heterogeneity rasters at multiple spatial scales (block and focal)
#         for all three sites, masked to deciduous-only pixels.
#
#         Sweeps:
#           block_facts = c(5, 10, 20, 50): aggregate(fact, fun="sd") on the
#             1 m source → block windows of 5, 10, 20, 50 m.
#           focal_Ks    = c(3, 5): focal(w=K, fun="sd") on the block_10_m base
#             → sliding windows of 30, 50 m (effective K × 10 m).
#             K=1 (focal_10_m) is excluded: terra::focal rejects w=1 as
#             non-meaningful. Its semantic equivalent is block_10_m (see
#             compute_heterogeneity roxygen). focal_Ks = c(3L, 5L) covers
#             effective windows of 30 m and 50 m.
#
#         Output files per (site × metric):
#           {metric}_sd_block_5_m.tif
#           {metric}_sd_block_10_m.tif
#           {metric}_sd_block_20_m.tif
#           {metric}_sd_block_50_m.tif
#           {metric}_sd_focal_30_m.tif
#           {metric}_sd_focal_50_m.tif
#           {metric}_sd_res_10_m.tif   ← legacy name, identical to block_10_m
#                                         (rétrocompatibilité SM6b pipeline)
#         Total: 7 files × 3 sites × 2 metrics = 42 writes.
#
#         Output root: revision/output/intermediate/sm6/{site}/
#
#         Mask used:
#           artifacts_deciduous_only_low_vegetation_majority_90_p_res_10_m.envi
#           (pre-computed per-site raster in
#           03_RESULTS/{site}/LiDAR/Heterogeneity_Masks/).
#           Combines BDForêt deciduous filter, low-vegetation exclusion
#           (fCover majority 90 %), and artifact removal.
#
#         CHM source: chm/res_1_m/chm.tif (standard, not pit-free).
#         Pit-free exists for Aigoual only; standard CHM used for all sites
#         to ensure inter-site consistency.
#
#         Addresses reviewer comment R2.8: CHM_std vs DSM_std comparison,
#         extended to multi-scale sweep for SM6a passe 2.
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/07_sm6_compute_heterogeneity.R")
# ---

library(here)
library(terra)
library(cli)

source(here::here("revision", "R", "paths.R"))
source(here::here("revision", "R", "sm6_heterogeneity.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites   <- c("Aigoual", "Blois", "Mormal")
metrics <- c("DSM", "CHM")
out_dir <- here::here("revision", "output", "intermediate", "sm6")

# Block window sizes in metres (source is 1 m → fact = size_m)
block_facts <- c(5L, 10L, 20L, 50L)

# Focal window sizes in 10 m pixels; effective window = K × 10 m
# K=1 excluded (terra::focal rejects w=1). block_10_m is its semantic
# equivalent — see header and compute_heterogeneity roxygen for details.
focal_Ks <- c(3L, 5L)

# ── Input path helpers ─────────────────────────────────────────────────────────

source_path <- function(site, metric) {
  if (metric == "DSM") {
    file.path(paths$ext_results, site, "LiDAR", "dsm", "res_1_m",
               "rasterize_canopy.vrt")
  } else {
    file.path(paths$ext_results, site, "LiDAR", "chm", "res_1_m", "chm.tif")
  }
}

mask_path <- function(site) {
  file.path(paths$ext_results, site, "LiDAR", "Heterogeneity_Masks",
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
    "SM6a: missing input file(s) — cannot proceed:\n",
    paste0("  ", missing_files, collapse = "\n")
  )
}

cli::cli_alert_success(
  "Pre-flight passed — all {length(sites) * length(metrics) + length(sites)} ",
  "input files found."
)

# ── Output directory ───────────────────────────────────────────────────────────

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# ── Compute and save ───────────────────────────────────────────────────────────

t0        <- proc.time()
out_paths <- character(0)   # new naming convention
leg_paths <- character(0)   # legacy naming (block_10_m duplicates)

n_block_expected <- length(block_facts) * length(sites) * length(metrics)
n_focal_expected <- length(focal_Ks)    * length(sites) * length(metrics)
n_leg_expected   <-                        length(sites) * length(metrics)

cli::cli_h2(
  "SM6a passe 2 — {n_block_expected} block + {n_focal_expected} focal + ",
  "{n_leg_expected} legacy = ",
  "{n_block_expected + n_focal_expected + n_leg_expected} total writes"
)

for (site in sites) {
  cli::cli_h3("Site: {site}")
  msk <- mask_path(site)

  for (metric in metrics) {
    cli::cli_alert_info("  Metric: {metric}")
    src <- source_path(site, metric)
    het_block10 <- NULL     # caches block_10_m for legacy + focal reuse

    # ── Block variants ─────────────────────────────────────────────────────────
    for (fact in block_facts) {
      size_m <- fact          # 1 m source → fact m aggregation window
      fn     <- build_output_filename(metric, "block", size_m)

      het <- compute_heterogeneity(
        source_raster_path = src,
        mask_raster_path   = msk,
        fact               = fact,
        method             = "block"
      )

      out_path  <- save_heterogeneity_raster(het, site, fn, out_dir)
      out_paths <- c(out_paths, out_path)
      cli::cli_alert_success("    {fn}")

      if (fact == 10L) het_block10 <- het   # save for legacy write + focal
    }

    # ── Legacy file: block_10_m with passe-1 naming (SM6b compatibility) ──────
    legacy_fn <- paste0(tolower(metric), "_sd_res_10_m.tif")
    leg_path  <- save_heterogeneity_raster(het_block10, site, legacy_fn, out_dir)
    leg_paths <- c(leg_paths, leg_path)
    cli::cli_alert_success("    {legacy_fn}  [legacy = block_10_m]")

    # ── Focal variants (sliding window SD on block_10_m base) ─────────────────
    for (K in focal_Ks) {
      size_m <- K * 10L      # K pixels on 10 m grid → effective window size
      fn     <- build_output_filename(metric, "focal", size_m)

      het <- compute_heterogeneity(
        source_raster_path = src,
        mask_raster_path   = msk,
        method             = "focal",
        K                  = K,
        base_raster        = het_block10   # reuse cached block_10_m
      )

      out_path  <- save_heterogeneity_raster(het, site, fn, out_dir)
      out_paths <- c(out_paths, out_path)
      cli::cli_alert_success("    {fn}")
    }
  }
}

elapsed <- round((proc.time() - t0)[["elapsed"]], 1)

# ── Console summary ────────────────────────────────────────────────────────────

all_written   <- c(out_paths, leg_paths)
file_sizes_kb <- round(file.info(all_written)[["size"]] / 1024, 1)

cli::cli_h1("SM6a passe 2 — résumé")
cli::cli_bullets(c(
  "v" = "Durée : {elapsed} s",
  "v" = "Nouveaux fichiers (convention passe 2) : {length(out_paths)}",
  "v" = "Fichiers legacy (rétrocompat SM6b)     : {length(leg_paths)}"
))
cli::cli_text("Détail par fichier :")
for (i in seq_along(all_written)) {
  label <- if (i > length(out_paths)) " [legacy]" else ""
  cli::cli_bullets(c(
    "i" = "{basename(all_written[i])}{label} — {file_sizes_kb[i]} KB"
  ))
}
