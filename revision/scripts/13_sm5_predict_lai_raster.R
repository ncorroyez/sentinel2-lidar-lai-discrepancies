# ---
# title:  06c_sm5_predict_lai_raster.R
# desc:   Orchestration — SM5 passe 3.
#         For each site and each d_opt_source scenario, applies the Pareto-
#         optimal SVR model (from prosail_opt.csv, script 06b) to the full
#         Sentinel-2 reflectance raster, producing spatially continuous
#         LAI_S2_opt maps used as inputs for SM6b.
#
#         Three scenarios are processed:
#           per_site  : site-specific d_opt and Column_opt (Pareto).
#                       Output: s2lai_summer_opt_per_site_res_10_m.tif
#           common    : all-sites Pareto d_opt applied to each site.
#                       (d_opt_source = "all_sites" in prosail_opt.csv)
#                       Output: s2lai_summer_opt_common_res_10_m.tif
#           fixed_4   : d_opt = 4 m for all sites (submitted-paper reference).
#                       Output: s2lai_summer_opt_fixed_4_res_10_m.tif
#
#         Both rasters per site are written to:
#           revision/output/intermediate/sm6/{site}/
#
#         Prerequisites:
#           - revision/output/intermediate/sm5/prosail_opt.csv  (script 06b)
#           - revision/output/intermediate/PROSAIL_Models/{site}/
#             LIDFa_lai_LMA_BROWN/{Column_opt}.rds              (script 02)
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/06c_sm5_predict_lai_raster.R")
# ---

library(here)
library(terra)
library(data.table)
library(cli)

source(here::here("revision", "R", "prosail_inversion.R"))

# ── Configuration ──────────────────────────────────────────────────────────────

norm_ref        <- "DSM_keepTrees"
name_strategy   <- "LIDFa_lai_LMA_BROWN"
band_prefixes   <- c("B03", "B04", "B08")

# Scenarios to process: named list of output suffix → d_opt_source label in prosail_opt.csv
scenarios <- list(
  per_site = "per_site",
  common   = "all_sites",
  fixed_4  = "fixed_4"
)

l2a_ids <- c(
  Aigoual = "L2A_T31TEJ_A031608_20210711T104217",
  Blois   = "L2A_T31TCN_A031222_20210614T105443",
  Mormal  = "L2A_T31UER_A031222_20210614T105443"
)

# ── Prerequisites ──────────────────────────────────────────────────────────────

prosail_opt_csv <- here::here(
  "revision", "output", "intermediate", "sm5", "prosail_opt.csv"
)
if (!file.exists(prosail_opt_csv)) {
  stop("prosail_opt.csv not found — run script 06b first.\n  ", prosail_opt_csv)
}

prosail_opt <- data.table::fread(prosail_opt_csv)

# ── Output directory ────────────────────────────────────────────────────────────

sm6_dir <- here::here("revision", "output", "intermediate", "sm6")

# ── Process each scenario × site ───────────────────────────────────────────────

t0 <- proc.time()

for (scenario in names(scenarios)) {
  source_label <- scenarios[[scenario]]
  cli::cli_h1("Scenario: {scenario}  (d_opt_source = {source_label})")

  opt_scen <- prosail_opt[d_opt_source == source_label & Norm == norm_ref &
                            Site %in% names(l2a_ids)]

  if (nrow(opt_scen) == 0L) {
    cli::cli_warn("No rows for d_opt_source={source_label} in prosail_opt.csv — skipping.")
    next
  }

  for (site in opt_scen$Site) {
    cli::cli_h2("Site: {site}")

    row        <- opt_scen[Site == site]
    column_opt <- row$Column_opt
    d_opt_val  <- row$d_opt

    cli::cli_alert_info("  Column_opt : {column_opt}")
    cli::cli_alert_info("  d_opt      : {d_opt_val} m")

    # ── SVR model ──────────────────────────────────────────────────────────────
    rds_path <- here::here(
      "revision", "output", "intermediate", "PROSAIL_Models",
      site, name_strategy, paste0(column_opt, ".rds")
    )
    if (!file.exists(rds_path)) {
      cli::cli_warn("  SVR model not found — skipping:\n  {rds_path}")
      next
    }
    cli::cli_alert_info("  Loading SVR ensemble ...")
    svr <- load_svr_ensemble(rds_path)

    # ── S2 reflectance raster ─────────────────────────────────────────────────
    l2a     <- l2a_ids[[site]]
    s2_path <- here::here(
      "03_RESULTS", site, "PROSAIL_Optimization", l2a,
      "Reflectance", paste0(l2a, "_Refl")
    )
    if (!file.exists(s2_path)) {
      cli::cli_warn("  S2 raster not found — skipping:\n  {s2_path}")
      next
    }
    cli::cli_alert_info("  Loading S2 raster ...")
    s2_rast <- terra::rast(s2_path)

    # ── Deciduous-only mask ───────────────────────────────────────────────────
    mask_path <- here::here(
      "03_RESULTS", site, "LiDAR", "Heterogeneity_Masks",
      "artifacts_deciduous_only_low_vegetation_majority_90_p_res_10_m.envi"
    )
    mask_rast <- if (file.exists(mask_path)) terra::rast(mask_path) else {
      cli::cli_warn("  Mask not found — output will not be masked for {site}.")
      NULL
    }

    # ── Predict ───────────────────────────────────────────────────────────────
    cli::cli_alert_info("  Predicting LAI on full raster ...")
    lai_opt <- predict_lai_full_raster(
      svr_ensemble  = svr,
      s2_raster     = s2_rast,
      band_prefixes = band_prefixes,
      mask_rast     = mask_rast
    )

    # ── Write output ─────────────────────────────────────────────────────────
    site_out_dir <- file.path(sm6_dir, site)
    dir.create(site_out_dir, recursive = TRUE, showWarnings = FALSE)

    out_fn   <- paste0("s2lai_summer_opt_", scenario, "_res_10_m.tif")
    out_path <- file.path(site_out_dir, out_fn)
    terra::writeRaster(lai_opt, out_path, overwrite = TRUE, gdal = "COMPRESS=LZW")

    n_valid <- sum(!is.na(terra::values(lai_opt)))
    cli::cli_alert_success("  Written: {out_fn} ({n_valid} valid pixels)")
  }
}

elapsed <- round((proc.time() - t0)[["elapsed"]], 1)

cli::cli_h1("SM5 passe 3 — résumé")
cli::cli_bullets(c(
  "v" = "Scénarios : {paste(names(scenarios), collapse=', ')}",
  "v" = "Sorties   : revision/output/intermediate/sm6/{site}/s2lai_summer_opt_{scenario}_res_10_m.tif",
  "i" = "Durée     : {elapsed} s"
))
