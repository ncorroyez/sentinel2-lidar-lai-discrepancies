# ---
# title:  07_compute_lai_als_dopt.R
# desc:   Orchestration — computes LAI_ALS_dopt rasters for each site under
#         two scenarios:
#
#           per_site : cumulative PAD from canopy top to each site's own
#                      Pareto d_opt (from dopt_reference.csv).
#           common   : cumulative PAD from canopy top to the combined-sites
#                      Pareto d_opt (same depth for all three sites).
#
#         Both scenarios apply the k/theta rescaling factor derived from
#         k_select and theta_select (must match settings used in 06).
#         LAI_ALS_dopt = scale_factor × sum(ladstack[[1:d_opt]]) per pixel,
#         where scale_factor = (k_ref / k_select) × cos(theta_select).
#
#         These rasters are the primary LAI input for SVR training in script 08.
#
# Prerequisite: dopt_reference.csv produced by script 06 must be present at
#   output/intermediate/sm5/dopt_reference.csv
#
# Output: one GeoTIFF per site × scenario at:
#   output/intermediate/lai_als_dopt/{site}/
#   LAI_ALS_dopt_{scenario}.tif
#   (scenario ∈ {per_site, common})
#
# Run from the project root (NC_Full/):
#   source("scripts/07_compute_lai_als_dopt.R")
# ---

library("here")
library("terra")
library("data.table")

source(here::here("R", "paths.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites         <- c("Aigoual", "Blois", "Mormal")
norm_select   <- "DSM_keepTrees"   # norm used in d_opt selection (06)
method_select <- "pareto"          # d_opt method (06)

# Must match the settings used when running 06.
k_ref        <- 0.5
k_select     <- 0.65
theta_select <- 0      # degrees from nadir (scan angle correction handled separately)

# "pool"    — combined d_opt from pooled-pixel row ("Sites combined")
# "average" — combined d_opt from site-averaged metrics row ("Sites averaged")
# Must match combined_mode used in 06.
combined_mode <- "average"

# TRUE  — always recompute and overwrite existing .tif
# FALSE — skip sites/scenarios whose output .tif already exists (resume mode)
overwrite <- TRUE

# ── Derived ────────────────────────────────────────────────────────────────────

combined_label <- if (combined_mode == "average") "Sites averaged" else "Sites combined"
scale_factor   <- (k_ref / k_select) * cos(theta_select * pi / 180)

cat(sprintf("k_select=%.2f  theta_select=%d°  scale_factor=%.4f\n",
            k_select, theta_select, scale_factor))

# ── Load d_opt reference ───────────────────────────────────────────────────────

dopt_csv <- file.path(paths$output, "intermediate", "sm5",
                        "dopt_reference.csv")
if (!file.exists(dopt_csv))
  stop("dopt_reference.csv not found — run script 06 first.\n  ", dopt_csv)

dopt_ref <- data.table::fread(dopt_csv)
dopt_ref <- dopt_ref[Norm == norm_select & method_dopt == method_select]

# Common d_opt: combined-sites row (Sites averaged or Sites combined)
d_opt_common <- dopt_ref[Site == combined_label, d_opt]
if (length(d_opt_common) != 1L)
  stop("Expected exactly one '", combined_label, "' row in dopt_reference.csv",
       " (Norm=", norm_select, ", method=", method_select, ").\n",
       "  Check that combined_mode in 06 and 07 match.")

cat("Common d_opt (", combined_label, ", pareto):", d_opt_common, "m\n")

# Per-site d_opt (individual sites only)
d_opt_per_site <- setNames(
  dopt_ref[Site %in% sites, d_opt],
  dopt_ref[Site %in% sites, Site]
)
cat("Per-site d_opt:", paste(names(d_opt_per_site),
                              d_opt_per_site, sep = "=", collapse = ", "), "m\n")

# ── Output directory ───────────────────────────────────────────────────────────

out_root <- file.path(paths$output, "intermediate", "lai_als_dopt")

# ── Helper: pre-computed PAD raster at depth d_opt, rescaled for k and theta ─
# The per-pixel canopy-top integration is done offline and stored as
#   testPADs/PAD_Profiles_DSM_keepTrees/PAD_X_40.tif  where X = 40 - d_opt + 0.5.
# Each pixel holds the PAR integrated over the top d_opt metres of canopy
# at k_ref. LAI_ALS_dopt(k, θ) = scale × PAD_raster.

pad_filename_dopt <- function(d_opt, canopy_max_m = 40) {
  sprintf("PAD_%.1f_%d.tif", canopy_max_m - d_opt + 0.5, canopy_max_m)
}

compute_lai_als_dopt_rast <- function(pad_dir, d_opt, scale = 1,
                                       canopy_max_m = 40) {
  pad_file <- file.path(pad_dir, pad_filename_dopt(d_opt, canopy_max_m))
  if (!file.exists(pad_file))
    stop("PAD raster not found for d_opt = ", d_opt, " :\n  ", pad_file)
  terra::rast(pad_file) * scale
}

# ── Processing loop ────────────────────────────────────────────────────────────

for (site in sites) {

  cat("\n── Processing site:", site, "──\n")

  pad_dir <- file.path(
    paths$prosail_lidar, site, "LiDAR",
    "testPADs", "PAD_Profiles_DSM_keepTrees"
  )

  out_dir <- file.path(out_root, site)
  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  scenarios <- list(
    per_site = d_opt_per_site[[site]],
    common   = d_opt_common
  )

  for (scenario in names(scenarios)) {

    d_opt   <- scenarios[[scenario]]
    out_tif <- file.path(out_dir, paste0("LAI_ALS_dopt_", scenario, ".tif"))

    cat("  scenario:", scenario, "| d_opt =", d_opt, "m",
        "| scale =", round(scale_factor, 4),
        "->", basename(out_tif), "\n")

    if (file.exists(out_tif) && !overwrite) {
      cat("    [SKIP] tif already exists\n")
      next
    }

    rast_dopt <- compute_lai_als_dopt_rast(pad_dir, d_opt,
                                            scale = scale_factor)
    terra::writeRaster(rast_dopt, out_tif, overwrite = TRUE, datatype = "FLT4S")
  }
}

cat("\n── 07 done — LAI_ALS_dopt rasters written to", out_root, "──\n")
