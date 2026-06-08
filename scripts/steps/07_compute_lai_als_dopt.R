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
#   revision/output/intermediate/sm5/dopt_reference.csv
#
# Output: one GeoTIFF per site × scenario at:
#   revision/output/intermediate/lai_als_dopt/{site}/
#   LAI_ALS_dopt_{scenario}.tif
#   (scenario ∈ {per_site, common})
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/07_compute_lai_als_dopt.R")
# ---

library("here")
library("terra")
library("data.table")

source(here::here("revision", "R", "paths.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites         <- c("Aigoual", "Blois", "Mormal")
norm_select   <- "DSM_keepTrees"   # norm used in d_opt selection (06)
method_select <- "pareto"          # d_opt method (06)

# Must match the settings used when running 06.
k_ref        <- 0.5
k_select     <- 0.6
theta_select <- 0      # degrees from nadir (scan angle correction handled separately)

# "pool"    — combined d_opt from pooled-pixel row ("Sites combined")
# "average" — combined d_opt from site-averaged metrics row ("Sites averaged")
# Must match combined_mode used in 06.
combined_mode <- "average"

# ── Derived ────────────────────────────────────────────────────────────────────

combined_label <- if (combined_mode == "average") "Sites averaged" else "Sites combined"
scale_factor   <- (k_ref / k_select) * cos(theta_select * pi / 180)

cat(sprintf("k_select=%.1f  theta_select=%d°  scale_factor=%.4f\n",
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

# ── Helper: cumulative PAD sum to depth d_opt, rescaled for k and theta ────────
# ladstack layer i = PAD in the i-th metre from canopy top (depth i) at k_ref.
# LAI_ALS_dopt(k, θ) = scale_factor × sum(ladstack[[1:d_opt]])

compute_lai_als_dopt_rast <- function(ladstack_path, d_opt, scale = 1) {
  lad <- terra::rast(ladstack_path)
  if (d_opt > terra::nlyr(lad))
    stop("d_opt (", d_opt, ") exceeds number of layers in ladstack (",
         terra::nlyr(lad), ")")
  sum(lad[[seq_len(d_opt)]], na.rm = FALSE) * scale
}

# ── Processing loop ────────────────────────────────────────────────────────────

for (site in sites) {

  cat("\n── Processing site:", site, "──\n")

  ladstack_path <- file.path(
    paths$prosail_lidar, site, "LiDAR",
    "PAD_Profiles_Classic", "ladstack.tif"
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

    if (file.exists(out_tif)) {
      cat("    [SKIP] tif already exists\n")
      next
    }

    rast_dopt <- compute_lai_als_dopt_rast(ladstack_path, d_opt,
                                            scale = scale_factor)
    terra::writeRaster(rast_dopt, out_tif, overwrite = TRUE, datatype = "FLT4S")
  }
}

cat("\n── 07 done — LAI_ALS_dopt rasters written to", out_root, "──\n")
