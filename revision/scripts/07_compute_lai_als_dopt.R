# ---
# title:  07_compute_lai_als_dopt.R
# desc:   Orchestration — computes LAI_ALS_dopt rasters for each site under
#         three scenarios:
#
#           per_site  : cumulative PAD from canopy top to each site's own
#                       Pareto d_opt (from dopt_reference.csv).
#           common    : cumulative PAD from canopy top to the All_sites
#                       Pareto d_opt (same depth for all three sites).
#           fixed_4   : cumulative PAD from canopy top to 4 m depth
#                       (reproduces the submitted paper).
#
#         LAI_ALS_dopt is computed as sum(ladstack[[1:d_opt]]) per pixel,
#         where layer 1 = topmost 1 m slice of the canopy.
#
#         These rasters are the primary LAI input for SVR training in script 08
#         (LiDAR_LAI_Best_Site_Depth variant of the simulation grid).
#
# Prerequisite: dopt_reference.csv produced by script 06 must be present at
#   revision/output/intermediate/sm5/dopt_reference.csv
#
# Output: one GeoTIFF per site × scenario at:
#   revision/output/intermediate/lai_als_dopt/{site}/
#   LAI_ALS_dopt_{scenario}.tif
#   (scenario ∈ {per_site, common, fixed_4})
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/07_compute_lai_als_dopt.R")
# ---

library("here")
library("terra")
library("data.table")

# ── Parameters ─────────────────────────────────────────────────────────────────

sites         <- c("Aigoual", "Blois", "Mormal")
norm_select   <- "DSM_keepTrees"   # norm used in d_opt selection (06)
method_select <- "pareto"          # d_opt method (06)
fixed_depth   <- 4L                # reference depth from submitted paper

# ── Load d_opt reference ───────────────────────────────────────────────────────

dopt_csv <- here::here("revision", "output", "intermediate", "sm5",
                        "dopt_reference.csv")
if (!file.exists(dopt_csv))
  stop("dopt_reference.csv not found — run script 06 first.\n  ", dopt_csv)

dopt_ref <- data.table::fread(dopt_csv)
dopt_ref <- dopt_ref[Norm == norm_select & method_dopt == method_select]

# Common d_opt = All_sites row
d_opt_common <- dopt_ref[Site == "All_sites", d_opt]
if (length(d_opt_common) != 1L)
  stop("Expected exactly one All_sites row in dopt_reference.csv")

cat("Common d_opt (All_sites, pareto):", d_opt_common, "m\n")

# Per-site d_opt (exclude All_sites row)
d_opt_per_site <- setNames(
  dopt_ref[Site != "All_sites", d_opt],
  dopt_ref[Site != "All_sites", Site]
)
cat("Per-site d_opt:", paste(names(d_opt_per_site),
                              d_opt_per_site, sep = "=", collapse = ", "), "m\n")

# ── Output directory ───────────────────────────────────────────────────────────

out_root <- here::here("revision", "output", "intermediate", "lai_als_dopt")

# ── Helper: cumulative PAD sum to depth d_opt ──────────────────────────────────
# ladstack layer i = PAD in the i-th metre from canopy top (depth i).
# LAI_ALS_dopt = sum of layers 1 through d_opt.

compute_lai_als_dopt_rast <- function(ladstack_path, d_opt) {
  lad <- terra::rast(ladstack_path)
  if (d_opt > terra::nlyr(lad))
    stop("d_opt (", d_opt, ") exceeds number of layers in ladstack (",
         terra::nlyr(lad), ")")
  sum(lad[[seq_len(d_opt)]], na.rm = FALSE)
}

# ── Processing loop ────────────────────────────────────────────────────────────

for (site in sites) {

  cat("\n── Processing site:", site, "──\n")

  ladstack_path <- here::here(
    "PROSAIL-Optimization", "01_DATA", site, "LiDAR",
    "PAD_Profiles_Classic", "ladstack.tif"
  )

  out_dir <- file.path(out_root, site)
  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  scenarios <- list(
    per_site = d_opt_per_site[[site]],
    common   = d_opt_common,
    fixed_4  = fixed_depth
  )

  for (scenario in names(scenarios)) {

    d_opt <- scenarios[[scenario]]
    out_tif <- file.path(out_dir, paste0("LAI_ALS_dopt_", scenario, ".tif"))

    cat("  scenario:", scenario, "| d_opt =", d_opt, "m ->", basename(out_tif), "\n")

    rast_dopt <- compute_lai_als_dopt_rast(ladstack_path, d_opt)
    terra::writeRaster(rast_dopt, out_tif, overwrite = TRUE, datatype = "FLT4S")
  }
}

cat("\n── 07 done — LAI_ALS_dopt rasters written to", out_root, "──\n")
