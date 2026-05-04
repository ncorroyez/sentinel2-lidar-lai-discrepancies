# ---
# title:  bundle_inputs.R
# desc:   Copies all input data required to run the full pipeline (phase 1→5)
#         from the SMB share (or any configured data_root) to revision/data/,
#         so the pipeline can be run locally with profile: local in config.yml.
#
#         Run once from the project root (NC_Full/) while the SMB is mounted:
#           source("revision/bundle_inputs.R")
#
#         The resulting revision/data/ tree (~300–500 Mo) can be zipped and
#         shared with anyone who has access to the code on GitHub.
#
# Recipient instructions:
#   1. Clone the code: git clone <github_url>
#   2. Extract the data bundle into revision/data/
#   3. Copy config.yml.template to config.yml
#   4. Set profile: local in config.yml and adjust data_root / output_root
#   5. source("revision/scripts/00_run_all.R")
# ---

library(here)
source(here::here("revision", "R", "paths.R"))

# ── Destination roots ─────────────────────────────────────────────────────────

dst_raw     <- here::here("revision", "data", "raw")
dst_results <- here::here("revision", "data", "results")
dst_lidar   <- here::here("revision", "data", "prosail_lidar")

sites <- c("Aigoual", "Blois", "Mormal")

l2a_ids <- c(
  Aigoual = "L2A_T31TEJ_A031608_20210711T104217",
  Blois   = "L2A_T31TCN_A031222_20210614T105443",
  Mormal  = "L2A_T31UER_A031222_20210614T105443"
)

# ── Helpers ───────────────────────────────────────────────────────────────────

#' Copy a directory tree, preserving relative structure
cp_dir <- function(src, dst) {
  if (!dir.exists(src)) {
    warning("Source directory not found, skipping: ", src)
    return(invisible(FALSE))
  }
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(src, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  for (f in files) {
    rel  <- sub(paste0("^", normalizePath(src), .Platform$file.sep), "", normalizePath(f))
    dest <- file.path(dst, rel)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    file.copy(f, dest, overwrite = TRUE)
  }
  invisible(TRUE)
}

#' Copy a single file, creating destination directory if needed
cp_file <- function(src, dst_dir) {
  if (!file.exists(src)) {
    warning("Source file not found, skipping: ", src)
    return(invisible(FALSE))
  }
  dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(src, file.path(dst_dir, basename(src)), overwrite = TRUE)
  invisible(TRUE)
}

# ── Per-site copy ─────────────────────────────────────────────────────────────

for (site in sites) {

  message("\n── Bundling site: ", site, " ──")

  l2a <- l2a_ids[[site]]

  # ── raw_data ─────────────────────────────────────────────────────────────────

  # Geo_Files — AOI bbox, UTM shapefile, field points
  cp_dir(
    file.path(paths$raw_data, site, "Geo_Files"),
    file.path(dst_raw, site, "Geo_Files")
  )

  # LiDAR metric rasters (mean, lskew) — used by 03a for sampling filter
  for (tif in c("mean_res_10_m.tif", "lskew_res_10_m.tif")) {
    cp_file(
      file.path(paths$raw_data, site, "LiDAR", tif),
      file.path(dst_raw, site, "LiDAR")
    )
  }

  # ── ext_results ──────────────────────────────────────────────────────────────

  # S2 geometry angles (used to build geom_s2 for SVR training)
  cp_dir(
    file.path(paths$ext_results, site, "PROSAIL_Optimization", "geomAcq_S2"),
    file.path(dst_results, site, "PROSAIL_Optimization", "geomAcq_S2")
  )

  # S2 L2A reflectances (used for sampling + raster prediction)
  cp_dir(
    file.path(paths$ext_results, site, "PROSAIL_Optimization", l2a),
    file.path(dst_results, site, "PROSAIL_Optimization", l2a)
  )

  # Metrics/Deciduous_Only — only the specific files needed by the pipeline
  # (not the full directory which contains ~8 GB of legacy outputs)
  dec_only_src <- file.path(paths$ext_results, site, "Metrics", "Deciduous_Only")
  dec_only_dst <- file.path(dst_results, site, "Metrics", "Deciduous_Only")
  for (tif in c("lidarlai_res_10_m.tif", "max_res_10_m.tif",
                "ladstack_classic.tif", "fCover_res_10_m.tif")) {
    cp_file(file.path(dec_only_src, tif), dec_only_dst)
  }
  # PAD profiles (DSM variant) — used by sm6_load_rasters for PAD figures
  cp_dir(
    file.path(dec_only_src, "PAD_Profiles_dsm_keepTrees"),
    file.path(dec_only_dst, "PAD_Profiles_dsm_keepTrees")
  )

  # LiDAR/dsm — 1 m DSM for heterogeneity (VRT + underlying tiles)
  cp_dir(
    file.path(paths$ext_results, site, "LiDAR", "dsm", "res_1_m"),
    file.path(dst_results, site, "LiDAR", "dsm", "res_1_m")
  )

  # LiDAR/chm — 1 m CHM
  cp_file(
    file.path(paths$ext_results, site, "LiDAR", "chm", "res_1_m", "chm.tif"),
    file.path(dst_results, site, "LiDAR", "chm", "res_1_m")
  )

  # Heterogeneity mask (deciduous only, 90th percentile)
  cp_dir(
    file.path(paths$ext_results, site, "LiDAR", "Heterogeneity_Masks"),
    file.path(dst_results, site, "LiDAR", "Heterogeneity_Masks")
  )

  # ── prosail_lidar ─────────────────────────────────────────────────────────────

  lidar_src <- file.path(paths$prosail_lidar, site, "LiDAR")

  # Ladstack (full LAI profile sum) — Phase A pool in 09_train_prosail_full
  cp_dir(
    file.path(lidar_src, "PAD_Profiles_Classic"),
    file.path(dst_lidar, site, "LiDAR", "PAD_Profiles_Classic")
  )

  # DSM-based PAD profiles — used by PAD figures + sweep
  cp_dir(
    file.path(lidar_src, "PAD_Profiles_dsm_keepTrees"),
    file.path(dst_lidar, site, "LiDAR", "PAD_Profiles_dsm_keepTrees")
  )

  # Max canopy height raster — coherence mask
  cp_file(
    file.path(lidar_src, "max_res_10_m.tif"),
    file.path(dst_lidar, site, "LiDAR")
  )

  message("  Done: ", site)
}

# ── Size summary ──────────────────────────────────────────────────────────────

dirs <- c(dst_raw, dst_results, dst_lidar)
sizes <- vapply(dirs, function(d) {
  fs <- list.files(d, recursive = TRUE, full.names = TRUE)
  sum(file.size(fs), na.rm = TRUE) / 1024^2  # MB
}, numeric(1))

message("\n── Bundle complete ──")
message(sprintf("  raw/           : %6.0f MB", sizes[[1]]))
message(sprintf("  results/       : %6.0f MB", sizes[[2]]))
message(sprintf("  prosail_lidar/ : %6.0f MB", sizes[[3]]))
message(sprintf("  TOTAL          : %6.0f MB", sum(sizes)))
message("\nShare revision/data/ as a zip — recipient extracts, sets profile: local in config.yml.")
