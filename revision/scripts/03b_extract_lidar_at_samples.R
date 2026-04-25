# ---
# title:  03b_extract_lidar_at_samples.R
# desc:   Orchestration — extracts LiDAR PAD profile values at the pixel
#         positions sampled by 03a_sample_s2_pixels.R, for all sites × PAD
#         normalisation types × depth layers.
#
#         For each (site × norm_type × depth) triplet, reads the GPKG of
#         sample positions produced by 03a and extracts PAD values from the
#         corresponding single-depth raster file, writing one CSV per triplet.
#         Total at execution: 3 sites × 4 norm types × 38 depths = 456 CSVs.
#
#         Mirrors PROSAIL-Optimization/02_CODES/Sentinel2/Main_s2_02B_sample_lidar.R
#         refactored with here::here() paths and revision/R/sentinel2_sampling.R.
#
# Prerequisite: 03a_sample_s2_pixels.R must have been run first (produces
#   the Sampling_stratified_uniform_nbSamples_5000.GPKG files).
#
# NOTE:   PAD rasters are read from 01_DATA/ (read-only).
#         CSVs are written to 03_RESULTS/{site}/PROSAIL_Optimization/sampling/
#         to remain compatible with the d_opt analysis pipeline (SM5).
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/03b_extract_lidar_at_samples.R")
# ---

library("here")
library("terra")
library("readr")

source(here::here("revision", "R", "paths.R"))
source(here::here("revision", "R", "sentinel2_sampling.R"))

# ── Global parameters — unchanged from Main_s2_02B_sample_lidar.R ─────────────

sites            <- c("Aigoual", "Blois", "Mormal")
nbSamplesS2Refl  <- 5000
Samplingmethod   <- "stratified_uniform"
h_min_values     <- c(10L, 15L, 20L)   # must match 03a

# PAD normalisation types and their sub-directories under
# 01_DATA/{site}/LiDAR/testPADs/
norm_types <- c(
  DSM          = "PAD_Profiles_DSM",
  DTM          = "PAD_Profiles_DTM",
  DSM_keepTrees = "PAD_Profiles_DSM_keepTrees",
  DTM_keepTrees = "PAD_Profiles_DTM_keepTrees"
)

# ── Processing loop ────────────────────────────────────────────────────────────

for (site in sites) {

  cat("Processing site:", site, "\n")

  for (h_min in h_min_values) {

  # Sample locations GPKG produced by 03a
  filename_sampling <- file.path(
    paths$ext_results, site, "PROSAIL_Optimization", "sampling",
    paste0("Sampling_", Samplingmethod,
           "_hmin", h_min,
           "_nbSamples_", nbSamplesS2Refl, ".GPKG")
  )

  out_sampling <- file.path(paths$ext_results, site, "PROSAIL_Optimization",
                             "sampling")

  for (norm_name in names(norm_types)) {

    pad_dir <- file.path(paths$prosail_lidar, site, "LiDAR",
                          "testPADs", norm_types[[norm_name]])

    # List all single-depth PAD rasters (pattern: PAD_X.5_40.tif)
    pad_files <- list.files(pad_dir,
                             pattern   = "PAD_.*\\.tif$",
                             full.names = TRUE)

    for (pad_file in pad_files) {

      # Extract depth from file name (format: PAD_X.5_40.tif)
      x_value <- floor(
        as.numeric(sub(".*_(\\d+\\.5)_40.tif", "\\1", basename(pad_file)))
      )
      depth <- 40 - x_value   # effective canopy depth (m from top)

      output_csv_path <- file.path(
        out_sampling,
        paste0("PAD_", norm_name, "_Depth_", depth,
               "_Samples_", Samplingmethod,
               "_hmin", h_min,
               "_nbSamples_", nbSamplesS2Refl, ".csv")
      )

      extract_pad_at_samples(
        sample_gpkg_path = filename_sampling,
        pad_rast_path    = pad_file,
        output_csv_path  = output_csv_path
      )
    }
  }

  }  # end h_min loop
}

cat("Done. CSVs written to 03_RESULTS/{site}/PROSAIL_Optimization/sampling/\n")
