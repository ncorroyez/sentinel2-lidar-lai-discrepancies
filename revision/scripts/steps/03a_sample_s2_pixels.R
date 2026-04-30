# ---
# title:  03a_sample_s2_pixels.R
# desc:   Orchestration — stratified-uniform sampling of Sentinel-2 pixels for
#         the three study sites (Aigoual, Blois, Mormal). Produces:
#           - Sampling_{method}_nbSamples_5000.GPKG  (sample positions)
#           - S2_reflectance_{method}_nbSamples_5000.csv  (S2 reflectances)
#         in 03_RESULTS/{site}/PROSAIL_Optimization/sampling/.
#
#         Mirrors PROSAIL-Optimization/02_CODES/Sentinel2/Main_s2_02A_extract.R
#         refactored with here::here() paths and revision/R/sentinel2_sampling.R.
#
# NOTE:   The AOI sampling preparation file (aoi_sampling_prep.gpkg) is written
#         to 03_RESULTS/{site}/PROSAIL_Optimization/sampling/ instead of
#         01_DATA/{site}/Geo_Files/ (original behaviour, violation of CLAUDE.md).
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/03a_sample_s2_pixels.R")
# ---

library("here")
library("terra")
library("sgsR")
library("dplyr")
library("readr")

source(here::here("revision", "R", "paths.R"))
source(here::here("revision", "R", "sentinel2_sampling.R"))

# ── Global parameters — unchanged from Main_s2_02A_extract.R ─────────────────

sites <- c("Aigoual", "Blois", "Mormal")

l2a_ids <- c(
  Aigoual = "L2A_T31TEJ_A031608_20210711T104217",
  Blois   = "L2A_T31TCN_A031222_20210614T105443",
  Mormal  = "L2A_T31UER_A031222_20210614T105443"
)

name_vect        <- "utm_init.shp"
nbSamples        <- 5000
Samplingmethods  <- c("stratified_uniform")
h_min_values     <- c(10L, 15L, 20L)   # dmin values (m); paper section 2.5.4

# ── Processing loop ────────────────────────────────────────────────────────────

for (site in sites) {

  l2a <- l2a_ids[[site]]

  # Input paths
  aoi_path          <- file.path(paths$raw_data, site, "Geo_Files", name_vect)
  S2_path           <- file.path(paths$ext_results, site, "PROSAIL_Optimization",
                                  l2a, "Reflectance", paste0(l2a, "_Refl"))
  lidar_dir         <- file.path(paths$prosail_lidar, site, "LiDAR")
  lidar_lai_path    <- file.path(lidar_dir, "PAD_Profiles_Classic", "ladstack.tif")
  mean_path         <- file.path(paths$raw_data, site, "LiDAR", "mean_res_10_m.tif")
  max_path          <- file.path(lidar_dir, "max_res_10_m.tif")
  lskew_path        <- file.path(paths$raw_data, site, "LiDAR", "lskew_res_10_m.tif")
  field_points_path <- file.path(paths$raw_data, site, "Geo_Files",
                                  "data_utm31n.geojson")

  # Output sampling directory
  out_sampling <- file.path(paths$ext_results, site, "PROSAIL_Optimization",
                             "sampling")
  dir.create(out_sampling, showWarnings = FALSE, recursive = TRUE)

  # Step 1 — clip AOI to valid S2 + LiDAR pixels
  # NOTE: output redirected to 03_RESULTS/ (original wrote to 01_DATA/).
  aoi_sampling_path <- file.path(out_sampling, "aoi_sampling_prep.gpkg")
  aoi_sampling_path <- align_and_remove_na_for_aoi(
    aoi_path       = aoi_path,
    S2_path        = S2_path,
    lidar_lai_path = lidar_lai_path,
    output_path    = aoi_sampling_path,
    save           = TRUE
  )

  # Step 2 — sample and extract S2 reflectances for each h_min
  for (Samplingmethod in Samplingmethods) {
    for (h_min in h_min_values) {

      set.seed(42)

      S2Refl <- get_s2_samples(
        aoi_path          = aoi_sampling_path,
        S2_path           = S2_path,
        mean_path         = mean_path,
        lskew_path        = lskew_path,
        lidar_lai_path    = lidar_lai_path,
        max_path          = max_path,
        field_points_path = field_points_path,
        nbSamples         = nbSamples,
        site              = site,
        method            = Samplingmethod,
        h_min             = h_min
      )

      # Save sample locations (GeoPackage) with sample_id
      filename_sampling <- file.path(out_sampling,
                                      paste0("Sampling_", Samplingmethod,
                                             "_hmin", h_min,
                                             "_nbSamples_", nbSamples, ".GPKG"))
      sample_locations <- S2Refl$sample_location
      sample_locations$sample_id <- seq_len(nrow(sample_locations))
      terra::writeVector(x        = sample_locations,
                          filename = filename_sampling,
                          filetype = "GPKG",
                          overwrite = TRUE)

      # Save S2 reflectances (TSV CSV)
      filename_refl <- file.path(out_sampling,
                                  paste0("S2_reflectance_", Samplingmethod,
                                         "_hmin", h_min,
                                         "_nbSamples_", nbSamples, ".csv"))
      s2_refl <- S2Refl$S2_refl
      readr::write_delim(x = s2_refl, file = filename_refl, delim = "\t")
    }
  }
}
