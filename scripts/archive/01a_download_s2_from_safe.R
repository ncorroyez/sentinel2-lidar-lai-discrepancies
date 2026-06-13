# ---
# title:  01a_download_s2_from_safe.R
# desc:   Orchestration — extracts Sentinel-2 L2A reflectances and cloud masks
#         from SAFE archives already present on disk for the three study sites
#         (Aigoual, Blois, Mormal).
#
#         Mirrors PROSAIL-Optimization/02_CODES/Sentinel2/Main_s2_00_download.R
#         refactored with here::here() paths and R/sentinel2_preprocessing.R.
#
# Prerequisite: SAFE archives must be present at
#   01_DATA/{site}/Sentinel-2/{date}/{safe_name}.SAFE
#
# Run from the project root (NC_Full/):
#   source("scripts/01a_download_s2_from_safe.R")
# ---

library("here")
library("preprocS2")
source(here::here("R", "paths.R"))
source(here::here("R", "sentinel2_preprocessing.R"))

# ── Site parameters — unchanged from Main_s2_00_download.R ────────────────────

sites <- c("Aigoual", "Blois", "Mormal")

dates <- c(
  Aigoual = "2021-07-11",
  Blois   = "2021-06-14",
  Mormal  = "2021-06-14"
)

safe_names <- c(
  Aigoual = "S2A_MSIL2A_20210711T104031_N0301_R008_T31TEJ_20210711T134956.SAFE",
  Blois   = "S2A_MSIL2A_20210614T105031_N0300_R051_T31TCN_20210614T140120.SAFE",
  Mormal  = "S2A_MSIL2A_20210614T105031_N0300_R051_T31UER_20210614T140120.SAFE"
)

resolution <- 10
s2source   <- "SAFE"

# ── Processing loop ────────────────────────────────────────────────────────────

for (site in sites) {

  safe_path  <- file.path(paths$raw_data, site, "Sentinel-2", dates[[site]],
                           safe_names[[site]])
  aoi_path   <- file.path(paths$raw_data, site, "Geo_Files", "utm_init.shp")
  out_root   <- file.path(paths$ext_results, site, "PROSAIL_Optimization")

  # Step 1 — extract reflectances and metadata from SAFE
  s2obj <- load_s2_from_safe(
    safe_path   = safe_path,
    path_vector = aoi_path,
    resolution  = resolution,
    s2source    = s2source
  )

  # Results land in a sub-directory named after the GRANULE identifier
  granule_dir <- file.path(out_root, basename(s2obj$S2_Bands$GRANULE))
  dir.create(granule_dir, recursive = TRUE, showWarnings = FALSE)

  # Step 2 — write cloud mask (SCL-derived binary mask via preprocS2)
  cloud_path <- file.path(granule_dir, "CloudMask")
  dir.create(cloud_path, showWarnings = FALSE, recursive = TRUE)
  generate_cloud_mask(
    S2_Stack   = s2obj$S2_Stack,
    Cloud_path = cloud_path,
    S2source   = s2source,
    saveRaw    = TRUE
  )

  # Step 3 — write reflectances as ENVI Int16 BIL
  refl_dir  <- file.path(granule_dir, "Reflectance")
  dir.create(refl_dir, showWarnings = FALSE, recursive = TRUE)
  refl_path <- file.path(refl_dir,
                          paste0(basename(s2obj$S2_Bands$GRANULE), "_Refl"))
  write_reflectance(
    S2_Stack  = s2obj$S2_Stack,
    S2_Bands  = s2obj$S2_Bands,
    refl_path = refl_path
  )
}
