# ---
# title:  01b_download_s2_from_cdse.R
# desc:   Orchestration — downloads Sentinel-2 L2A imagery from Copernicus
#         Data Space (CDSE) for the three study sites (Aigoual, Blois, Mormal).
#         preprocS2::get_s2_raster() handles download, reflectance extraction,
#         cloud mask, and geometry-of-acquisition angles in a single call.
#
#         Mirrors PROSAIL-Optimization/02_CODES/Sentinel2/Main_s2_01_download.R
#         refactored with here::here() paths and R/sentinel2_preprocessing.R.
#
# Prerequisite: CDSE credentials must be set in ~/.Renviron before running:
#   CDSE_ID="sh-XXXXXXXX-XXXX-XXXX-XXXX-XXXXXXXXXXXX"
#   CDSE_SECRET="XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
#   CURL_SSL_BACKEND="openssl"
#   Then restart R or call readRenviron("~/.Renviron").
#
# WARNING: do NOT read, display, or copy 01_DATA/copernicus-credentials*.yml.
#
# Run from the project root (NC_Full/):
#   source("scripts/01b_download_s2_from_cdse.R")
# ---

library("here")
library("preprocS2")
source(here::here("R", "paths.R"))
source(here::here("R", "sentinel2_preprocessing.R"))

# ── Site parameters — unchanged from Main_s2_01_download.R ────────────────────

sites <- c("Aigoual", "Blois", "Mormal")

dates <- c(
  Aigoual = "2021-07-11",
  Blois   = "2021-06-14",
  Mormal  = "2021-06-14"
)

name_vect <- "utm_init.shp"

# ── Processing loop ────────────────────────────────────────────────────────────

for (site in sites) {

  aoi_path   <- file.path(paths$raw_data, site, "Geo_Files", name_vect)
  output_dir <- file.path(paths$ext_results, site, "PROSAIL_Optimization")

  # download_s2_cdse() reads CDSE_ID and CDSE_SECRET from .Renviron at runtime
  download_s2_cdse(
    site       = site,
    date       = dates[[site]],
    aoi_path   = aoi_path,
    output_dir = output_dir
    # credentials_env defaults to c("CDSE_ID", "CDSE_SECRET")
  )
}
