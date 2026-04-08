# ================================================================
# rstactheia
# ================================================================

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

library(rhdf5)
library(data.table)
library(sf)
library(stringr)
library(lubridate)
library(dplyr)
library(terra)
library(rstac)
library(rstactheia)
library(lidR)
library(mapview)

# ---------------------
# Parameters
# ---------------------
# data_dir <- "../../01_DATA"
data_dir <- "/media/corroyez/MyPassport/01_DATA"
results_dir <- "../../03_RESULTS"
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- c("Aigoual")

for (site in sites){
  my_roi <- st_read(file.path(data_dir, site,
                              "Geo_Files/utm_init.shp"))
  roi_wgs84 <- st_transform(my_roi, crs = 4326)
  
  stac_items <- stac('https://api.stac.teledetection.fr') |>
    stac_search(
      collections = "geogedi",
      datetime = "2021-05-01T00:00:00Z/2022-12-31T23:59:59Z",
      bbox = st_bbox(roi_wgs84)
    ) |> 
    # ext_filter(`s2:mgrs_tile`=="30TXP") |>
    get_request() |>
    items_fetch()
  
  # create_api_key(save=T)
  
  stac_items <- stac_items |>
    items_sign_theia()
  # url_check <- stac_items$features[[1]]$assets$copc$href
  stac_items |> assets_download(output_dir = file.path(data_dir,
                                                       site,
                                                       "GEDI"))
}

las <- readLAS("//media/corroyez/MyPassport/01_DATA/Mormal/GEDI/sm1-gdc-geogedi/geogedi/GeoGEDI02_A_20220918T223850_O21340/GeoGEDI02_A_20220918T223850_O21340.copc.laz")
las_boundary <- st_bbox(las) %>% st_as_sfc()
mapview(las_boundary, 
        map.types = "Esri.WorldImagery", 
        alpha.regions = 0.2, 
        col.regions = "red")
