# ---
# title: "debug_preprocS2.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2025-01-13"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}
library(preprocS2)
library(sf)

# 1- define input & output directories
output_dir <- './RESULTS'
datetime <- as.Date('2021-06-14')

# 2- define area of interest
bbox <- sf::st_bbox(obj = c('xmin' = -54.1, 'ymin' = 5.0,
                            'xmax' = -54.0, 'ymax' = 5.1))
aoi_path <- 'utm_init.shp'
sf::st_write(obj = bbox_to_poly(x = bbox), dsn = aoi_path, overwrite = T)

# get tiling grid kml from https://sentiwiki.copernicus.eu/web/s2-products
path_S2tilinggrid <- "S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.kml"

# 3- get S2 acquisition and geometry of acquisition for study area
# https://shapps.dataspace.copernicus.eu/dashboard/#/account/settings
id <- "sh-3699814c-729b-4d29-abde-91163718cf5a"
secret <- "r8Z7znbv6A3uQHiFTFgnlg0mfvg4Y0Vy"
authentication <- list('id' = id, 'pwd' = secret)

get_s2_raster(
  aoi_path = aoi_path,
  # bbox = bbox,
  datetime = datetime, output_dir = output_dir, 
  path_S2tilinggrid = path_S2tilinggrid, 
  overwrite = T, geomAcq = T, authentication = authentication
)