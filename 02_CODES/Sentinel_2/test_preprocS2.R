# ---
# title: "test_preprocS2.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2025-01-11"
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
library(ggplot2)

# Import useful functions
source("../libraries/functions_create_masks.R")
source("../libraries/functions_sentinel_2.R")
source("../libraries/functions_plots.R")

# Pre-processing Parameters
# data_dir <- '../../01_DATA'
# data_dir <- '/media/corroyez/My Passport/01_DATA'
# data_dir <- '../../01_DATA'
# sites <- c("Aigoual", "Blois", "Mormal")
# sites <- "Aigoual" # Mormal Blois Aigoual
sites <- c("Aigoual")
results_dir <- "../../03_RESULTS"

for (site in sites){
  # 1- define input & output directories
  if (site == "Aigoual") {
    datetime <- as.Date('2021-07-11')
  } else if (site == "Blois") {
    datetime <- as.Date('2021-06-14')
  } else if (site == "Mormal") {
    datetime <- as.Date('2021-06-14')
  } else {
    stop("Unknown site. Please provide a valid site name.")
  }
  output_dir <- file.path(results_dir, site, 'preprocS2')
  
  # 2- define area of interest
  # bbox <- sf::st_bbox(obj = c('xmin' = -54.1, 'ymin' = 5.0,
  #                             'xmax' = -54.0, 'ymax' = 5.1))
  shp_path <- file.path(output_dir, 'utm_init.shp')
  shp_read <- st_read(shp_path)
  bbox <- st_bbox(shp_read)
  bbox_poly <- st_as_sfc(bbox)
  bbox_poly_deg <- st_transform(bbox_poly, crs = 4326)
  bbox_deg <- st_bbox(bbox_poly_deg)
  
  plot <- ggplot() +
    geom_sf(data = shp_read, fill = "lightblue", color = "darkblue") +
    geom_sf(data = bbox_poly_deg, fill = NA, color = "red", linetype = "dashed") +
    theme_minimal() +
    ggtitle("Shapefile and Bounding Box")
  plot(plot)
  
  aoi_path <- file.path(output_dir, 'S2_aoi.GPKG')
  sf::st_write(obj = bbox_to_poly(x = bbox_deg), dsn = aoi_path, driver = 'GPKG', 
               append = F, overwrite = T)
  
  vectdata <- sf::read_sf(file.path(output_dir, 'S2_aoi.GPKG'))[1, ]
  
  vect_bbox <- vectdata |>
    sf::st_transform(4326) |>
    sf::st_bbox() 
  
  # get tiling grid kml from https://sentiwiki.copernicus.eu/web/s2-products
  path_S2tilinggrid <- file.path("../../Sentinel-2_tiling_grid.kml")
  
  # 3- get S2 acquisition and geometry of acquisition for study area
  # https://shapps.dataspace.copernicus.eu/dashboard/#/account/settings
  id <- "sh-3699814c-729b-4d29-abde-91163718cf5a"
  secret <- "r8Z7znbv6A3uQHiFTFgnlg0mfvg4Y0Vy"
  authentication <- list('id' = id, 'pwd' = secret)
  
  get_s2_raster(
                # aoi_path = aoi_path,
                bbox = vect_bbox,
                datetime = datetime, output_dir = output_dir, 
                path_S2tilinggrid = path_S2tilinggrid, 
                overwrite = T, geomAcq = T, authentication = authentication
                )
}




# 1- define input & output directories
output_dir <- './RESULTS_bbox'
dir.create(path = output_dir, showWarnings = F, recursive = T)
datetime <- as.Date('2024-10-05')

# 2- define area of interest
bbox <- sf::st_bbox(obj = c('xmin' = -54.1, 'ymin' = 5.0, 
                            'xmax' = -54.0, 'ymax' = 5.1))

# get tiling grid kml from https://sentiwiki.copernicus.eu/web/s2-products
path_S2tilinggrid <- '../../Sentinel-2_tiling_grid.kml'

# 3- get S2 acquisition and geometry of acquisition for study area
# https://shapps.dataspace.copernicus.eu/dashboard/#/account/settings
id <- "sh-3699814c-729b-4d29-abde-91163718cf5a"
secret <- "r8Z7znbv6A3uQHiFTFgnlg0mfvg4Y0Vy"
authentication <- list('id' = id, 'pwd' = secret)

get_s2_raster(bbox = bbox, datetime = datetime, output_dir = output_dir, 
              path_S2tilinggrid = path_S2tilinggrid, siteName = 'FrenchGuiana',
              overwrite = F, geomAcq = T, authentication = authentication)




# 1- define input & output directories
output_dir <- './RESULTS_vector'
dir.create(path = output_dir, showWarnings = F, recursive = T)
datetime <- as.Date('2024-10-05')

# 2- define area of interest
bbox <- sf::st_bbox(obj = c('xmin' = -54.1, 'ymin' = 5.0, 
                            'xmax' = -54.0, 'ymax' = 5.1))
aoi_path <- file.path(output_dir, 'S2_FrenchGuiana_aoi.GPKG')
sf::st_write(obj = bbox_to_poly(x = bbox), dsn = aoi_path, driver = 'GPKG',append = F, overwrite = T)
# get tiling grid kml from https://sentiwiki.copernicus.eu/web/s2-products
path_S2tilinggrid <- '../../Sentinel-2_tiling_grid.kml'

# 3- get S2 acquisition and geometry of acquisition for study area
# https://shapps.dataspace.copernicus.eu/dashboard/#/account/settings
id <- "sh-3699814c-729b-4d29-abde-91163718cf5a"
secret <- "r8Z7znbv6A3uQHiFTFgnlg0mfvg4Y0Vy"
authentication <- list('id' = id, 'pwd' = secret)

get_s2_raster(aoi_path = aoi_path, datetime = datetime, output_dir = output_dir, 
              path_S2tilinggrid = path_S2tilinggrid, siteName = 'FrenchGuiana',
              overwrite = T, geomAcq = T, authentication = authentication)