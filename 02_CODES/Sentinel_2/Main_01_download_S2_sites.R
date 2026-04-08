# clean environment
rm(list = ls(all=TRUE)); gc()
if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(rstactheia)
library(preprocS2)

geomAcq = T
mask_path = NULL
resolution = 10
fraction_vegetation = 5
doublecheckColl = F
offset = 1000
offset_B2 = F
RadiometricFilter = NULL
clean_bands = T
overwrite = F
corr_BRF = F
cloudcover = 95
bbox <- NULL
library(parallel)
library(future)
library(future.apply)
library(progressr)


# 1- define input & output directories
input_dir <- '../01_DATA/sites'
output_dir <- '../03_RESULTS'
sites <- c('Aigoual', 'Blois', 'Mormal')
dates <- list('Aigoual' = '2021-07-11', 
              'Blois' = '2021-06-14', 
              'Mormal' = '2021-06-14')
name_vect <- 'utm_init.shp'
# siteName = site

for (site in sites){
  # 2- define area of interest
  aoi_path <- file.path(input_dir, site, name_vect)
  output_subdir <- file.path('../03_RESULTS', site)
  
  # get tiling grid kml from https://sentiwiki.copernicus.eu/web/s2-products
  path_S2tilinggrid <- 'Sentinel-2_tiling_grid.kml'
  
  # 3- get S2 acquisition and geometry of acquisition for study area
  id <- "sh-3699814c-729b-4d29-abde-91163718cf5a"
  secret <- "r8Z7znbv6A3uQHiFTFgnlg0mfvg4Y0Vy"
  authentication <- list('id' = id, 'pwd' = secret)
  
  get_s2_raster(aoi_path = aoi_path, datetime = dates[[site]], 
                output_dir = output_subdir, path_S2tilinggrid = path_S2tilinggrid, 
                overwrite = overwrite, geomAcq = geomAcq, 
                authentication = authentication, siteName = site)
}
