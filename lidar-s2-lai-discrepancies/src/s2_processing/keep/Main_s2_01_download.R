# clean environment
rm(list = ls(all=TRUE)); gc()
if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(preprocS2)

# 1- Directories & input data
input_dir <- '../../01_DATA'
output_dir <- '../../03_RESULTS'
sites <- c('Aigoual', 
           'Blois', 
           'Mormal')
dates <- list('Aigoual' = '2021-07-11', 
              'Blois' = '2021-06-14', 
              'Mormal' = '2021-06-14')
name_vect <- 'utm_init.shp'
path_S2tilinggrid <- '../../01_DATA/Sentinel-2_tiling_grid.kml'

# 2- define authentication for CDSE
# run the following command 
# file.edit("~/.Renviron")
# add ID and password for CDSE as well as CURL_SSL_BACKEND using this template
# CDSE_ID="sh-XXXXXXXX-XXXX-XXXX-XXXX-XXXXXXXXXXXX"
# CDSE_SECRET="XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
# CURL_SSL_BACKEND="openssl"
# restart session to make sure .Renviron is updated
id <- Sys.getenv("CDSE_ID")
secret <- Sys.getenv("CDSE_SECRET")
authentication <- list('id' = id, 'pwd' = secret)

# 3- download S2 data
for (site in sites){
  # 3.1- get site vector
  aoi_path <- file.path(input_dir, site, "Geo_Files", name_vect)
  output_subdir <- file.path(output_dir, site, "PROSAIL_Optimization")
  if (!dir.exists(output_subdir)) 
    dir.create(dirname(output_subdir), recursive = TRUE, showWarnings = F)

  # 3.2- Get S2 acquisition & Geometry of acquisition for study area
  get_s2_raster(aoi_path = aoi_path, 
                datetime = dates[[site]], 
                output_dir = output_subdir, 
                path_S2tilinggrid = NULL,
                overwrite = T, 
                geomAcq = T, 
                authentication = authentication, 
                siteName = site,
                keepCRS = T)
}
