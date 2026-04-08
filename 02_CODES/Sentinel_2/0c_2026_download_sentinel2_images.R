# ---
# title: "0_download_sentinel2_images.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-07-22"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# Import useful functions
library(readr)
library(readODS)
library(preprocS2)
source("../libraries/functions_sentinel_2.R")
source("../libraries/functions_general_tools.R")

# Directories & input data
input_dir <- '/media/corroyez/MyPassport/01_DATA'
output_dir <- '../../03_RESULTS'
sites <- c('Aigoual', 'Blois', 'Mormal')
# sites <- c('Hayes', 'Reine')
# sites <- c('Mormal')
# sites <- c('Hayes')
s2_dates <- read_ods(file.path(input_dir, "S2_Dates_Summer2021.ods"))
name_vect <- 'utm_from_lidar.gpkg'

# Define authentication for CDSE
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

options <- set_options_preprocS2(fun = 'get_s2_raster')
options$overwrite <- F
options$geom_acq <- T

# Download S2 data
for (site in sites){
  print(site)
  dateAcqs <- get_site_full_dates(site, s2_dates)
  # stop()
  for (dateAcq in as.character(dateAcqs)){
    print(dateAcq)
    # Get site vector
    aoi_path <- file.path(input_dir, site, "Geo_Files", name_vect)
    output_subdir <- file.path(input_dir, site, "Sentinel-2", dateAcq)
    if (!dir.exists(output_subdir)) 
      dir.create(dirname(output_subdir), recursive = TRUE, showWarnings = F)
    # stop()
    # Get S2 acquisition & Geometry of acquisition for study area
    a <- get_s2_raster(aoi_path = aoi_path, 
                       datetime = dateAcq, 
                       output_dir = output_subdir, 
                       options = options,
                       site_name = site
    )
  }
}

# apply radiometric filter
# plot 53 LAI TS
# lai s2 pas corrigé dans musica
# hobo en scénario
# lai_2_4_6_8
# gapfilled / corr lai lidar total / corr prof / corr RF