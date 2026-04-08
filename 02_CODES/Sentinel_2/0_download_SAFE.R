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
library(preprocS2)
source("../libraries/functions_sentinel_2.R")
source("../libraries/functions_general_tools.R")

# 1- Directories & input data
input_dir <- '/media/corroyez/MyPassport/01_DATA'
output_dir <- '../../03_RESULTS'
# sites <- c('Aigoual', 
#            'Blois', 
#            'Mormal')
sites <- c('Hayes', 'Reine')
# dates <- list('Aigoual' = '2021-07-11', 
#               'Blois' = '2021-06-14', 
#               'Mormal' = '2021-06-14')
data <- read_csv(file.path(input_dir, "S2_Dates.csv"))
name_vect <- 'utm_from_lidar.gpkg'

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
for (site in sites) {
  
  # Site AOI
  aoi_path <- file.path(input_dir, site, "Geo_Files", name_vect)
  aoi <- st_read(aoi_path, quiet = TRUE)
  
  # Dates for the site
  dateAcqs <- get_site_dates(site, data)
  
  for (dateAcq in dateAcqs) {
    output_subdir <- file.path(input_dir, site, "Sentinel-2", dateAcq)
    dir.create(output_subdir, recursive = TRUE, showWarnings = FALSE)
    sen2r(
      gui = FALSE,
      extent = aoi,
      timewindow = c(dateAcq, dateAcq),
      step_atmcorr = "l2a",
      list_prods = NULL,
      path_l2a = output_subdir,
      overwrite = FALSE,
      use_python = FALSE,
      preprocess = FALSE
    )
  }
}
