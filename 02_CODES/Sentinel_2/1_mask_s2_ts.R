# ---
# title: "1_mask_s2_ts.R"
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
library(terra)
library(zoo)
source("../libraries/functions_sentinel_2.R")
source("../libraries/functions_general_tools.R")

# Directories & input data
input_dir <- '/media/corroyez/MyPassport/01_DATA'
output_dir <- '../../03_RESULTS'
sites <- c('Aigoual', 'Blois', 'Mormal')
# sites <- c('Mormal')
data <- read_ods(file.path(input_dir, "S2_Dates_Summer2021.ods"))
nb_bands <- 10

for (site in sites){
  print(paste("Processing site:", site))
  dateAcqs <- get_site_full_dates(site, data)
  for (dateAcq in dateAcqs){
    dateAcq <- as.Date(dateAcq)
    f_rast <- file.path(input_dir, site, "Sentinel-2", dateAcq, "raster_samples", 
                        paste0(site, "_001_", dateAcq, ".tiff"))
    f_mask <- file.path(input_dir, site, "Sentinel-2", dateAcq, "raster_samples", 
                        paste0(site, "_001_", dateAcq, "_BIN_v2.tiff"))
    if (file.exists(f_rast) && file.exists(f_mask)) {
      r <- rast(f_rast)
      m <- rast(f_mask)
      r_masked <- mask(r, m, maskvalues = 0, updatevalue = NA)
      out_path <- file.path(input_dir, site, "Sentinel-2", dateAcq)
      if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
      out_name <- file.path(out_path, paste0(site, "_", dateAcq, "_masked.tif"))
      print(paste("Saving:", out_name))
      writeRaster(r_masked, out_name, overwrite = TRUE)
    } else {
      warning(paste("File missing for site:", site, " Date:", dateAcq))
    }
  }
}