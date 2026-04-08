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
library("readr")
source("../libraries/functions_sentinel_2.R")
source("../libraries/functions_general_tools.R")

# Pre-processing Parameters
data_dir <- '/media/corroyez/My Passport/01_DATA'
# data_dir <- '../../01_DATA'
sites <- c("Aigoual", "Blois", "Mormal")
sites <- "Aigoual" # Mormal Blois Aigoual

for (site in sites){
  
  # Setup 
  results_dir <- '../../03_RESULTS'
  save_dir <- file.path(data_dir, site, "S2_Images")
  
  # Masks
  masks_dir <- file.path(results_dir, site, "LiDAR/Heterogeneity_Masks")
  metrics_dir <- file.path(results_dir, site, "Metrics")
  if (!dir.exists(metrics_dir)) {
    dir.create(metrics_dir, showWarnings = FALSE, recursive = TRUE)
  }
  raw_dir <- file.path(metrics_dir, "Raw")
  if (!dir.exists(raw_dir)) {
    dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Extract all dates for the current site
  dateAcqs <- '2021-11-18'
  
  # Dates: yyyy-mm-dd format mandatory
  # if (site == "Mormal") {
  #   dateAcqs <- c('2021-06-14', '2021-12-21')
  # } else if (site == "Blois") {
  #   dateAcqs <- c('2021-06-14', '2021-12-21')
  # } else if (site == "Aigoual") {
  #   dateAcqs <- c('2021-07-11', '2021-12-21')
  # } else {
  #   stop("Unknown site. Please provide a valid site name.")
  # }
  
  # resolutions <- c(10,20)
  resolutions <- c(10) # 10 20
  
  for (resolution in resolutions){
    for (dateAcq in dateAcqs){
      if (resolution == 10){
        cat("Chosen resolution", resolution, "\n")
      } else if (resolution == 20){
        cat("Chosen resolution", resolution, "\n")
      } else{
        stop("Error: Resolution must be 10 or 20.\n")
      }
      s2_creation_directory <- paste("../../01_DATA", 
                                     site, 
                                     'S2_Images', 
                                     dateAcq,
                                     sep = "/")
      if (!dir.exists(s2_creation_directory)) {
        dir.create(s2_creation_directory,
                   showWarnings = FALSE,
                   recursive = TRUE)
      }
      results_path <- paste(results_dir, site, sep = '/')
      S2source <- 'SAFE'
      saveRaw <- TRUE
      
      # S-2 Pre-Processing: Cloud Mask, Reflectance
      results <- preprocess_S2_no_dl(path_vector = paste("../../01_DATA",
                                                         site,
                                                         "Shapefiles/utm_init.shp",
                                                         sep = "/"),
                                     s2_creation_directory,
                                     results_path,
                                     resolution = resolution,
                                     S2source = S2source, # SAFE THEIA
                                     saveRaw = TRUE)
      # copy_and_move_folder(source = s2_creation_directory,
      #                      destination = save_dir)
    }
  }
}
