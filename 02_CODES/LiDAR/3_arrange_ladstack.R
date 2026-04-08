# ---
# title: "3b.arrange_ladstack.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2025-03-10"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

library("lidR")
library("lidRmetrics")
library("raster")
library("plotly")
library("terra")
library("viridis")
library("future")
library("rgl")
library("lmom")
library("data.table")
library("solaR")

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

# Import useful functions
source("../libraries/functions_lidar.R")
source("../libraries/functions_create_masks.R")
source("../libraries/functions_plots.R")

# Directories
data <- "../../01_DATA"
# sites <- "Aigoual" # Mormal Blois Aigoual
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- c("Blois", "Mormal")
# sites <- c('Aigoual', 'Blois')
# sites <- c("Mormal")
# sites <- c("Blois")
forest_composition <- "Metrics/Deciduous_Only"

for (site in sites){
  results_path <- file.path("../../03_RESULTS", site)
  ladstack_dsm_rast <- terra::rast(file.path(results_path, 
                                             forest_composition, 
                                             "ladstack_dsm.tif"))
                                         # "Metrics/Deciduous_Only/20resladstack_dsm.tif"))
  ladstack_dtm_rast <- terra::rast(file.path(results_path, 
                                             forest_composition, 
                                             "ladstack_dtm.tif"))
  lidarlai_rast <- terra::rast(file.path(results_path, 
                                         forest_composition, 
                                         "lidarlai_res_10_m.tif"))
  max_rast <- terra::rast(file.path(results_path, 
                                    forest_composition, 
                                    "max_res_10_m.tif"))
  max_rast[max_rast > 40] <- 40
  max_values <- values(max_rast)
  
  lidarlai_rast[max_values < 20] <- NA
  writeRaster(lidarlai_rast, file.path(results_path, 
                                       "lidarlai_h_inf_20_res_10_m.tif"),
              overwrite = T)
  
  dsm_rast <- ladstack_dsm_rast[[1]] # * 0
  dsm_keepTrees_rast <- ladstack_dsm_rast[[1]] # * 0
  
  dtm_rast <- ladstack_dtm_rast[[1]] # * 0
  dtm_keepTrees_rast <- ladstack_dtm_rast[[1]] # * 0
  
  # DSM
  for (i in 1:nlyr(ladstack_dsm_rast)) {
    
    # current_layer <- ifel(is.na(ladstack_dsm_rast[[i]]), NA, ladstack_dsm_rast[[i]])
    current_layer <- ladstack_dsm_rast[[i]]
    
    if (i == 1){
      current_layer <- 0
    }
    
    # DSM Logic: If max value is not reached, accumulate (sum); else set to NA for higher layers
    dsm_rast <- ifel(ceiling(max_rast) >= i, dsm_rast + current_layer, NA)
    # dsm_rast[current_layer >= max_rast & is.na(dsm_rast)] <- NA  # For higher layers, set NA
    
    # DSM_keepTrees Logic: If max value is not reached, accumulate (sum); else add 0 for higher layers
    dsm_keepTrees_rast <- ifel(ceiling(max_rast) >= i, dsm_keepTrees_rast + current_layer, dsm_keepTrees_rast)
    # dsm_keepTrees_rast[current_layer >= max_rast & is.na(dsm_keepTrees_rast)] <- 0  # For higher layers, add 0
    
    # if (i == 20){
    #   stop()
    # }
    d <- 40 - i + 0.5
    # DSM_keepTrees
    subraster_filename <- file.path(results_path,
                                    forest_composition,
                                    "testPADs",
                                    "PAD_Profiles_DSM_keepTrees",
                                    paste0("PAD_", d, "_40.tif"))
    dir.create(path = subraster_filename,
               showWarnings = FALSE,
               recursive = TRUE)
    writeRaster(dsm_keepTrees_rast, subraster_filename, overwrite=TRUE)
    
    # DSM_keepTreesAbove20
    subraster_filename20 <- file.path(results_path,
                                      forest_composition,
                                      "testPADs",
                                      "PAD_Profiles_DSM_keepTreesAbove20",
                                      paste0("PAD_", d, "_40.tif"))
    writeRaster(ifel(max_rast < 20, NA, dsm_keepTrees_rast),
                subraster_filename20, overwrite=TRUE)
    
    # DSM
    subraster_filename_dsm <- file.path(results_path,
                                        forest_composition,
                                        "testPADs",
                                        "PAD_Profiles_DSM",
                                        paste0("PAD_", d, "_40.tif"))
    writeRaster(dsm_rast, subraster_filename_dsm, overwrite=TRUE)
    
    # DSM_Above20
    subraster_filename2 <- file.path(results_path,
                                     forest_composition,
                                     "testPADs",
                                     "PAD_Profiles_DSM_Above20",
                                     paste0("PAD_", d, "_40.tif"))
    writeRaster(ifel(max_rast < 20, NA, dsm_rast),
                subraster_filename2, overwrite=TRUE)
  }
  
  # DTM
  for (i in 1:nlyr(ladstack_dtm_rast)) {
    
    # current_layer <- ifel(is.na(ladstack_rast[[i]]), NA, ladstack_rast[[i]])
    current_layer <- ladstack_dtm_rast[[i]]
    
    if (i == 1){
      current_layer <- 0
    }
    
    # DSM Logic: If max value is not reached, accumulate (sum); else set to NA for higher layers
    dtm_rast <- ifel(ceiling(max_rast) >= i, dtm_rast + current_layer, NA)
    # dsm_rast[current_layer >= max_rast & is.na(dsm_rast)] <- NA  # For higher layers, set NA
    
    # DSM_keepTrees Logic: If max value is not reached, accumulate (sum); else add 0 for higher layers
    dtm_keepTrees_rast <- ifel(ceiling(max_rast) >= i, dtm_keepTrees_rast + current_layer, dtm_keepTrees_rast)
    # dsm_keepTrees_rast[current_layer >= max_rast & is.na(dsm_keepTrees_rast)] <- 0  # For higher layers, add 0
    
    # if (i == 20){
    #   stop()
    # }
    d <- 40 - i + 0.5
    # DSM_keepTrees
    subraster_filename <- file.path(results_path,
                                    forest_composition,
                                    "testPADs",
                                    "PAD_Profiles_DTM_keepTrees",
                                    paste0("PAD_", d, "_40.tif"))
    writeRaster(dtm_keepTrees_rast, subraster_filename, overwrite=TRUE)
    
    # DSM_keepTreesAbove20
    subraster_filename20 <- file.path(results_path,
                                      forest_composition,
                                      "testPADs",
                                      "PAD_Profiles_DTM_keepTreesAbove20",
                                      paste0("PAD_", d, "_40.tif"))
    writeRaster(ifel(max_rast < 20, NA, dtm_keepTrees_rast),
                subraster_filename20, overwrite=TRUE)
    
    # DSM
    subraster_filename_dtm <- file.path(results_path,
                                        forest_composition,
                                        "testPADs",
                                        "PAD_Profiles_DTM",
                                        paste0("PAD_", d, "_40.tif"))
    writeRaster(dtm_rast, subraster_filename_dtm, overwrite=TRUE)
    
    # DSM_Above20
    subraster_filename2 <- file.path(results_path,
                                     forest_composition,
                                     "testPADs",
                                     "PAD_Profiles_DTM_Above20",
                                     paste0("PAD_", d, "_40.tif"))
    writeRaster(ifel(max_rast < 20, NA, dtm_rast),
                subraster_filename2, overwrite=TRUE)
  }
}

# cumulative_sum <- cumulative_sum + current_layer
# cumulative_sum2 <- cumulative_sum2 + ladstack_rast[[i]]