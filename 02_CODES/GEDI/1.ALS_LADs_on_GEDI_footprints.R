# ---
# title: "1.ALS_LADs_on_GEDI_footprints.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2025-12-02"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

library(dplyr)
library(lubridate)
library(sf)
library(terra)
library(exactextractr)
library(randomForest)
library(caret)
library(broom)
library(stringr)
library(ggplot2)
library(tidyr)
library(lidR)
source("../libraries/GEDI/functions_GEDI.R")

set.seed(42)
data_dir <- "/media/corroyez/MyPassport/01_DATA"
metrics_dir <- "Deciduous_Only"
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- c("Blois")  # uncomment to test only one site

bin_breaks <- seq(0, 40, by = 5)
gedi_radius = 12.5

for (site in sites) {
  start_time <- Sys.time()
  message("Processing site: ", site)
  las_dir <- file.path(data_dir, site, "LiDAR/3-las_normalized_utm")
  results_dir <- file.path("../../03_RESULTS", site, "Metrics", metrics_dir)
  save_dir <- file.path("../../01_DATA", site, "GEDI")
  
  gedi_data <- load_gedi_and_geogedi_data(site, data_dir, results_dir,
                                          bin_breaks, gedi_radius)
  gedi_buf  <- gedi_data$buffer
  
  # Create LiDAR Catalog
  ctg <- readLAScatalog(las_dir)
  opt_select(ctg) <- "xyzcr" 
  opt_progress(ctg) <- FALSE
  message("Processing ", nrow(gedi_buf), " GEDI footprints...")
  
  # Iterate and Clip
  lad_results_list <- vector("list", nrow(gedi_buf))
  for (i in 1:nrow(gedi_buf)) {
    
    if (i %% 10 == 0) message(paste0("  Processing ", i, " / ",
                                     nrow(gedi_buf), "\n"))
    
    #   pt_geom <- st_centroid(gedi_buf$geometry[i])
    #   coords  <- st_coordinates(pt_geom)
    #   pt_x <- coords[1]
    #   pt_y <- coords[2]
    #   las_circle <- clip_circle(ctg, xcenter = pt_x, ycenter = pt_y, 
    #                             radius = gedi_radius)
    #   
    #   # Calculate metrics if data exists
    #   if (!is.null(las_circle) && !is.empty(las_circle)) {
    #     # Safety check
    #     metrics <- tryCatch({
    #       calculate_lad_metrics(las_circle, bin_size = 5)
    #     }, error = function(e) return(NULL))
    #     
    #     lad_results_list[[i]] <- metrics
    #   }
    # }
    
    # -------------------------------------------------------
    # A. PROCESS ORIGINAL COORDINATES (X_utm, Y_utm)
    # -------------------------------------------------------
    x_org <- gedi_buf$x_utm[i]
    y_org <- gedi_buf$y_utm[i]
    
    metrics_org <- NULL
    
    # Clip & Calculate
    las_org <- clip_circle(ctg, xcenter = x_org, ycenter = y_org, radius = gedi_radius)
    if (!is.null(las_org) && !is.empty(las_org)) {
      metrics_org <- tryCatch(calculate_lad_metrics(las_org, bin_size = 5), error=function(e) NULL)
    }
    
    # Rename with suffix (e.g., LAI_lidar_org)
    metrics_org <- add_suffix(metrics_org, "_org")
    
    # -------------------------------------------------------
    # B. PROCESS CORRECTED COORDINATES (x_corrected, y_corrected)
    # -------------------------------------------------------
    x_cor <- gedi_buf$x_corrected[i]
    y_cor <- gedi_buf$y_corrected[i]
    
    metrics_cor <- NULL
    
    # Only process if coordinates exist (not NA)
    if (!is.na(x_cor) && !is.na(y_cor)) {
      las_cor <- clip_circle(ctg, xcenter = x_cor, ycenter = y_cor, radius = gedi_radius)
      
      if (!is.null(las_cor) && !is.empty(las_cor)) {
        metrics_cor <- tryCatch(calculate_lad_metrics(las_cor, bin_size = 5), error=function(e) NULL)
      }
    }
    
    # Rename with suffix (e.g., LAI_lidar_cor)
    metrics_cor <- add_suffix(metrics_cor, "_cor")
    
    # -------------------------------------------------------
    # C. MERGE BOTH
    # -------------------------------------------------------
    # c(list, list) combines them into one flat list
    combined_metrics <- c(metrics_org, metrics_cor)
    
    if (length(combined_metrics) > 0) {
      lad_results_list[[i]] <- combined_metrics
    }
  }
  
  # Merge Results
  lad_df <- bind_rows(lad_results_list)
  
  if(nrow(lad_df) > 0) {
    # Combine with original GEDI metadata
    # Ensure row alignment matches (list index matches row index)
    gedi_final <- cbind(gedi_buf, lad_df)
    print(head(gedi_final))
    saveRDS(gedi_final, file.path("../../01_DATA", paste0(site, "_GEDI_LiDAR_Metrics.rds")))
    message("Saved results to: ", paste0(site, "_GEDI_LiDAR_Metrics_Both.rds"))
  } else {
    warning("No results generated.")
  }
  duration <- Sys.time() - start_time
  message("Duration: ", round(duration, 2), " ", units(duration))
}
