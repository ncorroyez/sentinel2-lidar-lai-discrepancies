# ---
# title: "3.calculate_lidar_metrics.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-07-22"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

library("lidR")
library("raster")
library("plotly")
library("terra")
library("viridis")
library("future")
library("rgl")
library("lmom")
library("data.table")
library("geojsonio")

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

# Import useful functions
source("../libraries/functions_lidar.R")
source("../libraries/functions_general_tools.R")
source("../libraries/functions_plots.R")

# Directories 
data <- "../../01_DATA"
results_dir <- "../../03_RESULTS"
figures_dir <- "../../04_FIGURES"
data_dir <- "/media/corroyez/My Passport/01_DATA"
# sites <- "Mormal" # Mormal Blois Aigoual
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- c("Blois", "Mormal")

# las <- readLAS(file.path(data_dir, "Aigoual/LiDAR/2-las_utm/LAS_741500_6336000.las"))
# las <- readLAS(file.path(data_dir, "Blois/LiDAR/2-las_utm/LAS_741500_6336000.las"))
# las <- readLAS(file.path(data_dir, "Mormal/LiDAR/2-las_utm/LAS_741500_6336000.las"))
# plot(las)

all_observations <- list()
all_predictors_filtered <- list()

for (site in sites){
  cat("Site:", site, "\n")
  data_site <- file.path(data, site)
  masks_dir <- file.path(results_dir, "LiDAR/Heterogeneity_Masks")
  
  # Samples
  site_data <- load_metrics_with_pad(site, results_dir, "Not_Masked")
  all_observations[[site]] <- site_data 
  
  geojson_file <- file.path(results_dir, site, "data_utm31n.geojson")
  geo_data <- geojson_read(geojson_file, what = "sp")
  coord_x <- geo_data@data$coord_x_utm31n
  coord_y <- geo_data@data$coord_y_utm31n
  coord_df <- data.frame(coord_x = coord_x, coord_y = coord_y)
  
  closest_matches <- list()
  for (j in 1:nrow(coord_df)) {
    target_coord <- coord_df[j, ]
    distances <- sqrt((site_data$lat - target_coord$coord_y)^2 +
                        (site_data$lon - target_coord$coord_x)^2)
    closest_index <- which.min(distances)
    closest_matches[[j]] <- cbind(site_data[closest_index, ], target_coord)
    closest_matches[[j]]$id_plot <- geo_data@data$id_plot[j]
  }
  closest_data <- do.call(rbind, closest_matches)
  # aaaaaaaa
  # closest_data <- closest_data[, !colnames(closest_data) %in% c("coord_x", "coord_y")]
  # closest_data <- site_data[site_data$lon %in% coord_df$coord_x & site_data$lat %in% coord_df$coord_y, ]
  closest_data <- closest_data[
  closest_data$stands %in% c("oak", "beech", "deciduous", "poplar") &
  # !is.na(closest_data$deltaLAI_dtm) &
  # !is.na(closest_data$deltaLAI_dsm) &
  closest_data$mean >= 2,
  ]
  all_predictors_filtered[[site]] <- closest_data
  # aaaaaaaaaa
  # LAS directory
  las_utm_dir <- file.path(data_dir,
                           site,
                           "LiDAR",
                           "3-las_normalized_utm"
  )
  las_utm_files <- list.files(las_utm_dir, pattern = "\\.las$", full.names = TRUE)
  
  # Catalog
  ctg <- readLAScatalog(las_utm_dir)
  
  # Resolution
  res <- 10
  
  # LiDR optimization
  opt_chunk_size(ctg) <- 0 # Processing by files
  opt_chunk_buffer(ctg) <- 10
  opt_filter(ctg) <- "-drop_class 2 -drop_z_below 2"
  # opt_select(ctg) <- "xyz" 
  
  for (i in 1:nrow(all_predictors_filtered[[site]])) {
    field_pt <- all_predictors_filtered[[site]][i, ]
    lon <- field_pt[["lon"]]
    lat <- field_pt[["lat"]]
    id_plot <- field_pt[["id_plot"]]
    lcv <- round(field_pt[["lcv"]], 2)
    lskew <- round(field_pt[["lskew"]], 2)
    mean <- round(field_pt[["mean"]], 2)
    max <- round(field_pt[["max"]], 2)
    vci <- round(field_pt[["vci"]], 2)
    lai <- round(field_pt[["lidar_lai"]], 2)
    las <- clip_circle(ctg, lon, lat, 5)
    # plot(las)
    # aaaaa
    output_file <- file.path(figures_dir, site, "clouds",
                             paste0("id_plot_", id_plot, "_x_", lon, "_y_", lat,
                                    "_lcv_", lcv, "_lskew_", lskew, 
                                    "_mean_", mean, "_max_", max,
                                    "_vci_", vci, "_lai_", lai,
                                    ".png"))
    # png(output_file, width = 800, height = 600)
    # open3d()
    plot(las)
    # par3d(windowRect = c(100, 100, 1600, 1200))
    # dev.off()
    rgl.snapshot(output_file)
    # ggsave(output_file, a)
    cat("Saved plot for site", site, "point", i, "to", output_file, "\n")
    # aaa
  }
}
combined_predictors_filtered <- do.call(rbind, all_predictors_filtered)
combined_predictors_all <- do.call(rbind, all_observations)