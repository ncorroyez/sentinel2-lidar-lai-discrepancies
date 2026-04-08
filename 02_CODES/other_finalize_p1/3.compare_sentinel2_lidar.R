# ---
# title: "3.compare_sentinel2_lidar.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-01-24"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

source("libraries/functions_plots.R")

# Setup
results_dir <- "../03_RESULTS"
metrics_dir <- "Metrics/Full_Composition"
sites <- c("Aigoual", "Blois", "Mormal")
atbd_dir <- "res_10m/PRO4SAIL_INVERSION_atbd/addmult/bands_3_4_8"
atbd_clumping_uniform_dir <- "res_10m/PRO4SAIL_INVERSION_atbd_clumping_uniform/addmult/bands_3_4_8"
atbd_clumping_gaussian_dir <- "res_10m/PRO4SAIL_INVERSION_atbd_clumping_gaussian/addmult/bands_3_4_8"

for (site in sites){
  cat(paste("Study Site:", site), "\n")
  if(site == "Aigoual"){
    image_name <- "L2A_T31TEJ_A031608_20210711T104217"
  }
  else if(site == "Blois"){
    image_name <- "L2A_T31TCN_A031222_20210614T105443"
  }
  else if(site == "Mormal"){
    image_name <- "L2A_T31UER_A031222_20210614T105443"
  }
  else{
    stop("Error: Site is wrong")
  }
  
  lidar_lai <- terra::rast(file.path(results_dir,
                                     site,
                                     metrics_dir,
                                     "lidarlai_res_10_m.tif"))
  lidar_lai_values <- values(lidar_lai)
  
  lidar_lai_top_canopy <- terra::rast(file.path(results_dir,
                                                site,
                                                metrics_dir,
                                                "lidarlaitopcanopy_res_10_m.tif"))
  lidar_lai_top_canopy_values <- values(lidar_lai_top_canopy)
  
  mask <- terra::rast(file.path(results_dir, 
                                site,
                                "LiDAR/Heterogeneity_Masks",
                                "artifacts_full_composition_low_vegetation_majority_90_p_res_10_m.envi"))
  
  s2_atbd_lai <- terra::rast(file.path(results_dir,
                                       site,
                                       image_name,
                                       atbd_dir,
                                       paste0(image_name, "_Refl_lai.envi")))
  s2_atbd_lai <- mask(s2_atbd_lai, mask)
  s2_atbd_lai_values <- values(s2_atbd_lai)
  
  s2_atbd_clumping_uniform_lai <- terra::rast(file.path(results_dir,
                                                        site,
                                                        image_name,
                                                        atbd_clumping_uniform_dir,
                                                        paste0(image_name, "_Refl_lai.envi")))
  s2_atbd_clumping_uniform_lai <- mask(s2_atbd_clumping_uniform_lai, mask)
  s2_atbd_clumping_uniform_lai_values <- values(s2_atbd_clumping_uniform_lai)
  
  s2_atbd_clumping_gaussian_lai <- terra::rast(file.path(results_dir,
                                                         site,
                                                         image_name,
                                                         atbd_clumping_gaussian_dir,
                                                         paste0(image_name, "_Refl_lai.envi")))
  s2_atbd_clumping_gaussian_lai <- mask(s2_atbd_clumping_gaussian_lai, mask)
  s2_atbd_clumping_gaussian_lai_values <- values(s2_atbd_clumping_gaussian_lai)
  
  cat(paste("Correlation LiDAR LAI - S2 LAI", 
              cor(lidar_lai_values, s2_atbd_lai_values, use = "pairwise"), "\n"))
  cat(paste("Correlation LiDAR LAI - S2 LAI Clumping Uniform", 
              cor(lidar_lai_values, s2_atbd_clumping_uniform_lai_values, use = "pairwise"), "\n"))
  cat(paste("Correlation LiDAR LAI - S2 LAI Clumping Gaussian",  
              cor(lidar_lai_values, s2_atbd_clumping_gaussian_lai_values, use = "pairwise"), "\n\n"))
  
  cat(paste("Correlation LiDAR LAI Top Canopy - S2 LAI", 
              cor(lidar_lai_top_canopy_values, s2_atbd_lai_values, use = "pairwise"), "\n"))
  cat(paste("Correlation LiDAR LAI Top Canopy - S2 LAI Clumping Uniform", 
              cor(lidar_lai_top_canopy_values, s2_atbd_clumping_uniform_lai_values, use = "pairwise"), "\n"))
  cat(paste("Correlation LiDAR LAI Top Canopy - S2 LAI Clumping Gaussian",  
              cor(lidar_lai_top_canopy_values, s2_atbd_clumping_gaussian_lai_values, use = "pairwise"), "\n\n"))
}
