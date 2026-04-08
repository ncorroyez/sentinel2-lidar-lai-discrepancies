# ---
# title: "2.prepare_masks"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-01-30"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

source("functions_plots.R")

# ---------------------------------------- Preparation: Forest Mask & Directories ------------------------------------

# Get one band only of Forest Mask
forest_path <- "../03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/Forest_Mask/forest.envi"
single_layer_forest <- subset(terra::rast(forest_path), 1)

# Directory to save
results_dir <- '../03_RESULTS/Mormal'
output_prosail <- file.path(results_dir, 'L2A_T31UER_A031222_20210614T105443/PRO4SAIL_INVERSION_')
output_dir <- file.path(results_dir,"LAI_PAI_Masks/LAI_S2")

# LiDAR
pai_lidar_fdb_path <- "../03_RESULTS/Mormal/LAI_PAI_Masks/PAI_LiDAR/fdeboissieu/pai_las_cleaned.tif"
lidar_mask_fdb_path <- "../03_RESULTS/Mormal/LAI_PAI_Masks/PAI_LiDAR/fdeboissieu"
lidar_mask <- generic_resample_and_mask(pai_lidar_fdb_path, 
                                  "LiDAR", 
                                  lidar_mask_fdb_path, 
                                  single_layer_forest)

# Low vegetation mask
combine_path <- "../03_RESULTS/Mormal/LiDAR_Masks/combine.tif"
envi_pai_path <- "../03_RESULTS/Mormal/LAI_PAI_Masks/PAI_LiDAR/fdeboissieu/lidar.envi"

combine <- open_envi_file_as_values(combine_path)
pai_matrix <- open_envi_file_as_values(envi_pai_path)
pai_ref <- open_envi_file(envi_pai_path)

# LiDAR
masked_pai <- mask_low_vegetation(pai_matrix, pai_ref, combine)
masked_pai_raster <- masked_pai$masked_raster
masked_pai_values <- masked_pai$masked_values

plot_histogram_2_vars(masked_pai_values, pai_matrix, 
                      "Histogram of Masked LiDAR PAI Values and LiDAR PAI Values", "PAI",
                      "Masked LiDAR", "LiDAR",
                      "../03_RESULTS/Mormal/Plots", "hist_lai_maskedlidar_lidar")

parameters <- c("brown", "chl", "ewt", "lai", "lidfa", "lma", "n", "psoil", "q")
ranges <- c("high", "low", "full")

noises <- c("addmult", "gaussian")
bands_res <- c("10m", "20m")

image <- 'L2A_T31UER_A031222_20210614T105443'
lai_name <- 'Refl_lai'

# Sentinel-2 ATBD
filename <- paste("atbd", "addmult", "10m", sep = "_")
lai_file <- 'L2A_T31UER_A031222_20210614T105443_Refl_lai.envi'
path_to_mask <- file.path(paste0(output_prosail, "atbd"), 
                          "addmult",
                          "10m",
                          lai_file)
output_masks <- file.path(output_dir, paste0("atbd", "_", "low_vegetation"))
rast <- low_vegetation_resample_and_mask(path_to_mask, 
                                         filename, 
                                         output_masks, 
                                         single_layer_forest,
                                         combine)

# Sentinel-2 Parameters-only
for (parameter in parameters){
  for (range in ranges){
    for (noise in noises){
      for (band_res in bands_res){
        filename <- paste("atbd", parameter, range, noise, band_res, sep = "_")
        lai_file <- 'L2A_T31UER_A031222_20210614T105443_Refl_lai.envi'
        path_to_mask <- file.path(paste0(output_prosail, "atbd", "_", 
                                         parameter, "_", range), 
                                  noise,
                                  band_res,
                                  lai_file)
        output_masks <- file.path(output_dir, paste0("atbd", "_", parameter,
                                                     "_", "low_vegetation"))
        rast <- low_vegetation_resample_and_mask(path_to_mask, 
                                          filename, 
                                          output_masks, 
                                          single_layer_forest,
                                          combine)
      }
    }
  }
}

distribs <- c("atbd",
              "q_zhang_et_al_2005", "brede_et_al_2020", 
              "hauser_et_al_2021", "sinha_et_al_2020",
              "verhoef_and_bach_2007", "shiklomanov_et_al_2016",
              "atbd_n_high", "atbd_n_low", "atbd_n_full",
              "atbd_chl_high", "atbd_chl_low", "atbd_chl_full",
              "atbd_brown_high", "atbd_brown_low", "atbd_brown_full",
              "atbd_ewt_high", "atbd_ewt_low", "atbd_ewt_full",
              "atbd_lma_high", "atbd_lma_low", "atbd_lma_full",
              "atbd_lidfa_high", "atbd_lidfa_low", "atbd_lidfa_full",
              "atbd_lai_high", "atbd_lai_low", "atbd_lai_full",
              "atbd_q_high", "atbd_q_low", "atbd_q_full",
              "atbd_psoil_high", "atbd_psoil_low", "atbd_psoil_full")

# Generic
# for (distrib in distribs){
#   for (noise in noises){
#     for (band_res in bands_res){
#       filename <- paste(distrib, noise, band_res, sep = "_")
#       lai_file <- 'L2A_T31UER_A031222_20210614T105443_Refl_lai.envi'
#       path_to_mask <- file.path(paste0(output_prosail, distrib), 
#                                 noise,
#                                 band_res,
#                                 lai_file)
#       output_masks <- file.path(output_dir, distrib)
#       rast <- generic_resample_and_mask(path_to_mask, 
#                                         filename, 
#                                         output_masks, 
#                                         single_layer_forest)
#     }
#   }
# }