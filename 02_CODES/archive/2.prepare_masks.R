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

# Study site
site <- "Mormal" # "Blois", "Aigoual"

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

parameters <- c("brown", "chl", "ewt", "lai", "lidfa", "lma", "n", "psoil", "q")
estimated_parameters <- c('lai') #lai, CHL, EWT, LMA, fCover, fAPAR, albedo
ranges <- c("high", "low", "full")
noises <- c("addmult", "gaussian")
bands_res <- c("10m", "20m")

image <- 'L2A_T31UER_A031222_20210614T105443'

# Parameters-only
for (parameter in parameters){
  for (estimated_parameter in estimated_parameters){
    for (range in ranges){
      for (noise in noises){
        for (band_res in bands_res){
          filename <- paste("atbd", "modif",
                            parameter, range, noise, band_res,
                            "estimated", estimated_parameter, sep = "_")
          envi_file <- sprintf('L2A_T31UER_A031222_20210614T105443_Refl_%s.envi',
                               estimated_parameter)
          path_to_mask <- file.path(paste0(output_prosail, "atbd", "_", 
                                         parameter, "_", range), 
                                    noise,
                                    band_res,
                                    envi_file)
          output_masks <- file.path(output_dir, paste0("atbd", "_", parameter))
          rast <- generic_resample_and_mask(path_to_mask, 
                                            filename, 
                                            output_masks, 
                                            single_layer_forest)
        }
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