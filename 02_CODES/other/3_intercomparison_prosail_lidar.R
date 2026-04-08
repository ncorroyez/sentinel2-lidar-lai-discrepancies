# ---
# title: "3.intercomparison_prosail_lidar.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-02-07"
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
source("functions_intercomparisons.R")

# ---------------------------------------- Preparation: Forest Mask & Directories ------------------------------------

# Directory to save
site <- "Mormal" # "Blois", "Aigoual
image <- "L2A_T31UER_A031222_20210614T105443"
resolution <- "res_20m"
results_path <- file.path("../03_RESULTS", site)

plots_dir <- file.path(results_path, "Plots")
# deciles_dir <- file.path(results_path, "LiDAR/Heterogeneity_Masks/Deciles")

# LiDAR and S-2 paths
image_path <- file.path(results_path, image, resolution)
lidar_output_dir <- file.path(results_path, "LiDAR/PAI/lidR")
s2_output_dir <- file.path(image_path)

# --------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------- ATBD vs LiDAR ----------------------------------------------------

# Mask
mask <- terra::rast(file.path(results_path, 
                              "LiDAR/Heterogeneity_Masks",
                              "artifacts_low_vegetation_majority_90_p_res_20_m.envi"))
# mask <- c("Refl, forest, low_vegetation")

estimated_var <- 'lai'
# estimated_vars <- c('lai', 'EWT', 'CHL', 'LMA')

# Bands and Noise Study

comparison_path <- file.path(plots_dir,
                             "@atbd_bands_noise_study")
dir.create(path = comparison_path, showWarnings = FALSE, recursive = TRUE)

bands <- c('bands_3_4_8', 'bands_3_4_5_6_7_8A_11_12')
noises <- c('addmult', 'mult')

row_names <- rep(noises, each = length(bands))
col_names <- rep(bands, length(noises))

# Intercomparisons

cat("Bands and Noise Study: Intercomparison", "\n")

correlation_atbd_bands_noise(bands,
                             noises,
                             row_names,
                             col_names,
                             mask,
                             plots_dir,
                             s2_output_dir,
                             image,
                             estimated_var,
                             comparison_path,
                             "atbd_bands_noise_correlation_matrix")

rmse_atbd_bands_noise(bands,
                      noises,
                      row_names,
                      col_names,
                      mask,
                      plots_dir,
                      s2_output_dir,
                      image,
                      estimated_var,
                      comparison_path,
                      "atbd_bands_noise_rmse_matrix")

bias_atbd_bands_noise(bands,
                      noises,
                      row_names,
                      col_names,
                      mask,
                      plots_dir,
                      s2_output_dir,
                      image,
                      estimated_var,
                      comparison_path,
                      "atbd_bands_noise_bias_matrix")

residuals_atbd_bands_noise(bands,
                           noises,
                           row_names,
                           col_names,
                           mask,
                           plots_dir,
                           s2_output_dir,
                           image,
                           estimated_var,
                           comparison_path)

error_atbd_bands_noise(bands,
                       noises,
                       row_names,
                       col_names,
                       mask,
                       plots_dir,
                       s2_output_dir,
                       image,
                       estimated_var,
                       comparison_path,
                       "atbd_bands_noise_error_matrix")

r_squared_atbd_bands_noise(bands,
                           noises,
                           row_names,
                           col_names,
                           mask,
                           plots_dir,
                           s2_output_dir,
                           image,
                           estimated_var,
                           comparison_path,
                           "atbd_bands_noise_r_squared_matrix")

std_atbd_bands_noise(bands,
                     noises,
                     row_names,
                     col_names,
                     mask,
                     plots_dir,
                     s2_output_dir,
                     image,
                     estimated_var,
                     comparison_path,
                     "atbd_bands_noise_std_matrix")


histograms_scatterplots_atbd_bands_noise(bands,
                                         noises, 
                                         row_names,
                                         col_names,
                                         mask,
                                         plots_dir,
                                         s2_output_dir,
                                         image,
                                         estimated_var,
                                         comparison_path)
# Refs

lidar_path <- file.path(lidar_output_dir, "lidar.envi")
lidar_raster <- open_raster_file(lidar_path)
lidar_raster <- terra::project(lidar_raster, mask)
lidar_raster <- lidar_raster * mask
lidar_raster_values <- values(lidar_raster)
valid_indices <- complete.cases(lidar_raster_values)

admd10m <- values(mask * terra::project(open_raster_file(file.path(s2_output_dir,
                                                                   "PRO4SAIL_INVERSION_atbd",
                                                                   "addmult",
                                                                   "bands_3_4_8",
                                                                   paste0(image, 
                                                                          '_', 
                                                                          'Refl',
                                                                          '_',
                                                                          estimated_var, 
                                                                          ".envi")) 
                                                         ),
                                        mask))
admd20m <- values(mask * terra::project(open_raster_file(file.path(s2_output_dir,
                                                                   "PRO4SAIL_INVERSION_atbd",
                                                                   "addmult",
                                                                   "bands_3_4_5_6_7_8A_11_12",
                                                                   paste0(image, 
                                                                          '_', 
                                                                          'Refl',
                                                                          '_',
                                                                          estimated_var, 
                                                                          ".envi")) 
                                                         ),
                                        mask))
md10m <- values(mask * terra::project(open_raster_file(file.path(s2_output_dir,
                                                                 "PRO4SAIL_INVERSION_atbd",
                                                                 "mult",
                                                                 "bands_3_4_8",
                                                                 paste0(image, 
                                                                        '_', 
                                                                        'Refl',
                                                                        '_',
                                                                        estimated_var, 
                                                                        ".envi")) 
                                                       ),
                                      mask))
md20m <- values(mask * terra::project(open_raster_file(file.path(s2_output_dir,
                                                                "PRO4SAIL_INVERSION_atbd",
                                                                "mult",
                                                                "bands_3_4_5_6_7_8A_11_12",
                                                                paste0(image, 
                                                                       '_', 
                                                                       'Refl',
                                                                       '_',
                                                                       estimated_var, 
                                                                       ".envi")) 
                                                       ),
                                      mask))

plot_histogram(vars_list = list(lidar_raster_values, admd10m, admd20m, md10m, md20m),
               var_labs = c("lidar", "admd10m", "admd20m", "md10m", "md20m"))
  
# Ref ATBD AddMult 3_4_8 
  
cat("Bands and Noise Study: Ref: ", "atbd_addmult_bands_3_4_8", "\n")

calculate_metrics_from_a_reference(file.path(s2_output_dir, 
                                             "PRO4SAIL_INVERSION_atbd",
                                             "addmult",
                                             "bands_3_4_8",
                                             paste0(image, 
                                                    '_', 
                                                    'Refl',
                                                    '_',
                                                    estimated_var, 
                                                    ".envi")),
                                   s2_output_dir,
                                   image,
                                   noises,
                                   bands,
                                   mask,
                                   valid_indices,
                                   estimated_var,
                                   "atbd_addmult_bands_3_4_8",
                                   comparison_path,
                                   "metrics_ref_atbd_addmult_bands_3_4_8")

# Ref LiDAR

cat("Bands and Noise Study: Ref: ", "lidar", "\n")

calculate_metrics_from_a_reference(lidar_path,
                                   s2_output_dir,
                                   image,
                                   noises,
                                   bands,
                                   mask,
                                   valid_indices,
                                   estimated_var,
                                   "lidar",
                                   comparison_path,
                                   "metrics_ref_lidar")

# Canopy, Leaf, and Literature
actuals <- c('atbd', 'lidar')
# actuals <- c('atbd')
predicteds <- c(
  "atbd_fixed_psi", "atbd_psi_low", "atbd_psi_high",
  "atbd_S2Geom", "atbd_JBGeom", "atbd",
  "atbd_fixed_tto", "atbd_tto_low", "atbd_tto_high",
  "atbd_fixed_tts", "atbd_tts_low", "atbd_tts_high",
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
  "atbd_psoil_high", "atbd_psoil_low", "atbd_psoil_full"
)

band_res <- "bands_3_4_8"
noise_applied <- "addmult"

estimated_vars <- c('lai')
# estimated_vars <- c('lai', 'EWT', 'CHL', 'LMA')

lidar_mask <- open_raster_file(file.path(lidar_output_dir, "lidar.envi"))
lidar_mask <- terra::project(lidar_mask, mask)
lidar_mask <- lidar_mask * mask

all_metrics <- list()

cat("Chosen Band Resolution", band_res, "\n")
cat("Chosen Noise", noise_applied, "\n")

for(estimated_var in estimated_vars){
  atbd_mask <- open_raster_file(
    file.path(s2_output_dir, 
              "PRO4SAIL_INVERSION_atbd",
              noise_applied,
              band_res,
              paste0(image, 
                     '_', 
                     'Refl',
                     '_',
                     estimated_var, 
                     ".envi")))
  atbd_mask <- terra::project(atbd_mask, mask)
  atbd_mask <- atbd_mask * mask 
  for(actual in actuals){
    if (actual == 'atbd') {
      reference <- atbd_mask
    } 
    else if (actual == 'lidar') {
      reference <- lidar_mask
    }
    else {
      stop("Error: Unknown mask type.\n")
    }
    for(predicted in predicteds){
      cat("Comparison between the reference", actual,
          "and the prediction", predicted, "\n")
      done_comparison <- paste("comparison",
                               "ref",
                               actual,
                               "pred",
                               predicted,
                               sep = "_")
      comparison_path <- file.path(plots_dir,
                                   done_comparison)
      dir.create(path = file.path(comparison_path),
                 showWarnings = FALSE, 
                 recursive = TRUE)
      predicted_mask <- open_raster_file(file.path(s2_output_dir, 
                                                   paste0("PRO4SAIL_INVERSION_",
                                                          predicted),
                                                   noise_applied,
                                                   band_res,
                                                   paste0(image,
                                                          "_",
                                                          "Refl",
                                                          "_",
                                                          estimated_var,
                                                          ".envi")))
      predicted_mask <- terra::project(predicted_mask, mask)
      predicted_mask <- predicted_mask * mask 
      metrics <- calculate_metrics(reference, 
                                   predicted_mask, 
                                   file.path(comparison_path,
                                             "metrics"))
      metrics$Actual <- actual
      metrics$Predicted <- predicted
      all_metrics[[length(all_metrics) + 1]] <- metrics
    }
  }
}

# Combine metrics from all iterations into a single data frame
combined_metrics <- do.call(rbind, all_metrics)

# Save combined metrics to a CSV file
combined_metrics_path <- file.path(plots_dir, "combined_metrics.csv")
write.csv(combined_metrics, file = combined_metrics_path, row.names = FALSE)
cat("Combined metrics saved to:", combined_metrics_path, "\n")
