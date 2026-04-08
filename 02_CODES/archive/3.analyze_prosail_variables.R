# ---
# title: "3.analyze_prosail_variables"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-01-31"
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
results_path <- "../03_RESULTS/Mormal"
plots_dir <- file.path(results_path, "Plots")
deciles_dir <- file.path(results_path, "/LiDAR_Masks/Deciles")

input_base_path <- file.path(results_path, "L2A_T31UER_A031222_20210614T105443")
lidar_output_dir <- file.path(results_path, "Biophysical_S2_LiDAR_Masks
                                             /PAI_LiDAR/fdeboissieu")
output_directory <- file.path(results_path, "Biophysical_S2_LiDAR_Masks/
                                             Biophysical_Variables_S2")

# --------------------------------------------------------------------------------------------------------------------

# ------------------------------------- Preparation: Get LiDAR & Sentinel-2 Mask Values ------------------------------

# LiDAR

lidar_mask <- open_raster_file_as_values(file.path(lidar_output_dir, "lidar.envi"))

# Sentinel-2

atbd_mask <- open_raster_file_as_values(file.path(output_directory, "atbd/atbd_addmult_10m.envi"))

plot_histogram(atbd_mask, "Histogram of Difference between ATBD 3 bands and ATBD Values",
                     "LAI", c("Difference ATBDs", "a", "rea", "raed"), 
                     plots_dir, "test")

# BROWN

analyze_one_parameter("brown",
                      output_directory,
                      plots_dir)

# CHL

analyze_one_parameter("chl",
                      output_directory,
                      plots_dir)

# EWT

analyze_one_parameter("ewt",
                      output_directory,
                      plots_dir)

# LAI

analyze_one_parameter("lai",
                      output_directory,
                      plots_dir)

# LIDFa

analyze_one_parameter("lidfa",
                      output_directory,
                      plots_dir)

# LMA

analyze_one_parameter("lma",
                      output_directory,
                      plots_dir)

# N

analyze_one_parameter("n",
                      output_directory,
                      plots_dir)

# psoil

analyze_one_parameter("psoil",
                      output_directory,
                      plots_dir)

# q

analyze_one_parameter("q",
                      output_directory,
                      plots_dir)


# --------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------- LAI Comparison --------------------------------------------------

path_atbd_addmult_20m <- '../03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/PRO4SAIL_INVERSION_atbd/addmult/20m/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi'
rast_atbd_addmult_20m <- terra::rast(path_atbd_addmult_20m)
rast_atbd_addmult_20m_resampled <- resample(rast_atbd_addmult_20m, single_layer_forest)
atbd_addmult_20m_mask <- mask(rast_atbd_addmult_20m_resampled, single_layer_forest)
save_envi_file(atbd_addmult_20m_mask, "atbd_addmult_20m", "../03_RESULTS/Mormal/LAI_PAI_Masks/LAI_S2", "ATBD")
atbd_addmult_20m_mask <- values(atbd_addmult_20m_mask)
atbd_addmult_20m_mask[atbd_addmult_20m_mask < 0] <- 0

path_atbd_addmult_10m <- '../03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/PRO4SAIL_INVERSION_atbd/addmult/10m/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi'
rast_atbd_addmult_10m <- terra::rast(path_atbd_addmult_10m)
rast_atbd_addmult_10m_resampled <- resample(rast_atbd_addmult_10m, single_layer_forest)
atbd_addmult_10m_mask <- mask(rast_atbd_addmult_10m_resampled, single_layer_forest)
save_envi_file(atbd_addmult_10m_mask, "atbd_addmult_10m", "../03_RESULTS/Mormal/LAI_PAI_Masks/LAI_S2", "ATBD")
atbd_addmult_10m_mask <- values(atbd_addmult_10m_mask)
atbd_addmult_10m_mask[atbd_addmult_10m_mask < 0] <- 0

path_atbd_gaussian_20m <- '../03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/PRO4SAIL_INVERSION_atbd/gaussian/20m/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi'
rast_atbd_gaussian_20m <- terra::rast(path_atbd_gaussian_20m)
rast_atbd_gaussian_20m_resampled <- resample(rast_atbd_gaussian_20m, single_layer_forest)
atbd_gaussian_20m_mask <- mask(rast_atbd_gaussian_20m_resampled, single_layer_forest)
save_envi_file(atbd_gaussian_20m_mask, "atbd_gaussian_20m", "../03_RESULTS/Mormal/LAI_PAI_Masks/LAI_S2", "ATBD")
atbd_gaussian_20m_mask <- values(atbd_gaussian_20m_mask)
atbd_gaussian_20m_mask[atbd_gaussian_20m_mask < 0] <- 0

path_atbd_gaussian_10m <- '../03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/PRO4SAIL_INVERSION_atbd/gaussian/10m/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi'
rast_atbd_gaussian_10m <- terra::rast(path_atbd_gaussian_10m)
rast_atbd_gaussian_10m_resampled <- resample(rast_atbd_gaussian_10m, single_layer_forest)
atbd_gaussian_10m_mask <- mask(rast_atbd_gaussian_10m_resampled, single_layer_forest)
save_envi_file(atbd_gaussian_10m_mask, "atbd_gaussian_10m", "../03_RESULTS/Mormal/LAI_PAI_Masks/LAI_S2", "ATBD")
atbd_gaussian_10m_mask <- values(atbd_gaussian_10m_mask)
atbd_gaussian_10m_mask[atbd_gaussian_10m_mask < 0] <- 0

path_atbd_test <- '../03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/PRO4SAIL_INVERSION_atbd/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi'
rast_atbd_test <- terra::rast(path_atbd_test)
rast_atbd_test_resampled <- resample(rast_atbd_test, combine)
atbd_test_mask <- mask(rast_atbd_test_resampled, combine)
save_envi_file(atbd_test_mask, "atbd_combine", "../03_RESULTS/Mormal/LAI_PAI_Masks/LAI_S2", "ATBD")
atbd_test_mask <- values(atbd_test_mask)
atbd_test_mask[atbd_test_mask < 0] <- 0

# --------------------------------------------------------------------------------------------------------------------

# ------------------------------------ Comparison: ATBD 8 bands vs ATBD 3 bands --------------------------------------

# Histogram
plot_histogram_1_var(atbd_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     plots_dir, "hist_lai_atbd")
plot_histogram_1_var(atbd_mask3b, "Histogram of ATBD 3 bands Values",
                     "LAI", "ATBD 3 bands", "green",
                     plots_dir, "hist_lai_atbd3b")
plot_histogram_1_var(atbd_mask3b - atbd_mask, "Histogram of Difference between ATBD 3 bands and ATBD Values",
                     "LAI", "Difference ATBDs", "green",
                     plots_dir, "hist_lai_diff_atbd")
plot_histogram_2_vars(atbd_mask3b, atbd_mask, 
                      "Histogram of ATBD 3 bands and ATBD LAI Values", "LAI",
                      "ATBD 3 bands", "ATBD",
                      plots_dir, "hist_lai_atbd3b_atbd")

# Density Scatterplot
plot_density_scatterplot(atbd_mask3b, atbd_mask,
                         "ATBD=f(ATBD 3 bands)", "ATBD 3 bands", "ATBD",
                         plots_dir, "scatterplot_atbd3b_atbd")

# Correlation test
correlation_test_function(atbd_mask3b, atbd_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# ------------------------------------------- Comparison: ATBD vs Raw LiDAR ------------------------------------------

# Histogram
plot_histogram_1_var(atbd_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     plots_dir, "hist_lai_atbd")
plot_histogram_1_var(lidar_mask, "Histogram of LiDAR PAI Values",
                     "PAI", "LiDAR", "green",
                     plots_dir, "hist_lai_lidar")
plot_histogram_2_vars(lidar_mask, atbd_mask, 
                      "Histogram of LiDAR and ATBD LAI Values", "PAI & LAI",
                      "LiDAR", "ATBD",
                      plots_dir, "hist_lai_lidar_atbd")
plot_histogram_2_vars(lidar_mask, atbd_mask3b, 
                      "Histogram of LiDAR and ATBD 3 bands LAI Values", "PAI & LAI",
                      "LiDAR", "ATBD 3 bands",
                      plots_dir, "hist_lai_lidar_atbd3b")

vars_list <- list(atbd_mask, atbd_mask3b, lidar_mask)
var_labs <- c("ATBD 8 bands", "ATBD 3 bands", "LiDAR")

plot_histogram_3_vars(vars_list, "LAI Distributions", "LAI", var_labs, plots_dir, "atbd_atbd3b_lidar")

# Density Scatterplot
plot_density_scatterplot(lidar_mask, atbd_mask,
                         "ATBD=f(LiDAR)", "LiDAR", "ATBD",
                         plots_dir, "scatterplot_lidar_atbd")
plot_density_scatterplot(atbd_mask3b, lidar_mask,
                         "LiDAR=f(ATBD 3 bands)", "ATBD 3 bands", "LiDAR",
                         plots_dir, "scatterplot_atbd3b_lidar")

# Correlation test
correlation_test_function(lidar_mask, atbd_mask, method = "pearson")
correlation_test_function(atbd_mask3b, lidar_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Modif N ------------------------------------------------

atbd_n_high <- resample_and_mask("n", "high", input_base_path, output_directory, single_layer_forest)
atbd_n_low <- resample_and_mask("n", "low", input_base_path, output_directory, single_layer_forest)
atbd_n_full <- resample_and_mask("n", "full", input_base_path, output_directory, single_layer_forest)

vars_list <- list(atbd_mask, atbd_n_high, atbd_n_low, atbd_n_full)
var_labs <- c("ATBD", "ATBD N high", "ATBD N low", "ATBD N full")

plot_histogram_4_vars(vars_list, "LAI Distributions", "LAI", var_labs, plots_dir, "modif_N")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Modif CHL ----------------------------------------------

atbd_chl_high <- resample_and_mask("chl", "high", input_base_path, output_directory, single_layer_forest)
atbd_chl_low <- resample_and_mask("chl", "low", input_base_path, output_directory, single_layer_forest)
atbd_chl_full <- resample_and_mask("chl", "full", input_base_path, output_directory, single_layer_forest)

vars_list <- list(atbd_mask, atbd_chl_high, atbd_chl_low, atbd_chl_full)
var_labs <- c("ATBD", "ATBD CHL high", "ATBD CHL low", "ATBD CHL full")

plot_histogram_4_vars(vars_list, "LAI Distributions", "LAI", var_labs, plots_dir, "modif_CHL")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Modif BROWN --------------------------------------------

atbd_brown_high <- resample_and_mask("brown", "high", input_base_path, output_directory, single_layer_forest)
atbd_brown_low <- resample_and_mask("brown", "low", input_base_path, output_directory, single_layer_forest)
atbd_brown_full <- resample_and_mask("brown", "full", input_base_path, output_directory, single_layer_forest)

vars_list <- list(atbd_mask, atbd_brown_high, atbd_brown_low, atbd_brown_full)
var_labs <- c("ATBD", "ATBD BROWN high", "ATBD BROWN low", "ATBD BROWN full")

plot_histogram_4_vars(vars_list, "LAI Distributions", "LAI", var_labs, plots_dir, "modif_BROWN")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Modif EWT ----------------------------------------------

atbd_ewt_high <- resample_and_mask("ewt", "high", input_base_path, output_directory, single_layer_forest)
atbd_ewt_low <- resample_and_mask("ewt", "low", input_base_path, output_directory, single_layer_forest)
atbd_ewt_full <- resample_and_mask("ewt", "full", input_base_path, output_directory, single_layer_forest)

vars_list <- list(atbd_mask, atbd_ewt_high, atbd_ewt_low, atbd_ewt_full)
var_labs <- c("ATBD", "ATBD EWT high", "ATBD EWT low", "ATBD EWT full")

plot_histogram_4_vars(vars_list, "LAI Distributions", "LAI", var_labs, plots_dir, "modif_EWT")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Modif LMA ----------------------------------------------

atbd_lma_high <- resample_and_mask("lma", "high", input_base_path, output_directory, single_layer_forest)
atbd_lma_low <- resample_and_mask("lma", "low", input_base_path, output_directory, single_layer_forest)
atbd_lma_full <- resample_and_mask("lma", "full", input_base_path, output_directory, single_layer_forest)

vars_list <- list(atbd_mask, atbd_lma_high, atbd_lma_low, atbd_lma_full)
var_labs <- c("ATBD", "ATBD LMA high", "ATBD LMA low", "ATBD LMA full")

plot_histogram_4_vars(vars_list, "LAI Distributions", "LAI", var_labs, plots_dir, "modif_LMA")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Modif LIDFa --------------------------------------------

atbd_lidfa_high <- resample_and_mask("lidfa", "high", input_base_path, output_directory, single_layer_forest)
atbd_lidfa_low <- resample_and_mask("lidfa", "low", input_base_path, output_directory, single_layer_forest)
atbd_lidfa_full <- resample_and_mask("lidfa", "full", input_base_path, output_directory, single_layer_forest)

vars_list <- list(atbd_mask, atbd_lidfa_high, atbd_lidfa_low, atbd_lidfa_full)
var_labs <- c("ATBD", "ATBD LIDFa high", "ATBD LIDFa low", "ATBD LIDFa full")

plot_histogram_4_vars(vars_list, "LAI Distributions", "LAI", var_labs, plots_dir, "modif_LIDFa")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Modif lai ----------------------------------------------

atbd_lai_high <- resample_and_mask("lai", "high", input_base_path, output_directory, single_layer_forest)
atbd_lai_low <- resample_and_mask("lai", "low", input_base_path, output_directory, single_layer_forest)
atbd_lai_full <- resample_and_mask("lai", "full", input_base_path, output_directory, single_layer_forest)

vars_list <- list(atbd_mask, atbd_lai_high, atbd_lai_low, atbd_lai_full)
var_labs <- c("ATBD", "ATBD lai high", "ATBD lai low", "ATBD lai full")

plot_histogram_4_vars(vars_list, "LAI Distributions", "LAI", var_labs, plots_dir, "modif_lai")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Modif q ------------------------------------------------

atbd_q_high <- resample_and_mask("q", "high", input_base_path, output_directory, single_layer_forest)
atbd_q_low <- resample_and_mask("q", "low", input_base_path, output_directory, single_layer_forest)
atbd_q_full <- resample_and_mask("q", "full", input_base_path, output_directory, single_layer_forest)

vars_list <- list(atbd_mask, atbd_q_high, atbd_q_low, atbd_q_full)
var_labs <- c("ATBD", "ATBD q high", "ATBD q low", "ATBD q full")

plot_histogram_4_vars(vars_list, "LAI Distributions", "LAI", var_labs, plots_dir, "modif_q")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Modif psoil --------------------------------------------

atbd_psoil_high <- resample_and_mask("psoil", "high", input_base_path, output_directory, single_layer_forest)
atbd_psoil_low <- resample_and_mask("psoil", "low", input_base_path, output_directory, single_layer_forest)
atbd_psoil_full <- resample_and_mask("psoil", "full", input_base_path, output_directory, single_layer_forest)

vars_list <- list(atbd_mask, atbd_psoil_high, atbd_psoil_low, atbd_psoil_full)
var_labs <- c("ATBD", "ATBD psoil high", "ATBD psoil low", "ATBD psoil full")

plot_histogram_4_vars(vars_list, "LAI Distributions", "LAI", var_labs, plots_dir, "modif_psoil")

# --------------------------------------------------------------------------------------------------------------------