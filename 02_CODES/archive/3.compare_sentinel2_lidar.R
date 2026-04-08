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
lidar_output_dir <- file.path(results_path, "LAI_PAI_Masks/PAI_LiDAR/lidR")
output_directory <- file.path(results_path, "LAI_PAI_Masks/LAI_S2")

# --------------------------------------------------------------------------------------------------------------------

# ------------------------------------- Preparation: Get LiDAR & Sentinel-2 Mask Values ------------------------------

# LiDAR
lidar_mask <- open_envi_file_as_values(file.path(lidar_output_dir, "lidar.envi"))

# Sentinel-2
atbd_10m_mask <- open_envi_file_as_values(file.path(output_directory, 
                                                    "atbd/atbd_addmult_10m.envi"))

atbd_20m_mask <- open_envi_file_as_values(file.path(output_directory, 
                                                    "atbd/atbd_addmult_20m.envi"))

param_list <- list(lidar_mask,
                   atbd_10m_mask,
                   atbd_20m_mask)

param_labs <- c("LiDAR", "ATBD 10m", "ATBD 20m")  

plot_histogram_1_var(atbd_10m_mask, "Histogram of ATBD 10m LAI Values",
                     "LAI", "ATBD 10m", "green",
                     plots_dir, "hist_lai_atbd_10m")

plot_histogram_1_var(atbd_20m_mask, "Histogram of ATBD 20m LAI Values",
                     "LAI", "ATBD 20m", "green",
                     plots_dir, "hist_lai_atbd_20m")

plot_histogram_1_var(lidar_mask, "Histogram of LiDAR PAI Values",
                     "PAI", "LiDAR", "green",
                     plots_dir, "hist_lai_lidar")

plot_histogram_3_vars(param_list, "LAI Distributions", "LAI", 
                      param_labs, plots_dir, "atbd_10m_20m_lidar_comp")

# Density Scatterplot
plot_density_scatterplot(lidar_mask, atbd_10m_mask,
                         "ATBD 10m=f(LiDAR)", "LiDAR", "ATBD 10m",
                         plots_dir, "scatterplot_lidar_atbd10m")

plot_density_scatterplot(lidar_mask, atbd_20m_mask,
                         "ATBD 20m=f(LiDAR)", "LiDAR", "ATBD 20m",
                         plots_dir, "scatterplot_lidar_atbd20m")

plot_density_scatterplot(atbd_10m_mask, atbd_20m_mask,
                         "ATBD 20m=f(ATBD 10m)", "ATBD 10m", "ATBD 20m",
                         plots_dir, "scatterplot_atbd10m_atbd20m")

# Correlation test
correlation_test_function(lidar_mask, atbd_10m_mask, method = "pearson")
correlation_test_function(lidar_mask, atbd_20m_mask, method = "pearson")
correlation_test_function(atbd_10m_mask, atbd_20m_mask, method = "pearson")


# --------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------- LAI Comparison --------------------------------------------------

# -------------------------------------- Comparison: ATBD 10m vs 20m vs Raw LiDAR ------------------------------------

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

# Density Scatterplot
plot_density_scatterplot(lidar_mask, atbd_mask,
                         "ATBD=f(LiDAR)", "LiDAR", "ATBD",
                         plots_dir, "scatterplot_lidar_atbd")

# Correlation test
correlation_test_function(lidar_mask, atbd_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# -------------------------------------- Comparison: ATBD 10m vs 20m vs AddMult vs Gauss -----------------------------

atbd_addmult_10m_mask <- open_envi_file_as_values(file.path(output_directory, "atbd/atbd_addmult_10m.envi"))
atbd_addmult_20m_mask <- open_envi_file_as_values(file.path(output_directory, "atbd/atbd_addmult_20m.envi"))
atbd_gaussian_10m_mask <- open_envi_file_as_values(file.path(output_directory, "atbd/atbd_gaussian_10m.envi"))
atbd_gaussian_20m_mask <- open_envi_file_as_values(file.path(output_directory, "atbd/atbd_gaussian_20m.envi"))

# Histogram
plot_histogram_1_var(atbd_addmult_10m_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     plots_dir, "hist_lai_atbd_10m_admd")
plot_histogram_1_var(atbd_addmult_20m_mask, "Histogram of LiDAR PAI Values",
                     "PAI", "LiDAR", "green",
                     plots_dir, "hist_lai_atbd_20m_admd")
plot_histogram_1_var(atbd_gaussian_10m_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     plots_dir, "hist_lai_atbd_10m_gauss")
plot_histogram_1_var(atbd_gaussian_20m_mask, "Histogram of LiDAR PAI Values",
                     "PAI", "LiDAR", "green",
                     plots_dir, "hist_lai_atbd_20m_gauss")
plot_histogram_2_vars(atbd_addmult_10m_mask, atbd_gaussian_10m_mask, 
                      "Histogram of LiDAR and ATBD LAI Values", "PAI & LAI",
                      "ATBD ADMD", "ATBD Gauss",
                      plots_dir, "hist_lai_atbd_10m")
plot_histogram_2_vars(atbd_addmult_20m_mask, atbd_gaussian_20m_mask, 
                      "Histogram of LiDAR and ATBD LAI Values", "PAI & LAI",
                      "ATBD ADMD", "ATBD Gauss",
                      plots_dir, "hist_lai_atbd_20m")
plot_histogram_4_vars(list(atbd_addmult_10m_mask, atbd_addmult_20m_mask,
                           atbd_gaussian_10m_mask, atbd_gaussian_20m_mask), 
                      "Histogram of LiDAR and ATBD LAI Values", "PAI & LAI",
                      c("ATBD ADMD 10m", "ATBD ADMD 20m", 
                        "ATBD Gaussian 10m", "ATBD Gaussian 20m"),
                      plots_dir, "hist_lai_atbd_bands_noise")

plot_density_scatterplot(atbd_addmult_20m_mask, lidar_mask,
                         "ATBD=f(LiDAR)", "LiDAR", "ATBD",
                         plots_dir, "scatterplot_20m_lidar")

# Density Scatterplot
plot_density_scatterplot(lidar_mask, atbd_mask,
                         "ATBD=f(LiDAR)", "LiDAR", "ATBD",
                         plots_dir, "scatterplot_lidar_atbd")

# Correlation test
correlation_test_function(lidar_mask, atbd_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# ------------------------------------- Comparison: Masked LiDAR vs Raw LiDAR ----------------------------------------

# Histogram
plot_histogram_1_var(lidar_mask_masked2m, "Histogram of ATBD LAI Values",
                     "LAI", "LiDAR Masked2m", "green",
                     plots_dir, "hist_lai_lidarmasked2m")
plot_histogram_1_var(lidar_mask, "Histogram of LiDAR PAI Values",
                     "LAI", "LiDAR", "green",
                     plots_dir, "hist_lai_lidar")
plot_histogram_2_vars(lidar_mask, lidar_mask_masked2m, 
                      "Histogram of LiDAR and LiDAR Masked2m LAI Values", "LAI",
                      "LiDAR", "LiDAR Masked2m",
                      plots_dir, "hist_lai_lidar_lidarmasked2m")

# Density Scatterplot
plot_density_scatterplot(lidar_mask, lidar_mask_masked2m,
                         "LiDAR Masked2m=f(LiDAR)", "LiDAR", "LiDAR Masked2m",
                         plots_dir, "scatterplot_lidar_lidarmasked2m")

# Correlation test
correlation_test_function(lidar_mask, lidar_mask_masked2m, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# -------------------------------------- Comparison: ATBD vs Q.Zhang et al. 2005 -------------------------------------

# Histogram
plot_histogram_1_var(zhang2005_mask, "Histogram of Zhang2005 LAI Values",
                     "LAI", "Zhang2005", "green",
                     plots_dir, "hist_lai_zhang2005")
plot_histogram_1_var(atbd_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     plots_dir, "zhang2005_mask")
plot_histogram_2_vars(atbd_mask, zhang2005_mask, 
                      "Histogram of Zhang2005 and ATBD LAI Values", "LAI",
                      "ATBD", "Zhang2005",
                      plots_dir, "hist_lai_zhang2005_atbd")

# Density Scatterplot
plot_density_scatterplot(zhang2005_mask, atbd_mask,
                         "ATBD=f(Zhang2005)", "Zhang2005", "ATBD",
                         plots_dir, "scatterplot_atbd_zhang2005")

# Correlation test
correlation_test_function(zhang2005_mask, atbd_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Brede et al. 2020 --------------------------------------

# Histogram
plot_histogram_1_var(brede2020_mask, "Histogram of Brede2020 LAI Values",
                     "LAI", "Brede2020", "green",
                     plots_dir, "hist_lai_brede2020")
plot_histogram_1_var(atbd_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     plots_dir, "hist_lai_atbd")
plot_histogram_2_vars(atbd_mask, brede2020_mask, 
                      "Histogram of Brede2020 and ATBD LAI Values", "LAI",
                      "ATBD", "Brede2020",
                      plots_dir, "hist_lai_brede2020_atbd")

# Density Scatterplot
plot_density_scatterplot(brede2020_mask, atbd_mask,
                         "ATBD=f(Brede2020)", "Brede2020", "ATBD",
                         plots_dir, "scatterplot_atbd_brede2020")

# Correlation test
correlation_test_function(brede2020_mask, atbd_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Hauser et al. 2021 --------------------------------------

# Histogram
plot_histogram_1_var(hauser2021_mask, "Histogram of Hauser2021 LAI Values",
                     "LAI", "Hauser2021", "green",
                     plots_dir, "hist_lai_hauser2021")
plot_histogram_1_var(atbd_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     plots_dir, "hist_lai_atbd")
plot_histogram_2_vars(atbd_mask, hauser2021_mask, 
                      "Histogram of Hauser2021 and ATBD LAI Values", "LAI",
                      "ATBD", "Hauser2021",
                      plots_dir, "hist_lai_hauser2021_atbd")

# Density Scatterplot
plot_density_scatterplot(hauser2021_mask, atbd_mask,
                         "ATBD=f(Hauser2021)", "Hauser2021", "ATBD",
                         plots_dir, "scatterplot_atbd_hauser2021")

# Correlation test
correlation_test_function(hauser2021_mask, atbd_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Sinha et al. 2020 --------------------------------------

# Histogram
plot_histogram_1_var(sinha2020_mask, "Histogram of Sinha2020 LAI Values",
                     "LAI", "Sinha2020", "green",
                     plots_dir, "hist_lai_sinha2020")
plot_histogram_1_var(atbd_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     plots_dir, "hist_lai_atbd")
plot_histogram_2_vars(atbd_mask, sinha2020_mask, 
                      "Histogram of Hauser2021 and ATBD LAI Values", "LAI",
                      "ATBD", "sinha2020",
                      plots_dir, "hist_lai_sinha2020_atbd")

# Density Scatterplot
plot_density_scatterplot(sinha2020_mask, atbd_mask,
                         "ATBD=f(Sinha2020)", "Sinha2020", "ATBD",
                         plots_dir, "scatterplot_atbd_sinha2020")

# Correlation test
correlation_test_function(sinha2020_mask, atbd_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Verhoef and Bach 2007-----------------------------------

# Histogram
plot_histogram_1_var(verhoef2007_mask, "Histogram of Verhoef2007 LAI Values",
                     "LAI", "Verhoef2007", "green",
                     plots_dir, "hist_lai_verhoef2007")
plot_histogram_1_var(lidar_mask - atbd_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     plots_dir, "hist_lai_atbd")
plot_histogram_2_vars(atbd_mask, verhoef2007_mask, 
                      "Histogram of Verhoef2007 and ATBD LAI Values", "LAI",
                      "ATBD", "verhoef2007",
                      plots_dir, "hist_lai_verhoef2007_atbd")

# Density Scatterplot
plot_density_scatterplot(verhoef2007_mask, atbd_mask,
                         "ATBD=f(Verhoef2007)", "Verhoef2007", "ATBD",
                         plots_dir, "scatterplot_atbd_verhoef2007")

# Correlation test
correlation_test_function(verhoef2007_mask, atbd_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Shiklomanov et al. 2016---------------------------------

# Histogram
plot_histogram_1_var(shiklomanov2016_mask, "Histogram of Shiklomanov2016 LAI Values",
                     "LAI", "Shiklomanov2016", "green",
                     plots_dir, "hist_lai_shiklomanov2016")
plot_histogram_1_var(atbd_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     plots_dir, "hist_lai_atbd")
plot_histogram_2_vars(atbd_mask, shiklomanov2016_mask, 
                      "Histogram of Shiklomanov2016 and ATBD LAI Values", "LAI",
                      "ATBD", "shiklomanov2016",
                      plots_dir, "hist_lai_shiklomanov2016_atbd")

# Density Scatterplot
plot_density_scatterplot(shiklomanov2016_mask, atbd_mask,
                         "ATBD=f(Shiklomanov2016)", "Shiklomanov2016", "ATBD",
                         plots_dir, "scatterplot_atbd_shiklomanov2016")

# Correlation test
correlation_test_function(shiklomanov2016_mask, atbd_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# ------------------------------ Comparison: Q.Zhang et al. 2005 vs Brede et al. 2020 --------------------------------

# Histogram
plot_histogram_1_var(brede2020_mask, "Histogram of Brede2020 LAI Values",
                     "LAI", "Zhang2005", "green",
                     plots_dir, "hist_lai_brede2020")
plot_histogram_1_var(zhang2005_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     plots_dir, "hist_lai_zhang2005")
plot_histogram_2_vars(zhang2005_mask, brede2020_mask, 
                      "Histogram of Zhang2005 and Brede2020 LAI Values", "LAI",
                      "Zhang2005", "Brede2020",
                      plots_dir, "hist_lai_zhang2005_brede2020")

# Density Scatterplot
plot_density_scatterplot(zhang2005_mask, brede2020_mask,
                         "ATBD=f(Brede2020)", "Brede2020", "ATBD",
                         plots_dir, "scatterplot_zhang2005_brede2020")

# Correlation test
correlation_test_function(zhang2005_mask, brede2020_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# ------------------------------ Comparison: Q.Zhang et al. 2005 vs Verhoef and Bach 2007 -----------------------------

# Histogram
plot_histogram_1_var(verhoef2007_mask, "Histogram of Verhoef2007 LAI Values",
                     "LAI", "Zhang2005", "green",
                     plots_dir, "hist_lai_verhoef2007")
plot_histogram_1_var(zhang2005_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     plots_dir, "hist_lai_zhang2005")
plot_histogram_2_vars(zhang2005_mask, verhoef2007_mask, 
                      "Histogram of Zhang2005 and Verhoef2007 LAI Values", "LAI",
                      "Zhang2005", "Verhoef2007",
                      plots_dir, "hist_lai_zhang2005_verhoef2007")

# Density Scatterplot
plot_density_scatterplot(zhang2005_mask, verhoef2007_mask,
                         "Zhang2005=f(Verhoef2007)", "Verhoef2007", "ATBD",
                         plots_dir, "scatterplot_zhang2005_verhoef2007")

# Correlation test
correlation_test_function(zhang2005_mask, verhoef2007_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

vars_list <- list(atbd_mask, lidar_mask, zhang2005_mask, brede2020_mask,
                  hauser2021_mask, sinha2020_mask,
                  verhoef2007_mask, shiklomanov2016_mask)

titles <- c("LAI Distributions")
xlab <- "LAI"
var_labs <- c("ATBD", "LiDAR", "Zhang2005", "Brede2020", "Hauser2021", "Sinha2020", 
              "Verhoef2007", "Shiklomanov2016")

plot_histogram_8_vars(vars_list, titles, xlab, var_labs, plots_dir, "8vars")

# --------------------------------------- Comparison: ATBD vs Modif N ------------------------------------------------

vars_list <- list(atbd_mask, atbd_n_high_mask, atbd_n_low_mask)
title <- "LAI Distributions"
xlab <- "LAI"
var_labs <- c("ATBD", "ATBD N high", "ATBD N low")

plot_histogram_3_vars(vars_list, title, xlab, var_labs, plots_dir, "Modif N")

# --------------------------------------------------------------------------------------------------------------------

# ------------------------------------------- Comparison: ATBD vs LiDAR Het ------------------------------------------

correlation_values <- c()
inc <- 0.1
for (low_quantile in seq(0.0, 0.99, by = inc)) {
  quantiles <- terra::rast(file.path(deciles_dir, 
                                     paste0("pai_masked_", 
                                            low_quantile, 
                                            "_", 
                                            low_quantile + inc,
                                            ".tif")))
  quantiles <- values(quantiles)
  correlation_value <- correlation_test_function(atbd_10m_mask, quantiles, method = "pearson")
  correlation_values <- c(correlation_values, correlation_value)
}

midpoints <- seq(0.05, 0.95, by = 0.1)
plot(1, type = "n", 
     xlab = "Decile", ylab = "Correlation Value", 
     xlim = c(0, 1), ylim = range(correlation_values),
     main = "Correlation Values for Deciles")
lines(midpoints, correlation_values, type = "l")
points(midpoints, correlation_values, pch = 16, col = "red")
axis(side = 1, at = seq(0.1, 0.9, by = 0.2), labels = seq(0.1, 0.9, by = 0.2))