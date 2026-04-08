# ---
# title: "3.compare_s2_distribs.R"
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
dirname <- "../03_RESULTS/Mormal/Plots"
results_path <- "../03_RESULTS/Mormal"
deciles_dir <- file.path(results_path, "/LiDAR_Masks/Deciles")

# --------------------------------------------------------------------------------------------------------------------

# ------------------------------------- Preparation: Get LiDAR & Sentinel-2 Mask Values ------------------------------

# LiDAR
pai_lidar_path <- "../03_RESULTS/Mormal/LAI_PAI_Masks/PAI_LiDAR/fdeboissieu/pai_las_cleaned.tif"
pai_lidar <- terra::rast(pai_lidar_path)

pai_lidar_resampled <- resample(pai_lidar, single_layer_forest)
lidar_mask <- mask(pai_lidar_resampled, single_layer_forest)

save_envi_file(lidar_mask, "lidar", "../03_RESULTS/Mormal/LAI_PAI_Masks/PAI_LiDAR/fdeboissieu", "LiDAR")

lidar_mask <- values(lidar_mask)
lidar_mask[lidar_mask < 0] <- 0

# Sentinel-2
path_atbd <- '../03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/PRO4SAIL_INVERSION_atbd/codist_LAI_TRUE/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi'
path_zhang2005 <- '../03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/PRO4SAIL_INVERSION_q_zhang_et_al_2005/codist_LAI_TRUE/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi'
path_brede2020 <- '../03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/PRO4SAIL_INVERSION_brede_et_al_2020/codist_LAI_TRUE/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi'
path_hauser2021 <- '../03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/PRO4SAIL_INVERSION_hauser_et_al_2021/codist_LAI_TRUE/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi'
path_sinha2020 <- '../03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/PRO4SAIL_INVERSION_sinha_et_al_2020/codist_LAI_TRUE/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi'
path_verhoef2007 <- '../03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/PRO4SAIL_INVERSION_verhoef_and_bach_2007/codist_LAI_TRUE/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi'
path_shiklomanov2016 <- '../03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/PRO4SAIL_INVERSION_shiklomanov_et_al_2016/codist_LAI_TRUE/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi'

# Raster
rast_atbd <- terra::rast(path_atbd)
rast_zhang2005 <- terra::rast(path_zhang2005)
rast_brede2020 <- terra::rast(path_brede2020)
rast_hauser2021 <- terra::rast(path_hauser2021)
rast_sinha2020 <- terra::rast(path_sinha2020)
rast_verhoef2007 <- terra::rast(path_verhoef2007)
rast_shiklomanov2016 <- terra::rast(path_shiklomanov2016)

# Resample to match dimensions
rast_atbd_resampled <- resample(rast_atbd, single_layer_forest)
rast_zhang2005_resampled <- resample(rast_zhang2005, single_layer_forest)
rast_brede2020_resampled <- resample(rast_brede2020, single_layer_forest)
rast_hauser2021_resampled <- resample(rast_hauser2021, single_layer_forest)
rast_sinha2020_resampled <- resample(rast_sinha2020, single_layer_forest)
rast_verhoef2007_resampled <- resample(rast_verhoef2007, single_layer_forest)
rast_shiklomanov2016_resampled <- resample(rast_shiklomanov2016, single_layer_forest)

# Mask
atbd_mask <- mask(rast_atbd_resampled, single_layer_forest)
zhang2005_mask <- mask(rast_zhang2005_resampled, single_layer_forest)
brede2020_mask <- mask(rast_brede2020_resampled, single_layer_forest)
hauser2021_mask <- mask(rast_hauser2021_resampled, single_layer_forest)
sinha2020_mask <- mask(rast_sinha2020_resampled, single_layer_forest)
verhoef2007_mask <- mask(rast_verhoef2007_resampled, single_layer_forest)
shiklomanov2016_mask <- mask(rast_shiklomanov2016_resampled, single_layer_forest)

# Write ENVI files
save_envi_file(atbd_mask, "atbd", "../03_RESULTS/Mormal/LAI_PAI_Masks/LAI_S2", "ATBD")
save_envi_file(zhang2005_mask, "zhang2005", "../03_RESULTS/Mormal/LAI_PAI_Masks/LAI_S2", "zhang2005_mask")
save_envi_file(brede2020_mask, "brede2020", "../03_RESULTS/Mormal/LAI_PAI_Masks/LAI_S2", "brede2020_mask")
save_envi_file(hauser2021_mask, "hauser2021", "../03_RESULTS/Mormal/LAI_PAI_Masks/LAI_S2", "hauser2021_mask")
save_envi_file(sinha2020_mask, "sinha2020", "../03_RESULTS/Mormal/LAI_PAI_Masks/LAI_S2", "sinha2020_mask")
save_envi_file(verhoef2007_mask, "verhoef2007", "../03_RESULTS/Mormal/LAI_PAI_Masks/LAI_S2", "verhoef2007_mask")
save_envi_file(shiklomanov2016_mask, "shiklomanov2016", "../03_RESULTS/Mormal/LAI_PAI_Masks/LAI_S2", "shiklomanov2016_mask")

# Keep only values
atbd_mask <- values(atbd_mask)
zhang2005_mask <- values(zhang2005_mask)
brede2020_mask <- values(brede2020_mask)
hauser2021_mask <- values(hauser2021_mask)
sinha2020_mask <- values(sinha2020_mask)
verhoef2007_mask <- values(verhoef2007_mask)
shiklomanov2016_mask <- values(shiklomanov2016_mask)

# Truncate negative values
atbd_mask[atbd_mask < 0] <- 0
zhang2005_mask[zhang2005_mask < 0] <- 0
brede2020_mask[brede2020_mask < 0] <- 0
hauser2021_mask[hauser2021_mask < 0] <- 0
sinha2020_mask[sinha2020_mask] <- 0
verhoef2007_mask[verhoef2007_mask < 0] <- 0
shiklomanov2016_mask[shiklomanov2016_mask < 0] <- 0

# --------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------- LAI Comparison --------------------------------------------------

# ------------------------------------------- Comparison: ATBD vs Raw LiDAR ------------------------------------------

# Histogram
plot_histogram_1_var(atbd_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     dirname, "hist_lai_atbd")
plot_histogram_1_var(lidar_mask, "Histogram of LiDAR PAI Values",
                     "PAI", "LiDAR", "green",
                     dirname, "hist_lai_lidar")
plot_histogram_2_vars(lidar_mask, atbd_mask, 
                      "Histogram of LiDAR and ATBD LAI Values", "PAI & LAI",
                      "LiDAR", "ATBD",
                      dirname, "hist_lai_lidar_atbd")

# Density Scatterplot
plot_density_scatterplot(lidar_mask, atbd_mask,
                         "ATBD=f(LiDAR)", "LiDAR", "ATBD",
                         dirname, "scatterplot_lidar_atbd")

# Correlation test
correlation_test_function(lidar_mask, atbd_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# -------------------------------------- Comparison: ATBD vs Q.Zhang et al. 2005 -------------------------------------

# Histogram
plot_histogram_1_var(zhang2005_mask, "Histogram of Zhang2005 LAI Values",
                     "LAI", "Zhang2005", "green",
                     dirname, "hist_lai_zhang2005")
plot_histogram_1_var(atbd_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     dirname, "zhang2005_mask")
plot_histogram_2_vars(atbd_mask, zhang2005_mask, 
                      "Histogram of Zhang2005 and ATBD LAI Values", "LAI",
                      "ATBD", "Zhang2005",
                      dirname, "hist_lai_zhang2005_atbd")

# Density Scatterplot
plot_density_scatterplot(zhang2005_mask, atbd_mask,
                         "ATBD=f(Zhang2005)", "Zhang2005", "ATBD",
                         dirname, "scatterplot_atbd_zhang2005")

# Correlation test
correlation_test_function(zhang2005_mask, atbd_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Brede et al. 2020 --------------------------------------

# Histogram
plot_histogram_1_var(brede2020_mask, "Histogram of Brede2020 LAI Values",
                     "LAI", "Brede2020", "green",
                     dirname, "hist_lai_brede2020")
plot_histogram_1_var(atbd_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     dirname, "hist_lai_atbd")
plot_histogram_2_vars(atbd_mask, brede2020_mask, 
                      "Histogram of Brede2020 and ATBD LAI Values", "LAI",
                      "ATBD", "Brede2020",
                      dirname, "hist_lai_brede2020_atbd")

# Density Scatterplot
plot_density_scatterplot(brede2020_mask, atbd_mask,
                         "ATBD=f(Brede2020)", "Brede2020", "ATBD",
                         dirname, "scatterplot_atbd_brede2020")

# Correlation test
correlation_test_function(brede2020_mask, atbd_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Hauser et al. 2021 --------------------------------------

# Histogram
plot_histogram_1_var(hauser2021_mask, "Histogram of Hauser2021 LAI Values",
                     "LAI", "Hauser2021", "green",
                     dirname, "hist_lai_hauser2021")
plot_histogram_1_var(atbd_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     dirname, "hist_lai_atbd")
plot_histogram_2_vars(atbd_mask, hauser2021_mask, 
                      "Histogram of Hauser2021 and ATBD LAI Values", "LAI",
                      "ATBD", "Hauser2021",
                      dirname, "hist_lai_hauser2021_atbd")

# Density Scatterplot
plot_density_scatterplot(hauser2021_mask, atbd_mask,
                         "ATBD=f(Hauser2021)", "Hauser2021", "ATBD",
                         dirname, "scatterplot_atbd_hauser2021")

# Correlation test
correlation_test_function(hauser2021_mask, atbd_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Sinha et al. 2020 --------------------------------------

# Histogram
plot_histogram_1_var(sinha2020_mask, "Histogram of Sinha2020 LAI Values",
                     "LAI", "Sinha2020", "green",
                     dirname, "hist_lai_sinha2020")
plot_histogram_1_var(atbd_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     dirname, "hist_lai_atbd")
plot_histogram_2_vars(atbd_mask, sinha2020_mask, 
                      "Histogram of Hauser2021 and ATBD LAI Values", "LAI",
                      "ATBD", "sinha2020",
                      dirname, "hist_lai_sinha2020_atbd")

# Density Scatterplot
plot_density_scatterplot(sinha2020_mask, atbd_mask,
                         "ATBD=f(Sinha2020)", "Sinha2020", "ATBD",
                         dirname, "scatterplot_atbd_sinha2020")

# Correlation test
correlation_test_function(sinha2020_mask, atbd_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Verhoef and Bach 2007-----------------------------------

# Histogram
plot_histogram_1_var(verhoef2007_mask, "Histogram of Verhoef2007 LAI Values",
                     "LAI", "Verhoef2007", "green",
                     dirname, "hist_lai_verhoef2007")
plot_histogram_1_var(lidar_mask - atbd_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     dirname, "hist_lai_atbd")
plot_histogram_2_vars(atbd_mask, verhoef2007_mask, 
                      "Histogram of Verhoef2007 and ATBD LAI Values", "LAI",
                      "ATBD", "verhoef2007",
                      dirname, "hist_lai_verhoef2007_atbd")

# Density Scatterplot
plot_density_scatterplot(verhoef2007_mask, atbd_mask,
                         "ATBD=f(Verhoef2007)", "Verhoef2007", "ATBD",
                         dirname, "scatterplot_atbd_verhoef2007")

# Correlation test
correlation_test_function(verhoef2007_mask, atbd_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------- Comparison: ATBD vs Shiklomanov et al. 2016---------------------------------

# Histogram
plot_histogram_1_var(shiklomanov2016_mask, "Histogram of Shiklomanov2016 LAI Values",
                     "LAI", "Shiklomanov2016", "green",
                     dirname, "hist_lai_shiklomanov2016")
plot_histogram_1_var(atbd_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     dirname, "hist_lai_atbd")
plot_histogram_2_vars(atbd_mask, shiklomanov2016_mask, 
                      "Histogram of Shiklomanov2016 and ATBD LAI Values", "LAI",
                      "ATBD", "shiklomanov2016",
                      dirname, "hist_lai_shiklomanov2016_atbd")

# Density Scatterplot
plot_density_scatterplot(shiklomanov2016_mask, atbd_mask,
                         "ATBD=f(Shiklomanov2016)", "Shiklomanov2016", "ATBD",
                         dirname, "scatterplot_atbd_shiklomanov2016")

# Correlation test
correlation_test_function(shiklomanov2016_mask, atbd_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# ------------------------------ Comparison: Q.Zhang et al. 2005 vs Brede et al. 2020 --------------------------------

# Histogram
plot_histogram_1_var(brede2020_mask, "Histogram of Brede2020 LAI Values",
                     "LAI", "Zhang2005", "green",
                     dirname, "hist_lai_brede2020")
plot_histogram_1_var(zhang2005_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     dirname, "hist_lai_zhang2005")
plot_histogram_2_vars(zhang2005_mask, brede2020_mask, 
                      "Histogram of Zhang2005 and Brede2020 LAI Values", "LAI",
                      "Zhang2005", "Brede2020",
                      dirname, "hist_lai_zhang2005_brede2020")

# Density Scatterplot
plot_density_scatterplot(zhang2005_mask, brede2020_mask,
                         "ATBD=f(Brede2020)", "Brede2020", "ATBD",
                         dirname, "scatterplot_zhang2005_brede2020")

# Correlation test
correlation_test_function(zhang2005_mask, brede2020_mask, method = "pearson")

# --------------------------------------------------------------------------------------------------------------------

# ------------------------------ Comparison: Q.Zhang et al. 2005 vs Verhoef and Bach 2007 -----------------------------

# Histogram
plot_histogram_1_var(verhoef2007_mask, "Histogram of Verhoef2007 LAI Values",
                     "LAI", "Zhang2005", "green",
                     dirname, "hist_lai_verhoef2007")
plot_histogram_1_var(zhang2005_mask, "Histogram of ATBD LAI Values",
                     "LAI", "ATBD", "green",
                     dirname, "hist_lai_zhang2005")
plot_histogram_2_vars(zhang2005_mask, verhoef2007_mask, 
                      "Histogram of Zhang2005 and Verhoef2007 LAI Values", "LAI",
                      "Zhang2005", "Verhoef2007",
                      dirname, "hist_lai_zhang2005_verhoef2007")

# Density Scatterplot
plot_density_scatterplot(zhang2005_mask, verhoef2007_mask,
                         "Zhang2005=f(Verhoef2007)", "Verhoef2007", "ATBD",
                         dirname, "scatterplot_zhang2005_verhoef2007")

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

plot_histogram_8_vars(vars_list, titles, xlab, var_labs, dirname, "8varsold")

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
  correlation_value <- correlation_test_function(brede2020_mask, quantiles, method = "pearson")
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