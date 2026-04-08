# ---
# title: "0.verif_atbd_biophysical_toolbox.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-02-20"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# Import useful functions
source("../libraries/functions_sentinel_2.R")
source("../libraries/functions_plots.R")
source("../libraries/functions_intercomparisons.R")
library(Metrics)

bands_choice_10m <- "10m"
bands_choice_20m <- "bands_3_4_5_6_7_8A_11_12"

biophysical_10m <- "_s2resampled_10m_biophysical10m.data"
biophysical_20m <- "_s2resampled_20m_biophysical.data"

biophysical_10m_20m <- "_s2resampled_10m_biophysical20m.data"

site <- "Mormal"
prosail_dir <- file.path("../../03_RESULTS", site) 
snap_dir <- "../../../Softwares/SNAP_Toolbox_Results"
snap_image <- "S2A_MSIL2A_20210614T105031_N0300_R051_T31UER_20210614T140120"
plots_dir <- file.path(prosail_dir, "Plots/verif_toolbox")

resample_method <- "bilinear" #nearest #bilinear

mask_10m <- terra::rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/LiDAR/Heterogeneity_Masks/artifacts_low_vegetation_majority_90_p_res_10_m.envi")
plot(mask_10m)

mask_20m <- terra::rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/LiDAR/Heterogeneity_Masks/artifacts_low_vegetation_majority_90_p_res_20_m.envi")
plot(mask_20m)

# 10m Img 10m Bands
# atbd_10m_prosail <- terra::rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/res_10m/PRO4SAIL_INVERSION_atbd/addmult/bands_3_4_8/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi")
# plot(atbd_10m_prosail)
# atbd_10m_prosail <- atbd_10m_prosail * mask_10m
# plot(atbd_10m_prosail)

# 20m Img 10m Bands
atbd_10m_prosail <- terra::rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/res_20m/PRO4SAIL_INVERSION_atbd/addmult/bands_3_4_8/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi")
plot(atbd_10m_prosail)
atbd_10m_prosail <- atbd_10m_prosail * mask_20m
plot(atbd_10m_prosail)

# 20m Img 20m Bands
atbd_20m_prosail <- terra::rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/res_20m/PRO4SAIL_INVERSION_atbd/addmult/bands_3_4_5_6_7_8A_11_12/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi")
atbd_20m_prosail <- atbd_20m_prosail * mask_20m
plot(atbd_20m_prosail)

# 10m Img 20m Bands
# atbd_1020m_prosail <- terra::rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/res_10m/PRO4SAIL_INVERSION_modifbandatbd/addmult/bands_3_4_5_6_7_8A_11_12/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi")
# atbd_1020m_prosail <- terra::project(atbd_1020m_prosail, mask_20m)
# atbd_1020m_prosail <- atbd_1020m_prosail * mask_20m
# atbd_1020m_prosail_values <- values(atbd_1020m_prosail)

atbd_10m_snap <- terra::rast(file.path(snap_dir,
                                       resample_method,
                                       paste0(snap_image, 
                                              biophysical_10m), 
                                       "lai.img"
))

# 10m Img 10m Bands
# atbd_10m_snap <- terra::project(atbd_10m_snap, mask_10m)
# atbd_10m_snap <- atbd_10m_snap * mask_10m

# 20m Img 10m Bands
atbd_10m_snap <- terra::project(atbd_10m_snap, mask_20m)
atbd_10m_snap <- atbd_10m_snap * mask_20m

# 20m Img 20m Bands
atbd_20m_snap <- terra::rast(file.path(snap_dir,
                                       resample_method,
                                       paste0(snap_image, 
                                              biophysical_20m),
                                       "lai.img"
))

atbd_20m_snap <- terra::project(atbd_20m_snap, mask_20m)
atbd_20m_snap <- atbd_20m_snap * mask_20m

# Snap 10m LAI BP 20m
atbd_10m20m_snap <- terra::rast(file.path(snap_dir,
                                          resample_method,
                                          paste0(snap_image, 
                                                 biophysical_10m_20m),
                                          "lai.img"
))

atbd_10m20m_snap <- terra::project(atbd_10m20m_snap, mask_20m)
atbd_10m20m_snap <- atbd_10m20m_snap * mask_20m

# LiDAR
lidar_lai <- terra::rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/Metrics/Deciduous_Only/lidarlai_res_10_m.tif")
lidar_lai <- terra::project(lidar_lai, mask_20m)
lidar_lai <- lidar_lai * mask_20m

# Values

atbd_10m_prosail_values <- values(atbd_10m_prosail)
atbd_20m_prosail_values <- values(atbd_20m_prosail)
atbd_10m_snap_values <- values(atbd_10m_snap)
atbd_20m_snap_values <- values(atbd_20m_snap)
atbd_10m20m_snap_values <- values(atbd_10m20m_snap)
lidar_lai_values <- values(lidar_lai)

all_values <- cbind(
  atbd_10m_prosail_values,
  atbd_20m_prosail_values,
  atbd_10m_snap_values,
  atbd_20m_snap_values,
  atbd_10m20m_snap_values,
  lidar_lai_values
)
complete_rows <- complete.cases(all_values)
atbd_10m_prosail_values <- atbd_10m_prosail_values[complete_rows, , drop = FALSE]
atbd_20m_prosail_values <- atbd_20m_prosail_values[complete_rows, , drop = FALSE]
atbd_10m_snap_values <- atbd_10m_snap_values[complete_rows, , drop = FALSE]
atbd_20m_snap_values <- atbd_20m_snap_values[complete_rows, , drop = FALSE]
atbd_10m20m_snap_values <- atbd_10m20m_snap_values[complete_rows, , drop = FALSE]
lidar_lai_values <- lidar_lai_values[complete_rows, , drop = FALSE]

# stop()

# Metrics

# 10m
calculate_metrics(atbd_10m_prosail, 
                  atbd_10m_snap, 
                  plots_dir,
                  "10m_prosail",
                  "10m_snap")

# 20m
calculate_metrics(atbd_20m_prosail, 
                  atbd_20m_snap, 
                  plots_dir,
                  "20m_prosail",
                  "20m_snap")

# PROSAIL
calculate_metrics(atbd_10m_prosail, 
                  atbd_20m_prosail, 
                  plots_dir,
                  "10m_prosail",
                  "20m_prosail")

# Toolbox
calculate_metrics(atbd_10m_snap, 
                  atbd_20m_snap, 
                  plots_dir,
                  "10m_snap",
                  "20m_snap")

# All
plot_histogram(vars_list = list(atbd_10m_prosail_values, atbd_20m_prosail_values,
                                atbd_10m_snap_values, atbd_20m_snap_values), 
               title = "LAI Distributions", 
               xlab = "LAI", 
               var_labs = c("10m_prosail",
                            "20m_prosail",
                            "10m_snap",
                            "20m_snap"), 
               dirname = plots_dir, 
               filename = ("hist_10m_prosail_20m_prosail_10m_snap_20m_snap"))

plot_histogram(vars_list = list(atbd_10m_prosail_values, atbd_20m_prosail_values,
                                atbd_10m_snap_values, atbd_20m_snap_values,
                                atbd_10m20m_snap_values, lidar_lai_values), 
               title = "LAI Distributions", 
               xlab = "LAI", 
               var_labs = c("10m_prosail",
                            "20m_prosail",
                            "10m_snap",
                            "20m_snap",
                            "10mBP20m",
                            "LiDAR")
)



# cor_10m <- cor(atbd_10m_prosail_values, atbd_10m_snap_values, 
#                use = "pairwise.complete.obs")
# 
# cor_20m <- cor(atbd_20m_prosail_values, atbd_20m_snap_values, 
#                use = "pairwise.complete.obs")
# 
# plot_histogram(vars_list = list(atbd_10m_prosail_values, atbd_10m_snap_values),
#                title = sprintf("LAI Distributions 10m SNAP resampled %s", 
#                                resample_method),
#                xlab = "LAI",
#                var_labs = c("10m PROSAIL", "10m SNAP"),
#                dirname = snap_dir,
#                filename = sprintf("hist_lai_10m_resampling_%s", 
#                                   resample_method),
#                limits = c(0.001, 8)
# )
# 
# plot_density_scatterplot(var_x = atbd_10m_prosail_values,
#                          var_y = atbd_10m_snap_values,
#                          xlab = "10m PROSAIL",
#                          ylab = sprintf("10m SNAP resampled %s", 
#                                         resample_method),
#                          dirname = snap_dir,
#                          filename = sprintf("scatter_lai_10m_resampling_%s", 
#                                             resample_method),
#                          )
# 
# plot_histogram(vars_list = list(atbd_20m_prosail_values, atbd_20m_snap_values),
#                title = sprintf("LAI Distributions 20m SNAP resampled %s", 
#                                resample_method),
#                xlab = "LAI",
#                var_labs = c("20m PROSAIL", "20m SNAP"),
#                dirname = snap_dir,
#                filename = sprintf("hist_lai_20m_resampling_%s", 
#                                   resample_method),
#                limits = c(0.001, 8)
# )
# 
# plot_density_scatterplot(var_x = atbd_20m_prosail_values,
#                          var_y = atbd_20m_snap_values,
#                          xlab = "20m PROSAIL",
#                          ylab = sprintf("20m SNAP resampled %s", 
#                                         resample_method),
#                          dirname = snap_dir,
#                          filename = sprintf("scatter_lai_20m_resampling_%s", 
#                                             resample_method),
# )

# cor_prosail <- cor(atbd_10m_prosail_values, 
#                    atbd_20m_prosail_values, 
#                    use = "pairwise.complete.obs")
# 
# cor_toolbox <- cor(atbd_10m_snap_values, 
#                    atbd_20m_snap_values, 
#                    use = "pairwise.complete.obs")
# 
# plot_histogram(vars_list = list(atbd_10m_prosail_values, 
#                                 atbd_20m_prosail_values),
#                title = "LAI Distributions PROSAIL",
#                xlab = "LAI",
#                var_labs = c("10m PROSAIL", "20m PROSAIL"),
#                dirname = snap_dir,
#                filename = sprintf("hist_lai_prosail_resampling_%s", 
#                                   resample_method),
#                limits = c(0.001, 8)
# )
# 
# plot_density_scatterplot(var_x = atbd_10m_prosail_values,
#                          var_y = atbd_20m_prosail_values,
#                          xlab = "10m PROSAIL",
#                          ylab = "20m PROSAIL",
#                          dirname = snap_dir,
#                          filename = sprintf("scatter_lai_prosail_resampling_%s", 
#                                             resample_method),
# )
# 
# plot_histogram(vars_list = list(atbd_10m_snap_values,
#                                 atbd_20m_snap_values),
#                title = sprintf("LAI Distributions Toolbox SNAP resampled %s", 
#                                resample_method),
#                xlab = "LAI",
#                var_labs = c("10m SNAP", "20m SNAP"),
#                dirname = snap_dir,
#                filename = sprintf("hist_lai_snap_resampling_%s", 
#                                   resample_method),
#                limits = c(0.001, 8)
# )
# 
# plot_density_scatterplot(var_x = atbd_10m_snap_values,
#                          var_y = atbd_20m_snap_values,
#                          xlab = sprintf("10m SNAP resampled %s", 
#                                         resample_method),
#                          ylab = sprintf("20m SNAP resampled %s", 
#                                         resample_method),
#                          dirname = snap_dir,
#                          filename = sprintf("scatter_lai_snap_resampling_%s", 
#                                             resample_method),
# )