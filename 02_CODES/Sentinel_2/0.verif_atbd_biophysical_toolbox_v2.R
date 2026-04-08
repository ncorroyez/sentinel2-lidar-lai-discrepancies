# 0.verif_atbd_biophysical_toolbox_10m_only.R
# Author: Nathan CORROYEZ, UMR TETIS
# Last update: 2025-07-30

# ---- Clear environment & set working directory ----
rm(list=ls(all=TRUE)); gc()
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

# ---- Libraries and functions ----
library(terra)
library(Metrics)
library(dplyr)
library(tidyr)
library(knitr)
library(ggplot2)
source("../libraries/functions_sentinel_2.R")         # if needed for reading
source("../libraries/functions_plots.R")              # for plot_histogram, etc.
source("../libraries/functions_intercomparisons.R")   # for calculate_metrics()

# ---- Paths & filenames ----
# site        <- "Mormal"
sites       <- c("Aigoual", "Blois", "Mormal")

# Init table to store all results
results <- data.frame(
  Site      = character(),
  Comparison = character(),
  r         = numeric(),
  RMSE      = numeric(),
  stringsAsFactors = FALSE
)

for (site in sites){
  mask_10m    <- rast(paste0("/home/corroyez/Documents/NC_Full/03_RESULTS/", site, "/LiDAR/Heterogeneity_Masks/artifacts_low_vegetation_majority_90_p_res_10_m.envi"))
  
  # Sentinel-2 ATBD products (10 m)
  s2_dir      <- file.path("../../03_RESULTS", site, "Metrics/Deciduous_Only")
  file_atbd   <- file.path(s2_dir, "s2lai_summer_atbd_res_10_m.tif")
  file_codist <- file.path(s2_dir, "s2lai_summer_atbd_codist_false_res_10_m.tif")
  
  # SNAP Toolbox LAI (10 m → then masked)
  snap_dir        <- "../../../Softwares/SNAP_Toolbox_Results" 
  snap_image      <- switch(site,
                            "Aigoual" = "S2A_MSIL2A_20210711T104031_N0301_R008_T31TEJ_20210711T134956",
                            "Blois"   = "S2A_MSIL2A_20210614T105031_N0300_R051_T31TCN_20210614T140120",
                            "Mormal"  = "S2A_MSIL2A_20210614T105031_N0300_R051_T31UER_20210614T140120")
  file_snap_lai   <- file.path(snap_dir, site,
                               paste0(snap_image, "_resampled_biophysical10m.data"),
                               "lai.img")
  
  plots_dir <- file.path(s2_dir, "Plots/verif_toolbox_10m")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive=TRUE)
  
  # ---- Load & mask rasters ----
  # ATBD standard
  atbd     <- rast(file_atbd)    |> project(mask_10m) |> mask(mask_10m)
  # ATBD codist=false
  codist   <- rast(file_codist)  |> project(mask_10m) |> mask(mask_10m)
  # SNAP toolbox
  snap_lai <- rast(file_snap_lai)|> project(mask_10m) |> mask(mask_10m)
  
  # ---- Extract values & keep only complete cases ----
  v_atbd   <- values(atbd)
  v_codist <- values(codist)
  v_snap   <- values(snap_lai)
  
  all_vals <- cbind(v_atbd, v_codist, v_snap)
  good     <- complete.cases(all_vals)
  v_atbd   <- v_atbd[good, , drop=FALSE]
  v_codist <- v_codist[good, , drop=FALSE]
  v_snap   <- v_snap[good, , drop=FALSE]
  
  # ---- Comparison metrics ----
  # Compare standard ATBD vs codist‐false
  # calculate_metrics(atbd, codist, plots_dir,
  #                   "ATBD_10m", "ATBD_codistfalse_10m")
  # 
  # # Compare standard ATBD vs SNAP
  # calculate_metrics(atbd, snap_lai, plots_dir,
  #                   "ATBD_10m", "SNAP_10m")
  # 
  # # Compare codist‐false vs SNAP
  # calculate_metrics(codist, snap_lai, plots_dir,
  #                   "ATBD_codistfalse_10m", "SNAP_10m")
  # 
  # ---- Histograms of LAI distributions ----
  # plot_histogram(
  #   vars_list = list(v_atbd, v_codist, v_snap),
  #   title     = "LAI Distributions (10 m)",
  #   xlab      = "LAI",
  #   var_labs  = c("ATBD", "ATBD_codist_false", "SNAP"),
  #   dirname   = plots_dir,
  #   filename  = "hist_lai_10m_comparison"
  # )
  # 
  # # ---- (Optional) Scatterplots ----
  # plot_density_scatterplot(
  #   var_x    = v_atbd, var_y = v_snap,
  #   xlab     = "ATBD LAI (10 m)",
  #   ylab     = "SNAP LAI (10 m)",
  #   dirname  = plots_dir,
  #   filename = "scatter_ATBD_vs_SNAP_10m"
  # )
  # 
  # plot_density_scatterplot(
  #   var_x    = v_codist, var_y = v_snap,
  #   xlab     = "ATBD_codistfalse LAI (10 m)",
  #   ylab     = "SNAP LAI (10 m)",
  #   dirname  = plots_dir,
  #   filename = "scatter_ATBDcodist_vs_SNAP_10m"
  # )
  
  
  # Append to results table
  results <- bind_rows(results,
                       data.frame(Site = site, 
                                  Comparison = "LAI[S2_ATBD_T]*\" vs \"*LAI[S2_ATBD_F]",
                                  r    = cor(v_atbd[,1], v_codist[,1]), 
                                  RMSE = Metrics::rmse(v_atbd[,1], v_codist[,1])),
                       
                       data.frame(Site = site, 
                                  Comparison = "LAI[S2_ATBD_T]*\" vs \"*LAI[SNAP]",
                                  r    = cor(v_atbd[,1], v_snap[,1]), 
                                  RMSE = Metrics::rmse(v_atbd[,1], v_snap[,1])),
                       
                       data.frame(Site = site, 
                                  Comparison = "LAI[S2_ATBD_F]*\" vs \"*LAI[SNAP]",
                                  r    = cor(v_codist[,1], v_snap[,1]), 
                                  RMSE = Metrics::rmse(v_codist[,1], v_snap[,1]))
  )
}
kable(results, digits = 3)

# Pivot longer for faceted heatmap
results_long <- results |>
  pivot_longer(cols = c(r, RMSE), names_to = "Metric", values_to = "Value")
results_long$Site <- factor(results_long$Site, levels = rev(c("Aigoual", "Blois", "Mormal")))

metric_labeller <- as_labeller(c(
  "r"    = "italic(r)",
  "RMSE" = "RMSE"
), default = label_parsed)

p <- ggplot(results_long, aes(x = Comparison, y = Site, fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Value, 2)), color = "black", size = 4) +
  # Use the custom labeller in facet_wrap
  facet_wrap(~Metric, nrow = 1, scales = "free", labeller = metric_labeller) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "Compared S2 ATBD Distributions",
       y = "Site",
       fill = "Value") +
  scale_x_discrete(labels = function(x) parse(text = x))

ggsave(
  filename = "comp_atbd.png",
  plot = p,
  width = 10,
  height = 10,
  dpi = 300
)
