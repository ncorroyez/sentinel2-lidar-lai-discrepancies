# ---
# title: "plot_initial_field_scatterplots.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-11-13"
# ---

# ----------- (Optional) Clear the environment and free memory -----------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --- Define working directory as the directory where the script is located ----
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# ------------------------------ Libraries -------------------------------------

library("lidR")
library("data.table")
library("raster")
library("plotly")
library("terra")
library("viridis")
library("future")
library("ggplot2")
library("Metrics")

# --------------------------- Import useful functions --------------------------

source("libraries/functions_plots.R")
source("libraries/functions_general_tools.R")

# --------------------------------- Setup --------------------------------------
results_dir <- "../03_RESULTS"
figures_dir <- "../04_FIGURES"
forest_composition <- "Not_Masked"
metrics_dir <- "Metrics/Not_Masked"
s2lai_filename <- "s2lai_res_10_m.tif"
lidarlai_filename <- "lidarlai_res_10_m.tif"
# sites <- c("Aigoual", "Blois", "Mormal")
sites <- c("Aigoual")
# -------------------------------- Process -------------------------------------
all_predictors_filtered <- list()
results <- data.frame(Site = character(), R2 = numeric(), RMSE = numeric(), stringsAsFactors = FALSE)
for (i in 1:length(sites)) {
  cat("Processing site:", sites[i], "\n")
  # site_data <- load_metrics_for_site(sites[i], results_dir, forest_composition)
  site_path <- file.path(results_dir, sites[i])
  # data_utm31n <- geojson_read(file.path(site_path, "data_utm31n.geojson"), what = "sp")
  # coord_x <- data_utm31n@data$coord_x_utm31n
  # coord_y <- data_utm31n@data$coord_y_utm31n
  # coord_df <- data.frame(coord_x = coord_x, coord_y = coord_y)
  #
  # closest_matches <- list()
  # for (j in 1:nrow(coord_df)) {
  #   target_coord <- coord_df[j, ]
  #   distances <- sqrt((site_data$lat - target_coord$coord_y)^2 +
  #                       (site_data$lon - target_coord$coord_x)^2)
  #   closest_index <- which.min(distances)
  #   closest_matches[[j]] <- cbind(site_data[closest_index, ], target_coord)
  # }
  # closest_data <- do.call(rbind, closest_matches)
  # closest_data <- closest_data[, !colnames(closest_data) %in% c("coord_x", "coord_y")]
  # all_predictors_filtered[[sites[i]]] <- closest_data

  data_utm31n <- vect(file.path(site_path, "data_utm31n.geojson"))
  s2lai_path <- file.path(results_dir, sites[i], metrics_dir, s2lai_filename)
  s2lai <- terra::rast(s2lai_path)
  lidarlai_path <- file.path(results_dir, sites[i], metrics_dir, lidarlai_filename)
  lidarlai <- terra::rast(lidarlai_path)
  s2lai_values <- extract(s2lai, data_utm31n)
  lidarlai_values <- extract(lidarlai, data_utm31n)
  aaaaaaaaaaaaaa
  lidarlai_values <- all_predictors_filtered[[sites[i]]]$lidar_lai
  s2lai_values <- all_predictors_filtered[[sites[i]]]$s2_lai
  # Calculate true R2 and RMSE
  model <- lm(lidarlai_values ~ s2lai_values)
  predictions <- predict(model)
  ss_res <- sum((lidarlai_values - predictions)^2)
  ss_tot <- sum((lidarlai_values - mean(lidarlai_values))^2)
  r2 <- 1 - (ss_res / ss_tot)
  cor <- cor(lidarlai_values, s2lai_values)
  # r2 <- summary(model)$r.squared
  rmse <- rmse(s2lai_values, lidarlai_values)

  # Store the results in a data frame
  results <- rbind(results, data.frame(Site = sites[i], R2 = r2, RMSE = rmse, R = cor))

  # Initial scatterplot
  plot_data <- data.frame(S2_LAI = s2lai_values, LIDAR_LAI = lidarlai_values)
  print(dim(plot_data))
  plot <- ggplot(plot_data, aes(x = S2_LAI, y = LIDAR_LAI)) +
    geom_point(alpha = 0.5) +
    labs(
      x = "Sentinel-2 LAI",
      y = "LiDAR LAI",
      title = paste(sites[i], "Sentinel-2 LAI vs LiDAR LAI Values")
    ) +
    theme_bw() +
    xlim(0, max(s2lai_values)) +
    ylim(0, max(lidarlai_values))
  ggsave(filename = file.path(figures_dir, sites[i], paste0(sites[i], "_initial_field_lai_scatterplot.png")), plot = plot)

  # Quantiles
  variables <- c("lcv", "mean", "lskew", "vci", "slope", "hillshade",
                 "rumple", "cv_lad", "cv_lad_dtm", "cv", "dtm")
  for (var in variables) {
    qt_var <- all_predictors_filtered[[sites[i]]] %>%
      mutate(class = cut(.[[var]], breaks = quantile(.[[var]], probs = c(0, 1/3, 2/3, 1)),
                         labels = c("Low", "Medium", "High"), include.lowest = TRUE))
    qt_plot <- ggplot(data = qt_var, aes(x = s2_lai, y = lidar_lai, color = class)) +
      geom_point(size = 2, alpha = 0.7) +  # Scatter plot points with color by class
      geom_abline(intercept = 0, slope = 1, color = "darkred", linetype = "dashed", linewidth = 1) +
      labs(
        title = paste("Scatterplot LiDAR LAI vs S2 LAI with", var, "Classes"),
        x = "S2 LAI",
        y = "LiDAR LAI",
        color = paste("Class (Quantiles of", var, ")")
      ) +
      scale_color_manual(values = c("Low" = "steelblue", "Medium" = "goldenrod", "High" = "tomato")) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "right"
      )
    ggsave(filename = file.path(figures_dir, sites[i], paste0(sites[i], "_", var, "_lai_scatterplot.png")), plot = qt_plot)
  }
}
print(results)
