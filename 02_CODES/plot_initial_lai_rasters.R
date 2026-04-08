# ---
# title: "plot_initial_lai_rasters.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-10-01"
# ---

# ----------- (Optional) Clear the environment and free memory -----------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --- Define working directory as the directory where the script is located ----
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# ------------------------------ Libraries -------------------------------------

library("tidyr")
library("plotly")
library("terra")
library("viridis")
library("dplyr")
library("ggplot2")
library("ggnewscale")

# --------------------------- Import useful functions --------------------------

source("libraries/functions_plots.R")

remove_outliers <- function(df, col) {
  Q1 <- quantile(df[[col]], 0.25, na.rm = TRUE)
  Q3 <- quantile(df[[col]], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  df <- df %>% filter(df[[col]] >= lower_bound & df[[col]] <= upper_bound)
  return(df)
}

# --------------------------------- Setup --------------------------------------
results_dir <- "../03_RESULTS"
figures_dir <- "../04_FIGURES"
metrics_dir <- "Metrics/Deciduous_Only"
s2lai_filename <- "s2lai_summer_v2_res_10_m.tif"
lidarlai_filename <- "lidarlai_res_10_m.tif"
sites <- c("Aigoual", "Blois", "Mormal")

# -------------------------------- Process -------------------------------------
all_sites_data <- list()
for (i in 1:length(sites)) {
  cat("Processing site:", sites[i], "\n")
  
  s2lai_path <- file.path(results_dir, sites[i], metrics_dir, s2lai_filename)
  s2lai <- terra::rast(s2lai_path)
  
  lidarlai_path <- file.path(results_dir, sites[i], metrics_dir, lidarlai_filename)
  lidarlai <- terra::rast(lidarlai_path)
  
  na_index_s2lai <- terra::countNA(s2lai)
  na_index_lidarlai <- terra::countNA(lidarlai)
  na_index <- na_index_s2lai + na_index_lidarlai
  s2lai[na_index > 0] <- NA
  lidarlai[na_index > 0] <- NA
  
  # Calculate Diff LAI
  delta_lai <- lidarlai - s2lai
  
  # Stack the rasters
  stacked_rasters <- c(lidarlai, s2lai, delta_lai)
  
  # Convert to data frame
  df <- as.data.frame(stacked_rasters, xy = TRUE)
  colnames(df) <- c("x", "y", "LiDAR_LAI", "Sentinel_2_LAI", "Delta_LAI")
  
  # Remove outliers for each LAI type
  df_cleaned <- df %>%
    remove_outliers("LiDAR_LAI") %>%
    remove_outliers("Sentinel_2_LAI") %>%
    remove_outliers("Delta_LAI")
  
  # Use pivot_longer instead of melt
  df_melted <- df_cleaned %>% # df_cleaned df
    pivot_longer(cols = c(LiDAR_LAI, Sentinel_2_LAI, Delta_LAI),
                 names_to = "Type",
                 values_to = "LAI")
  df_melted$Type <- factor(df_melted$Type, levels = c("LiDAR_LAI", "Sentinel_2_LAI", "Delta_LAI"))
  
  # Add site column for faceting
  df_melted$Site <- sites[i]
  all_sites_data[[i]] <- df_melted
  
  # Create output directory if it doesn't exist
  output_dir <- file.path(figures_dir, sites[i])
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create the plot
  p <- ggplot(df_melted, aes(x = x, y = y)) +
    # Layer for LiDAR_LAI and Sentinel_2_LAI with the first fill scale
    geom_raster(data = subset(df_melted, Type %in% c("LiDAR_LAI", "Sentinel_2_LAI")), 
                aes(fill = LAI)) +
    scale_fill_viridis_c(option = "C", name = "LiDAR and Sentinel-2 LAI") +  # Color scale for LAI values
    
    # Add a new scale for Diff_LAI
    new_scale_fill() +
    
    # Layer for Diff_LAI with a different fill scale
    geom_raster(data = subset(df_melted, Type == "Delta_LAI"), 
                aes(fill = LAI)) +
    scale_fill_viridis_c(option = "B", name = "Delta LAI") +  # Color scale for difference values
    
    facet_wrap(~ Type, nrow = 1) +
    labs(title = paste(sites[i], "LAI Comparisons"),
         x = "Longitude",
         y = "Latitude") +
    theme_bw() +
    theme(legend.position = "bottom")
  # plot(p)
  
  # Save the plot
  ggsave(filename = file.path(output_dir, paste0(sites[i], "_lai_comparison.png")),
         plot = p, width = 1920, height = 1080, units = "px", dpi = 150)
  
  # Create the boxplot for each LAI type
  p_boxplot <- ggplot(df_melted, aes(x = Type, y = LAI, fill = Type)) +
    geom_boxplot(outlier.shape = NA) +  # Hide outliers
    labs(title = paste(sites[i], "LAI Boxplot Comparison"),
         x = "LAI Type", 
         y = "LAI Values") +
    scale_fill_viridis_d() +  # Discrete viridis scale for boxplots
    theme_bw() +
    theme(legend.position = "none") + # No legend for boxplots
    coord_cartesian(ylim = range(df_melted$LAI, na.rm = TRUE))  # Set y-axis to the data range
  
  # Save the boxplot
  ggsave(filename = file.path(output_dir, paste0(sites[i], "_lai_boxplot_comparison.png")),
         plot = p_boxplot, width = 8, height = 6)
  
  cat("Finished processing site:", sites[i], "\n")
}

# Combine data for all sites into one data frame
all_sites_df <- do.call(rbind, all_sites_data)
all_sites_df$Type <- factor(all_sites_df$Type, levels = c("LiDAR_LAI", "Sentinel_2_LAI", "Delta_LAI"))

# Plot the combined boxplot with facets for each site
p_boxplot <- ggplot(all_sites_df, aes(x = Type, y = LAI, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "LAI Boxplot Comparison Across Sites",
       x = "LAI Type", 
       y = "LAI Values") +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~ Site) +  # Facet by site
  coord_cartesian(ylim = range(all_sites_df$LAI, na.rm = TRUE))

# Save the combined boxplot
ggsave(filename = file.path(figures_dir, "combined_lai_boxplot_comparison.png"),
       plot = p_boxplot, width = 12, height = 6)

cat("Boxplot for all sites saved successfully.\n")