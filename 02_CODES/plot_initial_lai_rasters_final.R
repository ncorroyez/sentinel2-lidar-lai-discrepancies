# ---
# title: "plot_initial_lai_rasters.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2025-09-07"
# ---

rm(list=ls(all=TRUE))
gc()

if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# ------------------------------ Libraries -------------------------------------
library("tidyr")
library("terra")
library("viridis")
library("dplyr")
library("ggplot2")
library("patchwork")
library("grid")

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
# sites <- c("Aigoual", "Blois", "Mormal")
sites <- c("Mormal", "Blois", "Aigoual")

# -------------------------------- Collect all data ----------------------------
all_sites_data <- list()

for (i in 1:length(sites)) {
  cat("Loading site:", sites[i], "\n")
  
  s2lai_path <- file.path(results_dir, sites[i], metrics_dir, s2lai_filename)
  s2lai <- terra::rast(s2lai_path)
  
  lidarlai_path <- file.path(results_dir, sites[i], metrics_dir, lidarlai_filename)
  lidarlai <- terra::rast(lidarlai_path)
  
  na_index <- terra::countNA(s2lai) + terra::countNA(lidarlai)
  s2lai[na_index > 0] <- NA
  lidarlai[na_index > 0] <- NA
  
  stacked_rasters <- c(lidarlai, s2lai)
  df <- as.data.frame(stacked_rasters, xy = TRUE)
  colnames(df) <- c("x", "y", "LiDAR LAI", "S2 LAI")
  
  df_cleaned <- df %>%
    remove_outliers("LiDAR LAI") %>%
    remove_outliers("S2 LAI")
  
  df_melted <- df_cleaned %>%
    pivot_longer(cols = c("LiDAR LAI", "S2 LAI"),
                 names_to = "Type",
                 values_to = "LAI")
  
  df_melted$Type <- factor(df_melted$Type,
                           levels = c("LiDAR LAI", "S2 LAI"),
                           labels = c("LAI[ALS]", "LAI[S2]"))
  df_melted$Site <- sites[i]
  
  all_sites_data[[i]] <- df_melted
}

# ------------------------ Compute global color scale --------------------------
all_sites_df <- do.call(rbind, all_sites_data)
lai_min <- min(all_sites_df$LAI, na.rm = TRUE)
lai_max <- max(all_sites_df$LAI, na.rm = TRUE)
cat("Global LAI range:", lai_min, "to", lai_max, "\n")

# ------------------------ Plot each site separately ---------------------------
site_plots <- list()

for (i in 1:length(sites)) {
  df_melted <- all_sites_data[[i]]
  
  p_raster <- ggplot(df_melted, aes(x = x, y = y, fill = LAI)) +
    geom_raster() +
    scale_fill_viridis_c(option = "C", name = "LAI",
                         limits = c(lai_min, lai_max)) +
    facet_wrap(~ Type, nrow = 1, labeller = label_parsed) +
    coord_equal() +
    labs(title = sites[i],
         x = NULL, # Explicitly remove labels here
         y = NULL) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 20, face = "bold"),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      
      # --- REQUESTED CHANGES ---
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  site_plots[[i]] <- p_raster
  
  # Save individual
  output_dir <- file.path(figures_dir, sites[i])
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  ggsave(filename = file.path(output_dir, paste0(sites[i], "_lai_raster_comparison.png")),
         plot = p_raster, width = 12, height = 6)
  
  cat("Raster plot saved for site:", sites[i], "\n")
}

# ------------------------ Assemble into one figure ----------------------------
combined_plot <- (site_plots[[1]] / site_plots[[2]] / site_plots[[3]]) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(filename = file.path(figures_dir, "all_sites_lai_raster_comparison.png"),
       plot = combined_plot, width = 12, height = 12, dpi = 300)

cat("Combined raster plot saved for all sites\n")
