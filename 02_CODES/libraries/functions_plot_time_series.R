# ---
# title: "functions_plot_time_series.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-09-26"
# ---

# Function to prepare LAI data from rasters and masks
prepare_lai_data <- function(sites, results_dir, dirs_to_keep, structural_metrics) {
  
  lai_data <- data.frame(site = character(),
                         composition = character(),
                         quantile = integer(),
                         layer = numeric(),
                         lai_mean = numeric(),
                         lai_std = numeric(),
                         structural_metric = character())
  
  for (site in sites){
    cat("Processing site: ", site, "\n")
    
    masks_dir <- file.path(results_dir, site, "LiDAR/Heterogeneity_Masks/Quantiles")
    metrics_dir <- file.path(results_dir, site, "Metrics")
    
    composition_dirs <- list.dirs(metrics_dir, recursive = FALSE)
    composition_dirs <- composition_dirs[grepl(paste(dirs_to_keep, collapse = "|"), composition_dirs)]
    
    for (composition_dir in composition_dirs) {
      
      composition <- basename(composition_dir)
      cat("Processing composition: ", composition, "\n")
      
      lai_rasters <- terra::rast(file.path(composition_dir, paste0(site, "_lai_stacked_raster.tif")))
      layer_names <- names(lai_rasters)
      
      for (structural_metric in structural_metrics) {
        cat("Processing metric: ", structural_metric, "\n")
        
        quantiles_dir <- file.path(masks_dir, composition, structural_metric, "Deciles", "3_Quantiles")
        quantiles_low <- c(0, 0.33, 0.67)
        quantiles_high <- c(0.33, 0.66, 1)
        
        for (quantile_idx in seq_along(quantiles_low)) {
          cat("Processing quantile: ", quantile_idx, "\n")
          
          qt <- terra::rast(file.path(quantiles_dir, paste0(structural_metric, "_heter_binary_mask_res_10_m_", 
                                                            quantiles_low[quantile_idx], "_", quantiles_high[quantile_idx], ".tif")))
          
          for (layer_idx in 1:nlyr(lai_rasters)) {
            lai_raster <- lai_rasters[[layer_idx]]
            lai_masked <- mask(lai_raster, qt)
            lai_values <- values(lai_masked, na.rm = TRUE)
            
            lai_mean <- mean(lai_values, na.rm = TRUE)
            lai_std <- sd(lai_values, na.rm = TRUE)
            
            lai_data <- rbind(lai_data, 
                              data.frame(site = site,
                                         composition = composition,
                                         quantile = quantile_idx,
                                         layer = layer_names[layer_idx],
                                         lai_mean = lai_mean,
                                         lai_std = lai_std,
                                         structural_metric = structural_metric))
          }
        }
      }
    }
  }
  
  # Convert 'layer' to YYYY-MM (year-month)
  lai_data$layer <- as.yearmon(substr(lai_data$layer, 1, 6), "%Y%m")
  
  return(lai_data)
}

# Function to plot LAI data faceted by sites
plot_lai_facet_sites <- function(lai_data, output_dir) {
  for (composition_name in unique(lai_data$composition)) {
    for (metric_name in unique(lai_data$structural_metric)) {
      
      composition_metric_data <- lai_data[lai_data$composition == composition_name & 
                                            lai_data$structural_metric == metric_name, ]
      
      if (nrow(composition_metric_data) > 0) {
        p <- ggplot(composition_metric_data, aes(x = layer, y = lai_mean, fill = factor(quantile))) + #color = factor(quantile))) +
          geom_boxplot() +
          # geom_line() +
          # geom_point() +
          # geom_errorbar(aes(ymin = lai_mean - lai_std, ymax = lai_mean + lai_std), width = 0.2) +
          labs(x = "Month", y = "Mean LAI", color = "Quantile Class",
               title = paste("Mean LAI with Standard Deviation for", composition_name, "using", metric_name)) +
          facet_wrap(~ site) + 
          theme_bw() +
          theme(legend.position = "top",
                axis.text.x = element_text(angle = 45, hjust = 1),
                plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
        
        print(p)
        
        ggsave(filename = paste0(output_dir, "/LAI_plot_3sites_", composition_name, "_", metric_name, ".png"), 
               plot = p, 
               width = 16,  # increase width
               height = 10,  # increase height
               units = "in",  # set units to inches
               dpi = 300)  # increase dpi for better resolution
      }
    }
  }
}

# Function to plot LAI data separately for each site, composition, and metric
plot_lai_separate <- function(lai_data, output_dir) {
  for (site_name in unique(lai_data$site)) {
    for (composition_name in unique(lai_data$composition)) {
      for (metric_name in unique(lai_data$structural_metric)) {
        
        site_composition_metric_data <- lai_data[lai_data$site == site_name &
                                                   lai_data$composition == composition_name & 
                                                   lai_data$structural_metric == metric_name, ]
        
        if (nrow(site_composition_metric_data) > 0) {
          p <- ggplot(site_composition_metric_data, aes(x = layer, y = lai_mean, fill = factor(quantile))) + # color = factor(quantile), group = quantile)) +
            geom_boxplot() +
            # geom_line() +
            # geom_point() +
            # geom_errorbar(aes(ymin = lai_mean - lai_std, ymax = lai_mean + lai_std), width = 0.2) +
            labs(x = "Month", y = "Mean LAI", color = "Quantile Class",
                 title = paste("Mean LAI with Standard Deviation for", site_name, 
                               "\nComposition:", composition_name, "| Metric:", metric_name)) +
            scale_x_yearmon(n.breaks = 20, format = "%b %Y") +
            theme_bw() +
            theme(legend.position = "top",
                  axis.text.x = element_text(angle = 45, hjust = 1),
                  plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
          
          print(p)
          
          ggsave(filename = paste0(output_dir, "/LAI_plot_1site_", site_name, "_", composition_name, "_", metric_name, ".png"), 
                 plot = p, 
                 width = 16,  # increase width
                 height = 10,  # increase height
                 units = "in",  # set units to inches
                 dpi = 300)  # increase dpi for better resolution
        }
      }
    }
  }
}

# Apply sliding window to the LAI raster using a focal function
# focal_fun <- function(values) {
#   if (all(is.na(values))) return(NA)
#   mean_val <- mean(values, na.rm = TRUE)
#   sd_val <- sd(values, na.rm = TRUE)
#   return(c(mean_val, sd_val))
# }

# Define a custom focal function for mean
focal_mean_fun <- function(values) {
  if (all(is.na(values))) return(NA)
  return(mean(values, na.rm = TRUE))
}

# Define a custom focal function for standard deviation
focal_sd_fun <- function(values) {
  if (all(is.na(values))) return(NA)
  return(sd(values, na.rm = TRUE))
}
