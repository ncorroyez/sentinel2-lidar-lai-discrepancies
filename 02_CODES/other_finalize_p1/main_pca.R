# ---
# title: "main_pca.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-09-27"
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
library(terra)
library(ggplot2)
library(dplyr)
library(factoextra)
library(vegan)
library(reshape2)
source("libraries/functions_general_tools.R")

# Pre-processing Parameters
# sites <- c("Aigoual", "Blois", "Mormal")
sites <- c("Aigoual", "Blois")
# sites <- "Aigoual" # Mormal Blois Aigoual
results_dir <- "../03_RESULTS"
forest_composition <- "Deciduous_Only" # Full_Composition Deciduous_Only Not_Masked
metrics_dir <- file.path("Metrics", forest_composition)

metric_files <- c("lidarlai_res_10_m.tif",
                  "s2lai_summer_res_10_m.tif",
                  "lcv_res_10_m.tif",
                  "mean_res_10_m.tif",
                  "max_res_10_m.tif",
                  "lskew_res_10_m.tif",
                  "hillshade_res_10_m.tif",
                  "slope_res_10_m.tif",
                  "aspect_res_10_m.tif",
                  "std_res_10_m.tif",
                  "cv_res_10_m.tif",
                  "rumple_res_10_m.tif",
                  "vci_res_10_m.tif",
                  "fCover_res_10_m.tif",
                  "ndvi_summer_res_10_m.tif",
                  "cv_lad_res_10_m.tif",
                  "cv_lad_dtm_res_10_m.tif"
)

for (site in sites){
  # metrics_files <- list.files(path = file.path(results_dir, site, metrics_dir), 
  #                             pattern = "\\.tif$", full.names = TRUE)
  metrics_files <- file.path(results_dir, site, metrics_dir, metric_files)
  
  # Read all the raster files and create a SpatRaster
  raster_stack <- rast(metrics_files)
  raster_values <- values(raster_stack)
  raster_values <- na.omit(raster_values)
  raster_df <- as.data.frame(raster_values)
  indices <- sample(1:nrow(raster_df), size = 120)
  raster_df <- raster_df[indices, ]
  colnames(raster_df) <- gsub("\\_res_10_m.tif$", "", basename(metrics_files))
  
  # Perform PCA using prcomp()
  pca_result <- prcomp(raster_df, center = TRUE, scale. = TRUE)
  
  # 1. Scree plot - Explained variance of each component
  fviz_eig(pca_result, addlabels = TRUE, barfill = "steelblue", 
           barcolor = "steelblue", linecolor = "red") + 
    labs(title = "Scree Plot", x = "Principal Components", y = "Explained Variance (%)")
  
  # 2. Biplot - Shows both variable loadings and PCA scores
  # The biplot combines the PCA scores (dots) and loadings (arrows)
  fviz_pca_biplot(pca_result, 
                  label = "var",           # Show variable labels (raster layers)
                  geom.ind = "point",      # Show PCA scores as points
                  repel = TRUE) +          # Avoid overlapping labels
    labs(title = paste(site, "Biplot of Principal Components"))
  
  # 3. Individual Component Loadings Plot (Contribution of variables to PCs)
  # This plot shows how much each original raster contributes to each component.
  # The higher the value, the more the variable contributes to the component.
  
  # Plot the loadings of raster layers on the first two principal components
  loadings <- pca_result$rotation
  loadings_df <- as.data.frame(loadings)
  
  # For PC1
  fviz_contrib(pca_result, choice = "var", axes = 1, top = ncol(loadings_df)) +
    labs(title = "Contributions of Variables to PC1")
  
  # For PC2
  fviz_contrib(pca_result, choice = "var", axes = 2, top = ncol(loadings_df)) +
    labs(title = "Contributions of Variables to PC2")
  
  fviz_contrib(pca_result, choice = "var", axes = 1:4, top = ncol(loadings_df)) +
    labs(title = "Contributions of Variables to PC1-4")
  
  # NMS using vegan's metaMDS function
  nms_result <- metaMDS(raster_df, distance = "bray", trymax = 100)
  
  # 1. Scree plot for NMS (Plot the stress of NMS)
  stressplot(nms_result)
  
  # 2. NMS Ordination Plot
  plot(nms_result, type = "n")  # Create a blank plot
  points(nms_result, display = "sites", col = "blue", pch = 16)  # Plot the sites (samples)
  # text(nms_result, display = "sites", col = "blue", pos = 3)  # Label the sites (optional)
  
  # 3. Biplot (Environmental Variables or Variables in the data)
  # This biplot will show the loadings of variables in NMS
  # biplot(nms_result, display = "species")
  
  # Optionally, plot the stress and explore the configuration
  # NMS results are stored in the "points" slot for site scores
  nms_scores <- as.data.frame(scores(nms_result))
  colnames(nms_scores) <- c("NMDS1", "NMDS2")
  
  # Plotting NMS1 vs NMS2 with ggplot2
  ggplot(nms_scores, aes(x = NMDS1, y = NMDS2)) +
    geom_point(color = "blue") +
    labs(title = paste(site, "NMS Ordination"), x = "NMDS1", y = "NMDS2")
  
  # Calculate the correlation between the variables and the NMDS axes
  correlations_nms <- cor(raster_df, nms_scores)
  
  # For each axis (NMDS1 and NMDS2), visualize the contribution of variables
  correlation_df <- as.data.frame(correlations_nms)
  correlation_df$variable <- rownames(correlation_df)
  
  # Plot the correlation of variables with NMDS1 and NMDS2
  ggplot(correlation_df, aes(x = variable)) +
    geom_bar(aes(y = NMDS1), stat = "identity", fill = "blue", alpha = 0.6) +
    geom_bar(aes(y = NMDS2), stat = "identity", fill = "red", alpha = 0.6) +
    coord_flip() +
    labs(title = paste(site, "Variable Contributions to NMDS1 and NMDS2"), 
         x = "Variable", y = "Correlation with NMDS Axes") +
    theme_minimal() +
    scale_y_continuous(limits = c(-1, 1))  # Limiting the y-axis for better visualization
}





# Create an empty list to store data from all sites
combined_raster_df <- list()

# Iterate over the sites to combine their data
for (site in sites){
  
  # Read in the metrics files for the current site
  metrics_files <- file.path(results_dir, site, metrics_dir, metric_files)
  
  # Read all the raster files and create a SpatRaster
  raster_stack <- rast(metrics_files)
  raster_values <- values(raster_stack)
  raster_values <- na.omit(raster_values)
  raster_df <- as.data.frame(raster_values)
  
  indices <- sample(1:nrow(raster_df), size = 120)
  raster_df <- raster_df[indices, ]
  
  colnames(raster_df) <- gsub("\\_res_10_m.tif$", "", basename(metrics_files))
  
  # Add site information for future identification
  raster_df$site <- site
  
  # Append the data to the list
  combined_raster_df[[site]] <- raster_df
}

# Combine all data from all sites into one data frame
combined_raster_df <- do.call(rbind, combined_raster_df)

# Remove site column before performing NMS
raster_data_combined <- combined_raster_df[, !names(combined_raster_df) %in% c("site")]

# Perform NMS using vegan's metaMDS function
nms_result_combined <- metaMDS(raster_data_combined, distance = "bray", trymax = 100)

# 1. Stress plot for NMS (Plot the stress of NMS)
stressplot(nms_result_combined)

# 2. NMS Ordination Plot
plot(nms_result_combined, type = "n")  # Create a blank plot
points(nms_result_combined, display = "sites", col = "blue", pch = 16)  # Plot the sites (samples)

# 3. Correlation of Variables with NMS Axes (NMDS1 and NMDS2)
nms_scores_combined <- as.data.frame(scores(nms_result_combined))
colnames(nms_scores_combined) <- c("NMDS1", "NMDS2")

# Calculate the correlation between the variables and the NMDS axes
correlations_nms_combined <- cor(raster_data_combined, nms_scores_combined)

# For each axis (NMDS1 and NMDS2), visualize the contribution of variables
correlation_nms_df <- as.data.frame(correlations_nms_combined)
correlation_nms_df$variable <- rownames(correlation_nms_df)

# Plot the correlation of variables with NMDS1 and NMDS2
ggplot(correlation_nms_df, aes(x = variable)) +
  geom_bar(aes(y = NMDS1), stat = "identity", fill = "blue", alpha = 0.6) +
  geom_bar(aes(y = NMDS2), stat = "identity", fill = "red", alpha = 0.6) +
  coord_flip() +
  labs(title = paste("All Sites Variable Contributions to NMDS1 and NMDS2"), 
       x = "Variable", y = "Correlation with NMDS Axes") +
  theme_minimal() +
  scale_y_continuous(limits = c(-1, 1))  # Limiting the y-axis for better visualization
