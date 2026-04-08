# ---
# title: "main_clustering.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-12-19"
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
library(cluster)
library(clusterSim)
library(geojsonio)
library(car)
source("libraries/functions_general_tools.R")

# Pre-processing Parameters
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- c("Blois", "Mormal")
# sites <- "Aigoual" # Mormal Blois Aigoual
results_dir <- "../03_RESULTS"
forest_composition <- "Not_Masked" # Full_Composition Deciduous_Only Not_Masked
metrics_dir <- file.path("Metrics", forest_composition)

# --------------------------------- Loading ------------------------------------
all_observations <- list()
all_predictors_filtered <- list()
for (site in sites) {
  cat("Site:", site, "\n")
  # Load metrics with padding and combine with spatial data
  site_data <- load_metrics_with_pad(site, results_dir, forest_composition)
  all_observations[[site]] <- site_data
  
  # Read GeoJSON file and extract coordinates
  geojson_file <- file.path(results_dir, site, "data_utm31n.geojson")
  geo_data <- geojson_read(geojson_file, what = "sp")
  coord_x <- geo_data@data$coord_x_utm31n
  coord_y <- geo_data@data$coord_y_utm31n
  coord_df <- data.frame(coord_x = coord_x, coord_y = coord_y)
  
  # Find closest matches
  closest_matches <- list()
  for (j in 1:nrow(coord_df)) {
    target_coord <- coord_df[j, ]
    distances <- sqrt((site_data$lat - target_coord$coord_y)^2 +
                        (site_data$lon - target_coord$coord_x)^2)
    closest_index <- which.min(distances)
    closest_matches[[j]] <- cbind(site_data[closest_index, ], target_coord)
  }
  closest_data <- do.call(rbind, closest_matches)
  closest_data <- closest_data[, !colnames(closest_data) %in% c("coord_x", "coord_y")]
  closest_data <- closest_data[
    closest_data$stands %in% c("oak", "beech", "deciduous", "poplar") &
      !is.na(closest_data$lidar_lai) &
      closest_data$mean >= 2,
  ]
  all_predictors_filtered[[site]] <- closest_data
}
combined_predictors_filtered <- do.call(rbind, all_predictors_filtered)

# -------------------------------- Variables -----------------------------------
selected_columns <- c(
                      "lcv",
                      "mean", 
                      # "max",
                      "vci",
                      "rumple",
                      # "std",
                      "lskew",
                      # "cv",
                      # "slope"
                      # "stands",
                      "fCover"
                      # "s2_lai",
                      # "hillshade"
                      # "deltaLAI",
                      # "lidar_lai"
)

# vif_data <- combined_predictors_filtered[, selected_columns]
# vif_data <- na.omit(vif_data)  # Remove rows with missing values
# lm_model <- lm(lskew ~ ., data = vif_data)
# vif_values <- vif(lm_model)
# print(vif_values)

# -------------------------------- All sites -----------------------------------
# Perform PCA on combined filtered predictors
pca_prep <- combined_predictors_filtered[, selected_columns]
correlation_matrix <- cor(combined_predictors_filtered[, selected_columns], use = "complete.obs")
print(correlation_matrix)
# pca_prep$stands <- as.numeric(pca_prep$stands)
pca_result <- prcomp(pca_prep, center = TRUE, scale. = TRUE)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Site <- rep(sites, sapply(all_predictors_filtered, nrow))

# Clustering (K-means)
k <- 2  # Number of clusters
set.seed(42)  # For reproducibility
pca_scores$Cluster <- kmeans(pca_scores[, 1:2], centers = k)$cluster
# pca_scores$Cluster <- pam(pca_scores[, 1:2], k)$clustering


# Final PCA Plot with Clusters and Sites
final_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = as.factor(Cluster), shape = Site)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1.5) +
  labs(
    title = "PCA Clustering Across Sites",
    x = "PC1",
    y = "PC2",
    color = "Cluster",
    shape = "Site"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
print(final_plot)

# Scree Plot - Explained Variance
scree_plot <- fviz_eig(pca_result, addlabels = TRUE, barfill = "steelblue", 
                       barcolor = "steelblue", linecolor = "red") +
  labs(title = "Scree Plot", x = "Principal Components", y = "Explained Variance (%)")
print(scree_plot)

# Biplot - Scores and Loadings
biplot <- fviz_pca_biplot(pca_result, 
                          label = "var",         # Show variable labels
                          repel = TRUE,          # Avoid overlapping labels
                          geom.ind = "point",    # Show PCA scores as points
                          pointshape = 21, 
                          pointsize = 2, 
                          fill.ind = as.factor(pca_scores$Cluster), # Color by cluster
                          palette = "jco",       # Cluster-based palette
                          shape.ind = as.factor(pca_scores$Site)   # Shape by site
) +
  labs(title = "Biplot of Principal Components", 
       fill = "Cluster", 
       shape = "Site")   # Add meaningful legends
print(biplot)

# Variable Contribution to PCs
contrib_pc1 <- fviz_contrib(pca_result, choice = "var", axes = 1, top = 10) +
  labs(title = "Variable Contributions to PC1")
print(contrib_pc1)
contrib_pc2 <- fviz_contrib(pca_result, choice = "var", axes = 2, top = 10) +
  labs(title = "Variable Contributions to PC2")
print(contrib_pc2)

# Analyze Cluster 2
cluster2_mormal_blois <- pca_scores %>%
  filter(Cluster == 2 & Site %in% c("Mormal", "Blois"))
original_filtered_data <- combined_predictors_filtered %>%
  filter(site %in% c("Mormal", "Blois")) %>%
  filter(row.names(.) %in% row.names(cluster2_mormal_blois))

# Metrics
silhouette_scores <- silhouette(pca_scores$Cluster, dist(pca_scores[, 1:2]))
# plot(silhouette_scores, col = 2:max(pca_scores$Cluster), border = NA)
mean_silhouette <- mean(silhouette_scores[, 3])
ch_index <- index.G1(pca_scores[, 1:2], pca_scores$Cluster)
db_index <- index.DB(pca_scores[, 1:2], pca_scores$Cluster)$DB

print(paste("Average Silhouette Score:", round(mean_silhouette, 3)))
print(paste("Calinski-Harabasz Index:", round(ch_index, 3)))
print(paste("Davies-Bouldin Index:", round(db_index, 3)))
# ------------------------------- Site by site ---------------------------------
# Perform PCA for each site individually
site_pca_results <- list()
site_pca_scores <- list()

for (site in sites) {
  # Filter data for the site
  site_data <- all_predictors_filtered[[site]]
  site_data <- site_data[, selected_columns]
  
  # Perform PCA
  site_pca <- prcomp(site_data, center = TRUE, scale. = TRUE)
  site_pca_results[[site]] <- site_pca
  
  # Convert PCA scores to a data frame and add site information
  site_scores <- as.data.frame(site_pca$x)
  site_scores$Site <- site
  site_pca_scores[[site]] <- site_scores
  
  # Plot site-specific scree plot
  scree_plot <- fviz_eig(site_pca, addlabels = TRUE, barfill = "steelblue", 
                         barcolor = "steelblue", linecolor = "red") +
    labs(title = paste("Scree Plot for", site), x = "Principal Components", y = "Explained Variance (%)")
  print(scree_plot)
  
  # Plot site-specific biplot
  biplot <- fviz_pca_biplot(site_pca, 
                            label = "var",         # Show variable labels
                            repel = TRUE,          # Avoid overlapping labels
                            geom.ind = "point",    # Show PCA scores as points
                            pointshape = 21, 
                            pointsize = 2,
                            palette = "jco") +
    labs(title = paste("Biplot of Principal Components for", site))
  print(biplot)
  
  # Variable Contribution to PCs
  contrib_pc1 <- fviz_contrib(site_pca, choice = "var", axes = 1, top = 10) +
    labs(title = paste("Variable Contributions to PC1 for", site))
  print(contrib_pc1)
  contrib_pc2 <- fviz_contrib(site_pca, choice = "var", axes = 2, top = 10) +
    labs(title = paste("Variable Contributions to PC2 for", site))
  print(contrib_pc2)
  
  # Final PCA Plot with Clusters and Sites
  site_scores$Cluster <- kmeans(site_scores[, 1:2], centers = k)$cluster
  final_plot <- ggplot(site_scores, aes(x = PC1, y = PC2, color = as.factor(Cluster))) +
    geom_point(size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1.5) +
    labs(
      title = paste("PCA Clustering for", site),
      x = "PC1",
      y = "PC2",
      color = "Cluster",
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  print(final_plot)
}
