# ---
# title: "plot_PAD_profiles.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-12-02"
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
library("fields")

# --------------------------- Import useful functions --------------------------

source("libraries/functions_plots.R")

# --------------------------------- Setup --------------------------------------
results_dir <- "../03_RESULTS"
figures_dir <- "../04_FIGURES"
metrics_dir <- "Metrics/Deciduous_Only"
sites <- c("Aigoual", "Blois", "Mormal")
pad_types <- c("PAD_Profiles", "PAD_Profiles_updated_modifminz")

# for (site in sites) {
#   # Predefine the pixel index for the current site
#   site_raster_files <- list.files(file.path(results_dir, site, metrics_dir, "PAD_Profiles_updated_modifminz"), 
#                                   pattern = "\\.tif$", 
#                                   full.names = TRUE, 
#                                   recursive = TRUE)
#   # Extract depths from file names
#   depths <- sapply(basename(site_raster_files), function(file) {
#     parts <- unlist(strsplit(file, "_|\\."))
#     as.numeric(parts[4]) - as.numeric(parts[2]) # Calculate depth
#   })
#   
#   # Combine files and depths into a data frame
#   file_depth_df <- data.frame(
#     file = site_raster_files,
#     depth = depths
#   )
#   
#   # Sort by descending depth
#   file_depth_df <- file_depth_df[order(file_depth_df$depth, decreasing = TRUE), ]
#   
#   # Load the sorted raster stack
#   site_raster_stack <- rast(file_depth_df$file)
#   
#   # Extract all pixel values into a matrix
#   site_values_matrix <- values(site_raster_stack)
#   
#   # Identify a valid pixel index
#   site_valid_pixel_index <- which(
#     rowSums(!is.na(site_values_matrix)) == nlyr(site_raster_stack) & # Valid in all layers
#       site_values_matrix[, 1] > site_values_matrix[, 2] + 0.1            # max_depth != max_depth - 1
#   )[1] # Select the first valid pixel
#   
#   # Extract pixel values for the selected pixel
#   selected_pixel_values <- site_values_matrix[site_valid_pixel_index, ]
#   
#   for (pad_type in pad_types) {
#     figures_path <- file.path(figures_dir, site, pad_type)
#     if (!dir.exists(figures_path)) {
#       dir.create(figures_path, showWarnings = FALSE, recursive = TRUE)
#     }
#     
#     pad_path <- file.path(results_dir, site, metrics_dir, pad_type)
#     raster_files <- list.files(pad_path, pattern = "\\.tif$", full.names = TRUE)
#     raster_stack <- rast(raster_files)
#     
#     depths <- sapply(basename(raster_files), function(file) {
#       parts <- unlist(strsplit(file, "_|\\."))
#       as.numeric(parts[4]) - as.numeric(parts[2])
#     })
#     
#     # Ensure depths are sorted along with the stack
#     sorted_indices <- order(depths, decreasing = TRUE)
#     depths <- depths[sorted_indices]
#     raster_stack <- raster_stack[[sorted_indices]]
#     
#     # Extract pixel values for the predefined pixel index of the current site
#     pixel_values <- values(raster_stack)[site_valid_pixel_index, ]
#     
#     plot(pixel_values, depths, type = "b", xlab = "Pixel Value", ylab = "Depth",
#          main = paste(site, "Pixel Value vs. Depth"))
#     
#   }
# }

allometry1 <- data.frame(
  height = seq(30, 0, by = -1),
  density = rev(c(0.03, 0.08, 0.15, 0.2, 0.27, 0.33, 0.36, 0.41, 0.45,
              0.43, 0.4, 0.42, 0.34, 0.3, 0.27, 0.31, 0.29,
              0.27, 0.24, 0.22, 0.15, 0.12, 0.1, 0.12, 
              0.14, 0.12, 0.1, 0.08, 0.05, 0.1, 0.4))
)

allometry2 <- data.frame(
  height = seq(30, 0, by = -1),
  density = rev(c(0.03, 0.1, 0.2, 0.3, 0.4, 0.38, 0.35, 0.4, 0.45,
              0.38, 0.43, 0.38, 0.33, 0.3, 0.35, 0.42,
              0.29, 0.22, 0.18, 0.14, 0.11, 0.16, 0.13,
              0.10, 0.08, 0.06, 0.04, 0.06, 0.08, 0.1, 0.3))
)

plot_allometry_with_colorbar(allometry1, 
                             title = "DTM Normalization",
                             key_point = "Key Point", 
                             x_position = 3, 
                             y_position = 8,
                             label_size = 2)
plot_allometry_with_colorbar(allometry2, 
                             title = "DSM Normalization", 
                             key_point = "Key Point", 
                             x_position = 3, 
                             y_position = 8,
                             label_size = 2)

# colors <- viridis(length(allometry1$height), option = "G")
# green_palette <- colorRampPalette(c("lightgreen", "green", "darkgreen"))(length(allometry1$height))
# plot(allometry1$density, allometry1$height, type = "n",
#      xlab = "density", ylab = "height", main = "Height vs Density",
#      ylim = c(0, max(allometry1$height)))
# # lines(allometry1$density, allometry1$height)
# for (i in 1:(length(allometry1$height)-1)) {
#   segments(allometry1$density[i], allometry1$height[i], 
#            allometry1$density[i+1], allometry1$height[i+1], 
#            col = green_palette[i], lwd = 2)  # Adjust line width and color
# }

# green_palette <- colorRampPalette(c("lightgreen", "green", "darkgreen"))(length(allometry2$height))
# green_palette <- rev(colorRampPalette(c("lightgreen", "yellowgreen", "darkgreen"))(length(allometry2$height)))
# 
# # Create the plot with no points, but with color for the line
# par(mar=c(5, 5, 1, 5))
# # par(mfrow=c(1,1))
# plot(allometry2$density, allometry2$height, type = "n",  # 'n' means no points
#      xlab = "Density", ylab = "Canopy Depth",
#      # ylim = c(0, max(allometry2$height)),
#      ylim = c(max(allometry2$height), 0),
#      axes = T,   # Remove default axes
#      bty = "n",      # Remove border
#      col.main = "black",   # Title color
#      col.lab = "black",    # Label color
#      cex.main = 1.2,       # Title size
#      cex.lab = 1.1)        # Label size
# text(3, 8, "Key Point", col = "blue", cex = 1.2)
# # Add the line with the custom green color gradient based on height
# for (i in 1:(length(allometry2$height)-1)) {
#   segments(allometry2$density[i], allometry2$height[i], 
#            allometry2$density[i+1], allometry2$height[i+1], 
#            col = rev(green_palette[i]), lwd = 2)  # Adjust line width and color
# }
# 
# image.plot(
#   x = c(1.1, 1.2),  # Position of the color bar
#   y = c(min(allometry2$height), max(allometry2$height)),  # Y range matching height
#   z = matrix(seq(min(allometry2$height), max(allometry2$height), length.out = 256), ncol = 1),  # Value range for the color scale
#   col = green_palette,  # Use the custom green color palette
#   legend.only = TRUE,  # Only show the legend (color bar)
#   legend.shrink = 0.9,  # Shrink the legend
#   legend.width = 1.2,   # Set width of the color bar
#   legend.mar = 4,       # Set margin for the color bar
#   horizontal = F,    # Set horizontal color bar
#   axis.args = list(cex.axis = 1, tck = -0.02, labels = F)
# )
# mtext("Cumulative Profiles", side=4)