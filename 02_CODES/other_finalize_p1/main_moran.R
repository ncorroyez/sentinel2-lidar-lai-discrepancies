# ---
# title: "main_moran.R"
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

# Import libraries
library(terra)
library(ggplot2)
library(dplyr)
library(spdep)

# Pre-processing Parameters
# sites <- "Aigoual" # Mormal Blois Aigoual
# results_dir <- "../03_RESULTS"
# forest_composition <- "Not_Masked" # Full_Composition Deciduous_Only Not_Masked
# metrics_dir <- file.path("Metrics", forest_composition)


sites <- c("Aigoual", "Blois", "Mormal")


# site <- "Aigoual"
for (site in sites){
  # lai_raster <- terra::rast(file.path(results_dir, site, metrics_dir, "lidarlai_res_10_m.tif"))
  lai_raster <- terra::rast(file.path(site, "lidarlai_res_10_m.tif"))
  
  # Randomly select a pixel
  # set.seed(42)  # For reproducibility
  # Ensure the pixel chosen has a valid value
  valid_pixel <- FALSE
  while (!valid_pixel) {
    # Randomly pick a pixel
    pixel_index <- sample(ncell(lai_raster), 1)
    pixel_value <- values(lai_raster)[pixel_index]
    
    # Check if the value is not NA
    if (!is.na(pixel_value)) {
      valid_pixel <- TRUE
      center_coords <- xyFromCell(lai_raster, pixel_index)  # Get spatial coordinates
    }
  }
  
  # Initialize parameters
  max_window_size <- 5  # Maximum size of the growing window (must be odd)
  moran_results <- data.frame(window_size = integer(), moran_i = numeric())
  
  # Loop through increasing window sizes
  for (window_size in seq(3, max_window_size, by = 2)) {
    half_window <- window_size / 2
    
    # Calculate spatial extent for the current window size
    res_x <- res(lai_raster)[1]
    res_y <- res(lai_raster)[2]
    ext_window <- ext(
      center_coords[1] - half_window * res_x,
      center_coords[1] + half_window * res_x,
      center_coords[2] - half_window * res_y,
      center_coords[2] + half_window * res_y
    )
    
    # Extract the window and its values using the calculated extent
    window <- crop(lai_raster, ext_window)
    values <- values(window, na.rm = TRUE)
    coords <- xyFromCell(window, seq_len(ncell(window)))
    
    # Plot window (optional)
    window_df <- data.frame(x = coords[, 1], y = coords[, 2], V1 = values)
    center_pixel <- data.frame(x = center_coords[1], y = center_coords[2])
    win_plot <- ggplot(window_df, aes(x = x, y = y, fill = V1)) +
      geom_tile() +
      geom_point(data = center_pixel, aes(x = x, y = y), 
                 color = "red", size = 8, shape = 21, fill = "red") +
      scale_fill_viridis_c(na.value = "transparent") +
      labs(title = paste(window_size, "x", window_size, 
                         "Window", "with Center Pixel Highlighted"), 
           x = "X Coordinate", y = "Y Coordinate", fill = "Value") +
      theme_minimal()
    print(win_plot)
    
    # Define neighbors
    # nb <- dnearneigh(coords, 0, max(dist(coords))) # problem when taking all neighbors
    
    # Definition simple du nombre de voisins
    # en attendant de reussir a se focaliser sur le voisin central uniquement
    k <- floor(window_size*window_size - 1) / 2
    nb <- knn2nb(knearneigh(coords, k = k))
    
    # Calculate associated spatial matrix weights
    lw <- nb2listw(nb, style = "W")
    
    # Look at lw's structure (all neighbors lists of all neighbors)
    # print(str(lw))
    
    # Calculate Moran's I
    moran_test <- moran.test(values, lw)
    moran_i <- moran_test
    print(moran_i)
    
    # Store Moran Evolution
    moran_results <- rbind(moran_results, 
                           data.frame(window_size = window_size, 
                                      moran_i = moran_i$estimate[1]))
  }
  
  # Plot Moran Evolution
  plot <- ggplot(moran_results, aes(x = window_size, y = moran_i)) +
    geom_line(color = "blue") +
    geom_point(color = "red") +
    labs(
      title = paste(site, "Evolution of Moran's I with Growing Window Size"),
      x = "Window Size (Edge Length)",
      y = "Moran's I"
    ) +
    theme_minimal()
  print(plot)
  
  
  # Optional: look at the chosen pixel
  # lai_raster_df <- as.data.frame(lai_raster, xy = TRUE, na.rm = TRUE)
  # names(lai_raster_df) <- c("x", "y", "value")
  # ggplot(lai_raster_df, aes(x = x, y = y, fill = value)) +
  #   geom_tile() +
  #   scale_fill_viridis_c(na.value = "transparent") +  # Use a color scale
  #   geom_point(aes(x = center_coords[1], y = center_coords[2]), 
  #              color = "red", size = 3) +  # Overlay the chosen pixel
  #   labs(title = "Chosen Pixel on Raster", fill = "Value") +
  #   theme_minimal()
}