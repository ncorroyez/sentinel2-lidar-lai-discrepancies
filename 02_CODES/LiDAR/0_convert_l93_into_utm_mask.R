# ---
# title: "0.convert_l93_into_utm_mask.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-05-21"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Load the required libraries
library(sf)

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# Directory setup
sites <- c("Aigoual", "Blois", "Mormal") # Mormal Blois Aigoual
# sites <- "Blois"
# sites <- c("Hayes", "Reine")

for (site in sites){
  input_dir <- file.path('/media/corroyez/MyPassport/01_DATA',
                         site,
                         "Geo_Files")
  input_shapefile <- file.path(input_dir, "bdforet_2.shp")
  output_shapefile <- file.path(input_dir, "bdforet_2_utm.gpkg")
  
  # Read the shapefile
  shp <- st_read(input_shapefile)
  
  # Check the current CRS of the shapefile
  original_crs <- st_crs(shp)
  cat("Original CRS:\n")
  print(original_crs)
  
  # Define the target CRS (UTM 31N)
  target_crs <- st_crs(32631) # EPSG code for UTM 31N
  
  # Transform the CRS from Lambert 93 to UTM 31N
  shp_transformed <- st_transform(shp, target_crs)
  
  # Check the new CRS of the transformed shapefile
  cat("Transformed CRS:\n")
  print(st_crs(shp_transformed))
  
  # Save the transformed shapefile
  st_write(shp_transformed, output_shapefile, append = F, overwrite = T)
  cat(site, "Shapefile successfully transformed and saved to UTM 31N.\n")
}
