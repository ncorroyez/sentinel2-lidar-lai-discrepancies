# ---
# title: "0_convert_l93_into_utm.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-03-26"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

library(tools)
library(lidR)

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

# l93_init <- readLAS("/home/corroyez/Downloads/LAStools/bin/LAS_747500_7006000.las")
# utm_init <- readLAS("/home/corroyez/Downloads/LAStools/bin/utm.las")
# utm_norm <- readLAS("/media/corroyez/MyPassport/01_DATA/Aigoual/LiDAR/2-las_utm/LAS_545486.27_4885544.09.las")

lastools_path <- "/home/corroyez/Downloads/LAStools/bin/"
new_epsg <- "-target_utm 31N"
actual_epsg <- "-epsg 2154"

# Chemin vers le dossier contenant les fichiers LAS
# datapath <- "/home/corroyez/Documents/NC_Full/01_DATA"
datapath <- "/media/corroyez/MyPassport/01_DATA"
site <- "Reine" # Mormal Blois Aigoual Hayes Reine
las_l93_dir <- file.path(datapath, site, "LiDAR/leaf_off/1-las_l93")
las_l93_files <- list.files(las_l93_dir, 
                            pattern = "\\.(las|laz)$", 
                            full.names = FALSE)
# las_l93_files <- las_l93_files[-(1:408)]
las_utm_dir <- file.path(datapath, site, "LiDAR/leaf_off/2-las_utm")
dir.create(path = las_utm_dir, showWarnings = FALSE, recursive = TRUE)

# setwd(lastools_path)
# Boucle sur les fichiers LAS du dossier
# for (file in las_l93_files) {
#   # print(file.info(file.path(las_l93_dir, file)))
#   # Chemin vers le fichier LAS reprojeté
#   # filename_out <- paste0(tools::file_path_sans_ext(file), "_utm.las")
#   filename_out <- paste0(tools::file_path_sans_ext(file), ".las")
#   
#   # LASTools reprojection
#   cmd <- paste0("LD_LIBRARY_PATH=/home/corroyez/Downloads/LAStools/bin/",
#                 " ",
#                 lastools_path,
#                 "las2las64 -i ", 
#                 # paste0(las_l93_dir, file),
#                 file.path(las_l93_dir, file), 
#                 " -o ", 
#                 # paste0(las_utm_dir, filename_out),
#                 file.path(las_utm_dir, filename_out), 
#                 " ", 
#                 actual_epsg, 
#                 " ", 
#                 new_epsg)
#   system(command = cmd, wait = TRUE)
#   print(cmd)
# }

# Loop over LAS files
for (file in las_l93_files) {
  
  # 1. Define paths
  input_path <- file.path(las_l93_dir, file)
  # Create a temporary output filename first
  temp_output_path <- file.path(las_utm_dir, paste0("temp_", file))
  
  # 2. Construct LAStools command
  cmd <- paste0("LD_LIBRARY_PATH=/home/corroyez/Downloads/LAStools/bin/ ",
                lastools_path, "las2las64 -i ", 
                input_path, 
                " -o ", 
                temp_output_path, 
                " ", 
                actual_epsg, 
                " ", 
                new_epsg)
  
  # 3. Execute Reprojection
  print(paste("Reprojecting:", file))
  system(command = cmd, wait = TRUE)
  
  # 4. Read header of the NEW file to get UTM coordinates
  # We use readLASheader because it is much faster than reading the whole LAS
  if (file.exists(temp_output_path)) {
    header <- lidR::readLASheader(temp_output_path)
    
    # Get the minimum X and Y (usually used for tiling names)
    # We round them to avoid decimals in filenames
    new_x <- header@PHB[["Min X"]]
    new_y <- header@PHB[["Min Y"]]
    
    # Construct new filename: e.g., LAS_350000_5600000.las
    new_filename <- paste0("LAS_", new_x, "_", new_y, ".las")
    final_output_path <- file.path(las_utm_dir, new_filename)
    
    # 5. Rename the file
    file.rename(from = temp_output_path, to = final_output_path)
    print(paste("--> Renamed to:", new_filename))
  } else {
    warning(paste("Failed to create", temp_output_path))
  }
}

# Renommer les fichiers UTM et supprimer le suffixe "_utm"
# utm_files <- list.files(las_utm_dir, pattern = "\\.las$", full.names = TRUE)
# for (file in utm_files) {
#   # Extraire les coordonnées du nom de fichier UTM
#   coords <- gsub(".las", "", basename(file))
#   coords <- unlist(strsplit(coords, "_"))
#   
#   x_coord <- coords[2]
#   y_coord <- coords[3]
#   
#   # Nouveau nom de fichier sans le suffixe "_utm"
#   new_filename <- paste0("LAS_", x_coord, "_", y_coord, ".las")
#   
#   # Renommer le fichier
#   file.rename(from = file, to = file.path(las_utm_dir, new_filename))
# }






# datapath <- "../../01_DATA"
# 
# las_dir <- file.path(datapath, "Mormal/LiDAR/1-las_l93")
# cat("lasDpath: ", las_dir, "\n")
# 
# las_utm_norm_dir <- file.path(datapath, "Mormal/LiDAR/3-las_normalize_utm")
# las_utm_norm_files <- list.files(path = las_utm_norm_dir, pattern = "\\.las$", 
#                                  full.names = TRUE)
# las_utm_norm_test <- readLAS(las_utm_norm_files[1])
# 
# # Define the output directory for the reprojected LAS files
# las_dir_utm <- file.path(datapath, "Mormal/LiDAR/2-las_utm")
# dir.create(las_dir_utm, showWarnings = FALSE)
# 
# # Define the transformation function
# convert_coordinates <- function(x, y) {
#   utm <- st_transform(st_point(cbind(x, y)), 32631)  # Transform to UTM Zone 31N
#   return(c(st_coordinates(utm)[1], st_coordinates(utm)[2]))  # Extract transformed coordinates
# }
# proj4string(las_l93) <- CRS("+init=epsg:2154")
# # Apply the transformation
# las_l93@data$XY <- convert_coordinates(las_l93@data$X, las_l93@data$Y)
# 
# # Update the header with the new coordinate system information
# las_l93@header@crs$wkt <- 'PROJCS["WGS 84 / UTM zone 31N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",3],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","32631"]]'
# 
# # Write the LAS file
# writeLAS(las_l93, "output_file_utm.las")
# 
# # Define the EPSG code for the UTM zone you want to project to
# # For example, for UTM zone 31N, the EPSG code is 32631
# utm_epsg <- 32631  # Change this to the appropriate EPSG code for your UTM zone
# 
# # Function to reproject a single LAS file
# reproject_las <- function(las_file, utm_epsg) {
#   las <- readLAS(las_file)
#   las <- filter_poi(las, !is.na(Z))  # Remove points with NA elevation
#   las@header@PHB$EPSG <- utm_epsg  # Update EPSG code in the header
#   return(las)
# }
# 
# # Get a list of LAS files in the directory
# las_files <- list.files(path = las_dir, pattern = "\\.las$", full.names = TRUE)
# 
# # Process each LAS file
# for (las_file in las_files) {
#   # Define output file name in the output directory
#   output_file <- file.path(las_dir_utm, basename(sub(".las$", "_utm.las", las_file)))
#   
#   # Reproject the LAS file
#   las_utm <- reproject_las(las_file, utm_epsg)
#   
#   # Write the reprojected LAS file
#   writeLAS(las_utm, output_file)
#   
#   cat("LAS file", las_file, "converted and saved as", output_file, "\n")
# }
