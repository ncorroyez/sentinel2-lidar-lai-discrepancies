# ---
# title: "main_stands.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-10-16"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

library(terra)

# Pre-processing Parameters
sites <- c("Aigoual", "Blois", "Mormal")
results_dir <- "../03_RESULTS"

# General Type
for (site in sites){
  masks_path <- file.path(results_dir, site, "LiDAR/Heterogeneity_Masks")
  metrics_dir <- file.path(results_dir, site, "Metrics/Not_Masked") # Not_Masked Deciduous_Only
  # metrics_dir <- file.path(results_dir, site, "Metrics/Deciduous_Only") # Not_Masked Deciduous_Only
  
  # Open the shapefile
  full_comp <- vect(file.path(masks_path, "site_full_composition.shp"))
  # full_comp <- vect(file.path(masks_path, "site_deciduous_only.shp"))
  raster_1m <- rast(file.path(masks_path, "chm.envi"))
  raster_10m <- rast(file.path(metrics_dir, "lidarlai_res_10_m.tif"))
  full_comp <- project(full_comp, crs(raster_10m))
  
  # Get the TFV_G11 attribute
  tfv_values <- full_comp$TFV_G11
  # Create a new column 'Stand_Type' with classified stand types
  full_comp$Stand_Type <- ifelse(grepl("feuillus|peupleraie", tfv_values, ignore.case = TRUE), "deciduous",
                                 ifelse(grepl("conifères", tfv_values, ignore.case = TRUE), "coniferous",
                                        ifelse(grepl("mixte", tfv_values, ignore.case = TRUE), "mixtures",
                                               ifelse(grepl("lande|formation herbacée|sans couvert", tfv_values, ignore.case = TRUE), "land", NA))))
  full_comp$Stand_Type <- factor(full_comp$Stand_Type, levels = c("deciduous", "coniferous", "mixtures", "land"))
  stand_type_raster <- rasterize(full_comp, raster_10m, field = "Stand_Type", touches=TRUE)
  # levels(stand_type_raster) <- data.frame(ID = 1:4, Class = c("deciduous", "coniferous", "mixtures", "land"))
  # Convert the raster to a factor with defined levels
  # levels(stand_type_raster) <- data.frame(ID = 1:4, Class = c("deciduous", "coniferous", "mixtures", "land"))
  # stand_type_raster <- project(stand_type_raster, raster_10m)
  
  # Define the output file path
  output_raster_path <- file.path(metrics_dir, "stand_type_raster.tif")
  # plot(stand_type_raster)
  # Save the raster as a GeoTIFF, keeping the categorical information
  writeRaster(stand_type_raster, output_raster_path, overwrite = TRUE, datatype = "INT1U")
  
  print(paste("Categorical stand type raster saved for site:", site))
}

# Define a named vector for stand type to number mapping
stand_type_mapping <- c(
  "oak" = 0,
  "deciduous" = 1,
  "beech" = 2,
  "poplar" = 3,
  "coniferous" = 4,
  "mixed" = 5,
  "douglas" = 6,
  "larch" = 7,
  "fir/spruce" = 8,
  "scots pine" = 9,
  "laricio pine/black pine" = 10,
  "mixed pines" = 11,
  "nc" = 12,
  "nr" = 13,
  "land" = 14
)

for (site in sites){
  masks_path <- file.path(results_dir, site, "LiDAR/Heterogeneity_Masks")
  metrics_dir <- file.path(results_dir, site, "Metrics/Not_Masked") # Not_Masked Deciduous_Only
  # metrics_dir <- file.path(results_dir, site, "Metrics/Deciduous_Only") # Not_Masked Deciduous_Only
  
  # Open the shapefile
  full_comp <- vect(file.path(masks_path, "site_full_composition.shp"))
  # full_comp <- vect(file.path(masks_path, "site_deciduous_only.shp"))
  raster_1m <- rast(file.path(masks_path, "chm.envi"))
  raster_10m <- rast(file.path(metrics_dir, "lidarlai_res_10_m.tif"))
  full_comp <- project(full_comp, crs(raster_10m))
  
  # Get the TFV_G11 attribute
  tfv_values <- full_comp$ESSENCE
  # print(tfv_values)
  # Create a new column 'Stand_Type' with classified stand types
  full_comp$Stand_Type <- ifelse(grepl("Chênes décidus", tfv_values, ignore.case = TRUE), "oak",
                                 ifelse(grepl("Feuillus", tfv_values, ignore.case = TRUE), "deciduous",
                                        ifelse(grepl("Hêtre", tfv_values, ignore.case = TRUE), "beech",
                                               ifelse(grepl("Peuplier", tfv_values, ignore.case = TRUE), "poplar",
                                                      ifelse(grepl("Conifères", tfv_values, ignore.case = TRUE), "coniferous",
                                                             ifelse(grepl("Mixte", tfv_values, ignore.case = TRUE), "mixed",
                                                                    ifelse(grepl("Douglas", tfv_values, ignore.case = TRUE), "douglas",
                                                                           ifelse(grepl("Mélèze", tfv_values, ignore.case = TRUE), "larch",
                                                                                  ifelse(grepl("Sapin, épicéa", tfv_values, ignore.case = TRUE), "fir/spruce",
                                                                                         ifelse(grepl("Pin sylvestre", tfv_values, ignore.case = TRUE), "scots pine",
                                                                                                ifelse(grepl("Pin laricio, pin noir", tfv_values, ignore.case = TRUE), "laricio pine/black pine",
                                                                                                       ifelse(grepl("Pins mélangés", tfv_values, ignore.case = TRUE), "mixed pines",
                                                                                                              ifelse(grepl("nc", tfv_values, ignore.case = TRUE), "nc",
                                                                                                                     ifelse(grepl("nr", tfv_values, ignore.case = TRUE), "nr",
                                                                                                                            ifelse(grepl("lande|formation herbacée|sans couvert", tfv_values, ignore.case = TRUE), "land", NA)))))))))))))))
  # full_comp$Stand_Type <- factor(full_comp$Stand_Type, 
  #                                levels = sort(c("oak", "deciduous", "beech", "poplar", "coniferous", 
  #                                                "mixed", "douglas", "larch", "fir/spruce", "scots pine", 
  #                                                "laricio pine/black pine", "mixed pines", "nc", "nr", "land")),
  #                                labels = sort(c("oak", "deciduous", "beech", "poplar", "coniferous", 
  #                                                "mixed", "douglas", "larch", "fir/spruce", "scots pine", 
  #                                                "laricio pine/black pine", "mixed pines", "nc", "nr", "land")))
  full_comp$Stand_Type_Number <- stand_type_mapping[full_comp$Stand_Type]
  stand_type_raster <- rasterize(full_comp, raster_10m, field = "Stand_Type_Number", touches=TRUE)
  # levels(stand_type_raster) <- data.frame(ID = 1:15, Class = c("oak", "deciduous", "beech", "poplar", "coniferous", "mixed", "douglas", "larch", 
                                                               # "fir/spruce", "scots pine", "laricio pine/black pine", "mixed pines", "NC", "NR", "land"))
  # stand_type_raster <- project(stand_type_raster, raster_10m)
  
  # Define the output file path
  output_raster_path <- file.path(metrics_dir, "stand_trees_raster.tif")
  
  # Save the raster as a GeoTIFF, keeping the categorical information
  writeRaster(stand_type_raster, output_raster_path, overwrite = TRUE, datatype = "INT1U")
  plot(stand_type_raster)
  print(paste("Precise stand type raster saved for site:", site))
}
# r <- rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/Metrics/Deciduous_Only/stand_trees_raster.tif")
