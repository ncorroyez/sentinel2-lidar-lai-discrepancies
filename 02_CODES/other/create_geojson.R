# ---
# title: "create_geojson.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2023-12-15"
# ---

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

library(sf)
library(tidyr)

# --------------------------------------------- data_plot.csv --------------------------------------------------------

data <- read.csv("../../../MaCCMic_Imprint/Data_placettes_Aigoual/data_plots_202402281200.csv", sep = ",")

data$coord_x_l93 <- as.numeric(gsub(",", ".", gsub("�", "", data$coord_x_l93)))
data$coord_y_l93 <- as.numeric(gsub(",", ".", gsub("�", "", data$coord_y_l93)))
data$long_wgs84 <- as.numeric(gsub(",", ".", data$long_wgs84))
data$lat_wgs84 <- as.numeric(gsub(",", ".", data$lat_wgs84))
data$elevation <- as.numeric(gsub(",", ".", data$elevation))

data_blois <- data[data$forest == "Blois", ]
data_mormal <- data[data$forest == "Mormal", ]
data_aigoual <- data[data$forest == "Aigoual", ]

coordinates_blois <- st_as_sf(data_blois, 
                              coords = c("long_wgs84", "lat_wgs84"), 
                              crs = 4326)
coordinates_mormal <- st_as_sf(data_mormal, 
                               coords = c("long_wgs84", "lat_wgs84"), 
                               crs = 4326) 
coordinates_aigoual <- st_as_sf(data_aigoual, 
                                coords = c("long_wgs84", "lat_wgs84"), 
                                crs = 4326) 

st_write(coordinates_blois, "data_blois_l93.geojson", driver = "GeoJSON", append = T)
st_write(coordinates_mormal, "data_mormal_l93.geojson", driver = "GeoJSON", append = T)
st_write(coordinates_aigoual, "data_aigoual_l93.geojson", driver = "GeoJSON", append = T)

# Transform to UTM Zone 31N (EPSG: 32631)
utm_blois <- st_transform(coordinates_blois, crs = 32631)
utm_mormal <- st_transform(coordinates_mormal, crs = 32631)
utm_aigoual <- st_transform(coordinates_aigoual, crs = 32631)

# Extract UTM coordinates and add new columns
utm_blois <- utm_blois %>%
  mutate(coord_x_utm31n = st_coordinates(.)[, 1],
         coord_y_utm31n = st_coordinates(.)[, 2]) %>%
  select(-coord_x_l93, -coord_y_l93)  # Remove L93 coordinates if they exist

utm_mormal <- utm_mormal %>%
  mutate(coord_x_utm31n = st_coordinates(.)[, 1],
         coord_y_utm31n = st_coordinates(.)[, 2]) %>%
  select(-coord_x_l93, -coord_y_l93)

utm_aigoual <- utm_aigoual %>%
  mutate(coord_x_utm31n = st_coordinates(.)[, 1],
         coord_y_utm31n = st_coordinates(.)[, 2]) %>%
  select(-coord_x_l93, -coord_y_l93)

# Save the transformed coordinates as GeoJSON
st_write(utm_blois, "data_blois_utm31n.geojson", driver = "GeoJSON", append = T)
st_write(utm_mormal, "data_mormal_utm31n.geojson", driver = "GeoJSON", append = T)
st_write(utm_aigoual, "data_aigoual_utm31n.geojson", driver = "GeoJSON", append = T)

# --------------------------------------------------------------------------------------------------------------------

# ------------------------------------------- data_lidar_pad.csv -----------------------------------------------------

# data <- read.csv("../../MaCCMic_Imprint/MaccMIC_dec2023/data_lidar_pad.csv", sep = "\t")
# 
# data$coord_x_l93 <- as.numeric(gsub(",", ".", gsub("�", "", data$coord_x_l93)))
# data$coord_y_l93 <- as.numeric(gsub(",", ".", gsub("�", "", data$coord_y_l93)))
# data$long_wgs84 <- as.numeric(gsub(",", ".", data$long_wgs84))
# data$lat_wgs84 <- as.numeric(gsub(",", ".", data$lat_wgs84))
# data$elevation <- as.numeric(gsub(",", ".", data$elevation))
# 
# data_blois <- data[data$forest == "Blois", ]
# data_mormal <- data[data$forest == "Mormal", ]
# 
# coordinates_blois <- st_as_sf(data_blois, 
#                               coords = c("long_wgs84", "lat_wgs84"), 
#                               crs = 4326)
# coordinates_mormal <- st_as_sf(data_mormal, 
#                                coords = c("long_wgs84", "lat_wgs84"), 
#                                crs = 4326) 
# 
# st_write(coordinates_blois, "data_blois.geojson", driver = "GeoJSON")
# st_write(coordinates_mormal, "data_mormal.geojson", driver = "GeoJSON")
