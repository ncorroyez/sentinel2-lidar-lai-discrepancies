# ------------------------------------------------------------------------------------------------ #
# Script to extract image information at point location for large point data sets (GEDI) 

# Example for Vosges data GEDI v2 (complete data of first mission session (2019-2023))

# ------------------------------------------------------------------------------------------------ #
# Authors: Sylvie Durrieu 
# Date creation : oct 2023
# ------------------------------------------------------------------------------------------------ #


#-------------------------------------------------------------------------------
# load libraries
#-------------------------------------------------------------------------------
library(terra)
library(sf)
library(bit64)
library(data.table)
library(purrr)
library(stringr)
library(tidyverse)


rm(list=ls())     # remove objects from the environment
gc()  


#------------------------------------------
# Images from which is information is to be extracted
#------------------------------------------

# image folder 

img_folder <- "D:Data\\images_S2\\"    # possible to put all the images in one folder; images can be renamed before processing

#------------------------------------------
# List of image files to be processed
#------------------------------------------
  
# if all the images are in the same folder 
imgs <- paste0(img_folder, list.files(img_folder, "*.tif"))

# optional; sort images by name
sort(imgs)
# optional; organize list of images according to the preferred order in the final data frame; 
# exclude images if necessary
imgs <- c(imgs[13:22], imgs[c(5,10:12)], imgs[7:9], imgs[1:4], imgs[6])

# gather images in a collection 

imgs_col <- sprc(imgs)

### if gedi is not a spatial object

gedi <- sf::st_as_sf(x = gedi, coords = c("lon1", "lat1"), crs = 32632)

#------------------------------------------
# For small point data sets, extract image information at point location with a loop 
# example with NFI data set
#------------------------------------------


#-------------------------------------
# load and prepare NFI data set
#-------------------------------------
# NFI file

nfi <- readRDS(file.choose())
# nfi <- readRDS(file.choose())\\data_Anouk\\old\\Grid30_sentinelIFN_202208_06_lai.rds")

# optional: store point coordinates to add them in the final table

lat <- nfi$lat1
lon <- nfi$lon1

# if the file is not a spatial object transform it to a spatial object 
nfi <- sf::st_as_sf(x = nfi, coords = c("lon1", "lat1"), crs = 32632)

#-----------------------------------------
# extract mean and sd at NFI plot with terra::zonal
#-----------------------------------------

nfi_v <- vect(nfi)  # transform the file in a terra vector object

size_buf <- 50      # define the size of the buffer

nfiv_buf <- buffer(nfi_v, width = size_buf)     # create the polygon file of the circular buffers 

data_nfi <- NULL            # initialize the output file
nom_col_nfi <- NULL         # initialize the vector of the column names for the output file

# loop on the images
t1 <- Sys.time()
for (i in 1:length(imgs_col)) {
  mean <- zonal(imgs_col[i], nfiv_buf, fun="mean", as.raster=FALSE, exact =TRUE, as.polygons = TRUE)   # extract mean
  std <- zonal(imgs_col[i], nfiv_buf, fun="sd", as.raster=FALSE, touche= TRUE)  # extract standard deviation
  std[[1]][is.na(std[[1]])] <- 0    # for Na sd values (when only one valeu to assess sd) 
  
  longchaine <- nchar(basename(sources(imgs_col[i])))    # image name
  #namefic <- substr(basename(sources(imgs_col[i])), min(longchaine, longchaine-8), longchaine-4)  # to reduce long names
  namefic <- substr(basename(sources(imgs_col[i])), 1, longchaine-4)          # image 
  nom <- c(paste0("Zonal",size_buf,"_", namefic, "_mean") , paste0("Zonal",size_buf,"_", namefic,"_sd") )     # name of the 2 columns
  
  data_nfi <- cbind(data_nfi, mean[[1]], std[[1]])
  nom_col_nfi <- c(nom_col_nfi, nom ) 
  mean <- NULL
  nom <- NULL
  print(i)
}

data_nfi <- as.data.frame(data_nfi)
colnames(data_nfi) <- nom_col_nfi
t2 <- Sys.time() -t1
t2

# add unique id to final table and point coordinates 

data_nfi$npp <- nfi$npp
data_nfi$lon1 <- lon
data_nfi$lat1 <- lat

# save file
saveRDS(data_nfi, "D:\\data\\data_Anouk\\calcul_mean_sd_S2_dif\\S2_dif_mean_sd_nfi_buf50.rds")


#------------------------------------------
# For large point data sets, extract image information at point location with parallelization 
# Example with NFI data set
#------------------------------------------


library(foreach)
library(doParallel)

# parallel::detectCores()

# function to process a group of data
# name of output columns can be changed in the function
# current names = Zonal + buf_size + imagefilename + function used in zonal stat

process_poly <- function(group, buf_size) {
  #-----------------------------------------
  # Libraries
  #-----------------------------------------
  library(terra)
  library(sf)
  library(bit64)
  library(data.table)
  library(purrr)
  library(stringr)
  library(tidyverse)
  
  #------------------------------------------
  # List of image files to be processed
  #------------------------------------------
  
  # if all the images are in the same folder 
  imgs <- paste0(img_folder, list.files(img_folder, "*.tif"))
  # optional: sort images by names 
  sort(imgs)
  #imgs <- c(imgs[13:22], imgs[c(5,10:12)], imgs[7:9], imgs[1:4], imgs[6])
  
  imgs_col <- sprc(imgs) 
  
  poly <- gedi[which(gedi$group==group),]    # select polygons of the group being processed
  
  gedi_v <- vect(poly)                            # transform to terra vector
  gediv_buf <- buffer(gedi_v, width = buf_size)   # create a vector of buffers (polygons)
  
  data <- NULL    # initialize output data frame
  nom_col <- NULL  # initialize column names for output data frame
  
  # loop on images
  for (i in 1:length(imgs_col)) {
      mean <- zonal(imgs_col[i], gediv_buf, fun="mean", as.raster=FALSE, exact =TRUE, as.polygons = TRUE)  # compute 
      std <- zonal(imgs_col[i], gediv_buf, fun="sd", as.raster=FALSE, touche= TRUE) 
      std[[1]][is.na(std[[1]])] <- 0
      longchaine <- nchar(basename(sources(imgs_col[i])))
      #namefic <- substr(basename(sources(imgs_col[i])), min(longchaine, longchaine-10), longchaine-4) # if long image names,  to take only part of the name
      namefic <- substr(basename(sources(imgs_col[i])), 1, longchaine-4)
      nom <- c(paste0("Zonal",buf_size,"_", namefic, "_mean") , paste0("Zonal",buf_size,"_", namefic,"_sd") )
      
      data <- cbind(data, mean[[1]], std[[1]])
      nom_col <- c(nom_col, nom ) 
      mean <- NULL
      nom <- NULL
      print(i)
    }
  colnames(data) <- nom_col
  data <- as.data.frame(data)
  data$shotnumber <- poly$shotnumber          # add unique id
  data$lon <- unlist(map(poly$geometry,1))    # add point coordinates
  data$lat <- unlist(map(poly$geometry,2))
  
  return(data)
}

# Image folder  

img_folder <- "D:Data\\images_all\\"     # folder containing the images to be processed

# laod GEDI data set with at least footprint location and shotnumber as unique id 
gedi <- readRDS(file.choose())       

## if gedi file is not a spatial object transform it to a spatial object
# gedi <- sf::st_as_sf(x = gedi, coords = c("lon", "lat"), crs = 32632)
# gedi <- sf::st_as_sf(x = gedi, coords = c("lonrecalc", "latrecalc"), crs = 32632)

# create groups of polygons

nb <- 1000    # number of polygons in each group
gedi$group <- floor(as.numeric(rownames(gedi))/ nb)+1     # create a variable "group" in gedi file


# Create a cluster for parallel execution (adjust the number of cores as needed)
cl <- makeCluster(4)  # Use 4 CPU cores, you can adjust this based on your system

# Register the cluster for parallel execution
registerDoParallel(cl)

# Initialize the output file  as a list
all_mean_sd_geodGedi_buf_50 <- NULL

groups <- 1:max(gedi$group)

t1 <- Sys.time()

all_mean_sd_geodGedi_buf_50 <- foreach(group = groups, .packages=c('terra','bit64','data.table','sf', 'purrr', 'stringr', 'tidyverse'), .combine = "rbind") %dopar% {
  
  data <- process_poly(group, 50)
  return(data)
  }

# Stop the cluster
stopCluster(cl)

# Convert the result to a data frame
all_mean_sd_geodGedi_buf_50 <- as.data.frame(all_mean_sd_geodGedi_buf_50) 

t2 <- Sys.time()
t2-t1

# save the final data frame 
saveRDS(all_mean_sd_geodGedi_buf_50, "D:\\data\\data_Anouk\\calcul_mean_sd_geogedi\\all_mean_sd_geodGedi_buf_50.rds")



