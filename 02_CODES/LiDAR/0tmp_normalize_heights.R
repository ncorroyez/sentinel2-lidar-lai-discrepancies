# ---
# title: "0.normalize_heights.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-03-15"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

library("lidR")
library("raster")
library("plotly")
library("terra")
library("viridis")
library("future")
library("future.apply")
library("rgl")
library("dplyr")
library("stringr")
library("rlang")

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

# Import useful functions
source("../libraries/functions_lidar.R")
source("../libraries/functions_plots.R")
source("../libraries/functions_chm.R")

# Parallel
# plan(multisession)

# Directories
data <- "/media/corroyez/My Passport/01_DATA"
sites <- c("Aigoual") # Mormal Blois Aigoual
# sites <- c("Aigoual", "Blois", "Mormal")
# sites <- c("Aigoual", "Blois")
# sites <- c("Mormal")
# sites <- c("Blois")
lidar <- "LiDAR"
results <- "../../03_RESULTS"
# max_chm <- 40.5

plan(sequential)
plan(multisession(workers = 4))
for (site in sites){
  
  s2_rast <- switch(site,
                    "Aigoual" = terra::rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Aigoual/L2A_T31TEJ_A031608_20210711T104217/Reflectance/res_10_m/L2A_T31TEJ_A031608_20210711T104217_Refl")[[1]],
                    "Blois" = terra::rast('/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/L2A_T31TCN_A031222_20210614T105443/Reflectance/res_10_m/L2A_T31TCN_A031222_20210614T105443_Refl')[[1]],
                    "Mormal" = terra::rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/Reflectance/res_10_m/L2A_T31UER_A031222_20210614T105443_Refl")[[1]]
  )
  
  # Directories setup
  data_site <- file.path(data, site)
  # UTM
  las_utm <- file.path(lidar, "2-las_utm")
  las_dir <- file.path(data_site, las_utm)
  las_utm_files <- list.files(las_dir, pattern = "\\.las$", full.names = TRUE)
  dir.create(path = las_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Normalized UTM by DTM
  las_norm_utm <- file.path(lidar, "3-las_normalized_utm")
  las_norm_dir <- file.path(data_site, las_norm_utm)
  las_norm_utm_files <- list.files(las_norm_dir, pattern = "\\.las$", full.names = TRUE)
  dir.create(path = las_norm_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Normalized UTM by DSM
  las_norm_utm_dsm <- file.path(lidar, "4-las_normalized_utm_dsm_chunked")
  las_norm_dsm_dir <- file.path(data_site, las_norm_utm_dsm)
  las_norm_utm_dsm_files <- list.files(las_norm_dsm_dir, pattern = "\\.las$", full.names = TRUE)
  dir.create(path = las_norm_dsm_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Test
  # testtt <- file.path(lidar, "b")
  # test_dir <- file.path(data_site, testtt)
  # test_files <- list.files(test_dir, pattern = "\\.las$", full.names = TRUE)
  # ctg_test <- readLAScatalog(test_dir)
  # opt_chunk_size(ctg_test) <- 0 # Processing by files
  # opt_chunk_buffer(ctg_test) <- 10
  # opt_select(ctg_test) <- "xyzc"
  # # opt_filter(ctg_test) <- "-keep_first -keep_last -drop_z_below_0"
  # opt_filter(ctg_test) <- "-drop_z_below_0 -keep_last -keep_first"
  # stop()
  
  # Catalog
  cat("lasDpath: ", las_dir, "\n")
  ctg <- readLAScatalog(las_dir)
  # stop()
  # normalized_ctg <- readLAScatalog(las_norm_dir)
  # normalized_dsm_ctg <- readLAScatalog(las_norm_dsm_dir)
  # plot(ctg)
  st_crs(ctg) <- st_crs(s2_rast)
  # st_crs(normalized_ctg) <- st_crs(s2_rast)
  # st_crs(normalized_dsm_ctg) <- st_crs(s2_rast)
  # stop()
  # DTM
  # dtm_dir <- file.path(results, site, lidar, "dtm/res_1_m")
  # dtm <- terra::rast(file.path(dtm_dir, "rasterize_terrain.vrt"))
  # dtm_files <- list.files(dtm_dir, pattern = "\\.tif$", full.names = TRUE)
  # 
  # # DSM
  # dsm_dir <- file.path(results, site, lidar, "dsm/res_1_m")
  # dsm <- terra::rast(file.path(dsm_dir, "rasterize_canopy.vrt"))
  # dsm_files <- list.files(dsm_dir, pattern = "\\.tif$", full.names = TRUE)
  
  # CHM
  # chm_dir <- file.path(results, site, lidar, "chm/res_1_m")
  # chm <- terra::rast(file.path(chm_dir, "chm.tif"))
  
  # LiDR optimization
  opt_chunk_size(ctg) <- 250 # Processing by files
  opt_chunk_buffer(ctg) <- 5
  opt_select(ctg) <- "xyzc"
  # opt_filter(ctg) <- "-keep_first -keep_last"
  # opt_filter(ctg) <- "-drop_z_below 0"
  set_lidr_threads(1)
  stop()
  # stop()
  # Test
  # mutate all of this %>% /
  # las <- readLAS(las_utm_files[[30]])
  # st_crs(las) <- st_crs(s2_rast)
  # las <- clip_circle(las, xcenter = median(las@data$X), ycenter = median(las@data$Y),
  #                    radius = 50)
  # stop()
  # las_dsm_test <- las
  # las_dsm_test$Z[las_dsm_test$Classification == 2] <- -100
  
  
  # plot(las, color = "Classification")
  
  # las_dtm <- normalize_height(las, algo = tin())
  # las_dtm <- classify_poi(las_dtm, class = int(2), poi = ~las_dtm@data$Z < 2)
  # las_dtm$Z[las_dtm$Classification == 2] <- 0
  # plot(las_dtm, color = "Classification")
  
  # las_dtm[las_dtm@data$Z < 2]@data$Classification <- 2
  # las_dtm_class2 <- classify_poi(las_dtm, class = int(2), poi = ~las_dtm@data$Z < 2)
  # plot(las_dtm_class2, color = "Classification")
  
  # lad_dtm <- cloud_metrics(las_dtm, ~myPAD_dtm(Z))
  # lad_dtm_class2 <- cloud_metrics(las_dtm_class2, ~myPAD_dtm(Z))
  
  # DSM Class 2
  # dsmm <- rasterize_canopy(las_dtm_class2,
  #                          res = 1,
  #                          algorithm = p2r()
  # )
  # las_dsm_class2 <- normalize_height(las_dtm_class2, algo = dsmm)
  # plot(las_dsm_class2, color = "Classification")
  # las_dsm_class2$Z <- las_dsm_class2$Z + 40
  
  # las_dsm_class2 <- classify_poi(las_dsm_class2, class = int(4), poi = ~las_dsm_class2@data$Z > 22)
  # las_dsm_class2$Z <- las_dsm_class2$Z + 40 - max_z
  # las_dsm_class2$Z <- las_dsm_class2$Z + (-las_dsm_class2@header$`Min Z`)
  
  # las_dsm_class2$Z[las_dsm_class2$Classification == 2] <- min(las_dsm_class2$Z[las_dsm_class2$Classification == 4])
  # las_dsm_class2$Z[las_dsm_class2$Classification == 2] <- 0
  # plot(las_dsm_class2, color = "Classification")
  
  # lad <- LAD(las_dsm_class2@data$Z)
  
  # lad_dsm_class2 <- cloud_metrics(las_dsm_class2, ~myPAD_dsm(Classification, Z))
  # lad_dsm_class2 <- cloud_metrics(las_dsm_class2, ~myPAD_dtm(Z))
  
  # DSM
  # dsmm <- rasterize_canopy(las,
  #                          res = 1,
  #                          algorithm = p2r()
  # )
  # las_dsm <- normalize_height(las, algo = dsmm)
  
  # plot(las, color = "Classification")
  # plot(las_dtm, color = "Classification")
  # plot(las_dsm, color = "Classification")
  # plot(las_dsm_class2, color = "Classification")
  # las_dsm$Z <- las_dsm$Z + 40
  # las_dsm <- classify_poi(las_dsm, class = int(2), poi = ~las_dsm@data$Z < 2)
  # las_dsm$Z[las_dsm$Classification == 2] <- 0
  # las_dsm$Z[las_dsm$Classification == 3] <- -las_dsm$Z[las_dsm$Classification == 3]
  # las_dsm$Z[las_dsm$Classification == 4] <- -las_dsm$Z[las_dsm$Classification == 4]
  # las_dsm$Z[las_dsm$Classification == 5] <- -las_dsm$Z[las_dsm$Classification == 5]
  # las_dsm$Z[las_dsm$Classification == 2] <- min(las_dsm$Z[las_dsm$Classification == 4])
  # las_dsm$Z[las_dsm$Classification == 2] <- 0
  # plot(las_dsm, color = "Classification")
  
  # max_z <- max(las_dtm$Z)
  # las_dtm$Z <- las_dtm$Z + 40 - max_z
  # lad_dtm <- LAD(las_dtm@data$Z, z0 = 2)
  # lad_dsm <- LAD(las_dsm@data$Z, z0 = 2)
  # plot(lad_dsm$lad, lad_dsm$z, type = 'l')
  # points(lad_dtm$lad, lad_dtm$z, type = 'l', col = 'blue')
  
  # lad_dsm <- cloud_metrics(las_dsm, ~myPAD_dsm(Classification, Z))
  # lad_dsm <- cloud_metrics(las_dsm, ~myPAD_dtm(Z))
  
  # las_dtm[las_dtm@data$Z < 2]@data$Classification <- 2
  # las_dsm_class2 <- classify_poi(las_dsm, class = int(2), poi = ~las_dsm@data$Z < 2)
  # plot(las_dsm_class2, color = "Classification")
  
  # writeRaster(lad_dtm, "lad/dtm.tif")
  # writeRaster(lad_dtm_class2, "lad/dtm_class2.tif")
  # writeRaster(lad_dsm, "lad/dsm.tif")
  # writeRaster(lad_dsm_class2, "lad/dsm_class2.tif")
  
  # lad_dtm_val <- unlist(lad_dtm)
  # lad_dtm_class2_val <- unlist(lad_dtm_class2)
  # lad_dsm_val <- unlist(lad_dsm)
  # lad_dsm_class2_val <- unlist(lad_dsm_class2)
  # 
  # df <- data.frame(dtm = lad_dtm_val,
  #                  dtm_class2 = lad_dtm_class2_val,
  #                  dsm = lad_dsm_val,
  #                  dsm_class2 = lad_dsm_class2_val
  # )
  # plot(df$dtm, type = "l", 
  #      col = "darkred", 
  #      ylim = c(0, 2)
  #      )
  # points(df$dtm_class2, type = "l", col = "red")
  # points(rev(df$dsm), type = "l", col = "blue")
  # points(rev(df$dsm_class2), type = "l", col = "lightblue")
  # points(df$dsm, type = "l", col = "blue")
  # points(df$dsm_class2, type = "l", col = "lightblue")
  
  # lad_dsm_class2 <- cloud_metrics(las_dsm_class2, ~myPAD_dsm(Classification, Z), res = 10)
  
  # las_filtered <- las_dtm[las_dtm@data$Z < 2]
  # las_filtered@data$Classification <- 2
  
  
  # stop()
  # Read
  # tile <- 60
  # las <- readLAS(las_utm_files[[tile]])
  # st_crs(las) <- st_crs(s2_rast)
  # las <- clip_circle(las, xcenter = median(las@data$X), ycenter = median(las@data$Y),
  # radius = 50)
  
  # thr <- 1
  # DTM
  # las_dtm <- normalize_height(las, algo = tin())
  # las_dtm <- classify_poi(las_dtm, class = int(2), poi = ~las_dtm@data$Z < thr)
  # las_dtm$Z[las_dtm$Classification == 2] <- 0
  
  # las_dtm2 <- end_dtm_norm(las)
  
  # DSM
  # dsmm <- rasterize_canopy(las,
  #                          res = 1,
  #                          pkg = "terra",
  #                          algorithm = p2r()
  # )
  # las_dsm <- normalize_height(las, algo = dsmm)
  # las_dsm$Z <- las_dsm$Z + 40
  # las_dsm <- classify_poi(las_dsm, class = int(2), poi = ~las_dsm@data$Z < thr)
  # las_dsm$Z[las_dsm$Classification == 2] <- 0
  
  # Calculation
  # las_dtm$Z <- las_dtm$Z + 40 - max(las_dtm$Z)
  # 
  # lad_dtm <- LAD(las_dtm@data$Z, z0 = thr)
  # lad_dsm <- LAD(las_dsm@data$Z, z0 = thr)
  # lad_dsm32 <- LAD(las_dsm3@data$Z)
  # lad_dtm2 <- LAD(las_dtm2@data$Z)
  
  # lad_dtm2 <- pixel_metrics(las_dtm, ~myPAD_dtm(Z))
  # lad_dsm2 <- pixel_metrics(las_dsm, ~myPAD_dsm(Classification, Z))
  
  # lad_dtm2 <- pixel_metrics(las_dtm, ~myPAD_dtm_new(Z))
  # lad_dsm2 <- pixel_metrics(las_dsm, ~myPAD_dsm_new(Z, 40), res = 10)
  
  # plot(lad_dsm$lad, lad_dsm$z, type = 'l', main = paste("Threshold =", thr,
  #                                                       "Site =", site,
  #                                                       "Tile =", tile))
  # points(lad_dtm$lad, lad_dtm$z, type = 'l', col = 'blue')
  # points(lad_dtm2$lad, lad_dtm2$z, type = 'l', col = 'green')
  # points(lad_dsm3$lad, lad_dsm3$z, type = 'l', col = 'red')
  
  # las_dsm3 <- catalog_map(ctg_test, end_dsm_norm)
  # las_dsm33 <- end_dsm_norm(las)
  # lad_dsm3 <- pixel_metrics(las_dsm3, ~myPAD_dsm_new(Z, 40), res = 10)
  
  # writeRaster(lad_dsm2, "test_dsm2.tif", overwrite=TRUE)
  # writeRaster(lad_dsm3, "test_dsm.tif", overwrite=TRUE)
  
  # stop()
  # las_dtm <- normalize_height(las, algo = tin())
  # plot(las_dtm)
  # print(summary(las_dtm@data$Z))
  # las_dtm <- NULL
  
  
  # las_dsm <- normalize_height(las, algo = dsmm)
  # print(summary(las_dsm@data$Z))
  # plot(las_dsm, color = "Classification")
  
  # las_inv <- las_dsm
  # las_inv@data$Z <- -las_inv@data$Z
  # print(summary(las_inv@data$Z))
  # plot(las_inv)
  # las_inv <- NULL
  
  # las_dsm_upped <- las_dsm
  # las_dsm_upped@data$Z <- las_dsm_upped@data$Z + 40
  # print(summary(las_dsm_upped@data$Z))
  # plot(las_dsm_upped)
  # las_dsm_upped <- NULL
  # las_dsm <- NULL
  
  # stop()
  
  # Normalization DTM
  # opt_output_files(ctg) <- paste0(file.path(las_norm_dir, "{*}_norm"))
  # normalized_ctg_dtm <- normalize_height(ctg, dtm)
  # normalized_ctg <- catalog_map(ctg, end_dtm_norm)
  # normalized_ctg_dtm <- classify_poi(normalized_ctg_dtm,
  #                                    class = int(2),
  #                                    poi = ~normalized_ctg_dtm@data$Z < 2)
  # normalized_ctg_dtm@data$Z[normalized_ctg_dtm@data$Classification == 2] <- 0
  
  # Normalization DSM
  
  # start_time <- Sys.time()
  
  # opt_output_files(ctg) <- paste0(file.path(las_norm_dsm_dir, "{*}_norm"))
  opt_output_files(ctg) <- file.path(las_norm_dsm_dir, "{XLEFT}_{YBOTTOM}_norm")
  # normalized_ctg <- normalize_height(ctg, algo = dsm)
  output_dsm <- catalog_map(ctg, end_dsm_norm)
  # las_utm_norm_dsm_files <- list.files(las_norm_dsm_dir, pattern = "\\.las$", full.names = TRUE)
  # chm_values <- values(chm, na.rm = TRUE)
  # chm_values[chm_values > 40] <- 40
  # max_chm <- max(values(chm), na.rm = TRUE)
  
  # total_values_count <- length(chm_values)
  # greater_than_40_count <- sum(chm_values > 40)
  # greater_than_45_count <- sum(chm_values > 45)
  # greater_than_50_count <- sum(chm_values > 50)
  #
  # # Print the counts
  # print(paste("Total number of values:", total_values_count))
  # print(paste("Number of values greater than 40:", greater_than_40_count))
  # print(paste("Number of values greater than 45:", greater_than_45_count))
  # print(paste("Number of values greater than 50:", greater_than_50_count))
  
  
  # stop()
  
  # max_chm <- 40
  # # for (i in seq_along(las_utm_norm_dsm_files)){
  # for (i in seq_along(las_utm_files)){
  #   print(i)
  #   # las_i <- readLAS(las_utm_norm_dsm_files[i])
  #   las_i <- readLAS(las_utm_files[i])
  #   dtm_i <- rast(dtm_files[i])
  #   # print(las_i$Z)
  #   las_i <- las_i - dtm_i
  #   las_i <- classify_poi(las_i, class = int(2), poi = ~las_i@data$Z < 2)
  #   las_i$Z[las_i$Classification == 2] <- 0
  #   las_i$Z <- las_i$Z + max_chm - max(las_i$Z)
  #   # print(las_i$Z)
  #   writeLAS(las_i, las_utm_norm_dsm_files[i])
  # }
  
  
  
  
  
  # max_chm <- 40
  # for (i in seq_along(las_utm_norm_dsm_files)){
  # # for (i in seq_along(las_utm_files)){
  #   print(i)
  #   las_i <- readLAS(las_utm_norm_dsm_files[i])
  #   # las_i <- readLAS(las_utm_files[i])
  #   # dsm_i <- rast(dsm_files[i])
  #   # print(las_i$Z)
  #   # las_i <- las_i - dsm_i
  #   las_i$Z <- las_i$Z + max_chm
  #   las_i <- classify_poi(las_i, class = int(2), poi = ~las_i@data$Z < 2)
  #   las_i$Z[las_i$Classification == 2] <- 0
  #   # print(las_i$Z)
  #   writeLAS(las_i, las_utm_norm_dsm_files[i])
  # }
  # end_time <- Sys.time()
  # cat(paste(site, "normalization time:", format_time(start_time, end_time)),
  #     "\n")
}
