# ---
# title: "3.calculate_lidar_metrics.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-07-22"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

library("lidR")
library("lidRmetrics")
library("raster")
library("plotly")
library("terra")
library("viridis")
library("future")
library("rgl")
library("lmom")
library("data.table")
library("solaR")

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

# Import useful functions
source("../libraries/functions_lidar.R")
source("../libraries/functions_create_masks.R")
source("../libraries/functions_plots.R")

# Directories
data <- "../../01_DATA"
# sites <- "Mormal" # Mormal Blois Aigoual
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- c("Aigoual", "Blois")
# sites <- c("Mormal")
data_dir <- "/media/corroyez/My Passport/01_DATA"
for (site in sites){
  s2_rast <- switch(site,
                    "Aigoual" = terra::rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Aigoual/L2A_T31TEJ_A031608_20210711T104217/Reflectance/res_10_m/L2A_T31TEJ_A031608_20210711T104217_Refl")[[1]],
                    "Blois" = terra::rast('/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/L2A_T31TCN_A031222_20210614T105443/Reflectance/res_10_m/L2A_T31TCN_A031222_20210614T105443_Refl')[[1]],
                    "Mormal" = terra::rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/Reflectance/res_10_m/L2A_T31UER_A031222_20210614T105443_Refl")[[1]]
  )
  
  data_site <- file.path(data, site)
  results_path <- file.path("../../03_RESULTS", site)
  masks_dir <- file.path(results_path, "LiDAR/Heterogeneity_Masks")
  # LAS directories
  # UTM: for all metrics except LAI
  # las_utm <- "LiDAR/2-las_utm/"
  # las_utm_dir <- file.path(data_site, las_utm)
  las_utm_dir <- file.path(data_dir,
                           site,
                           "LiDAR",
                           "2-las_utm"
  )
  las_utm_files <- list.files(las_utm_dir, pattern = "\\.las$", full.names = TRUE)
  
  # Normalized UTM: for LAI only
  # las_norm_utm <- "LiDAR/3-las_normalized_utm/"
  # las_norm_utm_dir <- file.path(data_site, las_norm_utm)
  las_norm_utm_dir <- file.path(data_dir,
                                site,
                                "LiDAR",
                                "3-las_normalized_utm"
  )
  las_norm_utm_files <- list.files(las_norm_utm_dir, pattern = "\\.las$", full.names = TRUE)
  
  # Normalized UTM DSM: for LAI only
  # las_norm_utm_dsm <- "LiDAR/4-las_normalized_utm_dsm/"
  # las_norm_utm_dsm_dir <- file.path(data_site, las_norm_utm_dsm)
  las_norm_utm_dsm_dir <- file.path(data_dir,
                                    site,
                                    "LiDAR",
                                    "4-las_normalized_utm_dsm" # 4-las_normalized_utm_dsm
  )
  las_norm_utm_dsm_files <- list.files(las_norm_utm_dsm_dir, pattern = "\\.las$", full.names = TRUE)
  
  # Output
  results_path <- file.path("../../03_RESULTS", site)
  metrics_dir <- file.path(results_path, "Metrics")
  raw_dir <- file.path(metrics_dir, "Raw")
  if (!dir.exists(raw_dir)) {
    dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Catalogs
  # normalized_ctg <- readLAScatalog(las_norm_utm_dir)
  # ctg <- readLAScatalog(las_utm_dir)
  # normalized_dsm_ctg <- readLAScatalog(las_norm_utm_dsm_files)
  # plot(ctg)
  
  # Resolution
  res <- 10
  
  # LiDR optimization
  # opt_chunk_size(normalized_ctg) <- 0 # Processing by files
  # opt_chunk_buffer(normalized_ctg) <- 10
  # opt_select(normalized_ctg) <- "xyz"
  # opt_chunk_size(ctg) <- 0 # Processing by files
  # opt_chunk_buffer(ctg) <- 10
  # opt_select(ctg) <- "xyz"
  # opt_chunk_size(normalized_dsm_ctg) <- 0 # Processing by files
  # opt_chunk_buffer(normalized_dsm_ctg) <- 10
  # opt_select(normalized_dsm_ctg) <- "xyzc"
  # set_lidr_threads(10)
  
  # CHM
  chm <- terra::rast(file.path(results_path, "LiDAR/chm/res_1_m/chm.tif"))
  
  # DTM
  dtm <- terra::rast(file.path(results_path, "LiDAR/dtm/res_10_m/rasterize_terrain.vrt"))
  
  # DSM
  dsm <- terra::rast(file.path(results_path, "LiDAR/dsm/res_1_m/rasterize_canopy.vrt"))
  
  # Load a raster to project metrics at its resolution and coordinates
  mask_10m <- terra::rast(file.path(results_path,
                                    "LiDAR/Heterogeneity_Masks/mnc_mean_heights_res_10_m.envi"))
  
  # dir <- file.path(metrics_dir, "Not_Masked/PAD_Profiles_updated")
  # files <- list.files(dir, full.names = T, pattern = "\\.tif$")
  # raster_stack <- rast(files)
  # cv_raster <- (app(raster_stack, sd) / app(raster_stack, mean)) * 100
  # cv_lad_filename <- "cv_lad_res_10_m.tif"
  # writeRaster(cv_raster,
  #             file.path(raw_dir, cv_lad_filename),
  #             overwrite = TRUE)
  # apply_and_save_masks(cv_raster,
  #                      cv_lad_filename,
  #                      masks_dir,
  #                      metrics_dir)
  # Canopy cover
  # chm_1m_masked <- mask_chm_thresh(chm, threshold = 2, masks_dir, res)
  # canopy_cover <- terra::project(chm_1m_masked, mask_10m, method = 'average')
  # canopy_cover_file <- "fCover_res_10_m.tif"
  # writeRaster(canopy_cover,
  #             filename = file.path(raw_dir, canopy_cover_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(canopy_cover,
  #                      canopy_cover_file,
  #                      masks_dir,
  #                      metrics_dir)
  
  # Gap Fraction
  # gap_fraction_raster <- lidR::pixel_metrics(normalized_ctg,
  #                                            ~wrapper_gap_fraction(Z),
  #                                            res = res)
  # gap_fraction_raster <- 1 - canopy_cover / terra::project(chm_1m_masked, mask_10m, method = 'max')
  # gap_fraction_file <- "gap_fraction_res_10_m.tif"
  # writeRaster(gap_fraction_raster,
  #             filename = file.path(raw_dir, gap_fraction_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(gap_fraction_raster,
  #                      gap_fraction_file,
  #                      masks_dir,
  #                      metrics_dir)
  
  # opt_filter(ctg) <- "-drop_class 2 -drop_z_below 2"
  # st_crs(ctg) <- st_crs(s2_rast)
  # st_crs(normalized_ctg) <- st_crs(s2_rast)
  # st_crs(normalized_dsm_ctg) <- st_crs(s2_rast)
  
  # Point clouds-related metrics
  # Lmoms (Lcv and Lskew)
  # lmoms_rasters <- lidR::pixel_metrics(ctg,
  #                                      ~as.list(lmom::samlmu(Z, nmom=3, ratios=F)),
  #                                      res = res)
  # 
  # lmom_1 <- lmoms_rasters[[1]]
  # lmom_2 <- lmoms_rasters[[2]]
  # lmom_3 <- lmoms_rasters[[3]]
  # lmom_1 <- terra::project(lmom_1, mask_10m)
  # lmom_2 <- terra::project(lmom_2, mask_10m)
  # lmom_3 <- terra::project(lmom_3, mask_10m)
  
  # Lcv
  # lcv <- lmom_2 / lmom_1
  # lcv_file <- "lcv_no_ground.tif"
  # writeRaster(lcv,
  #             filename = file.path(raw_dir, lcv_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(lcv,
  #                      lcv_file,
  #                      masks_dir,
  #                      metrics_dir)
  # Lskew
  # lskew <- lmom_3 / lmom_2
  # lskew_file <- "lskew_no_ground.tif"
  # writeRaster(lskew,
  #             filename = file.path(raw_dir, lskew_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(lskew,
  #                      lskew_file,
  #                      masks_dir,
  #                      metrics_dir)
  
  # VCI
  # vci <- lidR::pixel_metrics(ctg, ~VCI_local(Z), res = res)
  # vci <- terra::project(vci, mask_10m)
  # vci_file <- "vci_no_ground.tif"
  # writeRaster(vci,
  #             filename = file.path(raw_dir, vci_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(vci,
  #                      vci_file,
  #                      masks_dir,
  #                      metrics_dir)
  
  # VDR
  # vdr <- lidR::pixel_metrics(ctg, ~VDR(Z), res = res)
  # vdr <- terra::project(vdr, mask_10m)
  # vdr_file <- "vdr_res_10_m.tif"
  # writeRaster(vdr,
  #             filename = file.path(raw_dir, vdr_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(vdr,
  #                      vdr_file,
  #                      masks_dir,
  #                      metrics_dir)
  
  # CV (of norm ctg)
  # cvnorm <- lidR::pixel_metrics(normalized_ctg, ~CV(Z), res = res)
  # cvnorm <- terra::project(cvnorm, mask_10m)
  # cvnorm_file <- "cvnorm_res_10_m.tif"
  # writeRaster(cvnorm,
  #             filename = file.path(raw_dir, cvnorm_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(cvnorm,
  #                      cvnorm_file,
  #                      masks_dir,
  #                      metrics_dir)
  #
  
  # PAI calculation
  # pai <- lidR::pixel_metrics(normalized_ctg, ~myPAI(Z, zmin=2), res=s2_rast)
  # pai[pai == Inf] <- NA
  # pai <- terra::project(pai, mask_10m)
  # pai_file <- "lidarlai_res_10_m.tif"
  # writeRaster(pai,
  #             filename = file.path(raw_dir, pai_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(pai,
  #                      pai_file,
  #                      masks_dir,
  #                      metrics_dir)
  
  # pad_rasters <- lidR::pixel_metrics(normalized_ctg,
  #                                    ~myPAD(Z
  #                                    ),
  #                                    res=res)
  # pad_rasters <- terra::project(pad_rasters, mask_10m)
  # 
  # if (!dir.exists(file.path(raw_dir, "PAD_Profiles"))) {
  #   dir.create(file.path(raw_dir, "PAD_Profiles"),
  #              showWarnings = FALSE,
  #              recursive = TRUE)
  # }
  # for (i in 1:nlyr(pad_rasters)) {
  #   layer <- subset(pad_rasters, i)
  #   filename <- paste0("PAD_",
  #                      names(pad_rasters)[i],
  #                      "_",
  #                      40,
  #                      ".tif")
  #   pad_bool <- TRUE
  #   if (i == nlyr(pad_rasters)){
  #     filename <- "cvlad_res_10_m.tif"
  #     pad_bool <- FALSE
  #   }
  #   writeRaster(layer,
  #               filename = file.path(raw_dir, filename),
  #               overwrite = TRUE)
  #   apply_and_save_masks(layer,
  #                        filename,
  #                        masks_dir,
  #                        metrics_dir,
  #                        pad_bool = pad_bool,
  #                        norm_bool = FALSE)
  # }
  
  # PAD Old Norm calculation
  # pad_rasters <- lidR::pixel_metrics(normalized_ctg,
  #                                    ~myPAD(Z
  #                                    ),
  #                                    res=s2_rast)
  # # pad_rasters <- terra::project(pad_rasters, mask_10m)
  # if (!dir.exists(file.path(raw_dir, "PAD_Profiles_NA"))) {
  #   dir.create(file.path(raw_dir, "PAD_Profiles_NA"),
  #              showWarnings = FALSE,
  #              recursive = TRUE)
  # }
  # # Get the names of all layers
  # layer_names <- names(pad_rasters)
  # 
  # # Extract, stack and save the first 38 LAD layers
  # lad_layers <- layer_names[grepl("LAD_Layer_", layer_names)]
  # lad_extracted_values <- as.numeric(sub("LAD_Layer_", "", lad_layers))
  # lad_sorted_indices <- order(lad_extracted_values, decreasing = TRUE)
  # lad_selected_layers <- lad_layers[lad_sorted_indices]
  # lad_raster <- subset(pad_rasters, lad_selected_layers)
  # lad_filename <- "ladstack.tif"
  # writeRaster(lad_raster,
  #             file.path(raw_dir, lad_filename),
  #             overwrite = TRUE)
  # apply_and_save_masks(lad_raster,
  #                      lad_filename,
  #                      masks_dir,
  #                      metrics_dir,
  #                      pad_bool = TRUE)
  # 
  # # Extract and save the next 38 PAI layers individually
  # pai_layers <- layer_names[grepl("PAI_", layer_names)]
  # pai_extracted_values <- as.numeric(sub("PAI_", "", pai_layers))
  # pai_sorted_indices <- order(pai_extracted_values)
  # pai_selected_layers <- pai_layers[pai_sorted_indices[1:38]]
  # 
  # for (i in seq_along(pai_selected_layers)) {
  #   layer <- subset(pad_rasters, pai_selected_layers[i])
  #   layer_name <- names(layer)
  #   filename <- paste0("PAD_", sub("PAI_", "", layer_name), "_40.tif")
  #   filepath <- file.path(raw_dir, filename)
  #   writeRaster(layer, filepath, overwrite = TRUE)
  #   apply_and_save_masks(layer,
  #                        filename,
  #                        masks_dir,
  #                        metrics_dir,
  #                        pad_bool = TRUE)
  # }
  # # a
  # #
  # # # Extract and save the last 38 CV_LAD layers individually
  # # cv_lad_layers <- layer_names[grepl("CV_LAD_", layer_names)]
  # # cv_lad_extracted_values <- as.numeric(sub("CV_LAD_", "", cv_lad_layers))
  # # cv_lad_sorted_indices <- order(cv_lad_extracted_values)
  # # cv_lad_selected_layers <- cv_lad_layers[cv_lad_sorted_indices[1:38]]
  # #
  # # if (!dir.exists(file.path(raw_dir, "CV_LAD"))) {
  # #   dir.create(file.path(raw_dir, "CV_LAD"), showWarnings = FALSE, recursive = TRUE)
  # # }
  # #
  # # for (i in seq_along(cv_lad_selected_layers)) {
  # #   layer <- subset(pad_rasters, cv_lad_selected_layers[i])
  # #   layer_name <- names(layer)
  # #   filename <- paste0("CV_LAD_", sub("CV_LAD_", "", layer_name), "_40.tif")
  # #   writeRaster(layer, file.path(raw_dir, filename), overwrite = TRUE)
  # #   apply_and_save_masks(layer,
  # #                        filename,
  # #                        masks_dir,
  # #                        metrics_dir,
  # #                        pad_bool = TRUE,
  # #                        norm_bool = FALSE,
  # #                        cvlad_bool = TRUE)
  # # }
  # 
  # # PAD New Norm calculation
  # pad_rasters <- lidR::pixel_metrics(normalized_dsm_ctg,
  #                                    ~myPAD_updated(Classification,
  #                                                   Z
  #                                    ),
  #                                    res=s2_rast)
  # pad_rasters <- terra::project(pad_rasters, mask_10m)
  # cv_lad <- pad_rasters$cv_lad
  # cv_lad_filename <- "cv_lad_res_10_m.tif"
  # writeRaster(cv_lad,
  #             file.path(raw_dir, cv_lad_filename),
  #             overwrite = TRUE)
  # apply_and_save_masks(cv_lad,
  #                      cv_lad_filename,
  #                      masks_dir,
  #                      metrics_dir)
  # 
  # # pad_rasters_old <- lidR::pixel_metrics(normalized_ctg,
  # #                                        ~myPAD(Z
  # #                                        ),
  # #                                        res=res)
  # # pad_rasters_old <- terra::project(pad_rasters_old, mask_10m)
  # # cv_lad <- pad_rasters_old$cv_lad
  # # cv_lad_filename <- "cv_lad_dtm_res_10_m.tif"
  # # writeRaster(cv_lad,
  # #             file.path(raw_dir, cv_lad_filename),
  # #             overwrite = TRUE)
  # # apply_and_save_masks(cv_lad,
  # #                      cv_lad_filename,
  # #                      masks_dir,
  # #                      metrics_dir)
  # 
  # if (!dir.exists(file.path(raw_dir, "PAD_Profiles_updated_modifminz_NA"))) {
  #   dir.create(file.path(raw_dir, "PAD_Profiles_updated_modifminz_NA"),
  #              showWarnings = FALSE,
  #              recursive = TRUE)
  # }
  # # # Get the names of all layers
  # layer_names <- names(pad_rasters)
  # #
  # # # Extract, stack and save the first 38 LAD layers
  # lad_layers <- layer_names[grepl("LAD_Layer_", layer_names)]
  # lad_extracted_values <- as.numeric(sub("LAD_Layer_", "", lad_layers))
  # lad_sorted_indices <- order(lad_extracted_values, decreasing = TRUE)
  # lad_selected_layers <- lad_layers[lad_sorted_indices]
  # lad_raster <- subset(pad_rasters, lad_selected_layers)
  # lad_filename <- "ladstack.tif"
  # writeRaster(lad_raster,
  #             file.path(raw_dir, lad_filename),
  #             overwrite = TRUE)
  # apply_and_save_masks(lad_raster,
  #                      lad_filename,
  #                      masks_dir,
  #                      metrics_dir,
  #                      pad_bool = TRUE,
  #                      norm_bool = TRUE)
  # #
  # # Extract and save the next 38 PAI layers individually
  # pai_layers <- layer_names[grepl("PAI_", layer_names)]
  # pai_extracted_values <- as.numeric(sub("PAI_", "", pai_layers))
  # pai_sorted_indices <- order(pai_extracted_values)
  # pai_selected_layers <- pai_layers[pai_sorted_indices[1:38]]
  # 
  # for (i in seq_along(pai_selected_layers)) {
  #   layer <- subset(pad_rasters, pai_selected_layers[i])
  #   layer_name <- names(layer)
  #   filename <- paste0("PAD_", sub("PAI_", "", layer_name), "_40.tif")
  #   filepath <- file.path(raw_dir, filename)
  #   writeRaster(layer, filepath, overwrite = TRUE)
  #   apply_and_save_masks(layer,
  #                        filename,
  #                        masks_dir,
  #                        metrics_dir,
  #                        pad_bool = TRUE,
  #                        norm_bool = TRUE)
  # }
  # a
  #
  # # Extract and save the last 38 CV_LAD layers individually
  # cv_lad_layers <- layer_names[grepl("CV_LAD_", layer_names)]
  # cv_lad_extracted_values <- as.numeric(sub("CV_LAD_", "", cv_lad_layers))
  # cv_lad_sorted_indices <- order(cv_lad_extracted_values)
  # cv_lad_selected_layers <- cv_lad_layers[cv_lad_sorted_indices[1:38]]
  #
  # if (!dir.exists(file.path(raw_dir, "CV_LAD"))) {
  #   dir.create(file.path(raw_dir, "CV_LAD"), showWarnings = FALSE, recursive = TRUE)
  # }
  #
  # for (i in seq_along(cv_lad_selected_layers)) {
  #   layer <- subset(pad_rasters, cv_lad_selected_layers[i])
  #   layer_name <- names(layer)
  #   filename <- paste0("CV_LAD_", sub("CV_LAD_", "", layer_name), "_40.tif")
  #   writeRaster(layer, file.path(raw_dir, filename), overwrite = TRUE)
  #   apply_and_save_masks(layer,
  #                        filename,
  #                        masks_dir,
  #                        metrics_dir,
  #                        pad_bool = TRUE,
  #                        norm_bool = TRUE,
  #                        cvlad_bool = TRUE)
  # }
  
  # Shadows
  # shade <- perform_shadows_analysis(dsm, metrics_dir)
  # shade_file <- "shade_res_10_m.tif"
  # writeRaster(shade,
  #             filename = file.path(raw_dir, shade_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(shade,
  #                      shade_file,
  #                      masks_dir,
  #                      metrics_dir)
  
  # HillShade
  # hillshade <- my_hillshade(dsm) # 1 -
  # hillshade_file <- "hillshade_res_10_m.tif"
  # writeRaster(hillshade,
  #             filename = file.path(raw_dir, hillshade_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(hillshade,
  #                      hillshade_file,
  #                      masks_dir,
  #                      metrics_dir)
  
  # Calculate slope using the DTM
  # slope <- terrain(dtm, v = "slope", unit = "degrees")
  # slope_file <- "slope_res_10_m.tif"
  # writeRaster(slope,
  #             filename = file.path(raw_dir, slope_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(slope,
  #                      slope_file,
  #                      masks_dir,
  #                      metrics_dir)
  
  # Calculate aspect (orientation) using the DTM
  # aspect <- terrain(dtm, v = "aspect", unit = "degrees")
  # aspect_file <- "aspect_res_10_m.tif"
  
  # Function to classify aspect into 8 directional classes
  # classify_aspect <- function(aspect) {
  #   ifelse((aspect >= 337.5 | aspect < 22.5), 1,  # North
  #          ifelse((aspect >= 22.5 & aspect < 67.5), 2,  # North-East
  #                 ifelse((aspect >= 67.5 & aspect < 112.5), 3,  # East
  #                        ifelse((aspect >= 112.5 & aspect < 157.5), 4,  # South-East
  #                               ifelse((aspect >= 157.5 & aspect < 202.5), 5,  # South
  #                                      ifelse((aspect >= 202.5 & aspect < 247.5), 6,  # South-West
  #                                             ifelse((aspect >= 247.5 & aspect < 292.5), 7,  # West
  #                                                    8)))))))  # North-West
  # }
  
  # Apply classification to the aspect raster
  # aspect_classes <- app(aspect, fun = classify_aspect)
  
  # Label the aspect classes (for later reference)
  # aspect_labels <- c("North", "North-East", "East", "South-East",
  #                    "South", "South-West", "West", "North-West")
  
  # # Plot the classified aspect raster
  # plot(aspect_classes, main = "Aspect Classified into 8 Directions", col = rainbow(8), legend = FALSE)
  # legend("topright", legend = aspect_labels, fill = rainbow(8))
  
  # writeRaster(aspect_classes,
  #             filename = file.path(raw_dir, aspect_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(aspect_classes,
  #                      aspect_file,
  #                      masks_dir,
  #                      metrics_dir)
  
  # Northness
  # northness <- cos(aspect * pi / 180)
  # northness_file <- "northness_res_10_m.tif"
  # writeRaster(northness,
  #             filename = file.path(raw_dir, northness_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(northness,
  #                      northness_file,
  #                      masks_dir,
  #                      metrics_dir)
  
  # Northness slope
  # northness_slope <- slope * northness
  # northness_slope_file <- "northness_slope_res_10_m.tif"
  # writeRaster(northness,
  #             filename = file.path(raw_dir, northness_slope_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(northness_slope,
  #                      northness_slope_file,
  #                      masks_dir,
  #                      metrics_dir)
  
  # TWI
  # flow_dir <- terrain(dtm, v = "flowdir")
  # flow_accum <- flowAccumulation(flow_dir)
  # twi <- log(flow_accum / tan(slope * pi / 180))
  # twi_file <- "twi_res_10_m.tif"
  # writeRaster(twi,
  #             filename = file.path(raw_dir, twi_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(twi,
  #                      twi_file,
  #                      masks_dir,
  #                      metrics_dir)
  
  # Convert aspect to radians
  # aspect_radians <- aspect * pi / 180
  
  # Calculate sin and cos of aspect
  # sin_aspect <- sin(aspect_radians)  # East-West component
  # cos_aspect <- cos(aspect_radians)  # North-South component
  #
  # sin_aspect_file <- "aspect_sin_res_10_m.tif"
  # cos_aspect_file <- "aspect_cos_res_10_m.tif"
  # writeRaster(sin_aspect,
  #             filename = file.path(raw_dir, sin_aspect_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(sin_aspect,
  #                      sin_aspect_file,
  #                      masks_dir,
  #                      metrics_dir)
  # writeRaster(cos_aspect,
  #             filename = file.path(raw_dir, cos_aspect_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(cos_aspect,
  #                      cos_aspect_file,
  #                      masks_dir,
  #                      metrics_dir)
  
  # Solar Incidence
  # if (site == "Aigoual"){
  #   t <- as.POSIXct("2021-07-11 10:40:31", tz = "UTC")
  #   lat <- 44.12
  #   lon <- 3.52
  # } else if (site == "Blois"){
  #   t <- as.POSIXct("2021-06-14 10:50:31", tz = "UTC")
  #   lat <- 47.57
  #   lon <- 1.29
  # } else if (site == "Mormal"){
  #   t <- as.POSIXct("2021-06-14 10:50:31", tz = "UTC")
  #   lat <- 50.20
  #   lon <- 3.74
  # } else{
  #   stop("Error: Site must be Mormal, Blois, or Aigoual.\n")
  # }
  
  # Get sun angle and altitude
  # df_sun <- oce::sunAngle(t, longitude = lon, latitude = lat)
  #
  # # Extract altitude and azimuth angles
  # solar_elevation <- df_sun$altitude  # in degrees
  # solar_azimuth <- df_sun$azimuth     # in degrees
  #
  # # Convert angles from degrees to radians
  # solar_elevation_rad <- solar_elevation * (pi / 180)
  # solar_azimuth_rad <- solar_azimuth * (pi / 180)
  # slope_rad <- slope * (pi / 180)
  # aspect_rad <- aspect * (pi / 180)
  #
  # declination <- df_sun$declination
  # declination_rad <- declination * (pi / 180)
  # latitude_rad <- lat * (pi / 180)
  # longitude_sun_rad <- df_sun$azimuth * (pi / 180)
  #
  # # Calculate solar incidence angle
  # incidence_angle <- acos(cos(solar_elevation_rad) * cos(slope_rad) +
  #   sin(solar_elevation_rad) * sin(slope_rad) * cos(solar_azimuth_rad - aspect_rad))
  #
  # zenith_angle <- asin(sin(declination_rad) * sin(latitude_rad) +
  #                        cos(declination_rad) * cos(latitude_rad) * cos(longitude_sun_rad - lon * (pi / 180)))
  
  # angles <- calculate_solar_angles(lat, long, t, slope, azimuth = 180)
  # ratio <- cos(angles$zenith_angle) / cos(angles$incidence_angle)
  #
  # incidence_angle_file <- "solar_incidence_angle_res_10_m.tif"
  # writeRaster(ratio,
  #             filename = file.path(raw_dir, incidence_angle_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(ratio,
  #                      incidence_angle_file,
  #                      masks_dir,
  #                      metrics_dir)
  
  # # CHM-related metrics
  #
  # # Max
  # max <- terra::project(chm, mask_10m, method = 'max')
  # max_file <- "max_res_10_m.tif"
  # writeRaster(max,
  #             filename = file.path(raw_dir, max_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(max,
  #                      max_file,
  #                      masks_dir,
  #                      metrics_dir)
  #
  #
  # Mean
  mean <- terra::project(chm, s2_rast, method = 'average')
  mean_file <- "mean_res_10_m.tif"
  writeRaster(mean,
              filename = file.path(raw_dir, mean_file),
              filetype = "GTiff",
              overwrite = TRUE)
  apply_and_save_masks(mean,
                       mean_file,
                       masks_dir,
                       metrics_dir)
  
  # Rumple index
  # rumple <- catalog_map(ctg,
  #                       rumple_index_surface,
  #                       res = res,
  #                       .options = list(raster_alignment = res))
  
  # Define a window size (e.g., 3x3, 5x5). Adjust based on your needs
  # window_size <- 3
  #
  # # Define a custom function to compute Rumple Index for each window
  # rumple_focal_function <- function(x) {
  #   # Remove NAs from the focal window
  #   x <- na.omit(x)
  #
  #   # Check if there are enough valid values to compute Rumple Index
  #   if (length(x) >= 4) {  # Ensure there are enough points for a meaningful computation
  #     # Create a temporary raster for the valid values
  #     temp_raster <- terra::rast(matrix(x, nrow = window_size, ncol = window_size, byrow = TRUE))
  #     # Compute and return the Rumple Index
  #     return(rumple_index(temp_raster))
  #   } else {
  #     return(NA)  # Return NA if not enough valid values
  #   }
  # }
  #
  # rumple_function <- function(x) {
  #   x <- as.vector(x)  # Convert matrix to vector
  #   x <- x[!is.na(x)]  # Remove NA values
  #
  #   # Ensure there are enough values to calculate the Rumple Index
  #   if (length(x) < 2) return(NA)  # Not enough data
  #
  #   return(lidR::rumple_index(x))  # Calculate Rumple Index
  # }
  #
  # # Apply the focal function to the raster
  # rumple_raster <- terra::focal(mean_noNA,
  #                               w = matrix(1, nrow = window_size, ncol = window_size),
  #                               fun = rumple_function,
  #                               na.policy = "omit")  # omit NAs in calculations
  #
  # # Plot the resulting Rumple Index raster
  # plot(rumple_raster)
  
  # Initialize list to store rumple index results for each file
  # rumple_indices <- list()
  #
  # # Iterate over all las_utm_files using the correct loop sequence
  # for (i in seq_along(las_utm_files)) {
  #   print(paste("Processing file:", i))
  #
  #   # Read the LAS file
  #   las_i <- readLAS(las_utm_files[i])
  #
  #   # Step 1: Rasterize the DTM
  #   dtm <- rasterize_terrain(las_i, res = 1, pkg = "terra",
  #                            algorithm = tin(extrapolate = kriging()))
  #
  #   # Step 2: Rasterize the DSM
  #   dsm <- rasterize_canopy(las_i, res = 1, pkg = "terra",
  #                           algorithm = pitfree(thresholds = c(0, 10, 20),
  #                                               max_edge = c(0, 1)))
  #
  #   # Step 3: Calculate the CHM (Canopy Height Model)
  #   chm <- dsm - dtm
  #
  #   # Step 4: Calculate the Rumple Index from the CHM
  #   rumple_index_raster_i <- rumple_index(chm)
  #
  #   # Store the result in the list
  #   rumple_indices[[i]] <- rumple_index_raster_i
  # }
  # rumple <- terra::project(rumple, mask_10m)
  # rumple_file <- "rumple_res_10_m.tif"
  # writeRaster(rumple,
  #             filename = file.path(raw_dir, rumple_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(rumple,
  #                      rumple_file,
  #                      masks_dir,
  #                      metrics_dir)
  
  # # Standard deviation (std)
  # std <- aggregate(chm, fact=res, fun="std")
  # std <- terra::project(std,
  #                       mask_10m,
  #                       method = 'bilinear')
  # std_file <- "std_res_10_m.tif"
  # writeRaster(std,
  #             filename = file.path(raw_dir, std_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(std,
  #                      std_file,
  #                      masks_dir,
  #                      metrics_dir)
  # #
  # # # Coefficient of variation (cv)
  # cv <- 100 * std / mean
  # cv_file <- "cv_res_10_m.tif"
  # writeRaster(cv,
  #             filename = file.path(raw_dir, cv_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(cv,
  #                      cv_file,
  #                      masks_dir,
  #                      metrics_dir)
  # #
  # # # Variance
  # variance <- std * std
  # variance_file <- "variance_res_10_m.tif"
  # writeRaster(variance,
  #             filename = file.path(raw_dir, variance_file),
  #             filetype = "GTiff",
  #             overwrite = TRUE)
  # apply_and_save_masks(variance,
  #                      variance_file,
  #                      masks_dir,
  #                      metrics_dir)
}
