# ---
# title: "functions_masks.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-05-27"
# ---

library("lidR")
library("raster")
library("plotly")
library("terra")
library("viridis")
library("future")
library("sf")
library("dplyr")

#' Mask Forest Composition Area
#'
#' This function masks the forest composition within the specified site edges,
#' creating separate masks for deciduous-only, deciduous-flex,
#'
#' @param site_edges A character string specifying the path to the shapefile
#' containing the edges of the site.
#' @param forest_composition A character string specifying the path to the
#' shapefile containing forest composition information.
#' @param masks_dir A character string specifying the directory where the
#' generated masks will be saved.
#' @param resolution The resolution of the raster masks to be created.
#' @return A list containing raster masks for full composition, deciduous-only,
#' deciduous-flex, coniferous-only and coniferous-flex areas.
#' @export
mask_forest_composition_area <- function(site_edges,
                                         forest_composition,
                                         masks_dir,
                                         resolution){
  # Read shapefiles
  site_edges_shp <- st_read(site_edges)
  forest_composition_shp <- st_read(forest_composition)

  # Convert forest_composition from Lambert-93 to UTM 31N
  if (st_crs(site_edges_shp) != st_crs(forest_composition_shp)) {
    forest_composition_shp <- st_transform(forest_composition_shp,
                                           st_crs(site_edges_shp))
  }

  # Perform intersection: forest composition on the site only
  merged_shp <- st_intersection(site_edges_shp, forest_composition_shp)
  st_write(merged_shp, file.path(masks_dir,"site_full_composition.shp"),
           append=FALSE)

  # Define deciduous and coniferous terms
  deciduous_terms <- c("chêne", "chênes",
                       "hêtre", "hêtres",
                       "feuillu", "feuillus", "peupleraie")
  coniferous_terms <- c("conifère", "conifères",
                        "douglas", "pin", "pins",
                        "sapin", "épicéa")
  oak_term <- c("chêne", "chênes")
  beech_term <- c("hêtre", "hêtres")

  # Create a logical vector for inclusion and exclusion
  deciduous_logical <- grepl(paste(deciduous_terms, collapse = "|"),
                             merged_shp$TFV,
                             ignore.case = TRUE)
  coniferous_logical <- grepl(paste(coniferous_terms, collapse = "|"),
                              merged_shp$TFV,
                              ignore.case = TRUE)
  # oak_logical <- grepl(paste(oak_term, collapse = "|"),
  #                      merged_shp$TFV,
  #                      ignore.case = TRUE)
  # beech_logical <- grepl(paste(beech_term, collapse = "|"),
  #                        merged_shp$TFV,
  #                        ignore.case = TRUE)

  # Filter to create composition masks
  deciduous_only_shp <- merged_shp[deciduous_logical & !coniferous_logical, ]
  deciduous_flex_shp <- merged_shp[deciduous_logical, ]
  coniferous_only_shp <- merged_shp[!deciduous_logical & coniferous_logical, ]
  coniferous_flex_shp <- merged_shp[coniferous_logical, ]
  # oak_only_shp <- merged_shp[oak_logical, ]
  # beech_only_shp <- merged_shp[beech_logical, ]

  st_write(deciduous_only_shp, file.path(masks_dir,"site_deciduous_only.shp"),
           append=F, overwrite = T)
  st_write(deciduous_flex_shp, file.path(masks_dir,"site_deciduous_flex.shp"),
           append=F, overwrite = T)
  st_write(coniferous_only_shp, file.path(masks_dir,"site_coniferous_only.shp"),
           append=F, overwrite = T)
  st_write(coniferous_flex_shp, file.path(masks_dir,"site_coniferous_flex.shp"),
           append=F, overwrite = T)
  # st_write(oak_only_shp, file.path(masks_dir,"site_oak_only.shp"),
  #          append=FALSE)
  # st_write(beech_only_shp, file.path(masks_dir,"site_beech_only.shp"),
  #          append=FALSE)

  # Rasterize full composition, deciduous-flex and deciduous-only mask
  raster_mask <- raster(extent(deciduous_only_shp), resolution = resolution)
  raster_mask[] <- 1  # Set all cells to a default value, e.g., 1

  full_composition_mask <- rasterize(merged_shp, raster_mask)
  full_composition_mask <- terra::rast(full_composition_mask)
  full_composition_mask[!is.na(full_composition_mask)] <- 1

  deciduous_only_mask <- rasterize(deciduous_only_shp, raster_mask)
  deciduous_only_mask <- terra::rast(deciduous_only_mask)
  deciduous_only_mask[!is.na(deciduous_only_mask)] <- 1

  deciduous_flex_mask <- rasterize(deciduous_flex_shp, raster_mask)
  deciduous_flex_mask <- terra::rast(deciduous_flex_mask)
  deciduous_flex_mask[!is.na(deciduous_flex_mask)] <- 1

  coniferous_only_mask <- rasterize(coniferous_only_shp, raster_mask)
  coniferous_only_mask <- terra::rast(coniferous_only_mask)
  coniferous_only_mask[!is.na(coniferous_only_mask)] <- 1

  coniferous_flex_mask <- rasterize(coniferous_flex_shp, raster_mask)
  coniferous_flex_mask <- terra::rast(coniferous_flex_mask)
  coniferous_flex_mask[!is.na(coniferous_flex_mask)] <- 1

  # if (length(oak_only_shp$TFV) == 0){
  #   print("Oak composition NULL:skip")
  #   oak_only_mask <- rasterize(deciduous_only_shp, raster_mask)
  #   oak_only_mask <- terra::rast(oak_only_mask)
  #   values(oak_only_mask) <- 0
  # } else {
  #   oak_only_mask <- rasterize(oak_only_shp, raster_mask)
  #   oak_only_mask <- terra::rast(oak_only_mask)
  #   oak_only_mask[!is.na(oak_only_mask)] <- 1
  #
  #   # Oak Only
  #   save_basic_plot(plot_to_save = oak_only_mask,
  #                   dirname = masks_dir,
  #                   filename = sprintf("oak_only_mask_res_%s_m.png",
  #                                      resolution),
  #                   title = sprintf("oak_only_mask_res_%s_m",
  #                                   resolution))
  #   plot(oak_only_mask, main=sprintf("oak_only_mask_res_%s_m",
  #                                    resolution))
  #   save_envi_file(oak_only_mask,
  #                  sprintf("oak_only_mask_res_%s_m",
  #                          resolution),
  #                  masks_dir)
  # }
  #
  # if (length(beech_only_shp$TFV) == 0){
  #   print("Beech Composition NULL:skip")
  #   beech_only_mask <- rasterize(deciduous_only_shp, raster_mask)
  #   beech_only_mask <- terra::rast(beech_only_mask)
  #   values(beech_only_mask) <- 0
  # } else {
  #   beech_only_mask <- rasterize(beech_only_shp, raster_mask)
  #   beech_only_mask <- terra::rast(beech_only_mask)
  #   beech_only_mask[!is.na(beech_only_mask)] <- 1
  #
  #   # Beech Only
  #   save_basic_plot(plot_to_save = beech_only_mask,
  #                   dirname = masks_dir,
  #                   filename = sprintf("beech_only_mask_res_%s_m.png",
  #                                      resolution),
  #                   title = sprintf("beech_only_mask_res_%s_m",
  #                                   resolution))
  #   plot(beech_only_mask, main=sprintf("beech_only_mask_res_%s_m",
  #                                      resolution))
  #   save_envi_file(beech_only_mask,
  #                  sprintf("beech_only_mask_res_%s_m",
  #                          resolution),
  #                  masks_dir)
  # }

  # Plot and save
  # Full composition
  save_basic_plot(plot_to_save = full_composition_mask,
                  dirname = masks_dir,
                  filename = sprintf("full_composition_mask_res_%s_m.png",
                                     resolution),
                  title = sprintf("full_composition_mask_res_%s_m",
                                  resolution))
  plot(full_composition_mask, main=sprintf("full_composition_mask_res_%s_m",
                                           resolution))
  save_envi_file(full_composition_mask,
                 sprintf("full_composition_mask_res_%s_m",
                         resolution),
                 masks_dir)

  # Deciduous Flex
  save_basic_plot(plot_to_save = deciduous_flex_mask,
                  dirname = masks_dir,
                  filename = sprintf("deciduous_flex_mask_res_%s_m.png",
                                     resolution),
                  title = sprintf("deciduous_flex_mask_res_%s_m",
                                  resolution))
  plot(deciduous_flex_mask, main=sprintf("deciduous_flex_mask_res_%s_m",
                                         resolution))
  save_envi_file(deciduous_flex_mask,
                 sprintf("deciduous_flex_mask_res_%s_m",
                         resolution),
                 masks_dir)

  # Deciduous Only
  save_basic_plot(plot_to_save = deciduous_only_mask,
                  dirname = masks_dir,
                  filename = sprintf("deciduous_only_mask_res_%s_m.png",
                                     resolution),
                  title = sprintf("deciduous_only_mask_res_%s_m",
                                  resolution))
  plot(deciduous_only_mask, main=sprintf("deciduous_only_mask_res_%s_m",
                                         resolution))
  save_envi_file(deciduous_only_mask,
                 sprintf("deciduous_only_mask_res_%s_m",
                         resolution),
                 masks_dir)

  cat(paste0("Full composition, deciduous only and flex,",
             "and coniferous only and flex masks",
             "have been successfully saved at :",
             masks_dir, "\n"))

  # Coniferous Flex
  save_basic_plot(plot_to_save = coniferous_flex_mask,
                  dirname = masks_dir,
                  filename = sprintf("coniferous_flex_mask_res_%s_m.png",
                                     resolution),
                  title = sprintf("coniferous_flex_mask_res_%s_m",
                                  resolution))
  plot(coniferous_flex_mask, main=sprintf("coniferous_flex_mask_res_%s_m",
                                          resolution))
  save_envi_file(coniferous_flex_mask,
                 sprintf("coniferous_flex_mask_res_%s_m",
                         resolution),
                 masks_dir)

  # Coniferous Only
  save_basic_plot(plot_to_save = coniferous_only_mask,
                  dirname = masks_dir,
                  filename = sprintf("coniferous_only_mask_res_%s_m.png",
                                     resolution),
                  title = sprintf("coniferous_only_mask_res_%s_m",
                                  resolution))
  plot(coniferous_only_mask, main=sprintf("coniferous_only_mask_res_%s_m",
                                          resolution))
  save_envi_file(coniferous_only_mask,
                 sprintf("coniferous_only_mask_res_%s_m",
                         resolution),
                 masks_dir)

  cat(paste0("Composition masks ",
             "have been successfully saved at :",
             masks_dir, "\n"))

  return(list(full_composition_mask = full_composition_mask,
              deciduous_only_mask = deciduous_only_mask,
              deciduous_flex_mask = deciduous_flex_mask,
              coniferous_only_mask = coniferous_only_mask,
              coniferous_flex_mask = coniferous_flex_mask
              # oak_only_mask = oak_only_mask,
              # beech_only_mask = beech_only_mask
  )
  )
}

#' Mask Site Edges
#'
#' This function masks a reflectance raster with the specified site edges.
#'
#' @param reflectance A raster object containing the reflectance data.
#' @param site_edges_path A character string specifying the path to the
#' shapefile containing the edges of the site in UTM format.
#' @param masks_dir A character string specifying the directory where the
#' generated masks will be saved.
#' @param resolution The resolution of the raster masks to be created.
#' @return A masked reflectance raster.
#' @export
mask_site_edges <- function(reflectance,
                            site_edges_path, # utm format
                            masks_dir,
                            resolution){
  site_edges_mask <- terra::vect(site_edges_path)
  site_edges_mask <- mask(reflectance, site_edges_mask)

  viz <- site_edges_mask
  viz[is.na(viz)] <- 0
  viz[viz == 1] <- NA

  save_basic_plot(plot_to_save = site_edges_mask,
                  dirname = masks_dir,
                  filename = sprintf("site_edges_mask_res_%s_m.png",
                                     resolution),
                  title = sprintf("site_edges_mask_res_%s_m", resolution))
  plot(site_edges_mask, main=sprintf("site_edges_mask_res_%s_m", resolution))

  save_envi_file(viz, sprintf("site_edges_mask_res_%s_m_viz", resolution),
                 masks_dir)
  save_envi_file(site_edges_mask, sprintf("site_edges_mask_res_%s_m",
                                          resolution),
                 masks_dir)

  cat("Site edges mask has been successfully saved at :",
      masks_dir, "\n")
  return(site_edges_mask)
}

#' Mask Clouds
#'
#' This function masks a reflectance raster with the specified cloud mask.
#'
#' @param reflectance A raster object containing the reflectance data.
#' @param cloud_path A character string specifying the cloud mask raster's path.
#' @param site_edges A character string specifying the path to the shapefile
#' containing the edges of the site in UTM format.
#' @param masks_dir A character string specifying the directory where the
#' generated masks will be saved.
#' @param resolution The resolution of the raster masks to be created.
#' @return A masked reflectance raster.
#' @export
mask_clouds <- function(reflectance,
                        cloud_path,
                        site_edges,
                        masks_dir,
                        resolution){
  cloud_mask <- terra::rast(cloud_path)


  cloud_mask <- terra::project(cloud_mask,
                               site_edges)

  cloud_mask <- reflectance * cloud_mask

  save_basic_plot(plot_to_save = cloud_mask,
                  dirname = masks_dir,
                  filename = sprintf("cloud_mask_res_%s_m.png", resolution),
                  title = sprintf("cloud_mask_res_%s_m", resolution))
  plot(cloud_mask, main=sprintf("cloud_mask_res_%s_m", resolution))

  save_envi_file(cloud_mask, sprintf("cloud_mask_res_%s_m", resolution),
                 masks_dir)

  cat("Clouds mask has been successfully saved at :",
      masks_dir, "\n")
  return(cloud_mask)
}

#' Save CHM
#'
#' This function saves the Canopy Height Model (CHM).
#'
#' @param chm A raster object containing the Canopy Height Model.
#' @param masks_dir A character string specifying the directory where
#' the CHM will be saved.
#' @return The saved CHM.
#' @export
save_chm <- function(chm,
                     masks_dir){

  # chm <- dts_charged - dtm_charged
  save_basic_plot(plot_to_save = chm,
                  dirname = masks_dir,
                  filename = "chm.png",
                  title = "chm")
  save_envi_file(chm, "chm", masks_dir)
  cat("CHM has been successfully saved at :",
      masks_dir, "\n")
  return(chm)
}

#' Mask CHM by Threshold
#'
#' This function masks the Canopy Height Model (CHM) based on
#' a specified threshold.
#'
#' @param chm A raster object containing the Canopy Height Model.
#' @param threshold The threshold value for masking the CHM.
#' @param masks_dir A character string specifying the directory where the
#' generated masks will be saved.
#' @param resolution The resolution of the raster masks to be created.
#' @return A masked CHM.
#' @export
mask_chm_thresh <- function(chm,
                            threshold,
                            masks_dir,
                            resolution){
  chm_thresh <- chm > threshold

  save_basic_plot(plot_to_save = chm_thresh,
                  dirname = masks_dir,
                  filename = sprintf("chm_thresholded_%d_m_res_%s_m.png",
                                     threshold,
                                     resolution),
                  title = sprintf("chm_thresholded_%d_m_res_%s_m",
                                  threshold,
                                  resolution))
  plot(chm_thresh,
       main=sprintf("chm_thresholded_%d_m_res_%s_m",
                    threshold,
                    resolution))

  save_envi_file(chm_thresh,
                 sprintf("chm_thresholded_%d_m_res_%s_m",
                         threshold,
                         resolution),
                 masks_dir)
  cat(sprintf("CHM thresholded %d m has been successfully saved at :",
              threshold),
      masks_dir, "\n")
  return(chm_thresh)
}

#' Project CHM by Threshold
#'
#' This function projects the Canopy Height Model (CHM) thresholded mask onto
#' a reference raster with correct coordinates.
#'
#' @param chm_thresh A raster object containing the thresholded CHM mask.
#' @param raster_with_right_coordinates A raster object with the correct
#' coordinates for projection.
#' @param masks_dir A character string specifying the directory where the
#' generated masks will be saved.
#' @param resolution The resolution of the raster masks to be created.
#' @return A projected CHM thresholded mask.
#' @export
project_chm_thresh <- function(chm_thresh,
                               raster_with_right_coordinates,
                               masks_dir,
                               resolution){
  average_chm_thresh <- terra::project(chm_thresh,
                                       raster_with_right_coordinates,
                                       method = 'average')
  writeRaster(average_chm_thresh,
              filename = file.path(masks_dir,
                                   "gap_fraction.tif"),
              filetype = "GTiff",
              overwrite = T)
  save_basic_plot(plot_to_save = average_chm_thresh,
                  dirname = masks_dir,
                  filename = sprintf(
                    "chm_1_m_thresholded_average_to_res_%s_m.png", resolution),
                  title = sprintf("chm_1_m_thresholded_average_to_res_%s_m",
                                  resolution))
  plot(average_chm_thresh,
       main=sprintf("chm_1_m_thresholded_average_to_res_%s_m", resolution))

  save_envi_file(average_chm_thresh,
                 sprintf("chm_1_m_thresholded_average_to_res_%s_m", resolution),
                 masks_dir)
  cat("CHM thresholded projected has been successfully saved at :",
      masks_dir, "\n")
  return(average_chm_thresh)
}

#' Mask CHM Thresholded by Majority Project
#'
#' This function masks the Canopy Height Model (CHM) thresholded mask by
#' majority projection and saves the result.
#'
#' @param chm A raster object containing the Canopy Height Model.
#' @param average_chm_thresh A raster object containing the
#' thresholded CHM mask.
#' @param percentage The percentage of pixels to keep based on
#' majority projection.
#' @param masks_dir A character string specifying the directory where the
#' generated masks will be saved.
#' @param resolution The resolution of the raster masks to be created.
#' @return The masked CHM.
#' @export
mask_majority_project_chm_thresh <- function(chm,
                                             average_chm_thresh,
                                             percentage,
                                             masks_dir,
                                             resolution){

  # Majority
  average_chm_thresh_percentage <- average_chm_thresh > percentage
  save_basic_plot(plot_to_save = average_chm_thresh_percentage,
                  dirname = masks_dir,
                  filename = sprintf(
                    "average_chm_thresholded_%s_p_kept_res_%s_m.png",
                    percentage*100, resolution),
                  title = sprintf("average_chm_thresholded_%s_p_kept_res_%s_m",
                                  percentage*100, resolution))
  plot(average_chm_thresh_percentage,
       main = sprintf("average_chm_thresholded_%s_p_kept_res_%s_m",
                      percentage*100, resolution))

  save_envi_file(average_chm_thresh_percentage,
                 sprintf("average_chm_thresholded_%s_p_kept_res_%s_m",
                         percentage*100, resolution),
                 masks_dir)


  # Project CHM
  average_chm <- terra::project(chm,
                                average_chm_thresh,
                                method = 'average')
  chm_val <- values(average_chm)
  chm_val[chm_val < 0] <- 0

  index <- average_chm_thresh_percentage == 1
  chm_final_mask <- average_chm
  chm_final_mask[!index] <- NA

  save_basic_plot(plot_to_save = chm_final_mask,
                  dirname = masks_dir,
                  filename = sprintf(
                    "chm_masked_with_average_chm_thresholded_%s_p_kept_res_%s_m.png",
                    percentage*100, resolution),
                  title = sprintf(
                    "chm_masked_with_average_chm_thresholded_%s_p_kept_res_%s_m",
                    percentage*100, resolution))
  plot(chm_final_mask,
       main = sprintf(
         "chm_masked_with_average_chm_thresholded_%s_p_kept_res_%s_m",
         percentage*100, resolution))

  save_envi_file(chm_final_mask,
                 sprintf(
                   "chm_masked_with_average_chm_thresholded_%s_p_kept_res_%s_m",
                   percentage*100, resolution),
                 masks_dir)

  # Comparison to assess the mask
  chm_final_mask_val <- values(chm_final_mask)
  chm_final_mask_val[chm_final_mask_val < 0] <- 0

  chm_list <- list(chm_val, chm_final_mask_val)
  labs <- c("chm", "chm_final_mask")
  plot_histogram(vars_list = chm_list,
                 title = sprintf(
                   "CHM vs CHM masked with threshold majorated %s p res %s m",
                   percentage*100, resolution),
                 xlab = "heights",
                 var_labs = labs,
                 dirname = masks_dir,
                 filename = sprintf(
                   "hist_chm_vs_chm_masked_threshold_majorated_%s_p_kept_res_%s_m",
                   percentage*100, resolution))
  return(average_chm_thresh_percentage)
}

#' Create Vegetation Forest Mask
#'
#' This function generates a vegetation forest mask based on LiDAR data,
#' Sentinel-2 images, and shapefiles.
#'
#' @param data_dir Directory containing data files.
#' @param site Name of the site where data is collected
#' (options: "Mormal", "Blois", "Aigoual").
#' @param shapefiles_dir Directory containing shapefiles.
#' @param masks_dir Directory to save generated masks.
#' @param results_path Path to the directory containing results.
#' @param resolution Resolution of the masks in meters.
#' @return A list containing the final vegetation forest mask.
#' @export
create_vegetation_forest_mask <- function(data_dir,
                                          site,
                                          shapefiles_dir,
                                          masks_dir,
                                          results_path,
                                          resolution){

  utm_init <- file.path(shapefiles_dir, "utm_init.shp") # init final
  bdforet_2 <- file.path(shapefiles_dir, "bdforet_2.shp")

  if (site == "Mormal"){
    dateAcq <- '2021-06-14' # yyyy-mm-dd format mandatory
  } else if (site == "Blois"){
    dateAcq <- '2021-06-14' # yyyy-mm-dd format mandatory
  } else if (site == "Aigoual"){
    dateAcq <- '2021-07-11' # yyyy-mm-dd format mandatory
  } else{
    stop("Error: Site must be Mormal, Blois, or Aigoual.\n")
  }

  result_path <- paste(results_dir, site, sep = '/')

  # Calculate CHM from DTS and DTM
  chm <- terra::rast(file.path(results_path,
                               "LiDAR/chm/res_1_m/chm.tif"))

  # Keep deciduous-only and deciduous-flex forest
  composition_masks <- mask_forest_composition_area(utm_init,
                                                    bdforet_2,
                                                    masks_dir,
                                                    resolution = 10)
  full_composition_mask <- composition_masks$full_composition_mask
  deciduous_only_mask <- composition_masks$deciduous_only_mask
  deciduous_flex_mask <- composition_masks$deciduous_flex_mask
  coniferous_only_mask <- composition_masks$coniferous_only_mask
  coniferous_flex_mask <- composition_masks$coniferous_flex_mask
  # oak_only_mask <- composition_masks$oak_only_mask
  # beech_only_mask <- composition_masks$beech_only_mask

  s2_creation_directory <- paste(data_dir,
                                 site,
                                 'S2_Images',
                                 dateAcq,
                                 sep = "/")
  dir.create(path = s2_creation_directory,
             showWarnings = FALSE,
             recursive = TRUE)

  # S-2 Pre-Processing: Cloud Mask, Reflectance
  results <- preprocess_S2_no_dl(
                                 # dateAcq,
                                 utm_init,
                                 s2_creation_directory,
                                 result_path,
                                 resolution = resolution,
                                 S2source = 'SAFE',
                                 saveRaw = TRUE)

  cloud_path <- results$Cloud_File
  refl_path <- results$refl_path # Reflectance

  reflectance <- terra::rast(refl_path)
  reflectance <- subset(reflectance, 1)
  reflectance[] <- 1

  # Site edges Mask
  site_edges_mask <- mask_site_edges(reflectance,
                                     utm_init,
                                     masks_dir,
                                     resolution = 10)

  # Cloud Mask

  cloud_mask <- mask_clouds(reflectance,
                            cloud_path,
                            site_edges_mask,
                            masks_dir,
                            resolution = 10)

  # CHM

  chm <- save_chm(chm, masks_dir)

  # Threshold

  chm_thresh <- mask_chm_thresh(chm,
                                threshold = 2,
                                masks_dir,
                                resolution = 10)

  # Project

  average_chm_thresh <- project_chm_thresh(chm_thresh,
                                           reflectance,
                                           masks_dir,
                                           resolution = 10)

  # Majority

  # percs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  percentage <- 0.9
  average_chm_thresh_perc <- mask_majority_project_chm_thresh(chm,
                                                              average_chm_thresh,
                                                              percentage,
                                                              masks_dir,
                                                              resolution = 10)

  # Final Mask

  # Initial artifacts
  artifacts <- cloud_mask * site_edges_mask
  save_envi_file(artifacts,
                 sprintf(paste("artifacts",
                               "res_%s_m_viz",
                               sep = "_"),
                         resolution),
                 masks_dir)
  save_basic_plot(plot_to_save = artifacts,
                  dirname = masks_dir,
                  filename = sprintf(paste("artifacts",
                                           "res_%s_m_viz.png",
                                           sep = "_"),
                                     resolution),
                  title = sprintf(paste("artifacts",
                                        "res_%s_m_viz",
                                        sep = "_"),
                                  resolution))

  # Assign unique names to each SpatRaster object
  names(full_composition_mask) <- "full_composition"
  names(deciduous_only_mask) <- "deciduous_only"
  names(deciduous_flex_mask) <- "deciduous_flex"
  names(coniferous_only_mask) <- "coniferous_only"
  names(coniferous_flex_mask) <- "coniferous_flex"
  # names(oak_only_mask) <- "oak_only"
  # names(beech_only_mask) <- "beech_only"

  # Composition list
  compositions <- list(full_composition_mask,
                       deciduous_only_mask,
                       deciduous_flex_mask,
                       coniferous_only_mask,
                       coniferous_flex_mask
                       # oak_only_mask,
                       # beech_only_mask
  )

  for (composition in compositions) {
    if (identical(values(composition), values(full_composition_mask))
        && names(composition) == "full_composition") {
      type <- "full_composition"
    } else if (identical(values(composition), values(deciduous_only_mask))
               && names(composition) == "deciduous_only") {
      type <- "deciduous_only"
    } else if (identical(values(composition), values(deciduous_flex_mask))
               && names(composition) == "deciduous_flex") {
      type <- "deciduous_flex"
    } else if (identical(values(composition), values(coniferous_only_mask))
               && names(composition) == "coniferous_only") {
      type <- "coniferous_only"
    } else if (identical(values(composition), values(coniferous_flex_mask))
               && names(composition) == "coniferous_flex") {
      type <- "coniferous_flex"
      # } else if (identical(values(composition), values(oak_only_mask))
      #            && names(composition) == "oak_only") {
      #   type <- "oak_only"
      # } else if (identical(values(composition), values(beech_only_mask))
      #            && names(composition) == "beech_only") {
      #   type <- "beech_only"
    } else {
      stop("Error: Composition choice is not amongst right ones.\n")
    }

    # Print the type determined for the current composition
    print(paste("Composition type:", type))
    composition_resampled <- terra::resample(composition,
                                             site_edges_mask,
                                             method="near")

    # Artifacts + Composition
    artifacts_composition <- cloud_mask * site_edges_mask * composition_resampled
    save_envi_file(artifacts_composition,
                   sprintf(paste("artifacts",
                                 type,
                                 "res_%s_m_viz",
                                 sep = "_"),
                           resolution),
                   masks_dir)
    save_basic_plot(plot_to_save = artifacts_composition,
                    dirname = masks_dir,
                    filename = sprintf(paste("artifacts",
                                             type,
                                             "res_%s_m_viz.png",
                                             sep = "_"),
                                       resolution),
                    title = sprintf(paste("artifacts",
                                          type,
                                          "res_%s_m_viz",
                                          sep = "_"),
                                    resolution))

    final_mask <- (cloud_mask * site_edges_mask
                   * composition_resampled * average_chm_thresh_perc)
    save_envi_file(final_mask,
                   sprintf(paste("artifacts",
                                 type,
                                 "low_vegetation_majority_%s_p_res_%s_m_viz",
                                 sep = "_"),
                           percentage*100,
                           resolution),
                   masks_dir)
    save_basic_plot(plot_to_save = final_mask,
                    dirname = masks_dir,
                    filename = sprintf(paste("artifacts",
                                             type,
                                             "low_vegetation_majority_%s_p_res_%s_m_viz.png",
                                             sep = "_"),
                                       percentage*100,
                                       resolution),
                    title = sprintf(paste("artifacts",
                                          type,
                                          "low_vegetation_majority_%s_p_res_%s_m_viz",
                                          sep = "_"),
                                    percentage*100,
                                    resolution))
    final_mask[final_mask == 0] <- NA
    plot(final_mask)
    save_envi_file(final_mask,
                   sprintf(paste("artifacts",
                                 type,
                                 "low_vegetation_majority_%s_p_res_%s_m",
                                 sep = "_"),
                           percentage*100,
                           resolution),
                   masks_dir)
    save_basic_plot(plot_to_save = final_mask,
                    dirname = masks_dir,
                    filename = sprintf(paste("artifacts",
                                             type,
                                             "low_vegetation_majority_%s_p_res_%s_m.png",
                                             sep = "_"),
                                       percentage*100,
                                       resolution),
                    title = sprintf(paste("artifacts",
                                          type,
                                          "low_vegetation_majority_%s_p_res_%s_m",
                                          sep = "_"),
                                    percentage*100,
                                    resolution))
  }
}

#' Create Heterogeneity Quantiles
#'
#' This function calculates heterogeneity quantiles based on error values
#' and creates masks accordingly.
#'
#' @param error Error values.
#' @param err Name of the error parameter.
#' @param site Name of the site.
#' @param masks_dir Directory to save generated masks.
#' @param composition_mask_dir Directory containing composition masks.
#' @param choice Choice of quantile calculation method
#' ("equal_intervals" or "deciles").
#' @param resolution Resolution of the masks in meters.
#' @return No direct return value.
#' @export
create_heterogeneity_quantiles <- function(error,
                                           err,
                                           site,
                                           masks_dir,
                                           composition_mask_dir,
                                           choice,
                                           inc = 0.1,
                                           resolution = 10){

  # Equal Intervals
  if (choice == "equal_intervals"){
    class_filename <- "Equal_Intervals"
    class_dir <- file.path(masks_dir, class_filename)
    dir.create(path = class_dir, showWarnings = F, recursive = T)
    choice_equal_deciles_for_plot <- "Interval"

    # Calculate the range of error values
    error_min <- min(values(error), na.rm = TRUE)
    error_max <- max(values(error), na.rm = TRUE)

    # Determine the range
    range <- c(0, 1-inc)
    midpoints <- seq(range[1] + inc/2, range[2] + inc/2, by = inc)

    # Adjust dir
    nb_filename <- sprintf("%s_Intervals/", num_intervals)
    class_dir <- file.path(class_dir, nb_filename)
    dir.create(path = class_dir, showWarnings = F, recursive = T)

    # Calculate interval width
    interval_width <- (error_max - error_min) / num_intervals

    for (i in 1:num_intervals) {
      low_value <- error_min + interval_width * (i - 1)
      high_value <- error_min + interval_width * i

      low_quantile <- inc * (i-1)

      cat("Interval", i, "\n")
      cat("Low Value", low_value, "\n")
      cat("High Value", high_value, "\n")

      heter_raster <- mask(error,
                           mask = error >= high_value | error < low_value,
                           maskvalue = 1)

      plot(heter_raster,
           main = paste(sprintf("%s heter raster res %s m",
                                err, resolution),
                        low_quantile, "to", low_quantile + inc))


      writeRaster(heter_raster,
                  paste(class_dir, sprintf("%s_heter_raster_res_%s_m_",
                                           err, resolution),
                        low_quantile, "_",
                        low_quantile + inc, ".tif", sep = ''),
                  overwrite = T)

      heter_binary_mask <- error >= high_value | error < low_value

      plot(heter_binary_mask,
           main = paste(sprintf("%s heter binary mask res %s m",
                                err, resolution), low_quantile,
                        "to", low_quantile + inc))

      writeRaster(heter_raster,
                  paste(class_dir, sprintf("%s_heter_binary_mask_res_%s_m_",
                                           err, resolution),
                        low_quantile, "_", low_quantile + inc,
                        ".tif", sep = ''),
                  overwrite = T)
    }
  }
  # Deciles
  else if (choice == "deciles"){
    class_filename <- "Deciles"
    class_dir <- file.path(masks_dir, class_filename)
    dir.create(path = class_dir, showWarnings = F, recursive = T)
    choice_equal_deciles_for_plot <- "Quantile"

    # Determine the number of quantiles
    range <- c(0, 1-inc)
    midpoints <- seq(range[1] + inc/2, range[2] + inc/2, by = inc)
    quantile_range <- seq(range[1], range[2], by = inc)
    quantiles_list <- list()

    # Adjust dir
    nb_filename <- sprintf("%s_Quantiles/", length(quantile_range))
    class_dir <- file.path(class_dir, nb_filename)
    dir.create(path = class_dir, showWarnings = F, recursive = T)

    for (low_quantile in quantile_range) {
      quantiles <- quantile(values(error),
                            probs = c(low_quantile, low_quantile + inc),
                            na.rm = TRUE)
      quantiles_list[[as.character(low_quantile)]] <- quantiles
      print(quantiles)
      print(summary(quantiles))

      heter_raster <- mask(error,
                           mask = error >= quantiles[2]
                           |error < quantiles[1], maskvalue = 1)

      plot(heter_raster,
           main = paste(sprintf("%s heter raster res %s m",
                                err, resolution),
                        low_quantile, "to", low_quantile + inc))

      heter_binary_mask <- error >= quantiles[2] | error < quantiles[1]

      plot(heter_binary_mask,
           main = paste(sprintf("%s heter binary mask res %s m",
                                err, resolution), low_quantile,
                        "to", low_quantile + inc))

      low_quantile <- round(low_quantile, digits = 2)
      high_quantile <- round(low_quantile + inc, digits = 2)

      writeRaster(heter_raster,
                  paste(class_dir, sprintf("%s_heter_raster_res_%s_m_",
                                           err, resolution),
                        low_quantile, "_",
                        high_quantile, ".tif", sep = ''),
                  overwrite = T)

      writeRaster(heter_raster,
                  paste(class_dir, sprintf("%s_heter_binary_mask_res_%s_m_",
                                           err, resolution),
                        low_quantile, "_",
                        high_quantile, ".tif", sep = ''),
                  overwrite = T)
    }
  }
  else {
    stop("choice_equal_deciles is not Equal Interval or Decile")
  }

  # Apply quantiles to LiDAR and Sentinel-2 LAI
  lai_lidar <- terra::rast(file.path(composition_mask_dir,
                                     "lidarlai_res_10_m.tif"))
  lai_s2 <- terra::rast(file.path(composition_mask_dir,
                                  "s2lai_res_10_m.tif"))

  for (low_value in seq(range[1], range[2], by = inc)) {
    # Apply deciles or equal intervals
    low_value <- round(low_value, digits = 2)
    high_value <- round(low_value + inc, digits = 2)
    intervals <- terra::rast(file.path(class_dir,
                                       paste0(sprintf("%s_heter_binary_mask_res_%s_m_",
                                                      err, resolution),
                                              low_value, "_",
                                              high_value, ".tif")))
    # LiDAR
    pai_resample <- terra::project(lai_lidar, intervals)
    pai_masked <- mask(pai_resample, intervals)
    plot(pai_masked, main = paste("pai lidar", low_value, "to", high_value))
    writeRaster(pai_masked,
                paste(class_dir, sprintf("%s_pai_masked_res_%s_m_",
                                         err,
                                         resolution),
                      low_value, "_", high_value,
                      ".tif", sep = ''),
                overwrite = T)

    # S-2
    lai_s2_resample <- terra::project(lai_s2, intervals)
    lai_s2_mask <- mask(lai_s2_resample, intervals)
    plot(lai_s2_mask, main = paste("lai s2", low_value, "to", high_value))
    writeRaster(lai_s2_mask,
                paste(class_dir, sprintf("%s_lai_s2_masked_res_%s_m_",
                                         err,
                                         resolution),
                      low_value, "_", high_value,
                      ".tif", sep = ''),
                overwrite = T)

  }
}







create_cv_lad_quantiles <- function(error,
                                    err,
                                    min_depth,
                                    site,
                                    masks_dir,
                                    composition_mask_dir,
                                    choice,
                                    inc = 0.1,
                                    resolution = 10){

  # Equal Intervals
  if (choice == "equal_intervals"){
    class_filename <- "Equal_Intervals"
    class_dir <- file.path(masks_dir, class_filename)
    dir.create(path = class_dir, showWarnings = F, recursive = T)
    choice_equal_deciles_for_plot <- "Interval"

    # Calculate the range of error values
    error_min <- min(values(error), na.rm = TRUE)
    error_max <- max(values(error), na.rm = TRUE)

    # Determine the range
    range <- c(0, 1-inc)
    midpoints <- seq(range[1] + inc/2, range[2] + inc/2, by = inc)

    # Adjust dir
    nb_filename <- sprintf("%s_Intervals/", num_intervals)
    class_dir <- file.path(class_dir, nb_filename)
    dir.create(path = class_dir, showWarnings = F, recursive = T)

    # Calculate interval width
    interval_width <- (error_max - error_min) / num_intervals

    for (i in 1:num_intervals) {
      low_value <- error_min + interval_width * (i - 1)
      high_value <- error_min + interval_width * i

      low_quantile <- inc * (i-1)

      cat("Interval", i, "\n")
      cat("Low Value", low_value, "\n")
      cat("High Value", high_value, "\n")

      heter_raster <- mask(error,
                           mask = error >= high_value | error < low_value,
                           maskvalue = 1)

      plot(heter_raster,
           main = paste(sprintf("%s heter raster res %s m",
                                err, resolution),
                        low_quantile, "to", low_quantile + inc))


      writeRaster(heter_raster,
                  paste(class_dir, sprintf("%s_heter_raster_res_%s_m_",
                                           err, resolution),
                        low_quantile, "_",
                        low_quantile + inc, ".tif", sep = ''),
                  overwrite = T)

      heter_binary_mask <- error >= high_value | error < low_value

      plot(heter_binary_mask,
           main = paste(sprintf("%s heter binary mask res %s m",
                                err, resolution), low_quantile,
                        "to", low_quantile + inc))

      writeRaster(heter_raster,
                  paste(class_dir, sprintf("%s_heter_binary_mask_res_%s_m_",
                                           err, resolution),
                        low_quantile, "_", low_quantile + inc,
                        ".tif", sep = ''),
                  overwrite = T)
    }
  }
  # Deciles
  else if (choice == "deciles"){
    class_filename <- "Deciles"
    class_dir <- file.path(masks_dir, class_filename)
    dir.create(path = class_dir, showWarnings = F, recursive = T)
    choice_equal_deciles_for_plot <- "Quantile"

    # Determine the number of quantiles
    range <- c(0, 1-inc)
    midpoints <- seq(range[1] + inc/2, range[2] + inc/2, by = inc)
    quantile_range <- seq(range[1], range[2], by = inc)
    quantiles_list <- list()

    # Adjust dir
    nb_filename <- sprintf("%s_Quantiles/", length(quantile_range))
    class_dir <- file.path(class_dir, nb_filename)
    dir.create(path = class_dir, showWarnings = F, recursive = T)

    for (low_quantile in quantile_range) {
      quantiles <- quantile(values(error),
                            probs = c(low_quantile, low_quantile + inc),
                            na.rm = TRUE)
      quantiles_list[[as.character(low_quantile)]] <- quantiles
      print(quantiles)
      print(summary(quantiles))

      heter_raster <- mask(error,
                           mask = error >= quantiles[2]
                           |error < quantiles[1], maskvalue = 1)

      plot(heter_raster,
           main = paste(sprintf("%s heter raster res %s m",
                                err, resolution),
                        low_quantile, "to", low_quantile + inc))

      heter_binary_mask <- error >= quantiles[2] | error < quantiles[1]

      plot(heter_binary_mask,
           main = paste(sprintf("%s heter binary mask res %s m",
                                err, resolution), low_quantile,
                        "to", low_quantile + inc))

      low_quantile <- round(low_quantile, digits = 2)
      high_quantile <- round(low_quantile + inc, digits = 2)

      writeRaster(heter_raster,
                  paste(class_dir, sprintf("%s_heter_raster_res_%s_m_",
                                           err, resolution),
                        low_quantile, "_",
                        high_quantile, ".tif", sep = ''),
                  overwrite = T)

      writeRaster(heter_raster,
                  paste(class_dir, sprintf("%s_heter_binary_mask_res_%s_m_",
                                           err, resolution),
                        low_quantile, "_",
                        high_quantile, ".tif", sep = ''),
                  overwrite = T)
    }
  }
  else {
    stop("choice_equal_deciles is not Equal Interval or Decile")
  }

  # Apply quantiles to LiDAR and Sentinel-2 LAI
  pad_filename <- paste0("PAD_",
                         min_depth,
                         "_40.tif")
  lai_lidar <- terra::rast(file.path(composition_mask_dir,
                                     "PAD_Profiles_updated",
                                     pad_filename))
  lai_s2 <- terra::rast(file.path(composition_mask_dir,
                                  "s2lai_res_10_m.tif"))

  for (low_value in seq(range[1], range[2], by = inc)) {
    # Apply deciles or equal intervals
    low_value <- round(low_value, digits = 2)
    high_value <- round(low_value + inc, digits = 2)
    intervals <- terra::rast(file.path(class_dir,
                                       paste0(sprintf("%s_heter_binary_mask_res_%s_m_",
                                                      err, resolution),
                                              low_value, "_",
                                              high_value, ".tif")))
    # LiDAR
    pai_resample <- terra::project(lai_lidar, intervals)
    pai_masked <- mask(pai_resample, intervals)
    plot(pai_masked, main = paste("pai lidar", low_value, "to", high_value))
    writeRaster(pai_masked,
                paste(class_dir, sprintf("%s_pai_masked_res_%s_m_",
                                         err,
                                         resolution),
                      low_value, "_", high_value,
                      ".tif", sep = ''),
                overwrite = T)

    # S-2
    lai_s2_resample <- terra::project(lai_s2, intervals)
    lai_s2_mask <- mask(lai_s2_resample, intervals)
    plot(lai_s2_mask, main = paste("lai s2", low_value, "to", high_value))
    writeRaster(lai_s2_mask,
                paste(class_dir, sprintf("%s_lai_s2_masked_res_%s_m_",
                                         err,
                                         resolution),
                      low_value, "_", high_value,
                      ".tif", sep = ''),
                overwrite = T)

  }
}












#' Apply and Save Masks
#'
#' This function applies masks to a raster and saves the masked rasters
#' in different directories based on mask types.
#'
#' @param raster Raster to be masked.
#' @param raster_basename Basename of the raster file.
#' @param masks_dir Directory containing masks.
#' @param metrics_dir Directory to save masked rasters based on mask types.
#' @param pad_bool Boolean indicating if studied metric is a PAD profile.
#' @return No direct return value.
#' @export
apply_and_save_masks <- function(raster,
                                 raster_basename,
                                 masks_dir,
                                 metrics_dir,
                                 pad_classic = FALSE,
                                 pad_dtm = FALSE,
                                 pad_dsm = FALSE
                                 # cvlad_bool = FALSE
                                 ){

  # Open masks
  full_comp_mask <- terra::rast(
    file.path(masks_dir, "artifacts_full_composition_low_vegetation_majority_90_p_res_10_m.envi"))
  deciduous_flex_mask <- terra::rast(
    file.path(masks_dir, "artifacts_deciduous_flex_low_vegetation_majority_90_p_res_10_m.envi"))
  deciduous_only_mask <- terra::rast(
    file.path(masks_dir, "artifacts_deciduous_only_low_vegetation_majority_90_p_res_10_m.envi"))
  coniferous_flex_mask <- terra::rast(
    file.path(masks_dir, "artifacts_coniferous_flex_low_vegetation_majority_90_p_res_10_m.envi"))
  coniferous_only_mask <- terra::rast(
    file.path(masks_dir, "artifacts_coniferous_only_low_vegetation_majority_90_p_res_10_m.envi"))
  # oak_only_mask <- terra::rast(
  #   file.path(masks_dir, "artifacts_oak_only_low_vegetation_majority_90_p_res_10_m.envi"))
  # beech_only_mask <- terra::rast(
  #   file.path(masks_dir, "artifacts_beech_only_low_vegetation_majority_90_p_res_10_m.envi"))
  no_mask <- terra::rast(
    file.path(masks_dir, "chm_1_m_thresholded_average_to_res_10_m.envi"))

  # Project the raster to right coordinates and resolution
  projected_raster <- terra::project(raster, full_comp_mask)

  # Masks list
  composition_masks <- c("Full_Composition",
                         "Deciduous_Flex",
                         "Deciduous_Only",
                         "Coniferous_Flex",
                         "Coniferous_Only",
                         # "Oak_Only",
                         # "Beech_Only",
                         "Not_Masked")

  # raster_basename <- if (pad_bool &! norm_bool) file.path("PAD_Profiles", raster_basename)
  # raster_basename <- if (pad_bool && norm_bool) file.path("PAD_Profiles_updated", raster_basename)
  # else raster_basename
  # Initialize the base directory for the raster basename
  base_dir <- ""

  if (pad_classic) {
    base_dir <- "PAD_Profiles_Classic"
  } else if (!pad_classic & pad_dtm) {
    base_dir <- "PAD_Profiles_dtm_pit"
  } else if (!pad_classic & !pad_dtm & pad_dsm) {
    # base_dir <- "PAD_Profiles_updated_modifminz_NA"
    base_dir <- "PAD_Profiles_dsm_pit"
  } else if (!pad_classic & !pad_dtm & !pad_dsm) {
    base_dir <- ""
    # base_dir <- file.path("PAD_Profiles_updated", "CV_LAD")
  } else {
    base_dir <- ""
  }

  # Construct the full path for the raster basename
  raster_basename <- file.path(base_dir, raster_basename)

  for (composition_mask in composition_masks){
    if (!dir.exists(file.path(metrics_dir, composition_mask))) {
      dir.create(file.path(metrics_dir, composition_mask),
                 showWarnings = FALSE,
                 recursive = TRUE)
    }
    if (!dir.exists(file.path(metrics_dir,
                              composition_mask,
                              raster_basename))) {
      dir.create(file.path(metrics_dir,
                           composition_mask,
                           raster_basename),
                 showWarnings = FALSE,
                 recursive = TRUE)
    }
  }



  # Apply and save
  writeRaster(mask(projected_raster, full_comp_mask),
              filename = file.path(metrics_dir,
                                   "Full_Composition",
                                   raster_basename),
              filetype = "GTiff",
              overwrite = T)
  writeRaster(mask(projected_raster, deciduous_flex_mask),
              filename = file.path(metrics_dir,
                                   "Deciduous_Flex",
                                   raster_basename),
              filetype = "GTiff",
              overwrite = T)
  writeRaster(mask(projected_raster, deciduous_only_mask),
              filename = file.path(metrics_dir,
                                   "Deciduous_Only",
                                   raster_basename),
              filetype = "GTiff",
              overwrite = T)
  writeRaster(mask(projected_raster, coniferous_flex_mask),
              filename = file.path(metrics_dir,
                                   "Coniferous_Flex",
                                   raster_basename),
              filetype = "GTiff",
              overwrite = T)
  writeRaster(mask(projected_raster, coniferous_only_mask),
              filename = file.path(metrics_dir,
                                   "Coniferous_Only",
                                   raster_basename),
              filetype = "GTiff",
              overwrite = T)
  # writeRaster(mask(projected_raster, oak_only_mask),
  #             filename = file.path(metrics_dir,
  #                                  "Oak_Only",
  #                                  raster_basename),
  #             filetype = "GTiff",
  #             overwrite = T)
  # writeRaster(mask(projected_raster, beech_only_mask),
  #             filename = file.path(metrics_dir,
  #                                  "Beech_Only",
  #                                  raster_basename),
  #             filetype = "GTiff",
  #             overwrite = T)
  writeRaster(mask(projected_raster, no_mask),
              filename = file.path(metrics_dir,
                                   "Not_Masked",
                                   raster_basename),
              filetype = "GTiff",
              overwrite = T)
}
