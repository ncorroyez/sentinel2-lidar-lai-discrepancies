# ---
# title: "function_sentinel_2.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, 
# CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2023-11-07"
# ---

library(preprocS2)
library(prosail)
library(spinR)
library(raster)
library(stars)
library(terra)
library(ggplot2)
library(truncnorm)
library(bigRaster)

#' Downloads Sentinel-2 satellite images for a given acquisition date 
#' and region of interest.
#'
#' @param dateAcq Date of image acquisition in "YYYY-MM-DD" format.
#' @param path_vector Path to the vector file defining the region of interest.
#' @param output_directory Directory where the downloaded images will be saved.
#' 
#' @details This function checks if the images for the specified date 
#' have already been downloaded. 
#' If they haven't, it downloads them using the get_S2_L2A_Image function 
#' from the preprocS2 package.
#' 
#' @return If images are successfully downloaded, it prints a success message. 
#' If images have already been downloaded, it prints a message indicating 
#' that the process is skipped.
#' 
#' @export
download_S2_images <- function(dateAcq, path_vector, output_directory) {
  
  # Check if already downloaded
  date_without_hyphens <- gsub("-", "", dateAcq)
  subdirectories <- list.dirs(output_directory, # list all S2_Images
                              full.names = FALSE, recursive = FALSE)
  
  # Skip if yes
  if (any(grepl(date_without_hyphens, subdirectories))) {
    cat("Skipping download for date", dateAcq, 
        "because it has been already downloaded in", output_directory, ".\n")
  } else {
    Path_S2 <- get_S2_L2A_Image(l2a_path = output_directory, 
                                spatial_extent = path_vector, 
                                dateAcq = dateAcq,
                                DeleteL1C = TRUE, 
                                Sen2Cor = TRUE) # Fastest way to download
    cat("S2 image acquired at", dateAcq, 
        "has been successfully created and saved at :", output_directory, ".\n")
  }
}

#' Processes the downloaded Sentinel-2 images, including extraction, 
#' resampling, and stacking of data.
#'
#' @param output_directory Directory containing the downloaded 
#' Sentinel-2 images.
#' @param path_vector Path to the vector file defining the region of interest.
#' @param result_path Directory where the processed data will be saved.
#' @param resolution Resolution of the processed data.
#' @param S2source Source of the Sentinel-2 data (e.g., 'SAFE' or 'L2A').
#'
#' @details This function extracts data from the downloaded images, 
#' resamples it to the specified resolution, 
#' and stacks the data. It then returns the processed data object.
#'
#' @return Processed Sentinel-2 data object.
#' 
#' @export
process_S2_data <- function(output_directory, 
                            path_vector, 
                            result_path, 
                            resolution, 
                            S2source) {
  # Define raster path
  # Path_S2 <- file.path(output_directory, list.files(output_directory, 
  #                                                   pattern = 'FRE')) # FRE .SAFE
  
  Path_S2 <- file.path(output_directory, list.files(output_directory))
  
  # Create the result directory if it doesn't exist
  dir.create(path = result_path, showWarnings = FALSE, recursive = TRUE)
  
  # Extract, resample, and stack data
  S2obj <- preprocS2::extract_from_S2_L2A(
    Path_dir_S2 = Path_S2,
    path_vector = path_vector,
    S2source = S2source, # S2source
    resolution = resolution
  )
  cat("This image has been successfully processed 
      and the results are saved at :", Path_S2, "\n")
  # Update shapefile path if needed (reprojection)
  # forest_mask <- S2obj$forest_mask
  return(list(S2obj = S2obj,
              path_vector = S2obj$path_vector))
}

#' Generates a cloud mask from the Sentinel-2 stack.
#'
#' @param S2_Stack Stack of Sentinel-2 images.
#' @param Cloud_path Directory where the cloud mask will be saved.
#' @param S2source Source of the Sentinel-2 data (e.g., 'SAFE' or 'L2A').
#' @param saveRaw Logical indicating whether to save the raw cloud mask.
#'
#' @details This function generates a cloud mask from the Sentinel-2 stack using the save_cloud_s2 function from preprocS2 package.
#'
#' @return Object containing information about the generated cloud mask.
#'
#' @export
generate_cloud_mask <- function(S2_Stack, Cloud_path, S2source, saveRaw = TRUE) {
  # Write cloud mask and get filenames
  cloudmasks <- preprocS2::save_cloud_s2(
    S2_stars = S2_Stack,
    Cloud_path = Cloud_path,
    S2source = S2source,
    SaveRaw = saveRaw
  )
  cat("Cloud mask has been successfully created and saved at :", 
      Cloud_path, "\n")
  
  return(cloudmasks)
}

#' Writes the reflectance data to a file.
#'
#' @param S2_Stack Stack of Sentinel-2 images.
#' @param S2_Bands Information about the bands in the Sentinel-2 data.
#' @param refl_path Path where the reflectance data will be saved.
#'
#' @details This function writes the reflectance data to a file in ENVI format 
#' using the save_reflectance_s2 function from preprocS2 package.
#'
#' @export
write_reflectance <- function(S2_Stack, S2_Bands, refl_path) {
  # Save Reflectance file as ENVI image with BIL interleaves
  tile_S2 <- preprocS2::get_tile(S2_Bands$GRANULE)
  dateAcq_S2 <- preprocS2::get_dateAcq(S2product = dirname(dirname(S2_Bands$GRANULE)))
  
  preprocS2::save_reflectance_s2(
    S2_stars = S2_Stack, 
    Refl_path = refl_path,
    S2Sat = NULL, 
    tile_S2 = tile_S2, 
    dateAcq_S2 = dateAcq_S2,
    Format = 'ENVI', 
    datatype = 'Int16', 
    MTD = S2_Bands$metadata, 
    MTD_MSI = S2_Bands$metadata_MSI
  )
  cat("Reflectance has been successfully created and saved at :", 
      refl_path, "\n")
}

#' Reads a raster from a file using the brick() function.
#'
#' @param refl_path Path to the raster file.
#'
#' @return Raster object read from the file.
#'
#' @export
read_raster <- function(refl_path) {
  Refl <- brick(refl_path)
  return(Refl)
}

#' Computes spectral indices from the reflectance data.
#'
#' @param Refl Raster object containing reflectance data.
#' @param SensorBands Information about the bands in the Sentinel-2 data.
#' @param IndexList List of spectral indices to compute.
#' @param ReflFactor Factor to scale reflectance values.
#'
#' @details This function computes spectral indices from the reflectance data 
#' using the compute_S2SI_Raster function from the spinR package.
#'
#' @return List containing computed spectral indices.
#'
#' @export
compute_spectral_indices <- function(Refl, SensorBands, IndexList, ReflFactor) {
  Indices <- spinR::compute_S2SI_Raster(
    Refl = Refl,
    SensorBands = SensorBands,
    Sel_Indices = IndexList,
    ReflFactor = ReflFactor,
    StackOut = FALSE
  )
  return(Indices)
}

#' Saves computed spectral indices to files.
#'
#' @param Indices List containing computed spectral indices.
#' @param SI_path Path where the spectral indices will be saved.
#' @param S2_Bands Information about the bands in the Sentinel-2 data.
#'
#' @details This function saves the computed spectral indices to files 
#' in ENVI format.
#'
#' @export
save_spectral_indices <- function(Indices, 
                                  SI_path, 
                                  S2_Bands) {
  for (SpIndx in names(Indices$SpectralIndices)) {
    Index_Path <- file.path(SI_path, paste(basename(S2_Bands$GRANULE), 
                                           '_', SpIndx, sep = ''))
    stars::write_stars(st_as_stars(Indices$SpectralIndices[[SpIndx]]), 
                       dsn = Index_Path, driver =  "ENVI", type = 'Float32')
    # Write band name in HDR
    HDR <- read_ENVI_header(get_HDR_name(Index_Path))
    HDR$`band names` <- SpIndx
    write_ENVI_header(HDR = HDR, HDRpath = get_HDR_name(Index_Path))
  }
}

#' Updates the cloud mask based on radiometric filtering.
#'
#' @param CloudPath Directory containing the cloud mask.
#' @param NDVI_Thresh Threshold value for NDVI.
#' @param Indices List containing computed spectral indices.
#' @param cloudmasks Object containing information about the cloud mask.
#'
#' @details This function updates the cloud mask based on radiometric filtering
#' using the update_cloud_mask function.
#'
#' @return Path to the updated cloud mask.
#'
#' @export
update_cloud_mask <- function(CloudPath, NDVI_Thresh, Indices, cloudmasks) {
  Elim <- which(values(Indices$SpectralIndices[['NDVI']]) < NDVI_Thresh)
  CloudInit <- stars::read_stars(cloudmasks$BinaryMask)
  CloudInit$CloudMask_Binary[Elim] <- 0
  # Save the updated cloud mask
  Cloud_File <- file.path(CloudPath, 'CloudMask_Binary_Update')
  stars::write_stars(CloudInit, dsn = Cloud_File, 
                     driver = "ENVI", type = 'Byte')
  return(Cloud_File)
}

#' Preprocesses Sentinel-2 data.
#'
#' @param dateAcq Date of image acquisition in "YYYY-MM-DD" format.
#' @param path_vector Path to the vector file defining the region of interest.
#' @param output_directory Directory where the downloaded images are saved.
#' @param result_path Directory where the processed data will be saved.
#' @param resolution Resolution of the processed data.
#' @param S2source Source of the Sentinel-2 data (e.g., 'SAFE').
#' @param saveRaw Logical indicating whether to save the raw cloud mask.
#'
#' @details This function preprocesses Sentinel-2 data by downloading images,
#' processing them, computing spectral indices, and updating the cloud mask.
#'
#' @return List containing information about the processed data.
#'
#' @export
preprocess_S2 <- function(dateAcq, 
                          path_vector,
                          output_directory,
                          result_path, 
                          resolution = resolution, 
                          S2source,
                          saveRaw = TRUE) {
  # download_S2_images(dateAcq, path_vector, output_directory)
  process_S2_data <- process_S2_data(output_directory, path_vector, result_path, 
                                     resolution, S2source)
  print("1")
  S2obj <- process_S2_data$S2obj
  path_vector <- process_S2_data$path_vector
  
  results_site_path <- file.path(result_path,basename(S2obj$S2_Bands$GRANULE))
  dir.create(path = results_site_path,showWarnings = FALSE,recursive = TRUE)
  print("2")
  cloud_path <- file.path(results_site_path,
                          'CloudMask',
                          sprintf("res_%d_m", resolution))
  dir.create(path = cloud_path, showWarnings = FALSE, recursive = TRUE)
  
  refl_dir <- file.path(results_site_path, 
                        'Reflectance',
                        sprintf("res_%d_m", resolution))
  dir.create(path = refl_dir, showWarnings = FALSE, recursive = TRUE)
  refl_path <- file.path(refl_dir, paste(basename(S2obj$S2_Bands$GRANULE), 
                                         '_Refl', sep = ''))
  
  cloudmasks <- generate_cloud_mask(S2obj$S2_Stack, cloud_path, 
                                    S2source, saveRaw)
  write_reflectance(S2obj$S2_Stack, S2obj$S2_Bands, refl_path)
  print("3")
  # Compute spectral indices
  IndexList <- c('NDVI', 'EVI', 'SAVI')
  ReflFactor <- 10000
  HDR_Refl <- read_ENVI_header(get_HDR_name(refl_path))
  SensorBands <- HDR_Refl$wavelength
  Refl <- read_raster(refl_path)
  Indices <- compute_spectral_indices(Refl, SensorBands, IndexList, ReflFactor)
  
  # Save spectral indices
  SI_path <- file.path(results_site_path, 
                       'SpectralIndices',
                       sprintf("res_%d_m", resolution))
  dir.create(path = SI_path,showWarnings = FALSE, recursive = TRUE)
  save_spectral_indices(Indices, SI_path, S2obj$S2_Bands)
  
  # Update Cloud mask based on radiometric filtering
  NDVI_Thresh <- 0.5
  # Cloud_File <- update_cloud_mask(cloud_path, NDVI_Thresh, Indices, cloudmasks)
  # TODO: Fix update_cloud_mask for this func
  
  return(list(results_site_path = results_site_path,
              refl_path = refl_path, 
              HDR_Refl = HDR_Refl
              # Cloud_File = Cloud_File
              ))
}

preprocess_S2_no_dl <- function(path_vector,
                                output_directory,
                                result_path, 
                                resolution = resolution, 
                                S2source,
                                saveRaw = TRUE) {
  process_S2_data <- process_S2_data(output_directory, path_vector, result_path, 
                                     resolution, S2source)
  
  S2obj <- process_S2_data$S2obj
  path_vector <- process_S2_data$path_vector
  
  results_site_path <- file.path(result_path,basename(S2obj$S2_Bands$GRANULE))
  dir.create(path = results_site_path,showWarnings = FALSE,recursive = TRUE)
  
  cloud_path <- file.path(results_site_path,
                          'CloudMask',
                          sprintf("res_%d_m", resolution))
  dir.create(path = cloud_path, showWarnings = FALSE, recursive = TRUE)
  
  refl_dir <- file.path(results_site_path, 
                        'Reflectance',
                        sprintf("res_%d_m", resolution))
  dir.create(path = refl_dir, showWarnings = FALSE, recursive = TRUE)
  refl_path <- file.path(refl_dir, paste(basename(S2obj$S2_Bands$GRANULE), 
                                         '_Refl', sep = ''))
  
  cloudmasks <- generate_cloud_mask(S2obj$S2_Stack, cloud_path, 
                                    S2source, saveRaw)
  write_reflectance(S2obj$S2_Stack, S2obj$S2_Bands, refl_path)
  
  # Compute spectral indices
  IndexList <- c('NDVI', 'EVI', 'SAVI')
  ReflFactor <- 10000
  HDR_Refl <- read_ENVI_header(get_HDR_name(refl_path))
  SensorBands <- HDR_Refl$wavelength
  Refl <- read_raster(refl_path)
  Indices <- compute_spectral_indices(Refl, SensorBands, IndexList, ReflFactor)
  
  # Save spectral indices
  SI_path <- file.path(results_site_path, 
                       'SpectralIndices',
                       sprintf("res_%d_m", resolution))
  dir.create(path = SI_path,showWarnings = FALSE, recursive = TRUE)
  save_spectral_indices(Indices, SI_path, S2obj$S2_Bands)
  
  # Update Cloud mask based on radiometric filtering
  NDVI_Thresh <- 0.5
  Cloud_File <- update_cloud_mask(cloud_path, NDVI_Thresh, Indices, cloudmasks)
  print("coucou")
  return(list(results_site_path = results_site_path,
              refl_path = refl_path, 
              HDR_Refl = HDR_Refl, 
              Cloud_File = Cloud_File))
}


#' Calculates the observation and solar angles.
#'
#' @param sza Solar zenith angle (degrees).
#' @param saa Solar azimuth angle (degrees).
#' @param vza View zenith angle (degrees).
#' @param vaa View azimuth angle (degrees).
#'
#' @return List containing the observation and solar angles.
#'
#' @export
calculate_angles <- function(sza, saa, vza, vaa) {
  # Convert angles to radians
  sza_rad <- sza * pi / 180
  saa_rad <- saa * pi / 180
  vza_rad <- vza * pi / 180
  vaa_rad <- vaa * pi / 180
  
  tto <- vza
  tts <- sza
  
  cos_psi <- -cos(sza_rad) * cos(vza_rad) + sin(sza_rad) * sin(vza_rad) * cos(saa_rad - vaa_rad)
  psi <- acos(cos_psi) * 180 / pi  # Convert back to degrees
  
  psi <- vza - sza
  
  return(list(tto = tto, tts = tts, psi = psi))
}

#' Retrieves Sentinel-2 geometry information.
#'
#' @param refl_path Path to the reflectance file.
#'
#' @return List containing the minimum and maximum values of observation 
#' and solar angles.
#'
#' @export
get_s2_geometry <- function(refl_path) {
  xmlfile <- file.path(dirname(refl_path), 'MTD_TL.xml')
  S2Geom <- get_S2geometry(MTD_TL_xml = xmlfile) 
  # print(length(S2Geom$SZA))
  # print(length(S2Geom$SAA))
  # print(length(S2Geom$VZA))
  # print(length(S2Geom$VAA))
  
  # Go from sza, saa, vza, vaa to tto, tts, psi
  S2Geom <- calculate_angles(S2Geom$SZA, S2Geom$SAA, S2Geom$VZA, S2Geom$VAA)
  
  # Create min and max lists
  GeomAcq <- list()
  GeomAcq$min <- GeomAcq$max <- list()
  GeomAcq$min$tto <- min(S2Geom$tto)
  GeomAcq$max$tto <- max(S2Geom$tto)
  GeomAcq$min$tts <- min(S2Geom$tts)
  GeomAcq$max$tts <- max(S2Geom$tts)
  GeomAcq$min$psi <- min(S2Geom$psi)
  GeomAcq$max$psi <- max(S2Geom$psi)
  # print(GeomAcq)
  return(GeomAcq)
}


#' Retrieves sensor response function.
#'
#' @param HDR_Refl ENVI header information for the reflectance data.
#' @param Path_SensorResponse Path to the sensor response function file.
#'
#' @return Sensor response function.
#'
#' @export
get_sensor_response <- function(HDR_Refl, Path_SensorResponse = NULL) {
  SensorName <- HDR_Refl$`sensor type` # Sentinel_2A
  SRF <- GetRadiometry(SensorName, Path_SensorResponse = Path_SensorResponse)
  return(SRF)
}

#' Adjusts optical constants.
#'
#' @param SRF Sensor response function.
#'
#' @return List containing adjusted optical constants.
#'
#' @export
adjust_optical_constants <- function(SRF) {
  if (is.null(SpecPROSPECT)){
    SpecPROSPECT <- prospect::SpecPROSPECT_FullRange
  }
  if (is.null(SpecSOIL)){
    SpecSOIL <- prosail::SpecSOIL
  }
  if (is.null(SpecPROSPECT)){
    SpecATM <- prosail::SpecATM
  }
  wvl <- SpecPROSPECT$lambda
  SpecSensor <- PrepareSensorSimulation(SpecPROSPECT, SpecSOIL, SpecATM, SRF)
  SpecPROSPECT_Sensor <- SpecSensor$SpecPROSPECT_Sensor
  SpecSOIL_Sensor <- SpecSensor$SpecSOIL_Sensor
  SpecATM_Sensor <- SpecSensor$SpecATM_Sensor
  check_SpectralSampling(SpecPROSPECT, SpecSOIL, SpecATM)
  
  return(list(SpecSensor = SpecSensor,
              SpecPROSPECT_Sensor = SpecPROSPECT_Sensor, 
              SpecSOIL_Sensor = SpecSOIL_Sensor, 
              SpecATM_Sensor = SpecATM_Sensor))
}

#' Defines spectral bands and variables.
#'
#' @param HDR_Refl ENVI header information for the reflectance data.
#' @param S2BandSel List of selected Sentinel-2 bands.
#'
#' @return List containing information about spectral bands and variables.
#'
#' @export
define_spectral_bands_and_variables <- function(HDR_Refl, S2BandSel) {
  ImgBandNames <- strsplit(HDR_Refl$`band names`, split = ',')[[1]]
  Bands2Select <- list()
  for (bpvar in names(S2BandSel)) {
    Bands2Select[[bpvar]] <- match(S2BandSel[[bpvar]], ImgBandNames)
  }
  for (parm in Parms2Estimate){
    S2BandSelect[[parm]] <- c('B3','B4','B5','B6','B7','B8','B11','B12')
    Bands2Select[[parm]] <- match(S2BandSelect[[parm]], SRF$Spectral_Bands)
  }
  return(list(S2BandSel = S2BandSel, 
              ImgBandNames = ImgBandNames,
              Bands2Select = Bands2Select))
}

#' Defines noise levels.
#'
#' @return List containing noise levels.
#'
#' @export
define_noise_levels <- function() {
  NoiseLevel <- list()
  NoiseLevel$lai <- 0.05
  return(NoiseLevel)
}         

#' Defines acquisition geometry.
#'
#' @return List containing minimum and maximum values of observation 
#' and solar angles.
#'
#' @export
define_GeomAcq <- function(){
  GeomAcq <- list()
  GeomAcq$min <- GeomAcq$max <- list()
  GeomAcq$min$tto <- 0
  GeomAcq$max$tto <- 10
  GeomAcq$min$tts <- 20
  GeomAcq$max$tts <- 30
  GeomAcq$min$psi <- 0
  GeomAcq$max$psi <- 360
  return(GeomAcq)
}

#' Defines ATBD (Algorithm Theoretical Basis Document) acquisition geometry.
#'
#' @return List containing minimum and maximum values of observation 
#' and solar angles for ATBD.
#'
#' @export
define_atbd_GeomAcq <- function(){
  GeomAcq <- list()
  GeomAcq$min <- GeomAcq$max <- list()
  GeomAcq$min$tto <- 5
  GeomAcq$max$tto <- 15
  GeomAcq$min$tts <- 25
  GeomAcq$max$tts <- 35
  GeomAcq$min$psi <- 180
  GeomAcq$max$psi <- 210
  return(GeomAcq)
}

#' Extracts summary statistics for each variable in the dataset.
#'
#' @param data Data frame containing the variables.
#' @param filepath Directory path where the statistics file will be saved.
#'
#' @return Data frame containing the summary statistics.
#'
#' @export
extract_summary_stats <- function(data, filepath) {
  # Initialize an empty data frame to store the statistics
  statistics_df <- data.frame(Variable = character(0), 
                              Min = numeric(0), 
                              Max = numeric(0), 
                              Mean = numeric(0), 
                              StdDev = numeric(0))
  
  # Loop through each variable in the data frame and compute summary statistics
  for (variable_name in names(data)) {
    variable <- data[[variable_name]]
    min_value <- min(variable, na.rm = TRUE)
    max_value <- max(variable, na.rm = TRUE)
    mean_value <- mean(variable, na.rm = TRUE)
    std_dev_value <- sd(variable, na.rm = TRUE)
    
    # Append the statistics to the data frame
    statistics_df <- rbind(statistics_df, 
                           data.frame(Variable = variable_name, 
                                      Min = min_value,
                                      Max = max_value,
                                      Mean = mean_value,
                                      StdDev = std_dev_value))
  }
  write.table(statistics_df,
              file = file.path(filepath, "PROSAIL_Stats.txt"), 
              sep = "\t", row.names = FALSE)
  
  return(statistics_df)
}

#' Generates random samples based on specified mean, standard deviation, 
#' and range.
#'
#' @param mean Mean value of the distribution.
#' @param std Standard deviation of the distribution.
#' @param min_val Minimum value of the range.
#' @param max_val Maximum value of the range.
#' @param size Number of samples to generate.
#'
#' @return Vector of generated samples.
#'
#' @export
generate_samples <- function(mean, std, min_val, max_val, size) {
  samples <- rtruncnorm(n = size, 
                        a = min_val, 
                        b = max_val,
                        mean = mean, 
                        sd = std)
  return(samples)
}

#' Modifies the distribution of a parameter in the PROSAIL input based on
#'  the desired distribution.
#'
#' @param InputPROSAIL List containing PROSAIL input parameters.
#' @param parameter_name Name of the parameter to be modified.
#' @param desired_distribution Desired distribution type 
#' ('Gaussian', 'Uniform', 'Unique').
#' @param minval Minimum value for the parameter (optional).
#' @param maxval Maximum value for the parameter (optional).
#' @param mean Mean value for the parameter (optional).
#' @param std Standard deviation for the parameter (optional).
#' @param unique_val Unique value for the parameter (optional).
#'
#' @return Modified PROSAIL input parameters.
#'
#' @export
# modify_parameter_distribution <- function(InputPROSAIL,
#                                           parameter_name,
#                                           desired_distribution,
#                                           minval = NULL,
#                                           maxval = NULL,
#                                           mean = NULL,
#                                           std = NULL,
#                                           unique_val = NULL,
#                                           lidar_lai = NULL) {
#   set.seed(42)
#   if (parameter_name %in% names(InputPROSAIL)
#       || parameter_name == "Cv" || parameter_name == "Zeta"
#       || parameter_name == "fraction_brown" || parameter_name == "diss") {
#     cat("Modify", parameter_name, "\n")
#     if (desired_distribution == 'Gaussian') {
#       cat("Desired distribution", desired_distribution, "\n")
#       mean_val <- ifelse(!is.null(mean), unname(mean), GaussianDistrib$Mean[[parameter_name]])
#       std_val <- ifelse(!is.null(std), unname(std), GaussianDistrib$Std[[parameter_name]])
#       min_val <- ifelse(!is.null(minval), unname(minval), -Inf)
#       max_val <- ifelse(!is.null(maxval), unname(maxval), Inf)
#       
#       new_values <- generate_samples(mean_val, std_val, min_val, max_val, 
#                                      length(InputPROSAIL[[parameter_name]]))
#     } else if (desired_distribution == 'Uniform') {
#       cat("Desired distribution", desired_distribution, "\n")
#       min_val <- ifelse(!is.null(minval), unname(minval), 0)
#       max_val <- ifelse(!is.null(maxval), unname(maxval), 1)
#       
#       new_values <- runif(length(InputPROSAIL[[parameter_name]]), 
#                           min_val, max_val)
#     } else if (desired_distribution == 'Unique') {
#       cat("Desired distribution", desired_distribution, "\n")
#       # Set a unique value for the parameter for all samples
#       unique_val <- ifelse(!is.null(unique_val), unname(unique_val), unique_val[[parameter_name]])
#       new_values <- rep(unique_val, length(InputPROSAIL[[parameter_name]]))
#     } else {
#       warning("Unknown distribution type. No changes made.")
#       new_values <- InputPROSAIL[[parameter_name]]
#     }
#     
#     # Update the specified parameter with new values
#     InputPROSAIL[[parameter_name]] <- new_values
#   } else {
#     warning(paste("Parameter", parameter_name, "not found in InputPROSAIL. 
#                   No changes made. Choose one parameter in this list:", 
#                   str(InputPROSAIL)))
#   }
#   
#   return(InputPROSAIL)
# }
modify_parameter_distribution <- function(InputPROSAIL,
                                          parameter_name,
                                          desired_distribution,
                                          minval = NULL,
                                          maxval = NULL,
                                          mean = NULL,
                                          std = NULL,
                                          unique_val = NULL,
                                          lidar_lai = NULL,
                                          lidar_lai_depth = NULL) {
  set.seed(42)
  
  if (parameter_name %in% names(InputPROSAIL)
      || parameter_name %in% c("Cv", "Zeta", "fraction_brown", "diss")) {
    
    cat("Modify", parameter_name, "\n")
    
    n <- length(InputPROSAIL[[1]])  # number of rows to generate
    
    if (desired_distribution == "Gaussian") {
      cat("Desired distribution: Gaussian\n")
      mean_val <- ifelse(!is.null(mean), unname(mean), GaussianDistrib$Mean[[parameter_name]])
      std_val <- ifelse(!is.null(std), unname(std), GaussianDistrib$Std[[parameter_name]])
      min_val <- ifelse(!is.null(minval), unname(minval), -Inf)
      max_val <- ifelse(!is.null(maxval), unname(maxval), Inf)
      
      new_values <- generate_samples(mean_val, std_val, min_val, max_val, n)
      
    } else if (desired_distribution == "Uniform") {
      cat("Desired distribution: Uniform\n")
      min_val <- ifelse(!is.null(minval), unname(minval), 0)
      max_val <- ifelse(!is.null(maxval), unname(maxval), 1)
      
      new_values <- runif(n, min_val, max_val)
      
    } else if (desired_distribution == "Unique") {
      cat("Desired distribution: Unique\n")
      unique_val <- ifelse(!is.null(unique_val), unname(unique_val), unique_val[[parameter_name]])
      new_values <- rep(unique_val, n)
      
    } else if (desired_distribution == "LiDAR_LAI") {
      cat("Desired distribution: LiDAR_LAI\n")
      if (is.null(lidar_lai)) {
        stop("lidar_lai must be provided for LiDAR_LAI distribution.")
      }
      
      new_values <- sample(lidar_lai, size = n, replace = F)
      
    } else if (desired_distribution == "LiDAR_LAI_Best_Site_Depth") {
      cat("Desired distribution: LiDAR_LAI_Best_Site_Depth\n")
      if (is.null(lidar_lai_depth)) {
        stop("lidar_lai_depth must be provided for LiDAR_LAI distribution.")
      }
      
      new_values <- sample(lidar_lai_depth, size = n, replace = F)
      
    } else {
      warning("Unknown distribution type. No changes made.")
      new_values <- InputPROSAIL[[parameter_name]]
    }
    
    InputPROSAIL[[parameter_name]] <- new_values
    
  } else {
    warning(paste("Parameter", parameter_name, "not found in InputPROSAIL.\n",
                  "No changes made. Choose one parameter in this list:", 
                  paste(names(InputPROSAIL), collapse = ", ")))
  }
  
  return(InputPROSAIL)
}

# ------------------------------------------- Distribution Functions -------------------------------------------------

#' Sets the distribution of PROSAIL parameters based on the specified option.
#'
#' This function modifies the distribution of PROSAIL parameters according to 
#' the specified distribution option. Different distribution options are 
#' available, each representing a different parameterization scheme for the 
#' PROSAIL model.
#'
#' @param distrib_option A character string specifying the distribution option 
#'   to be used. Available options include:
#'   \itemize{
#'     \item{"atbd"}{ATBD (Algorithm Theoretical Basis Document) parameterization.}
#'     \item{"atbd_S2Geom"}{ATBD parameterization with S2 geometry.}
#'     \item{"atbd_JBGeom"}{ATBD parameterization with JB geometry.}
#'     \item{"atbd_fixed_psi"}{ATBD parameterization with fixed PSI.}
#'     \item{"atbd_psi_low"}{ATBD parameterization with low PSI.}
#'     \item{"atbd_psi_high"}{ATBD parameterization with high PSI.}
#'     \item{"atbd_fixed_tto"}{ATBD parameterization with fixed TTO.}
#'     \item{"atbd_tto_low"}{ATBD parameterization with low TTO.}
#'     \item{"atbd_tto_high"}{ATBD parameterization with high TTO.}
#'     \item{"atbd_fixed_tts"}{ATBD parameterization with fixed TTS.}
#'     \item{"atbd_tts_low"}{ATBD parameterization with low TTS.}
#'     \item{"atbd_tts_high"}{ATBD parameterization with high TTS.}
#'     \item{"q_zhang_et_al_2005"}{Parameter distribution based on Zhang et al. (2005).}
#'     \item{"brede_et_al_2020"}{Parameter distribution based on Brede et al. (2020).}
#'     \item{"hauser_et_al_2021"}{Parameter distribution based on Hauser et al. (2021).}
#'     \item{"sinha_et_al_2020"}{Parameter distribution based on Sinha et al. (2020).}
#'     \item{"verhoef_and_bach_2007"}{Parameter distribution based on Verhoef and Bach (2007).}
#'     \item{"shiklomanov_et_al_2016"}{Parameter distribution based on Shiklomanov et al. (2016).}
#'     \item{"atbd_n_high"}{ATBD parameterization with high nitrogen.}
#'     \item{"atbd_n_low"}{ATBD parameterization with low nitrogen.}
#'     \item{"atbd_n_full"}{ATBD parameterization with full range nitrogen.}
#'     \item{"atbd_chl_high"}{ATBD parameterization with high chlorophyll content.}
#'     \item{"atbd_chl_low"}{ATBD parameterization with low chlorophyll content.}
#'     \item{"atbd_chl_full"}{ATBD parameterization with full range chlorophyll content.}
#'     \item{"atbd_brown_high"}{ATBD parameterization with high brown pigment content.}
#'     \item{"atbd_brown_low"}{ATBD parameterization with low brown pigment content.}
#'     \item{"atbd_brown_full"}{ATBD parameterization with full range brown pigment content.}
#'     \item{"atbd_ewt_high"}{ATBD parameterization with high equivalent water thickness.}
#'     \item{"atbd_ewt_low"}{ATBD parameterization with low equivalent water thickness.}
#'     \item{"atbd_ewt_full"}{ATBD parameterization with full range equivalent water thickness.}
#'     \item{"atbd_lma_high"}{ATBD parameterization with high leaf mass per area.}
#'     \item{"atbd_lma_low"}{ATBD parameterization with low leaf mass per area.}
#'     \item{"atbd_lma_full"}{ATBD parameterization with full range leaf mass per area.}
#'     \item{"atbd_lidfa_high"}{ATBD parameterization with high leaf angle distribution factor A.}
#'     \item{"atbd_lidfa_low"}{ATBD parameterization with low leaf angle distribution factor A.}
#'     \item{"atbd_lidfa_full"}{ATBD parameterization with full range leaf angle distribution factor A.}
#'     \item{"atbd_lai_high"}{ATBD parameterization with high leaf area index.}
#'     \item{"atbd_lai_low"}{ATBD parameterization with low leaf area index.}
#'     \item{"atbd_lai_full"}{ATBD parameterization with full range leaf area index.}
#'     \item{"atbd_q_high"}{ATBD parameterization with high diffuse radiation factor.}
#'     \item{"atbd_q_low"}{ATBD parameterization with low diffuse radiation factor.}
#'     \item{"atbd_q_full"}{ATBD parameterization with full range diffuse radiation factor.}
#'     \item{"atbd_psoil_high"}{ATBD parameterization with high soil reflectance.}
#'     \item{"atbd_psoil_low"}{ATBD parameterization with low soil reflectance.}
#'     \item{"atbd_psoil_full"}{ATBD parameterization with full range soil reflectance.}
#'   }
#' @param refl_path Path to the reflectance file used for PROSAIL simulation.
#'
#' @return Modified PROSAIL input parameters.
#'
#' @export
set_prosail_distribution <- function(distrib_option, refl_path,
                                     lidar_lai = NULL,
                                     lidar_lai_depth = NULL) {
  switch(
    distrib_option,
    "q_zhang_et_al_2005" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "Uniform",
                                                    minval = 1,
                                                    maxval = 7.5)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 80)
      InputPROSAIL$CAR <- 0.25*InputPROSAIL$CHL # need to modify CAR too if CHL is
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Uniform",
                                                    minval = 1,
                                                    maxval = 4.5)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "EWT",
                                                    "Uniform",
                                                    minval = 0.001,
                                                    maxval = 0.15)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LMA",
                                                    "Uniform",
                                                    minval = 0.001,
                                                    maxval = 0.04)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Uniform",
                                                    minval = 0.00001,
                                                    maxval = 8)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "TypeLidf",
                                                    "Unique",
                                                    unique_val = 2)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Uniform",
                                                    minval = 10,
                                                    maxval = 89)
    },
    "brede_et_al_2020" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 8)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 80)
      InputPROSAIL$CAR <- 0.25*InputPROSAIL$CHL # need to modify CAR too if CHL is
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Uniform",
                                                    minval = 1,
                                                    maxval = 2.5)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "EWT",
                                                    "Uniform",
                                                    minval = 0.002,
                                                    maxval = 0.025)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LMA",
                                                    "Uniform",
                                                    minval = 0.001,
                                                    maxval = 0.025)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 1)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "tto",
                                                    "Unique",
                                                    unique_val = 0)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "tts",
                                                    "Uniform",
                                                    minval = 27.5,
                                                    maxval = 80)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "psi",
                                                    "Unique",
                                                    unique_val = 0)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "q",
                                                    "Unique",
                                                    unique_val = 0)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "psoil",
                                                    "Unique",
                                                    unique_val = 0)
    },
    "hauser_et_al_2021" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "Gaussian",
                                                    mean = 2,
                                                    std = 1,
                                                    minval = 0.01,
                                                    maxval = 3.5)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Gaussian",
                                                    mean = 30,
                                                    std = 20,
                                                    minval = 10,
                                                    maxval = 60)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CAR",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 15)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "ANT",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 10)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Uniform",
                                                    minval = 1.4,
                                                    maxval = 1.7)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "EWT",
                                                    "Uniform",
                                                    minval = 0.001,
                                                    maxval = 0.045)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LMA",
                                                    "Uniform",
                                                    minval = 0.001,
                                                    maxval = 0.040)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Unique",
                                                    unique_val = 0.01)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "TypeLidf",
                                                    "Unique",
                                                    unique_val = 2)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Uniform",
                                                    minval = 30,
                                                    maxval = 70)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "q",
                                                    "Unique",
                                                    unique_val = 0.01)
    },
    "sinha_et_al_2020" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "Gaussian",
                                                    mean = 1.5,
                                                    std = 0.7,
                                                    minval = 0,
                                                    maxval = 2.8)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Gaussian",
                                                    mean = 30,
                                                    std = 20,
                                                    minval = 0,
                                                    maxval = 60)
      InputPROSAIL$CAR <- 0.25*InputPROSAIL$CHL # need to modify CAR too if CHL is
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Uniform",
                                                    minval = 1.5,
                                                    maxval = 2.5)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "EWT",
                                                    "Uniform",
                                                    minval = 0.030,
                                                    maxval = 0.050)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LMA",
                                                    "Uniform",
                                                    minval = 0.012,
                                                    maxval = 0.030)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Uniform",
                                                    minval = 20,
                                                    maxval = 50)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "q",
                                                    "Unique",
                                                    unique_val = 0.01)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "tts",
                                                    "Unique",
                                                    unique_val = 47)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "psoil",
                                                    "Uniform",
                                                    minval = 0.2,
                                                    maxval = 1)
    },
    "verhoef_and_bach_2007" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- -0.3
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Unique",
                                                    unique_val = 2)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Unique",
                                                    unique_val = 60)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "EWT",
                                                    "Unique",
                                                    unique_val = 0.009)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LMA",
                                                    "Unique",
                                                    unique_val = 0.005)
      InputPROSAIL$CAR <- 0.25*InputPROSAIL$CHL # need to modify CAR too if CHL is
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Unique",
                                                    unique_val = -0.2)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "q",
                                                    "Unique",
                                                    unique_val = 0.05)
    },
    "shiklomanov_et_al_2016" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- -0.3
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Uniform",
                                                    minval = 1.09,
                                                    maxval = 3)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Uniform",
                                                    minval = 0.78,
                                                    maxval = 106.72)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CAR",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 25.3)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "EWT",
                                                    "Uniform",
                                                    minval = 0.0043,
                                                    maxval = 0.0439)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LMA",
                                                    "Uniform",
                                                    minval = 0.0017,
                                                    maxval = 0.0152)
    },
    "atbd" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE,
                                       GeomAcq = define_atbd_GeomAcq()
                                       )
      InputPROSAIL$LIDFb <- 0
    },
    "atbd_codist_false" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE,
                                       Codist_LAI = FALSE,
                                       GeomAcq = define_atbd_GeomAcq()
      )
      InputPROSAIL$LIDFb <- 0
    },
    "atbd_clumping_uniform" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE,
                                       GeomAcq = define_atbd_GeomAcq())
      InputPROSAIL$Cv <- 0
      InputPROSAIL$Zeta <- 0
      InputPROSAIL$fraction_brown <- 0
      InputPROSAIL$diss <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "Cv",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 1)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "Zeta",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 1)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "fraction_brown",
      #                                               "Unique",
      #                                               unique_val = 0.5)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "diss",
      #                                               "Unique",
      #                                               unique_val = 0.5)
      InputPROSAIL$LIDFb <- 0
    },
    "atbd_clumping_gaussian" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE,
                                       GeomAcq = define_atbd_GeomAcq())
      InputPROSAIL$Cv <- 0
      InputPROSAIL$Zeta <- 0
      InputPROSAIL$fraction_brown <- 0
      InputPROSAIL$diss <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "Cv",
                                                    "Gaussian",
                                                    minval = 0,
                                                    maxval = 1,
                                                    mean = 0.5,
                                                    std = 0.1)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "Zeta",
                                                    "Gaussian",
                                                    minval = 0,
                                                    maxval = 1,
                                                    mean = 0.5,
                                                    std = 0.1)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "fraction_brown",
      #                                               "Unique",
      #                                               unique_val = 0.5)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "diss",
      #                                               "Unique",
      #                                               unique_val = 0.5)
      InputPROSAIL$LIDFb <- 0
    },
    "atbd_JBGeom" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE,
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
    },
    "atbd_n_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Uniform",
                                                    minval = 1,
                                                    maxval = 3)
    },
    "atbd_n_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Uniform",
                                                    minval = 1,
                                                    maxval = 1.2)
    },
    "atbd_n_full" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Uniform",
                                                    minval = 1,
                                                    maxval = 5)
    },
    "atbd_chl_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Uniform",
                                                    minval = 40,
                                                    maxval = 100)
      InputPROSAIL$CAR <- 0.25*InputPROSAIL$CHL # need to modify CAR too if CHL is
    },
    "atbd_chl_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 60)
      InputPROSAIL$CAR <- 0.25*InputPROSAIL$CHL # need to modify CAR too if CHL is
    },
    "atbd_chl_full" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 100)
      InputPROSAIL$CAR <- 0.25*InputPROSAIL$CHL # need to modify CAR too if CHL is
    },
    "atbd_brown_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Uniform",
                                                    minval = 0.00001,
                                                    maxval = 3)
    },
    "atbd_brown_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Uniform",
                                                    minval = 0.01,
                                                    maxval = 1)
    },
    "atbd_brown_full" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Uniform",
                                                    minval = 0.00001,
                                                    maxval = 8)
    },
    "atbd_ewt_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "EWT",
                                                    "Uniform",
                                                    minval = 0.02,
                                                    maxval = 0.05)
    },
    "atbd_ewt_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "EWT",
                                                    "Uniform",
                                                    minval = 0.001,
                                                    maxval = 0.02)
    },
    "atbd_ewt_full" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "EWT",
                                                    "Uniform",
                                                    minval = 0.001,
                                                    maxval = 0.05)
    },
    "atbd_lma_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LMA",
                                                    "Uniform",
                                                    minval = 0.02,
                                                    maxval = 0.04)
    },
    "atbd_lma_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LMA",
                                                    "Uniform",
                                                    minval = 0.001,
                                                    maxval = 0.02)
    },
    "atbd_lma_full" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LMA",
                                                    "Uniform",
                                                    minval = 0.001,
                                                    maxval = 0.04)
    },
    "atbd_lidfa_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL$TypeLidf <- 2
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Uniform",
                                                    minval = 50,
                                                    maxval = 90)
    },
    "atbd_lidfa_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL$TypeLidf <- 2
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Uniform",
                                                    minval = 10,
                                                    maxval = 50)
    },
    "atbd_lidfa_full" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL$TypeLidf <- 2
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Uniform",
                                                    minval = 10,
                                                    maxval = 90)
    },
    "atbd_lai_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "Uniform",
                                                    minval = 1,
                                                    maxval = 7)
    },
    "atbd_lai_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 3)
    },
    "atbd_lai_full" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 10)
    },
    "atbd_q_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "q",
                                                    "Uniform",
                                                    minval = 0.25,
                                                    maxval = 0.5)
    },
    "atbd_q_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "q",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 0.25)
    },
    "atbd_q_full" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "q",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 0.5)
    },
    "atbd_psoil_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "psoil",
                                                    "Uniform",
                                                    minval = 0.5,
                                                    maxval = 1)
    },
    "atbd_psoil_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "psoil",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 0.5)
    },
    "atbd_psoil_full" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "psoil",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 1)
    },
    "atbd_S2Geom" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = get_s2_geometry(refl_path))
      InputPROSAIL$LIDFb <- 0
    },
    "atbd_fixed_tto" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL$tto <- 5
    },
    "atbd_tto_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "tto",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 5)
    },
    "atbd_tto_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "tto",
                                                    "Uniform",
                                                    minval = 5,
                                                    maxval = 10)
    },
    "atbd_fixed_tts" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL$tts <- 28
    },
    "atbd_tts_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "tts",
                                                    "Uniform",
                                                    minval = 20,
                                                    maxval = 25)
    },
    "atbd_tts_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "tts",
                                                    "Uniform",
                                                    minval = 25,
                                                    maxval = 30)
    },
    "atbd_fixed_psi" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL$psi <- 150
    },
    "atbd_psi_low" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "psi",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 180)
    },
    "atbd_psi_high" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "psi",
                                                    "Uniform",
                                                    minval = 180,
                                                    maxval = 360)
    },
    "atbd_tts_surprising" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = define_GeomAcq())
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "tts",
                                                    "Uniform",
                                                    minval = 27.5,
                                                    maxval = 80)
    },
    "atbd_optim_6m_try" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = get_s2_geometry(refl_path),
                                       Codist_LAI = F)
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Uniform",
                                                    minval = 30,
                                                    maxval = 50)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "Gaussian",
                                                    minval = 0.1,
                                                    maxval = 12,
                                                    mean = 5,
                                                    std = 3)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "LMA",
      #                                               "Uniform",
      #                                               minval = 27.5,
      #                                               maxval = 80)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Uniform",
                                                    minval = 0.01,
                                                    maxval = 0.02)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "N",
      #                                               "Uniform",
      #                                               minval = 27.5,
      #                                               maxval = 80)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 80)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "psoil",
                                                    "Uniform",
                                                    minval = 0.2,
                                                    maxval = 1)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "q",
      #                                               "Uniform",
      #                                               minval = 27.5,
      #                                               maxval = 80)
    },
    "atbd_optim_common" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       # GeomAcq = get_s2_geometry(refl_path),
                                       Codist_LAI = F)
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Gaussian",
                                                    minval = 30,
                                                    maxval = 50,
                                                    mean = 40,
                                                    std = 20)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 0)
    },
    "atbd_optim_common_brownmodif" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       # GeomAcq = get_s2_geometry(refl_path),
                                       Codist_LAI = F)
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Gaussian",
                                                    minval = 30,
                                                    maxval = 50,
                                                    mean = 40,
                                                    std = 20)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Gaussian",
                                                    minval = 0,
                                                    maxval = 1,
                                                    mean = 0,
                                                    std = 0.1)
    },
    "optim_Aigoual" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = get_s2_geometry(refl_path),
                                       Codist_LAI = F)
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Gaussian",
                                                    minval = 30,
                                                    maxval = 50,
                                                    mean = 40,
                                                    std = 20)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "Gaussian",
                                                    minval = 0,
                                                    maxval = 11,
                                                    mean = 5,
                                                    std = 3)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "LMA",
      #                                               "Uniform",
      #                                               minval = 27.5,
      #                                               maxval = 80)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Gaussian",
                                                    minval = 0,
                                                    maxval = 1,
                                                    mean = 0,
                                                    std = 0.1)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "N",
      #                                               "Uniform",
      #                                               minval = 1,
      #                                               maxval = 2)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "CHL",
      #                                               "Uniform",
      #                                               minval = 0,
      #                                               maxval = 80)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "psoil",
      #                                               "Uniform",
      #                                               minval = 0.2,
      #                                               maxval = 1)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "q",
      #                                               "Uniform",
      #                                               minval = 0.05,
      #                                               maxval = 0.6)
    },
    "optim_Blois" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       # GeomAcq = get_s2_geometry(refl_path),
                                       Codist_LAI = F)
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Gaussian",
                                                    minval = 20,
                                                    maxval = 50,
                                                    mean = 35,
                                                    std = 20)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "LiDAR_LAI_Best_Site_Depth",
                                                    lidar_lai_depth = lidar_lai_depth)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 0)
    },
    "optim_Blois_brownmodif" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       # GeomAcq = get_s2_geometry(refl_path),
                                       Codist_LAI = F)
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Gaussian",
                                                    minval = 20,
                                                    maxval = 50,
                                                    mean = 35,
                                                    std = 20)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "LiDAR_LAI_Best_Site_Depth",
                                                    lidar_lai_depth = lidar_lai_depth)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Gaussian",
                                                    minval = 0,
                                                    maxval = 1,
                                                    mean = 0,
                                                    std = 0.1)
    },
    "optim_Mormal" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = get_s2_geometry(refl_path),
                                       Codist_LAI = F)
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Gaussian",
                                                    minval = 30,
                                                    maxval = 50,
                                                    mean = 40,
                                                    std = 20)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "lai",
      #                                               "Gaussian",
      #                                               minval = 0.1,
      #                                               maxval = 11,
      #                                               mean = 5,
      #                                               std = 3)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "LiDAR_LAI_Best_Site_Depth",
                                                    lidar_lai_depth = lidar_lai_depth)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "LMA",
      #                                               "Uniform",
      #                                               minval = 27.5,
      #                                               maxval = 80)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 0)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "N",
      #                                               "Uniform",
      #                                               minval = 1,
      #                                               maxval = 2)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "CHL",
      #                                               "Uniform",
      #                                               minval = 0,
      #                                               maxval = 80)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "psoil",
      #                                               "Uniform",
      #                                               minval = 0.2,
      #                                               maxval = 1)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "q",
      #                                               "Uniform",
      #                                               minval = 0.05,
      #                                               maxval = 0.3)
    },
    "optim_lidarD_Aigoual" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = get_s2_geometry(refl_path),
                                       Codist_LAI = F)
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Gaussian",
                                                    minval = 30,
                                                    maxval = 60,
                                                    mean = 45,
                                                    std = 20)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "LiDAR_LAI",
                                                    lidar_lai = lidar_lai)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "LMA",
      #                                               "Uniform",
      #                                               minval = 27.5,
      #                                               maxval = 80)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Uniform",
                                                    minval = 0.01,
                                                    maxval = 0.02)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "N",
      #                                               "Uniform",
      #                                               minval = 1,
      #                                               maxval = 2)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 80)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "psoil",
                                                    "Uniform",
                                                    minval = 0.2,
                                                    maxval = 1)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "q",
      #                                               "Uniform",
      #                                               minval = 0.05,
      #                                               maxval = 0.6)
    },
    "optim_lidarD_Blois" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = get_s2_geometry(refl_path),
                                       Codist_LAI = F)
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Gaussian",
                                                    minval = 30,
                                                    maxval = 50,
                                                    mean = 40,
                                                    std = 20)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "LiDAR_LAI",
                                                    lidar_lai = lidar_lai)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "LMA",
      #                                               "Uniform",
      #                                               minval = 27.5,
      #                                               maxval = 80)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "BROWN",
                                                    "Uniform",
                                                    minval = 0.01,
                                                    maxval = 0.02)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "N",
      #                                               "Uniform",
      #                                               minval = 1,
      #                                               maxval = 2)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "CHL",
                                                    "Uniform",
                                                    minval = 0,
                                                    maxval = 80)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "psoil",
      #                                               "Uniform",
      #                                               minval = 0.2,
      #                                               maxval = 1)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "q",
      #                                               "Uniform",
      #                                               minval = 0.05,
      #                                               maxval = 0.3)
    },
    "optim_lidarD_Mormal" = {
      InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
                                       GeomAcq = get_s2_geometry(refl_path),
                                       Codist_LAI = F)
      InputPROSAIL$LIDFb <- 0
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "LIDFa",
                                                    "Gaussian",
                                                    minval = 30,
                                                    maxval = 60,
                                                    mean = 45,
                                                    std = 20)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "lai",
                                                    "LiDAR_LAI",
                                                    lidar_lai = lidar_lai)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "LMA",
      #                                               "Uniform",
      #                                               minval = 27.5,
      #                                               maxval = 80)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "BROWN",
      #                                               "Uniform",
      #                                               minval = 0.01,
      #                                               maxval = 0.02)
      InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
                                                    "N",
                                                    "Uniform",
                                                    minval = 1,
                                                    maxval = 2)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "CHL",
      #                                               "Uniform",
      #                                               minval = 0,
      #                                               maxval = 80)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "psoil",
      #                                               "Uniform",
      #                                               minval = 0.2,
      #                                               maxval = 1)
      # InputPROSAIL <- modify_parameter_distribution(InputPROSAIL,
      #                                               "q",
      #                                               "Uniform",
      #                                               minval = 0.05,
      #                                               maxval = 0.3)
    },
    {
      stop(paste("Error:", 
                 distrib_option,
                 "is an invalid distribution option.",
                 "Please refer to the function code to see the possibilities"))
    }
  )
  
  return(InputPROSAIL)
}
