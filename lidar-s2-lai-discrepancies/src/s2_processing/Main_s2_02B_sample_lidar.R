# clean environment
rm(list = ls(all=TRUE)); gc()
if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getSourceEditorContext()$path))
lapply(X = as.list(list.files(path = '../libraries', full.names = T)), FUN = source, verbose = F)
library(terra)
library(sgsR)

# 1- Directories & input data
set.seed(42)
input_dir <- '../../01_DATA'
output_dir <- '../../03_RESULTS'
sites <- c('Aigoual', 'Blois', 'Mormal')
# sites <- c('Blois', 'Mormal')
# sites <- c('Aigoual')

output_dir_samples <- file.path(output_dir, sites, 'PROSAIL_Optimization', 'sampling')
names(output_dir_samples) <- sites
nbSamplesS2Refl <- 5000
# Samplingmethods <- c('random', 'stratified', 'stratified_uniform')
Samplingmethods <- c('stratified_uniform')
lapply(X = output_dir_samples, FUN = dir.create, showWarnings = F, recursive = T)

# Loop through each site and sampling method
for (site in sites) {
  cat("Processing site:", site, "\n")
  
  # LiDAR LAI file path
  lidar_lai_path <- file.path(input_dir, site, 'LiDAR', 'lidarlai_res_10_m.tif')
  lidar_lai_path <- file.path(input_dir, site, 'LiDAR/PAD_Profiles_Classic',
                              'ladstack.tif')
  
  # PAD directories
  dsm_dir <- file.path(input_dir, site, 'LiDAR/testPADs', 'PAD_Profiles_DSM')
  dsm_Above20_dir <- file.path(input_dir, site, 'LiDAR/testPADs', 'PAD_Profiles_DSM_Above20')
  dsm_keepTrees_dir <- file.path(input_dir, site, 'LiDAR/testPADs', 'PAD_Profiles_DSM_keepTrees')
  dsm_keepTreesAbove20_dir <- file.path(input_dir, site, 'LiDAR/testPADs', 'PAD_Profiles_DSM_keepTreesAbove20')
  dtm_dir <- file.path(input_dir, site, 'LiDAR/testPADs', 'PAD_Profiles_DTM')
  dtm_keepTrees_dir <- file.path(input_dir, site, 'LiDAR/testPADs', 'PAD_Profiles_DTM_keepTrees')
  
  # Loop through sampling methods
  for (i in seq_along(Samplingmethods)) {
    Samplingmethod <- Samplingmethods[i]
    
    # Open sample location
    filename_sampling <- file.path(output_dir_samples[site], 
                                   paste0('Sampling_', Samplingmethod, 
                                          "_nbSamples_", nbSamplesS2Refl,'.GPKG'))
    sample_locations <- terra::vect(filename_sampling)
    samples_id <- sample_locations$sample_id
    
    # Open sample location
    # filename_sampling_above20 <- file.path(output_dir_samples[site], 
    #                                paste0('Sampling_', Samplingmethod, 
    #                                       "_Above20_nbSamples_5000.GPKG"))
    # sample_locations_above20 <- terra::vect(filename_sampling_above20)
    
    # 1. Sampling LiDAR LAI
    lidar_values <- unlist(terra::extract(sum(terra::rast(lidar_lai_path), 
                                              na.rm = T), sample_locations)[-1])
    # 
    # # Save LiDAR LAI values for sample points to CSV
    filename_lidar_lai <- file.path(output_dir_samples[site],
                                    paste0('LiDAR_LAI_Samples_', Samplingmethod,
                                           "_nbSamples_", nbSamplesS2Refl, '.csv'))
    readr::write_delim(x = as.data.frame(cbind(lidar_values, samples_id)),
                       file = filename_lidar_lai,
                       delim = '\t')

    # 2. Sampling PAD_Profiles_DTM files
    dtm_files <- list.files(dtm_dir, pattern = "PAD_.*\\.tif$", full.names = TRUE)

    for (dtm_file in dtm_files) {
      # Extract X.5 value from file name, assuming file format is PAD_X.5_40.tif
      x_value <- floor(as.numeric(sub(".*_(\\d+\\.5)_40.tif", "\\1", basename(dtm_file))))
      depth <- 40 - x_value # Retrieve effective depth

      # Sample and save from DTM
      dtm_values <- unlist(terra::extract(terra::rast(dtm_file), sample_locations)[-1])
      lidar_values <- dtm_values

      filename_dtm <- file.path(output_dir_samples[site],
                                paste0('PAD_DTM_Depth_', depth,
                                       '_Samples_', Samplingmethod,
                                       "_nbSamples_", nbSamplesS2Refl, '.csv'))
      readr::write_delim(x = as.data.frame(cbind(lidar_values, samples_id)),
                         file = filename_dtm,
                         delim = '\t')
    }

    # 3. Sampling PAD_Profiles_DSM files
    dsm_files <- list.files(dsm_dir, pattern = "PAD_.*\\.tif$", full.names = TRUE)

    for (dsm_file in dsm_files) {
      # Extract X.5 value from file name, assuming file format is PAD_X.5_40.tif
      x_value <- floor(as.numeric(sub(".*_(\\d+\\.5)_40.tif", "\\1", basename(dsm_file))))
      depth <- 40 - x_value # Retrieve effective depth

      # Sample and save from DSM
      dsm_values <- unlist(terra::extract(terra::rast(dsm_file), sample_locations)[-1])
      lidar_values <- dsm_values

      filename_dsm <- file.path(output_dir_samples[site],
                                paste0('PAD_DSM_Depth_', depth,
                                       '_Samples_', Samplingmethod,
                                       "_nbSamples_", nbSamplesS2Refl, '.csv'))
      readr::write_delim(x = as.data.frame(cbind(lidar_values, samples_id)),
                         file = filename_dsm,
                         delim = '\t')
    }
    
    # 4. Sampling PAD_Profiles_DSM_keepTrees files
    dsm_keepTrees_files <- list.files(dsm_keepTrees_dir, pattern = "PAD_.*\\.tif$", full.names = TRUE)
    
    for (dsm_keepTrees_file in dsm_keepTrees_files) {
      # Extract X.5 value from file name, assuming file format is PAD_X.5_40.tif
      x_value <- floor(as.numeric(sub(".*_(\\d+\\.5)_40.tif", "\\1", basename(dsm_keepTrees_file))))
      depth <- 40 - x_value # Retrieve effective depth
      
      # Sample and save from DSM
      dsm_keepTrees_values <- unlist(terra::extract(terra::rast(dsm_keepTrees_file), sample_locations)[-1])
      lidar_values <- dsm_keepTrees_values
      filename_dsm_keepTrees <- file.path(output_dir_samples[site], 
                                          paste0('PAD_DSM_keepTrees_Depth_', depth, 
                                                 '_Samples_', Samplingmethod, 
                                                 "_nbSamples_", nbSamplesS2Refl, '.csv'))
      readr::write_delim(x = as.data.frame(cbind(lidar_values, samples_id)),  
                         file = filename_dsm_keepTrees, 
                         delim = '\t')
    }
    
    # 5. Sampling PAD_Profiles_DTM_keepTrees files
    dtm_keepTrees_files <- list.files(dtm_keepTrees_dir, pattern = "PAD_.*\\.tif$", full.names = TRUE)
    
    for (dtm_keepTrees_file in dtm_keepTrees_files) {
      # Extract X.5 value from file name, assuming file format is PAD_X.5_40.tif
      x_value <- floor(as.numeric(sub(".*_(\\d+\\.5)_40.tif", "\\1", basename(dtm_keepTrees_file))))
      depth <- 40 - x_value # Retrieve effective depth
      
      # Sample and save from DTM
      dtm_keepTrees_values <- unlist(terra::extract(terra::rast(dtm_keepTrees_file), sample_locations)[-1])
      lidar_values <- dtm_keepTrees_values
      filename_dtm_keepTrees <- file.path(output_dir_samples[site], 
                                          paste0('PAD_DTM_keepTrees_Depth_', depth, 
                                                 '_Samples_', Samplingmethod, 
                                                 "_nbSamples_", nbSamplesS2Refl, '.csv'))
      readr::write_delim(x = as.data.frame(cbind(lidar_values, samples_id)),  
                         file = filename_dtm_keepTrees, 
                         delim = '\t')
    }
    
    # Actualize samples_id
    # samples_id <- sample_locations_above20$sample_id
    # 
    # # 6. Sampling PAD_Profiles_DSM_Above20 files
    # dsm_Above20_files <- list.files(dsm_Above20_dir, pattern = "PAD_.*\\.tif$", full.names = TRUE)
    # 
    # for (dsm_Above20_file in dsm_Above20_files) {
    #   # Extract X.5 value from file name, assuming file format is PAD_X.5_40.tif
    #   x_value <- floor(as.numeric(sub(".*_(\\d+\\.5)_40.tif", "\\1", basename(dsm_Above20_file))))
    #   depth <- 40 - x_value # Retrieve effective depth
    #   
    #   # Sample and save from DSM
    #   dsm_Above20_values <- unlist(terra::extract(terra::rast(dsm_Above20_file), sample_locations_above20)[-1])
    #   lidar_values <- dsm_Above20_values
    #   filename_dsm_Above20 <- file.path(output_dir_samples[site], 
    #                                       paste0('PAD_DSM_Above20_Depth_', depth, 
    #                                              '_Samples_', Samplingmethod, 
    #                                              "_nbSamples_", nbSamplesS2Refl, '.csv'))
    #   readr::write_delim(x = as.data.frame(cbind(lidar_values, samples_id)), 
    #                      file = filename_dsm_Above20, 
    #                      delim = '\t')
    # }
    # 
    # # 7. Sampling PAD_Profiles_DSM_keepTreesAbove20 files
    # dsm_keepTreesAbove20_files <- list.files(dsm_keepTreesAbove20_dir, pattern = "PAD_.*\\.tif$", full.names = TRUE)
    # 
    # for (dsm_keepTreesAbove20_file in dsm_keepTreesAbove20_files) {
    #   # Extract X.5 value from file name, assuming file format is PAD_X.5_40.tif
    #   x_value <- floor(as.numeric(sub(".*_(\\d+\\.5)_40.tif", "\\1", basename(dsm_keepTreesAbove20_file))))
    #   depth <- 40 - x_value # Retrieve effective depth
    #   
    #   # Sample and save from DTM
    #   dsm_keepTreesAbove20_values <- unlist(terra::extract(terra::rast(dsm_keepTreesAbove20_file), sample_locations_above20)[-1])
    #   lidar_values <- dsm_keepTreesAbove20_values
    #   filename_dsm_keepTreesAbove20 <- file.path(output_dir_samples[site], 
    #                             paste0('PAD_DSM_keepTreesAbove20_Depth_', depth, 
    #                                    '_Samples_', Samplingmethod, 
    #                                    "_nbSamples_", nbSamplesS2Refl, '.csv'))
    #   readr::write_delim(x = as.data.frame(cbind(lidar_values, samples_id)), 
    #                      file = filename_dsm_keepTreesAbove20, 
    #                      delim = '\t')
    # }
  }
}