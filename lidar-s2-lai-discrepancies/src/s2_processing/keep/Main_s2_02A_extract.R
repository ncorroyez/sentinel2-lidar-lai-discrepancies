# clean environment
rm(list = ls(all=TRUE)); gc()
if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getSourceEditorContext()$path))
lapply(X = as.list(list.files(path = '../libraries', full.names = T)), FUN = source, verbose = F)
library(terra)
library(sgsR)
library(dplyr)

# 1- Directories & input data
input_dir <- '../../01_DATA'
output_dir <- '../../03_RESULTS'
sites <- c('Aigoual', 'Blois', 'Mormal')
# sites <- c('Blois', 'Mormal')
# sites <- c('Aigoual')

output_dir_samples <- file.path(output_dir, sites, 'PROSAIL_Optimization', 'sampling')
names(output_dir_samples) <- sites
datesAcq <- list('Aigoual' = '2021-07-11',
                 'Blois' = '2021-06-14',
                 'Mormal' = '2021-06-14')
name_vect <- 'utm_init.shp'

# 2- define number of samples to extract
nbSamples <- 5000
Samplingmethods <- c('random', 'stratified', 'stratified_uniform')
Samplingmethods <- c('stratified_uniform')
lapply(X = output_dir_samples, FUN = dir.create, showWarnings = F, recursive = T)

# 3- Sample points for each site and method
for (site in sites) {
  l2a <- switch(site,
                "Aigoual" = "L2A_T31TEJ_A031608_20210711T104217",
                "Blois" = "L2A_T31TCN_A031222_20210614T105443",
                "Mormal" = "L2A_T31UER_A031222_20210614T105443")
  
  # Define paths for LiDAR metrics, AOI, and corresponding S2 reflectance data
  aoi_path <- file.path(input_dir, site, 'Geo_Files', name_vect)
  S2_path <- file.path(output_dir, site, 'PROSAIL_Optimization', l2a, 
                       'Reflectance', paste0(l2a, "_Refl"))
  sampling_path <- file.path(input_dir, site, "Geo_Files", "aoi_sampling_prep.gpkg")
  field_points_path <- file.path(input_dir, site, "Geo_Files", "data_utm31n.geojson")
  
  # lidar_lai_path <- file.path(input_dir, site, 'LiDAR', 'lidarlai_res_10_m.tif')
  lidar_lai_path <- file.path(input_dir, site, 'LiDAR/PAD_Profiles_Classic',
                              'ladstack.tif')
  # ladstack <- sum(rast(lidar_lai_path), na.rm = T)
  
  mean_path <- file.path(input_dir, site, 'LiDAR', 'mean_res_10_m.tif')
  max_path <- file.path(input_dir, site, 'LiDAR', 'max_res_10_m.tif')
  lskew_path <- file.path(input_dir, site, 'LiDAR', 'lskew_res_10_m.tif')
  
  lidar_lai_rast <- rast(lidar_lai_path)
  max_rast <- rast(max_path)
  
  # Ensure alignment of LiDAR and S-2 and remove NA
  aoi_sampling_path <- align_and_remove_na_for_aoi(aoi_path = aoi_path, 
                                                   S2_path = S2_path, 
                                                   lidar_lai_path = lidar_lai_path,
                                                   sampling_path = sampling_path,
                                                   save = T)
  
  # Loop through sampling methods
  for (Samplingmethod in Samplingmethods) {
    set.seed(42)
    S2Refl <- get_s2_samples(aoi_path = aoi_sampling_path, 
                             S2_path = S2_path,
                             mean_path = mean_path,
                             lskew_path = lskew_path,
                             lidar_lai_path = lidar_lai_path,
                             max_path = max_path,
                             field_points_path = field_points_path,
                             nbSamples = nbSamples, 
                             site = site,
                             method = Samplingmethod)
    # stop()
    # Save sample location
    filename_sampling <- file.path(output_dir_samples[site], 
                                   paste0('Sampling_', Samplingmethod,
                                          "_nbSamples_", nbSamples, '.GPKG'))
    sample_locations <- S2Refl$sample_location
    sample_locations$sample_id <- seq_len(nrow(sample_locations))
    terra::writeVector(x = sample_locations, 
                       filename = filename_sampling, 
                       filetype = 'GPKG', 
                       overwrite = TRUE)
    
    # Save Sentinel-2 reflectance corresponding to samples
    filename_refl <- file.path(output_dir_samples[site], 
                               paste0('S2_reflectance_', Samplingmethod,
                                      "_nbSamples_", nbSamples, '.csv'))
    s2_refl <- S2Refl$S2_refl
    readr::write_delim(x = s2_refl, 
                       file = filename_refl, 
                       delim = '\t')
  }
}
