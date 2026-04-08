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
sites <- c('Blois', 'Mormal')
# sites <- c('Aigoual')

output_dir_samples <- file.path(output_dir, sites, 'PROSAIL_Optimization', 'sampling')
names(output_dir_samples) <- sites
datesAcq <- list('Aigoual' = '2021-07-11',
                 'Blois' = '2021-06-14',
                 'Mormal' = '2021-06-14')
name_vect <- 'utm_init.shp'

# 2- define number of samples to extract
nbSamples <- 10000
Samplingmethods <- c('random', 'stratified', 'stratified_uniform')
# Samplingmethods <- c('stratified_uniform')
lapply(X = output_dir_samples, FUN = dir.create, showWarnings = F, recursive = T)

# 3- Sample points for each site and method
for (site in sites) {
  l2a <- switch(site,
                "Aigoual" = "L2A_T31TEJ_A031608_20210711T104217",
                "Blois" = "L2A_T31TCN_A031222_20210614T105443",
                "Mormal" = "L2A_T31UER_A031222_20210614T105443")
  # Define paths for LiDAR metrics, AOI, and corresponding S2 reflectance data
  aoi_path <- file.path(input_dir, site, 'Geo_Files', name_vect)
  # S2_path <- file.path(output_dir, site, 'PROSAIL_Optimization', 'raster_samples', 
  #                      paste0(site, '_', datesAcq[[site]], '.tiff'))
  S2_path <- file.path(output_dir, site, 'PROSAIL_Optimization', l2a, 
                       'Reflectance', paste0(l2a, "_Refl"))
  lidar_lai_path <- file.path(input_dir, site, 'LiDAR', 'lidarlai_res_10_m.tif')
  mean_path <- file.path(input_dir, site, 'LiDAR', 'mean_res_10_m.tif')
  lskew_path <- file.path(input_dir, site, 'LiDAR', 'lskew_res_10_m.tif')
  sampling_path <- file.path(input_dir, site, "Geo_Files", "aoi_sampling_prep.gpkg")
  field_points_path <- file.path(input_dir, site, "Geo_Files", "data_utm31n.geojson")
  
  # Ensure alignment of LiDAR and S-2 and remove NA
  aoi_sampling_path <- align_and_remove_na_for_aoi(aoi_path = aoi_path, 
                                                   S2_path = S2_path, 
                                                   lidar_lai_path = lidar_lai_path,
                                                   sampling_path = sampling_path,
                                                   save = T)
  
  # Loop through sampling methods
  for (Samplingmethod in Samplingmethods) {
    # Sample for S2 using the current method
    S2Refl <- get_s2_samples(aoi_path = aoi_sampling_path, 
                             S2_path = S2_path,
                             mean_path = mean_path,
                             lskew_path = lskew_path,
                             lidar_lai_path = lidar_lai_path,
                             field_points_path = field_points_path,
                             nbSamples = nbSamples, 
                             site = site,
                             method = Samplingmethod)
    
    # Save sample location
    filename_sampling <- file.path(output_dir_samples[site], 
                                   paste0('Sampling_', Samplingmethod, '.GPKG'))
    terra::writeVector(x = S2Refl$sample_location, 
                       filename = filename_sampling, 
                       filetype = 'GPKG', 
                       overwrite = TRUE)
    
    # Save Sentinel-2 reflectance corresponding to samples
    filename_refl <- file.path(output_dir_samples[site], 
                               paste0('S2_reflectance_', Samplingmethod, '.csv'))
    readr::write_delim(x = S2Refl$S2_refl, 
                       file = filename_refl, 
                       delim = '\t')
    
    # Apply S2 sampling to LiDAR
    sample_locations <- terra::vect(filename_sampling)
    lidar_values <- unlist(terra::extract(terra::rast(lidar_lai_path),
                                          sample_locations)[-1])
    hist(lidar_values)
    
    # Save LiDAR LAI values for sample points to CSV
    filename_lidar_lai <- file.path(output_dir_samples[site], 
                                    paste0('LiDAR_LAI_Samples_', Samplingmethod, '.csv'))
    readr::write_delim(x = as.data.frame(lidar_values), 
                       file = filename_lidar_lai, 
                       delim = '\t')
  }
}
