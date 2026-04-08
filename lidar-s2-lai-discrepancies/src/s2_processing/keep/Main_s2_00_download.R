# clean environment
rm(list = ls(all=TRUE)); gc()
if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(preprocS2)

# 1- Directories & input data
input_dir <- '../../01_DATA'
output_dir <- '../../03_RESULTS'
sites <- c('Aigoual', 'Blois', 'Mormal')

# define resolution & source of data
resolution <- 10
S2source <- 'SAFE'

for (site in sites){
  date <- switch(site,
                 'Aigoual' = '2021-07-11', 
                 'Blois' = '2021-06-14', 
                 'Mormal' = '2021-06-14')
  safe <- switch(site,
                'Aigoual' = 'S2A_MSIL2A_20210711T104031_N0301_R008_T31TEJ_20210711T134956.SAFE', 
                'Blois' = 'S2A_MSIL2A_20210614T105031_N0300_R051_T31TCN_20210614T140120.SAFE', 
                'Mormal' = 'S2A_MSIL2A_20210614T105031_N0300_R051_T31UER_20210614T140120.SAFE')
  safe_path <- file.path(input_dir, site, "Sentinel-2", date, safe)
  path_vector <- file.path(input_dir, site, "Geo_Files", "utm_init.shp")
  output_subdir <- file.path(output_dir, site, "PROSAIL_Optimization")
  if (!dir.exists(output_subdir)) 
    dir.create(dirname(output_subdir), recursive = TRUE, showWarnings = F)
  
  S2obj <- extract_from_S2_L2A(Path_dir_S2 = safe_path,
                               path_vector = path_vector,
                               S2source = S2source,
                               resolution = resolution)
  
  # create specific result directory corresponding to granule name
  results_site_path <- file.path(output_subdir, basename(S2obj$S2_Bands$GRANULE))
  dir.create(path = results_site_path,showWarnings = FALSE,recursive = TRUE)
  ##____________________________________________________________________##
  ##                        Write CLOUD MASK                            ##
  ##--------------------------------------------------------------------##
  # directory for cloud mask
  Cloud_path <- file.path(results_site_path,'CloudMask')
  dir.create(path = Cloud_path,showWarnings = FALSE,recursive = TRUE)
  # Filename for cloud mask
  cloudmasks <- save_cloud_s2(S2_stars = S2obj$S2_Stack,
                              Cloud_path = Cloud_path,
                              S2source = S2source, SaveRaw = T)
  ##____________________________________________________________________##
  ##                        Write REFLECTANCE                           ##
  ##--------------------------------------------------------------------##
  # directory for Reflectance
  Refl_dir <- file.path(results_site_path,'Reflectance')
  dir.create(path = Refl_dir,showWarnings = FALSE,recursive = TRUE)
  # filename for Reflectance
  Refl_path <- file.path(Refl_dir,paste(basename(S2obj$S2_Bands$GRANULE),'_Refl',sep = ''))
  
  # Save Reflectance file as ENVI image with BIL interleaves
  # metadata files are important to account for offset applied on S2 L2A products 
  tile_S2 <- preprocS2::get_tile(S2obj$S2_Bands$GRANULE)
  dateAcq_S2 <- preprocS2::get_dateAcq(S2product = dirname(dirname(S2obj$S2_Bands$GRANULE)))
  save_reflectance_s2(S2_stars = S2obj$S2_Stack, 
                      Refl_path = Refl_path,
                      S2Sat = NULL, 
                      tile_S2 = tile_S2, 
                      dateAcq_S2 = dateAcq_S2,
                      Format = 'ENVI', 
                      datatype = 'Int16', 
                      MTD = S2obj$S2_Bands$metadata, 
                      MTD_MSI = S2obj$S2_Bands$metadata_MSI)
}