# ---
# title: "1.calculate_lad_profiles.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-03-12"
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

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

# Import useful functions
source("functions_lidar.R")
source("../functions_plots.R")

# LiDR optimization

# opt_chunk_size(ctg) <- 0 # Processing by files
# opt_chunk_buffer(ctg) <- 10

set_lidr_threads(0)
get_lidr_threads()
# memory.limit(4000000000)

# Directories 
data <- "../../01_DATA"
site <- "Blois"
data_site <- file.path(data, site)
output_dir <- "../../03_RESULTS/Blois/LiDAR/PAI/lidR/"

# LAS
las_norm_utm <- "LiDAR/3-las_normalized_utm/"
las_norm_utm_dir <- file.path(data_site, las_norm_utm)
las_norm_utm_files <- list.files(las_norm_utm_dir, pattern = "\\.las$", 
                                 full.names = TRUE)

# DTM 
# dtm_path <- file.path(las_norm_utm_dir, "dtm/res_1m/rasterize_terrain.vrt")
# dtm_charged <- terra::rast(dtm_path)

# Initialization
z0 <- 2
# las_norm_utm_files_subset <- las_norm_utm_files[69:71]

# Catalog
cat("lasDpath: ", las_norm_utm_dir, "\n")
ctg <- readLAScatalog(las_norm_utm_files)
plot(ctg)

# Test
# las_test <- readLAS(las_norm_utm_files[72])
# lad_test <- LAD(las_test@data$Z, z0=z0)



# Calculations

# Calculate the maximum height (zmax) across all LAS files

t1 <- Sys.time()

# zmax <- max(sapply(las_norm_utm_files, 
#                    function(file) max(readLAS(file)$Z)))
# zmax <- round(zmax * 2) / 2

zmax <- max(ctg@data$Max.Z)

# Preallocate vectors
metrics_list <- vector("list", length(las_norm_utm_files))
names_list <- character(length(las_norm_utm_files))

# Parallel computation using future package
# plan(multisession, workers = 2)

for (i in 1: length(las_norm_utm_files)) {
  
  las_i <- readLAS(las_norm_utm_files[i])
  
  # Metrics function
  metrics_i <- calculate_lad_profiles(las_i, z0, zmax)
  
  # Export : save results in a list
  metrics_list[[i]] <- metrics_i
  names_list[i] = paste0("las_", i)
  
  cat(paste0("File: ", i, "\n"))
  flush.console()
}

# Parallelized computation of LAD profiles
# metrics_list <- future_lapply(las_norm_utm_files, function(file, z0, zmax) {
#   las_i <- readLAS(file)
#   metrics_i <- calculate_lad_profiles(las_i, z0, zmax)
#   # cat(paste0("File: ", i, "\n"))
#   # flush.console()
#   print(1)
#   return(metrics_i)
# }, z0 = z0, zmax = zmax)

names(metrics_list) <- names_list

for (i in seq_along(metrics_list)) {
  new_z <- seq(z0+0.5, zmax)
  missing_z <- setdiff(new_z, metrics_list[[i]]$PAD_Profile$z)
  missing_data <- data.frame(z = missing_z, lad = rep(0, length(missing_z)))
  
  metrics_list[[i]]$PAD_Profile <- rbind(
    metrics_list[[i]]$PAD_Profile, missing_data)
  
  metrics_list[[i]]$PAD_Profile <- metrics_list[[i]]$PAD_Profile[order(
    metrics_list[[i]]$PAD_Profile$z), ]
  
  metrics_list[[i]]$PAD_Profile$lad[is.na(
    metrics_list[[i]]$PAD_Profile$lad)] <- 0
  
  row.names(metrics_list[[i]]$PAD_Profile) <- NULL
}




# Initialize a list to store the aggregated LAD sums for varying heights
lad_sums <- list()
zmin <- z0+0.5
zmax <- round(zmax * 2) / 2
zvec <- seq(zmin, zmax, by = 1)
k <- 1

# Iterate over the height range and calculate the sums
for (i in zmin:zmax) {
  # Subset the metrics_list to get LAD values for heights up to i
  subset_metrics <- lapply(metrics_list, function(metrics) {
    subset_data <- subset(metrics$PAD_Profile, z >= i)
    return(subset_data)
  })
  
  # Calculate the sum of LAD values for the current height
  sum_lad <- numeric(length(subset_metrics[[1]]$lad))
  for (j in seq_along(subset_metrics)) {
    sum_lad <- sum_lad + subset_metrics[[j]]$lad
  }
  
  # Store the height (z) along with the sum in the list
  lad_sums[[k]] <- list(z = zvec[1:length(sum_lad)], 
                        sum_lad = sum_lad / 463)
  
  # Remove the first element from zvec for the next iteration
  zvec <- zvec[-1]
  
  k <- k + 1
}

# Create an empty plot
# plot(NA, 
#      NA,
#      xlim = c(zmin, zmax),
#      ylim = c(0, max(unlist(lapply(lad_sums, function(x) max(x$sum_lad))))), 
#      xlab = "Height (z)", 
#      ylab = "Sum of LAD values",
#      main = "Sum of LAD values for varying heights")
# 
# # Iterate over lad_sums list and plot each curve
# for (i in seq_along(lad_sums)) {
#   lines(lad_sums[[i]]$z, lad_sums[[i]]$sum_lad, 
#         xlab = "Height (z)", ylab = "Sum of LAD values", 
#         col = i)
# }
# 
# # Add legend
# legend("topright", 
#        legend = as.character(seq_along(lad_sums)), 
#        col = seq_along(lad_sums), 
#        lty = 1, 
#        title = "Curve", 
#        cex = 0.8)

lai_list <- list()
for (i in seq_along(lad_sums)) {
  lai_sum <- sum(lad_sums[[i]]$sum_lad, 
                 na.rm = TRUE)
  lai_list[[i]] <- list(z = lad_sums[[i]]$z, 
                        sum_lad = lai_sum)
}

# # Initialize a vector to store the aggregated LAD values
# sum_lad <- numeric(length(metrics_list[[i]]$PAD_Profile$z))
# 
# # Iterate over the vectors in lad_list and aggregate LAD values for each z value
# for (i in seq_along(metrics_list)) {
#   sum_lad <- sum_lad + metrics_list[[i]]$PAD_Profile$lad
# }
# 
# # Plot the sum of LAD values against z
# plot(sum_lad, 
#      metrics_list[[i]]$PAD_Profile$z, 
#      type = "l",
#      xlab = "Sum of LAD",
#      ylab = "z")


# Convert the list to a data frame
las_df <- flatten_and_merge(metrics_list, length(new_z))

# Define column names
col_names <- c("ID_LAS", 
               # "Rumple_Index",
               "Hmax", 
               "VCI", 
               "CV_PAD", 
               "PAI", 
               "Max_PAD", 
               "H_maxPAD", 
               paste0("PAD_", min(new_z):max(new_z)))

# Set column names
colnames(las_df) <- col_names

# Write the data frame to a CSV file
write.csv(las_df, 
          file = file.path(output_dir, paste0(site, "_las_data.csv")),
          row.names = FALSE)

t2 <- Sys.time()
time <- t2-t1
time








#Calculate the maximum height (zmax) across all LAS files
las_norm_utm_files_subset <- las_norm_utm_files[70:71]
zmax <- max(sapply(las_norm_utm_files_subset, 
                   function(file) max(readLAS(file)$Z)))
zmax <- round(zmax * 2) / 2

# Iterate over each LAS file and calculate LAD profiles for each height stratum
lad_profiles <- list()
max_height_less_than_desired <- FALSE
file_number <- 1

for (file in las_norm_utm_files) {
  las <- readLAS(file)
  las$Z <- las$Z + zmax - max(las$Z)
  test_lad_data <- list()
  test_lad_data$z <- 0
  
  if (max(las$Z) < zmax) {
    max_height_less_than_desired <- TRUE
  }
  
  for (z_start in seq(zmax - 1, by = -1)) {
    z_end <- zmax
    if (z_start >= z0) {
      # Filter points within the desired height range
      las_subset <- filter_poi(las, Z >= z_start & Z <= z_end)
      # Calculate LAD profile for the current height stratum
      lad_data <- LAD(las_subset$Z)
      # If identical, exit
      if (length(lad_data$z) == length(test_lad_data$z)) {
        break  # Exit the inner loop if lad_data is equal to test_lad_data
      }
      else {
        test_lad_data <- lad_data
      }
      
      # Handle cases where the height range exceeds the actual range of heights
      if (length(lad_data$lad) == 0) {
        lad_data <- data.frame(z = seq(z_start, z_end), 
                               lad = rep(0, z_end - z_start + 1))
      }
      # Store the LAD profile
      lad_profiles[[paste0("las_", file_number)]] <- list(
        z = lad_data$z, 
        lad = lad_data$lad
      )
      # lad_profiles[[paste(z_start, "_", z_end)]] <- sum(lad_data$lad, 
      #                                                   na.rm = TRUE)
    }
  }
  file_number <- file_number + 1
}

if (max_height_less_than_desired) {
  for (z_start in seq(2.5, by = 1)) {  # Start from the minimum height
    z_end <- zmax
    if (z_start >= z0) {
      # Create a LAD profile with all zero values
      lad_data <- data.frame(z = seq(z_start, z_end), 
                             lad = rep(0, z_end - z_start + 1))
      # Store the LAD profile
      lad_profiles[[paste(z_start, "_", z_end)]] <- lad_data
    }
  }
}


# Loop

list_metrics <- list()
names_list<- NULL

files <- las_norm_utm_files #_subset

t1 <- Sys.time()

zmax <- max(sapply(files, 
                   function(file) max(readLAS(file)$Z)))
# zmax <- round(zmax * 2) / 2 

# for (j in 1:10) # for all buffers (10 to 100m)
# {
for ( i in 1: length(files)) # For all tiles
{
  las_i <- readLAS(files[i])
  
  # Metrics function
  metrics_i <- calculate_lidar_metrics(las_i, z0, zmax)
  
  # Export : save results in a list
  list_metrics[[i]] <- metrics_i
  names_list[i] = paste0("las_", i)
  
}

# Final file setup
names(list_metrics) <- names_list

# Execution time
t2 <- Sys.time()
time <- t2-t1
print(paste0("Execution time = ", round(time, 2), " mins"))

# Convert list into DF

# PAD Profiles

# Search for lmax of profil_pad 

lad <- lapply(list_metrics, "[", "Profil_pad" )

l_max <- 0
index <- 0
vec_z <- 0
for (i in 1: length(lad)) { 
  
  l <- length(lad[i][[1]]$Profil_pad$z)
  if ( l > l_max) {
    index <- i
    vec_z <- lad[i][[1]]$Profil_pad$z
  }
  l_max <- max(l, l_max)
}

# Prepare new list

names_pad <- paste0("Pad_", vec_z)
names_pad_norm <- paste0("Pad_norm_", vec_z)

metrics_names=c("Rumple_index",
                "VCI", 
                "Hmax",
                "CV_PAD",
                "PAI", 
                "Max_PAD", 
                "H_maxPAD",
                # "Gap_fraction_raster", 
                # "Gap_fraction_pts",
                names_pad,
                names_pad_norm)    

# Metrics Table

Nb_col <- length(metrics_names)
Nb_lig <- length(list_metrics)

metrics_tab <- data.frame(matrix(NA,ncol= Nb_col,nrow=Nb_lig))

colnames(metrics_tab) <- metrics_names
rownames(metrics_tab)<- names(list_metrics)

for (i in 1: length(list_metrics)){
  
  vect <-NULL
  vec_non_pad <-NULL
  vec_z_pad <- NULL
  
  vect <- unlist(list_metrics[[i]])
  vec_non_pad <- vect[1:7]
  
  metrics_tab[i, 1:7] <- vec_non_pad
  
  
  vec_z_pad <- vect[- c(1:7)]
  # print(length(vec_z_pad))
  metrics_tab[i, 7:(7+(length(vec_z_pad)/2))] <- vec_z_pad[(length(
    vec_z_pad)/2 +1):length(vec_z_pad)]
  
}

# Add ID Plot as first column
metrics_tab$id_plot <- row.names(metrics_tab) 
metrics_tab <- metrics_tab[c(dim(metrics_tab)[2], 1:(dim(metrics_tab)[2]-1))]
head(metrics_tab)

write.table(metrics_tab,
            paste(output_dir,
                  "metrics",
                  "_",
                  site,
                  ".csv"), 
            row.names = FALSE, 
            col.names = TRUE, 
            dec=".", 
            sep = ";")

# }







# Calculations
results_norm_utm <- calculate_lad_profiles(las_norm_utm_files, z0)
lad_profiles_norm_utm <- results_norm_utm$lad
cv_lad_values_norm_utm <- results_norm_utm$cv_lad
lai_values_norm_utm <- results_norm_utm$lai

# for (i in 1:length(las_files)){
#   las_i <- readLAS(las_files[i])
#   lad_i <- LAD(las_i@data$Z, z0=z0)
#   lad_list[[i]] <- lad_i
#   CV_LAD_i <- 100 * sd(lad_i$lad, na.rm = TRUE) / mean(lad_i$lad, na.rm = TRUE)    
#   PAI_i <- sum(lad_i$lad, na.rm = TRUE)
# }
# print(lad_list)

# Initialize a vector to store the aggregated LAD values
sum_lad <- numeric(length(lad_list_modif[[1]]$z))

# Iterate over the vectors in lad_list and aggregate LAD values for each z value
for (i in seq_along(lad_list_modif)) {
  sum_lad <- sum_lad + lad_list_modif[[i]]$lad
}

cv_lad <- cv(sum_lad)
print(cv_lad)

# Plot the sum of LAD values against z
plot(sum_lad, 
     lad_list_modif[[1]]$z, 
     type = "l",
     xlab = "z",
     ylab = "Sum of LAD")

# Catalog
cat("lasDpath: ", las_dir, "\n")
ctg <- readLAScatalog(las_norm_utm_dir)
plot(ctg)
# plot(readLAS(las_files[100]))

cv <- function(x) {
  sd(x) / mean(x)
}

# Test on a small catalog
n = 5
ctg_2 = head(ctg, n = 100)
ctg_3 <- tail(ctg_2, n = 5)

# Function to calculate LAD for each file in the catalog
ladcv <- pixel_metrics(ctg, ~cv(LAD(Z)$lad), res = 10)
plot(ladcv)

# Export raster to tif
ladcv_file <- "ladcv.tif"
ladcv_path <- file.path(output_dir, ladcv_file)
writeRaster(ladcv, filename = ladcv_path, overwrite = TRUE)

# Calculate LAD
lad <- pixel_metrics(ctg_3, ~lad_custom(ctg_3), res = 10)
writeRaster(lad, filename = file.path(output_dir, "lad_10m.tif"), 
            overwrite = TRUE)

# Define custom function to compute CV
cv <- function(x) sd(x) / mean(x)

# Calculate CV of LAD
ladcv <- grid_metrics(lad, ~cv(LAD.Z.), res = 10)
writeRaster(lad, filename = file.path(output_dir, "cv_lad_10m.tif"), 
            overwrite = TRUE)

# Plot the result
plot(ladcv, main = "Coefficient of Variation (CV) of LAD")


calculate_lad_profiles <- function(las_i, z0, zmax){
  
  las_i <- filter_poi(las_i, Z >= z0)
  Hmax_i <- max(las_i$Z)
  if (Hmax_i < z0) {
    # Set all values to 0 or NA or NULL
    rumple_i <- NA
    VCI_i <- NA
    CV_PAD_i <- NA
    # Hmax_i <- NA
    PAI_i <- NA
    max_PAD_i  <- NA
    H_max_PAD_i <- NA
    pad_i_translated <- NULL
    pad_i_translated$z <- 0
    pad_i_translated$lad <- NA
  } else {
    
    # CHM
    # chm_i <- grid_canopy(las_i, res = 1,
    #                      pitfree(thresholds = c(0, 2, 5, 10, 15),
    #                              max_edge = c(0, 1),
    #                              subcircle = 0.5))
    
    # DTM
    # dtm_i <- rasterize_terrain(las_i, res = 1, algorithm = tin())
    
    # DTS
    # dts_i <- rasterize_canopy(las_i, res = 1, p2r())
    
    # CHM
    # chm_i <- dts_i # - dtm_i
    
    # Calculations before translation
    # rumple_i <- rumple_index(chm_i)
    VCI_i <- VCI(las_i$Z, zmax = Hmax_i)
    
    # Translation
    las_i$Z <- las_i$Z + zmax - max(las_i$Z)
    pad_i_translated <- LAD(las_i@data$Z, z0=0)
    
    # Calculations after translation
    CV_PAD_i <- 100 * sd(pad_i_translated$lad, 
                         na.rm = TRUE) / mean(pad_i_translated$lad, 
                                              na.rm = TRUE)
    
    PAI_i <- sum(pad_i_translated$lad, na.rm = TRUE)
    max_PAD_i <- max(pad_i_translated$lad, na.rm = TRUE)
    H_max_PAD_i <- pad_i_translated$z[which(pad_i_translated$lad==max_PAD_i)]
  }
  
  # Output
  metrics = list(
    # rumple_i,
    Hmax_i,
    VCI_i,
    CV_PAD_i, 
    PAI_i, 
    max_PAD_i, 
    H_max_PAD_i, 
    pad_i_translated)
  
  names(metrics)=c(
    # "Rumple_Index",
    "Hmax",
    "VCI",
    "CV_PAD",
    "PAI",
    "Max_PAD",
    "H_maxPAD",
    "PAD_Profile")
  
  return(metrics)
}

flatten_and_merge <- function(metrics_list, nb_z) {
  flattened_list <- lapply(seq_along(metrics_list), function(i) {
    las <- metrics_list[[i]]
    pad_profile <- as.data.frame(las$PAD_Profile)
    
    # Flatten the PAD_Profile data frame before combining
    pad_values <- unlist(pad_profile$lad)
    
    # Ensure that the number of columns in pad_profile matches nb_z
    if (length(pad_values) < nb_z) {
      pad_values <- c(pad_values, rep(0, nb_z - length(pad_values)))
    } else if (length(pad_values) > nb_z) {
      pad_values <- pad_values[1:nb_z]  # Trim excess values
    }
    
    las_num <- paste0("las_", i)
    
    return(c(las_num, 
             as.numeric(unlist(las[c(
               # "Rumple_Index", 
               "Hmax", 
               "VCI", 
               "CV_PAD", 
               "PAI", 
               "Max_PAD", 
               "H_maxPAD")], 
               use.names = FALSE)), 
             pad_values))
  })
  
  return(do.call(rbind, flattened_list))
}

calculate_lidar_metrics <- function(las_i, z0, zmax){
  
  # Rumple
  chm_i <- grid_canopy(las_i, res = 1, 
                       pitfree(thresholds = c(0, 2, 5, 10, 15), 
                               max_edge = c(0, 1), 
                               subcircle = 0.5))
  
  rumple_i <- rumple_index(chm_i)
  
  # Z
  las_i$Z <- las_i$Z + zmax - max(las_i$Z)
  
  lad_profiles <- list()
  for (z_start in seq(zmax - 1, by = -1)) {
    z_end <- zmax
    if (z_start >= z0) {
      # Calculate LAD profile for the current height stratum
      pad_i_norm <- LAD(las_i$Z)
      
      # Calculate the total LAD within the interval
      total_lad <- sum(pad_i_norm$lad[
        pad_i_norm$z >= z_end & pad_i_norm$z <= zmax])
      
      # pad_i_norm_sum <- data.frame(z = paste0(z_start, "-", z_end), 
      #                          lad = sum(pad_i_norm$lad, na.rm = TRUE))
      # pad_i_norm_sum <- sum(pad_i_norm$lad, na.rm = TRUE)
      
      # Store the LAD profile
      lad_profiles[[paste0(z_start, "_", z_end)]] <- pad_i_norm$lad
    }
  }
  
  las_subset <- filter_poi(las_i, Z >= z0)
  if (length(las_subset$Z) == 0) {   
    Hmax_i <- max(las_i$Z)
    VCI_i <-NA
    
    pad_i <-NULL
    pad_i$z <-0
    pad_i$lad <- NA 
    
    pad_i_norm <- data.frame(z = seq(z_start, z_end), 
                             lad = rep(0, z_end - z_start + 1)) 
    
    CV_PAD_i <- NA
    PAI_i <- NA
    max_PAD_i <- NA
    H_max_PAD_i <- NA
  } 
  else {
    
    Hmax_i <- max(las_i$Z)
    VCI_i <- VCI(las_subset$Z, by = 1, zmax= Hmax_i)
    
    step <- 1
    pad_i <-LAD(las_i$Z, dz=step, k=0.5, z0=z0)
    
    if (dim(pad_i)[1] == 0) {
      
      pad_i <- NULL
      nb_ind_z <- ceiling(Hmax_i - (z0 + 0.5*step )) 
      
      pad_i <- data.frame(matrix(NA,ncol= 2, nrow=nb_ind_z)) 
      colnames(pad_i) = c("z", "lad")
      pad_i$z <- seq(from = (z0+0.5*step),
                     to = (z0+ (nb_indice_z - 0.5)*step), by =step) 
      pad_i$lad <- rep(NA, nb_ind_z)
      
      CV_PAD_i <- NA
      PAI_i <- NA
      max_PAD_i <- NA
      H_max_PAD_i <- NA
      
    } else {
      
      CV_PAD_i <- 100 * sd(pad_i$lad, 
                           na.rm = TRUE) / mean(pad_i$lad, 
                                                na.rm = TRUE) 
      PAI_i <- sum(pad_i$lad, na.rm = TRUE)  
      
      max_PAD_i <- max(pad_i$lad, na.rm = TRUE)
      H_max_PAD_i <- pad_i$z[which(pad_i$lad==max_PAD_i)]
      
      # Norm
      CV_PAD_norm_i <- 100 * sd(pad_i$lad, 
                                na.rm = TRUE) / mean(pad_i$lad, 
                                                     na.rm = TRUE) 
      PAI_norm_i <- sum(pad_i$lad, na.rm = TRUE)  
      
      max_PAD_norm_i <- max(pad_i$lad, na.rm = TRUE)
      H_max_PAD_norm_i <- pad_i$z[which(pad_i$lad==max_PAD_i)]
      
    }
  }
  
  # # gap fraction raster ave seuil de hauteur pour definir les trouees
  # # raster moins lisse que pour rumple, restitue mieux les trouees 
  # chm_tr_i <- grid_canopy(las_i_h, res=0.5, p2r() )    
  # 
  # gapf_raster_i <- 100 * (length(
  #   chm_tr_i@data@values[!is.na(chm_tr_i@data@values) 
  #                        & chm_tr_i@data@values < seuil_gap ]))
  # / ( length( chm_tr_i@data@values[!is.na(chm_tr_i@data@values)] ) )
  # 
  # # gap fraction issu du ratio des premiers retours parvenant sous le seuil 
  # # des hauteurs par rapport au nbre total de premier retours
  # 
  # gapf_pts_i <- 100 * dim(filter_poi(
  #   las_i_h, ReturnNumber ==1L & Z <seuil_gap)@data)[1] /  
  #   dim(filter_first(las_i_h)@data)[1]
  
  
  metrics = list(rumple_i,
                 VCI_i, 
                 Hmax_i,
                 CV_PAD_i, 
                 PAI_i, 
                 max_PAD_i, 
                 H_max_PAD_i, 
                 # gapf_raster_i, 
                 # gapf_pts_i,
                 pad_i,
                 pad_i_norm_sum
  )
  names(metrics)=c("Rumple_index",
                   "VCI",
                   "Hmax", 
                   "CV_PAD",
                   "PAI",
                   "Max_PAD",
                   "H_maxPAD",
                   # "Gap_fraction_raster",
                   # "Gap_fraction_pts",
                   "Profil_pad",
                   "Profil_pad_norm_from_top")
  
  return(metrics)
  
}