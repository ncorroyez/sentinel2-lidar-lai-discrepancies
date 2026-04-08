# clean environment
rm(list = ls(all=TRUE)); gc()
if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getSourceEditorContext()$path))
lapply(X = as.list(list.files(path = '../libraries', full.names = T)), FUN = source, verbose = F)
library(terra)
library(dplyr)
library(ggplot2)
library(fmsb)
library(tidyr)
library(rlang)
library(data.table)
library(stringr)

# 1- Directories & input data
output_dir      <- '../../03_RESULTS'
sites           <- c('Aigoual', 'Blois', 'Mormal')
name_strategy   <- "LIDFa_lai_LMA_BROWN_N_CHL_psoil_q_Final"
parms2test      <- c("LIDFa", "lai", "LMA", "BROWN", "N", "CHL", "psoil", "q")
nbSamplesS2Refl <- 5000

sampling_methods <- c('stratified_uniform')
norm_methods     <- c('DSM_keepTrees')
# norm_methods     <- c('DSM_keepTrees', 'DTM_keepTrees')
# norm_methods     <- c(
#                       'DSM_keepTrees',
#                       'DTM_keepTrees', 'DSM', 'DTM')
# norm_methods <- c("DSM")
depths           <- 1:38

# Function to compute metrics given lidar & S2 vectors
compute_metrics <- function(lidar_vals, s2_vals) {
  
  valid <- complete.cases(lidar_vals, s2_vals)
  lidar_vals <- lidar_vals[valid]
  s2_vals <- s2_vals[valid]
  
  if (length(lidar_vals) < 2) {
    return(rep(NA, 7))
  }
  
  lm_mod <- lm(s2_vals ~ lidar_vals)
  r       <- cor(s2_vals, lidar_vals, use = "complete.obs")
  r2      <- summary(lm_mod)$r.squared
  rmse    <- sqrt(mean((s2_vals - lidar_vals)^2, na.rm = TRUE))
  nrmse   <- rmse / IQR(lidar_vals, na.rm = TRUE)
  slope   <- coef(lm_mod)[2]
  intercept <- coef(lm_mod)[1]
  bias    <- mean(s2_vals - lidar_vals, na.rm = TRUE)
  mae     <- mean(abs(s2_vals - lidar_vals), na.rm = TRUE)
  names_out <- c("R", "R2", "RMSE", "NRMSE", "Slope", "Intercept", "Bias", "MAE")
  return(setNames(c(r, r2, rmse, nrmse, slope, intercept, bias, mae), names_out))
}

# Initialize combined results container
combined_results <- data.table()

# ------------------------------------------------------------------------------
# PART 1: per-site processing
# ------------------------------------------------------------------------------
for (site in sites) {
  message("Processing site: ", site)
  
  for (norm in norm_methods) {
    for (depth in depths) {
      print(depth)
      for (sampling_method in sampling_methods) {
        
        # Build file paths
        lidar_csv <- file.path(output_dir, site, 'PROSAIL_Optimization',
                               'sampling',
                               paste0("PAD_", norm, "_Depth_", depth,
                                      "_Samples_", sampling_method,
                                      "_nbSamples_", nbSamplesS2Refl, ".csv"))
        s2_csv    <- file.path(output_dir, site, 'PROSAIL_Optimization',
                               'PROSAIL_Models', name_strategy,
                               paste0("LAI_estimated_", sampling_method,
                                      "_nbSamples_", nbSamplesS2Refl, ".csv"))
        # Skip if missing
        if (!file.exists(lidar_csv) || !file.exists(s2_csv)) next
        
        # Load data
        lidar_dt <- fread(lidar_csv, header = TRUE, sep = "\t")
        s2_dt <- fread(s2_csv, sep = "\t")
        # s2_dt$samples_id <- 1:nrow(s2_dt)
        # s2_dt <- s2_dt[samples_id %in% lidar_dt$samples_id]
        # s2_dt <- setdiff(colnames(s2_dt), "samples_id")
        
        if (nrow(lidar_dt) != nrow(s2_dt)) {
          stop("Row mismatch for site ", site,
               ", norm=", norm, ", depth=", depth)
        }
        
        # Compute metrics for all columns in s2_dt at once
        metric_list <- lapply(colnames(s2_dt), function(colname) {
          vals_s2 <- s2_dt[[colname]]
          m <- compute_metrics(lidar_dt$lidar_values, vals_s2)
          
          atbd_flag <- str_count(colname, "\\.1(\\D|$)") == length(parms2test)
          
          data.frame(
            Site      = site,
            Norm      = norm,
            Depth     = depth,
            Method    = sampling_method,
            Column    = colname,
            R         = round(m["R"], 3),
            R2        = round(m["R2"], 3),
            RMSE      = round(m["RMSE"], 3),
            NRMSE     = round(m["NRMSE"], 3),
            Slope     = round(m["Slope"], 3),
            Intercept = round(m["Intercept"], 3),
            Bias      = round(m["Bias"], 3),
            MAE       = round(m["MAE"], 3),
            ATBD      = atbd_flag,
            stringsAsFactors = FALSE
          )
        })
        
        # Bind all metrics together at once
        combined_results <- rbindlist(list(combined_results, rbindlist(metric_list)))
        # }
        
      }  # sampling_method
    }    # depth
  }      # norm
}        # site
first_combined <- combined_results

out_csv <- file.path(output_dir,
                     paste0("all_results_combined_", name_strategy, "13_05.csv"))
if (!file.exists(out_csv)) {
  write.csv(combined_results, out_csv, row.names = FALSE)
} else {
  write.table(combined_results, out_csv, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
}
# ------------------------------------------------------------------------------
# PART 2: “All sites” aggregated processing
# ------------------------------------------------------------------------------
message("Processing ALL SITES aggregated case")

for (norm in norm_methods) {
  for (depth in depths) {
    print(depth)
    for (sampling_method in sampling_methods) {
      
      # Initialize empty data.tables for aggregation
      agg_lidar <- data.table()
      agg_s2    <- data.table()
      
      # Aggregate each site's data using a loop
      for (site in sites) {
        lidar_csv <- file.path(output_dir, site, 'PROSAIL_Optimization',
                               'sampling',
                               paste0("PAD_", norm, "_Depth_", depth,
                                      "_Samples_", sampling_method,
                                      "_nbSamples_", nbSamplesS2Refl, ".csv"))
        s2_csv    <- file.path(output_dir, site, 'PROSAIL_Optimization',
                               'PROSAIL_Models', name_strategy,
                               paste0("LAI_estimated_", sampling_method,
                                      "_nbSamples_", nbSamplesS2Refl, ".csv"))
        if (!file.exists(lidar_csv) || !file.exists(s2_csv)) next
        
        # Load data
        lidar_dt <- fread(lidar_csv, header = TRUE, sep = "\t")
        s2_dt <- fread(s2_csv, sep = "\t")
        
        # Align the samples by `samples_id`
        s2_dt$samples_id <- 1:nrow(s2_dt)
        s2_dt <- s2_dt[samples_id %in% lidar_dt$samples_id]
        
        # Append the current site's data to the aggregation tables
        agg_lidar <- rbind(agg_lidar, lidar_dt)
        agg_s2    <- rbind(agg_s2, s2_dt)
      }
      
      # If no data was aggregated, skip
      if (nrow(agg_lidar) == 0) next
      
      # Compute metrics on the aggregated data for all columns in agg_s2 at once
      metric_list <- lapply(colnames(agg_s2), function(colname) {
        vals_s2 <- agg_s2[[colname]]
        
        # Only keep valid data (complete cases)
        valid <- complete.cases(agg_lidar$lidar_values, vals_s2)
        
        # Compute metrics for the valid data
        m <- compute_metrics(agg_lidar$lidar_values[valid], vals_s2[valid])
        
        # Determine the ATBD flag based on column name
        atbd_flag <- str_count(colname, "\\.1(\\D|$)") == length(parms2test)
        
        # Return the metrics in a data.table format
        data.table(
          Site      = "All_sites",
          Norm      = norm,
          Depth     = depth,
          Method    = sampling_method,
          Column    = colname,
          R         = round(m["R"], 3),
          R2        = round(m["R2"], 3),
          RMSE      = round(m["RMSE"], 3),
          NRMSE     = round(m["NRMSE"], 3),
          Slope     = round(m["Slope"], 3),
          Intercept = round(m["Intercept"], 3),
          Bias      = round(m["Bias"], 3),
          MAE       = round(m["MAE"], 3),
          ATBD      = atbd_flag,
          stringsAsFactors = FALSE
        )
      })
      
      # Combine all results for this iteration
      combined_results <- rbindlist(list(combined_results, rbindlist(metric_list)))
    } # sampling_method
  }   # depth
}     # norm

# ------------------------------------------------------------------------------
# Save combined results for both per-site and All_sites
# ------------------------------------------------------------------------------
# out_csv <- file.path(output_dir,
#                      paste0("all_results_combined_", name_strategy, "28_04.csv"))
# fwrite(combined_results, out_csv)
# message("Saved combined results to: ", out_csv)

out_csv <- file.path(output_dir,
                     paste0("all_results_combined_", name_strategy, "13_05.csv"))
if (!file.exists(out_csv)) {
  write.csv(combined_results, out_csv, row.names = FALSE)
} else {
  write.table(combined_results, out_csv, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
}