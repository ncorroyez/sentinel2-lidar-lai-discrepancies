# clean environment
rm(list = ls(all=TRUE)); gc()
if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getSourceEditorContext()$path))
lapply(X = as.list(list.files(path = '../libraries', full.names = T)), FUN = source, verbose = F)
library(terra)

# 1- Directories & input data
output_dir <- '../../03_RESULTS'
sites <- c('Aigoual', 'Blois', 'Mormal')
# sites <- c('Aigoual', 'Blois')
sites <- c("Aigoual")
datesAcq <- list('Aigoual' = '2021-07-11', 
                 'Blois' = '2021-06-14', 
                 'Mormal' = '2021-06-14')

# name_strategy <- 'LIDFa_only'
# name_strategy <- 'LAI_only'
# name_strategy <- 'LMA_only'
# name_strategy <- 'BROWN_only'
# name_strategy <- 'LIDFa_BROWN'
# name_strategy <- 'LAI_LMA'
name_strategy <- "LIDFa_lai_LMA"
nbSamplesS2Refl <- 5000

strategies <- readRDS(file.path(output_dir, "Simulations_Strategy",
                                name_strategy, "simulation_strategy.rds"))
prosail_models_dir <- file.path(output_dir, sites, 
                                'PROSAIL_Optimization', 'PROSAIL_Models',
                                name_strategy)
names(prosail_models_dir) <- sites
lapply(X = prosail_models_dir, dir.create, showWarnings = F, recursive = T)
Samplingmethods <- c('random', 'stratified', 'stratified_uniform')
Samplingmethods <- 'random'

for (site in sites){
  for (Samplingmethod in Samplingmethods){
    s2_csv_path <- file.path(prosail_models_dir[site], 
                             paste0('LAI_estimated_', Samplingmethod,
                                    '_nbSamples_', nbSamplesS2Refl, '.csv'))
    s2_lai <- read.csv(s2_csv_path, header = T, sep = "\t")
    lidar_csv_path <- file.path(output_dir, site, 'PROSAIL_Optimization', 
                                'sampling', paste0('LiDAR_LAI_Samples_', 
                                                   Samplingmethod,
                                                   '_nbSamples_', 
                                                   nbSamplesS2Refl, '.csv'))
    lidar_lai <- read.csv(lidar_csv_path, header = T, sep = "\t")
    
    if (nrow(lidar_lai) != nrow(s2_lai)) {
      stop("Number of rows in Lidar LAI and Sentinel-2 LAI data do not match for site: ", site)
    }
    
    # Initialize results container
    metrics <- data.frame(Column = colnames(s2_lai), 
                          R = NA, 
                          R2 = NA, 
                          RMSE = NA, 
                          Slope = NA, 
                          Intercept = NA,
                          Bias = NA)
    
    # Calculate metrics for each column
    for (i in seq_along(colnames(s2_lai))) {
      s2_column <- s2_lai[[i]]
      
      # Extract parameters (e.g., LIDFa, lai) from the column name
      colname_split <- strsplit(colnames(s2_lai)[i], "_")[[1]]
      param_list <- list()
      for (param in colname_split) {
        if (grepl("\\.", param)) {
          param_parts <- strsplit(param, "\\.")[[1]]
          param_name <- param_parts[1]
          param_value <- as.numeric(param_parts[2])
          param_list[[param_name]] <- param_value
        }
      }
      
      # Find corresponding strategies for the parameters
      strategy_info <- lapply(names(param_list), function(param_name) {
        if (param_name %in% names(strategies$combination)) {
          param_value <- param_list[[param_name]]
          strategy_info <- strategies$combination[[param_name]]
          return(strategy_info[[param_value]])
        }
        return(NULL)  # If no strategy found for the parameter
      })
      
      # Exclude NA pairs
      valid_idx <- complete.cases(lidar_lai$lidar_values, s2_column)
      lidar_values <- lidar_lai$lidar_values[valid_idx]
      s2_values <- s2_column[valid_idx]
      
      # Compute metrics
      lm_model <- lm(s2_values ~ lidar_values)
      r <- round(cor(lidar_values, s2_values, use = "complete.obs", method = "pearson"), 2)
      r2 <- round(summary(lm_model)$r.squared, 2)
      rmse <- round(sqrt(mean((lidar_values - s2_values)^2)), 2)
      slope <- round(coef(lm_model)[2], 2)
      intercept <- round(coef(lm_model)[1], 2)
      bias <- round(mean(lidar_values - s2_values, na.rm = TRUE), 2)
      
      # Store results
      metrics[i, 1] <- as.character(colnames(s2_lai)[i])
      metrics[i, 2:7] <- as.numeric(c(r, r2, rmse, slope, intercept, bias))
      
      # Add parameter strategies to the column name for printing
      # strategy_text <- paste(sapply(strategy_info, function(x) paste(names(x), "=", x, collapse = ", ")), collapse = " | ")
      # metrics$Column[i] <- paste(metrics$Column[i], "(", strategy_text, ")", sep = " ")
    }
    
    # Analyse top 5% according to a criterion
    n <- ceiling(0.05 * nrow(metrics))
    criterion <- metrics$R2 # Choose between R, R2, RMSE, Slope, Intercept, Bias
    criterion <- metrics$RMSE
    criterion <- metrics$Slope
    top_5_percent <- head(metrics[order(-criterion), ], n)
    
    # Print results for the site
    print(paste("Site:", site))
    # print("All Metrics:")
    # print(metrics)
    print("Top 5% Metrics:")
    print(top_5_percent)
  }
}

# ---------------------------- Visualizations ----------------------------------
# Bar plot for R² or RMSE
# bar_plot <- ggplot(data = metrics, aes(x = Column, y = R2)) +
#   geom_bar(stat = "identity", fill = "skyblue") +
#   labs(title = paste("R² by Parameter Combination for", site),
#        x = "Parameter Combination",
#        y = "R²") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(file.path(output_dir, paste0("barplot_", name_strategy, ".png")), 
#        bar_plot, width = 10, height = 8, dpi = 300)