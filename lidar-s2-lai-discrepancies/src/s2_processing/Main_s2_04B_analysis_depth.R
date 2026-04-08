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

# 1- Directories & input data
output_dir <- '../../03_RESULTS'
sites <- c('Aigoual', 'Blois', 'Mormal')
# sites <- c('Blois', 'Mormal')
# sites <- c("Aigoual")
datesAcq <- list('Aigoual' = '2021-07-11', 
                 'Blois' = '2021-06-14', 
                 'Mormal' = '2021-06-14')

# name_strategy <- 'LIDFa_only'
# name_strategy <- 'LIDFa_BROWN'
# name_strategy <- 'LAI_LMA'
# name_strategy <- 'LIDFa_lai'
# name_strategy <- "LIDFa_lai_LMA"
# parms2test <- c("LIDFa", "lai", "LMA")
# name_strategy <- c("BROWN_CHL_N_psoil_q")
# parms2test <- c("BROWN", "CHL", "N", "psoil", "q")
name_strategy <- c("LIDFa_lai_LMA_BROWN_N_CHL_psoil_q")
parms2test <- c("LIDFa", "lai", "LMA", "BROWN", "N", "CHL", "psoil", "q")
nbSamplesS2Refl <- 5000

strategies <- readRDS(file.path(output_dir, "Simulations_Strategy",
                                name_strategy, "simulation_strategy.rds"))
prosail_models_dir <- file.path(output_dir, sites, 
                                'PROSAIL_Optimization', 'PROSAIL_Models',
                                name_strategy)
names(prosail_models_dir) <- sites
lapply(X = prosail_models_dir, dir.create, showWarnings = F, recursive = T)

# sampling_methods <- c('random', 'stratified', 'stratified_uniform')
# norm_methods <- c('DTM', 'DSM')

sampling_methods <- c('stratified_uniform')
# norm_methods <- c('DSM_Above20','DSM_keepTreesAbove20', 'DTM', 'DSM_keepTrees', 'DSM')
# norm_methods <- c('DTM', 'DTM_keepTrees', 'DSM')
norm_methods <- c('DSM_keepTrees', 'DTM_keepTrees')
# norm_methods <- 'DSM_keepTrees'
depths <- 1:38
# depths <- c(2, 5, 8, 10, 12, 15, 18, 20, 22, 25, 28, 30)

# Initialize combined results container
combined_results <- data.frame()
all_sites_data <- list()

# Loop through sites, normalization methods, depths, and sampling methods
for (site in c(sites, "All_sites")) {
  print(site)
  
  if (site == "All_sites") {
    aggregated_lidar_lai <- aggregated_s2_lai <- data.frame()
  }
  
  for (norm in norm_methods) {
    # depths <- 1:35
    # depths <- 31:35
    # if (norm == "DSM_Above20"){
    #   # depths <- 1:20
    #   next
    # }
    for (depth in depths) {
      print(depth)
      for (sampling_method in sampling_methods) {
        
        # File paths
        lidar_csv_path <- file.path(output_dir, site, 'PROSAIL_Optimization', 
                                    'sampling', paste0("PAD_", norm, "_Depth_", 
                                                       depth, "_Samples_", sampling_method,
                                                       '_nbSamples_', nbSamplesS2Refl, ".csv"))
        s2_csv_path <- file.path(output_dir, site, 'PROSAIL_Optimization', 
                                 'PROSAIL_Models', name_strategy, 
                                 paste0('LAI_estimated_', sampling_method,
                                        '_nbSamples_', nbSamplesS2Refl, '.csv'))
        
        # Check if both files exist
        if (!file.exists(lidar_csv_path) || !file.exists(s2_csv_path)) {
          next  # Skip to the next iteration if files are missing
        }
        
        # Load data
        lidar_lai <- fread(lidar_csv_path, header = TRUE, sep = "\t")
        s2_lai_init <- fread(s2_csv_path, header = TRUE, sep = "\t")
        s2_lai <- s2_lai_init[lidar_lai$samples_id, ]
        
        if (nrow(lidar_lai) != nrow(s2_lai)) {
          stop("Number of rows in LiDAR and Sentinel-2 data do not match for site: ", site)
        }
        
        if (site == "All_sites") {
          aggregated_lidar_lai <- rbind(aggregated_lidar_lai, lidar_lai)
          aggregated_s2_lai <- rbind(aggregated_s2_lai, s2_lai)
        }
        
        # Initialize results container
        metrics <- data.frame(Column = colnames(s2_lai), 
                              R = NA, 
                              R2 = NA, 
                              RMSE = NA,
                              NRMSE = NA,  # Updated to NRMSE
                              Slope = NA, 
                              Intercept = NA,
                              Bias = NA,
                              ATBD = NA,
                              stringsAsFactors = FALSE)  
        
        # Calculate metrics for each column
        for (i in seq_along(colnames(s2_lai))) {
          # print(i)
          s2_column <- s2_lai[[i]]
          
          # Exclude NA pairs
          valid_idx <- complete.cases(lidar_lai$lidar_values, s2_column)
          lidar_values <- lidar_lai$lidar_values[valid_idx]
          s2_values <- s2_column[valid_idx]
          
          # Compute metrics
          lm_model <- lm(lidar_values ~ s2_values)
          r <- cor(lidar_values, s2_values, use = "complete.obs", method = "pearson")
          r2 <- summary(lm_model)$r.squared
          rmse <- sqrt(mean((s2_values - lidar_values)^2, na.rm = TRUE))
          nrmse <- sqrt(mean((s2_values - lidar_values)^2, na.rm = TRUE)) / IQR(s2_values, na.rm = TRUE)
          slope <- coef(lm_model)[2]
          intercept <- coef(lm_model)[1]
          bias <- mean(s2_values - lidar_values, na.rm = TRUE)
          
          # Check if ATBD configuration
          is_atbd <- stringr::str_count(colnames(s2_lai)[i], "\\.1(\\D|$)") == length(parms2test)
          
          # Store results
          metrics[i, 1] <- as.character(colnames(s2_lai)[i])
          metrics[i, 2:8] <- as.numeric(c(r, r2, rmse, nrmse, slope, intercept, bias))
          metrics[i, 9] <- is_atbd  # Store ATBD status
        }
        
        # Find best R2 models (handling ties)
        best_r2_candidates <- metrics[metrics$R == max(metrics$R), ]
        if (nrow(best_r2_candidates) > 1) {
          best_r2_candidates <- best_r2_candidates[which.min(best_r2_candidates$NRMSE), ]
        }
        if (nrow(best_r2_candidates) > 1) {
          best_r2_row <- best_r2_candidates[which.min(abs(best_r2_candidates$Bias)), ]
        } else {
          best_r2_row <- best_r2_candidates
        }
        best_r2_row$Criterion <- "Best R2"
        best_r2_row$Depth <- depth
        
        # Find best NRMSE models (handling ties)
        best_nrmse_candidates <- metrics[metrics$NRMSE == min(metrics$NRMSE), ]
        if (nrow(best_nrmse_candidates) > 1) {
          best_nrmse_candidates <- best_nrmse_candidates[which.max(best_nrmse_candidates$R), ]
        }
        if (nrow(best_nrmse_candidates) > 1) {
          best_nrmse_row <- best_nrmse_candidates[which.min(abs(best_nrmse_candidates$Bias)), ]
        } else {
          best_nrmse_row <- best_nrmse_candidates
        }
        best_nrmse_row$Criterion <- "Best NRMSE"
        best_nrmse_row$Depth <- depth
        
        # Combine best models
        top_models <- rbind(best_r2_row, best_nrmse_row)
        numeric_columns <- c("R", "R2", 'RMSE', "NRMSE", "Slope", "Intercept", "Bias")
        metrics[, numeric_columns] <- round(metrics[, numeric_columns], 2)
        top_models[, numeric_columns] <- round(top_models[, numeric_columns], 2)
        
        # Save metrics and top models
        result_dir <- file.path(output_dir, site, 'PROSAIL_Optimization', 
                                'sampling_results', name_strategy,
                                norm, paste0("Depth_", depth))
        dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
        
        write.csv(metrics, file.path(result_dir, paste0(sampling_method, "_metrics.csv")), row.names = FALSE)
        write.csv(top_models, file.path(result_dir, paste0(sampling_method, "_top_models.csv")), row.names = FALSE)
        
        # Generate scatter plots for best models
        for (i in 1:nrow(top_models)) {
          model_info <- top_models[i, ]
          model_column <- model_info$Column
          criterion <- model_info$Criterion
          
          valid_idx <- complete.cases(lidar_lai$lidar_values, s2_lai[[model_column]])
          lidar_values <- lidar_lai$lidar_values[valid_idx]
          s2_values <- s2_lai[[model_column]][valid_idx]
          
          plot_title <- paste(criterion, "-", model_column, 
                              "\nSite:", site, "Norm:", norm, 
                              "Depth:", depth, "Method:", sampling_method,
                              "\nATBD:", ifelse(model_info$ATBD, "Yes", "No"))
          
          plot_filename <- file.path(result_dir, 
                                     paste0(sampling_method, "_", tolower(gsub(" ", "_", criterion)), 
                                            "_scatter_", model_column, ".png"))
          
          create_scatterplot(s2_values, lidar_values, plot_title, plot_filename)
        }
        
        # Generate scatter plots for ATBD configurations (even if not best)
        # atbd_models <- metrics[metrics$ATBD == TRUE, ]
        # for (i in 1:nrow(atbd_models)) {
        #   model_info <- atbd_models[i, ]
        #   model_column <- model_info$Column
        #   
        #   valid_idx <- complete.cases(lidar_lai$lidar_values, s2_lai[[model_column]])
        #   lidar_values <- lidar_lai$lidar_values[valid_idx]
        #   s2_values <- s2_lai[[model_column]][valid_idx]
        #   
        #   plot_title <- paste("ATBD Configuration -", model_column, 
        #                       "\nSite:", site, "Norm:", norm, 
        #                       "Depth:", depth, "Method:", sampling_method,
        #                       "\nATBD: Yes")
        #   
        #   plot_filename <- file.path(result_dir, 
        #                              paste0(sampling_method, "_atbd_scatter_", model_column, ".png"))
        #   
        #   create_scatterplot(s2_values, lidar_values, plot_title, plot_filename)
        # }
        
        # Append top performers to the log
        log_file <- file.path(output_dir, paste0("top_performers_log_", name_strategy, ".csv"))
        top_metrics <- data.frame(Site = site, Norm = norm, Depth = depth, Method = sampling_method, top_models)
        if (!file.exists(log_file)) {
          write.csv(top_metrics, log_file, row.names = FALSE)
        } else {
          write.table(top_metrics, log_file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
        }
        
        # Append metrics to combined results
        combined_results <- rbind(combined_results, 
                                  cbind(Site = site, Norm = norm, 
                                        Method = sampling_method, Depth = depth, metrics))
        
        # Save metadata
        metadata <- paste("Date:", Sys.Date(), 
                          "\nSite:", site, 
                          "\nNorm:", norm, 
                          "\nDepth:", depth, 
                          "\nMethod:", sampling_method, 
                          "\nCode Version:", system("git rev-parse HEAD", intern = TRUE))
        writeLines(metadata, file.path(result_dir, "metadata.txt"))
      }
    }
  }
}

# Plot NRMSE and R²
top_performers <- read.csv(log_file)
top_performers_filtered <- top_performers[, c("R", "NRMSE", "Site", "Norm", "Depth", "Method")]
for (site in sites) {
  site_data <- top_performers_filtered[top_performers_filtered$Site == site, ]
  site_data$logNRMSE <- log10(site_data$NRMSE)
  
  scatter_plot <- ggplot(site_data, aes(x = R, y = logNRMSE, color = Depth, shape = Norm)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(title = paste("R vs log(NRMSE) for Top Performers - Site:", site),
         x = "R",
         y = "log(NRMSE)",
         color = "Depth",
         shape = "Norm") +
    theme_bw(base_size = 14)
  plot_filename <- file.path(output_dir, paste0("R_vs_logNRMSE_top_performers_", 
                                                name_strategy, "_", site, "24_04.png"))
  ggsave(plot_filename, scatter_plot, width = 10, height = 10, dpi = 300)
}

# Save combined results
file_save <- file.path(output_dir, paste0("all_results_combined_", 
                                          name_strategy, "24_04.csv"))
# write.csv(combined_results, file_save, row.names = FALSE)
if (!file.exists(file_save)) {
  write.csv(combined_results, file_save, row.names = FALSE)
} else {
  write.table(combined_results, file_save, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
}
stop()

# ---------------------------- Visualizations ----------------------------------
# combined_results top_metrics
top_metrics <- read.csv(log_file, header = TRUE)
global_min <- min(top_metrics$R2, na.rm = TRUE)
global_max <- max(top_metrics$R2, na.rm = TRUE)
metrics <- c("R2", "NRMSE", "Bias")

#### HEATMAP
# Aggregate metrics (e.g., mean R²) across sites for each combination
for (site in unique(top_metrics$Site)) {
  # Filter data for the current site
  site_data <- top_metrics %>% filter(Site == site)
  
  for (metric in metrics) {
    # Aggregate the chosen metric for each combination of Norm, Depth, and Method
    heatmap_data <- site_data %>%
      group_by(Norm, Depth, Method) %>%
      summarise(Mean_Metric = mean(!!sym(metric), na.rm = TRUE), .groups = 'drop')
    
    # Convert Depth to a factor for proper ordering on the heatmap
    heatmap_data$Depth <- factor(heatmap_data$Depth, levels = unique(heatmap_data$Depth))
    
    # Define the heatmap plot
    heatmap_plot <- ggplot(heatmap_data, aes(x = Depth, y = Norm, fill = Mean_Metric)) +
      geom_tile(color = "white", linewidth = 0.5) +  # Add tiles with white borders
      facet_wrap(~ Method, ncol = 1) +              # Facet by sampling method
      scale_fill_gradient(low = "yellow", high = "red", name = paste("Mean", metric)) +
      labs(
        title = paste("Performance Heatmap for Site:", site, "-", metric),
        x = "Depth",
        y = "Normalization Method",
        fill = paste("Mean", metric)
      ) +
      theme_bw() +
      theme(
        text = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
        strip.text = element_text(size = 12, face = "bold"), # Facet label styling
        panel.grid.major = element_blank(),                 # Remove grid lines
        panel.grid.minor = element_blank()
      )
    
    # Save the plot for the current site and metric
    ggsave(file.path(output_dir, site, paste0("performance_heatmap_", metric, "_", name_strategy, ".png")), 
           heatmap_plot, width = 10, height = 8, dpi = 300)
  }
}

# Aggregate metrics (e.g., mean R²) across sites for each combination
for (site in unique(top_metrics$Site)) {
  # Filter data for the current site
  site_data <- top_metrics %>% filter(Site == site)
  
  # Aggregate metrics (e.g., mean R²) for each combination of Norm, Depth, and Method
  curve_data <- site_data %>%
    group_by(Norm, Depth, Method) %>%
    summarise(Mean_R2 = mean(R2, na.rm = TRUE), .groups = 'drop')
  
  # Convert Depth and Norm to factors for proper ordering
  curve_data$Depth <- factor(curve_data$Depth, levels = unique(curve_data$Depth))
  curve_data$Norm <- factor(curve_data$Norm, levels = unique(curve_data$Norm))
  
  # Create the curve plot
  curve_plot <- ggplot(curve_data, aes(x = Depth, y = Mean_R2, color = Norm, group = Norm)) +
    geom_line(size = 1) +                               # Add lines for the curves
    geom_point(size = 3) +                              # Add points for visibility
    facet_wrap(~ Method, ncol = 1) +                    # Facet by sampling method
    scale_color_manual(values = c("blue", "green", "red"), name = "Normalization Method") +
    labs(
      title = paste("Performance Curve for Site:", site),
      x = "Depth",
      y = "Mean R²",
      color = "Normalization Method"
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 12, family = "serif"),
      axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
      strip.text = element_text(size = 12, face = "bold"), # Facet label styling
      panel.grid.major = element_blank(),                 # Remove grid lines
      panel.grid.minor = element_blank()
    )
  
  # Save the plot for the current site
  ggsave(file.path(output_dir, site, paste0("performance_curve_", name_strategy, ".png")), 
         curve_plot, width = 10, height = 8, dpi = 300)
}

#### IMPORTANCE
for (site in unique(top_metrics$Site)) {
  # Filter data for the current site
  site_data <- top_metrics %>% filter(Site == site)
  
  # Aggregate parameter importance for the current site
  parameter_importance <- site_data %>%
    group_by(Column) %>%  # Group by parameter name
    summarise(
      Frequency = n(),  # Count how often the parameter appears in the top 5%
      Mean_R2 = mean(R2, na.rm = TRUE),  # Mean R² for the parameter
      .groups = 'drop'
    ) %>%
    arrange(desc(Frequency))  # Sort by frequency in descending order
  parameter_importance <- parameter_importance %>% filter(Frequency > 1)
  
  # Create the parameter importance ranking plot
  importance_plot <- ggplot(parameter_importance, aes(x = reorder(Column, Frequency), y = Frequency, fill = Mean_R2)) +
    geom_bar(stat = "identity", color = "black") +  # Bar plot
    scale_fill_gradient(low = "yellow", high = "red", name = "Mean R²", 
                        limits = c(global_min, global_max)) +
    labs(
      title = paste("Parameter Importance Ranking for Site:", site),
      x = "Parameter",
      y = "Frequency in Top 5% Models",
      fill = "Mean R²"
    ) +
    coord_flip() +  # Flip axes for better readability
    theme_bw() +
    theme(
      text = element_text(size = 12, family = "serif"),
      axis.text.y = element_text(face = "bold"),  # Bold parameter names
      panel.grid.major.y = element_blank(),  # Remove horizontal grid lines
      panel.grid.minor.y = element_blank(),
      legend.position = "bottom"
    )
  
  # Save the plot for the current site
  ggsave(file.path(output_dir, site, paste0("parameter_importance_ranking_", name_strategy, ".png")), 
         importance_plot, width = 10, height = 6, dpi = 300)
}

# Prepare data for radar plot
radar_data <- combined_results %>%
  group_by(Method, Norm) %>%  # Group by sampling method and normalization
  summarise(
    Slope = mean(Slope, na.rm = TRUE),
    Intercept = mean(Intercept, na.rm = TRUE),
    Bias = mean(Bias, na.rm = TRUE),
    RMSE = mean(RMSE, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  unite("Group", Method, Norm, sep = " + ") %>%  # Combine method and norm into group names
  select(Group, Slope, Intercept, Bias, RMSE) %>%
  mutate(across(-Group, scale)) %>%  # Normalize metrics (Z-score)
  as.data.frame()

# Define axis ranges (customize based on your data ranges)
max_values <- c(1.5, 3, 2, 5)  # Custom max for Slope, Intercept, Bias, RMSE
min_values <- c(0.5, -1, -1, 0)  # Custom min values

# Create radar chart base
radar_df <- rbind(max_values, min_values, radar_data[-1])  # First 2 rows define ranges

# Set colors for groups
colors <- viridis::viridis(nrow(radar_data))
par(mar = c(1, 1, 1, 1))

# Create radar plot
radarchart(
  radar_df,
  axistype = 1,
  pcol = colors,
  plwd = 2,
  plty = 1,
  cglcol = "grey",
  cglty = 1,
  axislabcol = "grey",
  caxislabels = seq(-2, 2, 1),  # Custom axis labels
  vlcex = 0.8,
  title = "Slope-Bias Radar Plot: Configuration Comparison"
)

# Add legend
legend(
  "topright",
  legend = radar_data$Group,
  bty = "n",
  pch = 20,
  col = colors,
  text.col = "black",
  cex = 0.8,
  pt.cex = 1.5
)

# Add reference lines for ideal values (slope=1, bias=0)
abline(h = 0, col = "blue", lty = 2)  # Bias=0 reference
abline(h = 1, col = "green", lty = 2)  # Slope=1 reference (if normalized)

# ------------------------------- Analysis -------------------------------------
# Main analysis
results_list <- list()

for (site in sites) {
  # Load combined results for the site
  site_results <- combined_results %>%
    filter(Site == site) %>%
    mutate(Is_ATBD = str_detect(col, "\\.1(\\D|$)"))
  
  # Find optimal depths for each normalization method
  optimal_depths <- site_results %>%
    group_by(Norm, Depth) %>%
    summarise(Mean_R2 = mean(R2, na.rm = TRUE), .groups = "drop") %>%
    group_by(Norm) %>%
    slice_max(Mean_R2, n = 1) %>%
    rename(Optimal_Depth = Depth)
  
  # ATBD Analysis
  atbd_results <- site_results %>%
    filter(Is_ATBD) %>%
    left_join(optimal_depths, by = "Norm") %>%
    filter(Depth == Optimal_Depth)
  
  # Optimal Configuration Analysis (non-ATBD)
  optimal_results <- site_results %>%
    filter(!Is_ATBD) %>%
    group_by(Norm, Method, Depth, Column) %>%
    summarise(Max_R2 = max(R2), .groups = "drop") %>%
    left_join(optimal_depths, by = "Norm") %>%
    filter(Depth == Optimal_Depth) %>%
    group_by(Norm) %>%
    slice_max(Max_R2, n = 1)
  
  # Generate plots
  plot_dir <- file.path(output_dir, site, "PROSAIL_Optimization", "Figures")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ATBD Plots
  for (norm in norm_methods) {
    current_data <- atbd_results %>% filter(Norm == norm)
    if (nrow(current_data) > 0) {
      depth <- unique(current_data$Optimal_Depth)
      
      # Load corresponding data
      lidar_path <- file.path(output_dir, site, "PROSAIL_Optimization", 
                              "sampling", paste0("PAD_", norm, "_Depth_", depth, "_Samples_random.csv"))
      s2_path <- file.path(output_dir, site, "PROSAIL_Optimization",
                           "PROSAIL_Models", name_strategy, "LAI_estimated_random.csv")
      
      if (file.exists(lidar_path) && file.exists(s2_path)) {
        lidar_data <- read.csv(lidar_path, sep = "\t")
        s2_data <- read.csv(s2_path, sep = "\t")
        
        # Find ATBD column
        atbd_col <- grep("\\.1$", colnames(s2_data), value = TRUE)[1]
        
        create_scatterplot(s2_data[[atbd_col]], lidar_data$lidar_values,
                           paste("ATBD vs LiDAR -", site, "-", norm),
                           file.path(plot_dir, paste0("ATBD_vs_LiDAR_", site, "_", norm, ".png")))
      }
    }
  }
  
  # Optimal Configuration Plots
  for (i in 1:nrow(optimal_results)) {
    row <- optimal_results[i,]
    lidar_path <- file.path(output_dir, site, "PROSAIL_Optimization", 
                            "sampling", paste0("PAD_", row$Norm, "_Depth_", 
                                               row$Optimal_Depth, "_Samples_", row$Method, ".csv"))
    s2_path <- file.path(output_dir, site, "PROSAIL_Optimization",
                         "PROSAIL_Models", name_strategy, 
                         paste0("LAI_estimated_", row$Method, ".csv"))
    
    if (file.exists(lidar_path) && file.exists(s2_path)) {
      lidar_data <- read.csv(lidar_path, sep = "\t")
      s2_data <- read.csv(s2_path, sep = "\t")
      
      create_scatterplot(lidar_data$lidar_values, s2_data[[row$Column]],
                         paste("Optimal vs LiDAR -", site, "-", row$Norm),
                         file.path(plot_dir, paste0("Optimal_vs_LiDAR_", site, "_", row$Norm, ".png")))
    }
  }
  
  # Save results
  results_list[[site]] <- list(
    ATBD = atbd_results,
    Optimal = optimal_results
  )
}

# Generate summary report
summary_report <- file.path(output_dir, paste0("LAI_Validation_Summary_", name_strategy, ".html"))
rmarkdown::render(
  input = system.file("rmd", "validation_report.Rmd"),
  output_file = summary_report,
  params = list(results = results_list)
)