# ---
# title: "plot_initial_scatterplots.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-07-09"
# ---

# ----------- (Optional) Clear the environment and free memory -----------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --- Define working directory as the directory where the script is located ----
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# ------------------------------ Libraries -------------------------------------

library(lidR)
library(data.table)
library(raster)
library(plotly)
library(terra)
library(viridis)
library(future)
library(rayshader)
library(broom)
library(Metrics)
library(ragg)
library(patchwork)
library(tidyr)

# --------------------------- Import useful functions --------------------------

source("libraries/functions_plots.R")

# --------------------------------- Setup --------------------------------------
results_dir   <- "../03_RESULTS"
figures_dir   <- "../04_FIGURES"
metrics_dir   <- "Metrics/Deciduous_Only"
s2_filename   <- "s2lai_summer_atbd_res_10_m.tif"
s2_common     <- "s2lai_summer_depth_study_common_res_10_m.tif"
s2_best_indiv <- "s2lai_summer_best_indiv_res_10_m.tif"
lidar_filename <- "PAD_Profiles_Classic/ladstack.tif"
pai_opt       <- "PAD_Profiles_dsm_keepTrees/PAD_34.5_40.tif"
mean_filename <- "mean_res_10_m.tif"
max_filename <- "max_res_10_m.tif"
# vci_filename  <- "vci_dsm_res_10_m.tif" # "vci_res_10_m.tif"
vci_filename  <- "vci_res_10_m.tif"
lcv_filename  <- "lcv_res_10_m.tif"
lskew_filename  <- "lskew_res_10_m.tif"
cvlad_filename  <- "cvlad_opt_depth_res_10_m.tif"  #"cvlad_opt_depth_res_10_m.tif"
dsm_sd_filename  <- "dsm_sd_res_10_m.tif"
sites         <- c("Aigoual", "Blois", "Mormal")
# sites         <- c("Aigoual", "Blois")
set.seed(42)
# -------------------------------- Process -------------------------------------
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)
#-----------------------------------------------------------------------------
# Generate scatterplots for each site (LAI comparison)
#-----------------------------------------------------------------------------
plots_lai <- lapply(sites, function(site) {
  message("Processing LAI for ", site)
  remove_y <- switch(site, "Aigoual" = F, "Blois" = T, "Mormal" = T)
  
  s2_path   <- file.path(results_dir, site, metrics_dir, s2_filename)
  lidar_path<- file.path(results_dir, site, metrics_dir, lidar_filename)
  mean_path <- file.path(results_dir, site, metrics_dir, mean_filename)
  pai_path <- file.path(results_dir, site, metrics_dir, pai_opt)
  
  s2    <- terra::rast(s2_path)
  lidar <- sum(terra::rast(lidar_path), na.rm = T)
  mean  <- terra::rast(mean_path)
  pai <- terra::rast(pai_path)
  
  na_index <- terra::countNA(s2) + terra::countNA(lidar) + terra::countNA(mean) + terra::countNA(pai)
  s2[na_index > 0]    <- NA
  lidar[na_index > 0] <- NA
  mean[na_index > 0]   <- NA
  pai[na_index > 0] <- NA
  
  plot_lai_scatter(values(lidar), values(s2), site,
                   xlim = c(0, 15), ylim = c(0, 15),
                   add_annotation = TRUE, annotation_pos = "topleft",
                   remove_y = remove_y)
})

combined_lai <- wrap_plots(plots_lai, nrow = 1, guides = "collect") &
  theme(legend.position = "none")

# Save combined LAI figure
ggsave(file.path(figures_dir, "All_Sites_LAI_Scatter.png"), combined_lai,
       dpi = 300, width = 18, height = 6, units = "in", type = "cairo")

#-----------------------------------------------------------------------------
# Generate scatterplots for ΔLAI vs Mean height
#-----------------------------------------------------------------------------
delta_plots <- lapply(sites, function(site) {
  message("Processing ΔLAI vs Mean for ", site)
  
  s2_path   <- file.path(results_dir, site, metrics_dir, s2_filename)
  lidar_path<- file.path(results_dir, site, metrics_dir, lidar_filename)
  mean_path <- file.path(results_dir, site, metrics_dir, mean_filename)
  pai_path <- file.path(results_dir, site, metrics_dir, pai_opt)
  
  s2    <- terra::rast(s2_path)
  lidar <- sum(terra::rast(lidar_path), na.rm = T)
  mean  <- terra::rast(mean_path)
  pai <- terra::rast(pai_path)
  
  delta <- lidar - s2
  na_index <- terra::countNA(s2) + terra::countNA(lidar) + terra::countNA(mean)
  delta[na_index > 0] <- NA
  mean[na_index > 0]  <- NA
  
  plot_lai_scatter(values(mean), values(delta), site,
                   xlab = "Mean Height (m)",
                   ylab = expression(Delta~LAI~"(LiDAR - S2)"),
                   xlim = c(0, 40), ylim = c(-15, 15),
                   add_annotation = TRUE, annotation_pos = "bottomright")
})

combined_delta <- wrap_plots(delta_plots, nrow = 1, guides = "collect") &
  theme(legend.position = "none")

# Save combined ΔLAI figure
ggsave(file.path(figures_dir, "All_Sites_DeltaLAI_vs_Mean.png"), combined_delta,
       dpi = 300, width = 18, height = 6, units = "in", type = "cairo")

# ------------------------------------------------------------------------------
# ΔLAI vs LCV class (quantile-based)
# ------------------------------------------------------------------------------
variance_stats_list <- list()
lcv_plots <- lapply(sites, function(site) {
  message("Processing ΔLAI vs LCV class (quantile) for ", site)
  
  depth <- switch(site,
                  "Aigoual" = "PAD_35.5_40.tif",
                  "Blois"   = "PAD_35.5_40.tif",
                  "Mormal"  = "PAD_31.5_40.tif")
  
  s2_best_path <- file.path(results_dir, site, metrics_dir, s2_best_indiv)
  # s2_best_path <- file.path(results_dir, site, metrics_dir, s2_filename)
  # lidar_path   <- file.path(results_dir, site, metrics_dir, lidar_filename)
  lidar_path   <- file.path(results_dir, site, metrics_dir,
                            "PAD_Profiles_dsm_keepTrees", depth)
  lcv_path     <- file.path(results_dir, site, metrics_dir, mean_filename)
  
  s2_best <- terra::rast(s2_best_path)
  # lidar   <- sum(terra::rast(lidar_path), na.rm = TRUE)
  lidar <- terra::rast(lidar_path)
  lcv     <- terra::rast(lcv_path)
  
  # Compute deltaLAI
  delta <- lidar - s2_best
  
  na_index <- terra::countNA(delta) + terra::countNA(lcv)
  delta[na_index > 0] <- NA
  lcv[na_index > 0]   <- NA
  
  # Build data.frame
  df <- data.frame(
    delta = as.numeric(values(delta)),
    s2    = as.numeric(values(s2_best)),
    lidar = as.numeric(values(lidar)),
    lcv   = as.numeric(values(lcv))
  )
  df <- df[complete.cases(df), ]
  
  if (nrow(df) < 10) return(NULL)  # Prevent failing on empty or tiny site
  
  # Create LCV classes based on quantiles (e.g., quartiles)
  quantiles <- quantile(df$lcv, probs = seq(0, 1, by = 0.25), na.rm = TRUE)
  quantiles <- unique(quantiles)
  if (length(quantiles) < 2) return(NULL)
  
  df$lcv_class <- cut(df$lcv, breaks = quantiles, include.lowest = TRUE, dig.lab = 2)
  
  # Pearson correlation, normalized RMSE, bias, and slope per LCV class
  stats_results <- by(df, df$lcv_class, function(sub_df) {
    # Filter NA
    valid <- complete.cases(sub_df$lidar, sub_df$s2)
    x <- sub_df$s2[valid]
    y <- sub_df$lidar[valid]
    
    cor_val   <- cor(y, x, method = "pearson")
    rmse_val  <- sqrt(mean((y - x)^2)) / IQR(y)
    bias_val  <- mean(y - x)
    slope_val <- coef(lm(y ~ x))[2]
    
    data.frame(
      class = as.character(unique(sub_df$lcv_class)),
      cor   = round(cor_val, 2),
      rmse  = round(rmse_val, 2),
      bias  = round(bias_val, 2),
      slope = round(slope_val, 2),
      n     = length(y)
    )
  })
  stats_results <- do.call(rbind, stats_results)
  print(stats_results)
  
  # Variance
  var_stats <- aggregate(df$delta, by = list(class = df$lcv_class), FUN = var)
  colnames(var_stats) <- c("class", "variance")
  var_stats$site <- site
  variance_stats_list[[site]] <<- var_stats
  
  # Plot
  p <- ggplot(df, aes(x = lcv_class, y = delta)) +
    geom_boxplot(outlier.shape = NA, fill = "seagreen3", alpha = 0.7) +
    labs(title = site,
         x = "LCV Quantile Class",
         y = expression(Delta~LAI~"(LiDAR - S2 best)")) +
    ylim(c(-10,15))+
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
})
lcv_plots <- Filter(Negate(is.null), lcv_plots)
combined_lcv <- wrap_plots(lcv_plots, nrow = 1)

# Plot variance per class and site
variance_df <- do.call(rbind, variance_stats_list)
p_var <- ggplot(variance_df, aes(x = class, y = variance, fill = site)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Variance of ΔLAI per LCV Class",
       x = "LCV Quantile Class",
       y = "Variance of ΔLAI") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_var)

# Save LCV boxplot figure
ggsave(file.path(figures_dir, "All_Sites_DeltaLAI_vs_VCIATBDclass_Quantile.png"), combined_lcv,
       dpi = 300, width = 18, height = 6, units = "in", type = "cairo")

# ------------------------------------------------------------------------------
ATBD   <- "s2lai_summer_atbd_res_10_m.tif"
Optimized <- "s2lai_summer_best_indiv_res_10_m.tif"
# Optimized     <- "s2lai_summer_depth_study_common_res_10_m.tif"

compare_all_sites <- function(metric = c("mean", "vci", "cvlad", "lcv", "dsm_sd")) {
  # LIMITER LE LAI LIDAR
  metric <- match.arg(metric)
  all_stats <- list()
  all_plot_data <- list()
  
  for (site in sites) {
    message("Processing ", site)
    
    depth <- switch(site,
                    "Aigoual" = "PAD_35.5_40.tif",
                    # "Blois"   = "PAD_35.5_40.tif",
                    "Blois"   = "PAD_25.5_40.tif",
                    "Mormal"  = "PAD_31.5_40.tif")
    
    lidar_path <- file.path(results_dir, site, metrics_dir,
                            "PAD_Profiles_dsm_keepTrees", depth)
    
    full_lidar_path <- file.path(results_dir, site, metrics_dir,
                            "lidarlai_res_10_m.tif")
    
    metric_file <- switch(metric,
                          mean = mean_filename,
                          vci  = vci_filename,
                          cvlad  = cvlad_filename,
                          lcv = lcv_filename,
                          dsm_sd = dsm_sd_filename)
    metric_path <- file.path(results_dir, site, metrics_dir, metric_file)
    
    lidar_rast  <- terra::rast(lidar_path)
    metric_rast <- terra::rast(metric_path)
    
    if (metric == "mean") {
      if (site == "Aigoual") {
        metric_rast[metric_rast > 30] <- 30
      } else {
        metric_rast[metric_rast > 40] <- 40
      }
    }
    
    for (s2_cfg in c("ATBD", "Optimized")) {
      s2_path <- file.path(results_dir, site, metrics_dir, get(s2_cfg))
      s2_rast <- terra::rast(s2_path)
      
      na1 <- terra::countNA(lidar_rast) + terra::countNA(s2_rast)
      lidar_rast[na1 > 0] <- NA
      s2_rast[na1 > 0] <- NA
      
      # Compute delta
      delta <- lidar_rast - s2_rast
      na_idx <- terra::countNA(delta) + terra::countNA(metric_rast)
      delta[na_idx > 0] <- NA
      metric_rast[na_idx > 0] <- NA
      
      df <- data.frame(
        delta  = as.numeric(values(delta)),
        s2     = as.numeric(values(s2_rast)),
        lidar  = as.numeric(values(lidar_rast)),
        metric = as.numeric(values(metric_rast))
      )
      df <- df[complete.cases(df), ]
      if (nrow(df) < 10) next
      
      # Manual class breaks
      breaks <- switch(metric,
                       mean = c(-Inf, 10, 20, 30, 40, Inf),
                       vci  = c(-Inf, 0.4, 0.6, Inf),
                       cvlad = c(-Inf, 0.2, 0.4, Inf),
                       lcv = c(-Inf, 0.005, 0.03, Inf),
                       dsm_sd = c(-Inf, 2.5, 5, Inf))
      labels <- switch(metric,
                       mean = c("<10", "10–20", "20–30", "30–40", ">40"),
                       vci  = c("<0.4", "0.4–0.6", ">0.6"),
                       cvlad = c("<0.2", "0.2–0.4", ">0.4"),
                       lcv = c("<0.005", "0.005–0.03", ">0.03"),
                       dsm_sd = c("<2.5", "2.5–5", ">5"))
      
      df$class  <- cut(df$metric, breaks = breaks, labels = labels, include.lowest = TRUE)
      df$class <- factor(df$class, levels = labels)
      df$site   <- site
      df$config <- s2_cfg
      
      # df <- do.call(rbind, lapply(split(df, df$class), function(sub_df) {
      #   if (nrow(sub_df) > 50000) {
      #     sub_df <- sub_df[sample(nrow(sub_df), 50000), ]
      #   }
      #   return(sub_df)
      # }))
      
      
      # ─── SUBSAMPLE AFTER CLASSING ────────────────────────────────────────
      # How many classes actually present?
      present_classes <- levels(df$class)[ tapply(df$class, df$class, length) > 0 ]
      NC <- length(present_classes)
      if (NC==0) next
      
      # Base quota per class, plus distribute remainder
      total_n   <- 15000
      base_n    <- floor(total_n/NC)
      rem       <- total_n - base_n*NC
      set.seed(42)
      
      sampled_idx <- unlist(lapply(seq_along(present_classes), function(i) {
        cls <- present_classes[i]
        idx <- which(df$class==cls)
        want <- base_n + ifelse(i <= rem, 1, 0)
        if (length(idx) <= want) return(idx)
        sample(idx, want)
      }))
      df <- df[sampled_idx, ]
      
      
      all_plot_data[[length(all_plot_data)+1]] <- df
      
      # Stats per class
      stats <- by(df, df$class, function(sub_df) {
        x <- sub_df$s2
        y <- sub_df$lidar
        data.frame(
          site   = site,
          config = s2_cfg,
          class  = as.character(unique(sub_df$class)),
          cor    = round(cor(y, x, method = "pearson"), 2),
          rmse   = round(sqrt(mean((y - x)^2)), 2),
          # nrmse  = round(sqrt(mean((y - x)^2)) / IQR(y), 2),
          bias   = round(mean(y - x), 2),
          slope  = round(coef(lm(y ~ x))[2], 2),
          n      = length(y)
        )
      })
      all_stats[[length(all_stats)+1]] <- do.call(rbind, stats)
    }
  }
  
  # Combine
  stats_df <- do.call(rbind, all_stats)
  plot_df  <- do.call(rbind, all_plot_data)
  plot_df$class <- factor(plot_df$class, levels = unique(unlist(switch(metric,
                                                                       mean = c("<10", "10–20", "20–30", "30–40", ">40"),
                                                                       vci  = c("<0.4", "0.4–0.6", ">0.6"),
                                                                       cvlad = c("<0.2", "0.2–0.4", ">0.4"),
                                                                       lcv = c("<0.005", "0.005–0.03", ">0.03"),
                                                                       dsm_sd = c("<2.5", "2.5–5", ">5")))))
  
  print(stats_df)
  write.csv(stats_df, file = paste0("delta_lai_stats_by_", metric, ".csv"), row.names = FALSE)
  
  # Plot: boxplot per class and config
  p <- ggplot(plot_df, aes(x = class, y = delta, fill = config)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, position = position_dodge(0.8), width = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
    facet_wrap(~site, scales = "free", ncol = 3) +
    labs(
      x = switch(metric,
                 mean = "Mean Height (meters)",
                 vci  = "Vertical Complexity Index (VCI)",
                 dsm_sd = expression(paste(CHM[SD], " (m)")),
                 paste("Class of", metric)),
      y = expression(paste(
        Delta,                      # Δ
        LAI,                        # LAI
        " (",                       # space + "("
        LAI[ALS_dopt],                   # LAI₍ₐₗₛ₎
        " − ",                     # minus sign with spaces
        LAI[S2],                    # LAI₍ₛ₂₎
        ")"                         # closing ")"
      )),
      fill = "Sentinel-2\nconfiguration"
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(-5, 5), breaks = seq(-5, 5, 1), expand = expansion(mult = c(0.02, 0.02))) +
    scale_fill_manual(values = c("ATBD" = "#1f78b4", "Optimized" = "#33a02c")) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      # legend.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      # strip.text = element_text(size = 13, face = "bold"),
      axis.title = element_text(size = 14)
      # plot.margin = margin(10, 10, 10, 10)
    )
  
  print(p)
  
  invisible(list(stats = stats_df, plot_data = plot_df, plot = p))
}

plot_summary_metrics_by_class <- function(stats_df, metric = c("mean", "vci", "cvlad", "lcv", "dsm_sd")) {
  metric <- match.arg(metric)
  
  # Ensure proper class order
  stats_df$class <- factor(
    stats_df$class,
    levels = switch(metric,
                    mean = c("<10", "10–20", "20–30", "30–40", ">40"),
                    vci  = c("<0.4", "0.4–0.6", ">0.6"),
                    cvlad = c("<0.2", "0.2–0.4", ">0.4"),
                    lcv = c("<0.005", "0.005–0.03", ">0.03"),
                    dsm_sd = c("<2.5", "2.5–5", ">5"))
  )
  stats_df <- as.data.frame(stats_df) %>%
    dplyr::rename(
      R     = cor,
      # NRMSE = nrmse,
      RMSE = rmse,
      Bias  = bias,
      Slope = slope
    ) %>%
    dplyr::select(site, config, class, R, RMSE, Bias, Slope)
  
  
  # Reshape for ggplot
  stats_long <- tidyr::pivot_longer(
    stats_df,
    cols = c(R, RMSE, Bias, Slope), #c(R, NRMSE, Bias, Slope),
    names_to = "metric_name",
    values_to = "value"
  )
  
  # Order metric_name factor levels
  stats_long$metric_name <- factor(stats_long$metric_name,
                                   # levels = c("R", "NRMSE", "Bias", "Slope"))
                                   levels = c("R", "RMSE", "Bias", "Slope"))
  
  # Label for x-axis
  x_label <- switch(metric,
                    mean = "Mean Height (m)",
                    vci  = "Vertical Complexity Index (VCI)",
                    dsm_sd = "CHM Standard Deviation")
  
  # Plot
  p <- ggplot(stats_long, aes(x = class, y = value, group = config, color = config)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    facet_grid(metric_name ~ site, scales = "free_y") +
    scale_color_manual(values = c("ATBD" = "#1f78b4", "Optimized" = "#33a02c")) +
    labs(
      x = x_label,
      y = NULL,
      color = "Sentinel-2 Configuration"
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 12)
    )
  return(p)
}

plot_density_scatter_by_class <- function(plot_data, metric = c("mean", "vci", "cvlad", "lcv", "dsm_sd")) {
  metric <- match.arg(metric)
  
  # Ensure proper factor order
  class_levels <- switch(metric,
                         mean = c("<10", "10–20", "20–30", "30–40", ">40"),
                         vci  = c("<0.4", "0.4–0.6", ">0.6"),
                         cvlad = c("<0.2", "0.2–0.4", ">0.4"),
                         lcv = c("<0.005", "0.005–0.03", ">0.03"),
                         dsm_sd = c("<2.5", "2.5–5", ">5")
  )
  
  # Ensure proper factor levels
  plot_data$class <- factor(plot_data$class, levels = class_levels)
  plot_data$config <- factor(plot_data$config)
  
  all_combos <- expand.grid(
    class = class_levels,
    config = levels(plot_data$config),
    site = unique(plot_data$site)
  )
  plot_data <- plot_data %>%
    right_join(all_combos, by = c("class", "config", "site"))
  
  # Get unique sites
  sites <- unique(plot_data$site)
  site_plots <- list()
  combined_plot <- NULL
  
  for (i in seq_along(sites)) {
    s <- sites[i]
    site_data <- subset(plot_data, site == s)
    
    # Determine axis labels based on position
    y_label <- if (i == 1) expression(LAI[S2]) else ""
    x_label <- if (i == 2) expression(LAI[ALS_dopt]) else ""
    
    p <- ggplot(site_data, aes(x = lidar, y = s2)) +
      geom_hex(bins = 200) +
      facet_grid(class ~ config) +
      scale_fill_viridis_c(option = "C", name = "Count", na.value = NA) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      coord_equal() +
      labs(
        x = x_label,
        y = y_label,
        title = s
      ) +
      xlim(c(0, 10)) + ylim(c(0, 10)) + 
      theme_bw(base_size = 20) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        # plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(t = 100)),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.margin = ggplot2::margin(t = 100, r = 5, b = 50, l = 5, unit = "pt")
      )
    
    combined_plot <- if (is.null(combined_plot)) p else combined_plot | p
  }
  
  return(combined_plot)
}

# Mean Height
result_mean <- compare_all_sites(metric = "mean")
plot_mean <- plot_summary_metrics_by_class(result_mean$stats, metric = "mean")
plot_mean_hex <- plot_density_scatter_by_class(result_mean$plot_data, metric = "mean")

# Cv_LAD
result_cvlad <- compare_all_sites(metric = "cvlad")
plot_cvlad <- plot_summary_metrics_by_class(result_cvlad$stats, metric = "cvlad")
plot_cvlad_hex <- plot_density_scatter_by_class(result_cvlad$plot_data, metric = "cvlad")

# VCI
result_vci <- compare_all_sites(metric = "vci")
plot_vci <- plot_summary_metrics_by_class(result_vci$stats, metric = "vci")
plot_vci_hex <- plot_density_scatter_by_class(result_vci$plot_data, metric = "vci")
# ggsave("vci_hex.png", plot_vci_hex, width = 12, height = 4)

# Lcv
result_lcv <- compare_all_sites(metric = "lcv")
plot_lcv <- plot_summary_metrics_by_class(result_lcv$stats, metric = "lcv")
plot_lcv_hex <- plot_density_scatter_by_class(result_lcv$plot_data, metric = "lcv")

# dsm_sd
result_dsm_sd <- compare_all_sites(metric = "dsm_sd")
plot_dsm_sd <- plot_summary_metrics_by_class(result_dsm_sd$stats, metric = "dsm_sd")
plot_dsm_sd_hex <- plot_density_scatter_by_class(result_dsm_sd$plot_data, metric = "dsm_sd")
ggsave("dsm_sd_hex.png", plot_dsm_sd_hex, width = 12, height = 4)

# for (site in sites){
#   cvlad <- rast(file.path("/home/corroyez/Documents/NC_Full/03_RESULTS", site, "Metrics/Deciduous_Only/cvlad_opt_depth_res_10_m.tif"))
#   vci <- rast(file.path("/home/corroyez/Documents/NC_Full/03_RESULTS", site, "Metrics/Deciduous_Only/vci_res_10_m.tif"))
#   
#   na_mask <- is.na(cvlad) | is.na(vci)
#   
#   # Apply the mask to both rasters
#   cvlad_matched <- mask(cvlad, na_mask, maskvalues = TRUE)
#   writeRaster(cvlad_matched, filename = file.path("/home/corroyez/Documents/NC_Full/03_RESULTS", site, "Metrics/Deciduous_Only/cvlad_opt_depth_res_10_m.tif"), 
#               overwrite = T)
# }

for (site in sites){
  ladstack <- rast(file.path("/home/corroyez/Documents/NC_Full/03_RESULTS", site, 
                             "Metrics/Deciduous_Only/PAD_Profiles_dsm_keepTrees/ladstack.tif"))
  vci <- rast(file.path("/home/corroyez/Documents/NC_Full/03_RESULTS", site, 
                        "Metrics/Deciduous_Only/vci_res_10_m.tif"))
  
  depth <- switch(site,
                  "Aigoual" = 5,
                  "Blois"   = 5,
                  "Mormal"  = 9)
  lad <- rast(file.path("/home/corroyez/Documents/NC_Full/03_RESULTS", site,
                        "Metrics/Deciduous_Only/lidarlai_optim_depth_res_10_m.tif"))
  
  na_mask <- is.na(lad) | is.na(vci)
  lad <- mask(lad, na_mask, maskvalues = TRUE)
  
  mean_rast <- focal(lad, w = 5, fun = mean, na.rm = TRUE)
  sd_rast <- focal(lad, w = 5, fun = sd, na.rm = TRUE)
  cvlad <- sd_rast / mean_rast
  
  # print(vci)
  # print("\n")
  # print(cvlad)
  # print("\n")
  
  cvlad <- mask(cvlad, is.na(cvlad) | is.na(vci), maskvalues = TRUE)
  
  # Optionally: save or plot
  # plot(cvlad, main = paste0("CV_LAD - ", site))
  writeRaster(cvlad, 
              filename = file.path("/home/corroyez/Documents/NC_Full/03_RESULTS", site, 
                                   "Metrics/Deciduous_Only/cvlad_opt_depth_2Jun_res_10_m.tif"), 
              overwrite = TRUE)
}

# ------------------------------------------------------------------------------

