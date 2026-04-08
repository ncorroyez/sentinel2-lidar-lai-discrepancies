# ---
# title: "three_factors_analysis.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2025-07-11"
# ---

# ----------- (Optional) Clear the environment and free memory -----------------
rm(list=ls(all=TRUE))
gc()

# --- Define working directory as the directory where the script is located ----
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# --------------------------------- Libraries ----------------------------------
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(purrr)
library(relaimpo)
library(ranger)
library(rpart)
library(rpart.plot)
library(hexbin)
library(patchwork)

# --------------------------------- Setup --------------------------------------
results_dir <- "../03_RESULTS"
metrics_dir <- "Metrics/Deciduous_Only"
# lidar_lai_file <- "lidarlai_res_10_m.tif"
lidar_lai_file <- "ladstack_classic.tif"
max_file <- "max_res_10_m.tif"
mns_sd_file <- "dsm_sd_res_10_m.tif"
s2_files <- list(
  LAI_S2_ATBD = "s2lai_summer_atbd_res_10_m.tif",
  LAI_S2_opt  = "s2lai_summer_best_indiv_res_10_m.tif"
)
s2_labels <- c(
  LAI_S2_ATBD = expression(LAI[S2_ATBD]),
  LAI_S2_opt  = expression(LAI[S2_opt])
)
lidar_labels <- c(
  LAI_ALS = expression(LAI[ALS]),
  LAI_ALS_dopt  = expression(LAI[ALS_dopt])
)
sites <- c("Aigoual", "Blois", "Mormal")
set.seed(42)

# ---------------------- 1. Load & Assemble Data ------------------------------
analyze_three_factor_relationships <- function() {
  all_data <- list()
  
  # Define the specific combinations to run
  combinations <- list(
    c(s2_config = "LAI_S2_ATBD", lidar_config = "LAI_ALS"),
    c(s2_config = "LAI_S2_ATBD", lidar_config = "LAI_ALS_dopt"),
    c(s2_config = "LAI_S2_opt", lidar_config = "LAI_ALS_dopt")
  )
  
  for (site in sites) {
    # optimal_fn <- switch(site,
    #                      "Aigoual" = "PAD_35.5_40.tif",
    #                      "Blois"= "PAD_35.5_40.tif",
    #                      "Mormal" = "PAD_31.5_40.tif")
    optimal_fn <- switch(site,
                         "Aigoual" = "PAD_36.5_40.tif",
                         "Blois"= "PAD_36.5_40.tif",
                         "Mormal" = "PAD_36.5_40.tif")
    # "Aigoual" = "PAD_34.5_40.tif",
    # "Blois" = "PAD_20.5_40.tif",
    # "Mormal" = "PAD_26.5_40.tif")
    base <- file.path(results_dir, site, metrics_dir)
    
    # load rasters
    rast_opt <- rast(file.path(base, "PAD_Profiles_dsm_keepTrees", optimal_fn))
    # rast_full <- rast(file.path(base, lidar_lai_file))
    rast_full <- sum(rast(file.path(base, lidar_lai_file)), na.rm = TRUE)
    rast_het <- rast(file.path(base, mns_sd_file))
    rast_maxh <- rast(file.path(base, max_file))
    
    # get global maximum LAI if you really need an an explicit upper bound:
    max_lai_value <- max(values(rast_full), na.rm = TRUE)
    
    # ----- Build the reference mask (using FULL lidar) -----
    ref_mask <- rast_full
    ref_mask[] <- is.na(rast_full[]) |
      is.na(rast_het[]) |
      is.na(rast_maxh[]) |
      (rast_full[] < 2) |# <-- threshold based on FULL
      (rast_full[] > max_lai_value) |
      (rast_maxh[] <= 10) |
      (rast_maxh[] >= 40)
    
    for (combo in combinations) {
      cfg <- combo["s2_config"]
      depth_tag <- combo["lidar_config"]
      
      # Use the normal S2 file for individual sites
      rast_s2 <- rast(file.path(base, s2_files[[cfg]]))
      rast_lid <- if (depth_tag == "LAI_ALS_dopt") rast_opt else rast_full
      
      # add S2 into mask
      na_mask <- ref_mask | is.na(rast_s2)
      
      # apply mask
      lid <- mask(rast_lid, na_mask, maskvalue=TRUE)
      s2  <- mask(rast_s2, na_mask, maskvalue=TRUE)
      het <- mask(rast_het, na_mask, maskvalue=TRUE)
      maxh <- mask(rast_maxh, na_mask, maskvalue=TRUE)
      
      # stack into a data.frame
      df <- data.frame(
        lidar = as.numeric(values(lid)),
        s2 = as.numeric(values(s2)),
        mns_sd = as.numeric(values(het)),
        height = as.numeric(values(maxh))
      ) %>%
        drop_na()
      
      # annotate
      df <- df %>%
        mutate(
          site = site,
          s2_config = cfg,
          lidar_config = depth_tag
        )
      
      all_data[[length(all_data) + 1]] <- df
      print(str(df))
    }
  }
  
  bind_rows(all_data)
}

# ---------------------- 2. Compute Metrics -----------------------------
compute_metrics <- function(df) {
  
  # --- Part 1: Binned Analysis (Sample once) ---
  # Create heterogeneity classes
  df_binned <- df %>%
    dplyr::mutate(
      mns_sd_class = cut(
        mns_sd,
        breaks = c(-Inf, 2.5, 5, Inf),
        labels = c("Low", "Medium", "High"),
        right = FALSE,
        include.lowest = TRUE
      )
    ) %>%
    dplyr::group_by(site, mns_sd_class, lidar_config, s2_config) %>%
    dplyr::mutate(row_id = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  # Sample IDs to ensure consistent representation across configs and sites
  set.seed(42)
  sampled_ids_binned <- df_binned %>%
    dplyr::group_by(site, mns_sd_class) %>%
    dplyr::distinct(row_id) %>%
    dplyr::slice_sample(n = 5000, replace = FALSE) %>%
    dplyr::ungroup()
  
  # Filter data using the sampled IDs
  sampled_data <- df_binned %>%
    dplyr::inner_join(sampled_ids_binned, 
                      by = c("site", "mns_sd_class", "row_id"))
  
  # --- Part 2: Create unified dataset with both binned + Total ---
  data_for_metrics <- sampled_data %>%
    dplyr::select(-row_id) %>%
    dplyr::bind_rows(
      sampled_data %>%
        dplyr::select(-row_id) %>%
        dplyr::mutate(mns_sd_class = "Total") %>%
        dplyr::distinct()
    )
  
  # --- Part 3: Calculate metrics ---
  calculate_summary_metrics <- function(data_frame) {
    data_frame %>%
      dplyr::group_by(site, s2_config, lidar_config, mns_sd_class) %>%
      dplyr::filter(dplyr::n() >= 10) %>%
      dplyr::summarise(
        R = cor(lidar, s2, use = "complete.obs"),
        RMSE = sqrt(mean((lidar - s2)^2, na.rm = TRUE)),
        Bias = mean(lidar - s2, na.rm = TRUE),
        Slope = coef(lm(lidar ~ s2))[2],
        .groups = "drop"
      ) %>%
      dplyr::ungroup()
  }
  
  metrics <- calculate_summary_metrics(data_for_metrics) %>%
    dplyr::mutate(
      across(c(R, RMSE, Bias, Slope), ~round(., 2)),
      mns_sd_class = factor(mns_sd_class,
                            levels = c("Low", "Medium", "High", "Total"))
    ) %>%
    dplyr::arrange(site, lidar_config, s2_config)
  
  # --- Part 4: Return ---
  return(list(data_for_metrics = data_for_metrics, 
              metrics = metrics))
}

compute_metrics <- function(df) {
  
  # --- Part 1: Binned Analysis (Sample once) ---
  # Create heterogeneity classes
  df_binned <- df %>%
    dplyr::mutate(
      mns_sd_class = cut(
        mns_sd,
        breaks = c(-Inf, 2.5, 5, Inf),
        labels = c("Low", "Medium", "High"),
        right = FALSE,
        include.lowest = TRUE
      )
    ) %>%
    dplyr::group_by(site, mns_sd_class, lidar_config, s2_config) %>%
    dplyr::mutate(row_id = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  # Sample IDs to ensure consistent representation across configs and sites
  set.seed(42)
  sampled_ids_binned <- df_binned %>%
    dplyr::group_by(site, mns_sd_class) %>%
    dplyr::distinct(row_id) %>%
    dplyr::slice_sample(n = 5000, replace = FALSE) %>%
    dplyr::ungroup()
  
  # Filter data using the sampled IDs
  sampled_data <- df_binned %>%
    dplyr::inner_join(sampled_ids_binned,
                      by = c("site", "mns_sd_class", "row_id"))
  
  # --- Part 2: Create unified dataset with both binned + Total ---
  data_for_metrics <- sampled_data %>%
    dplyr::select(-row_id) %>%
    dplyr::bind_rows(
      sampled_data %>%
        dplyr::select(-row_id) %>%
        dplyr::mutate(mns_sd_class = "Total") %>%
        dplyr::distinct()
    )
  
  # --- Part 3: Calculate metrics and p-values ---
  calculate_summary_metrics <- function(data_frame) {
    data_frame %>%
      dplyr::group_by(site, s2_config, lidar_config, mns_sd_class) %>%
      dplyr::filter(dplyr::n() >= 10) %>%
      dplyr::summarise(
        # Correlation and RMSE
        R = cor(lidar, s2, use = "complete.obs"),
        RMSE = sqrt(mean((lidar - s2)^2, na.rm = TRUE)),
        
        # Bias and p-value
        Bias = mean(lidar - s2, na.rm = TRUE),
        Bias_p_value = if (sum(!is.na(lidar - s2)) > 1) {
          t.test(lidar - s2)$p.value
        } else {
          NA
        },
        
        # Slope and p-value
        Slope = coef(lm(lidar ~ s2))[2],
        Slope_p_value = if (sum(!is.na(lidar) & !is.na(s2)) > 2) {
          model <- lm(lidar ~ s2)
          # Test if slope is significantly different from 1
          t_stat <- (coef(model)[2] - 1) / summary(model)$coefficients[2, 2]
          2 * pt(abs(t_stat), df = df.residual(model), lower.tail = FALSE)
        } else {
          NA
        },
        .groups = "drop"
      ) %>%
      dplyr::ungroup()
  }
  
  metrics <- calculate_summary_metrics(data_for_metrics) %>%
    dplyr::mutate(
      across(c(R, RMSE, Bias, Slope), ~round(., 2)),
      across(c(Bias_p_value, Slope_p_value), ~round(., 3)),
      mns_sd_class = factor(mns_sd_class,
                            levels = c("Low", "Medium", "High", "Total"))
    ) %>%
    dplyr::arrange(site, lidar_config, s2_config)
  
  # --- Part 4: Return ---
  return(list(data_for_metrics = data_for_metrics,
              metrics = metrics))
}

# -------------------------- 3. Main Execution ---------------------------------
# 1. Load data for the three individual sites using standard files
raw_data_individual_sites <- analyze_three_factor_relationships()

# 2. Prepare aggregated data for "Sites combined" for all combinations
# Use the raw data from the individual sites to create the aggregated data frames
agg_atbd_als <- raw_data_individual_sites %>%
  filter(s2_config == "LAI_S2_ATBD" & lidar_config == "LAI_ALS") %>%
  ungroup() %>%
  mutate(site = "Sites combined")

agg_atbd_dopt <- raw_data_individual_sites %>%
  filter(s2_config == "LAI_S2_ATBD" & lidar_config == "LAI_ALS_dopt") %>%
  ungroup() %>%
  mutate(site = "Sites combined")

# Aggregated data for the special LAI_S2_opt vs LAI_ALS_dopt case
# This requires loading the specific S2 file for each of the three sites
special_s2_file <- "s2lai_summer_depth_study_common_res_10_m.tif"
special_lidar_file <- "PAD_36.5_40.tif"

special_data_list <- list()

for (site in sites) {
  base <- file.path(results_dir, site, metrics_dir)
  
  rast_opt_combined <- rast(file.path(base, "PAD_Profiles_dsm_keepTrees", special_lidar_file))
  rast_full_combined <- sum(rast(file.path(base, lidar_lai_file)), na.rm = TRUE)
  rast_het_combined <- rast(file.path(base, mns_sd_file))
  rast_maxh_combined <- rast(file.path(base, max_file))
  rast_s2_combined <- rast(file.path(base, special_s2_file))
  
  # Re-apply the mask to the special combined data
  max_lai_value_combined <- max(values(rast_full_combined), na.rm = TRUE)
  ref_mask_combined <- rast_full_combined
  ref_mask_combined[] <- is.na(rast_full_combined[]) |
    is.na(rast_het_combined[]) |
    is.na(rast_maxh_combined[]) |
    (rast_full_combined[] < 2) |
    (rast_full_combined[] > max_lai_value_combined) |
    (rast_maxh_combined[] <= 10) |
    (rast_maxh_combined[] >= 40)
  na_mask_combined <- ref_mask_combined | is.na(rast_s2_combined)
  
  lid_combined <- mask(rast_opt_combined, na_mask_combined, maskvalue = TRUE)
  s2_combined <- mask(rast_s2_combined, na_mask_combined, maskvalue = TRUE)
  het_combined <- mask(rast_het_combined, na_mask_combined, maskvalue = TRUE)
  maxh_combined <- mask(rast_maxh_combined, na_mask_combined, maskvalue = TRUE)
  
  # Create the special "Sites combined" data frame for the current site
  df_site_special <- data.frame(
    lidar = as.numeric(values(lid_combined)),
    s2 = as.numeric(values(s2_combined)),
    mns_sd = as.numeric(values(het_combined)),
    height = as.numeric(values(maxh_combined))
  ) %>%
    drop_na() %>%
    mutate(
      site = site, # Temporarily keep the original site label
      s2_config = "LAI_S2_opt",
      lidar_config = "LAI_ALS_dopt"
    )
  special_data_list[[site]] <- df_site_special
}
agg_opt_dopt_special <- bind_rows(special_data_list) %>%
  mutate(site = "Sites combined")

# 3. Combine all datasets
raw_data <- bind_rows(
  raw_data_individual_sites,
  agg_atbd_als,
  agg_atbd_dopt,
  agg_opt_dopt_special
)

# 4. Compute metrics on the complete dataset
results <- compute_metrics(raw_data)
final_metrics <- results$metrics
data_for_metrics <- results$data_for_metrics

# stop()

final_metrics <- final_metrics %>%
  # filter(!(s2_config == "LAI_S2_opt" & lidar_config == "LAI_ALS")) %>%
  mutate(Relationship = paste(s2_config, "vs", lidar_config)) %>%
  dplyr::select(site, Relationship, mns_sd_class, R, RMSE, Bias, Slope, Slope_p_value)

final_metrics_R <- final_metrics %>%
  dplyr::select(site, Relationship, mns_sd_class, R)

final_metrics_RMSE <- final_metrics %>%
  dplyr::select(site, Relationship, mns_sd_class, RMSE)

final_metrics_Bias <- final_metrics %>%
  dplyr::select(site, Relationship, mns_sd_class, Bias)

final_metrics_Slope <- final_metrics %>%
  dplyr::select(site, Relationship, mns_sd_class, Slope)

plot_metric_by_het <- function(df, metric_col, y_label = NULL) {
  
  # Default y_label if not provided
  if (is.null(y_label)) y_label <- metric_col
  ylim <- switch(metric_col,
                 "R" = c(0, 0.9),
                 "RMSE" = c(0, 3.5),
                 "Slope" = c(0, 2),
                 "Bias" = c(-3.5, 2),
  )
  
  # Prepare shape group
  plot_data <- df %>%
    mutate(
      shape_group = ifelse(mns_sd_class == "Total", 
                           "Total", "By Heterogeneity Class")
    )
  
  # Shift All_Heterogeneity stars (3 relationships per site)
  stars_data <- plot_data %>%
    filter(shape_group == "Total") %>%
    group_by(site) %>%
    arrange(Relationship) %>%
    mutate(x_pos = c(1.7, 2, 2.3)) %>%
    ungroup()
  
  # Build plot
  p <- ggplot() +
    # Normal heterogeneity classes (points + lines)
    geom_point(
      data = plot_data %>% filter(shape_group != "Total"),
      aes(
        x = mns_sd_class, 
        y = .data[[metric_col]], 
        color = Relationship, 
        shape = shape_group, 
        group = Relationship
      ),
      size = 3
    ) +
    geom_line(
      data = plot_data %>% filter(shape_group != "Total"),
      aes(
        x = mns_sd_class, 
        y = .data[[metric_col]], 
        color = Relationship, 
        group = Relationship
      ),
      size = 1.2
    ) +
    
    # Stars for All_Heterogeneity
    geom_point(
      data = stars_data,
      aes(
        x = x_pos, 
        y = .data[[metric_col]], 
        color = Relationship, 
        shape = shape_group,
        group = Relationship
      ),
      size = 5, stroke = 2
    ) +
    
    facet_wrap(~ site) +
    scale_shape_manual(
      values = c("By Heterogeneity Class" = 16, "Total" = 8)
    ) +
    scale_color_discrete(
      name = "Relationship",
      labels = c(
        "LAI_S2_ATBD vs LAI_ALS"     = expression(LAI[S2_ATBD] ~ "vs" ~ LAI[ALS]),
        "LAI_S2_ATBD vs LAI_ALS_dopt"= expression(LAI[S2_ATBD] ~ "vs" ~ LAI[ALS_dopt]),
        "LAI_S2_opt vs LAI_ALS_dopt"  = expression(LAI[S2_opt] ~ "vs" ~ LAI[ALS_dopt])
      )
    ) +
    labs(
      x = expression(paste("Heterogeneity Class (", CHM[SD], ")")),
      y = y_label,
      shape = "Type"
    ) +
    guides(
      color = guide_legend(order = 1),  # Relationship legend first
      shape = guide_legend(order = 2)   # Type legend second, below
    ) +
    ylim(ylim) +  
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",    # stack legends vertically
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(filename = paste0(metric_col, ".png"),
         plot = p, width = 7, height = 7, dpi = 300)
}

# plot_metric_by_het(final_metrics_R, "R", "Pearson's R")
# plot_metric_by_het(final_metrics_RMSE, "RMSE")
# plot_metric_by_het(final_metrics_Bias, "Bias")
# plot_metric_by_het(final_metrics_Slope, "Slope")

# ------------------------------------------------------------------------------

plot_all_metrics_grid <- function(final_metrics) {
  
  # Gather into long format
  long_df <- final_metrics %>%
    tidyr::pivot_longer(
      cols = c(R, RMSE, Bias, Slope),
      names_to = "metric", values_to = "value"
    ) %>%
    dplyr::mutate(
      shape_group = mns_sd_class,
      slope_significant = ifelse(!is.na(Slope_p_value) & Slope_p_value < 0.05,
                                 "significant", "ns"),
      # Use *strings* that can later be parsed as expressions
      metric = factor(metric,
                      levels = c("R", "RMSE", "Bias", "Slope"),
                      labels = c("italic(r)", "RMSE", "Bias", "Slope"))
    )
  
  # Extract Total rows for star points
  stars_data <- long_df %>%
    dplyr::filter(shape_group == "Total") %>%
    dplyr::mutate(
      mns_sd_class = "Total",
      shape_group = "Total"
    )
  
  # Non-total points
  points_data <- long_df %>%
    dplyr::filter(shape_group != "Total") %>%
    dplyr::mutate(
      shape_group = "By Heterogeneity Class"
    )
  
  # Combine points and stars
  plot_data <- dplyr::bind_rows(points_data, stars_data)
  
  # Define a custom dark green to light green color gradient
  green_gradient_colors <- c(
    "LAI_S2_ATBD vs LAI_ALS"      = "#004d00", # Darkest green
    "LAI_S2_ATBD vs LAI_ALS_dopt" = "#008000", # Medium green
    "LAI_S2_opt vs LAI_ALS_dopt"  = "#4CBB17"  # Lightest green
  )
  
  # Reference lines: only for Bias (y=0) and Slope (y=1)
  ref_lines <- tibble::tibble(
    metric = factor(c("Bias", "Slope"), 
                    levels = c("italic(r)", "RMSE", "Bias", "Slope")),
    yintercept = c(0, 1)
  )
  
  # Custom Y limits for each metric
  ylim_data <- tibble::tibble(
    metric = factor(c("italic(r)", "RMSE", "Bias", "Slope"),
                    levels = c("italic(r)", "RMSE", "Bias", "Slope")),
    ymin = c(0, 0, NA, NA)
  )
  
  # Plot
  p <- ggplot2::ggplot() +
    # Reference lines
    ggplot2::geom_hline(
      data = ref_lines,
      ggplot2::aes(yintercept = yintercept),
      linetype = "dashed", color = "grey40"
    ) +
    # Invisible layer for setting lower y limits
    ggplot2::geom_blank(
      data = ylim_data,
      ggplot2::aes(y = ymin)
    ) +
    # Points and lines for Low/Medium/High
    # ggplot2::geom_point(
    #   data = plot_data %>% dplyr::filter(shape_group != "Total"),
    #   ggplot2::aes(
    #     x = mns_sd_class, y = value,
    #     color = Relationship, group = Relationship
    #   ),
    #   size = 2.5
    # ) +
    ggplot2::geom_point(
      data = plot_data %>% dplyr::filter(shape_group != "Total"),
      ggplot2::aes(
        x = mns_sd_class, y = value,
        color = Relationship, 
        group = Relationship,
        fill = dplyr::if_else(metric == "Slope" & slope_significant == "ns", "white", Relationship)
      ),
      shape = 21, size = 3, stroke = 1
    ) +
    ggplot2::geom_line(
      data = plot_data %>% dplyr::filter(shape_group != "Total"),
      ggplot2::aes(
        x = mns_sd_class, y = value,
        color = Relationship, group = Relationship
      ),
      size = 1
    ) +
    # Star points for Total
    ggplot2::geom_point(
      data = plot_data %>% dplyr::filter(shape_group == "Total"),
      ggplot2::aes(
        x = mns_sd_class, y = value,
        color = Relationship
      ),
      shape = 8, size = 5, stroke = 1.2
    ) +
    # ggplot2::facet_grid(metric ~ site, scales = "free_y") +
    ggplot2::facet_grid(metric ~ site, scales = "free_y", 
                        labeller = ggplot2::labeller(metric = ggplot2::label_parsed)) +
    ggplot2::scale_x_discrete(
      limits = c("Low", "Medium", "High", "Total")
    ) +
    # Use the manually defined green colors
    ggplot2::scale_color_manual(
      name = "Relationship",
      values = green_gradient_colors,
      labels = c(
        "LAI_S2_ATBD vs LAI_ALS"      = expression(LAI[S2_ATBD] ~ "vs" ~ LAI[ALS]),
        "LAI_S2_ATBD vs LAI_ALS_dopt" = expression(LAI[S2_ATBD] ~ "vs" ~ LAI[ALS_dopt]),
        "LAI_S2_opt vs LAI_ALS_dopt"  = expression(LAI[S2_opt] ~ "vs" ~ LAI[ALS_dopt])
      )
    ) +
    ggplot2::scale_fill_manual(
      name = "Relationship",
      values = c(green_gradient_colors, white = "white"),
      guide = "none"
    ) +
    ggplot2::labs(
      x = expression(paste("Horizontal Heterogeneity (", CHM[SD], ")")),
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey90")
    )
  
  # Save
  ggplot2::ggsave("All_Metrics_Grid.png", plot = p, width = 9, height = 9, dpi = 300)
  
  return(p)
}

# plot_all_metrics_grid(final_metrics)
# scatterplots pour illustrer
# lais, max, mns

stop()

hb <- hexbin(
  x = data_for_metrics$s2,
  y = data_for_metrics$lidar,
  xbins = 100
)
count_limits <- range(hb@count[hb@count > 0])

plot_hex_scatter <- function(data_for_metrics, site_name, count_limits) {
  
  # Filter for one site
  df_site <- data_for_metrics %>%
    filter(site == site_name)
  
  # Relationship labels for facet (plotmath kept for facet labels)
  relationship_labels <- c(
    "LAI_S2_ATBD.LAI_ALS"      = "LAI[S2_ATBD] * ' vs ' * LAI[ALS]",
    "LAI_S2_ATBD.LAI_ALS_dopt" = "LAI[S2_ATBD] * ' vs ' * LAI[ALS[dopt]]",
    "LAI_S2_opt.LAI_ALS_dopt"  = "LAI[S2_opt] * ' vs ' * LAI[ALS[dopt]]"
  )
  
  # Create column for facetting
  df_site <- df_site %>%
    mutate(
      relationship = paste(s2_config, lidar_config, sep = "."),
      relationship = factor(relationship, levels = names(relationship_labels)),
      mns_sd_class = factor(mns_sd_class, levels = c("Low", "Medium", "High", "Total"))
    )
  
  # Compute metrics per relationship × mns_sd_class
  metrics_df <- df_site %>%
    group_by(relationship, mns_sd_class) %>%
    summarise(
      R = round(cor(lidar, s2, use = "complete.obs"), 2),
      RMSE = round(sqrt(mean((lidar - s2)^2, na.rm = TRUE)), 2),
      Bias = round(mean(lidar - s2, na.rm = TRUE), 2),
      Slope = round(coef(lm(lidar ~ s2))[2], 2),
      .groups = "drop"
    ) %>%
    mutate(
      # HTML label: italic r, then line breaks with <br>
      label_html = paste0(
        "<i>r</i> = ", R, "<br/>",
        "RMSE = ", RMSE, "<br/>",
        "Bias = ", Bias, "<br/>",
        "Slope = ", Slope
      ),
      # adjust x/y to place the label nicely inside each panel
      x = 0.2,   # small positive to avoid clipping on left
      y = 15   # near the top (adjust to your y limits)
    )
  
  # Plot using geom_richtext for consistent sizes and HTML formatting
  p <- ggplot(df_site, aes(x = s2, y = lidar)) +
    geom_hex(bins = 100) +
    scale_fill_viridis_c(
      trans = "log10",
      limits = count_limits,
      oob = scales::squish
    ) + 
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    facet_grid(
      mns_sd_class ~ relationship,
      labeller = labeller(
        relationship = as_labeller(relationship_labels, default = label_parsed)
      )
    ) +
    # geom_richtext from ggtext: no parse needed; uses HTML tags and consistent font size
    geom_richtext(
      data = metrics_df,
      aes(x = x, y = y, label = label_html),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      lineheight = 0.2,
      size = 3.5,
      fill = NA,
      label.color = NA,
      label.padding = grid::unit(rep(0, 4), "pt")
    ) +
    labs(
      x = expression(LAI[S2]),
      y = expression(LAI[ALS]),
      fill = "Count (log scale)",
      title = site_name
    ) +
    xlim(c(0, 8)) +
    ylim(c(0, 15)) +
    theme_bw(base_size = 16) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  ggsave(paste0("hex_scatter_", site_name, ".png"), p, 
         width = 8, height = 8, dpi = 300)
  return(p)
}

plots <- lapply(unique(data_for_metrics$site), function(s) {
  plot_hex_scatter(data_for_metrics, s, count_limits)
})

names(plots) <- unique(data_for_metrics$site)
p_combined <- (plots[["Aigoual"]] | plots[["Blois"]]) /
  (plots[["Mormal"]] | plots[["Sites combined"]])
p_combined <- p_combined +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(
  "hex_scatter_all_sites_2x2.png",
  p_combined,
  width = 15,
  height = 15,
  dpi = 300
)
