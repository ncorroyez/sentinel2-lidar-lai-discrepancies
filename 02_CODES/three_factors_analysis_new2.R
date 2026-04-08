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
library(purrr)
library(relaimpo)
library(ranger)
library(rpart)
library(rpart.plot)

# --------------------------------- Setup --------------------------------------
results_dir <- "../03_RESULTS"
metrics_dir <- "Metrics/Deciduous_Only"
# lidar_lai_file <- "lidarlai_res_10_m.tif"
lidar_lai_file <- "ladstack_classic.tif"
max_file <- "max_res_10_m.tif"
mns_sd_file <- "dsm_sd_res_10_m.tif"
s2_files <- list(
  LAI_S2_ATBD = "s2lai_summer_atbd_res_10_m.tif",
  LAI_S2_opt = "s2lai_summer_best_indiv_res_10_m.tif"
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
    
    print(sum(!is.na(values(rast_opt))))
    print(sum(!is.na(values(rast_full))))
    
    # get global maximum LAI if you really need an explicit upper bound:
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
    
    for (cfg in names(s2_files)) {
      rast_s2 <- rast(file.path(base, s2_files[[cfg]]))
      
      for (depth_tag in c("LAI_ALS_dopt", "LAI_ALS")) {
        rast_lid <- if (depth_tag == "LAI_ALS_dopt") rast_opt else rast_full
        
        # add S2 into mask
        # print(ref_mask)
        # print(rast_s2)
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
          drop_na() # %>%
        # filter to lidar ≥2, height >10, and (optionally) ≤ max_lai_value
        # filter(
        #lidar >= 2,
        #height > 10,
        #lidar <= max_lai_value
        # )
        
        # skip if too few pixels
        # if (nrow(df) < 10) next
        
        # uniform sample up to 5000 points
        # if (nrow(df) > 50000) {
        #df <- df %>% slice_sample(n = 50000)
        # }
        
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
        # labels = c("1: < 2.5", "2: 2.5-5", "3: > 5"),
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
    dplyr::inner_join(sampled_ids_binned, by = c("site", "mns_sd_class", "row_id"))

  # --- Part 2: Prepare datasets for both analyses from the sampled data ---
  # Binned data, with the "Sites combined" site
  binned_data_for_metrics <- bind_rows(
    sampled_data,
    sampled_data %>% dplyr::mutate(site = "Sites combined")
  )
  
  # Full site data, with the "Sites combined" site, by wrapping the binned classes
  full_site_data_for_metrics <- bind_rows(
    sampled_data,
    sampled_data %>% dplyr::mutate(site = "Sites combined")
  ) %>%
    dplyr::mutate(mns_sd_class = "Total") %>%
    # Drop the `row_id` as it's no longer relevant for this grouping
    dplyr::select(-row_id)
  print(str(full_site_data_for_metrics))
  # --- Part 3: Calculate metrics for both analyses ---
  # Function to calculate metrics
  calculate_summary_metrics <- function(data_frame) {
    data_frame %>%
      dplyr::group_by(site, s2_config, lidar_config, mns_sd_class) %>%
      dplyr::filter(dplyr::n() >= 10) %>%
      dplyr::summarise(
        # mns_sd_center = mean(mns_sd, na.rm = TRUE),
        R = cor(s2, lidar, use = "complete.obs"),
        RMSE = sqrt(mean((s2 - lidar)^2, na.rm = TRUE)),
        Bias = mean(s2 - lidar, na.rm = TRUE),
        Slope = coef(lm(s2 ~ lidar))[2],
        .groups = "drop"
      ) %>%
      dplyr::ungroup()
  }
  
  # Calculate metrics for both datasets
  metrics_binned <- calculate_summary_metrics(binned_data_for_metrics)
  metrics_full <- calculate_summary_metrics(full_site_data_for_metrics)
  
  # --- Part 4: Combine and return ---
  a <- bind_rows(metrics_binned, metrics_full) %>%
    mutate(
      across(c(R, RMSE, Bias, Slope), ~round(., 2)),
      mns_sd_class = factor(
        mns_sd_class,
        levels = c("Low", "Medium", "High", "Total")
      )
    ) %>%
    arrange(site, lidar_config, s2_config)
  
  return(list(full_site_data_for_metrics, binned_data_for_metrics, a))
}

# -------------------------- 3. Main Execution ---------------------------------
raw_data <- analyze_three_factor_relationships()
final_metrics <- compute_metrics(raw_data)

stop()

final_metrics <- final_metrics %>%
  filter(!(s2_config == "LAI_S2_opt" & lidar_config == "LAI_ALS")) %>%
  mutate(Relationship = paste(s2_config, "vs", lidar_config)) %>%
  dplyr::select(site, Relationship, mns_sd_class, R, RMSE, Bias, Slope)

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

plot_metric_by_het(final_metrics_R, "R", "Pearson's R")
plot_metric_by_het(final_metrics_RMSE, "RMSE")
plot_metric_by_het(final_metrics_Bias, "Bias")
plot_metric_by_het(final_metrics_Slope, "Slope")

# ------------------------------------------------------------------------------

plot_all_metrics_grid <- function(final_metrics) {
  
  # Gather into long format
  long_df <- final_metrics %>%
    tidyr::pivot_longer(cols = c(R, RMSE, Bias, Slope),
                        names_to = "metric", values_to = "value") %>%
    dplyr::mutate(
      shape_group = mns_sd_class,
      metric = factor(metric, levels = c("R", "RMSE", "Bias", "Slope"))
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
    metric = factor(c("Bias", "Slope"), levels = c("R", "RMSE", "Bias", "Slope")),
    yintercept = c(0, 1)
  )
  
  # Plot
  p <- ggplot2::ggplot() +
    # Reference lines
    ggplot2::geom_hline(
      data = ref_lines,
      ggplot2::aes(yintercept = yintercept),
      linetype = "dashed", color = "grey40"
    ) +
    # Points and lines for Low/Medium/High
    ggplot2::geom_point(
      data = plot_data %>% dplyr::filter(shape_group != "Total"),
      ggplot2::aes(
        x = mns_sd_class, y = value,
        color = Relationship, group = Relationship
      ),
      size = 2.5
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
    ggplot2::facet_grid(metric ~ site, scales = "free_y") +
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
    ggplot2::labs(
      x = expression(paste("Horizontal Heterogeneity (", CHM[SD], ")")),
      y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey90")
    )
  
  # Save
  ggplot2::ggsave("All_Metrics_Grid.png", plot = p, width = 10, height = 10, dpi = 300)
  
  return(p)
}

plot_all_metrics_grid(final_metrics)
# scatterplots pour illustrer
# lais, max, mns
