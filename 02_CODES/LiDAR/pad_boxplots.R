# ---
# title: "pad_boxplots.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2025-07-23"
# ---

# ----------- (Optional) Clear the environment and free memory -----------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --- Define working directory as the directory where the script is located ----
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# --------------------------------- Libraries ----------------------------------
library(terra)
library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)

make_pad_value_vs_depth_boxplot <- function(site,
                                            results_dir    = "../../03_RESULTS",
                                            metrics_subdir = "Metrics/Deciduous_Only/PAD_Profiles_dsm_keepTrees",
                                            # metrics_subdir = "Metrics/Deciduous_Only/PAD_Profiles_DTM_keepTrees",
                                            pad_pattern    = "^PAD_([0-9]+\\.?[0-9]*)_40\\.tif$",
                                            max_height     = 40) {
  
  # Build path & list files
  pad_dir   <- file.path(results_dir, site, metrics_subdir)
  pad_files <- list.files(pad_dir, pattern = pad_pattern, full.names = TRUE)
  
  # Read all rasters, pull PAD values, compute depth
  pad_df <- map_dfr(pad_files, function(f) {
    th       <- as.numeric(str_match(basename(f), pad_pattern)[,2])
    depth_m  <- max_height - th + 0.5
    vals     <- values(rast(f))
    
    tibble(
      PAD_value  = vals,
      Depth      = depth_m
    ) %>%
      filter(!is.na(PAD_value))
  })
  
  # Prepare reversed depth levels
  depth_levels <- sort(unique(pad_df$Depth), decreasing = TRUE)
  
  # Plot
  ggplot(pad_df, aes(
    x = PAD_value,
    y = factor(Depth, levels = depth_levels)
  )) +
    geom_boxplot(outlier.shape = NA) + 
    scale_y_discrete(drop = FALSE) +
    labs(
      title = paste0("PAD‑value distributions by depth — ", site),
      x     = "PAD value",
      y     = "Depth (m)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
}

# for (s in c("Aigoual")) {
#   print(make_pad_value_vs_depth_boxplot(s))
# }
# stop()

# ---- Example usage for all three sites ----
# for (s in c("Aigoual","Blois","Mormal")) {
#   print(make_pad_value_vs_depth_boxplot(s))
# }

# ------------------------------------------------------------------------------
make_all_pad_value_vs_depth_boxplot_dual <- function(sites = c("Aigoual", "Blois", "Mormal"),
                                                     results_dir = "../../03_RESULTS",
                                                     metrics_subdirs = list(
                                                       CHM = "Metrics/Deciduous_Only/PAD_Profiles_DSM_keepTrees2",
                                                       DTM = "Metrics/Deciduous_Only/PAD_Profiles_DTM_keepTrees2"
                                                     ),
                                                     pad_pattern = "^PAD_([0-9]+\\.?[0-9]*)_40\\.tif$",
                                                     max_height = 40) {
  # Helper to read PAD values for one site and normalization
  read_pad_data <- function(site, subdir, norm_type) {
    pad_dir <- file.path(results_dir, site, subdir)
    pad_files <- list.files(pad_dir, pattern = pad_pattern, full.names = TRUE)
    
    map_dfr(pad_files, function(f) {
      th <- as.numeric(str_match(basename(f), pad_pattern)[, 2])
      depth_m <- max_height - th + 0.5
      vals <- values(rast(f))
      
      tibble(
        PAD_value = vals,
        Depth = depth_m,
        NormType = norm_type,
        Site = site
      ) %>% filter(!is.na(PAD_value))
    })
  }
  
  # Read data for all sites and both normalization types
  all_pad_df <- map_dfr(sites, function(site) {
    bind_rows(
      read_pad_data(site, metrics_subdirs$CHM, "CHM"),
      read_pad_data(site, metrics_subdirs$DTM, "DTM")
    )
  })
  
  # Set reversed depth levels
  depth_levels <- sort(unique(all_pad_df$Depth), decreasing = TRUE)
  
  # Plot
  ggplot(all_pad_df, aes(
    x = PAD_value,
    y = factor(Depth, levels = depth_levels),
    fill = NormType
  )) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
    facet_wrap(~Site, scales = "fixed") +
    scale_y_discrete(drop = FALSE) +
    scale_fill_discrete(
      labels = c(
        expression(CHM[ALS]),
        expression(DTM[ALS])
      )
    ) +
    labs(
      # title = "PAD‑value distributions by depth and normalization type",
      x = "PAD value",
      y = "Depth (m)",
      fill = "Normalization"
    ) +
    xlim(c(0, 15)) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.text = element_text(face = "bold", size = 13)
    )
}
pad_plot <- make_all_pad_value_vs_depth_boxplot_dual()
# print(pad_plot)
ggsave(
  filename = "PAD_boxplot_dual_sites.png",
  plot = pad_plot,
  width = 8,
  height = 8,
  dpi = 300
)
