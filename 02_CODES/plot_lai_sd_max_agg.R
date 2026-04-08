# ---
# title: "plot_lai_sd_max_agg.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-12-02"
# ---

# ----------- (Optional) Clear the environment and free memory -----------------
rm(list=ls(all=TRUE)) 
gc() 

if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path)); getwd()
}

library("tidyr")
library("terra")
library("viridis")
library("dplyr")
library("ggplot2")
if (!require("patchwork")) install.packages("patchwork")
library("patchwork")
if (!require("scales")) install.packages("scales") # Needed for log breaks

# --------------------------------- Setup --------------------------------------
results_dir <- "../03_RESULTS"
figures_dir <- "../04_FIGURES"
metrics_dir <- "Metrics/Deciduous_Only"

if(!dir.exists(figures_dir)) dir.create(figures_dir)

sites <- c("Aigoual", "Blois", "Mormal")
lai_file <- "lidarlai_res_10_m.tif"
sd_file <- "dsm_sd_res_10_m.tif"
max_file <- "max_res_10_m.tif"

# --------------------------- Data Preparation ---------------------------------

prepare_data <- function(x_file, y_file){
  all_df <- list()
  for(site in sites){
    x_path <- file.path(results_dir, site, metrics_dir, x_file)
    y_path <- file.path(results_dir, site, metrics_dir, y_file)
    
    if(file.exists(x_path) & file.exists(y_path)){
      r_stack <- c(rast(x_path), rast(y_path))
      df <- as.data.frame(r_stack, xy = FALSE, na.rm = TRUE)
      colnames(df) <- c("X", "Y")
      df$Site <- site
      all_df[[site]] <- df
    }
  }
  bind_rows(all_df)
}

# --------------------------- Unified Scale Logic ------------------------------

combinations <- list(
  list(x = sd_file, y = max_file, xlab = expression(CHM[SD]), ylab = "Max Height (m)"),
  list(x = lai_file, y = max_file, xlab = expression(LAI[ALS]), ylab = "Max Height (m)"),
  list(x = sd_file, y = lai_file, xlab = expression(CHM[SD]), ylab = expression(LAI[ALS]))
)

# 1. Pre-calculate Global Max Count for unified legend
message("Calculating global max count...")
max_count_global <- 0

for(combo in combinations){
  df_temp <- prepare_data(combo$x, combo$y)
  p_temp <- ggplot(df_temp, aes(x=X, y=Y)) + geom_hex(bins=100)
  build <- ggplot_build(p_temp)
  current_max <- max(build$data[[1]]$count, na.rm=TRUE)
  max_count_global <- max(max_count_global, current_max)
}
message("Global max count set to: ", max_count_global)

# --------------------------- Plotting Loop ------------------------------------

plot_list <- list()

for(i in seq_along(combinations)){
  
  combo <- combinations[[i]]
  df <- prepare_data(combo$x, combo$y)
  
  p <- ggplot(df, aes(x = X, y = Y)) +
    geom_hex(bins = 100) + 
    geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 1.2) +
    
    # Unified Color Scale
    scale_fill_viridis(
      option = "viridis", 
      trans = "log", 
      name = "Count", 
      limits = c(1, max_count_global), 
      breaks = scales::breaks_log(n = 5)
    ) +
    
    facet_grid(~Site) + 
    labs(x = combo$xlab, y = combo$ylab) +
    
    theme_minimal(base_size = 20) + 
    theme(
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 18, color = "black"),
      strip.text = element_text(size = 24, face = "bold"),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 16),
      legend.key.height = unit(1.5, "cm"),
      
      # Ensure plot tag (a, b, c) is large and bold
      plot.tag = element_text(size = 30) 
    )
  
  # Remove Site headers for rows 2 and 3
  if(i != 1) {
    p <- p + theme(strip.text = element_blank())
  }
  
  plot_list[[i]] <- p
}

# --------------------------- Combine with Tags --------------------------------

final_plot <- (plot_list[[1]] / plot_list[[2]] / plot_list[[3]]) + 
  plot_layout(guides = "collect") + 
  
  # Add (a), (b), (c) tags
  plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)"))) &
  
  # Theme adjustments for final output
  theme(plot.tag.position = c(0.01, 1)) # Aligns tag to top-left corner

# --------------------------- Save ---------------------------------------------

ggsave(
  filename = file.path(figures_dir, "Combined_Hex_Grid_ABC_300dpi.png"), 
  plot = final_plot, 
  width = 14, 
  height = 14, 
  dpi = 300,
  bg = "white"
)

message("Saved with tags (a), (b), (c).")