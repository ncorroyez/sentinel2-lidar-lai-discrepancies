# ---
# title: "plot_lai_sd_max_sep.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-12-02"
# ---

# ----------- (Optional) Clear the environment and free memory -----------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --- Define working directory as the directory where the script is located ----
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# ------------------------------ Libraries -------------------------------------

library("tidyr")
library("plotly")
library("terra")
library("viridis")
library("dplyr")
library("ggplot2")
library("ggnewscale")
library("fields")

# --------------------------- Import useful functions --------------------------

source("libraries/functions_plots.R")

# --------------------------------- Setup --------------------------------------
results_dir <- "../03_RESULTS"
figures_dir <- "../04_FIGURES"
metrics_dir <- "Metrics/Deciduous_Only"
sites <- c("Aigoual", "Blois", "Mormal")
lai_file <- "lidarlai_res_10_m.tif"
sd_file <- "dsm_sd_res_10_m.tif"
max_file <- "max_res_10_m.tif"


# Function to prepare data for all sites
prepare_data <- function(x_file, y_file, xlab, ylab){
  all_df <- list()
  
  for(site in sites){
    x_r <- rast(file.path(results_dir, site, metrics_dir, x_file))
    y_r <- rast(file.path(results_dir, site, metrics_dir, y_file))
    
    # Ensure rasters are aligned
    rasters <- c(x_r, y_r)
    rasters <- crop(rasters, ext(rasters[[1]]))
    rasters <- mask(rasters, rasters[[1]])
    
    df <- as.data.frame(c(x_r, y_r), xy = FALSE, na.rm = TRUE)
    colnames(df) <- c("X", "Y")
    df$Site <- site
    
    all_df[[site]] <- df
  }
  
  bind_rows(all_df)
}

# List of combinations
combinations <- list(
  list(x = sd_file, y = max_file, xlab = expression(CHM[SD]), ylab = "Max Height (m)", fname = "1sd_max_hex.png"),
  list(x = lai_file, y = max_file, xlab = expression(LAI[ALS]), ylab = "Max Height (m)", fname = "2lai_max_hex.png"),
  list(x = sd_file, y = lai_file, xlab = expression(CHM[SD]), ylab = expression(LAI[ALS]), fname = "3sd_lai_hex.png")
)

# Loop through combinations and plot faceted by site
for(combo in combinations){
  df <- prepare_data(combo$x, combo$y, combo$xlab, combo$ylab)
  
  p <- ggplot(df, aes(x = X, y = Y)) +
    geom_hex(bins = 200) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    scale_fill_viridis(option = "viridis", trans = "log", name = "Count") +
    facet_wrap(~Site) +
    labs(x = combo$xlab, 
         y = combo$ylab
         # title = paste(combo$xlab, "vs", combo$ylab)
         ) +
    theme_minimal()
  
  ggsave(filename = file.path(figures_dir, combo$fname), plot = p, width = 12, height = 4)
}
