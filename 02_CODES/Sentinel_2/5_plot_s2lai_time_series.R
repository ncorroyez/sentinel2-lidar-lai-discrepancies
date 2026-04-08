# ---
# title: "5_plot_s2lai_time_series.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-09-24"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# Import useful functions
library("readr")
library("terra")
library("ggplot2")
library("dplyr")
library("zoo")
source("../libraries/functions_plot_time_series.R")
source("../libraries/functions_general_tools.R")

# Pre-processing Parameters
# data_dir <- '/media/corroyez/MyPassport/01_DATA'
data_dir <- '../../01_DATA'
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- "Aigoual" # Mormal Blois Aigoual
results_dir <- "../../03_RESULTS"
output_dir <- "../../04_FIGURES/lai_time_series"

# Read S2 Acquisition Dates (1/month)
file_path <- file.path(data_dir, "S2_Dates.csv")
data <- read_csv(file_path)

# List of directories to keep
dirs_to_keep <- c(
  # "Full_Composition",
  # "Not_Masked",
  "Deciduous_Only"
  # "Deciduous_Flex",
  # "Coniferous_Only",
  # "Coniferous_Flex"
)

lai_data <- data.frame(
  site = character(),
  composition = character(),
  layer = character(),
  lai_mean = numeric(),
  lai_sd   = numeric()
)

# --- MAIN LOOP ---
for (site in sites) {
  
  cat("\nProcessing site:", site, "\n")
  
  metrics_dir <- file.path(results_dir, site, "Metrics")
  composition_dirs <- list.dirs(metrics_dir, recursive = FALSE)
  composition_dirs <- composition_dirs[grepl(paste(dirs_to_keep, collapse = "|"), composition_dirs)]
  
  for (composition_dir in composition_dirs) {
    
    composition <- basename(composition_dir)
    cat("  Composition:", composition, "\n")
    
    # Load LAI stack
    lai_stack <- terra::rast(
      file.path(composition_dir, "stack", paste0(site, "_lai_stacked_raster.tif"))
    )
    
    layer_names <- names(lai_stack)
    
    # Compute LAI stats for each layer
    for (i in 1:nlyr(lai_stack)) {
      
      r <- lai_stack[[i]]
      lai_vals <- values(r, na.rm = TRUE)
      
      lai_data <- rbind(
        lai_data,
        data.frame(
          site = site,
          composition = composition,
          layer = layer_names[i],
          lai_mean = mean(lai_vals, na.rm = TRUE),
          lai_sd   = sd(lai_vals, na.rm = TRUE)
        )
      )
    }
  }
}

# --- DATE PROCESSING ---
# Convert layer names YYYYMMDD → YYYY-MM
lai_data$layer <- as.yearmon(substr(lai_data$layer, 1, 6), "%Y%m")

# --- PLOT 1: Facet by site (per composition) ---
for (composition in unique(lai_data$composition)) {
  
  df <- lai_data %>% filter(composition == !!composition)
  
  p <- ggplot(df, aes(x = layer, y = lai_mean, group = site, color = site)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lai_mean - lai_sd, ymax = lai_mean + lai_sd), width = 0.1) +
    facet_wrap(~ site) +
    labs(
      title = paste("Temporal LAI -", composition),
      x = "Date",
      y = "LAI (Mean ± SD)"
    ) +
    scale_x_yearmon(
      breaks = sort(unique(df$layer)),
      format = "%m/%y"   # e.g., July 2021 -> 07/21
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
  
  ggsave(
    filename = file.path(output_dir, paste0("LAI_time_series_", composition, ".png")),
    plot = p,
    width = 10,
    height = 6
  )
}

# --- PLOT 2: Facet by composition (per site) ---
for (site in unique(lai_data$site)) {
  
  df <- lai_data %>% filter(site == !!site)
  
  p <- ggplot(df, aes(x = layer, y = lai_mean, group = composition, color = composition)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lai_mean - lai_sd, ymax = lai_mean + lai_sd), width = 0.1) +
    facet_wrap(~ composition) +
    labs(
      title = paste("Temporal LAI -", site),
      x = "Date",
      y = "LAI (Mean ± SD)"
    ) +
    scale_x_yearmon(
      breaks = sort(unique(df$layer)),
      format = "%m/%y"   # e.g., July 2021 -> 07/21
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
  
  ggsave(
    filename = file.path(output_dir, paste0("LAI_time_series_", site, ".png")),
    plot = p,
    width = 10,
    height = 6
  )
}
