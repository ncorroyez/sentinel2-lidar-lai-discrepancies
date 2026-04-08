# ---
# title: "0.validation_data.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-02-27"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

library("lidR")
library("data.table")
library("raster")
library("plotly")
library("terra")
library("viridis")
library("dplyr")
library("future")

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# Import useful functions
source("Sentinel_2/functions_sentinel_2.R")
source("functions_plots.R")
source("LiDAR/functions_lidar.R")

# Pre-processing Parameters & Directories
data_dir <- '../01_DATA/'
copernicus_dir <- file.path(data_dir, 'COPERNICUS_GBOV_RM7_20240227101352/RM07')
deciduous_dir <- file.path(data_dir, 
                           'COPERNICUS_GBOV_RM7_20240227101352/Deciduous_Forests')

# Create the Deciduous Forests directory if it doesn't exist
if (!dir.exists(deciduous_dir)) {
  dir.create(deciduous_dir)
}

# List all files in the RM07 directory
files <- list.files(copernicus_dir, full.names = TRUE)

# Filter only csv files
csv_files <- files[endsWith(files, ".csv")]

# Initialize an empty list to store dataframes
dfs <- list()

# Loop through each csv file
for (file in csv_files) {
  # Read the csv file
  data <- read.csv(file, sep = ";", header = TRUE, stringsAsFactors = FALSE)
  
  # Check if "IGBP_class" column exists
  if ("IGBP_class" %in% colnames(data)) {
    # Check if "Deciduous" word is present in the "IGBP_class" column
    if (any(grepl("Deciduous", data$IGBP_class))) {
      # Copy the csv file to Deciduous Forests directory
      # file.copy(file, file.path(deciduous_dir, basename(file)))
      
      # Filter out lines containing -999
      data <- data[!apply(data, 1, function(x) any(x == -999)), ]
      
      # Append dataframe to the list
      dfs[[length(dfs) + 1]] <- data
    }
  }
}

# Check if there are any dataframes in the list
if (length(dfs) > 0) {
  # Combine all dataframes into one
  combined_df <- bind_rows(dfs)
  
  # Export combined dataframe to CSV
  write.csv(combined_df, 
            file = file.path(deciduous_dir, "combined_deciduous_data.csv"),
            row.names = FALSE)
  
  # Optional: If you want to view the combined dataframe
  # print(combined_df)
} else {
  cat("No files with 'Deciduous' data found in the Copernicus directory.")
}

# Read the combined CSV file
combined_df <- read.csv(file.path(deciduous_dir, "combined_deciduous_data.csv"),
                        sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Filter lines where the date is in June
combined_df <- combined_df[substr(combined_df$TIME_IS, 5, 6) %in% c("06", "07"), ]

# Keep only unique PLOT_ID lines
combined_df <- combined_df %>%
  distinct(PLOT_ID, .keep_all = TRUE)

# Export the filtered dataframe to a new CSV file
write.csv(combined_df, 
          file = file.path(deciduous_dir, "june_july_unique_plot_data.csv"), 
          row.names = FALSE)