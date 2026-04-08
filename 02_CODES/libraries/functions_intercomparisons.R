# ---
# title: "functions_intercomparisons.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-02-06"
# ---

library(terra)
library(ggplot2)
library(reshape2)
library(cowplot)
library(dplyr)
library(hrbrthemes)
library(raster)
library(rasterVis)
library(latticeExtra)
library(RColorBrewer)

#' Create Data Frame
#'
#' This function creates a data frame with specified row and column names.
#'
#' @param bands A vector of band names.
#' @param noises A vector of noise types.
#' @param row_names A vector of row names.
#' @param col_names A vector of column names.
#'
#' @return A data frame with specified row and column names.
#'
#' @examples
#' bands <- c("B1", "B2")
#' noises <- c("N1", "N2")
#' row_names <- c("Row1", "Row2")
#' col_names <- c("Col1", "Col2")
#' create_dataframe(bands, noises, row_names, col_names)
#'
#' @export
create_dataframe <- function(bands, 
                             noises,
                             row_names,
                             col_names){
  
  data_matrix <- matrix(NA, nrow = length(row_names), ncol = length(col_names))
  df <- data.frame(data_matrix)
  
  rownames(df) <- paste(row_names, col_names, sep = "_")
  colnames(df) <- paste(row_names, col_names, sep = "_")
  return(df)
}

#' Create List of Variables for ATBD Bands and Noises
#'
#' This function creates a list of data matrices for different bands and noises.
#'
#' @param bands A vector of band names.
#' @param noises A vector of noise types.
#' @param mask A spatial mask for the data.
#' @param plots_dir Directory for storing plots.
#' @param s2_output_dir Directory containing S2 output data.
#' @param image The image identifier.
#' @param estimated_var The variable being estimated.
#'
#' @return A list of data matrices for different bands and noises.
#'
#' @examples
#' bands <- c("B1", "B2")
#' noises <- c("N1", "N2")
#' mask <- terra::rast("path/to/mask")
#' create_atbd_bands_noises_vars_list(bands, noises, mask, "plots_dir", 
#' "s2_output_dir", "image", "estimated_var")
#'
#' @export
create_atbd_bands_noises_vars_list <- function(bands,
                                               noises,
                                               mask,
                                               plots_dir,
                                               s2_output_dir,
                                               image,
                                               estimated_var){
  
  # Create a list to store data matrices
  data_list <- list()
  
  # Load your data into the list, assuming you have them named as: 
  # '10m_addmult', '20m_addmult', etc.
  for (noise in noises) {
    for (band in bands) {
      data_path <- file.path(s2_output_dir,
                             'PRO4SAIL_INVERSION_atbd',
                             noise,
                             band,
                             paste0(image, 
                                    '_', 
                                    'Refl',
                                    '_',
                                    estimated_var, 
                                    ".envi"))
      data_name <- paste0(noise, "_", band)
      data <- terra::rast(data_path) 
      data <- terra::project(data, mask)
      data <- data * mask
      data <- values(data)
      data_list[[data_name]] <- data
    }
  }
  return(data_list)
}

#' Create List of Rasters for ATBD Bands and Noises
#'
#' This function creates a list of raster data for different bands and noises.
#'
#' @param bands A vector of band names.
#' @param noises A vector of noise types.
#' @param mask A spatial mask for the data.
#' @param plots_dir Directory for storing plots.
#' @param s2_output_dir Directory containing S2 output data.
#' @param image The image identifier.
#' @param estimated_var The variable being estimated.
#'
#' @return A list of raster data for different bands and noises.
#'
#' @examples
#' bands <- c("B1", "B2")
#' noises <- c("N1", "N2")
#' mask <- terra::rast("path/to/mask")
#' create_atbd_bands_noises_rasters_list(bands, noises, mask, "plots_dir", 
#' "s2_output_dir", "image", "estimated_var")
#'
#' @export
create_atbd_bands_noises_rasters_list <- function(bands,
                                                  noises,
                                                  mask,
                                                  plots_dir,
                                                  s2_output_dir,
                                                  image,
                                                  estimated_var){
  
  # Create a list to store data matrices
  data_list <- list()
  
  # Load your data into the list, assuming you have them named as: 
  # '10m_addmult', '20m_addmult', etc.
  for (noise in noises) {
    for (band in bands) {
      data_path <- file.path(s2_output_dir,
                             'PRO4SAIL_INVERSION_atbd',
                             noise,
                             band,
                             paste0(image, 
                                    '_', 
                                    'Refl',
                                    '_',
                                    estimated_var, 
                                    ".envi"))
      data_name <- paste0(noise, "_", band)
      data <- terra::rast(data_path)
      data <- terra::project(data, mask)
      data <- data * mask
      data_list[[data_name]] <- data
    }
  }
  return(data_list)
}

#' Save Data Frame to CSV
#'
#' This function saves a data frame to a CSV file.
#'
#' @param df The data frame to save.
#' @param df_path_to_save The path where the data frame should be saved.
#' @param df_name_to_save The name of the CSV file to save.
#'
#' @return None. The data frame is saved as a CSV file.
#'
#' @examples
#' df <- data.frame(matrix(rnorm(10), nrow = 5))
#' save_df(df, "path/to/save", "my_dataframe")
#'
#' @export
save_df <- function(df, df_path_to_save, df_name_to_save){
  if (is.null(df_path_to_save) || is.null(df_name_to_save)) {
    stop("DF Path or Name cannot be NULL.")
  }
  write.csv(format(df, digits = 3), 
            file.path(df_path_to_save, 
                      paste0(df_name_to_save, ".csv")),
            row.names = TRUE)
  cat("DF has been successfully saved at", 
      paste0(df_name_to_save, ".csv"), "\n")
}

#' Correlation Analysis for ATBD Bands and Noises
#'
#' This function performs a correlation analysis for different bands and noises,
#'  and saves the result as a CSV file.
#'
#' @param bands A vector of band names.
#' @param noises A vector of noise types.
#' @param row_names A vector of row names.
#' @param col_names A vector of column names.
#' @param mask A spatial mask for the data.
#' @param plots_dir Directory for storing plots.
#' @param s2_output_dir Directory containing S2 output data.
#' @param image The image identifier.
#' @param estimated_var The variable being estimated.
#' @param comparison_path The path where the comparison result should be saved.
#' @param csv_name The name of the CSV file to save the comparison result.
#'
#' @return None. The correlation analysis result is saved as a CSV file.
#'
#' @examples
#' bands <- c("B1", "B2")
#' noises <- c("N1", "N2")
#' row_names <- c("Row1", "Row2")
#' col_names <- c("Col1", "Col2")
#' mask <- terra::rast("path/to/mask")
#' correlation_atbd_bands_noise(bands, noises, row_names, col_names, mask,
#' "plots_dir", "s2_output_dir", "image", "estimated_var", "comparison_path", 
#' "csv_name")
#'
#' @export
correlation_atbd_bands_noise <- function(bands,
                                         noises,
                                         row_names,
                                         col_names,
                                         mask,
                                         plots_dir,
                                         s2_output_dir,
                                         image,
                                         estimated_var,
                                         comparison_path,
                                         csv_name){
  
  df <- create_dataframe(bands, noises, row_names, col_names)
  
  data_list <- create_atbd_bands_noises_vars_list(bands,
                                                  noises,
                                                  mask,
                                                  plots_dir,
                                                  s2_output_dir,
                                                  image,
                                                  estimated_var)
  
  cor_matrix <- cor(do.call(cbind, data_list), use = "pairwise.complete.obs")
  
  for (i in 1:length(row_names)) {
    for (j in 1:length(col_names)) {
      df[i, j] <- cor_matrix[i, j]
    }
  }
  save_df(df, 
          comparison_path,
          csv_name)
  df <- df[2:3, c(1, 4)]
  save_df(df, 
          comparison_path, 
          paste0(csv_name, "_subset"))
}

#' RMSE Analysis for ATBD Bands and Noises
#'
#' This function performs an RMSE analysis for different bands and noises, and saves the result as a CSV file.
#'
#' @param bands A vector of band names.
#' @param noises A vector of noise types.
#' @param row_names A vector of row names.
#' @param col_names A vector of column names.
#' @param mask A spatial mask for the data.
#' @param plots_dir Directory for storing plots.
#' @param s2_output_dir Directory containing S2 output data.
#' @param image The image identifier.
#' @param estimated_var The variable being estimated.
#' @param comparison_path The path where the comparison result should be saved.
#' @param csv_name The name of the CSV file to save the comparison result.
#'
#' @return None. The RMSE analysis result is saved as a CSV file.
#'
#' @examples
#' bands <- c("B1", "B2")
#' noises <- c("N1", "N2")
#' row_names <- c("Row1", "Row2")
#' col_names <- c("Col1", "Col2")
#' mask <- terra::rast("path/to/mask")
#' rmse_atbd_bands_noise(bands, noises, row_names, col_names, mask, "plots_dir", "s2_output_dir", "image", "estimated_var", "comparison_path", "csv_name")
#'
#' @export
rmse_atbd_bands_noise <- function(bands,
                                  noises,
                                  row_names,
                                  col_names,
                                  mask,
                                  plots_dir,
                                  s2_output_dir,
                                  image,
                                  estimated_var,
                                  comparison_path,
                                  csv_name){
  
  df <- create_dataframe(bands, noises, row_names, col_names)
  
  data_list <- create_atbd_bands_noises_vars_list(bands,
                                                  noises,
                                                  mask,
                                                  plots_dir,
                                                  s2_output_dir,
                                                  image,
                                                  estimated_var)
  
  for(i in seq_along(row_names)){
    for(j in seq_along(col_names)){
      var1 <- data_list[[i]]  # Use double brackets to access list elements
      var1 <- var1[!is.nan(var1)]
      var2 <- data_list[[j]]  # Use double brackets to access list elements
      var2 <- var2[!is.nan(var2)]
      rmse <- Metrics::rmse(var1, var2)
      df[i, j] <- rmse
    }
  }
  save_df(df, 
          comparison_path,
          csv_name)
  df <- df[2:3, c(1, 4)]
  save_df(df, 
          comparison_path, 
          paste0(csv_name, "_subset"))
}

#' Bias Analysis for ATBD Bands and Noises
#'
#' This function performs a bias analysis for different bands and noises, 
#' and saves the result as a CSV file.
#'
#' @param bands A vector of band names.
#' @param noises A vector of noise types.
#' @param row_names A vector of row names.
#' @param col_names A vector of column names.
#' @param mask A spatial mask for the data.
#' @param plots_dir Directory for storing plots.
#' @param s2_output_dir Directory containing S2 output data.
#' @param image The image identifier.
#' @param estimated_var The variable being estimated.
#' @param comparison_path The path where the comparison result should be saved.
#' @param csv_name The name of the CSV file to save the comparison result.
#'
#' @return None. The bias analysis result is saved as a CSV file.
#'
#' @examples
#' bands <- c("B1", "B2")
#' noises <- c("N1", "N2")
#' row_names <- c("Row1", "Row2")
#' col_names <- c("Col1", "Col2")
#' mask <- terra::rast("path/to/mask")
#' bias_atbd_bands_noise(bands, noises, row_names, col_names, mask, "plots_dir",
#'  "s2_output_dir", "image", "estimated_var", "comparison_path", "csv_name")
#'
#' @export
bias_atbd_bands_noise <- function(bands,
                                  noises,
                                  row_names,
                                  col_names,
                                  mask,
                                  plots_dir,
                                  s2_output_dir,
                                  image,
                                  estimated_var,
                                  comparison_path,
                                  csv_name){
  
  df <- create_dataframe(bands, noises, row_names, col_names)
  
  data_list <- create_atbd_bands_noises_vars_list(bands,
                                                  noises,
                                                  mask,
                                                  plots_dir,
                                                  s2_output_dir,
                                                  image,
                                                  estimated_var)
  
  for(i in seq_along(row_names)){
    for(j in seq_along(col_names)){
      var1 <- data_list[[i]]  # Use double brackets to access list elements
      var1 <- var1[!is.nan(var1)]
      var2 <- data_list[[j]]  # Use double brackets to access list elements
      var2 <- var2[!is.nan(var2)]
      bias <- Metrics::bias(var1, var2)
      df[i, j] <- bias
    }
  }
  save_df(df, 
          comparison_path,
          csv_name)
  df <- df[2:3, c(1, 4)]
  save_df(df, 
          comparison_path, 
          paste0(csv_name, "_subset"))
}

#' Residual Analysis for ATBD Bands and Noises
#'
#' This function calculates residuals between different bands and noises, 
#' and saves the residuals as raster files.
#'
#' @param bands A vector of band names.
#' @param noises A vector of noise types.
#' @param row_names A vector of row names.
#' @param col_names A vector of column names.
#' @param mask A spatial mask for the data.
#' @param plots_dir Directory for storing plots.
#' @param s2_output_dir Directory containing S2 output data.
#' @param image The image identifier.
#' @param estimated_var The variable being estimated.
#' @param comparison_path The path where the residuals should be saved.
#'
#' @return None. The residuals are saved as raster files.
#'
#' @examples
#' bands <- c("B1", "B2")
#' noises <- c("N1", "N2")
#' row_names <- c("Row1", "Row2")
#' col_names <- c("Col1", "Col2")
#' mask <- terra::rast("path/to/mask")
#' residuals_atbd_bands_noise(bands, noises, row_names, col_names, mask, 
#' "plots_dir", "s2_output_dir", "image", "estimated_var", "comparison_path")
#'
#' @export
residuals_atbd_bands_noise <- function(bands,
                                       noises,
                                       row_names,
                                       col_names,
                                       mask,
                                       plots_dir,
                                       s2_output_dir,
                                       image,
                                       estimated_var,
                                       comparison_path){
  
  df <- create_dataframe(bands, noises, row_names, col_names)
  
  data_list <- create_atbd_bands_noises_rasters_list(bands,
                                                     noises,
                                                     mask,
                                                     plots_dir,
                                                     s2_output_dir,
                                                     image,
                                                     estimated_var)
  
  for(i in seq_along(row_names)){
    for(j in seq_along(col_names)){
      if(!identical(data_list[[i]], data_list[[j]])){
        var1 <- data_list[[i]]  # Use double brackets to access list elements
        var2 <- data_list[[j]]  # Use double brackets to access list elements
        residuals <- var1 - var2
        save_envi_file(residuals,
                       paste0(rownames(df)[i], 
                              "_", 
                              colnames(df)[j], 
                              "_residuals"),
                       file.path(comparison_path, 
                                 "residuals",
                                 paste0(rownames(df)[i], 
                                        "_", 
                                        colnames(df)[j])))
      }
    }
  }
}

#' Error Analysis for ATBD Bands and Noises
#'
#' This function calculates the mean absolute error for different bands 
#' and noises, and saves the result as a CSV file.
#'
#' @param bands A vector of band names.
#' @param noises A vector of noise types.
#' @param row_names A vector of row names.
#' @param col_names A vector of column names.
#' @param mask A spatial mask for the data.
#' @param plots_dir Directory for storing plots.
#' @param s2_output_dir Directory containing S2 output data.
#' @param image The image identifier.
#' @param estimated_var The variable being estimated.
#' @param comparison_path The path where the comparison result should be saved.
#' @param csv_name The name of the CSV file to save the comparison result.
#'
#' @return None. The error analysis result is saved as a CSV file.
#'
#' @examples
#' bands <- c("B1", "B2")
#' noises <- c("N1", "N2")
#' row_names <- c("Row1", "Row2")
#' col_names <- c("Col1", "Col2")
#' mask <- terra::rast("path/to/mask")
#' error_atbd_bands_noise(bands, noises, row_names, col_names, mask,
#' "plots_dir", "s2_output_dir", "image", "estimated_var", "comparison_path",
#' "csv_name")
#'
#' @export
error_atbd_bands_noise <- function(bands,
                                   noises,
                                   row_names,
                                   col_names,
                                   mask,
                                   plots_dir,
                                   s2_output_dir,
                                   image,
                                   estimated_var,
                                   comparison_path,
                                   csv_name){
  
  df <- create_dataframe(bands, noises, row_names, col_names)
  
  data_list <- create_atbd_bands_noises_vars_list(bands,
                                                  noises,
                                                  mask,
                                                  plots_dir,
                                                  s2_output_dir,
                                                  image,
                                                  estimated_var)
  
  for(i in seq_along(row_names)){
    for(j in seq_along(col_names)){
      var1 <- data_list[[i]]  # Use double brackets to access list elements
      var1 <- var1[!is.nan(var1)]
      var2 <- data_list[[j]]  # Use double brackets to access list elements
      var2 <- var2[!is.nan(var2)]
      residuals <- var1 - var2
      error <- mean(abs(residuals)) 
      df[i, j] <- error
    }
  }
  save_df(df, 
          comparison_path,
          csv_name)
  df <- df[2:3, c(1, 4)]
  save_df(df, 
          comparison_path, 
          paste0(csv_name, "_subset"))
}

#' R-Squared Analysis for ATBD Bands and Noises
#'
#' This function calculates the R-squared values for different bands and noises,
#'  and saves the result as a CSV file.
#'
#' @param bands A vector of band names.
#' @param noises A vector of noise types.
#' @param row_names A vector of row names.
#' @param col_names A vector of column names.
#' @param mask A spatial mask for the data.
#' @param plots_dir Directory for storing plots.
#' @param s2_output_dir Directory containing S2 output data.
#' @param image The image identifier.
#' @param estimated_var The variable being estimated.
#' @param comparison_path The path where the comparison result should be saved.
#' @param csv_name The name of the CSV file to save the comparison result.
#'
#' @return None. The R-squared analysis result is saved as a CSV file.
#'
#' @examples
#' bands <- c("B1", "B2")
#' noises <- c("N1", "N2")
#' row_names <- c("Row1", "Row2")
#' col_names <- c("Col1", "Col2")
#' mask <- terra::rast("path/to/mask")
#' r_squared_atbd_bands_noise(bands, noises, row_names, col_names, mask, 
#' "plots_dir", "s2_output_dir", "image", "estimated_var", "comparison_path", 
#' "csv_name")
#'
#' @export
r_squared_atbd_bands_noise <- function(bands,
                                       noises,
                                       row_names,
                                       col_names,
                                       mask,
                                       plots_dir,
                                       s2_output_dir,
                                       image,
                                       estimated_var,
                                       comparison_path,
                                       csv_name){
  
  df <- create_dataframe(bands, noises, row_names, col_names)
  
  data_list <- create_atbd_bands_noises_vars_list(bands,
                                                  noises,
                                                  mask,
                                                  plots_dir,
                                                  s2_output_dir,
                                                  image,
                                                  estimated_var)
  
  for(i in seq_along(row_names)){
    for(j in seq_along(col_names)){
      var1 <- data_list[[i]]  # Use double brackets to access list elements
      var1 <- var1[!is.nan(var1)]
      var2 <- data_list[[j]]  # Use double brackets to access list elements
      var2 <- var2[!is.nan(var2)]
      r_squared <- cor(var1, var2)^2
      df[i, j] <- r_squared
    }
  }
  save_df(df, 
          comparison_path,
          csv_name)
  df <- df[2:3, c(1, 4)]
  save_df(df, 
          comparison_path, 
          paste0(csv_name, "_subset"))
}

#' Standard Deviation Analysis for ATBD Bands and Noises
#'
#' This function calculates the standard deviation of residuals for different 
#' bands and noises, and saves the result as a CSV file.
#'
#' @param bands A vector of band names.
#' @param noises A vector of noise types.
#' @param row_names A vector of row names.
#' @param col_names A vector of column names.
#' @param mask A spatial mask for the data.
#' @param plots_dir Directory for storing plots.
#' @param s2_output_dir Directory containing S2 output data.
#' @param image The image identifier.
#' @param estimated_var The variable being estimated.
#' @param comparison_path The path where the comparison result should be saved.
#' @param csv_name The name of the CSV file to save the comparison result.
#'
#' @return None. The standard deviation analysis result is saved as a CSV file.
#'
#' @examples
#' bands <- c("B1", "B2")
#' noises <- c("N1", "N2")
#' row_names <- c("Row1", "Row2")
#' col_names <- c("Col1", "Col2")
#' mask <- terra::rast("path/to/mask")
#' std_atbd_bands_noise(bands, noises, row_names, col_names, mask, "plots_dir",
#' "s2_output_dir", "image", "estimated_var", "comparison_path", "csv_name")
#'
#' @export
std_atbd_bands_noise <- function(bands,
                                 noises,
                                 row_names,
                                 col_names,
                                 mask,
                                 plots_dir,
                                 s2_output_dir,
                                 image,
                                 estimated_var,
                                 comparison_path,
                                 csv_name){
  
  df <- create_dataframe(bands, noises, row_names, col_names)
  
  data_list <- create_atbd_bands_noises_vars_list(bands,
                                                  noises,
                                                  mask,
                                                  plots_dir,
                                                  s2_output_dir,
                                                  image,
                                                  estimated_var)
  
  for(i in seq_along(row_names)){
    for(j in seq_along(col_names)){
      var1 <- data_list[[i]]  # Use double brackets to access list elements
      var1 <- var1[!is.nan(var1)]
      var2 <- data_list[[j]]  # Use double brackets to access list elements
      var2 <- var2[!is.nan(var2)]
      residuals <- var1 - var2
      std <- sd(residuals)
      df[i, j] <- std
    }
  }
  save_df(df, 
          comparison_path,
          csv_name)
  df <- df[2:3, c(1, 4)]
  save_df(df, 
          comparison_path, 
          paste0(csv_name, "_subset"))
}

#' Histogram and Scatterplot Analysis for ATBD Bands and Noises
#'
#' This function generates histograms and scatterplots for different bands 
#' and noises, and saves the plots in the specified directory.
#'
#' @param bands A vector of band names.
#' @param noises A vector of noise types.
#' @param row_names A vector of row names.
#' @param col_names A vector of column names.
#' @param mask A spatial mask for the data.
#' @param plots_dir Directory for storing plots.
#' @param s2_output_dir Directory containing S2 output data.
#' @param image The image identifier.
#' @param estimated_var The variable being estimated.
#' @param comparison_path The path where the plots should be saved.
#'
#' @return None. The histograms and scatterplots are saved in 
#' the specified directory.
#'
#' @examples
#' bands <- c("B1", "B2")
#' noises <- c("N1", "N2")
#' row_names <- c("Row1", "Row2")
#' col_names <- c("Col1", "Col2")
#' mask <- terra::rast("path/to/mask")
#' histograms_scatterplots_atbd_bands_noise(bands, noises, row_names, 
#' col_names, mask, "plots_dir", "s2_output_dir", "image", "estimated_var", 
#' "comparison_path")
#'
#' @export
histograms_scatterplots_atbd_bands_noise <- function(bands,
                                                     noises,
                                                     row_names,
                                                     col_names,
                                                     mask,
                                                     plots_dir,
                                                     s2_output_dir,
                                                     image,
                                                     estimated_var,
                                                     comparison_path){
  
  df <- create_dataframe(bands, noises, row_names, col_names)
  
  data_list <- create_atbd_bands_noises_vars_list(bands,
                                                  noises,
                                                  mask,
                                                  plots_dir,
                                                  s2_output_dir,
                                                  image,
                                                  estimated_var)
  
  for(i in seq_along(row_names)){
    for(j in seq_along(col_names)){
      if(!identical(data_list[[i]], data_list[[j]])){
        if((names(data_list)[i] == "addmult_10m" 
            && names(data_list)[j] == "addmult_20m")
           || (names(data_list)[i] == "mult_10m" 
               && names(data_list)[j] == "mult_20m")
           || (names(data_list)[i] == "addmult_10m" 
               && names(data_list)[j] == "mult_10m")
           || (names(data_list)[i] == "addmult_20m" 
               && names(data_list)[j] == "mult_20m")
        ) {
          residuals <- data_list[[i]] - data_list[[j]]
          vars_list <- list(data_list[[i]], data_list[[j]], residuals)
          var_labs <- c(names(data_list)[i], names(data_list)[j], "residuals")
          plot_histogram(vars_list, 
                         "LAI Distributions", 
                         "LAI",
                         var_labs, 
                         comparison_path, 
                         paste0(
                           "hist",
                           "_",
                           names(data_list)[i],
                           "_",
                           names(data_list)[j],
                           "_",
                           "residuals"
                         ))
          plot_density_scatterplot(var_x = data_list[[i]],
                                   var_y = data_list[[j]],
                                   xlab = names(data_list)[j],
                                   ylab = names(data_list)[i],
                                   dirname = comparison_path, 
                                   filename = paste0(
                                     "scatterplot",
                                     "_",
                                     names(data_list)[i],
                                     "_",
                                     names(data_list)[j],
                                     "_"))
        }
      }
    }
  }
}

#' Calculate Metrics for Actual and Predicted Values
#'
#' This function calculates various metrics (correlation, R-squared, RMSE,
#' bias, error, standard deviation)
#' between actual and predicted values, generates histogram and scatterplot
#' visualizations, and saves the results.
#'
#' @param actual RasterLayer of actual values.
#' @param predicted RasterLayer of predicted values.
#' @param save_path Directory to save the results.
#' @param actual_lab Label for the actual values. Default is "actual".
#' @param predicted_lab Label for the predicted values. Default is "predicted".
#'
#' @return A data frame containing the calculated metrics.
#'
#' @examples
#' actual <- terra::rast("path/to/actual")
#' predicted <- terra::rast("path/to/predicted")
#' save_path <- "path/to/save/results"
#' calculate_metrics(actual, predicted, save_path)
#'
#' @export
calculate_metrics <- function(actual, 
                              predicted, 
                              save_path,
                              actual_lab = NULL,
                              predicted_lab = NULL
) {
  if (is.null(actual_lab)){
    actual_lab <- "actual"
  }
  if (is.null(predicted_lab)){
    predicted_lab <- "predicted"
  }
  
  # Take values
  actual_values <- values(actual)
  predicted_values <- values(predicted)
  
  # Residuals
  residuals_raster <- actual - predicted
  residuals <- values(residuals_raster)
  
  # Check if the lengths of actual and predicted are the same
  if (length(actual) != length(predicted)) {
    stop("Lengths of actual and predicted values must be the same.")
  }
  
  hist_list <- list(actual_values, predicted_values, residuals)
  
  # Plot histogram: actual, predicted, residuals
  plot_histogram(vars_list = hist_list, 
                 title = "LAI Distributions", 
                 xlab = "LAI", 
                 var_labs = sprintf(c(actual_lab, predicted_lab, "Residuals")), 
                 dirname = save_path, 
                 filename = sprintf("hist_%s_%s_residuals",
                                    actual_lab,
                                    predicted_lab))
  
  plot_density_scatterplot(var_x = predicted_values,
                           var_y = actual_values,
                           xlab = sprintf("%s", predicted_lab), 
                           ylab = sprintf("%s", actual_lab),
                           dirname = save_path, 
                           filename = sprintf("scatterplot_%s_%s",
                                              actual_lab,
                                              predicted_lab),
                           xlimits = c(0, max(predicted_values)),
                           ylimits = c(0, max(actual_values)))
  
  # Remove NAs from both actual and predicted
  valid_indices <- complete.cases(actual_values, predicted_values)
  actual_values <- actual_values[valid_indices]
  predicted_values <- predicted_values[valid_indices]
  residuals <- residuals[valid_indices]
  
  # Check if there are remaining values after removing NAs
  if (length(actual) == 0) {
    stop("No valid data points remaining after removing NAs.")
  }
  
  # Calculate metrics
  correlation <- cor(actual_values, predicted_values)
  ssres <- sum((residuals)^2)
  sstot <- sum((actual_values - mean(actual_values))^2)
  print(ssres)
  print(sstot)
  r_squared <- 1 - (ssres / sstot)
  rmse <- sqrt(mean(residuals^2))
  bias <- mean(residuals)
  error <- mean(abs(residuals))
  std_dev <- sd(residuals)
  
  # Store metrics in a data frame
  metrics_df <- data.frame(
    R = correlation,
    R_Squared = r_squared,
    RMSE = rmse,
    Bias = bias,
    Error = error,
    Std_Dev = std_dev
  )
  
  # Save residuals map
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  writeRaster(residuals_raster, 
              file.path(save_path, sprintf("residuals_%s_%s.envi",
                                           actual_lab,
                                           predicted_lab)),
              overwrite = TRUE)
  
  # Save metrics to a CSV file
  write.csv(format(metrics_df, digits = 3), 
            file = file.path(save_path, sprintf("metrics_%s_%s.csv",
                                                actual_lab,
                                                predicted_lab)),
            row.names = FALSE)
  
  # Print the metrics data frame
  print(metrics_df)
  
  # Return the metrics data frame
  return(metrics_df)
}

#' Calculate Metrics from a Reference Raster
#'
#' This function calculates various metrics (correlation, R-squared, RMSE, bias,
#'  error, standard deviation)
#' between a reference raster and a set of predicted rasters, 
#' generates histogram and scatterplot visualizations,
#' and saves the results.
#'
#' @param reference_raster_path Path to the reference raster file.
#' @param rasters_dir Directory containing the predicted raster files.
#' @param image Image identifier.
#' @param noises Vector of noise types.
#' @param bands Vector of band names.
#' @param mask Spatial mask for the data.
#' @param valid_indices Indices of valid data points.
#' @param estimated_var Variable being estimated.
#' @param reference_type Type of the reference data.
#' @param comparison_path Directory to save the results.
#' @param csv_name Name of the CSV file to save the comparison result.
#'
#' @return A data frame containing the calculated metrics.
#'
#' @examples
#' reference_raster_path <- "path/to/reference_raster"
#' rasters_dir <- "path/to/rasters"
#' image <- "image_identifier"
#' noises <- c("noise1", "noise2")
#' bands <- c("band1", "band2")
#' mask <- terra::rast("path/to/mask")
#' valid_indices <- complete.cases(values(mask))
#' estimated_var <- "estimated_variable"
#' reference_type <- "reference_type"
#' comparison_path <- "path/to/save/results"
#' csv_name <- "comparison_results.csv"
#' calculate_metrics_from_a_reference(reference_raster_path, rasters_dir, image,
#' noises, bands, mask, valid_indices, estimated_var, reference_type, 
#' comparison_path, csv_name)
#'
#' @export
calculate_metrics_from_a_reference <- function(reference_raster_path,
                                               rasters_dir,
                                               image,
                                               noises,
                                               bands,
                                               mask,
                                               valid_indices,
                                               estimated_var,
                                               reference_type,
                                               comparison_path,
                                               csv_name) {
  # Create an empty list to store the metrics data frames
  metrics_list <- list()
  
  reference <- open_raster_file(reference_raster_path)
  reference <- terra::project(reference, mask)
  reference <- reference * mask
  reference_values <- values(reference)
  
  for (noise in noises){
    for (band in bands){
      # Raster to compare with reference
      raster <- open_raster_file(
        file.path(rasters_dir, 
                  "PRO4SAIL_INVERSION_atbd",
                  noise,
                  band,
                  paste0(image, 
                         '_', 
                         'Refl',
                         '_',
                         estimated_var, 
                         ".envi")))
      raster <- terra::project(raster, mask)
      raster <- raster * mask
      raster_values <- values(raster)
      
      print(length(reference_values))
      print(length(raster_values))
      
      # Reference
      residuals_raster <- reference - raster
      residuals <- values(residuals_raster)
      
      # Remove NAs
      # valid_indices <- complete.cases(reference_values, raster_values)
      # valid_indices <- complete.cases(raster_values)
      reference_valuess <- reference_values[valid_indices]
      raster_values <- raster_values[valid_indices]
      residuals <- residuals[valid_indices]
      
      print(length(reference_valuess))
      print(length(raster_values))
      
      # Calculate metrics
      correlation <- cor(reference_valuess, raster_values)
      ssres <- sum((residuals)^2)
      sstot <- sum((reference_valuess - mean(reference_valuess))^2)
      # print(ssres)
      # print(sstot)
      r_squared <- 1 - (ssres / sstot)
      rmse <- sqrt(mean(residuals^2))
      bias <- mean(residuals)
      error <- mean(abs(residuals))
      std_dev <- sd(residuals)
      
      # Plots
      plot_histogram(list(reference_valuess, raster_values, residuals),
                     "LAI Distributions",
                     "LAI",
                     c(reference_type, 
                       sprintf("atbd_%s_%s", noise, band), 
                       "Residuals"),
                     comparison_path,
                     paste0(
                       "hist",
                       "_",
                       reference_type,
                       "_",
                       sprintf("atbd_%s_%s", noise, band),
                       "_",
                       "residuals"
                     ),
                     c(-5,10))
      plot_density_scatterplot(var_x = raster_values,
                               var_y = reference_valuess,
                               xlab = sprintf("atbd_%s_%s", noise, band),
                               ylab = reference_type,
                               dirname = comparison_path,
                               filename = paste0(
                                 "scatterplot",
                                 "_",
                                 reference_type,
                                 "_",
                                 sprintf("atbd_%s_%s", noise, band),
                                 "_"),
                               xlimits = c(0,10),
                               ylimits = c(0,10))
      save_envi_file(residuals_raster,
                     paste0(reference_type, 
                            "_", 
                            sprintf("atbd_%s_%s", noise, band), 
                            "_residuals"),
                     file.path(comparison_path, 
                               "residuals",
                               paste0(reference_type, 
                                      "_", 
                                      sprintf("atbd_%s_%s", noise, band))))
      
      # Create a data frame for the metrics
      metrics_df <- data.frame(
        Ref_Metrics = paste0(reference_type, "_atbd_", noise, "_", band),
        R = correlation,
        R_squared = r_squared,
        RMSE = rmse,
        Bias = bias,
        Error = error,
        Std_Dev = std_dev,
        stringsAsFactors = FALSE
      )
      
      # Append the data frame to the list
      metrics_list[[noise]] <- rbind(metrics_list[[noise]], metrics_df)
    }
  }
  
  # Combine all data frames into one
  metrics_df <- do.call(rbind, metrics_list)
  
  # Print the metrics data frame
  print(metrics_df)
  
  # Export the data frame to a tabular format (e.g., CSV)
  save_df(metrics_df, 
          comparison_path,
          csv_name)
}
