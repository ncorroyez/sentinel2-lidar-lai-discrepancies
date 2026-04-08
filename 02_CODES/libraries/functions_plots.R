# ---
# title: "function_plots.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2023-11-01"
# ---

library("terra")
library("ggplot2")
library("reshape2")
library("cowplot")
library("dplyr")
library("hrbrthemes")
# library("raster")
# library("rasterVis")
library("latticeExtra")
library("RColorBrewer")
library("gridExtra")

#' Save a raster file in ENVI format
#' 
#' This function saves a raster object in ENVI format.
#'
#' @param raster_to_save_in_envi Raster object to be saved.
#' @param output_name Name of the output file without extension.
#' @param results_path Directory where the file will be saved.
save_envi_file <- function(raster_to_save_in_envi, output_name, results_path) {
  if (!dir.exists(results_path)) {
    dir.create(results_path, showWarnings = FALSE, recursive = TRUE)
  }
  filename <- file.path(results_path, paste0(output_name, ".envi"))
  writeRaster(raster_to_save_in_envi, filename = filename, overwrite = TRUE)
  cat("ENVI file has been successfully created and saved at:", filename, "\n")
}

#' Save a raster file in TIF format
#' 
#' This function saves a raster object in TIF format.
#'
#' @param raster_to_save_in_tif Raster object to be saved.
#' @param output_name Name of the output file without extension.
#' @param results_path Directory where the file will be saved.
save_tif_file <- function(raster_to_save_in_tif, 
                          output_name,
                          results_path){
  if (!dir.exists(results_path)) {
    dir.create(results_path, showWarnings = FALSE, recursive = TRUE)
  }
  filename <- file.path(results_path, paste0(output_name, "tif"))
  writeRaster(raster_to_save_in_tif, 
              filename = filename, 
              overwrite=TRUE)
  cat("TIF file has been successfully created and saved at :", 
      filename, "\n")
}

#' Open a raster file
#' 
#' This function opens a raster file with supported extensions
#' ("tif", "tiff", "envi").
#'
#' @param raster_path Path to the raster file.
#' @return Raster object.
open_raster_file <- function(raster_path){
  file_extension <- tools::file_ext(raster_path)
  valid_extensions <- c("tif", "tiff", "envi")
  if (length(file_extension) == 1 && file_extension %in% valid_extensions) {
    return(terra::rast(raster_path))
  }
  else {
    stop("Unsupported file format. 
         Please provide a file with one of these extensions: ",  
         valid_extensions)
  }
}

#' Keep positive values in a raster
#' 
#' This function replaces negative values in a raster with zero.
#'
#' @param raster Raster object.
#' @return Raster object with non-negative values.
keep_positive_values <- function(raster) {
  raster_positive_values <- values(raster)
  raster_positive_values[raster_positive_values < 0] <- 0
  return(raster_positive_values)
}

#' Open raster file and return values
#' 
#' This function opens a raster file and returns its values.
#'
#' @param raster_path Path to the raster file.
#' @return Values of the raster.
open_raster_file_as_values <- function(raster_path) {
  raster_mask <- open_raster_file(raster_path)
  return(values(raster_mask))
}

#' Open raster file and return non-negative values
#' 
#' This function opens a raster file and returns its non-negative values.
#'
#' @param raster_path Path to the raster file.
#' @return Non-negative values of the raster.
open_raster_file_as_positive_values <- function(raster_path) {
  raster_mask <- open_raster_file(raster_path)
  return(keep_positive_values(raster_mask))
}

#' Mask low vegetation in a raster
#' 
#' This function masks low vegetation areas in a raster based on a mask.
#'
#' @param var_matrix Matrix of raster values.
#' @param var_raster Raster object to be masked.
#' @param low_vegetation Logical mask for low vegetation.
#' @return List containing masked raster and masked values.
mask_low_vegetation <- function(var_matrix, var_raster, low_vegetation){
  # Mask
  index <- low_vegetation == 1 # Index where vegetation is high
  masked_values <- var_matrix
  masked_values[!index] <- NA # Mask the low vegetation (index=0)
  
  # Take the reference raster and replace its values by the masked raster 
  masked_raster <- var_raster
  values(masked_raster) <- masked_values
  setMinMax(masked_raster)
  
  return(list(masked_raster = masked_raster,
              masked_values = masked_values)
  )
}

#' Save a basic plot
#' 
#' This function saves a basic plot to a PNG file.
#'
#' @param plot_to_save Plot object to be saved.
#' @param dirname Directory where the file will be saved.
#' @param filename Name of the output file.
#' @param title Title of the plot.
#' @param width_pixels Width of the plot in pixels.
#' @param height_pixels Height of the plot in pixels.
#' @param res Resolution of the plot.
save_basic_plot <- function(plot_to_save,
                            dirname,
                            filename,
                            title = NULL,
                            width_pixels = 1920,
                            height_pixels = 1080,
                            res = 200) {
  png(file.path(dirname, filename), 
      width = width_pixels, 
      height = height_pixels, 
      units = "px", res = res)
  plot(plot_to_save, main = title)
  dev.off()
}

#' Save an x-y plot
#' 
#' This function saves an x-y plot to a PNG file.
#'
#' @param xvar X values for the plot.
#' @param yvar Y values for the plot.
#' @param dirname Directory where the file will be saved.
#' @param filename Name of the output file.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param title Title of the plot.
#' @param width_pixels Width of the plot in pixels.
#' @param height_pixels Height of the plot in pixels.
#' @param res Resolution of the plot.
save_x_y_plot <- function(xvar,
                          yvar,
                          dirname,
                          filename,
                          xlab = NULL,
                          ylab = NULL,
                          title = NULL,
                          width_pixels = 1920,
                          height_pixels = 1080,
                          res = 200) {
  png(file.path(dirname, filename), 
      width = width_pixels, 
      height = height_pixels, 
      units = "px", res = res)
  plot(xvar, yvar, main = title, xlab=xlab, ylab=ylab)
  dev.off()
}

#' Plot density scatterplot
#' 
#' This function plots a density scatterplot and optionally saves it.
#'
#' @param var_x X values for the plot.
#' @param var_y Y values for the plot.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param title Title of the plot.
#' @param dirname Directory where the file will be saved.
#' @param filename Name of the output file.
#' @param xlimits Limits for the x-axis.
#' @param ylimits Limits for the y-axis.
plot_density_scatterplot <- function(var_x,
                                     var_y,
                                     xlab, 
                                     ylab,
                                     title = NULL,
                                     dirname = NULL, 
                                     filename = NULL,
                                     xlimits = NULL,
                                     ylimits = NULL) {
  
  if (!is.vector(var_x)) {
    var_x <- as.vector(var_x)
  }
  if (!is.vector(var_y)) {
    var_y <- as.vector(var_y)
  }
  
  data <- data.frame(
    my_x = var_x,
    my_y = var_y
  )
  
  if (is.null(xlimits)){
    xlimits <- c(min(data$my_x - 0.5), max(data$my_x + 0.5))
  }
  if (is.null(ylimits)){
    ylimits <- c(min(data$my_y - 0.5), max(data$my_y + 0.5))
  }
  
  lm_model <- lm(data$my_y ~ data$my_x, data = data)
  # print(coef(lm_model))
  
  # equation <- sprintf("italic(y) == %.3f * italic(x) + %.3f ",
  #                     coef(lm_model)[2],
  #                     coef(lm_model)[1])
  equation <- sprintf("y2 = %.3f x + %.3f",
                      round(coef(lm_model)[2],2),
                      round(coef(lm_model)[1],2))
  
  statsReg1 <- cor.test(data$my_x, data$my_y)$estimate
  statsReg21 <- summary(lm_model)$r.squared
  
  ssres <- sum((data$my_y - data$my_x)^2)
  sstot <- sum((data$my_y - mean(data$my_y))^2)
  # print(ssres)
  # print(sstot)
  r_squared <- 1 - (ssres / sstot)
  
  statsReg2 <- r_squared
  statsReg3 <- Metrics::rmse(na.omit(data$my_y), na.omit(data$my_x))
  
  # Set the text size parameters
  text_size <- 24
  
  p3 <- ggplot(data, aes(x=var_x, y=var_y)) +
    # ggtitle(paste(sprintf("R: %.3f,", statsReg1), 
    #               sprintf("R² (model): %.3f,", statsReg21),
    #               sprintf("R² (formula): %.3f,", statsReg2),
    #               sprintf("RMSE: %.3f,", statsReg3),
    #               sprintf("Equation: %s", equation))) +
    ggtitle(title) +
    xlab(xlab) + ylab(ylab) +
    scale_fill_continuous(type = "viridis") +
    scale_x_continuous(limits = xlimits) +
    scale_y_continuous(limits = ylimits) +
    # geom_bin2d(bins = 400) +
    # geom_smooth(method="lm" , color="blue", se=TRUE, lwd=1.5) +
    geom_bin2d(bins = 400) +  # Fill aesthetic mapped to bin count
    geom_smooth(method="lm" , color="blue", se=TRUE, lwd=1.5) +
    geom_abline(slope = 1, intercept = 0, color="red", lwd=1.5) + 
    theme(
      # Adjust main title size
      plot.title = element_text(size = text_size * 1.2, face = "bold", hjust = 0.5),
      
      # Adjust axis label size
      axis.title.x = element_text(size = text_size * 1.2, color = "black"),
      axis.title.y = element_text(size = text_size * 1.2, color = "black"),
      axis.text = element_text(size = text_size * 1.2, color = "black"),
      legend.text = element_text(size = text_size * 1.2, color = "black")
      
      # Adjust text annotations size
      # Uncomment the lines below if you want to add text annotations
      # axis.text = element_text(size = text_size * 0.8), 
      # axis.text.x = element_text(size = text_size * 0.8),
      # axis.text.y = element_text(size = text_size * 0.8)
    ) +
    annotate("text", x = 8.5, y = 5,
             label = sprintf("%s", equation), size = text_size * 0.5, color = "black") +
    annotate("text", x = 9, y = 8,
             label = "y1 = x", size = text_size * 0.5, color = "black") +
    annotate("text", x = 9, y = 2,
             label = sprintf("R: %.2f", statsReg1), size = text_size * 0.5, color = "black")
  # annotate("text", x = 1.0685, y = 9,
  #          label = sprintf("R²: %.2f", statsReg21), size = text_size * 0.5, color = "black") +
  # coord_fixed(ratio = 1)
  plot(p3)
  
  if (!is.null(dirname) && !is.null(filename)){
    file_path <- file.path(dirname, filename)
    ggsave(filename = paste0(file_path, ".png"),
           plot = p3, device = "png", scale = 1, 
           width = 1920, height = 1080, units = "px", dpi = 100)
    cat("Plots saved at: ", paste0(file_path, ".png"), "\n")
  }
  return(p3)
}

#' Plot Histogram
#'
#' This function creates and optionally saves a histogram for one or more
#' sets of variables.
#'
#' @param vars_list A list or vector of numeric values to be plotted.
#' @param title A character string specifying the title of the plot. 
#' Default is "LAI Distributions".
#' @param xlab A character string specifying the x-axis label. Default is "LAI".
#' @param var_labs A vector of character strings specifying the labels 
#' for the variables.
#' @param dirname A character string specifying the directory name 
#' where the plot will be saved.
#' @param filename A character string specifying the name of the file 
#' to save the plot.
#' @param limits A numeric vector of length 2 specifying the limits of 
#' the x-axis. Default is the range of the data.
#'
#' @return A histogram plot. Optionally saves the plot as a PNG file 
#' if `dirname` and `filename` are provided.
#'
#' @examples
#' vars_list <- list(rnorm(1000), rnorm(1000, mean = 3))
#' var_labs <- c("Group 1", "Group 2")
#' plot_histogram(vars_list, title = "My Histogram", var_labs = var_labs, 
#' dirname = "plots", filename = "histogram")
#'
#' @export
plot_histogram <- function(vars_list, 
                           title = NULL, 
                           xlab = NULL, 
                           var_labs, 
                           dirname = NULL, 
                           filename = NULL,
                           limits = NULL # c(0.1, 10)
) {
  
  if (length(vars_list) == 0) {
    cat("No variables provided. Exiting.\n")
    return(NULL)
  }
  if (is.null(title)) {
    title <- "LAI Distributions"
  }
  if (is.null(xlab)) {
    xlab <- "LAI"
  }
  if (is.matrix(vars_list) 
      || is.array(vars_list)
      || is.double(vars_list)) { # A single variable (i.e. not a list) was provided
    vars_list <- as.vector(vars_list)
    fill_colors = brewer.pal(3, "Set2")
    data <- data.frame(
      type = rep(as.character(var_labs), each = length(vars_list)),
      value = vars_list
    )
    if (is.null(limits)){
      limits = c(min(data$value, na.rm = TRUE), max(data$value, na.rm = TRUE))
    }
    
    h <- ggplot(data, aes(x = value, fill = type)) +
      geom_histogram(bins = 500, alpha = 0.4, position = "identity") +
      scale_fill_manual(values = fill_colors[1:length(vars_list)]) +
      labs(x = xlab, y = "Frequency", fill = "") +
      ggtitle(title) +
      scale_x_continuous(limits = limits)
  }
  if (is.list(vars_list)) { # A list is provided
    fill_colors = brewer.pal(length(vars_list), "Set2")
    data <- data.frame(
      type = rep(var_labs, each = length(vars_list[[1]])),
      value = do.call(c, vars_list) # do.call(c, vars_list) # unlist(vars_list)
    )
    if (is.null(limits)){
      limits = c(min(data$value, na.rm = TRUE), max(data$value, na.rm = TRUE))
    }
    
    h <- ggplot(data, aes(x = value, fill = type)) +
      geom_histogram(bins = 500, alpha = 0.6, position = "identity") +
      scale_fill_manual(values = fill_colors[1:length(vars_list)]) +
      labs(x = xlab, y = "Frequency", fill = "") +
      ggtitle(title) +
      scale_x_continuous(limits = limits)
  }
  print(h)
  
  if (!is.null(dirname) && !is.null(filename)) {
    file_path <- file.path(dirname, filename)
    ggsave(filename = paste0(file_path, ".png"),
           plot = h, device = "png", scale = 1, 
           width = 1920, height = 1080, units = "px", dpi = 200)
    cat("Plot is successfully saved at: ", paste0(file_path, ".png"), "\n")
  } else {
    cat("Plot is not saved because", "dirname =", dirname, 
        "and filename = ", filename, "\n")
  }
}

#' Correlation Test Function
#'
#' This function performs a correlation test between two numeric variables 
#' and prints the result.
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param method A character string specifying the method for the correlation
#'  test. Default is "pearson". Other options include "kendall" and "spearman".
#'
#' @return The estimated correlation coefficient.
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' correlation <- correlation_test_function(x, y)
#'
#' @export
correlation_test_function <- function(x, y, method = "pearson") {
  correlation_test <- cor.test(x, y, method = method)
  print(correlation_test)
  
  p_value <- correlation_test$p.value
  if (p_value < 0.05) {
    cat("Correlation is statistically significant (p < 0.05)\n")
  } else {
    cat("Correlation is not statistically significant (p >= 0.05)\n")
  }
  return(correlation_test$estimate)
}

#' Plot Distributions
#'
#' This function creates and optionally saves histograms 
#' for each variable in a data frame.
#'
#' @param InputPROSAIL A data frame containing the variables to be plotted.
#' @param save_plots A logical value indicating whether to save the plots. 
#' Default is \code{FALSE}.
#' @param dirname A character string specifying the directory name 
#' where the plots will be saved.
#' @param filename A character string specifying the name of the file 
#' to save the plots.
#'
#' @return A grid of histogram plots. Optionally saves the plots as a PNG file
#'  if \code{save_plots} is \code{TRUE} and 
#'  both \code{dirname} and \code{filename} are provided.
#'
#' @examples
#' InputPROSAIL <- data.frame(var1 = rnorm(1000), var2 = rnorm(1000, mean = 3))
#' plot_distributions(InputPROSAIL, save_plots = TRUE, dirname = "plots", 
#' filename = "distributions")
#'
#' @export
plot_distributions <- function(InputPROSAIL,
                               save_plots = FALSE,
                               dirname = NULL,
                               filename = NULL) {
  melted_df <- melt(InputPROSAIL)
  custom_green <- rgb(0, 200, 30, maxColorValue = 255)
  
  plot_grid <- lapply(unique(melted_df$variable), function(var) {
    ggplot(data = subset(melted_df, variable == var), aes(value)) +
      geom_histogram(bins = 30, fill = custom_green, color = "black") +
      ggtitle(paste("Histogram of", var)) +
      xlab(var) + ylab("Frequency")
  })
  
  final_plot <- plot_grid(plotlist = plot_grid, ncol = 4)
  
  if (save_plots) {
    if (!is.null(dirname) & !is.null(filename)) {
      if (!dir.exists(dirname)) {
        dir.create(dirname)
      }
      file_path <- file.path(dirname, filename)
      ggsave(filename = paste0(file_path, ".png"), final_plot, 
             width = 1920, height = 1080, units = "px", dpi = 100)
      cat("Plots saved at: ", paste0(file_path, ".png"), "\n")
    } else {
      warning("Please provide both directory name and file name for saving plots.")
    }
  } else {
    print(final_plot)
  }
}

create_data_one_site <- function(site,
                                 all_observations,
                                 field_data,
                                 spatial_sampled_data) {
  # Identify common columns across all data frames for the given site
  if (is.null(all_observations) || 
      is.null(field_data) || 
      is.null(spatial_sampled_data)) {
    stop(paste("Data for site", site, "is missing in one or more sources."))
  }
  
  common_cols <- dplyr::intersect(
    dplyr::intersect(names(all_observations), 
                     names(field_data)),
    names(spatial_sampled_data)
  )
  # Remove unwanted variables
  excluded_vars <- c("stands", "lat", "lon")
  common_cols <- setdiff(common_cols, excluded_vars)
  common_cols <- c("lcv", "s2_lai", "mean")
  print(common_cols)
  
  # Select only the common columns
  all <- all_observations %>%
    dplyr::select(all_of(common_cols)) %>%
    mutate(source = "All")
  
  field <- field_data %>%
    dplyr::select(all_of(common_cols)) %>%
    mutate(source = "Field")
  
  sampled <- spatial_sampled_data %>%
    dplyr::select(all_of(common_cols)) %>%
    mutate(source = "SYS Sampling")
  
  # Combine all data
  combined_data <- bind_rows(all, field, sampled)
  
  # Reshape into long format for plotting
  long_data <- combined_data %>%
    tidyr::pivot_longer(cols = common_cols,
                        names_to = "variable",
                        values_to = "value")
  
  # Optionally order variables for consistent plotting
  long_data$variable <- factor(long_data$variable, levels = common_cols)
  
  return(long_data)
}

plot_data_one_site <- function(site,
                               all_observations,
                               field_data,
                               spatial_sampled_data,
                               output_dir,
                               save = FALSE) {
  long_data <- create_data_one_site(site,
                                    all_observations,
                                    field_data,
                                    spatial_sampled_data)
  
  plot <- ggplot(long_data, aes(x = source, y = value, fill = variable)) +
    # Histogram with density
    ggdist::stat_halfeye(
      # adjust bandwidth
      adjust = 1,
      # move to the right
      justification = -0.2,
      # remove the slub interval
      .width = 0,
      point_colour = NA
    ) +
    geom_boxplot(aes(x = source, y = value),
                 width = 0.2,
                 fill = "white",
                 color = "black",
                 outlier.shape = NA,
                 alpha = 0.5) +  # Boxplot without outliers
    # Individual points (jittered) for the raincloud effect
    # ggdist::stat_dots(
    #   ## orientation to the left
    #   side = "left",
    #   ## move geom to the left
    #   justification = 1.2
    #   ## adjust grouping (binning) of observations
    #   # binwidth = .25
    # ) +
    ## remove white space on the sides
    # coord_cartesian(xlim = c(1.3, 2.9)) +
    facet_wrap(~ variable, scales = "free") +
    labs(title = paste("Raincloud Plots for", site),
         x = "Data Source",
         y = "Value") +
    # theme_bw() +
    tidyquant::scale_fill_tq() +
    tidyquant::theme_tq() +
    theme(legend.position = "none") +
    coord_flip()
  if (save){
    ggsave(filename = file.path(output_dir, paste0(site, "_all_raincloud_plots.png")),
           plot = plot)
  }
  else {
    print(plot)
  }
}

create_violin_data_one_site <- function(site,
                                        sys_data, # all_observations,
                                        field_data,
                                        spatial_sampled_data,
                                        variables = c("lidar_lai", "s2_lai")) {
  sys <- sys_data %>%
    dplyr::select(all_of(variables)) %>%
    mutate(source = "SYS Sampling")
  
  field <- field_data %>%
    dplyr::select(all_of(variables)) %>%
    mutate(source = "Field")
  
  sampled <- spatial_sampled_data %>%
    dplyr::select(all_of(variables)) %>%
    mutate(source = "Proposed Stratified Sampling")
  
  combined_data <- bind_rows(field, sys,
                             sampled
                             )
  # combined_data <- combined_data %>%
  #   mutate(deltaLAI = lidar_lai - s2_lai)
  long_data <- combined_data %>%
    tidyr::pivot_longer(cols = all_of(variables),
                        names_to = "variable",
                        values_to = "value")
  long_data$variable <- factor(long_data$variable, levels = variables)
  long_data$source <- factor(long_data$source, levels = c("Field", "SYS Sampling",
                                                          "Proposed Stratified Sampling"
                                                          ))
  return(long_data)
}

plot_violin_one_site <- function(site,
                                 all_observations,
                                 field_data,
                                 spatial_sampled_data,
                                 output_dir,
                                 save = FALSE) {
  long_data <- create_violin_data_one_site(site,
                                           all_observations,
                                           field_data,
                                           spatial_sampled_data)
  
  plot <- ggplot(long_data, aes(x = source, y = value, fill = variable)) +
    geom_violin(trim = FALSE) +
    facet_wrap(~ variable, scales = "free_y") +
    labs(title = paste("Violin Plots for", site),
         x = "Data Source",
         y = "Value") +
    theme_bw() +
    theme(legend.position = "none") # Hide legend for cleaner plot
  if (save){
    ggsave(filename = file.path(output_dir, paste0(site, "_violin_plot.png")),
           plot = plot)
  }
  else {
    print(plot)
  }
}

plot_raincloud_one_site <- function(site,
                                    sys_data, # all_observations
                                    field_data,
                                    spatial_sampled_data,
                                    output_dir,
                                    variables,
                                    save = FALSE) {
  long_data <- create_violin_data_one_site(site,
                                           sys_data,
                                           field_data,
                                           spatial_sampled_data,
                                           variables)
  long_data <- long_data %>%
    group_by(variable) %>%
    mutate(density = (value - min(value)) / (max(value) - min(value))) %>%
    ungroup()
  
  plot <- ggplot(long_data, aes(x = source, y = value, fill = variable)) +
    # Histogram with density
    ggdist::stat_halfeye(
      # data = long_data,
      # aes(x = source, y = density, fill = variable),
      # adjust bandwidth
      adjust = 0.4,
      # scale = 1,
      # move to the right
      # justification = -0.2,
      # remove the slub interval
      .width = 0,
      point_colour = NA
    ) +
    geom_boxplot(aes(x = source, y = value),
                 width = 0.2,
                 fill = "white",
                 color = "black",
                 outlier.shape = NA,
                 # justification = 0.1,
                 alpha = 0.5) +  # Boxplot without outliers
    # Individual points (jittered) for the raincloud effect
    # ggdist::stat_dots(
    #   ## orientation to the left
    #   side = "left",
    #   ## move geom to the left
    #   justification = 1.2,
    #   ## adjust grouping (binning) of observations
    #   binwidth = NA
    # ) +
    ## remove white space on the sides
    # coord_cartesian(xlim = c(1.3, 2.9)) +
    facet_wrap(~ variable, scales = "free") +
    labs(title = paste("Raincloud Plots for", site),
         x = "Data Source",
         y = "Value") +
    # theme_bw() +
    tidyquant::scale_fill_tq() +
    tidyquant::theme_tq() +
    theme(legend.position = "none")
    # coord_flip()
  if (save){
    ggsave(filename = file.path(output_dir, paste0(site, "_raincloud_plot.png")),
           plot = plot,
           height = 8,  # Increase height
           width = 16,  # Increase width
           units = "in")  # Set units as inches
  }
  else {
    print(plot)
  }
}

create_violin_data <- function(all_observations,
                               field_data,
                               spatial_sampled_data,
                               sites) {
  combined_data_list <- list()
  for (site in sites) {
    # obs_data <- all_observations[[site]] %>%
    #   dplyr::select(lidar_lai, s2_lai) %>%
    #   mutate(source = "All", site = site)
    
    filtered_data <- field_data[[site]] %>%
      dplyr:: select(lidar_lai, s2_lai) %>%
      mutate(source = "Field", site = site)
    
    variogram_data <- spatial_sampled_data[[site]] %>%
      dplyr::select(lidar_lai, s2_lai) %>%
      mutate(source = "Variogram 500m", site = site)
    
    combined_site_data <- bind_rows(
      # obs_data,
      filtered_data,
      variogram_data)
    combined_site_data <- combined_site_data %>%
      mutate(deltaLAI = lidar_lai - s2_lai)
    combined_data_list[[site]] <- combined_site_data
  }
  
  combined_data <- bind_rows(combined_data_list)
  long_data <- combined_data %>%
    tidyr::pivot_longer(cols = c("lidar_lai", "s2_lai", "deltaLAI"),
                 names_to = "variable",
                 values_to = "value")
  long_data$variable <- factor(long_data$variable,
                               levels = c("lidar_lai", "s2_lai", "deltaLAI"))
  
  return(long_data)
}

plot_violin <- function(all_observations,
                        field_data,
                        spatial_sampled_data,
                        sites,
                        output_dir,
                        save = FALSE) {
  long_data <- create_violin_data(all_observations,
                                  field_data,
                                  spatial_sampled_data,
                                  sites)
  
  plot <- ggplot(long_data, aes(x = source, y = value, fill = variable)) +
    geom_violin(trim = FALSE) +
    facet_wrap(~ site, scales = "free_y") +
    labs(title = "Violin Plots for All Sites",
         x = "Data Source",
         y = "Value") +
    theme_bw() # +
  # theme(legend.position = "none")
  if (save){
    ggsave(filename = file.path(output_dir, "all_sites_violin_plot.png"),
           plot = plot,
           height = 8,  # Increase height
           width = 16,  # Increase width
           units = "in")  # Set units as inches
  }
  else {
    print(plot)
  }
}

plot_raincloud_all_sites <- function(all_observations,
                                     field_data,
                                     spatial_sampled_data,
                                     sites,
                                     output_dir,
                                     save = FALSE) {
  # Generate the long-format data for all sites
  long_data <- create_violin_data(all_observations,
                                  field_data,
                                  spatial_sampled_data,
                                  sites)
  
  # Create the raincloud plot
  plot <- ggplot(long_data, aes(x = source, y = value, fill = variable)) +
    # Raincloud density
    stat_halfeye(
      adjust = 0.4,  # Adjust bandwidth
      justification = -0.2,  # Move density slightly left
      .width = 0,  # Disable slab intervals
      point_colour = NA  # Remove point indicators on density
    ) +
    # Boxplot in the center of the raincloud
    geom_boxplot(
      width = 0.2,
      fill = "white",
      color = "black",
      outlier.shape = NA,
      alpha = 0.5
    ) +
    # Jittered dots for individual observations
    ggdist::stat_dots(
      side = "left",  # Align dots to the left
      justification = 1.2  # Position them to the left of the boxplot
    ) +
    # Facet for different sites
    facet_wrap(~ site, scales = "free_y") +
    # Add titles and labels
    labs(
      title = "Raincloud Plots for All Sites",
      x = "Data Source",
      y = "Value"
    ) +
    # Aesthetic adjustments
    theme_tq() +
    scale_fill_tq() +
    theme(
      legend.position = "none"
    ) +
    coord_flip()
  
  # Save the plot
  if (save){
    ggsave(
      filename = file.path(output_dir, "all_sites_raincloud_plot.png"),
      plot = plot,
      height = 8,  # Increase height
      width = 16,  # Increase width
      units = "in"  # Set units as inches
    )
  }
  else {
    print(plot)
  }
}

plot_allometry_with_colorbar <- function(allometry_data, 
                                         title = "Plot",
                                         key_point = NULL, 
                                         x_position = 3, 
                                         y_position = 8,
                                         label_size = 1.5) {
  
  # Generate the custom green color palette (reverse gradient)
  # green_palette <- rev(colorRampPalette(c("lightgreen", "yellowgreen", "darkgreen"))(length(allometry_data$height)))
  # green_palette <- rev(colorRampPalette(c("darkgreen", "mediumseagreen", "lightgreen"))(length(allometry_data$height)))
  green_palette <- rev(colorRampPalette(c("darkgreen", "forestgreen", "limegreen", "yellowgreen", "lightgreen"))(length(allometry_data$height)))
  # Set up plot margins to leave space for the color bar
  par(mar = c(5, 5, 2, 8), xpd = TRUE)
  
  # Create the plot with no points but with color for the line
  plot(allometry_data$density, allometry_data$height, type = "n",  # 'n' means no points
       xlab = "LAD (m²/m²)", ylab = "Canopy Depth (m)",
       ylim = c(max(allometry_data$height), 0),  # Reverse y-axis so 0 is at the top
       main = title,
       axes = TRUE,   # Show default axes
       bty = "n",      # Remove border
       col.main = "black",   # Title color
       col.lab = "black",    # Label color
       cex.main = label_size * 1.2,       # Title size
       cex.lab = label_size,        # Label size
       cex.axis = label_size)        # Label size
  
  # Add optional key point annotation if provided
  if (!is.null(key_point)) {
    text(x_position, y_position, key_point, col = "blue", cex = label_size)
  }
  
  # Add the line with the custom green color gradient based on height
  for (i in 1:(length(allometry_data$height)-1)) {
    segments(allometry_data$density[i], allometry_data$height[i], 
             allometry_data$density[i+1], allometry_data$height[i+1], 
             col = green_palette[i], lwd = 4)  # Adjust line width and color
  }
  
  # Draw the shaded area under the curve
  # polygon(
  #   x = c(allometry_data$density, min(allometry_data$density)),
  #   y = c(allometry_data$height, max(allometry_data$height)),
  #   col = adjustcolor("green", alpha.f = 0.2),
  #   border = NA  # No border for the filled area
  # )
  
  # Add the color bar (image plot)
  image.plot(
    x = c(1.1, 1.2),  # Position of the color bar
    y = c(min(allometry_data$height), max(allometry_data$height)),  # Y range matching height
    z = matrix(seq(min(allometry_data$height), max(allometry_data$height), length.out = 256), ncol = 1),  # Value range for the color scale
    col = green_palette,  # Use the custom green color palette
    legend.only = TRUE,  # Only show the legend (color bar)
    legend.shrink = 0.9,  # Shrink the legend
    legend.width = 1.2,   # Set width of the color bar
    legend.mar = 5,       # Set margin for the color bar
    horizontal = FALSE,   # Set vertical color bar
    axis.args = list(cex.axis = label_size, tck = -0.02, labels = FALSE) # Remove axis labels
  )
  
  # Arrow
  arrow_x <- 0.48
  arrows(x0 = arrow_x, y0 = min(allometry_data$height), 
         x1 = arrow_x, y1 = max(allometry_data$height), 
         length = 0.1, col = "black", lwd = 4)
  # Min / Max LAD
  text(x = arrow_x, y = min(allometry_data$height), 
       labels = min(allometry_data$height), 
       pos = 3, cex = label_size)
  text(x = arrow_x, y = max(allometry_data$height), 
       labels = sum(allometry_data$density), 
       pos = 1, cex = label_size)
  
  # Add a title to the color bar (optional)
  mtext("Cumulative Profiles", side = 4, line = 2, cex = label_size)  # Title for the color bar
  
  # Add a custom title for the plot
  # mtext(title, side = 3, line = 2, cex = 1.5)  # Add plot title
}

plot_lai_scatter <- function(x, y, site,
                             xlab = "LiDAR LAI",
                             ylab = "Sentinel-2 LAI",
                             xlim = c(0, 15),
                             ylim = c(0, 15),
                             add_annotation = TRUE,
                             annotation_pos = c("topleft", "bottomright"),
                             remove_y = FALSE) {
  
  set.seed(42)
  # Prepare data
  df <- na.omit(data.frame(x = as.numeric(x), y = as.numeric(y)))
  df <- sample_n(df, 5000)
  
  # Fit model and get stats
  fit   <- lm(y ~ x, data = df)
  stats <- broom::glance(fit)
  coefs <- coef(fit)
  rmse_val <- Metrics::rmse(df$y, predict(fit))
  r_val <- cor(df$x, df$y, method = "pearson")
  
  # Prepare annotation text
  eqn_text  <- sprintf("y = %.2f x + %.2f", coefs[2], coefs[1])
  r2_text   <- sprintf("R² = %.2f", stats$r.squared)
  rmse_text <- sprintf("RMSE = %.2f", rmse_val)
  r_text   <- sprintf("r = %.2f", r_val)
  
  if (remove_y){
    ylab <- NULL
  }
  
  # Base plot
  p <- ggplot(df, aes(x = x, y = y)) +
    # geom_bin2d(bins = 400) +
    # scale_fill_viridis_c(name = "Count", option = "D") +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", color = "blue", se = TRUE, size = 1) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    coord_fixed(ratio = 1) +
    xlim(xlim) + ylim(ylim) +
    labs(title = NULL, x = xlab, y = ylab) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "right")
  
  # Add annotations if requested
  if (isTRUE(add_annotation)) {
    pos <- match.arg(annotation_pos)
    if (pos == "topleft") {
      x_annot <- xlim[1] + 0.02 * diff(xlim)
      y_vals  <- rev(ylim[2] - c(0.05, 0.12, 0.19) * diff(ylim))
      hjust   <- 0
    } else {
      x_annot <- xlim[2] - 0.02 * diff(xlim)
      y_vals  <- ylim[1] + c(0.19, 0.12, 0.05) * diff(ylim)
      hjust   <- 1
    }
    # p <- p +
    # annotate("text", x = x_annot, y = y_vals[1], label = eqn_text,  hjust = hjust, size = 5) +
    # annotate("text", x = x_annot, y = y_vals[2], label = r2_text,   hjust = hjust, size = 5) +
    # annotate("text", x = x_annot, y = y_vals[3], label = rmse_text, hjust = hjust, size = 5)
    p <- p +
      annotate("text", x = x_annot, y = y_vals[1], label = eqn_text,  hjust = hjust, size = 5) +
      annotate("text", x = x_annot, y = y_vals[2], label = rmse_text, hjust = hjust, size = 5) +
      annotate("text", x = x_annot, y = y_vals[3], label = r_text,    hjust = hjust, size = 5)
  }
  # ggsave(filename, width = 10, height = 10, dpi = 300)
  # print(p)
  
  return(p)
}
