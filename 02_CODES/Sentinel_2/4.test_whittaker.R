# ---
# title: "4.test_whittaker.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-02-13"
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
library(terra)
library(zoo)
library(readODS)
library(stringr)
library(signal)
library(ptw)
source("../libraries/functions_create_masks.R")
source("../libraries/functions_sentinel_2.R")
source("../libraries/functions_plots.R")
source("../libraries/functions_general_tools.R")

# 2. PARAMETRES ---------------------------------------------------------------
set.seed(42)
data_dir <- '/media/corroyez/MyPassport/01_DATA'
results_dir <- '../../03_RESULTS'
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- c("Mormal")
distribs <- c("atbd", "atbd_optim_common")

# Grille temporelle commune pour le Machine Learning
common_dates <- seq(from = as.Date("2021-06-01"), 
                    to   = as.Date("2021-09-30"), 
                    by   = "10 days")

# Dates d'acquisition LiDAR par site (Format: Named Vector)
lidar_dates <- c(
  "Aigoual" = "2021-07-18",
  "Blois"   = "2021-06-17",
  "Mormal"  = "2021-06-16"
)

# Lecture du fichier de dates S2
dates_file <- file.path(data_dir, "S2_Dates_Summer2021.ods")
data_dates <- read_ods(dates_file)

# 3. BOUCLE PRINCIPALE --------------------------------------------------------

for (site in sites) {
  cat("\n========================================\n")
  cat("PROCESSING SITE:", site, "\n")
  
  # Chemins et Dates
  res_path <- file.path(results_dir, site, "Metrics/Deciduous_Only")
  dateAcqs <- get_site_full_dates(site, data_dates) # Votre fonction
  lidar_date_site <- as.Date(lidar_dates[site])
  
  # Chargement du LiDAR de référence
  path_lidar <- file.path(res_path, "lidarlai_res_10_m.tif")
  
  if(!file.exists(path_lidar)) {
    cat("[ERREUR] LiDAR non trouvé pour", site, ":", path_lidar, "\n")
    next
  }
  r_lidar <- rast(path_lidar)
  
  for (distrib in distribs) {
    cat("\n  -> Distribution:", distrib, "\n")
    
    # A. Chargement et Tri des images S2
    # ----------------------------------
    pattern <- paste0("^s2lai_(", paste(dateAcqs, collapse = "|"), ")_", distrib, "_res_10_m\\.tif$")
    files <- list.files(res_path, pattern = pattern, full.names = TRUE)
    
    if(length(files) < 2) { cat("     [SKIP] Pas assez de fichiers S2.\n"); next }
    
    # Tri chronologique (Sécurité)
    dates_found <- as.Date(str_extract(basename(files), "\\d{4}-\\d{2}-\\d{2}"))
    ord <- order(dates_found)
    files <- files[ord]
    current_dates <- dates_found[ord]
    
    r_stack <- rast(files)
    terra::time(r_stack) <- current_dates
    
    # B. Gap-Filling (Interpolation Linéaire)
    # ---------------------------------------
    cat("     1. Gap-Filling (Interpolation)...\n")
    # On remplit les trous avant tout lissage ou calibration
    r_filled <- terra::approximate(r_stack, method = "linear", rule = 2)
    
    
    # C. Calibration du Lambda Whittaker avec le LiDAR
    # ------------------------------------------------
    cat("     2. Calibration du Lambda (Whittaker vs LiDAR)...\n")
    
    # On extrait un échantillon aléatoire de 2000 pixels (LiDAR + Série S2)
    # pour tester quel lambda donne la courbe passant le plus près du LiDAR.
    # pts_sample <- spatSample(r_lidar, size=10000, method="random", na.rm=TRUE, as.points=TRUE, values=TRUE)
    
    val_p98 <- global(r_lidar, 
                      fun = function(x) quantile(x, 0.98, na.rm=TRUE))[[1]]
    breaks <- seq(from = 2, to = floor(val_p98), by = 1)
    breaks <- c(breaks, floor(val_p98) + 1)
    breaks <- unique(breaks)
    r_strata <- classify(r_lidar, rcl = breaks, include.lowest = TRUE)
    r_strata <- mask(r_strata, r_filled[[1]])
    n_bins <- length(breaks) - 1
    target_total <- 5000
    pts_per_bin <- ceiling(target_total / n_bins)
    pts_sample <- spatSample(r_strata, 
                             size = pts_per_bin, 
                             method = "stratified", 
                             xy = TRUE, na.rm = TRUE, values = FALSE)
    # stop()
    if(length(pts_sample) > 0) {
      # vals_lidar_sample <- pts_sample$V1
      vals_lidar_sample <- extract(r_lidar, pts_sample)
      vals_s2_sample <- extract(r_filled, pts_sample)#[,-1] # Séries S2 (sans ID)
      # stop()
      # vals_lidar_sample <- values(r_lidar, mat = FALSE) 
      # vals_s2_sample    <- values(r_filled, mat = TRUE)
      
      # Lambdas à tester
      lambdas_test <- c(0.01, 0.1, 1, 2, 10, 15, 50, 100, 500)
      rmse_results <- numeric(length(lambdas_test))
      
      # Conversion dates en numérique pour approx()
      s2_dates_num <- as.numeric(current_dates)
      lidar_date_num <- as.numeric(lidar_date_site)
      
      for(i in seq_along(lambdas_test)) {
        lam <- lambdas_test[i]
        
        # Fonction calcul erreur pour un pixel
        calc_pixel_error <- function(y_s2, target_lidar) {
          if(any(is.na(y_s2))) return(NA)
          # Lissage
          y_smooth <- tryCatch(ptw::whit2(y_s2, lambda = lam), error=function(e) y_s2)
          # Valeur prédite à la date LiDAR
          pred_at_lidar <- approx(s2_dates_num, y_smooth, xout=lidar_date_num)$y
          return((pred_at_lidar - target_lidar)^2)
        }
        
        # Appliquer sur l'échantillon (mapply est plus rapide qu'une boucle)
        # as.list(as.data.frame(t(vals_s2_sample))) transforme le data.frame en liste de vecteurs (lignes)
        errors <- mapply(calc_pixel_error, as.list(as.data.frame(t(vals_s2_sample))), vals_lidar_sample[[1]])
        rmse_results[i] <- sqrt(mean(errors, na.rm=TRUE))
      }
      
      best_lambda <- lambdas_test[which.min(rmse_results)]
      cat("        => Lambda Optimal trouvé :", best_lambda, "(RMSE:", round(min(rmse_results), 2), ")\n")
      
    } else {
      cat("        [WARN] Pas de points communs S2/LiDAR trouvés. Lambda par défaut = 15.\n")
      best_lambda <- 15
    }
    # stop()
    
    # D. Application du Lissage (Whittaker Optimal)
    # ---------------------------------------------
    cat("     3. Application du Lissage (Whit2)...\n")
    
    whittaker_fun <- function(y) {
      if (all(is.na(y))) return(y)
      # On utilise le lambda calibré
      return(ptw::whit2(y, lambda = best_lambda)) 
    }
    
    r_smoothed <- app(r_filled, fun = whittaker_fun, cores = 1)
    names(r_smoothed) <- paste0("LAI_Smooth_", current_dates)
    
    # Sauvegarde intermédiaire (Série lissée sur dates d'origine)
    out_name_smooth <- file.path(res_path, paste0(site, "_LAI_SmoothedW_", distrib, ".tif"))
    # writeRaster(r_smoothed, out_name_smooth, overwrite = TRUE)
    
    
    # E. Harmonisation (Grille Commune)
    # ---------------------------------
    cat("     4. Harmonisation temporelle (10 jours)...\n")
    
    # Il est important d'assigner le temps au raster lissé pour que l'interpolation fonctionne
    terra::time(r_smoothed) <- current_dates
    
    resample_fun <- function(y) {
      if (all(is.na(y))) return(rep(NA, length(common_dates)))
      return(approx(x = as.numeric(current_dates), 
                    y = y, 
                    xout = as.numeric(common_dates), 
                    rule = 2)$y)
    }
    
    r_harmonized <- app(r_smoothed, fun = resample_fun, cores = 1)
    names(r_harmonized) <- paste0("LAI_", common_dates)
    terra::time(r_harmonized) <- common_dates
    
    # Sauvegarde Série Harmonisée (Non corrigée LiDAR)
    out_name_harm <- file.path(res_path, paste0(site, "_LAI_HarmonizedW_10d_", distrib, ".tif"))
    # writeRaster(r_harmonized, out_name_harm, overwrite = TRUE)
    
    
    # F. Correction au Prorata (LiDAR Fusion)
    # ---------------------------------------
    cat("     5. Correction LiDAR (Ratio/Prorata)...\n")
    
    # 1. Calculer la valeur S2 EXACTE à la date du LiDAR (pixel par pixel)
    # On utilise la série lissée (r_smoothed) car elle contient la dynamique propre
    get_s2_at_lidar_date <- function(y) {
      if (all(is.na(y))) return(NA)
      return(approx(x = as.numeric(current_dates), 
                    y = y, 
                    xout = as.numeric(lidar_date_site), 
                    rule = 2)$y)
    }
    
    r_s2_at_lidar_date <- app(r_smoothed, fun = get_s2_at_lidar_date, cores = 1)
    
    # 2. Alignement spatial (sécurité)
    if (!compareGeom(r_lidar, r_s2_at_lidar_date, stopOnError = FALSE)) {
      r_lidar_proj <- project(r_lidar, r_s2_at_lidar_date, method = "bilinear")
    } else {
      r_lidar_proj <- r_lidar
    }
    
    # 3. Calcul du Ratio (LiDAR / S2)
    # Si S2 est ~0, on clamp pour éviter l'explosion
    r_ratio <- r_lidar_proj / r_s2_at_lidar_date
    r_ratio <- clamp(r_ratio, lower = 0.1, upper = 10, values = FALSE)
    r_ratio <- ifel(is.finite(r_ratio), r_ratio, 1) # Remplacer NA/Inf par 1
    
    # 4. Application du ratio à toute la série harmonisée
    r_corrected <- r_harmonized * r_ratio
    names(r_corrected) <- paste0("LAI_Corr_", common_dates)
    
    # Sauvegarde Finale (Corrigée)
    out_name_corr <- file.path(res_path, paste0(site, "_LAI_HarmonizedW_10d_LiDAR_Corrected_", distrib, ".tif"))
    # writeRaster(r_corrected, out_name_corr, overwrite = TRUE)
    
    cat("     [OK] Final saved:", basename(out_name_corr), "\n")
    
  } # End Loop Distrib
} # End Loop Site