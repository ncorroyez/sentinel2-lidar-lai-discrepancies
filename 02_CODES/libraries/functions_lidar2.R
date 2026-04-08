# ---
# title: "functions_lidar2.R"
# ---

library("lidR")

# Calculate Plant Area Index (PAI) / LAI with a specific starting height
myPAI <- function(z, zmin = 1, k = 0.5) {
  Nout <- sum(z < zmin)
  Nin <- length(z)
  if (Nin == 0) return(NA)
  
  PAI <- -log(Nout / Nin) / k
  return(PAI)
}

# Calculate Plant Area Density (PAD) Layers starting from z0 = 1
myPAD_z1 <- function(z) {
  # Calculate LAD with a base height of 1m
  lad_i <- LAD(z, z0 = 1)
  lad_i$z <- floor(lad_i$z) + 0.5 
  
  # Ensure all targeted vertical bins are included
  new_z <- seq(1.5, 39.5)
  missing_z <- setdiff(new_z, lad_i$z)
  missing_data <- data.frame(z = missing_z, lad = rep(NA, length(missing_z)))
  
  lad <- rbind(lad_i, missing_data)
  lad <- lad[order(lad$z), ]
  row.names(lad) <- NULL
  
  # Format as a list of numeric values for pixel_metrics
  LADs <- numeric(length(new_z))
  k <- 1
  for (i in new_z) {
    if (!is.na(lad$lad[lad$z == i])) {
      LADs[k] <- lad$lad[lad$z == i]
    } else {
      LADs[k] <- NA
    }
    k <- k + 1
  }
  
  LADs_list <- setNames(as.list(LADs), paste0("LAD_Layer_", new_z))
  return(LADs_list)
}