# ==============================================================================
# Chapter 3 — by-archetype (P1-P4) scores replacing the dense/open LAI-3.86 split.
#
# For the four forcing scenarios kept in the chapter (LiDAR fixed, S2 alone,
# S2 opt, LiDAR x S2 combination static + dynamic) plus the raw-bands random
# forest correction, score the time-matched summer daytime offset (dTmax) per
# canopy archetype P1-P4 (Chapter 1 clustering) and pooled over all 53 plots:
# between-plot R2 (with a plot-resampling bootstrap CI95), bias, RMSE, n.
#
# Inputs (all pre-computed; no MuSICA run):
#  - out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv  (rmusica repo)
#      archetype P + observed dTmax per plot (n = 53)
#  - out_files/Chapter3/tables/perplot_dtmax_conventions.csv (rmusica repo)
#      per-scenario per-plot dTmax; column d_new = time-matched convention
#      (written by c3_tm_extract.R)
#  - out_files/Chapter3/lai_prep/df_plots_real53.rds + the Blois 2021-06-14
#      L2A reflectance raster (NC_Full/03_RESULTS, read-only): features for the
#      LOO random-forest correction on the 10 raw S2 bands (recipe of
#      c3_tm_rf.R / c3_rf_rawbands.R: ntree = 500, mtry = 2, seed = 42; the
#      target here is the archetype-frame observation dTmax_obs).
#
# Output:
#  - manuscripts/ch3/tables/Table_archetype_scores.csv (NC_Full repo)
#
# No per-plot simulated coupling slope is persisted anywhere, so this table
# covers dTmax only. The slope would need re-extraction from the NetCDFs in
# out_files/Chapter3/nc_genuine53_windcorr/<SCENARIO>/HOBO_<id>.nc against the
# macro series of in_files/musica_in_Blois_pblh.nc (see c3_stratified_wc.R).
#
#   Rscript scripts/ch3/c3_archetype_scores.R
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(terra)
  library(randomForest)
})

# ---- path constants (inputs span two project roots, so no here::here()) ------
RMUSICA <- "/home/corroyez/Documents/z_Example_rmusica_31012025"
NC_FULL <- "/home/corroyez/Documents/NC_Full"
P_CLUSTER <- file.path(RMUSICA, "out_files/Chapter1/tables",
                       "tab_hobo_perplot_cluster.csv")
P_PERPLOT <- file.path(RMUSICA, "out_files/Chapter3/tables",
                       "perplot_dtmax_conventions.csv")
P_PLOTS   <- file.path(RMUSICA, "out_files/Chapter3/lai_prep",
                       "df_plots_real53.rds")
P_REFL    <- file.path(NC_FULL, "03_RESULTS/Blois",
                       "L2A_T31TCN_A031222_20210614T105443/Reflectance/res_10_m",
                       "L2A_T31TCN_A031222_20210614T105443_Refl")
P_OUT     <- file.path(NC_FULL, "manuscripts/ch3/tables",
                       "Table_archetype_scores.csv")

N_BOOT    <- 2000L
BOOT_SEED <- 123L
RF_SEED   <- 42L

# ---- chapter scenarios: label -> name in perplot_dtmax_conventions.csv -------
SCEN <- c(
  "LiDAR fixed"                 = "STATIC_ALS",
  "Sentinel-2 alone"            = "STATIC_S2_ATBD",
  "Sentinel-2 opt"              = "STATIC_S2_OPT",
  "LiDAR x S2 (static rescale)" = "STATIC_S2_RESCALED",
  "LiDAR x S2 (dynamic annual)" = "DYN_S2_ANNUAL",
  # Robustness variants of the LiDAR forcing, scored by archetype so that the
  # robustness appendix rests on the archetype reading rather than on the
  # retired two-stratum split.
  "LiDAR at k = 0.5"            = "STATIC_ALS_K05",
  "LiDAR at h_min = 3 m"        = "STATIC_ALS_HMIN3",
  "LiDAR at h_min = 5 m"        = "STATIC_ALS_HMIN5"
)

# ---- archetypes + observations -----------------------------------------------
cl <- fread(P_CLUSTER)
stopifnot(nrow(cl) == 53, all(c("id_plot", "dTmax_obs", "P") %in% names(cl)))
cnt <- cl[, .N, by = P][order(P)]
exp_cnt <- c(P1 = 8L, P2 = 12L, P3 = 13L, P4 = 20L)
stopifnot(identical(setNames(cnt$N, cnt$P), exp_cnt))
cat("archetype counts OK:", paste(cnt$P, cnt$N, collapse = " | "), "\n")

# ---- per-plot simulated dTmax (time-matched: d_new) --------------------------
pp <- fread(P_PERPLOT)
missing_scen <- setdiff(SCEN, unique(pp$scenario))
if (length(missing_scen))
  stop("scenarios missing from per-plot file: ", paste(missing_scen, collapse = ", "))
sim <- pp[scenario %in% SCEN, .(scenario, id_plot, sim = d_new)]
sim[, label := names(SCEN)[match(scenario, SCEN)]]

# ---- LOO random forest on the 10 raw S2 bands (correction model) -------------
# Recipe follows c3_tm_rf.R (ntree = 500, mtry = 2, seed = 42, LOO over the 53
# plots); the target is dTmax_obs so the correction is scored in the same
# frame as the four forcing scenarios.
rf_row <- NULL
if (file.exists(P_PLOTS) && file.exists(paste0(P_REFL, ".hdr"))) {
  df <- as.data.table(readRDS(P_PLOTS))
  r  <- rast(P_REFL)
  names(r) <- c("B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A",
                "B11", "B12")
  pts <- vect(as.data.frame(df[, .(x, y)]), geom = c("x", "y"),
              crs = "EPSG:32631")
  X <- as.data.frame(as.data.table(terra::extract(r, pts))[, -1])
  M <- merge(cl[, .(id_plot, y = dTmax_obs)],
             cbind(df[, .(id_plot)], X), by = "id_plot")
  stopifnot(nrow(M) == 53)
  Xm <- as.data.frame(M[, -c("id_plot", "y")])
  set.seed(RF_SEED)
  pred <- rep(NA_real_, nrow(M))
  for (i in seq_len(nrow(M))) {
    rf <- randomForest(x = Xm[-i, , drop = FALSE], y = M$y[-i],
                       ntree = 500, mtry = 2)
    pred[i] <- as.numeric(predict(rf, Xm[i, , drop = FALSE]))
  }
  rf_row <- data.table(scenario = "RF_RAWBANDS_LOO", id_plot = M$id_plot,
                       sim = pred,
                       label = "Correction (RF raw S2 bands)")
  cat("RF correction: LOO predictions computed (seed", RF_SEED, ")\n")
} else {
  cat("RF inputs not found -> correction row skipped (NA)\n")
}
sim <- rbind(sim, rf_row, fill = TRUE)

# ---- merge with archetypes ---------------------------------------------------
D <- merge(sim, cl[, .(id_plot, obs = dTmax_obs, P)], by = "id_plot")

# ---- scoring: R2 (+ bootstrap CI95), bias, RMSE, n ---------------------------
score <- function(obs, sim_v) {
  ok <- is.finite(obs) & is.finite(sim_v)
  o <- obs[ok]; s <- sim_v[ok]; n <- length(o)
  if (n < 3)
    return(list(n = n, R2 = NA_real_, R2_lo = NA_real_, R2_hi = NA_real_,
                bias = NA_real_, RMSE = NA_real_))
  r2 <- cor(o, s)^2
  set.seed(BOOT_SEED)
  bs <- replicate(N_BOOT, {
    i <- sample.int(n, n, replace = TRUE)
    if (sd(o[i]) == 0 || sd(s[i]) == 0) NA_real_ else cor(o[i], s[i])^2
  })
  ci <- unname(quantile(bs, c(.025, .975), na.rm = TRUE))
  list(n = n, R2 = r2, R2_lo = ci[1], R2_hi = ci[2],
       bias = mean(s - o), RMSE = sqrt(mean((s - o)^2)))
}

by_arch <- D[, score(obs, sim), by = .(label, scenario, group = P)]
pooled  <- D[, score(obs, sim), by = .(label, scenario)][, group := "All (53)"]
res <- rbind(by_arch, pooled)
lev <- c(names(SCEN), "Correction (RF raw S2 bands)")
res[, label := factor(label, levels = lev)]
setorder(res, label, group)
setcolorder(res, c("label", "scenario", "group", "n", "R2", "R2_lo", "R2_hi",
                   "bias", "RMSE"))

fwrite(res, P_OUT)
cat("wrote", P_OUT, sprintf("(%d rows)\n", nrow(res)))
print(res[, .(label, group, n, R2 = round(R2, 3),
              CI = sprintf("[%.2f, %.2f]", R2_lo, R2_hi),
              bias = round(bias, 2), RMSE = round(RMSE, 2))], nrow = 100)
