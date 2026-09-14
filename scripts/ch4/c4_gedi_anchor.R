# ==============================================================================
# c4_gedi_anchor.R — THE Chapter 4 test: can a wall-to-wall LAI product anchored
# on GEDI (zero ALS in features AND in training target) recover the ALS truth?
#
# Design: RF(target = GEDI PAI at power-beam footprints; features = S2 LAI +
# FORMS-H, both available wall-to-wall) -> predict on the canonical pixel pool
# -> score against LAI_ALS (k = 0.65). Baselines: raw S2, GEDI-ratio-corrected
# S2, and the ALS-supervised LOSO RF of Ch3 (Table3b reference).
#
# Footprint filter (from c4 quality sweep): power == "full", sensitivity >= 0.90
# (Blois r jumps 0.49 -> 0.67 with power beams; sensitivity beyond 0.90 adds
# little once power-only). Pixel pool: same recipe/seed as c3_dopt_crosssite.R.
# Read-only on 01_DATA / 03_RESULTS. Run from NC_Full root.
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(data.table); library(randomForest); library(ggplot2)
  library(here)
})
source(here::here("scripts", "_article_style.R"))
sites   <- c("Aigoual", "Blois", "Mormal")
lai_min <- 2.0
n_samp  <- 5000L
k_scale <- 0.5 / 0.65
set.seed(42)
dir_tabs <- here::here("manuscripts/ch4", "tables")
dir_figs <- here::here("manuscripts/ch4", "figures")

# --- footprints: GEDI anchor pool + FORMS-H feature --------------------------
fp <- readRDS(here::here("output", "intermediate", "c4", "gedi_matchup_3site.rds"))
fp <- fp[power == "full" & !is.na(lai_s2)]
cat("anchor footprints (power beams):\n"); print(fp[, .N, by = site])

formsh <- rast(here::here("01_DATA", "FORMS-H_Height_10m_cm.tif"))  # EPSG:2154, cm
pts <- vect(fp[, .(x, y)], geom = c("x", "y"), crs = "EPSG:32631")
fp[, FORMS_H := terra::extract(formsh, project(pts, "EPSG:2154"))[, 2] / 100]
fp <- fp[is.finite(FORMS_H)]
cat("with FORMS-H:", nrow(fp), "\n")

# --- pixel evaluation pool (c3_dopt_crosssite recipe, no d_opt needed) -------
als_fn <- function(s) here::here("03_RESULTS", s, "Metrics", "Deciduous_Only",
                                 "ladstack_classic.tif")
s2_fn  <- function(s) here::here("output", "intermediate", "sm6", s,
                                 "s2lai_summer_atbd_T_res_10_m.tif")
build_site <- function(s) {
  lai_als <- sum(rast(als_fn(s)), na.rm = TRUE)
  lai_s2  <- rast(s2_fn(s))
  if (!compareGeom(lai_als, lai_s2, stopOnError = FALSE))
    lai_s2 <- resample(lai_s2, lai_als)
  ext154 <- project(as.polygons(ext(lai_als), crs = crs(lai_als)), "EPSG:2154")
  fh <- project(crop(formsh, ext154), lai_als) / 100
  d <- data.table(site = s,
                  LAI_ALS = as.numeric(values(lai_als)),
                  LAI_S2  = as.numeric(values(lai_s2)),
                  FORMS_H = as.numeric(values(fh)))
  d <- d[is.finite(LAI_ALS) & is.finite(LAI_S2) & is.finite(FORMS_H) &
           LAI_ALS >= lai_min]
  d[, LAI_ALS := LAI_ALS * k_scale]
  d[sample(.N, min(.N, n_samp))]
}
px <- rbindlist(lapply(sites, build_site))
cat("evaluation pixels:\n"); print(px[, .N, by = site])

# --- models ------------------------------------------------------------------
metr <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs)
  list(r2   = cor(pred[ok], obs[ok])^2,
       rmse = sqrt(mean((pred[ok] - obs[ok])^2)),
       bias = mean(pred[ok] - obs[ok]))
}
res <- rbindlist(lapply(sites, function(s) {
  te      <- px[site == s]
  fp_loc  <- fp[site == s]
  fp_oth  <- fp[site != s]
  px_oth  <- px[site != s]

  # GEDI-anchored ratio corrections (magnitude only, no ALS anywhere)
  ratio_loc  <- fp_loc[, mean(pai_gedi) / mean(lai_s2)]
  ratio_loso <- fp_oth[, mean(pai_gedi) / mean(lai_s2)]

  # GEDI-anchored RF (pattern + magnitude, no ALS anywhere)
  rf_gedi_loc  <- randomForest(x = fp_loc[, .(LAI_S2 = lai_s2, FORMS_H)],
                               y = fp_loc$pai_gedi, ntree = 400)
  rf_gedi_loso <- randomForest(x = fp_oth[, .(LAI_S2 = lai_s2, FORMS_H)],
                               y = fp_oth$pai_gedi, ntree = 400)

  # ALS-supervised LOSO reference (needs ALS at the training sites; Ch3 Table3b)
  rf_als_loso <- randomForest(x = px_oth[, .(LAI_S2, FORMS_H)],
                              y = px_oth$LAI_ALS, ntree = 400)

  preds <- list(
    RAW_S2          = te$LAI_S2,
    RATIO_GEDI_loc  = te$LAI_S2 * ratio_loc,
    RATIO_GEDI_loso = te$LAI_S2 * ratio_loso,
    RF_GEDI_loc     = predict(rf_gedi_loc,  te[, .(LAI_S2, FORMS_H)]),
    RF_GEDI_loso    = predict(rf_gedi_loso, te[, .(LAI_S2, FORMS_H)]),
    RF_ALS_loso     = predict(rf_als_loso,  te[, .(LAI_S2, FORMS_H)]))
  rbindlist(lapply(names(preds), function(m)
    c(list(site = s, model = m, n_train_fp = if (grepl("loc", m)) nrow(fp_loc)
           else if (grepl("GEDI_loso", m)) nrow(fp_oth) else NA_integer_),
      metr(preds[[m]], te$LAI_ALS))))
}))
res[, `:=`(r2 = round(r2, 3), rmse = round(rmse, 2), bias = round(bias, 2))]
cat("\n=== GEDI-anchored wall-to-wall recovery vs ALS truth (pixels) ===\n")
print(dcast(res, model ~ site, value.var = c("r2", "rmse", "bias")))
fwrite(res, file.path(dir_tabs, "Table5_c4_gedi_anchor.csv"))

# --- figure ------------------------------------------------------------------
lab <- c(RAW_S2 = "Raw S2", RATIO_GEDI_loc = "S2 × GEDI ratio (local)",
         RATIO_GEDI_loso = "S2 × GEDI ratio (LOSO)",
         RF_GEDI_loc = "RF GEDI-trained (local)",
         RF_GEDI_loso = "RF GEDI-trained (LOSO)",
         RF_ALS_loso = "RF ALS-supervised (LOSO ref)")
rl <- melt(res, id.vars = c("site", "model"), measure.vars = c("r2", "rmse"),
           variable.name = "metric")
rl[, model := factor(lab[model], levels = lab)]
rl[, metric := factor(metric, c("r2", "rmse"),
                      c("R² vs ALS LAI", "RMSE (m²/m²) vs ALS LAI"))]
p <- ggplot(rl, aes(model, value, fill = grepl("GEDI", model))) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.3, size = 2.8) +
  facet_grid(metric ~ site, scales = "free_y") +
  scale_fill_manual(values = c(`FALSE` = "grey55", `TRUE` = "#7570b3"),
                    guide = "none") +
  labs(x = NULL, y = NULL,
       title = "Wall-to-wall LAI anchored on GEDI (no ALS) vs the ALS truth",
       subtitle = paste("Purple = GEDI-anchored (ALS-free). Grey references:",
                        "raw S2 and the ALS-supervised LOSO RF.")) +
  theme_article() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7.5))
ggsave(file.path(dir_figs, "Fig4_c4_gedi_anchor.png"), p,
       width = 10, height = 6, dpi = 200)
cat("\nDONE\n")
