# ==============================================================================
# Chapter 3 — RECONSTRUCTED durable generator for the cross-site d_opt transfer
# (Table 3 / Fig 4). Leave-one-site-out RF recovery of full LiDAR LAI (LAI_ALS),
# four feature sets, showing the ORACLE COLLAPSE:
#   S2  <  S2+FORMS-H  <  +d_opt(TRUE, ALS)  ,  but  +d_opt(PREDICTED, non-ALS) ≈ S2+FORMS-H
# "d_opt feature" = per-pixel LAI_ALS_DOPT (d_opt-truncated LiDAR LAI, common 7 m).
# "predicted d_opt" = RF(LAI_ALS_DOPT ~ S2 + FORMS-H) trained on the train sites
#   (a non-ALS reconstruction) — redundant with S2+FORMS-H by construction.
# Read-only on 01_DATA / 03_RESULTS. Run from NC_Full root.
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(data.table); library(randomForest); library(ggplot2); library(here)
})
sites   <- c("Aigoual", "Blois", "Mormal")
lai_min <- 2.0
n_samp  <- 5000L
set.seed(42)

als_fn  <- function(s) here("03_RESULTS", s, "Metrics", "Deciduous_Only", "ladstack_classic.tif")
s2_fn   <- function(s) here("output","intermediate","sm6", s, "s2lai_summer_atbd_T_res_10_m.tif")
dopt_fn <- function(s) here("output","intermediate","lai_als_dopt", s, "LAI_ALS_dopt_common.tif")
formsh  <- rast(here("01_DATA", "FORMS-H_Height_10m_cm.tif"))   # EPSG:2154, cm

k_scale <- 0.5 / 0.65        # full LiDAR LAI to study k=0.65 (R² is scale-invariant)
build_site <- function(s) {
  lai_als <- sum(rast(als_fn(s)), na.rm = TRUE)   # raw k=0.5 for the LAI>=2 pool; scaled below
  lai_s2  <- rast(s2_fn(s)); dopt <- rast(dopt_fn(s))
  if (!compareGeom(lai_als, lai_s2,  stopOnError=FALSE)) lai_s2 <- resample(lai_s2, lai_als)
  if (!compareGeom(lai_als, dopt,    stopOnError=FALSE)) dopt   <- resample(dopt,   lai_als)
  # FORMS-H: crop the France raster to this site (in Lambert-93), then project to grid
  ext154 <- project(as.polygons(ext(lai_als), crs = crs(lai_als)), "EPSG:2154")
  fh <- project(crop(formsh, ext154), lai_als) / 100      # cm -> m
  d <- data.table(site = s,
                  LAI_ALS   = as.numeric(values(lai_als)),
                  LAI_S2    = as.numeric(values(lai_s2)),
                  FORMS_H   = as.numeric(values(fh)),
                  LAI_DOPT  = as.numeric(values(dopt)))
  d <- d[is.finite(LAI_ALS) & is.finite(LAI_S2) & is.finite(FORMS_H) &
         is.finite(LAI_DOPT) & LAI_ALS >= lai_min]   # pool on raw k=0.5 LAI>=2
  d[, LAI_ALS := LAI_ALS * k_scale]                  # report/target at k=0.65 (R²-invariant)
  d[sample(.N, min(.N, n_samp))]
}
dt <- rbindlist(lapply(sites, build_site))
cat("sampled pixels per site:\n"); print(dt[, .N, by = site])

r2 <- function(pred, obs) { ok <- is.finite(pred) & is.finite(obs); if (sum(ok) < 5) return(NA);
  cor(pred[ok], obs[ok])^2 }
featsets <- list(
  S2                  = c("LAI_S2"),
  `S2+FORMSH`         = c("LAI_S2","FORMS_H"),
  `S2+FORMSH+dopt_TRUE`= c("LAI_S2","FORMS_H","LAI_DOPT"),
  `S2+FORMSH+dopt_PRED`= c("LAI_S2","FORMS_H","LAI_DOPT_PRED"))

res <- rbindlist(lapply(sites, function(test) {
  tr <- dt[site != test]; te <- copy(dt[site == test])
  # non-ALS predicted d_opt: RF(LAI_DOPT ~ S2 + FORMS_H) on train, predict on both
  m_pred <- randomForest(LAI_DOPT ~ LAI_S2 + FORMS_H, data = tr, ntree = 400)
  tr$LAI_DOPT_PRED <- predict(m_pred, tr); te$LAI_DOPT_PRED <- predict(m_pred, te)
  out <- data.table(site = test)
  for (fn in names(featsets)) {
    f  <- featsets[[fn]]
    m  <- randomForest(x = tr[, ..f], y = tr$LAI_ALS, ntree = 400)
    out[[fn]] <- r2(predict(m, te[, ..f]), te$LAI_ALS)
  }
  out
}))
cat("\n=== LOSO R² (predicting LAI_ALS on held-out site) ===\n")
print(res[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)], row.names = FALSE)
fwrite(res, here("output","tables","Table3b_dopt_crosssite_recovery.csv"))
fwrite(res, here("manuscripts/ch3","tables","Table3b_dopt_crosssite_recovery.csv"))

# ── Fig 4 regenerated from the reconstructed table ──────────────────────────
m <- melt(res, id.vars="site", variable.name="feat", value.name="r2")
lab <- c(S2="S2 only", `S2+FORMSH`="S2 + FORMS-H",
         `S2+FORMSH+dopt_TRUE`="S2 + FORMS-H + d_opt (TRUE, ALS)",
         `S2+FORMSH+dopt_PRED`="S2 + FORMS-H + d_opt (predicted, non-ALS)")
m[, feat := factor(lab[as.character(feat)], levels = rev(lab))]
cols <- c("S2 only"="#F08C7A","S2 + FORMS-H"="#1FA3A3",
          "S2 + FORMS-H + d_opt (TRUE, ALS)"="#9B5DE5",
          "S2 + FORMS-H + d_opt (predicted, non-ALS)"="#C9B8E8")
g <- ggplot(m, aes(r2, feat, fill = feat)) + geom_col(width=0.7) +
  geom_vline(xintercept=0, colour="grey60") + facet_wrap(~ site) +
  scale_fill_manual(values=cols, guide="none") +
  labs(title="Cross-site LAI transfer (leave-one-site-out)",
       subtitle="R² predicting LiDAR LAI on the held-out site. TRUE (ALS) d_opt lifts it; a satellite-predicted d_opt collapses back to baseline.",
       x=expression(R^2), y=NULL) +
  theme_minimal(base_size=11) + theme(plot.subtitle=element_text(size=8.5))
for (d in c(here("output","figures"), here("manuscripts/ch3","figures")))
  ggsave(file.path(d, "Fig4_dopt_transfer.png"), g, width=9.5, height=4, dpi=150)
cat("\nDONE\n")
