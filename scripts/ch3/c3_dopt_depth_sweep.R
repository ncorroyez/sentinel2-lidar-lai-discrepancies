# ==============================================================================
# Chapter 3 — cross-site LiDAR-LAI recovery vs d_opt DEPTH choice.
# Tests how the d_opt feature (LAI integrated over the TOP d metres of canopy)
# affects held-out R²: shallow (what S2 senses) -> deep (-> full LiDAR LAI).
# Depths: 5, common(7), per_site(8/6/10), 10, 15, 20 m. Uniform-LAI n=15000.
# LAI_ALS_dopt(d) = scale_factor x PAD_(40-d+0.5)_40.tif  (k=0.65). Run from NC_Full root.
# ==============================================================================
suppressPackageStartupMessages({ library(terra); library(data.table); library(randomForest); library(here) })
sites <- c("Aigoual","Blois","Mormal"); lai_min <- 2; k_scale <- 0.5/0.65
als <- function(s) here("03_RESULTS", s, "Metrics", "Deciduous_Only", "ladstack_classic.tif")
s2a <- function(s) here("output","intermediate","sm6", s, "s2lai_summer_atbd_T_res_10_m.tif")
fh0 <- rast(here("01_DATA","FORMS-H_Height_10m_cm.tif"))
# Explicit PAD dir used by step 07 (top-d_opt-metres integration, k_ref=0.5).
pad_d <- function(s, d) here("03_RESULTS", s, "Metrics", "Not_Masked", "testPADs",
                             "PAD_Profiles_DSM_keepTrees", sprintf("PAD_%.1f_40.tif", 40 - d + 0.5))
d_persite <- c(Aigoual = 8, Blois = 6, Mormal = 10)   # per-site d_opt (Ch2)

# ---- base pool once (lai_raw, S2, FORMS_H) + all d_opt variants as columns ---
full <- rbindlist(lapply(sites, function(s) {
  laR <- sum(rast(als(s)), na.rm=TRUE); g <- laR
  ld <- function(r) { if (!compareGeom(g, r, stopOnError=FALSE)) r <- resample(r, g); as.numeric(values(r)) }
  sa <- rast(s2a(s)); e <- project(as.polygons(ext(g), crs=crs(g)), "EPSG:2154"); fh <- crop(fh0, e)
  d <- data.table(site=s, lai_raw=as.numeric(values(g)), LAI_S2=ld(sa),
                  FORMS_H=ld(project(fh, g))/100)
  d[, dop_5        := ld(rast(pad_d(s,5)))  * k_scale]
  d[, dop_common7  := ld(rast(pad_d(s,7)))  * k_scale]
  d[, dop_per_site := ld(rast(pad_d(s,d_persite[[s]]))) * k_scale]
  d[, dop_10       := ld(rast(pad_d(s,10))) * k_scale]
  d[, dop_15       := ld(rast(pad_d(s,15))) * k_scale]
  d[, dop_20       := ld(rast(pad_d(s,20))) * k_scale]
  d[is.finite(lai_raw) & is.finite(LAI_S2) & is.finite(FORMS_H) & lai_raw >= lai_min]
}))
cat("pool per site:\n"); print(full[, .N, by=site])
cat("mean LAI_dopt by depth (Blois):\n")
print(full[site=="Blois", lapply(.SD, mean, na.rm=TRUE),
           .SDcols=patterns("^dop_")][, lapply(.SD, round, 2)])

samp_uniform <- function(d, n) { hi <- quantile(d$lai_raw, .98, names=FALSE)
  br <- seq(lai_min, hi, by=1); if (tail(br,1)<hi) br <- c(br, hi)
  d <- d[lai_raw<=hi]; d[, bin:=cut(lai_raw, br, include.lowest=TRUE)]
  per <- ceiling(n/nlevels(d$bin)); d[, .SD[sample(.N, min(.N, per))], by=bin][, bin:=NULL][] }
r2 <- function(p,o){ok<-is.finite(p)&is.finite(o); if(sum(ok)<5) NA else cor(p[ok],o[ok])^2}

set.seed(42)
S <- full[, samp_uniform(.SD, 15000), by=site]; S[, LAI_ALS := lai_raw * k_scale]
dop_cols <- c("dop_5","dop_common7","dop_per_site","dop_10","dop_15","dop_20")

loso_feat <- function(dt, feats) rbindlist(lapply(sites, function(test) {
  tr <- dt[site!=test]; te <- dt[site==test]
  okt <- complete.cases(tr[, c("LAI_ALS", feats), with=FALSE])
  oke <- complete.cases(te[, c("LAI_ALS", feats), with=FALSE])
  m <- randomForest(x = tr[okt, ..feats], y = tr$LAI_ALS[okt], ntree=300)
  data.table(site=test, r2=r2(predict(m, te[oke, ..feats]), te$LAI_ALS[oke])) }))

base <- loso_feat(S, c("LAI_S2","FORMS_H"))[, .(site, S2_FORMSH=r2)]
out  <- base
for (dc in dop_cols) {
  rr <- loso_feat(S, c("LAI_S2","FORMS_H", dc))
  out[[dc]] <- rr$r2
}
cat("\n=== held-out R² predicting LiDAR LAI: S2+FORMS-H + d_opt(depth) ===\n")
print(out[, lapply(.SD, function(x) if(is.numeric(x)) round(x,3) else x)], row.names=FALSE)
cat("\n=== d_opt LIFT (depth - S2+FORMSH baseline) ===\n")
lift <- copy(out); for (dc in dop_cols) lift[[dc]] <- round(out[[dc]] - out$S2_FORMSH, 3)
print(lift[, c("site", dop_cols), with=FALSE], row.names=FALSE)
fwrite(out, here("output","tables","Table3c_dopt_depth_sweep.csv"))
fwrite(out, here("manuscripts/ch3","tables","Table3c_dopt_depth_sweep.csv"))
cat("\nDONE\n")
