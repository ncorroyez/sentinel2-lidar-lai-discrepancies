# ==============================================================================
# Chapter 3 — cross-site d_opt recovery under DIFFERENT sampling schemes, to test
# whether the d_opt benefit strengthens with Ch2-style uniform-LAI sampling
# (LAI in [2, p98], bins of 1) and/or more training points — vs plain random.
# Same LOSO RF design as c3_dopt_crosssite.R. Target LAI_ALS at k=0.65.
# Run from NC_Full root:  Rscript scripts/ch3/c3_dopt_crosssite_sampling.R
# ==============================================================================
suppressPackageStartupMessages({ library(terra); library(data.table); library(randomForest); library(here) })
sites   <- c("Aigoual", "Blois", "Mormal")
lai_min <- 2.0; k_scale <- 0.5/0.65
als  <- function(s) here("03_RESULTS", s, "Metrics", "Deciduous_Only", "ladstack_classic.tif")
s2a  <- function(s) here("output","intermediate","sm6", s, "s2lai_summer_atbd_T_res_10_m.tif")
dpc  <- function(s) here("output","intermediate","lai_als_dopt", s, "LAI_ALS_dopt_common.tif")
formsh <- rast(here("01_DATA", "FORMS-H_Height_10m_cm.tif"))

# ---- full pool per site (all pixels passing the canonical filter) -----------
full <- rbindlist(lapply(sites, function(s) {
  laR <- sum(rast(als(s)), na.rm = TRUE); sa <- rast(s2a(s)); dp <- rast(dpc(s))
  for (nm in c("sa","dp")) if (!compareGeom(laR, get(nm), stopOnError=FALSE)) assign(nm, resample(get(nm), laR))
  ext154 <- project(as.polygons(ext(laR), crs = crs(laR)), "EPSG:2154")
  fh <- project(crop(formsh, ext154), laR) / 100
  d <- data.table(site = s, lai_raw = as.numeric(values(laR)),
                  LAI_S2 = as.numeric(values(sa)), FORMS_H = as.numeric(values(fh)),
                  LAI_DOPT = as.numeric(values(dp)))
  d[is.finite(lai_raw) & is.finite(LAI_S2) & is.finite(FORMS_H) & is.finite(LAI_DOPT) & lai_raw >= lai_min]
}))
cat("full pool per site:\n"); print(full[, .N, by = site])

# ---- samplers (return row indices within a per-site dt) ---------------------
samp_random  <- function(d, n) d[sample(.N, min(.N, n))]
samp_uniform <- function(d, n) {                 # Ch2-style: uniform over LAI in [2, p98]
  hi <- quantile(d$lai_raw, 0.98, names = FALSE)
  br <- seq(lai_min, hi, by = 1.0); if (tail(br,1) < hi) br <- c(br, hi)
  d <- d[lai_raw <= hi]; d[, bin := cut(lai_raw, br, include.lowest = TRUE)]
  per <- ceiling(n / nlevels(d$bin))
  d[, .SD[sample(.N, min(.N, per))], by = bin][, bin := NULL][]
}

featsets <- list(S2 = "LAI_S2", `S2+FH` = c("LAI_S2","FORMS_H"),
                 `+dopt_TRUE` = c("LAI_S2","FORMS_H","LAI_DOPT"),
                 `+dopt_PRED` = c("LAI_S2","FORMS_H","LAI_DOPT_PRED"))
r2 <- function(p,o){ok<-is.finite(p)&is.finite(o); if(sum(ok)<5) NA else cor(p[ok],o[ok])^2}

run_loso <- function(dt) {
  rbindlist(lapply(sites, function(test) {
    tr <- dt[site != test]; te <- copy(dt[site == test])
    mp <- randomForest(LAI_DOPT ~ LAI_S2 + FORMS_H, data = tr, ntree = 300)
    tr$LAI_DOPT_PRED <- predict(mp, tr); te$LAI_DOPT_PRED <- predict(mp, te)
    o <- data.table(site = test)
    for (fn in names(featsets)) { f <- featsets[[fn]]
      m <- randomForest(x = tr[, ..f], y = tr$LAI_ALS, ntree = 300)
      o[[fn]] <- r2(predict(m, te[, ..f]), te$LAI_ALS) }
    o }))
}

configs <- list(
  list(lab="random n=5000",       fn=function(d) samp_random(d,5000)),
  list(lab="uniform-LAI n=5000",  fn=function(d) samp_uniform(d,5000)),
  list(lab="uniform-LAI n=15000", fn=function(d) samp_uniform(d,15000)))

set.seed(42)
for (cf in configs) {
  dt <- full[, cf$fn(.SD), by = site]
  dt[, LAI_ALS := lai_raw * k_scale]
  res <- run_loso(dt)
  cat(sprintf("\n===== %-20s (n/site: %s) =====\n", cf$lab,
              paste(dt[, .N, by=site]$N, collapse="/")))
  print(res[, lapply(.SD, function(x) if(is.numeric(x)) round(x,3) else x)], row.names = FALSE)
  cat(sprintf("  d_opt LIFT (TRUE - S2+FH): %s | PRED-S2+FH: %s\n",
      paste(round(res$`+dopt_TRUE`-res$`S2+FH`,3),collapse="/"),
      paste(round(res$`+dopt_PRED`-res$`S2+FH`,3),collapse="/")))
}
cat("\nDONE\n")
