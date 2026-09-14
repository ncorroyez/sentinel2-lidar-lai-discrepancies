# ==============================================================================
# Chapter 3 — d_opt depth sweep with the OPERATIONAL substitution:
# for a pixel where canopy h_max < d_opt, the whole canopy is within the optical
# depth (S2 sees all of it) -> use LAI_S2 instead of LAI_ALS_dopt as the feature.
# This removes the circular "full-LiDAR-LAI injection" on short pixels and isolates
# the real d_opt signal (tall pixels, where d_opt actually truncates).
# Uniform-LAI n=5000 (speed). Compares vs the plain (no-substitution) feature.
# Run from NC_Full root:  Rscript scripts/ch3/c3_dopt_depth_hmax.R
# ==============================================================================
suppressPackageStartupMessages({ library(terra); library(data.table); library(randomForest); library(here) })
sites <- c("Aigoual","Blois","Mormal"); lai_min <- 2; k_scale <- 0.5/0.65
als <- function(s) here("03_RESULTS", s, "Metrics", "Deciduous_Only", "ladstack_classic.tif")
hmx <- function(s) here("03_RESULTS", s, "Metrics", "Deciduous_Only", "max_res_10_m.tif")
s2a <- function(s) here("output","intermediate","sm6", s, "s2lai_summer_atbd_T_res_10_m.tif")
fh0 <- rast(here("01_DATA","FORMS-H_Height_10m_cm.tif"))
pad_d <- function(s, d) here("03_RESULTS", s, "Metrics", "Not_Masked", "testPADs",
                             "PAD_Profiles_DSM_keepTrees", sprintf("PAD_%.1f_40.tif", 40 - d + 0.5))
d_persite <- c(Aigoual = 8, Blois = 6, Mormal = 10)
depths <- c(5, 7, 10, 15, 20)

full <- rbindlist(lapply(sites, function(s) {
  g <- sum(rast(als(s)), na.rm=TRUE)
  ld <- function(r){ if(!compareGeom(g,r,stopOnError=FALSE)) r<-resample(r,g); as.numeric(values(r)) }
  e <- project(as.polygons(ext(g), crs=crs(g)), "EPSG:2154")
  d <- data.table(site=s, lai_raw=as.numeric(values(g)), LAI_S2=ld(rast(s2a(s))),
                  FORMS_H=ld(project(crop(fh0,e), g))/100, Hmax=ld(rast(hmx(s))))
  for (dd in depths) d[[paste0("dop_",dd)]] <- ld(rast(pad_d(s,dd))) * k_scale
  d[["dop_ps"]] <- ld(rast(pad_d(s, d_persite[[s]]))) * k_scale
  d[["d_ps"]]   <- d_persite[[s]]
  d[is.finite(lai_raw)&is.finite(LAI_S2)&is.finite(FORMS_H)&is.finite(Hmax)&lai_raw>=lai_min]
}))
cat("pool:\n"); print(full[, .N, by=site])

samp_uniform <- function(d,n){hi<-quantile(d$lai_raw,.98,names=FALSE);br<-seq(lai_min,hi,by=1)
  if(tail(br,1)<hi)br<-c(br,hi);d<-d[lai_raw<=hi];d[,bin:=cut(lai_raw,br,include.lowest=TRUE)]
  per<-ceiling(n/nlevels(d$bin));d[,.SD[sample(.N,min(.N,per))],by=bin][,bin:=NULL][]}
r2<-function(p,o){ok<-is.finite(p)&is.finite(o);if(sum(ok)<5)NA else cor(p[ok],o[ok])^2}
set.seed(42); S <- full[, samp_uniform(.SD,5000), by=site]; S[, LAI_ALS := lai_raw*k_scale]

# share of pixels substituted (Hmax < d)
cat("\n% pixels Hmax<d (substituted to LAI_S2):\n")
for (dd in c(depths)) cat(sprintf("  d=%2d: %s\n", dd,
   paste(S[, sprintf("%s=%d%%", site, round(100*mean(Hmax<dd))), by=site]$V1, collapse=" ")))

loso <- function(dt, fc) rbindlist(lapply(sites, function(test){
  tr<-dt[site!=test];te<-dt[site==test]
  ok<-complete.cases(tr[,c("LAI_ALS",fc),with=FALSE]);oe<-complete.cases(te[,c("LAI_ALS",fc),with=FALSE])
  m<-randomForest(x=tr[ok,..fc],y=tr$LAI_ALS[ok],ntree=300)
  data.table(site=test,r2=r2(predict(m,te[oe,..fc]),te$LAI_ALS[oe]))}))

# build feature variants: plain dop(d), and substituted dop(d)|Hmax<d -> LAI_S2
res <- data.table(site=sites)
res$S2_FH <- loso(S, c("LAI_S2","FORMS_H"))$r2
for (dd in depths) {
  S[[paste0("plain_",dd)]] <- S[[paste0("dop_",dd)]]
  S[[paste0("subs_",dd)]]  <- fifelse(S$Hmax < dd, S$LAI_S2, S[[paste0("dop_",dd)]])
  res[[paste0("plain_d",dd)]] <- loso(S, c("LAI_S2","FORMS_H", paste0("plain_",dd)))$r2
  res[[paste0("subs_d",dd)]]  <- loso(S, c("LAI_S2","FORMS_H", paste0("subs_",dd)))$r2
}
# per_site
S$plain_ps <- S$dop_ps
S$subs_ps  <- fifelse(S$Hmax < S$d_ps, S$LAI_S2, S$dop_ps)
res$plain_ps <- loso(S, c("LAI_S2","FORMS_H","plain_ps"))$r2
res$subs_ps  <- loso(S, c("LAI_S2","FORMS_H","subs_ps"))$r2

cat("\n=== held-out R²: PLAIN dop(d) vs SUBSTITUTED (Hmax<d -> LAI_S2) ===\n")
print(res[, lapply(.SD, function(x) if(is.numeric(x)) round(x,3) else x)], row.names=FALSE)
fwrite(res, here("output","tables","Table3d_dopt_hmax_substitution.csv"))
fwrite(res, here("manuscripts/ch3","tables","Table3d_dopt_hmax_substitution.csv"))
cat("\nDONE\n")
