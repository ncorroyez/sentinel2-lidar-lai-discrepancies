# ==============================================================================
# Chapter 3 — LAYERED HYBRID LAI: give each sensor the layer it measures best.
#   tall pixel (Hmax > d_opt): hybrid = S2(top, calibrated) + LiDAR(bottom, exact)
#       LAI_below = LAI_ALS_full - LAI_ALS_dopt(d)   (ground -> Hmax-d, from LiDAR)
#       S2 calibrated to the top scale: S2 x b_top  (b_top = mean(LAI_dopt)/mean(S2))
#   short pixel (Hmax <= d_opt): S2 sees the whole canopy -> hybrid = S2 x b_full
# Evaluated vs full LiDAR LAI (k=0.65), leave-one-site-out, depths 7/10/15 m.
# (Diagnostic/fusion: the bottom uses LiDAR, so NOT operational — it tests how well
#  S2's top + LiDAR's bottom reconstruct the total, i.e. the d_opt decomposition.)
# Run from NC_Full root:  Rscript scripts/ch3/c3_hybrid_lai.R
# ==============================================================================
suppressPackageStartupMessages({ library(terra); library(data.table); library(here) })
sites <- c("Aigoual","Blois","Mormal"); lai_min <- 2; k_scale <- 0.5/0.65
als <- function(s) here("03_RESULTS", s, "Metrics", "Deciduous_Only", "ladstack_classic.tif")
hmx <- function(s) here("03_RESULTS", s, "Metrics", "Deciduous_Only", "max_res_10_m.tif")
s2a <- function(s) here("output","intermediate","sm6", s, "s2lai_summer_atbd_T_res_10_m.tif")
pad_d <- function(s, d) here("03_RESULTS", s, "Metrics", "Not_Masked", "testPADs",
                            "PAD_Profiles_DSM_keepTrees", sprintf("PAD_%.1f_40.tif", 40 - d + 0.5))
depths <- c(7, 10, 15)

full <- rbindlist(lapply(sites, function(s) {
  g <- sum(rast(als(s)), na.rm=TRUE)
  ld <- function(r){ if(!compareGeom(g,r,stopOnError=FALSE)) r<-resample(r,g); as.numeric(values(r)) }
  d <- data.table(site=s, LAI_full=as.numeric(values(g))*k_scale, S2=ld(rast(s2a(s))), Hmax=ld(rast(hmx(s))))
  for (dd in depths) d[[paste0("dopt_",dd)]] <- ld(rast(pad_d(s,dd))) * k_scale
  d[is.finite(LAI_full)&is.finite(S2)&is.finite(Hmax)&(LAI_full/k_scale)>=lai_min]
}))
samp <- function(d,n){hi<-quantile(d$LAI_full,.98,names=FALSE);lo<-lai_min*k_scale
  br<-seq(lo,hi,length.out=8);d<-d[LAI_full<=hi];d[,b:=cut(LAI_full,br,include.lowest=TRUE)]
  per<-ceiling(n/nlevels(d$b));d[,.SD[sample(.N,min(.N,per))],by=b][,b:=NULL][]}
set.seed(42); S <- full[, samp(.SD,5000), by=site]
met <- function(p,o){ok<-is.finite(p)&is.finite(o);c(r2=cor(p[ok],o[ok])^2,
  rmse=sqrt(mean((p[ok]-o[ok])^2)),bias=mean(p[ok]-o[ok]))}

baseline <- rbindlist(lapply(c("raw_S2","ratioA"), function(rn) {
  m<-rbindlist(lapply(sites,function(test){tr<-S[site!=test];te<-S[site==test]
    p<-if(rn=="raw_S2") te$S2 else te$S2*(mean(tr$LAI_full)/mean(tr$S2))
    as.data.table(as.list(met(p,te$LAI_full)))}))
  data.table(method=rn, d=NA, R2=round(mean(m$r2),3),RMSE=round(mean(m$rmse),3),Bias=round(mean(m$bias),3))}))

hyb <- rbindlist(lapply(depths, function(dd) {
  dc <- paste0("dopt_",dd)
  m <- rbindlist(lapply(sites, function(test){
    tr<-S[site!=test]; te<-copy(S[site==test])
    tall_tr <- tr$Hmax>dd; short_tr <- !tall_tr
    b_top  <- mean(tr[[dc]][tall_tr])/mean(tr$S2[tall_tr])
    b_full <- if(any(short_tr)) mean(tr$LAI_full[short_tr])/mean(tr$S2[short_tr]) else mean(tr$LAI_full)/mean(tr$S2)
    below  <- pmax(te$LAI_full - te[[dc]], 0)
    p <- ifelse(te$Hmax>dd, te$S2*b_top + below, te$S2*b_full)
    as.data.table(as.list(met(p, te$LAI_full)))[, site:=test] }))
  data.table(method="hybrid_layered", d=dd,
             R2=round(mean(m$r2),3), RMSE=round(mean(m$rmse),3), Bias=round(mean(m$bias),3),
             pct_tall=round(100*mean(S$Hmax>dd))) }), fill=TRUE)

cat("=== baseline corrections vs full LiDAR LAI (LOSO) ===\n"); print(baseline, row.names=FALSE)
cat("\n=== LAYERED HYBRID (S2 top + LiDAR bottom) vs full LiDAR LAI (LOSO) ===\n")
print(hyb, row.names=FALSE)
# per-site for the best depth
cat("\n=== hybrid per-site (d=10) ===\n")
ps <- rbindlist(lapply(sites, function(test){tr<-S[site!=test];te<-copy(S[site==test])
  tall<-tr$Hmax>10; b_top<-mean(tr$dopt_10[tall])/mean(tr$S2[tall])
  b_full<-mean(tr$LAI_full[!tall])/mean(tr$S2[!tall]); below<-pmax(te$LAI_full-te$dopt_10,0)
  p<-ifelse(te$Hmax>10, te$S2*b_top+below, te$S2*b_full)
  as.data.table(as.list(met(p,te$LAI_full)))[,site:=test]}))
print(ps[, .(site, r2=round(r2,3), rmse=round(rmse,3), bias=round(bias,3))], row.names=FALSE)
fwrite(rbind(baseline, hyb, fill=TRUE), here("manuscripts/ch3","tables","Table3g_hybrid_layered_lai.csv"))
cat("\nDONE\n")
