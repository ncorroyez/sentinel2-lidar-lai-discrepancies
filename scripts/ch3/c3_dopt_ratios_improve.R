# ==============================================================================
# Chapter 3 — do the d_opt-based pairings / opt LUT improve the S2 correction?
#  (1) Correlation matrix: {S2_ATBD, S2_opt} x {LAI_ALS_full, LAI_ALS_dopt(7m)}
#      — does S2 track the TOP it senses (LAI_ALS_dopt) better than the FULL LAI?
#  (2) Correction recipes vs full LiDAR LAI (LOSO, k=0.65):
#      ratioA (full/ATBD), ratioD (full/opt), raw opt, and the d_opt two-step
#      (S2 -> top via b, top -> full via g) which should collapse to ratioA.
# Uniform-LAI n=5000. Run from NC_Full root:  Rscript scripts/ch3/c3_dopt_ratios_improve.R
# ==============================================================================
suppressPackageStartupMessages({ library(terra); library(data.table); library(here) })
sites <- c("Aigoual","Blois","Mormal"); lai_min <- 2; k_scale <- 0.5/0.65
als  <- function(s) here("03_RESULTS", s, "Metrics", "Deciduous_Only", "ladstack_classic.tif")
hmx  <- function(s) here("03_RESULTS", s, "Metrics", "Deciduous_Only", "max_res_10_m.tif")
s2a  <- function(s) here("output","intermediate","sm6", s, "s2lai_summer_atbd_T_res_10_m.tif")
s2o  <- function(s) here("output","intermediate","sm6", s, "s2lai_summer_opt_res_10_m.tif")
dpc  <- function(s) here("output","intermediate","lai_als_dopt", s, "LAI_ALS_dopt_common.tif")

full <- rbindlist(lapply(sites, function(s) {
  g <- sum(rast(als(s)), na.rm=TRUE)
  ld <- function(r){ if(!compareGeom(g,r,stopOnError=FALSE)) r<-resample(r,g); as.numeric(values(r)) }
  data.table(site=s, LAI_full=as.numeric(values(g))*k_scale, LAI_dopt=ld(rast(dpc(s))),
             S2_ATBD=ld(rast(s2a(s))), S2_opt=ld(rast(s2o(s))), Hmax=ld(rast(hmx(s))))[
    is.finite(LAI_full)&is.finite(S2_ATBD)&is.finite(S2_opt)&is.finite(LAI_dopt)&(LAI_full/k_scale)>=lai_min]
}))
samp <- function(d,n){hi<-quantile(d$LAI_full,.98,names=FALSE);lo<-lai_min*k_scale
  br<-seq(lo,hi,length.out=8);d<-d[LAI_full<=hi];d[,b:=cut(LAI_full,br,include.lowest=TRUE)]
  per<-ceiling(n/nlevels(d$b));d[,.SD[sample(.N,min(.N,per))],by=b][,b:=NULL][]}
set.seed(42); S <- full[, samp(.SD,5000), by=site]

# ---- (1) correlation matrix (R²) ----
cat("=== R² of S2 vs LiDAR targets (does S2 track the TOP it senses?) ===\n")
cm <- S[, .(`S2_ATBD~full`=cor(S2_ATBD,LAI_full)^2, `S2_ATBD~dopt7`=cor(S2_ATBD,LAI_dopt)^2,
            `S2_opt~full`=cor(S2_opt,LAI_full)^2,   `S2_opt~dopt7`=cor(S2_opt,LAI_dopt)^2), by=site]
print(cm[, lapply(.SD, function(x) if(is.numeric(x)) round(x,3) else x)], row.names=FALSE)
cat(sprintf("POOLED: S2_ATBD~full=%.3f S2_ATBD~dopt=%.3f | S2_opt~full=%.3f S2_opt~dopt=%.3f\n",
  cor(S$S2_ATBD,S$LAI_full)^2, cor(S$S2_ATBD,S$LAI_dopt)^2,
  cor(S$S2_opt,S$LAI_full)^2, cor(S$S2_opt,S$LAI_dopt)^2))

# ---- (2) correction recipes vs FULL LAI, LOSO ----
met <- function(p,o){ok<-is.finite(p)&is.finite(o);c(r2=cor(p[ok],o[ok])^2,
  rmse=sqrt(mean((p[ok]-o[ok])^2)),bias=mean(p[ok]-o[ok]))}
recipes <- list(
  raw_ATBD     = function(tr,te) te$S2_ATBD,
  ratioA_ATBD  = function(tr,te){ a<-mean(tr$LAI_full)/mean(tr$S2_ATBD); te$S2_ATBD*a },
  raw_opt      = function(tr,te) te$S2_opt,
  ratioD_opt   = function(tr,te){ a<-mean(tr$LAI_full)/mean(tr$S2_opt); te$S2_opt*a },
  # d_opt two-step: S2 -> top (b), top -> full (g). b*g = ratioA (should match).
  dopt_2step   = function(tr,te){ b<-mean(tr$LAI_dopt)/mean(tr$S2_ATBD)
                   g<-mean(tr$LAI_full)/mean(tr$LAI_dopt); te$S2_ATBD*b*g },
  # opt paired to the TOP (consistent pair), then scaled to full by structural g
  opt_dopt_2step = function(tr,te){ b<-mean(tr$LAI_dopt)/mean(tr$S2_opt)
                   g<-mean(tr$LAI_full)/mean(tr$LAI_dopt); te$S2_opt*b*g })
cat("\n=== correction recipes vs FULL LiDAR LAI (LOSO mean over 3 held-out sites) ===\n")
out <- rbindlist(lapply(names(recipes), function(rn){
  m<-rbindlist(lapply(sites, function(test){tr<-S[site!=test];te<-S[site==test]
    as.data.table(as.list(met(recipes[[rn]](tr,te), te$LAI_full)))}))
  data.table(recipe=rn, R2=round(mean(m$r2),3), RMSE=round(mean(m$rmse),3), Bias=round(mean(m$bias),3))}))
print(out, row.names=FALSE)
fwrite(cm, here("manuscripts/ch3","tables","Table3f_S2_vs_dopt_correlations.csv"))
fwrite(out, here("manuscripts/ch3","tables","Table3f_dopt_ratio_recipes.csv"))
cat("\nDONE\n")
