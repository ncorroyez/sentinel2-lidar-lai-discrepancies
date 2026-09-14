# ==============================================================================
# Chapter 3 — IMPROVING the S2 LAI (explicit corrections), with the h_max-aware
# d_opt operationalisation: where canopy height < d_opt, S2 sees the whole canopy
# -> keep raw S2; where height >= d_opt, S2 misses the sub-d_opt layer -> correct.
# Evaluated vs LiDAR LAI (k=0.65), leave-one-site-out (calibrate on 2 sites, test
# on the 3rd). Height from FORMS-H (operational, non-ALS) vs LiDAR Hmax (oracle).
# d_opt = 7 m (operational common) and 10 m. Uniform-LAI n=5000.
# Run from NC_Full root:  Rscript scripts/ch3/c3_s2_improvement.R
# ==============================================================================
suppressPackageStartupMessages({ library(terra); library(data.table); library(here) })
sites <- c("Aigoual","Blois","Mormal"); lai_min <- 2; k_scale <- 0.5/0.65
als <- function(s) here("03_RESULTS", s, "Metrics", "Deciduous_Only", "ladstack_classic.tif")
hmx <- function(s) here("03_RESULTS", s, "Metrics", "Deciduous_Only", "max_res_10_m.tif")
s2a <- function(s) here("output","intermediate","sm6", s, "s2lai_summer_atbd_T_res_10_m.tif")
fh0 <- rast(here("01_DATA","FORMS-H_Height_10m_cm.tif"))

full <- rbindlist(lapply(sites, function(s) {
  g <- sum(rast(als(s)), na.rm=TRUE)
  ld <- function(r){ if(!compareGeom(g,r,stopOnError=FALSE)) r<-resample(r,g); as.numeric(values(r)) }
  e <- project(as.polygons(ext(g), crs=crs(g)), "EPSG:2154")
  data.table(site=s, LAI_ALS=as.numeric(values(g))*k_scale, LAI_S2=ld(rast(s2a(s))),
             FORMS_H=ld(project(crop(fh0,e), g))/100, Hmax=ld(rast(hmx(s))))[
    is.finite(LAI_ALS)&is.finite(LAI_S2)&is.finite(FORMS_H)&is.finite(Hmax)&(LAI_ALS/k_scale)>=lai_min]
}))
samp <- function(d,n){hi<-quantile(d$LAI_ALS,.98,names=FALSE);lo<-lai_min*k_scale
  br<-seq(lo,hi,length.out=8);d<-d[LAI_ALS<=hi];d[,b:=cut(LAI_ALS,br,include.lowest=TRUE)]
  per<-ceiling(n/nlevels(d$b));d[,.SD[sample(.N,min(.N,per))],by=b][,b:=NULL][]}
set.seed(42); S <- full[, samp(.SD,5000), by=site]

met <- function(pred, obs){ ok<-is.finite(pred)&is.finite(obs)
  c(r2=cor(pred[ok],obs[ok])^2, rmse=sqrt(mean((pred[ok]-obs[ok])^2)), bias=mean(pred[ok]-obs[ok])) }

# correction recipes: return predicted LAI given (train, test rows, d, height col)
recipes <- list(
  raw           = function(tr,te,d,H) te$LAI_S2,
  ratioA        = function(tr,te,d,H){ a<-mean(tr$LAI_ALS)/mean(tr$LAI_S2); te$LAI_S2*a },
  # plain physics: amplify by height/d_opt for ALL pixels (shrinks short ones)
  physics_plain = function(tr,te,d,H){ amp<-tr[[H]]/d; c<-mean(tr$LAI_ALS)/mean(tr$LAI_S2*amp)
                   (te[[H]]/d)*te$LAI_S2*c },
  # h_max-aware physics: keep raw S2 where height < d_opt, amplify only where >= d_opt
  physics_hawa  = function(tr,te,d,H){ amp<-ifelse(tr[[H]]>=d, tr[[H]]/d, 1); c<-mean(tr$LAI_ALS)/mean(tr$LAI_S2*amp)
                   ifelse(te[[H]]>=d, te[[H]]/d, 1)*te$LAI_S2*c },
  # ratio applied only where height >= d_opt (else keep raw S2)
  ratioA_hawa   = function(tr,te,d,H){ sel<-tr[[H]]>=d; a<-mean(tr$LAI_ALS[sel])/mean(tr$LAI_S2[sel])
                   ifelse(te[[H]]>=d, te$LAI_S2*a, te$LAI_S2) })

run <- function(d, H) {
  rbindlist(lapply(names(recipes), function(rn) {
    m <- rbindlist(lapply(sites, function(test){
      tr<-S[site!=test]; te<-S[site==test]; p<-recipes[[rn]](tr,te,d,H)
      as.data.table(as.list(met(p, te$LAI_ALS)))[, site:=test] }))
    data.table(recipe=rn, R2=mean(m$r2), RMSE=mean(m$rmse), Bias=mean(m$bias)) }))
}

for (cfg in list(list(d=7,H="FORMS_H"), list(d=10,H="FORMS_H"), list(d=15,H="FORMS_H"), list(d=10,H="Hmax"))) {
  cat(sprintf("\n===== d_opt=%dm | height=%s%s | %% pixels H<d: %s =====\n", cfg$d, cfg$H,
              if(cfg$H=="Hmax")" (LiDAR oracle)" else " (operational, FORMS-H)",
              paste(S[, sprintf("%s=%d%%", site, round(100*mean(get(cfg$H)<cfg$d))), by=site]$V1, collapse=" ")))
  out <- run(cfg$d, cfg$H)
  print(out[, lapply(.SD, function(x) if(is.numeric(x)) round(x,3) else x)], row.names=FALSE)
  fwrite(out, here("manuscripts/ch3","tables",
                   sprintf("Table3e_s2improve_d%d_%s.csv", cfg$d, cfg$H)))
}
cat("\nDONE\n")
