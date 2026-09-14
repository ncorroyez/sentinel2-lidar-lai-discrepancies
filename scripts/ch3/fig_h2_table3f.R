# ==============================================================================
# Chapter 3 — Fig H2 regenerated on the SAME sample as Table 3f, so its annotated
# R² matches the text (Blois, LAI ≥ 2, k_scale = 0.5/0.65, 5000-pixel uniform-LAI
# sample). Two panels: Sentinel-2 ATBD LAI vs full LiDAR LAI (R² ≈ 0.03) and vs
# top-d_opt LiDAR LAI (R² ≈ 0.27). Styled to _article_style.R.
# Run from NC_Full root:  Rscript scripts/fig_h2_table3f.R
# ==============================================================================
suppressPackageStartupMessages({ library(terra); library(data.table); library(ggplot2) })
NC <- "/home/corroyez/Documents/NC_Full"
source(file.path(NC, "scripts/_article_style.R"))
sites <- c("Aigoual","Blois","Mormal"); lai_min <- 2; k_scale <- 0.5/0.65
als <- function(s) file.path(NC,"03_RESULTS",s,"Metrics","Deciduous_Only","ladstack_classic.tif")
s2a <- function(s) file.path(NC,"output","intermediate","sm6",s,"s2lai_summer_atbd_T_res_10_m.tif")
dpc <- function(s) file.path(NC,"output","intermediate","lai_als_dopt",s,"LAI_ALS_dopt_common.tif")
# Build all three sites and sample exactly as Table3f (set.seed(42); by=site), so
# the Blois draw is byte-identical to the table.
full <- rbindlist(lapply(sites, function(s){
  g <- sum(rast(als(s)), na.rm=TRUE)
  ld <- function(r){ if(!compareGeom(g,r,stopOnError=FALSE)) r<-resample(r,g); as.numeric(values(r)) }
  data.table(site=s, LAI_full=as.numeric(values(g))*k_scale, LAI_dopt=ld(rast(dpc(s))), S2_ATBD=ld(rast(s2a(s))))[
    is.finite(LAI_full)&is.finite(S2_ATBD)&is.finite(LAI_dopt)&(LAI_full/k_scale)>=lai_min] }))
samp <- function(d,n){hi<-quantile(d$LAI_full,.98,names=FALSE);lo<-lai_min*k_scale
  br<-seq(lo,hi,length.out=8);d<-d[LAI_full<=hi];d[,b:=cut(LAI_full,br,include.lowest=TRUE)]
  per<-ceiling(n/nlevels(d$b));d[,.SD[sample(.N,min(.N,per))],by=b][,b:=NULL][]}
set.seed(42); Sall <- full[, samp(.SD,5000), by=site]; S <- Sall[site=="Blois"]
r2f <- cor(S$S2_ATBD,S$LAI_full)^2; r2d <- cor(S$S2_ATBD,S$LAI_dopt)^2
cat(sprintf("Blois (Table3f sample n=%d): S2_ATBD~full R²=%.3f | S2_ATBD~dopt7 R²=%.3f\n", nrow(S), r2f, r2d))

D <- rbindlist(list(
  data.table(target=sprintf("full LiDAR LAI (R² = %.2f)", r2f), x=S$LAI_full, y=S$S2_ATBD),
  data.table(target=sprintf("top-d_opt LiDAR LAI (R² = %.2f)", r2d), x=S$LAI_dopt, y=S$S2_ATBD)))
g2 <- ggplot(D, aes(x,y)) +
  geom_point(colour="#1A9850", size=0.5, alpha=0.25) +
  geom_smooth(method="lm", colour="black", linewidth=0.8, se=FALSE) +
  facet_wrap(~target, scales="free_x") +
  labs(x=expression("LiDAR LAI (one-sided,"~m^2/m^2*")"),
       y=expression(Sentinel*"-2 LAI ("*m^2/m^2*")")) +
  theme_article(11)
ggsave_article(file.path(NC,"outputs/figures_article_ch3/FigH2_s2_tracks_top"), g2, 6.9, 3.4)
cat("DONE -> FigH2_s2_tracks_top (Table3f sample)\n")
