# ==============================================================================
# Chapter 3 — Fig "saturation" (S2 vs LiDAR LAI) in the ARTICLE style.
# Target journal: Remote Sensing of Environment (double-column 170 mm; vector PDF).
# Sources the shared style so theme/palette/sizing match every figure.
# Run from NC_Full root:  Rscript scripts/make_fig_saturation.R
# ==============================================================================
suppressPackageStartupMessages({ library(terra); library(data.table); library(ggplot2); library(here) })
source(here("scripts", "_article_style.R"))

sites <- c("Aigoual","Blois","Mormal"); lai_min <- 2.0
k_scale <- 0.5/0.65
als_fn <- function(s) here("03_RESULTS",s,"Metrics","Deciduous_Only","ladstack_classic.tif")
s2_fn  <- function(s) here("output","intermediate","sm6",s,"s2lai_summer_atbd_T_res_10_m.tif")
dt <- rbindlist(lapply(sites, function(s){
  lai_raw <- sum(rast(als_fn(s)), na.rm=TRUE); lai_s2 <- rast(s2_fn(s))
  if(!compareGeom(lai_raw,lai_s2,stopOnError=FALSE)) lai_s2 <- resample(lai_s2,lai_raw,method="bilinear")
  d <- data.table(site=s, lai_raw=as.numeric(values(lai_raw)), lai_s2=as.numeric(values(lai_s2)))
  d <- d[is.finite(lai_raw)&is.finite(lai_s2)&lai_raw>=lai_min]; d[,lai_als:=lai_raw*k_scale][,lai_raw:=NULL]; d }))
outdir <- here("outputs","figures_article_ch3"); dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
saveRDS(dt, file.path(outdir,"fig_saturation_data.rds"))                     # auditable underlying data

# goodness-of-fit (sensor intercomparison), aligned/identical format
r <- cor(dt$lai_s2, dt$lai_als); rmse <- sqrt(mean((dt$lai_s2-dt$lai_als)^2)); mae <- mean(abs(dt$lai_s2-dt$lai_als))
lab <- sprintf("italic(r)==%.2f", r); lab2 <- sprintf("RMSE==%.2f", rmse); lab3 <- sprintf("MAE==%.2f", mae)
set.seed(1); dpl <- dt[sample(.N, min(.N, 18000))]

# (a) 1:1 scatter, square, equal limits, dashed 1:1, GoF in fixed corner, no title
ga <- ggplot(dpl, aes(lai_als, lai_s2, colour=site)) +
  geom_point(alpha=0.06, size=0.4) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey40") +
  geom_smooth(method="lm", se=FALSE, linewidth=0.9) +
  annotate("text", x=0.2, y=c(7.6,7.1,6.6), hjust=0, size=3.2, parse=TRUE, label=c(lab,lab2,lab3)) +
  annotate("text", x=-Inf, y=Inf, hjust=-0.2, vjust=1.4, label="(a)", fontface="bold", size=4) +
  scale_colour_manual(values=PAL_SITE, name=NULL) +
  coord_equal(xlim=c(0,8), ylim=c(0,8), expand=FALSE) +
  labs(x=expression(LAI[LiDAR]~(m^2~m^-2)), y=expression(LAI[S2]~(m^2~m^-2))) +
  theme_article(12) + legend_corner(0.98,0.02) +
  guides(colour=guide_legend(override.aes=list(alpha=1,size=1.5)))

# (b) compressed distribution, no title
dlong <- rbind(data.table(source="LiDAR (ALS)", lai=dt$lai_als),
               data.table(source="Sentinel-2 (PROSAIL)", lai=dt$lai_s2))
gb <- ggplot(dlong, aes(lai, fill=source, colour=source)) +
  geom_density(alpha=0.35, linewidth=0.6) +
  annotate("text", x=-Inf, y=Inf, hjust=-0.2, vjust=1.4, label="(b)", fontface="bold", size=4) +
  scale_fill_manual(values=c("LiDAR (ALS)"="#E8746A","Sentinel-2 (PROSAIL)"="#1FA3A3"), name=NULL) +
  scale_colour_manual(values=c("LiDAR (ALS)"="#E8746A","Sentinel-2 (PROSAIL)"="#1FA3A3"), name=NULL) +
  labs(x=expression(LAI~(m^2~m^-2)), y="Density") +
  theme_article(12) + legend_corner(0.98,0.98)

library(patchwork); g <- ga + gb + plot_layout(widths=c(1,1.15))
ggsave_article(file.path(outdir,"Fig1_saturation"), g, width=9.0, height=4.4)
ggsave_article(here("manuscripts/ch3","figures","Fig1_saturation_styled"), g, width=9.0, height=4.4)
cat(sprintf("DONE  r=%.2f RMSE=%.2f MAE=%.2f  -> %s\n", r, rmse, mae, file.path(outdir,"Fig1_saturation.{png,pdf}")))
