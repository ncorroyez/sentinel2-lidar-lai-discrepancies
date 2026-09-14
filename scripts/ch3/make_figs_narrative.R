# Chapter 3 — Fig 3 / S1 / S2 in ARTICLE style. Target: RSE. Sources shared style.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(here) })
source(here("scripts","_article_style.R"))
outdir <- here("outputs","figures_article_ch3"); dir.create(outdir,showWarnings=FALSE,recursive=TRUE)
sv <- function(g,nm,w,h){ ggsave_article(file.path(outdir,nm),g,w,h)
  ggsave_article(here("manuscripts/ch3","figures",paste0(nm,"_styled")),g,w,h) }

# Fig 3 — magnitude recovery (Table6 LOO)
t6 <- fread(here("manuscripts/ch3","tables","Table6_LAIrecovery_LOO.csv"))[method!="dopt_lai_rescale"]
lab <- c(raw_s2="Raw S2", global_rescale="S2 x mean ratio (ratio-A)", rf_s2only="RF (S2 only)",
         dopt_physics="S2 x H/d_opt", rf_struct_only="RF (LiDAR structure)", rf_full="RF (S2 + LiDAR structure)")
t6[, lab:=reorder(lab[method], r2)]
t6[, grp:=ifelse(applicable_beyond_lidar,"S2-only (works beyond LiDAR)","needs LiDAR-derived inputs (oracle)")]
g3 <- ggplot(t6, aes(r2, lab, fill=grp)) + geom_col(width=0.7) +
  geom_text(aes(label=sprintf("%.2f",r2)), hjust=-0.15, size=3) +
  scale_fill_manual(values=c("S2-only (works beyond LiDAR)"="#E8746A","needs LiDAR-derived inputs (oracle)"="#9B5DE5"), name=NULL) +
  scale_x_continuous(limits=c(0,0.82), expand=c(0,0)) +
  labs(x=expression(R^2~"(predicting LiDAR LAI, leave-one-plot-out)"), y=NULL) +
  theme_article(11) + theme(legend.position="bottom")
sv(g3,"Fig3_magnitude_recovery",6.7,3.8)

# Fig S1 — annual S2 LAI series + smoother
ts <- as.data.table(readRDS(here("output","intermediate","blois_s2_lai_ts_2021.rds")))
cor <- as.data.table(readRDS(here("output","intermediate","blois_s2_lai_corrected_daily_2021.rds")))
ts_summ <- ts[, .(lai=mean(LAI_S2_ATBD,na.rm=TRUE)), by=doy]
gS1 <- ggplot() + geom_point(data=ts_summ, aes(doy,lai), colour="#1FA3A3", size=1.8) +
  geom_line(data=cor, aes(doy, LAI_S2_ATBD__raw), colour="#0B6E6E", linewidth=0.9) +
  labs(x="Day of year", y=expression(LAI[S2]~(m^2~m^-2))) + theme_article(11)
sv(gS1,"FigS1_annual_series",6.7,3.4)

# Fig S2 — per-plot corrected trajectories
pp <- as.data.table(readRDS(here("output","intermediate","blois_perplot_daily_v2_2021.rds")))
ppl <- melt(pp, id.vars=c("id","doy"), measure.vars=c("rescaled_ALS","rescaled_dopt","ML_rf"),
            variable.name="variant", value.name="lai")
vlab <- c(rescaled_ALS="x full LiDAR LAI", rescaled_dopt="x LiDAR LAI (d_opt)", ML_rf="RF (structure + S2)")
ppl[, variant:=vlab[as.character(variant)]]
gS2 <- ggplot(ppl, aes(doy,lai,group=id)) + geom_line(alpha=0.18, colour="#1A9850") + facet_wrap(~variant) +
  labs(x="Day of year", y=expression(LAI~(m^2~m^-2))) + theme_article(11)
sv(gS2,"FigS2_perplot_corrections",9,3.2)
cat("DONE narrative Fig3/S1/S2\n")
