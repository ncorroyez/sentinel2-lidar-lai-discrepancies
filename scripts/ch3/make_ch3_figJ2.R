# ==============================================================================
# Fig J2 — split-threshold robustness. The dense/open reversal is not an artefact
# of the median: LiDAR wins the dense stratum and S2 wins the open stratum across
# a wide range of LAI split thresholds. Vertical line = the split actually used (3.86),
# which sits near the median LiDAR LAI of the 53 plots (median 3.80). The sweep now
# spans LAI 2.0-5.5, so it covers the physical buffering-regime breakpoint (observed
# slope & ΔTmax ~ LAI) = 2.89, marked on the figure. The reversal holds from 2.0 to 4.5.
# Caveat: below 2.4 the open stratum has only n=8 plots (LAI gap between 1.07 and 2.37,
# so thresholds 1.1-2.35 all give the same 8/45 partition); n_open is annotated per point.
# Data: manuscripts/ch3/tables/figJ_threshold_sweep.csv (c3_split_explore_wc.R).
#   Rscript scripts/ch3/make_ch3_figJ2.R
# ==============================================================================
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
source("scripts/_article_style.R")
TAB<-"manuscripts/ch3/tables"; FIG<-"manuscripts/ch3/figures"
W<-fread(file.path(TAB,"figJ_threshold_sweep.csv"))
Wm<-rbind(
  data.table(thr=W$thr,R2=W$LiDAR_dense_dT,series="LiDAR in dense stratum"),
  data.table(thr=W$thr,R2=W$S2_open_dT,   series="Sentinel-2 in open stratum"),
  data.table(thr=W$thr,R2=W$S2_dense_dT,  series="Sentinel-2 in dense stratum"),
  data.table(thr=W$thr,R2=W$LiDAR_open_dT,series="LiDAR in open stratum"))
Wm[,series:=factor(series,levels=c("LiDAR in dense stratum","Sentinel-2 in open stratum","Sentinel-2 in dense stratum","LiDAR in open stratum"))]
pal<-c("LiDAR in dense stratum"="#0072B2","Sentinel-2 in open stratum"="#D55E00",
       "Sentinel-2 in dense stratum"="#9ECAE1","LiDAR in open stratum"="#FDBE85")
pJ2<-ggplot(Wm,aes(thr,R2,colour=series))+
  geom_vline(xintercept=3.86,linetype="dashed",colour="grey45")+
  annotate("text",x=3.86,y=0.02,label="split used (3.86)",angle=90,vjust=-0.4,hjust=0,size=3,colour="grey35")+
  geom_vline(xintercept=2.89,linetype="dotted",colour="grey45")+
  annotate("text",x=2.89,y=0.02,label="buffering breakpoint (2.89)",angle=90,vjust=-0.4,hjust=0,size=3,colour="grey35")+
  geom_line(linewidth=0.9)+geom_point(size=1.6)+
  geom_text(data=W,aes(thr,0.86,label=n_open),inherit.aes=FALSE,size=2.5,colour="grey40")+
  geom_text(data=W,aes(thr,0.82,label=n_dense),inherit.aes=FALSE,size=2.5,colour="grey40")+
  annotate("text",x=min(W$thr)-0.15,y=0.86,label="n open",hjust=1,size=2.5,colour="grey40")+
  annotate("text",x=min(W$thr)-0.15,y=0.82,label="n dense",hjust=1,size=2.5,colour="grey40")+
  scale_colour_manual(values=pal,name=NULL)+
  scale_y_continuous(limits=c(0,0.88),breaks=seq(0,0.8,0.2))+
  coord_cartesian(clip="off")+
  labs(x=expression("LAI split threshold"~(m^2~m^-2)),y="Between-plot ΔTmax R²",
       title="The dense/open reversal is robust to the split threshold",
       subtitle="LiDAR wins the dense stratum and Sentinel-2 the open stratum across LAI 2.0–4.5, covering the\nbuffering breakpoint (2.89); below 2.4 the open stratum holds only n = 8 plots")+
  theme_article(12)+
  theme(legend.position="bottom",legend.direction="horizontal",
        plot.margin=margin(6,10,6,42))+
  guides(colour=guide_legend(nrow=2,byrow=TRUE))
ggsave_article(file.path(FIG,"FigJ2_threshold_sweep"),pJ2,8.2,5.0)
cat("DONE: FigJ2_threshold_sweep -> manuscripts/ch3/figures/\n")
