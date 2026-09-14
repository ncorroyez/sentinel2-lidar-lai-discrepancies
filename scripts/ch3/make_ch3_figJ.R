# ==============================================================================
# Fig J — stratified diagnostic (§3.4 centerpiece, windcorr genuine-53).
# The pooled Sentinel-2 R2 advantage is an OPEN-canopy effect; in dense (saturated)
# canopy the LiDAR wins the between-plot ranking. Split at median LiDAR LAI.
# Data: manuscripts/ch3/tables/figJ_stratified.csv (c3_stratified_wc.R).
#   Rscript scripts/ch3/make_ch3_figJ.R   (from NC_Full root)
# ==============================================================================
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
source("scripts/_article_style.R")
source("scripts/_figure_label_check.R")
TAB<-"manuscripts/ch3/tables"; FIG<-"manuscripts/ch3/figures"
J<-fread(file.path(TAB,"figJ_stratified.csv"))
Jm<-melt(J,id.vars=c("stratum","sensor","n"),measure.vars=c("dT_R2","sl_R2"),variable.name="metric",value.name="R2")
Jm[,metric:=factor(metric,levels=c("dT_R2","sl_R2"),labels=c("Amplitude ΔTmax (R²)","Coupling slope (R²)"))]
Jm[,stratum:=factor(stratum,levels=c("Open (LAI<med)","Dense (LAI>=med)","Pooled (all 53)"),
   labels=c("Open\ncanopy","Dense\ncanopy","Pooled\n(all 53)"))]
Jm[,sensor:=factor(sensor,levels=c("LiDAR","Sentinel-2"))]
check_labels(Jm$R2, "manuscripts/ch3/Chapter3_article_standalone_EN.md", what="Fig J — R2 stratifies")

pJ<-ggplot(Jm,aes(stratum,R2,fill=sensor))+
  geom_col(width=0.7,position=position_dodge(0.72))+
  geom_text(aes(label=fmt2(R2)),position=position_dodge(0.72),vjust=-0.35,size=2.9,colour="grey25")+
  facet_wrap(~metric)+scale_y_continuous(limits=c(0,0.82),expand=expansion(mult=c(0,0.04)))+
  scale_fill_manual(values=c("LiDAR"=PAL_SENSOR[["LiDAR"]],"Sentinel-2"=PAL_SENSOR[["Sentinel-2"]]),name=NULL)+
  labs(x=NULL,y="Between-plot R² (sim vs observed)",
       title="The Sentinel-2 R² advantage is an open-canopy artefact",
       subtitle="In dense (saturated) canopy the LiDAR wins the ranking; pooling hides the reversal (split at LAI 3.86, near the median)")+
  theme_article(12)+theme(legend.position="top")
ggsave_article(file.path(FIG,"FigJ_stratified"),pJ,8.6,4.4)
cat("DONE: FigJ_stratified -> manuscripts/ch3/figures/\n")
