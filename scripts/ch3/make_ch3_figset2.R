# ==============================================================================
# Chapter 3 additional figures E/F/G/H/I (v3.2.3, genuine-53).
#   E premise (buffering ~ LiDAR LAI, obs & MuSICA).  F why S2 fails.
#   G operational-vs-oracle ablation.  H validation scatter (amplitude collapse).
#   I seasonal windows (dynamic does not help at shoulders).
#   Rscript scripts/ch3/make_ch3_figset2.R   (from NC_Full root)
# ==============================================================================
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(scales)})
source("scripts/_article_style.R")
TAB<-"manuscripts/ch3/tables"; FIG<-"manuscripts/ch3/figures"
FAM<-c("LiDAR"="#0072B2","Sentinel-2"="#D55E00","Fusion"="#CC79A7","d_opt truncated"="#666666")

## ---- Fig E : premise -------------------------------------------------------
E<-fread(file.path(TAB,"figE_premise.csv"))
Em<-rbind(
  data.table(LAI=E$LAI_ALS,val=E$do,metric="Daytime offset ΔTmax (°C)",src="Observed (HOBO)"),
  data.table(LAI=E$LAI_ALS,val=E$sim_do,metric="Daytime offset ΔTmax (°C)",src="MuSICA (LiDAR LAI)"),
  data.table(LAI=E$LAI_ALS,val=E$so,metric="Coupling slope",src="Observed (HOBO)"),
  data.table(LAI=E$LAI_ALS,val=E$sim_so,metric="Coupling slope",src="MuSICA (LiDAR LAI)"))
rlab<-E[,.(r_obs=cor(LAI_ALS,do),r_sim=cor(LAI_ALS,sim_do))]
pE<-ggplot(Em,aes(LAI,val,colour=src))+geom_point(size=1.8,alpha=0.8)+
  geom_smooth(method="lm",se=FALSE,linewidth=0.8)+facet_wrap(~metric,scales="free_y")+
  scale_colour_manual(values=c("Observed (HOBO)"="grey20","MuSICA (LiDAR LAI)"="#0072B2"),name=NULL)+
  labs(x=expression("LiDAR LAI"~(m^2~m^-2)),y=NULL,
       title="Understory buffering scales with canopy LAI — and MuSICA reproduces it",
       subtitle=sprintf("Observed r = %.2f; MuSICA r = %.2f (ΔTmax)",rlab$r_obs,rlab$r_sim))+
  theme_article(12)+legend_corner(0.99,0.99)
ggsave_article(file.path(FIG,"FigE_premise"),pE,8.4,4.4)

## ---- Fig F : why S2 fails ---------------------------------------------------
F<-fread(file.path(TAB,"figF_s2_vs_lidar.csv"))
pF<-ggplot(F,aes(LAI_ALS,S2opt,colour=fCover))+
  geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey55")+
  geom_point(size=3)+scale_colour_viridis_c(option="D",name="fCover")+
  annotate("text",x=0.5,y=6.4,label="S2 over-reads open plots",hjust=0,size=3.3,colour="grey25")+
  annotate("text",x=5.6,y=2.4,label="S2 saturates in\ndense canopy",hjust=1,size=3.3,colour="grey25")+
  labs(x=expression("LiDAR LAI (reference)"~(m^2~m^-2)),y=expression("Sentinel-2 optical LAI"~(m^2~m^-2)),
       title="Sentinel-2 fails at both ends of the LAI range",
       subtitle="Over-estimates open/sparse canopy, saturates in dense canopy (dashed = 1:1)")+
  theme_article(12)+legend_corner(0.99,0.02)
ggsave_article(file.path(FIG,"FigF_s2_vs_lidar"),pF,6.8,4.8)

## ---- Fig G : operational-vs-oracle ablation --------------------------------
G<-fread(file.path(TAB,"figG_ablation.csv"))
G[,tier:=fifelse(feature_set=="LiDAR LAI (direct)","ceiling",fifelse(operational,"operational (no LiDAR)","oracle (needs LiDAR)"))]
G[,feature_set:=factor(feature_set,levels=G[order(LOO_R2_LAI)]$feature_set)]
Gm<-melt(G,id.vars=c("feature_set","tier"),measure.vars=c("LOO_R2_LAI","cor_dTmax"),variable.name="metric",value.name="v")
Gm[,v:=ifelse(metric=="cor_dTmax",abs(v),v)]
Gm[,metric:=factor(metric,levels=c("LOO_R2_LAI","cor_dTmax"),labels=c("Recovers the LiDAR LAI (R²)","Predicts the observed buffering (|r|)"))]
tpal<-c("operational (no LiDAR)"="#D55E00","oracle (needs LiDAR)"="#0072B2","ceiling"="grey40")
pG<-ggplot(Gm,aes(v,feature_set,fill=tier))+geom_col(width=0.7)+facet_wrap(~metric)+
  scale_fill_manual(values=tpal,guide="none")+scale_x_continuous(limits=c(0,1))+
  labs(x=NULL,y=NULL,title="Operational S2+height recovers ~60–80%; adding LiDAR reaches the ceiling",
       subtitle="But S2+FORMS-H ≈ FORMS-H alone: canopy height, not optical S2, carries the signal")+
  theme_article(11)+theme(legend.position="none")
ggsave_article(file.path(FIG,"FigG_ablation"),pG,8.6,4.2)

## ---- Fig H : validation scatter (amplitude collapse) -----------------------
H<-fread(file.path(TAB,"figH_scatter.csv"))
Hm<-rbind(data.table(obs=H$do,sim=H$LiDAR,src="LiDAR full"),data.table(obs=H$do,sim=H$S2,src="Sentinel-2 (dynamic)"))
lim<-range(c(Hm$obs,Hm$sim)); st<-Hm[,.(r=cor(sim,obs),rmse=sqrt(mean((sim-obs)^2)),sdr=sd(sim)/sd(obs)),by=src]
st[,lab:=sprintf("r = %.2f\nRMSE = %.2f °C\nSD rec. = %.2f",r,rmse,sdr)]
pH<-ggplot(Hm,aes(obs,sim,colour=src))+geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey55")+
  geom_point(size=2.2,alpha=0.85)+facet_wrap(~src)+coord_equal(xlim=lim,ylim=lim)+
  geom_text(data=st,aes(x=lim[1]+0.2,y=lim[2]-0.3,label=lab),hjust=0,vjust=1,size=3,colour="grey20",inherit.aes=FALSE)+
  scale_colour_manual(values=c("LiDAR full"="#0072B2","Sentinel-2 (dynamic)"="#D55E00"),guide="none")+
  labs(x="Observed ΔTmax offset (°C)",y="Simulated ΔTmax offset (°C)",
       title="The Sentinel-2 field collapses to a flat range",
       subtitle="Observations span −2.6 to +5 °C; simulations recover a fraction of that spread")+
  theme_article(12)
ggsave_article(file.path(FIG,"FigH_validation_scatter"),pH,8.0,4.6)

## ---- Fig I : seasonal windows ----------------------------------------------
I<-fread(file.path(TAB,"Table5_shoulders_windows_genuine.csv"))
keep<-c("CONST (LAI fixe an)","STATIC (pheno param.)","DYN S2 (mag+timing)","FUSION par couche")
lab<-c("CONST (LAI fixe an)"="Fixed (constant)","STATIC (pheno param.)"="Parametric phenology","DYN S2 (mag+timing)"="Dynamic S2","FUSION par couche"="Layered fusion")
I<-I[scenario%in%keep]; I[,scen:=factor(lab[scenario],levels=lab)]; I[,window:=factor(window,levels=c("leafout","summer","autumn"),labels=c("Leaf-out (Apr–May)","Summer (JJAS)","Autumn (Oct–Nov)"))]
pI<-ggplot(I,aes(scen,R2,fill=scen))+geom_col(width=0.72)+facet_wrap(~window)+
  scale_fill_manual(values=c("Fixed (constant)"="#666666","Parametric phenology"="#0072B2","Dynamic S2"="#D55E00","Layered fusion"="#CC79A7"),guide="none")+
  labs(x=NULL,y="Between-plot R² (ΔTmax)",title="Time-varying LAI does not beat fixed LAI at the shoulders",
       subtitle="At leaf-out the parametric curve wins; dynamic Sentinel-2 is weakest (v3.2.3)")+
  theme_article(11)+theme(axis.text.x=element_text(angle=30,hjust=1))
ggsave_article(file.path(FIG,"FigI_seasonal_windows"),pI,8.4,4.6)

cat("DONE: FigE/F/G/H/I -> manuscripts/ch3/figures/\n")
