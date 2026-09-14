# Chapter 3 — Fig H1/H2/H3 (d_opt decomposition / layered hybrid) ARTICLE style.
suppressPackageStartupMessages({ library(terra); library(data.table); library(ggplot2); library(here) })
source(here("scripts","_article_style.R"))
sites <- c("Aigoual","Blois","Mormal"); lai_min<-2; k_scale<-0.5/0.65; d<-7
als<-function(s) here("03_RESULTS",s,"Metrics","Deciduous_Only","ladstack_classic.tif")
hmx<-function(s) here("03_RESULTS",s,"Metrics","Deciduous_Only","max_res_10_m.tif")
s2a<-function(s) here("output","intermediate","sm6",s,"s2lai_summer_atbd_T_res_10_m.tif")
dpc<-function(s) here("output","intermediate","lai_als_dopt",s,"LAI_ALS_dopt_common.tif")
outdir<-here("outputs","figures_article_ch3"); sv<-function(g,nm,w,h){ ggsave_article(file.path(outdir,nm),g,w,h)
  ggsave_article(here("manuscripts/ch3","figures",paste0(nm,"_styled")),g,w,h) }
pool<-rbindlist(lapply(sites,function(s){ g<-sum(rast(als(s)),na.rm=TRUE)
  ld<-function(r){if(!compareGeom(g,r,stopOnError=FALSE))r<-resample(r,g);as.numeric(values(r))}
  data.table(site=s, full=as.numeric(values(g))*k_scale, S2=ld(rast(s2a(s))),
             dopt=ld(rast(dpc(s))), Hmax=ld(rast(hmx(s))))[
    is.finite(full)&is.finite(S2)&is.finite(dopt)&is.finite(Hmax)&(full/k_scale)>=lai_min] }))
set.seed(1); P<-pool[, .SD[sample(.N,min(.N,6000))], by=site]
P[, hybrid:=NA_real_]
for(test in sites){ tr<-P[site!=test]; tall<-tr$Hmax>d
  b_top<-mean(tr$dopt[tall])/mean(tr$S2[tall]); b_full<-mean(tr$full[!tall])/mean(tr$S2[!tall])
  P[site==test, hybrid:=ifelse(Hmax>d, S2*b_top+pmax(full-dopt,0), S2*b_full)] }
r2d<-P[, .(r2=cor(hybrid,full)^2), by=site]; r2d[, lbl:=sprintf("italic(R)^2==%.2f",r2)]
saveRDS(P, file.path(outdir,"figH_hybrid_pool.rds"))
# H1
g1<-ggplot(P, aes(full,hybrid,colour=site)) + geom_point(alpha=.06,size=.5) +
  geom_abline(slope=1,linetype="dashed",colour="grey40") + facet_wrap(~site) +
  geom_text(data=r2d, aes(x=0.4,y=6.6,label=lbl), parse=TRUE, inherit.aes=FALSE, hjust=0, size=3) +
  scale_colour_manual(values=PAL_SITE, guide="none") + coord_equal(xlim=c(0,7),ylim=c(0,7),expand=FALSE) +
  labs(x=expression(full~LiDAR~LAI~(m^2~m^-2)), y=expression(hybrid~LAI~(m^2~m^-2))) + theme_article(11)
sv(g1,"FigH1_hybrid_vs_lidar",9,3.4)
# H2
B<-P[site=="Blois"]
dd<-rbind(data.table(target="full LiDAR LAI", x=B$full,y=B$S2,r2=cor(B$full,B$S2)^2),
          data.table(target="top-d_opt LiDAR LAI", x=B$dopt,y=B$S2,r2=cor(B$dopt,B$S2)^2))
dd[, lab:=sprintf("%s  (R² = %.2f)",target,r2)]
g2<-ggplot(dd, aes(x,y)) + geom_point(alpha=.08,size=.5,colour="#1A9850") +
  geom_smooth(method="lm",se=FALSE,colour="black",linewidth=.8) + facet_wrap(~lab, scales="free_x") +
  labs(x=expression(LiDAR~LAI~(m^2~m^-2)), y=expression(LAI[S2]~(m^2~m^-2))) + theme_article(11)
sv(g2,"FigH2_s2_tracks_top",8,3.4)
# H3
ts<-as.data.table(readRDS(here("output","intermediate","blois_s2_lai_ts_2021.rds")))
sm<-ts[, .(S2=mean(LAI_S2_ATBD,na.rm=TRUE)), by=doy][order(doy)]
bott<-mean(B$full)-mean(B$dopt); btop<-mean(B$dopt)/mean(B$S2); afull<-mean(B$full)/mean(B$S2)
q05<-quantile(sm$S2,.05); smr<-max(sm$S2); sm[, frac:=pmax((S2-q05)/(smr-q05),0)]
sm[, `:=`(raw_S2=S2, rescaled=S2*afull, hybrid_static=bott+S2*btop, hybrid_pheno=bott*frac+S2*btop)]
ml<-melt(sm, id.vars="doy", measure.vars=c("raw_S2","rescaled","hybrid_static","hybrid_pheno"),
         variable.name="series", value.name="LAI")
labs_s<-c(raw_S2="raw S2 (saturated)", rescaled="S2 x global ratio", hybrid_static="hybrid, static LiDAR bottom", hybrid_pheno="hybrid, phenology-scaled bottom")
ml[, series:=factor(labs_s[as.character(series)], levels=labs_s)]
g3<-ggplot(ml, aes(doy,LAI,colour=series)) +
  annotate("rect",xmin=152,xmax=273,ymin=-Inf,ymax=Inf,fill="grey88",alpha=.6) +
  annotate("text",x=212,y=0.3,label="summer (LiDAR date)",size=2.7,colour="grey35") +
  geom_line(linewidth=0.9) + scale_colour_manual(values=c("#9ecae1","#3182bd","#e6550d","#1A9850"), name=NULL) +
  labs(x="Day of year", y=expression(LAI~(m^2~m^-2))) + theme_article(11) + theme(legend.position="bottom")
sv(g3,"FigH3_hybrid_seasonal",8.5,4.6)
cat("DONE FigH1/H2/H3\n")
