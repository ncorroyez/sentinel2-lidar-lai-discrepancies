# ==============================================================================
# Chapter 3 — illustration figures for the d_opt decomposition / layered hybrid,
# and the SEASONAL behaviour of the fusion (Blois annual S2 series).
#   FigH1  hybrid LAI vs full LiDAR LAI, summer, per site (R²=0.70)
#   FigH2  diagnostic: S2 tracks the TOP (d_opt) not the FULL LAI (Blois)
#   FigH3  seasonal hybrid trajectory: static vs phenology-scaled bottom — the
#          LiDAR bottom is a summer snapshot, so a static bottom over-estimates winter.
# Run from NC_Full root:  Rscript scripts/ch3/c3_hybrid_figs.R
# ==============================================================================
suppressPackageStartupMessages({ library(terra); library(data.table); library(ggplot2); library(patchwork); library(here) })
sites <- c("Aigoual","Blois","Mormal"); lai_min <- 2; k_scale <- 0.5/0.65; d <- 7
als <- function(s) here("03_RESULTS", s, "Metrics", "Deciduous_Only", "ladstack_classic.tif")
hmx <- function(s) here("03_RESULTS", s, "Metrics", "Deciduous_Only", "max_res_10_m.tif")
s2a <- function(s) here("output","intermediate","sm6", s, "s2lai_summer_atbd_T_res_10_m.tif")
dpc <- function(s) here("output","intermediate","lai_als_dopt", s, "LAI_ALS_dopt_common.tif")
site_cols <- c(Aigoual="#E8746A", Blois="#2EA02E", Mormal="#5B8FF0")

pool <- rbindlist(lapply(sites, function(s){ g<-sum(rast(als(s)),na.rm=TRUE)
  ld<-function(r){if(!compareGeom(g,r,stopOnError=FALSE))r<-resample(r,g);as.numeric(values(r))}
  data.table(site=s, full=as.numeric(values(g))*k_scale, S2=ld(rast(s2a(s))),
             dopt=ld(rast(dpc(s))), Hmax=ld(rast(hmx(s))))[
    is.finite(full)&is.finite(S2)&is.finite(dopt)&is.finite(Hmax)&(full/k_scale)>=lai_min] }))
set.seed(1); P <- pool[, .SD[sample(.N, min(.N, 6000))], by=site]

# hybrid (LOSO-calibrated top), summer
P[, hybrid := NA_real_]
for (test in sites) { tr<-P[site!=test]; tall<-tr$Hmax>d
  b_top<-mean(tr$dopt[tall])/mean(tr$S2[tall]); b_full<-mean(tr$full[!tall])/mean(tr$S2[!tall])
  P[site==test, hybrid := ifelse(Hmax>d, S2*b_top + pmax(full-dopt,0), S2*b_full)] }
r2lab <- P[, sprintf("%s: R²=%.2f", site, cor(hybrid,full)^2), by=site]$V1

# ---- FigH1: hybrid vs full LiDAR, per site ----
g1 <- ggplot(P, aes(full, hybrid, colour=site)) +
  geom_point(alpha=.06, size=.5) + geom_abline(slope=1, linetype="dashed") +
  facet_wrap(~site) + scale_colour_manual(values=site_cols, guide="none") +
  coord_equal(xlim=c(0,7), ylim=c(0,7)) +
  labs(title="Layered hybrid (S2 top + LiDAR bottom) reconstructs total LAI",
       subtitle=paste("Summer, per-pixel vs full LiDAR LAI.", paste(r2lab, collapse="  ·  ")),
       x=expression(full~LiDAR~LAI~(m^2/m^2)), y="hybrid LAI") +
  theme_minimal(base_size=11)
ggsave(here("manuscripts/ch3","figures","FigH1_hybrid_vs_lidar.png"), g1, width=9, height=3.6, dpi=150)

# ---- FigH2: S2 tracks the TOP, not the FULL (Blois) ----
B <- P[site=="Blois"]
dd <- rbind(data.table(target="full LiDAR LAI", x=B$full, y=B$S2, r2=cor(B$full,B$S2)^2),
            data.table(target="top-d_opt LiDAR LAI", x=B$dopt, y=B$S2, r2=cor(B$dopt,B$S2)^2))
dd[, lab := sprintf("%s (R²=%.2f)", target, r2)]
g2 <- ggplot(dd, aes(x, y)) + geom_point(alpha=.08, size=.5, colour="#2EA02E") +
  geom_smooth(method="lm", se=FALSE, colour="black", linewidth=.8) + facet_wrap(~lab, scales="free_x") +
  labs(title="Blois: Sentinel-2 tracks the canopy TOP it senses, not the full LAI",
       x=expression(LiDAR~LAI~(m^2/m^2)), y=expression(LAI[S2])) + theme_minimal(base_size=11)
ggsave(here("manuscripts/ch3","figures","FigH2_s2_tracks_top.png"), g2, width=8, height=3.6, dpi=150)

# ---- FigH3: seasonal hybrid (Blois annual S2 series) ----
ts <- as.data.table(readRDS(here("output","intermediate","blois_s2_lai_ts_2021.rds")))
sm <- ts[, .(S2=mean(LAI_S2_ATBD, na.rm=TRUE)), by=doy][order(doy)]
bott <- mean(B$full) - mean(B$dopt)                       # LiDAR bottom (summer)
btop <- mean(B$dopt) / mean(B$S2)                         # top calibration
afull<- mean(B$full) / mean(B$S2)
q05 <- quantile(sm$S2, .05); smr <- max(sm$S2)
sm[, frac := pmax((S2 - q05)/(smr - q05), 0)]             # greenness fraction (0 winter,1 summer)
sm[, `:=`(raw_S2 = S2,
          rescaled = S2 * afull,                          # global ratio to full magnitude
          hybrid_static = bott + S2*btop,                 # bottom fixed (summer snapshot)
          hybrid_pheno  = bott*frac + S2*btop)]           # bottom follows phenology
ml <- melt(sm, id.vars="doy", measure.vars=c("raw_S2","rescaled","hybrid_static","hybrid_pheno"),
           variable.name="series", value.name="LAI")
labs_s <- c(raw_S2="raw S2 (saturated)", rescaled="S2 × global ratio",
            hybrid_static="hybrid, static LiDAR bottom", hybrid_pheno="hybrid, phenology-scaled bottom")
ml[, series := factor(labs_s[as.character(series)], levels=labs_s)]
g3 <- ggplot(ml, aes(doy, LAI, colour=series)) + geom_line(linewidth=1) +
  annotate("rect", xmin=152, xmax=273, ymin=-Inf, ymax=Inf, fill="grey85", alpha=.5) +
  annotate("text", x=212, y=0.3, label="summer (LiDAR date)", size=2.8, colour="grey35") +
  scale_colour_manual(values=c("#9ecae1","#3182bd","#e6550d","#31a354"), name=NULL) +
  labs(title="Seasonal behaviour of the layered fusion (Blois 2021)",
       subtitle="The LiDAR bottom is a summer snapshot: a STATIC bottom over-estimates winter LAI;\nscaling the bottom by phenology is needed off-summer. Valid by construction only near the LiDAR date.",
       x="Day of year", y=expression(LAI~(m^2/m^2))) +
  theme_minimal(base_size=11) + theme(legend.position="bottom", plot.subtitle=element_text(size=8.5))
ggsave(here("manuscripts/ch3","figures","FigH3_hybrid_seasonal.png"), g3, width=8.5, height=5, dpi=150)

cat("means: bottom=",round(bott,2)," b_top=",round(btop,2)," a_full=",round(afull,2),"\n")
cat("DONE FigH1, FigH2, FigH3\n")
