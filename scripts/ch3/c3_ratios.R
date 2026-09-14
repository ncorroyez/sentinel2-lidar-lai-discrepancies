# ==============================================================================
# Chapter 3 — DURABLE reconstruction of the correction-ratio table (Table_ratios),
# at the study convention k_select = 0.65. Replaces the lost /tmp script.
# Ratios of per-site MEAN LAI over the canonical SM6a pool (deciduous, fCover>=90%,
# raw-k0.5 LAI >= 2; full LiDAR LAI then reported at k=0.65 via k_scale=0.5/0.65).
#   A = LAI_ALS        / LAI_S2_ATBD     (full LiDAR / ATBD S2)  <- the transferable one
#   B = LAI_ALS_dopt   / LAI_S2_ATBD     (d_opt-truncated LiDAR / ATBD S2), common 7m & per-site
#   C = LAI_ALS_dopt   / LAI_S2_opt      (top / forest-tuned S2)
#   D = LAI_ALS        / LAI_S2_opt      (full / forest-tuned S2)
# NB LAI_ALS_dopt rasters are ALREADY at k=0.65 (step 07 scale_factor). Only the
# full-canopy sum(ladstack) needs k_scale here.
# Run from NC_Full root:  Rscript scripts/ch3/c3_ratios.R
# ==============================================================================
suppressPackageStartupMessages({ library(terra); library(data.table); library(here) })
sites   <- c("Aigoual", "Blois", "Mormal")
lai_min <- 2.0
k_scale <- 0.5 / 0.65                                  # PAD k=0.5 -> study k=0.65
als  <- function(s) here("03_RESULTS", s, "Metrics", "Deciduous_Only", "ladstack_classic.tif")
s2a  <- function(s) here("output","intermediate","sm6", s, "s2lai_summer_atbd_T_res_10_m.tif")
s2o  <- function(s) here("output","intermediate","sm6", s, "s2lai_summer_opt_res_10_m.tif")
dpc  <- function(s) here("output","intermediate","lai_als_dopt", s, "LAI_ALS_dopt_common.tif")
dps  <- function(s) here("output","intermediate","lai_als_dopt", s, "LAI_ALS_dopt_per_site.tif")

M <- rbindlist(lapply(sites, function(s) {
  laR <- sum(rast(als(s)), na.rm = TRUE)              # raw k=0.5
  sa  <- rast(s2a(s)); so <- rast(s2o(s)); dc <- rast(dpc(s)); dp <- rast(dps(s))
  for (nm in c("sa","so","dc","dp"))
    if (!compareGeom(laR, get(nm), stopOnError=FALSE)) assign(nm, resample(get(nm), laR))
  d <- data.table(laR = as.numeric(values(laR)), sa = as.numeric(values(sa)),
                  so = as.numeric(values(so)), dc = as.numeric(values(dc)),
                  dp = as.numeric(values(dp)))
  d <- d[is.finite(laR) & is.finite(sa) & laR >= lai_min]        # canonical pool
  data.table(site = s,
             LAI_ALS       = mean(d$laR) * k_scale,             # k=0.65
             LAI_ALS_dopt_common  = mean(d$dc, na.rm=TRUE),     # already k=0.65
             LAI_ALS_dopt_persite = mean(d$dp, na.rm=TRUE),
             LAI_S2_ATBD   = mean(d$sa),
             LAI_S2_opt    = mean(d$so, na.rm=TRUE))
}))
M[, `:=`(
  A_full_ATBD        = LAI_ALS / LAI_S2_ATBD,
  B_common_ATBD      = LAI_ALS_dopt_common  / LAI_S2_ATBD,
  B_persite_ATBD     = LAI_ALS_dopt_persite / LAI_S2_ATBD,
  C_common_opt       = LAI_ALS_dopt_common  / LAI_S2_opt,
  D_full_opt         = LAI_ALS / LAI_S2_opt)]
cv <- function(x) sd(x) / mean(x)
ratios <- c("A_full_ATBD","B_common_ATBD","B_persite_ATBD","C_common_opt","D_full_opt")
cvrow  <- as.list(setNames(sapply(ratios, function(r) cv(M[[r]])), ratios))
out <- rbind(M, data.table(site = "CV_crosssite",
             LAI_ALS=NA, LAI_ALS_dopt_common=NA, LAI_ALS_dopt_persite=NA,
             LAI_S2_ATBD=NA, LAI_S2_opt=NA, as.data.table(cvrow)), fill = TRUE)

cat("=== Per-site mean LAI (k=0.65) ===\n")
print(M[, .(site, LAI_ALS, LAI_ALS_dopt_common, LAI_S2_ATBD, LAI_S2_opt)], digits = 3)
cat("\n=== Correction ratios (k=0.65) + cross-site CV ===\n")
print(out[, c("site", ratios), with = FALSE], digits = 3)

fwrite(out, here("output","tables","Table_ratios_correction.csv"))
fwrite(out, here("manuscripts/ch3","tables","Table_ratios_correction.csv"))
cat("\nDONE\n")
