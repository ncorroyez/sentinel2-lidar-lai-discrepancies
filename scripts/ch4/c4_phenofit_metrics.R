# ==============================================================================
# c4_phenofit_metrics.R — Chapter 4: phenometrics with the STANDARD method
# (double-logistic Elmore fit + TRS50 threshold via the `phenofit` package),
# i.e. the exact metric of Cotrina-Sanchez et al. (2026) — makes our S2-vs-GEDI
# EOS gap directly comparable to their continental ~35 days.
#
# Series used (all durable CSV/RDS produced by earlier c4 scripts):
#  - Blois S2 ATBD annual series (19 dates, site mean)          -> SOS/EOS
#  - Blois GEDI per-acquisition adjusted composites (Table14)   -> SOS/EOS
#  - Mormal GEDI monthly adjusted composites (Table12)          -> SOS/EOS
#  - Mormal GEDI per-layer monthly composites (Table17)         -> SOS/EOS
# Output: manuscripts/ch4/tables/Table19_c4_phenofit_trs50.csv
# ==============================================================================
suppressMessages({
  library(phenofit); library(data.table); library(here)
})
dir_tabs <- here::here("manuscripts/ch4", "tables")

# TRS50 on the weekly-interpolated, wWHIT-smoothed series (threshold crossing,
# phenofit::whit2 = the smoother of the Cotrina-Sanchez chain). The Elmore
# double-logistic is kept as a control column: with ONE pooled season of sparse
# points it symmetrises the long autumn shoulder and lands too early, so the
# smoothed-crossing is the primary metric here (Ferrara-style TRS50).
trs50 <- function(doy, val, label) {
  ok <- is.finite(doy) & is.finite(val)
  wk <- seq(min(doy[ok]), max(doy[ok]), by = 7)
  vw <- approx(doy[ok], val[ok], wk)$y
  sm <- as.numeric(whit2(vw, lambda = 10))
  half <- min(sm) + (max(sm) - min(sm)) / 2
  up <- which(sm >= half)
  sos <- if (length(up)) wk[up[1]] else NA_real_
  eos <- if (length(up)) wk[rev(up)[1]] else NA_real_
  elm <- tryCatch({
    fit <- curvefit(vw, wk, tout = 1:365, methods = "Elmore")
    PhenoTrs(fit$model$Elmore, t = 1:365, approach = "Trs", trs = 0.5,
             IsPlot = FALSE)
  }, error = function(e) c(sos = NA_real_, eos = NA_real_))
  data.table(series = label, n_pts = sum(ok),
             SOS_trs50 = round(sos), EOS_trs50 = round(eos),
             SOS_elmore = round(elm[["sos"]]), EOS_elmore = round(elm[["eos"]]))
}

# Blois S2 (site-mean ATBD, 19 dates 2021)
s2 <- as.data.table(readRDS(here::here("output", "intermediate",
                                       "blois_s2_lai_ts_2021.rds")))
s2 <- s2[, .(lai = mean(LAI_S2_ATBD, na.rm = TRUE)), by = doy][order(doy)]

# Blois GEDI per-acquisition, structure-adjusted (power pool of Table14)
t14 <- fread(file.path(dir_tabs, "Table14_c4_blois_dates.csv"))
t14 <- t14[pool == "power"][, doy := as.integer(format(as.Date(date), "%j"))]
# drop only May 29 2021 (n = 5 and inconsistent); Apr 11 2022 (n = 5) is kept:
# it is the sole pre-leaf-out baseline that anchors the seasonal amplitude
t14 <- t14[as.Date(date) != as.Date("2021-05-29")]

# Mormal monthly adjusted composites (Table12) + per-layer (Table17)
t12 <- fread(file.path(dir_tabs, "Table12_c4_gedi_monthly.csv"))[site == "Mormal"]
t17 <- fread(file.path(dir_tabs, "Table17_c4_layer_monthly.csv"))[site == "Mormal"]

res <- rbind(
  trs50(s2$doy, s2$lai,                       "Blois S2 ATBD series"),
  trs50(t14$doy, t14$pai_adj,                 "Blois GEDI (adjusted, 9 dates)"),
  trs50(30.4 * (t12$m - 0.5), t12$pai_adj,    "Mormal GEDI total (monthly)"),
  t17[, trs50(30.4 * (m - 0.5), val, paste("Mormal GEDI", layer[1])), by = layer][, -1])
print(res)
cat("\nEOS gap Blois (GEDI - S2):",
    res[series == "Blois GEDI (adjusted, 9 dates)", EOS_trs50] -
      res[series == "Blois S2 ATBD series", EOS_trs50], "days",
    "(Cotrina-Sanchez 2026 continental mean: ~35 days)\n")
fwrite(res, file.path(dir_tabs, "Table19_c4_phenofit_trs50.csv"))
cat("DONE\n")
