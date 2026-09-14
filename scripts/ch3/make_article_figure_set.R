# Chapter 3 — assemble the article figure set with narrative numbering + manifest.
# Edit MAIN / ANNEX to reorder; never rename source files. Run from NC_Full root.
suppressPackageStartupMessages({ library(here); library(data.table) })
src <- here("outputs","figures_article_ch3"); dst <- file.path(src,"numbered")
dir.create(dst, showWarnings=FALSE, recursive=TRUE)
# narrative number -> produced figure stem (Fig 2 merged into Fig 1; Fig 5 from z_Example)
# Body Fig 0-8 -> submission Figure 1-9 (0-indexed body figure becomes Figure_1, etc.)
MAIN <- list(
  "Figure_1"="Fig0_lai_sensitivity",            # body Fig 0
  "Figure_2"="Fig1_saturation_scatter",         # body Fig 1
  "Figure_3"="Fig2_saturation_distribution",    # body Fig 2
  "Figure_4"="Fig3_magnitude_recovery",         # body Fig 3
  "Figure_5"="Fig4_dopt_transfer",              # body Fig 4
  "Figure_6"="Fig5_premise_slope_LAI",          # body Fig 5
  "Figure_7"="Fig_stratified_diagnostic",       # body Fig 6
  "Figure_8"="Fig_error_vs_structure_abs",      # body Fig 7
  "Figure_9"="Fig_monthly_metrics")             # body Fig 8
ANNEX <- list(
  "Figure_A1"="FigH1_hybrid_vs_lidar", "Figure_A2"="FigH2_s2_tracks_top",
  "Figure_A3"="FigH3_hybrid_seasonal", "Figure_A4"="Fig4e_summer_lai_scenarios",
  "Figure_A5"="FigSh_shapley_clhs_clusters", "Figure_A6"="FigS1_annual_series",
  "Figure_A7"="FigS2_perplot_corrections", "Figure_A8"="Fig_monthly_obs_regime")
man <- list(); missing <- character()
for (set in list(MAIN, ANNEX)) for (num in names(set)) {
  stem <- set[[num]]; ok <- FALSE
  for (ext in c("png","pdf")) {
    f <- file.path(src, paste0(stem,".",ext))
    if (file.exists(f)) { file.copy(f, file.path(dst, paste0(num,".",ext)), overwrite=TRUE); ok <- TRUE }
  }
  if (ok) man[[length(man)+1]] <- data.table(number=num, source=stem)
  else { missing <- c(missing, sprintf("%s <- %s", num, stem)) }
}
fwrite(rbindlist(man), file.path(dst,"manifest.csv"))
writeLines(c("# MISSING figures (not produced in NC_Full):", missing), file.path(dst,"MISSING.txt"))
cat("numbered set ->", dst, "\n"); print(rbindlist(man))
cat("\nMISSING:\n"); cat(missing, sep="\n"); cat("\n")
