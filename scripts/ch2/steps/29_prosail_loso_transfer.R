# ---
# title:  29_prosail_loso_transfer.R
# desc:   Leave-one-site-out transfer test of the PROSAIL LUT configurations,
#         and dispersion of RMSE over the full 225-configuration ensemble.
#
#         Answers the R4 comment on revision 3 (RSE-D-25-04417): the Pareto
#         configuration is currently selected per site, with nothing held out,
#         and scoring against LAI_ALS_dopt lowers RMSE on its own. This step
#         re-scores the existing 225-configuration ensemble; no PROSAIL
#         training and no inversion are re-run.
#
#         Three analyses:
#           (A) validation — reproduces the four Table A.16 rows per site.
#           (B) LOSO on the configuration only. The LUT is chosen on the two
#               training sites among the 135 configurations whose LAI prior is
#               not derived from LiDAR (lai in {ATBD, OPT#1, OPT#2}), then
#               applied unchanged at the held-out site, at that site's own
#               Pareto d_opt.
#           (C) LOSO on the configuration AND the depth. Both the integration
#               depth and the LUT come from the two training sites, so nothing
#               about the held-out site enters the selection.
#           (D) dispersion of RMSE over the 225 configurations at each site.
#
#         Writes:
#           output/tables/Table_A17_loso_transfer.csv
#           output/tables/Table_A18_rmse_spread_225.csv
#           output/tables/Table_A16_partial_attribution_check.csv
#
# Run from the project root (NC_Full/):
#   Rscript scripts/ch2/steps/29_prosail_loso_transfer.R
# ---

suppressMessages({
  library(here)
  library(data.table)
  library(cli)
})

source(here::here("R", "paths.R"))
source(here::here("R", "sm5_dopt.R"))

# ── Configuration — must match scripts/ch2/steps/12_sm5_select_prosail_opt.R ──────

norm_select         <- "DSM_keepTrees"
h_min_select        <- 10L
lai_scenario_select <- "common"
max_depth           <- 10L      # d_thr; matches scripts/ch2/steps/06_sm5_select_dopt.R
sites               <- c("Aigoual", "Blois", "Mormal")
criteria            <- c("R", "RMSE", "Bias", "Slope")

# LAI prior levels that do not sample the LiDAR LAI raster (R/prosail_lut.R):
#   1 = ATBD, 2 = OPT#1, 3 = OPT#2, 4 = LAI_ALS, 5 = LAI_ALS_dopt
lai_free_levels <- 1:3

sm5_dir <- file.path(paths$output, "intermediate", "sm5")
out_dir <- file.path(paths$output, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load ──────────────────────────────────────────────────────────────────────

metrics <- data.table::fread(
  file.path(sm5_dir, "all_results_combined_LIDFa_lai_LMA_BROWN.csv"))
metrics <- metrics[Norm == norm_select & h_min == h_min_select &
                     lai_scenario == lai_scenario_select & Site %in% sites]

atbd <- data.table::fread(
  file.path(sm5_dir, "all_results_atbd_LIDFa_lai_LMA_BROWN.csv"))
atbd <- atbd[Norm == norm_select & h_min == h_min_select & Site %in% sites]

dopt_ref <- data.table::fread(file.path(sm5_dir, "dopt_reference.csv"))
dopt_ref <- dopt_ref[Norm == norm_select & method_dopt == "pareto" & Site %in% sites]
d_opt_site <- setNames(as.integer(dopt_ref$d_opt), dopt_ref$Site)

pro_opt <- data.table::fread(file.path(sm5_dir, "prosail_opt.csv"))
pro_opt <- pro_opt[d_opt_source == "per_site" & Norm == norm_select]
col_pareto <- setNames(pro_opt$Column_opt, pro_opt$Site)

metrics[, lai_lvl := as.integer(sub(".*_lai=([0-9]+)_.*", "\\1", Column))]
ATBD_COL <- "LIDFa=1_lai=1_LMA=1_BROWN=1"

# ── Helpers ───────────────────────────────────────────────────────────────────

# Minimisable form of the four criteria, rounded to 2 dp exactly as
# compute_pareto_front() does, so scores are comparable with the paper's.
minimisable <- function(dt) {
  cbind(R     = 1 - round(dt$R, 2),
        RMSE  = round(dt$RMSE, 2),
        Bias  = abs(round(dt$Bias, 2)),
        Slope = abs(round(dt$Slope, 2) - 1))[, criteria, drop = FALSE]
}

# Avg_score in the sense of compute_pareto_front(norm_method = "minmax"):
# higher is better. Normalisation is over the whole candidate set here
# (norm_scope = "all"), because a Front-1 scope is undefined once scores from
# two sites have to be averaged.
avg_score <- function(dt) {
  m <- minimisable(dt)
  m_norm <- apply(m, 2L, function(x) {
    lo <- min(x, na.rm = TRUE); hi <- max(x, na.rm = TRUE)
    if (hi == lo) rep(0, length(x)) else (x - lo) / (hi - lo)
  })
  round(1 - rowMeans(m_norm), 4)
}

pct_retained <- function(v_cfg, v_atbd, v_ref) {
  round(100 * (v_atbd - v_cfg) / (v_atbd - v_ref), 0)
}

row_at <- function(site, depth, column) {
  metrics[Site == site & Depth == depth & Column == column]
}

# ── (A) Validation — reproduce Table A.16 ─────────────────────────────────────

cli::cli_h1("(A) Reproduction of Table A.16")

a16 <- rbindlist(lapply(sites, function(s) {
  d <- d_opt_site[[s]]
  cp <- col_pareto[[s]]
  lv <- as.integer(sub(".*_lai=([0-9]+)_.*", "\\1", cp))
  # (iii) same non-LAI priors, LAI prior reverted to ATBD
  cp_rev <- sub("_lai=[0-9]+_", "_lai=1_", cp)
  # (iv) best fully ALS-free configuration by RMSE
  free <- metrics[Site == s & Depth == d & lai_lvl %in% lai_free_levels]
  cp_free <- free[which.min(RMSE), Column]

  cfgs <- list(
    list(lab = "ATBD baseline",             col = ATBD_COL),
    list(lab = "Full Pareto site-opt",      col = cp),
    list(lab = "Same non-LAI + lai = ATBD", col = cp_rev),
    list(lab = "Best ALS-free",             col = cp_free))

  ref <- row_at(s, d, cp)
  bas <- row_at(s, d, ATBD_COL)

  rbindlist(lapply(cfgs, function(cfg) {
    r <- row_at(s, d, cfg$col)
    data.table(
      Site = s, d_opt = d, Configuration = cfg$lab, Column = cfg$col,
      r = round(r$R, 2), RMSE = round(r$RMSE, 2), Bias = round(r$Bias, 2),
      Slope = round(r$Slope, 2),
      pct_RMSE  = pct_retained(round(r$RMSE, 2), round(bas$RMSE, 2), round(ref$RMSE, 2)),
      pct_Bias  = pct_retained(abs(round(r$Bias, 2)), abs(round(bas$Bias, 2)),
                               abs(round(ref$Bias, 2))),
      lai_prior_from_LiDAR = as.integer(sub(".*_lai=([0-9]+)_.*", "\\1", cfg$col)) %in% 4:5)
  }))
}))
print(a16)
fwrite(a16, file.path(out_dir, "Table_A16_partial_attribution_check.csv"))

# ── (B) and (C) Leave-one-site-out ────────────────────────────────────────────

cli::cli_h1("(B)/(C) Leave-one-site-out transfer")

loso <- rbindlist(lapply(sites, function(s_out) {
  s_tr <- setdiff(sites, s_out)
  d_out <- d_opt_site[[s_out]]

  # Reference improvement reported in the paper for the held-out site
  bas <- row_at(s_out, d_out, ATBD_COL)
  ref <- row_at(s_out, d_out, col_pareto[[s_out]])
  rmse_b <- round(bas$RMSE, 2); rmse_r <- round(ref$RMSE, 2)
  bias_b <- abs(round(bas$Bias, 2)); bias_r <- abs(round(ref$Bias, 2))

  emit <- function(mode, column, depth, note) {
    r <- row_at(s_out, depth, column)
    data.table(
      Held_out = s_out, Trained_on = paste(s_tr, collapse = " + "),
      Mode = mode, Depth_applied = depth, Column = column,
      r = round(r$R, 2), RMSE = round(r$RMSE, 2), Bias = round(r$Bias, 2),
      Slope = round(r$Slope, 2),
      pct_RMSE = pct_retained(round(r$RMSE, 2), rmse_b, rmse_r),
      pct_Bias = pct_retained(abs(round(r$Bias, 2)), bias_b, bias_r),
      Note = note)
  }

  # ---- (B) configuration held out, each site keeps its own Pareto depth -----
  sc_tr <- lapply(s_tr, function(s) {
    x <- metrics[Site == s & Depth == d_opt_site[[s]]]
    setorder(x, Column)
    data.table(Column = x$Column, score = avg_score(x))
  })
  tr <- Reduce(function(a, b) merge(a, b, by = "Column"), sc_tr)
  tr[, score_mean := rowMeans(.SD), .SDcols = patterns("^score")]
  tr[, lai_lvl := as.integer(sub(".*_lai=([0-9]+)_.*", "\\1", Column))]

  cfg_free <- tr[lai_lvl %in% lai_free_levels][which.max(score_mean), Column]
  cfg_any  <- tr[which.max(score_mean), Column]

  out <- list(
    emit("B — LUT from training sites, ALS-free priors", cfg_free, d_out,
         "config selected on the two other sites; depth = held-out site's d_opt"),
    emit("B — LUT from training sites, all 225 priors", cfg_any, d_out,
         "same, LiDAR-derived LAI priors allowed as candidates"))

  # ---- (C) depth AND configuration held out --------------------------------
  # Two rules for transferring the depth. (C1) is the primary one: the mean of
  # the two training sites' published d_opt, each obtained with the paper's own
  # Pareto rule (scripts/ch2/steps/06_sm5_select_dopt.R). (C2) re-runs that Pareto
  # rule jointly on the two training sites, with the min-max normalisation taken
  # over all candidate depths, since a Front-1 scope is undefined once two sites
  # have to be averaged.
  d_c1 <- as.integer(round(mean(unlist(d_opt_site[s_tr]))))

  sc_d <- lapply(s_tr, function(s) {
    x <- atbd[Site == s & Depth <= max_depth & Depth >= 1L]
    setorder(x, Depth)
    data.table(Depth = x$Depth, score = avg_score(x))
  })
  trd <- Reduce(function(a, b) merge(a, b, by = "Depth"), sc_d)
  trd[, score_mean := rowMeans(.SD), .SDcols = patterns("^score")]
  d_c2 <- trd[which.max(score_mean), Depth]

  cfg_at_depth <- function(dd) {
    sc <- lapply(s_tr, function(s) {
      x <- metrics[Site == s & Depth == dd]
      setorder(x, Column)
      data.table(Column = x$Column, score = avg_score(x))
    })
    t2 <- Reduce(function(a, b) merge(a, b, by = "Column"), sc)
    t2[, score_mean := rowMeans(.SD), .SDcols = patterns("^score")]
    t2[, lai_lvl := as.integer(sub(".*_lai=([0-9]+)_.*", "\\1", Column))]
    t2[lai_lvl %in% lai_free_levels][which.max(score_mean), Column]
  }

  out <- c(out, list(
    emit("C — depth and LUT from training sites", cfg_at_depth(d_c1), d_c1,
         sprintf("d = %d m (mean of the two training sites' d_opt) and LUT both selected on %s",
                 d_c1, paste(s_tr, collapse = " + "))),
    emit("C — ATBD LUT at the transferred depth", ATBD_COL, d_c1,
         sprintf("baseline at the transferred depth d = %d m", d_c1)),
    emit("C2 — depth from a joint Pareto over the training sites",
         cfg_at_depth(d_c2), d_c2,
         sprintf("robustness check: d = %d m from the joint depth Pareto", d_c2))))

  rbindlist(out)
}))
print(loso)
fwrite(loso, file.path(out_dir, "Table_A17_loso_transfer.csv"))

# ── (D) Dispersion of RMSE over the 225 configurations ────────────────────────

cli::cli_h1("(D) RMSE dispersion over the 225-configuration ensemble")

spread <- rbindlist(lapply(sites, function(s) {
  d <- d_opt_site[[s]]
  x <- metrics[Site == s & Depth == d]
  x225 <- round(x$RMSE, 2)
  x135 <- round(x[lai_lvl %in% lai_free_levels]$RMSE, 2)
  rmse_atbd   <- round(row_at(s, d, ATBD_COL)$RMSE, 2)
  rmse_pareto <- round(row_at(s, d, col_pareto[[s]])$RMSE, 2)
  data.table(
    Site = s, d_opt = d,
    n_all = length(x225), n_ALSfree = length(x135),
    RMSE_min = min(x225), RMSE_q1 = round(quantile(x225, .25), 2),
    RMSE_median = round(median(x225), 2), RMSE_q3 = round(quantile(x225, .75), 2),
    RMSE_max = max(x225),
    RMSE_min_ALSfree = min(x135), RMSE_median_ALSfree = round(median(x135), 2),
    RMSE_ATBD = rmse_atbd, ATBD_percentile = round(100 * mean(x225 <= rmse_atbd), 0),
    RMSE_Pareto = rmse_pareto,
    Pareto_percentile = round(100 * mean(x225 <= rmse_pareto), 0),
    n_better_than_Pareto = sum(x225 < rmse_pareto),
    n_ALSfree_better_than_ATBD = sum(x135 < rmse_atbd))
}))
print(spread)
fwrite(spread, file.path(out_dir, "Table_A18_rmse_spread_225.csv"))

cli::cli_alert_success("Written to {out_dir}")

# ── (E) Robustness of the training criterion, and cross-site aggregates ──────

cli::cli_h1("(E) Alternative training criterion and aggregates")

# Alternative to the paper's Pareto score: mean over the two training sites of
# the RMSE reduction relative to that site's own ATBD baseline. Expressed as a
# ratio so that a site with a large baseline error does not dominate.
alt <- rbindlist(lapply(sites, function(s_out) {
  s_tr <- setdiff(sites, s_out)
  red <- lapply(s_tr, function(s) {
    d <- d_opt_site[[s]]
    x <- metrics[Site == s & Depth == d]
    setorder(x, Column)
    base <- round(row_at(s, d, ATBD_COL)$RMSE, 2)
    data.table(Column = x$Column, red = (base - round(x$RMSE, 2)) / base)
  })
  tr <- Reduce(function(a, b) merge(a, b, by = "Column"), red)
  tr[, red_mean := rowMeans(.SD), .SDcols = patterns("^red")]
  tr[, lai_lvl := as.integer(sub(".*_lai=([0-9]+)_.*", "\\1", Column))]
  cfg <- tr[lai_lvl %in% lai_free_levels][which.max(red_mean), Column]
  d <- d_opt_site[[s_out]]
  r <- row_at(s_out, d, cfg)
  data.table(Held_out = s_out, Criterion = "mean relative RMSE reduction",
             Column = cfg, RMSE = round(r$RMSE, 2), Bias = round(r$Bias, 2))
}))
print(alt)

agg <- function(cols, depths) {
  vals <- mapply(function(s, cc, dd) {
    r <- row_at(s, dd, cc); c(RMSE = round(r$RMSE, 2), Bias = abs(round(r$Bias, 2)))
  }, sites, cols, depths)
  rowMeans(vals)
}
d_own    <- unname(d_opt_site[sites])
m_atbd   <- agg(rep(ATBD_COL, 3L), d_own)
m_pareto <- agg(unname(col_pareto[sites]), d_own)
m_losoB  <- agg(loso[grepl("^B .*ALS-free", Mode)][match(sites, Held_out), Column], d_own)
m_losoC  <- agg(loso[grepl("^C — depth", Mode)][match(sites, Held_out), Column],
                loso[grepl("^C — depth", Mode)][match(sites, Held_out), Depth_applied])
m_losoC2 <- agg(loso[grepl("^C2", Mode)][match(sites, Held_out), Column],
                loso[grepl("^C2", Mode)][match(sites, Held_out), Depth_applied])

summary_dt <- data.table(
  Configuration = c("ATBD baseline", "Full Pareto site-opt (paper)",
                    "LOSO, ALS-free LUT, own d_opt",
                    "LOSO, ALS-free LUT and depth (mean of training d_opt)",
                    "LOSO, ALS-free LUT and depth (joint depth Pareto)"),
  mean_RMSE = round(c(m_atbd[["RMSE"]], m_pareto[["RMSE"]],
                      m_losoB[["RMSE"]], m_losoC[["RMSE"]], m_losoC2[["RMSE"]]), 2),
  mean_absBias = round(c(m_atbd[["Bias"]], m_pareto[["Bias"]],
                         m_losoB[["Bias"]], m_losoC[["Bias"]], m_losoC2[["Bias"]]), 2))
summary_dt[, `:=`(
  dRMSE_vs_ATBD = round(m_atbd[["RMSE"]] - mean_RMSE, 2),
  dBias_vs_ATBD = round(m_atbd[["Bias"]] - mean_absBias, 2))]
summary_dt[, `:=`(
  pct_RMSE_retained = round(100 * dRMSE_vs_ATBD / (m_atbd[["RMSE"]] - m_pareto[["RMSE"]]), 0),
  pct_Bias_retained = round(100 * dBias_vs_ATBD / (m_atbd[["Bias"]] - m_pareto[["Bias"]]), 0))]
print(summary_dt)
fwrite(summary_dt, file.path(out_dir, "Table_A17b_loso_summary.csv"))
fwrite(alt, file.path(out_dir, "Table_A17c_criterion_robustness.csv"))
