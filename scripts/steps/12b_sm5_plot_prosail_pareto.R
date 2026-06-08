# ---
# title:  12b_sm5_plot_prosail_pareto.R
# desc:   Visualisation — Pareto front plot for optimal PROSAIL configuration
#         selection (SM5 passe 2). For each of two d_opt modes:
#
#           per_site  : each site's own Pareto d_opt (from prosail_opt.csv)
#           all_sites : d_opt from "Sites averaged" in dopt_reference.csv,
#                       applied uniformly to all sites
#
#         Shows all PROSAIL configurations in R² vs RMSE space at the
#         selected depth. Non-dominated (Pareto front) configurations are
#         highlighted in blue; the selected Column_opt is marked with a red
#         star; the ATBD reference config is marked in orange.
#
#         Layout: facet_grid(d_opt_mode ~ site) — 2 rows × 3 columns,
#         norm = DSM_keepTrees.
#
#         Reads:
#           revision/output/intermediate/sm5/prosail_opt.csv
#           revision/output/intermediate/sm5/dopt_reference.csv
#           revision/output/intermediate/sm5/
#             all_results_combined_LIDFa_lai_LMA_BROWN.csv
#
#         Figures written to:
#           revision/output/figures/sm5/prosail_pareto_DSM.pdf
#           revision/output/figures/sm5/prosail_pareto_DSM.png
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/12b_sm5_plot_prosail_pareto.R")
# ---

library("here")
library("data.table")
library("ggplot2")
library("patchwork")
library("scales")
library("cli")

# ── Parameters ─────────────────────────────────────────────────────────────────

norm_select      <- "DSM_keepTrees"
h_min_select     <- 10L
lai_scenario_sel <- "per_site"
sites            <- c("Aigoual", "Blois", "Mormal")

dopt_mode_labels <- c(
  per_site  = "d_opt — per site",
  all_sites = "d_opt — sites averaged"
)

cat_colours <- c(
  dominated      = "grey75",
  `non-dominated` = "#4393c3",
  ATBD           = "#e08214",
  selected       = "#d6604d"
)
cat_sizes <- c(
  dominated      = 1.2,
  `non-dominated` = 2.2,
  ATBD           = 2.8,
  selected       = 2.8
)
cat_alpha <- c(
  dominated      = 0.30,
  `non-dominated` = 0.85,
  ATBD           = 0.90,
  selected       = 1.00
)

# ── Load data ──────────────────────────────────────────────────────────────────

opt_csv  <- file.path(paths$output, "intermediate", "sm5",
                       "prosail_opt.csv")
dref_csv <- file.path(paths$output, "intermediate", "sm5",
                       "dopt_reference.csv")
comb_csv <- file.path(paths$output, "intermediate", "sm5",
                       "all_results_combined_LIDFa_lai_LMA_BROWN.csv")

for (f in c(opt_csv, dref_csv, comb_csv)) {
  if (!file.exists(f)) stop("Missing required file:\n  ", f)
}

cli::cli_alert_info("Loading prosail_opt.csv ...")
opt_dt <- data.table::fread(opt_csv)

cli::cli_alert_info("Loading dopt_reference.csv ...")
dref   <- data.table::fread(dref_csv)

cli::cli_alert_info("Loading all_results_combined (~1M rows, please wait) ...")
combined <- data.table::fread(comb_csv)
combined <- combined[
  Norm         == norm_select &
  h_min        == h_min_select &
  lai_scenario == lai_scenario_sel
]
cli::cli_alert_success("Filtered to {nrow(combined)} rows.")

# ── Pareto helpers ─────────────────────────────────────────────────────────────

identify_pareto_front <- function(m_mat) {
  n         <- nrow(m_mat)
  dominated <- logical(n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      if (all(m_mat[j, ] <= m_mat[i, ]) && any(m_mat[j, ] < m_mat[i, ])) {
        dominated[i] <- TRUE
        break
      }
    }
  }
  !dominated
}

compute_dist_utopia <- function(sub) {
  m1 <- 1 - sub$R;   m2 <- sub$RMSE
  m3 <- abs(sub$Bias); m4 <- abs(sub$Slope - 1)
  valid <- !is.na(m1) & !is.na(m2) & !is.na(m3) & !is.na(m4)
  out   <- rep(NA_real_, nrow(sub))
  if (sum(valid) < 2L) return(out)
  m_v <- cbind(m1[valid], m2[valid], m3[valid], m4[valid])
  m_norm <- apply(m_v, 2L, function(x) {
    lo <- min(x); hi <- max(x)
    if (hi == lo) return(rep(0, length(x)))
    (x - lo) / (hi - lo)
  })
  out[valid] <- sqrt(rowSums(m_norm^2))
  out
}

pareto_select <- function(sub) {
  m_mat <- cbind(1 - sub$R, sub$RMSE, abs(sub$Bias), abs(sub$Slope - 1))
  valid <- complete.cases(m_mat)
  is_front <- logical(nrow(sub))
  if (sum(valid) >= 2L) {
    front_idx <- identify_pareto_front(m_mat[valid, , drop = FALSE])
    is_front[which(valid)[front_idx]] <- TRUE
  }
  dist <- compute_dist_utopia(sub)
  best_idx <- which.min(dist)
  list(
    is_front = is_front,
    dist     = dist,
    col_opt  = if (length(best_idx) > 0L) sub$Column[best_idx] else NA_character_
  )
}

# ── Build combined selection table ─────────────────────────────────────────────

# Per-site mode: d_opt and Column_opt from prosail_opt.csv
per_site_dt <- opt_dt[
  d_opt_source == "per_site" & Norm == norm_select,
  .(Site, d_opt, Column_opt)
]

# All-sites mode: d_opt from "Sites averaged" row in dopt_reference.csv.
# Normalise combined labels to "Sites averaged" (the label used by script 06).
combined_labels <- c("Sites averaged", "Sites_averaged",
                     "Sites combined",  "Sites_combined", "All_sites")

all_sites_dref <- dref[
  Site %in% combined_labels & Norm == norm_select & method_dopt == "pareto",
  .(d_opt = d_opt[[1L]])
]
if (nrow(all_sites_dref) == 0L) {
  cli::cli_warn("No 'Sites averaged' d_opt found in dopt_reference.csv — ",
                "skipping all_sites mode.")
  all_sites_d <- NA_integer_
} else {
  all_sites_d <- all_sites_dref$d_opt
  cli::cli_alert_info("All-sites d_opt ({norm_select}) = {all_sites_d}")
}

# ── Build panel data for each (mode × site) ────────────────────────────────────

build_panel_dt <- function(site_name, d_val, col_opt_ref = NULL) {
  sub <- combined[Site == site_name & Depth == d_val]
  if (nrow(sub) == 0L) {
    cli::cli_warn("No data for Site={site_name}, Depth={d_val} — skipping.")
    return(NULL)
  }
  pf <- pareto_select(sub)

  # If no pre-selected column given, use the Pareto minimum
  col_opt <- if (!is.null(col_opt_ref)) col_opt_ref else pf$col_opt

  sub[, dist_utopia := pf$dist]
  sub[, on_front    := pf$is_front]
  sub[, category    := data.table::fcase(
    Column == col_opt,                   "selected",
    ATBD   == TRUE & Column != col_opt,  "ATBD",
    on_front & Column != col_opt,        "non-dominated",
    default                              = "dominated"
  )]
  sub[, category := factor(
    category, levels = c("dominated", "non-dominated", "ATBD", "selected")
  )]
  data.table::setorder(sub, category)
  sub[, site_name := site_name]
  sub[, col_opt   := col_opt]
  sub
}

panels <- list()

for (s in sites) {
  # per_site
  d_ps  <- per_site_dt[Site == s, d_opt]
  co_ps <- per_site_dt[Site == s, Column_opt]
  if (length(d_ps) > 0L)
    panels[[paste0("per_site_", s)]] <- build_panel_dt(s, d_ps, co_ps)

  # all_sites
  if (!is.na(all_sites_d))
    panels[[paste0("all_sites_", s)]] <- build_panel_dt(s, all_sites_d)
}

panels <- Filter(Negate(is.null), panels)

if (length(panels) == 0L) stop("No panel data — check inputs.")

# ── Single-panel plot function ─────────────────────────────────────────────────

make_pareto_panel <- function(pdata, title = NULL, show_x = TRUE,
                               show_y = TRUE, show_legend = FALSE) {
  col_opt_row <- pdata[category == "selected"]

  p <- ggplot2::ggplot(
    pdata,
    ggplot2::aes(x = R2, y = RMSE,
                 colour = category, size = category, alpha = category)
  ) +
    ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::scale_colour_manual(
      values = cat_colours,
      labels = c(
        dominated       = "Dominated",
        `non-dominated` = "Pareto front",
        ATBD            = "ATBD config",
        selected        = "Selected (opt)"
      ),
      name = NULL
    ) +
    ggplot2::scale_size_manual(values = cat_sizes, guide = "none") +
    ggplot2::scale_alpha_manual(values = cat_alpha, guide = "none") +
    ggplot2::labs(
      title = title,
      x     = if (show_x) expression(R^2) else NULL,
      y     = if (show_y) "RMSE (m²/m²)" else NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5,
                                                size = 11),
      legend.position  = if (show_legend) "right" else "none",
      panel.grid.minor = ggplot2::element_blank()
    )

  # Overlay star for selected config
  if (nrow(col_opt_row) > 0L) {
    p <- p + ggplot2::annotate(
      "point",
      x = col_opt_row$R2[[1L]], y = col_opt_row$RMSE[[1L]],
      shape = 8, size = 5.5, colour = cat_colours[["selected"]], stroke = 1.4
    )
  }

  p
}

# ── Assemble grid ──────────────────────────────────────────────────────────────

dopt_modes <- c("per_site", "all_sites")

row_panels <- list()

for (mode in dopt_modes) {
  mode_panels <- list()
  for (si in seq_along(sites)) {
    s   <- sites[si]
    key <- paste0(mode, "_", s)
    if (!key %in% names(panels)) next
    pd  <- panels[[key]]
    d_v <- unique(pd$Depth)
    ttl <- paste0(s, "  (d_opt = ", d_v, " m)")
    p_i <- make_pareto_panel(
      pd,
      title       = ttl,
      show_x      = (mode == dopt_modes[length(dopt_modes)]),
      show_y      = (si == 1L),
      show_legend = (si == length(sites))
    )
    mode_panels[[s]] <- p_i
  }
  if (length(mode_panels) > 0L)
    row_panels[[mode]] <- mode_panels
}

# Stack rows (per_site top, all_sites bottom)
assembled_rows <- lapply(seq_along(row_panels), function(ri) {
  mode   <- names(row_panels)[ri]
  row_ps <- row_panels[[mode]]
  strip  <- patchwork::wrap_plots(row_ps, nrow = 1L, guides = "collect") &
    ggplot2::theme(legend.position = "none")
  strip + patchwork::plot_annotation(
    title = dopt_mode_labels[[mode]],
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0)
    )
  )
})

fig <- patchwork::wrap_plots(assembled_rows, ncol = 1L) +
  patchwork::plot_annotation(
    title    = paste0(
      "Pareto PROSAIL configuration selection - Norm: ",
      norm_select
    ),
    subtitle = paste0(
      "h_min = ", h_min_select,
      " m | LAI scenario = ", lai_scenario_sel,
      " | Criteria: R2, RMSE, |Bias|, |Slope-1|"
    ),
    caption  = "Red * = selected config  |  Blue = Pareto front  |  Orange = ATBD"
  )

# ── Save ───────────────────────────────────────────────────────────────────────

out_dir <- file.path(paths$output, "figures", "sm5")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

n_rows <- length(assembled_rows)
h_cm   <- 12 * n_rows + 4

out_pdf <- file.path(out_dir, "prosail_pareto_DSM.pdf")
out_png <- file.path(out_dir, "prosail_pareto_DSM.png")

ggplot2::ggsave(out_pdf, fig,
                width = 30, height = h_cm, units = "cm", device = cairo_pdf)
ggplot2::ggsave(out_png, fig,
                width = 30, height = h_cm, units = "cm",
                device = "png", dpi = 600, bg = "white")

cli::cli_alert_success("Written: {out_pdf}")
cli::cli_alert_success("Written: {out_png}")

# ── Parameter extraction helper ─────────────────────────────────────────────────
# Column format: "LIDFa=X_lai=Y_LMA=Z_BROWN=W"

extract_col_params <- function(col_vec) {
  data.table::data.table(
    LIDFa = as.numeric(regmatches(col_vec,
               regexpr("(?<=LIDFa=)[^_]+", col_vec, perl = TRUE))),
    lai   = as.numeric(regmatches(col_vec,
               regexpr("(?<=_lai=)[^_]+",  col_vec, perl = TRUE))),
    LMA   = as.numeric(regmatches(col_vec,
               regexpr("(?<=LMA=)[^_]+",   col_vec, perl = TRUE))),
    BROWN = as.numeric(regmatches(col_vec,
               regexpr("(?<=BROWN=)[^_]+", col_vec, perl = TRUE)))
  )
}

# ── Collect Front-1 rows for per_site mode ─────────────────────────────────────

front1_list <- vector("list", length(sites))
atbd_list   <- vector("list", length(sites))
names(front1_list) <- names(atbd_list) <- sites

for (s in sites) {
  d_ps  <- per_site_dt[Site == s, d_opt]
  co_ps <- per_site_dt[Site == s, Column_opt]
  if (length(d_ps) == 0L) next

  pd <- build_panel_dt(s, d_ps, co_ps)
  if (is.null(pd)) next

  front1_list[[s]] <- pd[on_front == TRUE]
  atbd_list[[s]]   <- pd[ATBD == TRUE]
}

front1_all <- data.table::rbindlist(front1_list, use.names = TRUE, fill = TRUE)
atbd_all   <- data.table::rbindlist(atbd_list,   use.names = TRUE, fill = TRUE)

if (nrow(front1_all) == 0L || !"site_name" %in% names(front1_all)) {
  cli::cli_warn(
    "front1_all is empty — run 12_sm5_select_prosail_opt.R first to populate prosail_opt.csv."
  )
  stop("Stopping: no Front-1 data available.", call. = FALSE)
}

# ── Console summary: Front-1 ranges and ATBD values ───────────────────────────

cli::cli_h2("Front-1 metric ranges vs ATBD")
metrics_disp <- c("R", "RMSE", "Bias", "Slope")

for (s in sites) {
  cli::cli_h3(s)
  f1   <- front1_all[site_name == s]
  atbd <- atbd_all[site_name == s]
  for (m in metrics_disp) {
    rng      <- range(f1[[m]], na.rm = TRUE)
    atbd_val <- if (nrow(atbd) > 0L) round(atbd[[m]][[1L]], 2) else NA
    cli::cli_text(
      "  {m}: Front-1 [{round(rng[1], 2)}, {round(rng[2], 2)}]  |  ATBD = {atbd_val}"
    )
  }
}

# ── Console: Pareto Front-1 ranked table ──────────────────────────────────────

cli::cli_h2("Pareto Front-1 -- ranked by avg normalised criteria")

front1_param_labels <- list(
  LIDFa = c("1"="ATBD",  "2"="OPT#1", "3"="OPT#2",
             "4"="OPT#3", "5"="OPT#4"),
  lai   = c("1"="ATBD",  "2"="OPT#1",         "3"="OPT#2",
             "4"="OPT#3", "5"="LAIALS",         "6"="LAIALS_dopt"),
  LMA   = c("1"="ATBD",  "2"="OPT#1", "3"="OPT#2"),
  BROWN = c("1"="ATBD",  "2"="OPT#1", "3"="OPT#2")
)

fwd_range <- function(x) {
  lo <- min(x, na.rm = TRUE); hi <- max(x, na.rm = TRUE)
  if (hi <= lo) return(rep(0, length(x)))
  (x - lo) / (hi - lo)
}
inv_range <- function(x) {
  lo <- min(x, na.rm = TRUE); hi <- max(x, na.rm = TRUE)
  if (hi <= lo) return(rep(0, length(x)))
  (hi - x) / (hi - lo)
}

front1_ranked_list <- list()

for (s in sites) {
  cli::cli_h3("{s}")
  f1 <- data.table::copy(front1_all[site_name == s])
  if (nrow(f1) == 0L) { message("  (no data)"); next }

  params <- extract_col_params(f1$Column)
  f1[, ALA_lbl   := front1_param_labels$LIDFa[as.character(params$LIDFa)]]
  f1[, LAI_lbl   := front1_param_labels$lai[  as.character(params$lai)]]
  f1[, LMA_lbl   := front1_param_labels$LMA[  as.character(params$LMA)]]
  f1[, BROWN_lbl := front1_param_labels$BROWN[as.character(params$BROWN)]]

  f1[, R_norm     := fwd_range(R)]
  f1[, RMSE_norm  := inv_range(RMSE)]
  f1[, Bias_norm  := inv_range(abs(Bias))]
  f1[, Slope_norm := inv_range(abs(Slope - 1))]
  f1[, Avg_norm   := rowMeans(.SD),
     .SDcols = c("R_norm", "RMSE_norm", "Bias_norm", "Slope_norm")]

  data.table::setorder(f1, -Avg_norm)
  f1[, Rank := .I]

  front1_ranked_list[[s]] <- f1

  hdr <- sprintf(
    "  %-8s %-12s %-6s %-6s  %5s %6s %7s %6s  %6s %9s %9s %10s  %8s  %s",
    "ALA", "LAI", "LMA", "BROWN", "R", "RMSE", "Bias", "Slope",
    "R_n", "RMSE_n", "Bias_n", "Slope_n", "Avg_n", "Rank"
  )
  message(hdr)
  message(strrep("-", nchar(hdr)))

  for (i in seq_len(nrow(f1))) {
    row <- f1[i]
    message(sprintf(
      "  %-8s %-12s %-6s %-6s  %5.2f %6.2f %7.2f %6.2f  %6.2f %9.2f %9.2f %10.2f  %8.2f  %d",
      row$ALA_lbl, row$LAI_lbl, row$LMA_lbl, row$BROWN_lbl,
      row$R, row$RMSE, row$Bias, row$Slope,
      row$R_norm, row$RMSE_norm, row$Bias_norm, row$Slope_norm,
      row$Avg_norm, row$Rank
    ))
  }
}

front1_ranked_all <- data.table::rbindlist(front1_ranked_list,
                                            use.names = TRUE, fill = TRUE)
ranked_csv <- here::here(
  "revision", "output", "intermediate", "sm5",
  "prosail_pareto_front1_ranked.csv"
)
data.table::fwrite(
  front1_ranked_all[, .(
    site = site_name, ALA = ALA_lbl, LAI = LAI_lbl,
    LMA = LMA_lbl, BROWN = BROWN_lbl,
    R        = round(R,        2),
    RMSE     = round(RMSE,     2),
    Bias     = round(Bias,     2),
    Slope    = round(Slope,    2),
    R_norm   = round(R_norm,   2),
    RMSE_norm  = round(RMSE_norm,  2),
    Bias_norm  = round(Bias_norm,  2),
    Slope_norm = round(Slope_norm, 2),
    Avg_norm   = round(Avg_norm,   2),
    Rank
  )],
  ranked_csv
)
cli::cli_alert_success("Written: {ranked_csv}  ({nrow(front1_ranked_all)} rows)")

# ── Figure 2: Front-1 metric boxplots ─────────────────────────────────────────
# Boxes = distribution of all Front-1 configs per site at per_site d_opt.
# Red ★ = ATBD config value.

cli::cli_h2("Figure 2 — Front-1 metric boxplots")

metrics_long <- data.table::melt(
  front1_all,
  id.vars       = "site_name",
  measure.vars  = c("R", "RMSE", "Bias", "Slope"),
  variable.name = "Metric",
  value.name    = "Value"
)
atbd_long <- data.table::melt(
  atbd_all,
  id.vars       = "site_name",
  measure.vars  = c("R", "RMSE", "Bias", "Slope"),
  variable.name = "Metric",
  value.name    = "Value"
)
metrics_long[, site_name := factor(site_name, levels = sites)]
atbd_long[,   site_name := factor(site_name, levels = sites)]

ylims_lower <- c(R = 0.4, RMSE = 0, Bias = -1, Slope = 0)
ylabels <- list(
  R     = expression(italic(r)),
  RMSE  = expression(RMSE~(m^2/m^2)),
  Bias  = expression(Bias~(m^2/m^2)),
  Slope = "Slope"
)

plots_box <- lapply(c("R", "RMSE", "Bias", "Slope"), function(m) {
  df_m   <- metrics_long[Metric == m]
  atbd_m <- atbd_long[Metric == m]

  ggplot2::ggplot(df_m, ggplot2::aes(x = site_name, y = Value)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_point(
      data  = atbd_m,
      ggplot2::aes(x = site_name, y = Value),
      shape = 8, size = 3.2, stroke = 1.2, colour = "red"
    ) +
    ggplot2::coord_cartesian(ylim = c(ylims_lower[[m]], NA)) +
    ggplot2::labs(x = NULL, y = ylabels[[m]]) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10))
})

fig_box <- (plots_box[[1L]] | plots_box[[2L]]) /
  (plots_box[[3L]] | plots_box[[4L]]) +
  patchwork::plot_annotation(
    title = "Pareto Front 1 - metric distributions by site  (* = ATBD config)",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0.5)
    )
  )

out_box_pdf <- file.path(out_dir, "prosail_front1_metrics_DSM.pdf")
out_box_png <- file.path(out_dir, "prosail_front1_metrics_DSM.png")
ggplot2::ggsave(out_box_pdf, fig_box,
                width = 22, height = 20, units = "cm", device = cairo_pdf)
ggplot2::ggsave(out_box_png, fig_box,
                width = 22, height = 20, units = "cm", device = "png", dpi = 600, bg = "white")
cli::cli_alert_success("Written: {out_box_pdf}")
cli::cli_alert_success("Written: {out_box_png}")

# ── Figure 3: Parameter frequency in Front-1 ───────────────────────────────────
# facet_grid(site ~ parameter): parameters as columns, sites as rows.
# Each panel shows proportion of each param value in Front-1 configs.
# Red star marks the value of the selected (opt) config for that site.
# Dashed line = uniform baseline (1 / n_distinct_values).

cli::cli_h2("Figure 3 — Front-1 parameter frequency")

# Internal extraction names match the Column string; display labels differ.
params_extract <- c("LIDFa", "lai",  "LMA",   "BROWN")
param_grid_raw <- list(LIDFa = 1:5,  lai = 1:6, LMA = 1:3, BROWN = 1:3)
param_labels   <- c(LIDFa = "ALA",   lai = "LAI", LMA = "LMA", BROWN = "BROWN")

# Display order: alphabetical on renamed labels
params_disp    <- sort(unique(param_labels))          # ALA, BROWN, LAI, LMA
param_grid_disp <- setNames(
  param_grid_raw[names(param_labels)[order(param_labels)]],
  params_disp
)

front1_pcols <- data.table::copy(front1_all[, .(site_name, Column)])
extracted     <- extract_col_params(front1_pcols$Column)
front1_pcols  <- cbind(front1_pcols, extracted)

freq_long <- data.table::melt(
  front1_pcols,
  id.vars       = "site_name",
  measure.vars  = params_extract,
  variable.name = "Parameter",
  value.name    = "Value"
)
freq_long[, Parameter := param_labels[as.character(Parameter)]]

freq_tbl <- freq_long[, .N, by = .(site_name, Parameter, Value)]
freq_tbl[, freq := N / sum(N), by = .(site_name, Parameter)]

# Full (site × parameter × value) grid — fill missing with 0
full_grid <- data.table::rbindlist(lapply(params_disp, function(p) {
  data.table::CJ(site_name = sites, Parameter = p, Value = param_grid_disp[[p]])
}))
freq_full <- merge(full_grid, freq_tbl, by = c("site_name", "Parameter", "Value"),
                   all.x = TRUE)
freq_full[is.na(freq), freq := 0]
freq_full[is.na(N),    N    := 0L]
freq_full[, site_name := factor(site_name, levels = sites)]
freq_full[, Parameter := factor(Parameter, levels = params_disp)]

uniform_dt <- data.table::data.table(
  Parameter    = factor(params_disp, levels = params_disp),
  uniform_freq = 1 / lengths(param_grid_disp)
)

# Selected config: one param value per (site, parameter)
selected_params <- data.table::rbindlist(lapply(sites, function(s) {
  co <- per_site_dt[Site == s, Column_opt]
  if (length(co) == 0L) return(NULL)
  ex <- extract_col_params(co)
  data.table::melt(
    data.table::data.table(site_name = s, ex),
    id.vars       = "site_name",
    variable.name = "Parameter",
    value.name    = "Value"
  )
}))
selected_params[, Parameter := param_labels[as.character(Parameter)]]
selected_params[, site_name := factor(site_name, levels = sites)]
selected_params[, Parameter := factor(Parameter, levels = params_disp)]

# Merge to get bar height for star positioning
star_dt <- merge(selected_params, freq_full,
                 by = c("site_name", "Parameter", "Value"))

# Per-parameter x-axis label map.
# Label strings are R plotmath expression syntax:
#   quoted literals ('ATBD') for plain text, LAI["ALS"] for subscripts.
# LIDFa (displayed as ALA) has 5 values → value 5 = OPT#4, not LAI.
x_label_map <- list(
  ALA   = c("1" = "'ATBD'",  "2" = "'OPT#1'", "3" = "'OPT#2'",
            "4" = "'OPT#3'", "5" = "'OPT#4'"),
  LAI   = c("1" = "'ATBD'",  "2" = "'OPT#1'", "3" = "'OPT#2'",
            "4" = "'OPT#3'", "5" = 'LAI["ALS"]', "6" = 'LAI["ALS_dopt"]'),
  LMA   = c("1" = "'ATBD'",  "2" = "'OPT#1'", "3" = "'OPT#2'"),
  BROWN = c("1" = "'ATBD'",  "2" = "'OPT#1'", "3" = "'OPT#2'")
)

all_x_levels <- c("'ATBD'",  "'OPT#1'", "'OPT#2'", "'OPT#3'", "'OPT#4'",
                  'LAI["ALS"]', 'LAI["ALS_dopt"]')

freq_full[, Value_label := x_label_map[[as.character(Parameter)]][as.character(Value)],
          by = Parameter]
freq_full[, Value_label := factor(Value_label, levels = all_x_levels)]

star_dt[,   Value_label := x_label_map[[as.character(Parameter)]][as.character(Value)],
          by = Parameter]
star_dt[,   Value_label := factor(Value_label, levels = all_x_levels)]

fig_freq <- ggplot2::ggplot(
    freq_full,
    ggplot2::aes(x = Value_label, y = freq)
  ) +
  ggplot2::geom_col(fill = "steelblue", width = 0.7) +
  ggplot2::geom_text(
    data        = freq_full[N > 0L],
    ggplot2::aes(x = Value_label, y = freq * 0.65, label = N),
    colour      = "white", size = 3.5, fontface = "bold", vjust = 0.3,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_point(
    data        = star_dt,
    ggplot2::aes(x = Value_label, y = freq),
    shape       = 8, size = 5.0, stroke = 1.5, colour = "red",
    inherit.aes = FALSE,
    position    = ggplot2::position_nudge(y = 0.085)
  ) +
  ggplot2::geom_hline(
    data        = uniform_dt,
    ggplot2::aes(yintercept = uniform_freq),
    linetype    = "dashed", linewidth = 0.7, colour = "black",
    inherit.aes = FALSE
  ) +
  ggplot2::scale_x_discrete(labels = scales::label_parse(), drop = TRUE) +
  ggplot2::facet_grid(site_name ~ Parameter, scales = "free_x") +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.15))
  ) +
  ggplot2::labs(
    x = "Parameter Configuration",
    y = "Proportion of Parameter Configurations"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
    strip.text       = ggplot2::element_text(face = "bold"),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 30, hjust = 1)
  )

out_freq_pdf <- file.path(out_dir, "prosail_front1_params_DSM.pdf")
out_freq_png <- file.path(out_dir, "prosail_front1_params_DSM.png")
ggplot2::ggsave(out_freq_pdf, fig_freq,
                width = 28, height = 20, units = "cm", device = cairo_pdf)
ggplot2::ggsave(out_freq_png, fig_freq,
                width = 28, height = 20, units = "cm", device = "png", dpi = 600, bg = "white")
cli::cli_alert_success("Written: {out_freq_pdf}")
cli::cli_alert_success("Written: {out_freq_png}")

cli::cli_h1("Done — SM5 Pareto PROSAIL selection figures")
