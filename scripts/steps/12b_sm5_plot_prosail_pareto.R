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
#           output/intermediate/sm5/prosail_opt.csv
#           output/intermediate/sm5/dopt_reference.csv
#           output/intermediate/sm5/
#             all_results_combined_LIDFa_lai_LMA_BROWN.csv
#
#         Figures written to:
#           output/figures/prosail_pareto_DSM.pdf
#           output/figures/prosail_pareto_DSM.png
#
# Run from the project root (NC_Full/):
#   source("scripts/12b_sm5_plot_prosail_pareto.R")
# ---

library("here")
library("data.table")
library("ggplot2")
library("patchwork")
library("scales")
library("cli")

source(here::here("R", "paths.R"))

# Single source of truth for Pareto front + utopian-distance selection.
# Ensures consistency with `prosail_opt.csv` produced by step 12.
source(here::here("R", "sm5_dopt.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

norm_select      <- "DSM_keepTrees"
h_min_select     <- 10L
lai_scenario_sel <- "common"

# Pareto criteria — must match scripts/steps/12_sm5_select_prosail_opt.R
pareto_criteria    <- c("R", "RMSE", "Bias", "Slope")
pareto_norm_method <- "minmax"
# pareto_norm_method <- "max"
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
# Pareto front identification and Avg_score (mean of Front-1-normalised
# criteria) delegate to compute_pareto_front() in R/sm5_dopt.R so that this
# figure script and prosail_opt.csv (step 12) stay strictly coherent (same
# front, same selection, same Avg_score).

pareto_select <- function(sub) {
  if (nrow(sub) < 2L) {
    return(list(is_front  = logical(nrow(sub)),
                avg_score = rep(NA_real_, nrow(sub)),
                col_opt   = NA_character_))
  }

  # compute_pareto_front() takes a metrics_subset with columns
  # Depth, R, RMSE, Bias, Slope. We use row index as a stand-in for Depth.
  subset_for_pf <- data.table::data.table(
    Depth = seq_len(nrow(sub)),
    R     = sub$R,
    RMSE  = sub$RMSE,
    Bias  = sub$Bias,
    Slope = sub$Slope
  )
  pf <- compute_pareto_front(subset_for_pf,
                              criteria    = pareto_criteria,
                              norm_method = pareto_norm_method,
                              norm_scope  = "front1")

  is_front <- logical(nrow(sub))
  if (!anyNA(pf$pareto_front_depths))
    is_front[pf$pareto_front_depths] <- TRUE

  # Map all_distances (Avg_score named numeric, NA for dominated configs)
  # back to a per-row vector.
  avg_v <- rep(NA_real_, nrow(sub))
  if (!is.null(pf$all_distances)) {
    idx_keep      <- as.integer(names(pf$all_distances))
    avg_v[idx_keep] <- as.numeric(pf$all_distances)
  }

  col_opt <- if (!is.na(pf$d_opt)) sub$Column[pf$d_opt] else NA_character_

  list(is_front = is_front, avg_score = avg_v, col_opt = col_opt)
}

# ── Build combined selection table ─────────────────────────────────────────────

# Per-site mode: d_opt and Column_opt from prosail_opt.csv
per_site_dt <- opt_dt[
  d_opt_source == "per_site" & Norm == norm_select,
  .(Site, d_opt, Column_opt)
]

# All-sites mode: d_opt and per-site Column_opt from prosail_opt.csv too
all_sites_dt <- opt_dt[
  d_opt_source == "all_sites" & Norm == norm_select,
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

## Per-criterion normalisation helper — driven by `pareto_norm_method`.
## "higher = better" convention: returns 1 - minmax(x) so that rowMeans
## of the per-criterion *_norm columns equals Avg_score (1 = best, 0 = worst).
## Applied on the FULL candidate set inside build_panel_dt below.
inv_range <- function(x) {
  lo <- min(x, na.rm = TRUE); hi <- max(x, na.rm = TRUE)
  if (pareto_norm_method == "max") {
    if (hi == 0) return(rep(1, length(x)))
    return(1 - x / hi)
  }
  if (hi <= lo) return(rep(1, length(x)))
  1 - (x - lo) / (hi - lo)
}

build_panel_dt <- function(site_name, d_val, col_opt_ref = NULL) {
  sub <- combined[Site == site_name & Depth == d_val]
  if (nrow(sub) == 0L) {
    cli::cli_warn("No data for Site={site_name}, Depth={d_val} — skipping.")
    return(NULL)
  }
  pf <- pareto_select(sub)

  # If no pre-selected column given, use the Pareto minimum
  col_opt <- if (!is.null(col_opt_ref)) col_opt_ref else pf$col_opt

  sub[, Avg_score := pf$avg_score]
  sub[, on_front  := pf$is_front]

  # Per-criterion min-max norms on the FRONT-1 subset only (mirrors
  # compute_pareto_front's norm_scope = "front1" for PROSAIL OPT, which
  # ignores absurd configurations among the ~225 grid points). Rows outside
  # Front-1 keep their norm columns as NA — the ranked CSV only retains
  # Front-1 rows.
  front_mask <- sub$on_front
  if (any(front_mask)) {
    R_min     <- 1 - round(sub$R[front_mask],     2)
    RMSE_min  <- round(sub$RMSE[front_mask],      2)
    Bias_min  <- round(abs(round(sub$Bias[front_mask],  2)),     2)
    Slope_min <- round(abs(round(sub$Slope[front_mask], 2) - 1), 2)

    sub[, R_norm     := NA_real_]
    sub[, RMSE_norm  := NA_real_]
    sub[, Bias_norm  := NA_real_]
    sub[, Slope_norm := NA_real_]
    sub[front_mask, R_norm     := inv_range(R_min)]
    sub[front_mask, RMSE_norm  := inv_range(RMSE_min)]
    sub[front_mask, Bias_norm  := inv_range(Bias_min)]
    sub[front_mask, Slope_norm := inv_range(Slope_min)]
  } else {
    sub[, R_norm     := NA_real_]
    sub[, RMSE_norm  := NA_real_]
    sub[, Bias_norm  := NA_real_]
    sub[, Slope_norm := NA_real_]
  }

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

out_dir <- file.path(paths$output, "figures")
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

# ── Collect Front-1 rows per mode (per_site, all_sites) ────────────────────────

collect_front1 <- function(mode_dt, mode_label) {
  fl <- vector("list", length(sites)); names(fl) <- sites
  al <- vector("list", length(sites)); names(al) <- sites
  for (s in sites) {
    d_v  <- mode_dt[Site == s, d_opt]
    co_v <- mode_dt[Site == s, Column_opt]
    if (length(d_v) == 0L) next
    pd <- build_panel_dt(s, d_v, co_v)
    if (is.null(pd)) next
    pd[, dopt_source := mode_label]
    fl[[s]] <- pd[on_front == TRUE]
    al[[s]] <- pd[ATBD == TRUE]
  }
  list(front1 = data.table::rbindlist(fl, use.names = TRUE, fill = TRUE),
       atbd   = data.table::rbindlist(al, use.names = TRUE, fill = TRUE))
}

mode_data <- list(
  per_site  = collect_front1(per_site_dt,  "per_site"),
  all_sites = collect_front1(all_sites_dt, "all_sites")
)

# Default front1_all / atbd_all = per_site (used by Figure 2 boxplots that
# rely on the per_site d_opt; the all_sites view is in the ranked CSV).
front1_all <- mode_data$per_site$front1
atbd_all   <- mode_data$per_site$atbd

if (nrow(front1_all) == 0L || !"site_name" %in% names(front1_all)) {
  cli::cli_warn(
    "front1_all is empty — run 12_sm5_select_prosail_opt.R first to populate prosail_opt.csv."
  )
  stop("Stopping: no Front-1 data available.", call. = FALSE)
}

# ── Console summary: Front-1 ranges and ATBD values ───────────────────────────
# Bias and Slope have two views:
#   - signed range  (raw)     : sign tells over/under estimation
#   - optimisation range      : |Bias|  and  |Slope - 1|  — these are the
#                                criteria Pareto actually minimises.
# We print both so the reader can interpret a negative-Bias member correctly.

cli::cli_h2("Front-1 metric ranges vs ATBD")

fmt_range <- function(x) sprintf("[%.2f, %.2f]", min(x, na.rm=TRUE),
                                  max(x, na.rm=TRUE))

for (s in sites) {
  cli::cli_h3(s)
  f1   <- front1_all[site_name == s]
  atbd <- atbd_all[site_name == s]
  atbd_R     <- if (nrow(atbd) > 0L) atbd$R[[1L]]     else NA_real_
  atbd_RMSE  <- if (nrow(atbd) > 0L) atbd$RMSE[[1L]]  else NA_real_
  atbd_Bias  <- if (nrow(atbd) > 0L) atbd$Bias[[1L]]  else NA_real_
  atbd_Slope <- if (nrow(atbd) > 0L) atbd$Slope[[1L]] else NA_real_

  cli::cli_text("  R         : Front-1 {fmt_range(f1$R)}    |  ATBD = {round(atbd_R, 2)}")
  cli::cli_text("  RMSE      : Front-1 {fmt_range(f1$RMSE)} |  ATBD = {round(atbd_RMSE, 2)}")
  cli::cli_text(
    "  Bias raw  : Front-1 {fmt_range(f1$Bias)}      |  ATBD = {round(atbd_Bias, 2)}"
  )
  cli::cli_text(
    "  |Bias|    : Front-1 {fmt_range(abs(f1$Bias))} |  ATBD = {round(abs(atbd_Bias), 2)}   (Pareto criterion)"
  )
  cli::cli_text(
    "  Slope raw : Front-1 {fmt_range(f1$Slope)}     |  ATBD = {round(atbd_Slope, 2)}"
  )
  cli::cli_text(
    "  |Slope-1| : Front-1 {fmt_range(abs(f1$Slope - 1))} |  ATBD = {round(abs(atbd_Slope - 1), 2)}   (Pareto criterion)"
  )
}

# ── Console: Pareto Front-1 ranked table ──────────────────────────────────────

cli::cli_h2("Pareto Front-1 -- ranked by avg normalised criteria")

front1_param_labels <- list(
  LIDFa = c("1"="ATBD",  "2"="OPT#1", "3"="OPT#2", "4"="OPT#3", "5"="OPT#4"),
  lai   = c("1"="ATBD",  "2"="OPT#1", "3"="OPT#2",
             "4"="LAIALS", "5"="LAIALS_dopt"),
  LMA   = c("1"="ATBD",  "2"="OPT#1", "3"="OPT#2"),
  BROWN = c("1"="ATBD",  "2"="OPT#1", "3"="OPT#2")
)

## Build the ranked CSV (one per mode) — same logic per mode.
## inv_range() is defined above (before build_panel_dt) so it is available
## both inside build_panel_dt (for the *_norm columns over ALL candidates)
## and here.
build_ranked <- function(front1_all) {
  rl <- list()
  for (s in sites) {
    f1 <- data.table::copy(front1_all[site_name == s])
    if (nrow(f1) == 0L) next

    params <- extract_col_params(f1$Column)
    f1[, ALA_lbl   := front1_param_labels$LIDFa[as.character(params$LIDFa)]]
    f1[, LAI_lbl   := front1_param_labels$lai[  as.character(params$lai)]]
    f1[, LMA_lbl   := front1_param_labels$LMA[  as.character(params$LMA)]]
    f1[, BROWN_lbl := front1_param_labels$BROWN[as.character(params$BROWN)]]

    # *_norm columns are inherited from build_panel_dt — they were computed
    # over the FULL candidate set (~225 configs at d_opt) so that
    # rowMeans(*_norm) reproduces Avg_score for the Front-1 subset.
    data.table::setorder(f1, -Avg_score)
    f1[, Rank := .I]
    rl[[s]] <- f1
  }
  data.table::rbindlist(rl, use.names = TRUE, fill = TRUE)
}

print_ranked <- function(ranked_dt, mode_label) {
  cli::cli_h2("Ranked Front-1 — d_opt source: {mode_label}")
  hdr <- sprintf(
    "  %-8s %-12s %-6s %-6s  %5s %6s %7s %6s  %6s %9s %9s %10s  %9s  %s",
    "ALA", "LAI", "LMA", "BROWN", "R", "RMSE", "Bias", "Slope",
    "R_n", "RMSE_n", "Bias_n", "Slope_n", "Avg_score", "Rank"
  )
  for (s in sites) {
    cli::cli_h3(s)
    sub <- ranked_dt[site_name == s]
    if (nrow(sub) == 0L) { message("  (no data)"); next }
    message(hdr); message(strrep("-", nchar(hdr)))
    for (i in seq_len(nrow(sub))) {
      row <- sub[i]
      message(sprintf(
        "  %-8s %-12s %-6s %-6s  %5.2f %6.2f %7.2f %6.2f  %6.2f %9.2f %9.2f %10.2f  %9.2f  %d",
        row$ALA_lbl, row$LAI_lbl, row$LMA_lbl, row$BROWN_lbl,
        row$R, row$RMSE, row$Bias, row$Slope,
        row$R_norm, row$RMSE_norm, row$Bias_norm, row$Slope_norm,
        row$Avg_score, row$Rank
      ))
    }
  }
}

write_ranked_csv <- function(ranked_dt, mode_label) {
  out_csv <- here::here(
    "output", "intermediate", "sm5",
    paste0("prosail_pareto_front1_ranked_", mode_label, ".csv")
  )
  data.table::fwrite(
    ranked_dt[, .(
      site = site_name, ALA = ALA_lbl, LAI = LAI_lbl,
      LMA = LMA_lbl, BROWN = BROWN_lbl,
      R          = round(R,         2),
      RMSE       = round(RMSE,      2),
      Bias       = round(Bias,      2),
      Slope      = round(Slope,     2),
      R_norm     = round(R_norm,    2),
      RMSE_norm  = round(RMSE_norm, 2),
      Bias_norm  = round(Bias_norm, 2),
      Slope_norm = round(Slope_norm, 2),
      Avg_score  = round(Avg_score, 2),
      Rank
    )],
    out_csv
  )
  cli::cli_alert_success("Written: {out_csv}  ({nrow(ranked_dt)} rows)")
}

for (current_mode in c("per_site", "all_sites")) {
  ranked <- build_ranked(mode_data[[current_mode]]$front1)
  print_ranked(ranked, current_mode)
  write_ranked_csv(ranked, current_mode)
}

# ── Figure 2: Front-1 metric boxplots (one per mode) ──────────────────────────
# Boxes = distribution of all Front-1 configs per site at that mode's d_opt.
# Red ★ = ATBD config value.

cli::cli_h2("Figure 2 — Front-1 metric boxplots (per mode)")

ylims_lower <- c(R = 0.4, RMSE = 0, Bias = -1, Slope = 0)
ylabels <- list(
  R     = expression(italic(r)),
  RMSE  = expression(RMSE~(m^2/m^2)),
  Bias  = expression(Bias~(m^2/m^2)),
  Slope = "Slope"
)

make_fig_box <- function(front1_dt, atbd_dt, mode_label) {
  metrics_long <- data.table::melt(
    front1_dt, id.vars = "site_name",
    measure.vars = c("R", "RMSE", "Bias", "Slope"),
    variable.name = "Metric", value.name = "Value"
  )
  atbd_long <- data.table::melt(
    atbd_dt, id.vars = "site_name",
    measure.vars = c("R", "RMSE", "Bias", "Slope"),
    variable.name = "Metric", value.name = "Value"
  )
  metrics_long[, site_name := factor(site_name, levels = sites)]
  atbd_long[,   site_name := factor(site_name, levels = sites)]

  plots <- lapply(c("R", "RMSE", "Bias", "Slope"), function(m) {
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

  (plots[[1L]] | plots[[2L]]) / (plots[[3L]] | plots[[4L]]) +
    patchwork::plot_annotation(
      title = paste0("Pareto Front 1 — metric distributions  (* = ATBD)  ",
                     "[d_opt: ", mode_label, "]"),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0.5)
      )
    )
}

for (current_mode in c("per_site", "all_sites")) {
  fig_box <- make_fig_box(mode_data[[current_mode]]$front1,
                          mode_data[[current_mode]]$atbd,
                          current_mode)
  out_box_pdf <- file.path(out_dir,
    paste0("prosail_front1_metrics_DSM_", current_mode, ".pdf"))
  out_box_png <- file.path(out_dir,
    paste0("prosail_front1_metrics_DSM_", current_mode, ".png"))
  ggplot2::ggsave(out_box_pdf, fig_box, width = 22, height = 20, units = "cm",
                  device = cairo_pdf)
  ggplot2::ggsave(out_box_png, fig_box, width = 22, height = 20, units = "cm",
                  device = "png", dpi = 600, bg = "white")
  cli::cli_alert_success("Written: {out_box_pdf}")
  cli::cli_alert_success("Written: {out_box_png}")
}

# ── Figure 3: Parameter frequency in Front-1 (one per mode) ───────────────────
# facet_grid(site ~ parameter): parameters as columns, sites as rows.
# Each panel shows proportion of each param value in Front-1 configs.
# Red star marks the value of the selected (opt) config for that site.
# Dashed line = uniform baseline (1 / n_distinct_values).

cli::cli_h2("Figure 3 — Front-1 parameter frequency (per mode)")

# Internal extraction names match the Column string; display labels differ.
params_extract <- c("LIDFa", "lai",  "LMA",   "BROWN")
param_grid_raw <- list(LIDFa = 1:5,  lai = 1:5, LMA = 1:3, BROWN = 1:3)
param_labels   <- c(LIDFa = "ALA",   lai = "LAI", LMA = "LMA", BROWN = "BROWN")

# Display order: alphabetical on renamed labels
params_disp    <- sort(unique(param_labels))          # ALA, BROWN, LAI, LMA
param_grid_disp <- setNames(
  param_grid_raw[names(param_labels)[order(param_labels)]],
  params_disp
)

make_fig_freq <- function(front1_dt, mode_dt, mode_label) {
  front1_pcols <- data.table::copy(front1_dt[, .(site_name, Column)])
  extracted    <- extract_col_params(front1_pcols$Column)
  front1_pcols <- cbind(front1_pcols, extracted)

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

  full_grid <- data.table::rbindlist(lapply(params_disp, function(p) {
    data.table::CJ(site_name = sites, Parameter = p,
                    Value = param_grid_disp[[p]])
  }))
  freq_full <- merge(full_grid, freq_tbl,
                     by = c("site_name", "Parameter", "Value"), all.x = TRUE)
  freq_full[is.na(freq), freq := 0]
  freq_full[is.na(N),    N    := 0L]
  freq_full[, site_name := factor(site_name, levels = sites)]
  freq_full[, Parameter := factor(Parameter, levels = params_disp)]

  uniform_dt <- data.table::data.table(
    Parameter    = factor(params_disp, levels = params_disp),
    uniform_freq = 1 / lengths(param_grid_disp)
  )

  # Selected config: from mode_dt (per_site_dt or all_sites_dt)
  selected_params <- data.table::rbindlist(lapply(sites, function(s) {
    co <- mode_dt[Site == s, Column_opt]
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

  star_dt <- merge(selected_params, freq_full,
                   by = c("site_name", "Parameter", "Value"))

  x_label_map <- list(
    ALA   = c("1" = "'ATBD'",  "2" = "'OPT#1'", "3" = "'OPT#2'",
              "4" = "'OPT#3'", "5" = "'OPT#4'"),
    LAI   = c("1" = "'ATBD'",  "2" = "'OPT#1'", "3" = "'OPT#2'",
              "4" = 'LAI["ALS"]', "5" = 'LAI["ALS_dopt"]'),
    LMA   = c("1" = "'ATBD'",  "2" = "'OPT#1'", "3" = "'OPT#2'"),
    BROWN = c("1" = "'ATBD'",  "2" = "'OPT#1'", "3" = "'OPT#2'")
  )
  all_x_levels <- c("'ATBD'", "'OPT#1'", "'OPT#2'", "'OPT#3'", "'OPT#4'",
                    'LAI["ALS"]', 'LAI["ALS_dopt"]')

  freq_full[, Value_label := x_label_map[[as.character(Parameter)]][as.character(Value)],
            by = Parameter]
  freq_full[, Value_label := factor(Value_label, levels = all_x_levels)]
  star_dt[,   Value_label := x_label_map[[as.character(Parameter)]][as.character(Value)],
            by = Parameter]
  star_dt[,   Value_label := factor(Value_label, levels = all_x_levels)]

  ggplot2::ggplot(freq_full, ggplot2::aes(x = Value_label, y = freq)) +
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
      y = "Proportion of Parameter Configurations",
      subtitle = paste0("d_opt source: ", mode_label)
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      strip.text       = ggplot2::element_text(face = "bold"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(angle = 30, hjust = 1)
    )
}

mode_dts <- list(per_site = per_site_dt, all_sites = all_sites_dt)
for (current_mode in c("per_site", "all_sites")) {
  fig_freq <- make_fig_freq(mode_data[[current_mode]]$front1,
                            mode_dts[[current_mode]],
                            current_mode)
  out_freq_pdf <- file.path(out_dir,
    paste0("prosail_front1_params_DSM_", current_mode, ".pdf"))
  out_freq_png <- file.path(out_dir,
    paste0("prosail_front1_params_DSM_", current_mode, ".png"))
  ggplot2::ggsave(out_freq_pdf, fig_freq, width = 28, height = 20, units = "cm",
                  device = cairo_pdf)
  ggplot2::ggsave(out_freq_png, fig_freq, width = 28, height = 20, units = "cm",
                  device = "png", dpi = 600, bg = "white")
  cli::cli_alert_success("Written: {out_freq_pdf}")
  cli::cli_alert_success("Written: {out_freq_png}")
}

cli::cli_h1("Done — SM5 Pareto PROSAIL selection figures")
