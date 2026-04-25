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

opt_csv  <- here::here("revision", "output", "intermediate", "sm5",
                       "prosail_opt.csv")
dref_csv <- here::here("revision", "output", "intermediate", "sm5",
                       "dopt_reference.csv")
comb_csv <- here::here("revision", "output", "intermediate", "sm5",
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
  m1 <- 1 - sub$R2;  m2 <- sub$RMSE
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
  m_mat <- cbind(1 - sub$R2, sub$RMSE, abs(sub$Bias), abs(sub$Slope - 1))
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
      "Pareto PROSAIL configuration selection — Norm: ",
      norm_select
    ),
    subtitle = paste0(
      "h_min = ", h_min_select,
      " m | LAI scenario = ", lai_scenario_sel,
      " | Criteria: R², RMSE, |Bias|, |Slope−1|"
    ),
    caption  = "Red ★ = selected config  |  Blue = Pareto front  |  Orange = ATBD"
  )

# ── Save ───────────────────────────────────────────────────────────────────────

out_dir <- here::here("revision", "output", "figures", "sm5")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

n_rows <- length(assembled_rows)
h_cm   <- 12 * n_rows + 4

out_pdf <- file.path(out_dir, "prosail_pareto_DSM.pdf")
out_png <- file.path(out_dir, "prosail_pareto_DSM.png")

ggplot2::ggsave(out_pdf, fig,
                width = 30, height = h_cm, units = "cm")
ggplot2::ggsave(out_png, fig,
                width = 30, height = h_cm, units = "cm",
                device = "png", dpi = 300)

cli::cli_alert_success("Written: {out_pdf}")
cli::cli_alert_success("Written: {out_png}")
cli::cli_h1("Done — SM5 Pareto PROSAIL selection figure")
