# ---
# title:  22_sm5_plot_prosail_param_distributions.R
# desc:   Figures — parameter distributions for the PROSAIL optimisation steps.
#
#           Panel (a): LAI configurations
#             - Theoretical: ATBD, OPT#1, OPT#2, OPT#3
#             - Empirical (LiDAR): LAIALS and LAIALS_dopt per site + pooled
#           Panel (b): ALA, BROWN, LMA configurations
#             - ATBD (from prosail::get_atbd_lut_input()) and OPT#1–4
#             - Sub-panels: ALA | BROWN | LMA
#
#         Reads:
#           output/intermediate/sm5/prosail_opt.csv
#           03_RESULTS/{site}/Metrics/Deciduous_Only/ladstack_classic.tif
#           03_RESULTS/{site}/Metrics/Deciduous_Only/max_res_10_m.tif
#
#         Outputs:
#           output/figures/sm_prosail_param_distribs.pdf/.png
#
# Run from the project root (NC_Full/):
#   source("scripts/22_sm5_plot_prosail_param_distributions.R")
# ---

library(here)
library(data.table)
library(ggplot2)
library(patchwork)
library(terra)
library(truncnorm)
library(prosail)
library(cli)

source(here::here("R", "paths.R"))

set.seed(42)
sites <- c("Aigoual", "Blois", "Mormal")

# ── I/O ────────────────────────────────────────────────────────────────────────

out_dir <- file.path(paths$output, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── d_opt per site ─────────────────────────────────────────────────────────────

opt_csv <- file.path(paths$output, "intermediate", "sm5", "prosail_opt.csv")
if (!file.exists(opt_csv))
  stop("Missing prosail_opt.csv — run script 12 first:\n  ", opt_csv)
opt_dt       <- data.table::fread(opt_csv)
dopt_by_site <- setNames(
  opt_dt[Norm == "DSM_keepTrees" & d_opt_source == "per_site", d_opt],
  opt_dt[Norm == "DSM_keepTrees" & d_opt_source == "per_site", Site]
)
cli::cli_alert_info(
  "d_opt per site: {paste(names(dopt_by_site), dopt_by_site, sep='=', collapse=', ')}"
)

# ── ATBD parameter samples ─────────────────────────────────────────────────────

cli::cli_alert_info("Sampling ATBD LUT parameters (n = 10 000)...")
atbd_ip <- prosail::get_atbd_lut_input(nb_samples = 10000L,
                                         codistribution_lai = FALSE)
# Column names: lai, lidf_a, lma, brown (all lower-case)
cli::cli_alert_success("ATBD sample columns: {paste(names(atbd_ip), collapse=', ')}")

# ── LiDAR LAI distributions per site ──────────────────────────────────────────

cli::cli_h2("Loading LiDAR LAI rasters...")

lai_full_vals <- list()    # LAIALS: full-depth cumulative PAD
lai_dopt_vals <- list()    # LAIALS_dopt: PAD to site-specific d_opt

for (site in sites) {
  base_ext <- file.path(paths$ext_results, site, "Metrics", "Deciduous_Only")
  lad_path  <- file.path(base_ext, "ladstack_classic.tif")
  max_path  <- file.path(base_ext, "max_res_10_m.tif")

  if (!file.exists(lad_path)) {
    cli::cli_warn("Missing ladstack_classic.tif for {site} — skipping")
    next
  }

  cli::cli_alert_info("{site}: reading ladstack ({lad_path}) ...")
  lad      <- terra::rast(lad_path)           # 38 layers (1 m depth each)
  nlyr_lad <- terra::nlyr(lad)

  # Full-depth LAI (sum all layers — absolute-height bins)
  lai_full <- terra::app(lad, fun = sum, na.rm = TRUE)

  # d_opt-truncated LAI: prefer the script-07 raster (per-pixel canopy-top
  # PAD), then raw PAD file. NEVER sum the bottom of `ladstack_classic.tif`:
  # its layers are absolute-height bins (LAD_Layer_2.5..39.5), not
  # depth-from-canopy-top.
  d_val          <- dopt_by_site[[site]]
  lai_dopt_07    <- file.path(paths$output, "intermediate", "lai_als_dopt",
                              site, "LAI_ALS_dopt_per_site.tif")
  pad_x          <- sprintf("PAD_%.1f_40.tif", 40 - d_val + 0.5)
  pad_raw_path   <- file.path(base_ext, "PAD_Profiles_dsm_keepTrees", pad_x)

  if (file.exists(lai_dopt_07)) {
    lai_dopt <- terra::rast(lai_dopt_07)
  } else if (file.exists(pad_raw_path)) {
    lai_dopt <- terra::rast(pad_raw_path)
  } else {
    cli::cli_warn(
      "LAI_ALS_dopt missing for {site} (d_opt = {d_val}). ",
      "Run scripts/steps/07_compute_lai_als_dopt.R first. Skipping site."
    )
    next
  }

  # Mask to pixels with max height > 10 m
  if (file.exists(max_path)) {
    max_r <- terra::rast(max_path)
    mask  <- max_r > 10
    lai_full <- terra::mask(lai_full, mask, maskvalues = 0)
    lai_dopt <- terra::mask(lai_dopt, mask, maskvalues = 0)
  }

  v_full <- terra::values(lai_full, na.rm = TRUE)[, 1L]
  v_dopt <- terra::values(lai_dopt, na.rm = TRUE)[, 1L]

  # Keep pixels where LAI > 2 (h_min filter equivalent)
  v_full <- v_full[v_full > 2]
  v_dopt <- v_dopt[v_dopt > 0]

  # Subsample if large
  n_max <- 50000L
  if (length(v_full) > n_max) v_full <- sample(v_full, n_max)
  if (length(v_dopt) > n_max) v_dopt <- sample(v_dopt, n_max)

  lai_full_vals[[site]] <- v_full
  lai_dopt_vals[[site]] <- v_dopt

  cli::cli_alert_success(
    "{site}: LAIALS n={length(v_full)}, LAIALS_dopt n={length(v_dopt)} (d_opt={d_val} m)"
  )
}


# ── Helper: build density rows ─────────────────────────────────────────────────

# Returns data.frame with (x, density, distribution, group)
build_dens <- function(x_grid, vals_or_pars, name, group) {
  if (length(vals_or_pars) > 10L) {
    # Empirical kernel density
    d <- density(vals_or_pars, from = min(x_grid), to = max(x_grid),
                 n = length(x_grid))
    data.frame(x = d$x, density = d$y, distribution = name, group = group)
  } else {
    # Truncated normal from named parameter vector
    p <- vals_or_pars
    y <- truncnorm::dtruncnorm(x_grid,
                                a    = p["min"], b    = p["max"],
                                mean = p["mean"], sd  = p["sd"])
    data.frame(x = x_grid, density = y, distribution = name, group = group)
  }
}

# ── Panel (a) — LAI data ───────────────────────────────────────────────────────

x_lai <- seq(0, 15, length.out = 5000L)

lai_theory <- list(
  "ATBD"  = atbd_ip[["lai"]],
  "OPT#1" = c(min = 0, max =  9, mean = 4, sd = 3),
  "OPT#2" = c(min = 0, max = 11, mean = 5, sd = 3)
)

site_labels <- sites

df_lai_theory <- data.table::rbindlist(lapply(names(lai_theory), function(nm)
  build_dens(x_lai, lai_theory[[nm]], nm, "Others")
))

df_lai_empir <- data.table::rbindlist(c(
  lapply(site_labels, function(s)
    build_dens(x_lai, lai_full_vals[[s]], paste0("LAIALS - ", s),      "LAI")),
  lapply(site_labels, function(s)
    build_dens(x_lai, lai_dopt_vals[[s]], paste0("LAIALS_dopt - ", s), "LAI"))
))

df_lai <- rbind(df_lai_theory, df_lai_empir)

# Colour and linetype
lai_levels <- c(
  "ATBD", "OPT#1", "OPT#2",
  paste0("LAIALS - ",      sites),
  paste0("LAIALS_dopt - ", sites)
)

lai_colours <- c(
  "ATBD"                  = "#1f78b4",
  "OPT#1"                 = "#6a3d9a",
  "OPT#2"                 = "#e31a1c",
  "LAIALS - Aigoual"      = "#33a02c",
  "LAIALS - Blois"        = "#B22222",
  "LAIALS - Mormal"       = "#ccaa00",
  "LAIALS_dopt - Aigoual" = "#b2df8a",
  "LAIALS_dopt - Blois"   = "#ff9a99",
  "LAIALS_dopt - Mormal"  = "#fdbf6f"
)

lai_linetypes <- c(
  "ATBD"  = "solid", "OPT#1" = "solid", "OPT#2" = "solid",
  "LAIALS - Aigoual"      = "solid",
  "LAIALS - Blois"        = "solid",
  "LAIALS - Mormal"       = "solid",
  "LAIALS_dopt - Aigoual" = "dashed",
  "LAIALS_dopt - Blois"   = "dashed",
  "LAIALS_dopt - Mormal"  = "dashed"
)

lai_legend <- c(
  "ATBD"   = "ATBD",
  "OPT#1"  = "OPT#1",
  "OPT#2"  = "OPT#2",
  "OPT#3"  = "OPT#3",
  "LAIALS - Aigoual"      = expression(LAI[ALS] ~ "- Aigoual"),
  "LAIALS - Blois"        = expression(LAI[ALS] ~ "- Blois"),
  "LAIALS - Mormal"       = expression(LAI[ALS] ~ "- Mormal"),
  "LAIALS_dopt - Aigoual" = expression(LAI[ALS_dopt] ~ "- Aigoual"),
  "LAIALS_dopt - Blois"   = expression(LAI[ALS_dopt] ~ "- Blois"),
  "LAIALS_dopt - Mormal"  = expression(LAI[ALS_dopt] ~ "- Mormal")
)

df_lai[, distribution := factor(distribution, levels = lai_levels)]
df_lai[, group        := factor(group, levels = c("Others", "LAI"))]

facet_lbl_lai <- ggplot2::as_labeller(
  c("Others" = "Others", "LAI" = "LAI"),
  default = ggplot2::label_value
)

# ── Panel (b) — ALA, BROWN, LMA data ──────────────────────────────────────────

other_dists <- list(
  ALA = list(
    "ATBD"  = atbd_ip[["lidf_a"]],
    "OPT#1" = c(min = 20, max = 50, mean = 35, sd = 20),
    "OPT#2" = c(min = 25, max = 55, mean = 40, sd = 20),
    "OPT#3" = c(min = 30, max = 55, mean = 40, sd = 20),
    "OPT#4" = c(min = 30, max = 60, mean = 45, sd = 20)
  ),
  BROWN = list(
    "ATBD"  = atbd_ip[["brown"]],
    "OPT#1" = c(min = 0, max = 1, mean = 0, sd = 0.1),
    "OPT#2" = NULL                             # point mass at BROWN = 0
  ),
  LMA = list(
    "ATBD"  = atbd_ip[["lma"]],
    "OPT#1" = c(min = 0.003, max = 0.02,  mean = 0.010, sd = 0.005),
    "OPT#2" = c(min = 0.003, max = 0.030, mean = 0.015, sd = 0.005)
  )
)

x_ranges <- list(
  ALA   = seq(  0,   90, length.out = 5000L),
  BROWN = seq(  0,  1.5, length.out = 5000L),
  LMA   = seq(0.0, 0.035, length.out = 5000L)
)

df_other_list <- lapply(names(other_dists), function(param) {
  x <- x_ranges[[param]]
  rows <- lapply(names(other_dists[[param]]), function(nm) {
    v <- other_dists[[param]][[nm]]
    if (is.null(v)) {
      # BROWN OPT#2: vertical spike at 0 (BROWN always 0)
      data.frame(x = c(0, 0), density = c(0, 1),
                 distribution = nm, group = param)
    } else {
      build_dens(x, v, nm, param)
    }
  })
  data.table::rbindlist(rows, fill = TRUE)
})

df_other <- data.table::rbindlist(df_other_list, fill = TRUE)

other_levels <- c("ATBD", "OPT#1", "OPT#2", "OPT#3", "OPT#4")
other_colours <- c(
  "ATBD"  = "#1f78b4",
  "OPT#1" = "#6a3d9a",
  "OPT#2" = "#e31a1c",
  "OPT#3" = "#ff7f00",
  "OPT#4" = "#654321"
)

df_other[, distribution := factor(distribution, levels = other_levels)]
df_other[, group        := factor(group, levels = c("ALA", "BROWN", "LMA"))]

# ── Shared theme ───────────────────────────────────────────────────────────────

theme_pub <- ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    legend.position   = "bottom",
    legend.title      = ggplot2::element_blank(),
    legend.text       = ggplot2::element_text(size = 9),
    legend.key.width  = ggplot2::unit(1.4, "cm"),
    strip.text        = ggplot2::element_text(face = "bold", size = 10),
    strip.background  = ggplot2::element_rect(fill = "grey92", colour = NA),
    axis.title        = ggplot2::element_text(size = 11),
    panel.grid.minor  = ggplot2::element_blank()
  )

# ── Panel (a) ─────────────────────────────────────────────────────────────────

p_lai <- ggplot2::ggplot(
  df_lai,
  ggplot2::aes(x = x, y = density,
               colour   = distribution,
               linetype = distribution)
) +
  ggplot2::geom_line(linewidth = 1.1, na.rm = TRUE) +
  ggplot2::facet_wrap(~ group, scales = "free_y", ncol = 1,
                       labeller = facet_lbl_lai) +
  ggplot2::scale_colour_manual(values = lai_colours, labels = lai_legend) +
  ggplot2::scale_linetype_manual(values = lai_linetypes, guide = "none") +
  ggplot2::labs(x = "LAI", y = "Density") +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2)) +
  theme_pub

# ── Panel (b) ─────────────────────────────────────────────────────────────────

p_other <- ggplot2::ggplot(
  df_other,
  ggplot2::aes(x = x, y = density, colour = distribution)
) +
  ggplot2::geom_line(linewidth = 1.1, na.rm = TRUE) +
  ggplot2::facet_wrap(~ group, scales = "free", ncol = 1) +
  ggplot2::scale_colour_manual(values = other_colours) +
  ggplot2::labs(x = "Value", y = "Density") +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1)) +
  theme_pub

# ── Combine ────────────────────────────────────────────────────────────────────

fig <- (p_lai | p_other) +
  patchwork::plot_annotation(tag_levels = list(c("(a)", "(b)"))) &
  ggplot2::theme(
    plot.tag          = ggplot2::element_text(size = 13, face = "bold"),
    plot.tag.position = c(0.02, 0.99)
  )

# ── Save ───────────────────────────────────────────────────────────────────────

out_pdf <- file.path(out_dir, "sm_prosail_param_distribs.pdf")
out_png <- file.path(out_dir, "sm_prosail_param_distribs.png")

ggplot2::ggsave(out_pdf, fig, width = 30, height = 22, units = "cm", device = cairo_pdf)
ggplot2::ggsave(out_png, fig, width = 30, height = 22, units = "cm",
                device = "png", dpi = 600, bg = "white")

cli::cli_alert_success("Written: {out_pdf}")
cli::cli_alert_success("Written: {out_png}")
cli::cli_h1("Done — PROSAIL parameter distribution figures")
