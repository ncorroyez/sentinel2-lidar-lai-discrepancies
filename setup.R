# ---
# title:  setup.R
# desc:   One-shot dependency installer for the LiDAR x Sentinel-2 LAI
#         pipeline. Installs all CRAN packages, the three jbferet GitHub
#         packages (prosail / prospect / preprocS2) and the archived
#         liquidSVM package, then verifies that every package loads.
#
#         Run once after cloning the repository:
#           source("setup.R")
#
#         The script is idempotent: packages already installed are left
#         alone. To force reinstall, set force = TRUE below.
# ---

force <- FALSE   # set to TRUE to reinstall packages already present

# ── CRAN packages ─────────────────────────────────────────────────────────────

cran_pkgs <- c(
  # Utilities
  "here", "yaml", "cli", "stringr", "remotes",
  # Geospatial
  "terra", "sf", "stars", "sgsR", "lidR",
  # Data
  "data.table", "readr", "dplyr",
  # Plotting
  "ggplot2", "patchwork", "scales",
  # Parallel
  "future", "future.apply", "progressr",
  # Stats / PROSAIL training
  "truncnorm"
)

installed <- rownames(utils::installed.packages())
missing   <- if (force) cran_pkgs else setdiff(cran_pkgs, installed)

if (length(missing)) {
  message("Installing ", length(missing), " CRAN package(s): ",
          paste(missing, collapse = ", "))
  utils::install.packages(missing)
} else {
  message("All CRAN packages already installed.")
}

# ── liquidSVM (archived from CRAN in 2020) ────────────────────────────────────
#
# liquidSVM is the SVR backend used by prosail::PROSAIL_Hybrid_Train. The
# package was archived from CRAN; install from the CRAN archive tarball.

if (force || !requireNamespace("liquidSVM", quietly = TRUE)) {
  message("Installing liquidSVM from CRAN archive...")
  liquidsvm_url <- paste0(
    "https://cran.r-project.org/src/contrib/Archive/liquidSVM/",
    "liquidSVM_1.2.4.tar.gz"
  )
  utils::install.packages(liquidsvm_url, repos = NULL, type = "source")
} else {
  message("liquidSVM already installed.")
}

# ── GitHub packages (jbferet) ─────────────────────────────────────────────────

gh_pkgs <- c(
  prospect  = "jbferet/prospect",
  prosail   = "jbferet/prosail",
  preprocS2 = "jbferet/preprocS2"
)

for (pkg in names(gh_pkgs)) {
  if (force || !requireNamespace(pkg, quietly = TRUE)) {
    message("Installing ", pkg, " from github::", gh_pkgs[[pkg]], " ...")
    remotes::install_github(gh_pkgs[[pkg]], upgrade = "never")
  } else {
    message(pkg, " already installed (", utils::packageVersion(pkg), ")")
  }
}

# ── Verification ──────────────────────────────────────────────────────────────

all_pkgs <- c(cran_pkgs, "liquidSVM", names(gh_pkgs))
ok <- vapply(all_pkgs, function(p) requireNamespace(p, quietly = TRUE),
             logical(1))

message("\n── Setup summary ──")
for (p in all_pkgs) {
  status <- if (ok[[p]]) "OK " else "FAIL"
  ver    <- if (ok[[p]]) as.character(utils::packageVersion(p)) else "-"
  message(sprintf("  [%s] %-15s %s", status, p, ver))
}

if (!all(ok)) {
  stop("Some packages failed to install. See messages above.")
}

# Minimum prosail version check — the pipeline relies on generate_lut_4sail()
# and spec_soil_atbd_v2, both introduced in prosail 3.0.0.
if (utils::packageVersion("prosail") < "3.0.0") {
  warning("prosail < 3.0.0 detected. The pipeline requires >= 3.0.0. ",
          "Re-run setup.R with force = TRUE.")
}

message("\nAll dependencies installed. You can now run the pipeline:")
message("  source(\"scripts/00_run_all.R\")")
