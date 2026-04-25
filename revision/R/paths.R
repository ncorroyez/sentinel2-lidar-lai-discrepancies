# ---
# title:  paths.R
# desc:   Reads config.yml at the project root and exports a `paths` list
#         with absolute paths to all external data directories. Source this
#         file at the top of any script that reads from 01_DATA, 03_RESULTS,
#         or PROSAIL-Optimization.
#
#         Usage:
#           source(here::here("revision", "R", "paths.R"))
#           # then use paths$ext_results, paths$raw_data, etc.
#
#         Fallback: if config.yml is absent, data_root defaults to the
#         project root (here::here()), preserving backward compatibility
#         for local setups where 01_DATA/ and 03_RESULTS/ live in the
#         project directory.
# ---

local({
  cfg_path <- here::here("config.yml")

  if (file.exists(cfg_path)) {
    if (!requireNamespace("yaml", quietly = TRUE))
      stop("Package 'yaml' is required to read config.yml.\n",
           "Install with: install.packages('yaml')")
    cfg       <- yaml::read_yaml(cfg_path)
    data_root <- cfg$data_root
    dirs      <- cfg$dirs
  } else {
    # No config.yml — assume data lives at the project root (default layout).
    data_root <- here::here()
    dirs <- list(
      raw_data      = "01_DATA",
      ext_results   = "03_RESULTS",
      prosail_lidar = "PROSAIL-Optimization/01_DATA",
      prosail_codes = "PROSAIL-Optimization/02_CODES"
    )
  }

  if (!dir.exists(data_root))
    stop(
      "data_root does not exist: ", data_root, "\n",
      "Check config.yml (or copy config.yml.template to config.yml ",
      "and set data_root to your local SMB mount path)."
    )

  paths <<- list(
    data_root     = data_root,
    raw_data      = file.path(data_root, dirs$raw_data),
    ext_results   = file.path(data_root, dirs$ext_results),
    prosail_lidar = file.path(data_root, dirs$prosail_lidar),
    prosail_codes = file.path(data_root, dirs$prosail_codes)
  )

  invisible(paths)
})
