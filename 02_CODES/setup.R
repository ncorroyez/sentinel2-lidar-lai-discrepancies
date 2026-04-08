# ---
# title: "setup.R"
# aim: Directories architecture & install dependencies
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, 
# CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-05-24"
# ---

# ON-GOING, NOT FINISHED

# Function to install CRAN packages
install_if_missing <- function(packages) {
  installed_packages <- rownames(installed.packages())
  for (pkg in packages) {
    if (pkg == "stars") {
      if (!("stars" %in% installed_packages) 
          || packageVersion("stars") != "0.6.4") {
        devtools::install_version("stars", 
                                  version = "0.6.4", 
                                  repos = "http://cran.us.r-project.org")
      } else {
        cat("'stars' package version 0.6.4 is already installed\n")
      }
    } else {
      if (!pkg %in% installed_packages) {
        install.packages(pkg)
      } else {
        cat(paste("Package", pkg, "is already installed.\n"))
      }
    }
  }
}

# Install devtools packages
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Read dependencies.txt and extract packages names
installed_pkgs <- read.table("dependencies.txt", header = FALSE, stringsAsFactors = FALSE)[,1]

# Install packages
install_if_missing(cran_packages)

# Install GitHub packages
for (pkg in github_packages) {
  devtools::install_github(pkg)
}

# Directories architecture
create_dir_if_missing <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    cat(paste("Created directory:", dir, "\n"))
  } else {
    cat(paste("Already existing directory:", dir, "\n"))
  }
}

create_dir_if_missing("01_DATA", showWarnings = FALSE)
create_dir_if_missing("02_CODES", showWarnings = FALSE)
create_dir_if_missing("03_RESULTS", showWarnings = FALSE)

cat("**Setup configuration: SUCCESS**.\n")
