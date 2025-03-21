# List of packages available on CRAN
cran_packages <- c("igraph", "miceadds", "MASS", "dti", "Rfast", "gtools", 
                   "pracma", "nimble", "ar.matrix", "rlist", 
                   "doParallel", "foreach", "ramcmc")

# Install and load CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

# Handle ts.extend separately from CRAN Archive
if (!requireNamespace("ts.extend", quietly = TRUE)) {
  message("Installing ts.extend from CRAN Archive...")
  ts_extend_url <- "https://cran.r-project.org/src/contrib/Archive/ts.extend/ts.extend_0.1.0.tar.gz"
  ts_extend_tar <- tempfile(fileext = ".tar.gz")
  download.file(ts_extend_url, ts_extend_tar, mode = "wb")
  install.packages(ts_extend_tar, repos = NULL, type = "source")
}
library(ts.extend)