# ================================================================
# Include global core packages
# ================================================================

# ================================================================
# -------------------- GLOBAL CORE PACKAGES ----------------------
# ================================================================

core_packages <- c(
  "phyloseq",
  "tidyverse",
  "vegan",
  "Hmisc",
  "svglite"
)

installed <- core_packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(core_packages[!installed])

lapply(core_packages, library, character.only = TRUE)

# ================================================================
# -------------------- LOADING FUNCTION --------------------------
# ================================================================

load_pkg <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
