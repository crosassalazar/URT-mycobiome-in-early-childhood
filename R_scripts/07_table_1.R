# ================================================================
# Create Table 1
# ================================================================

# ================================================================
# -------------------- CLEAN ENVIRONMENT -------------------------
# ================================================================

rm(list = ls())
gc()

# ================================================================
# -------------------- SET WORKING DIRECTORY ---------------------
# ================================================================

setwd("~/Desktop/ITS")

# ================================================================
# -------------------- LOAD SOURCES AND PACKAGES -----------------
# ================================================================

source("./R_scripts/02_global_configuration.R")
source("./R_scripts/03_global_packages.R") # phyloseq, tidyverse, vegan, Hmisc, svglite
source("./R_scripts/04_global_labels.R")

load_pkg("arsenal")

# ================================================================
# -------------------- LOAD PHYLOSEQ -----------------------------
# ================================================================

ps <- readRDS(file.path(DIR_OUT, "ps_counts.rds"))

# ================================================================
# -------------------- EXTRACT METADATA AND RECODE VARIABLES -----
# ================================================================

metadata <- data.frame(sample_data(ps))

metadata <- metadata %>%
  mutate(
    efdelivery = ifelse(efdelivery == 2, 0, 1)
  )

# ================================================================
# -------------------- FACTOR VARIABLES --------------------------
# ================================================================

cols <- c(
  "efsex", "race_white_dev", "efdelivery", "any_breastfd_dev",
  "mat_asthma_dev", "young_sib_dev", "early_life_anti_dev",
  "early_life_smoke_dev", "insurance_dev", "petany_dev",
  "efrural_dev", "any_lrti_dev", "mold_dev", "season_collection"
)

# Keep only variables that exist

vars_present <- intersect(cols, names(metadata))

# Convert to factor

metadata[vars_present] <- lapply(metadata[vars_present], factor)

# ================================================================
# -------------------- LABEL VARIABLES ---------------------------
# ================================================================

metadata <- apply_labels(metadata)

# ================================================================
# -------------------- TABLE 1 SETTINGS --------------------------
# ================================================================

tbl_control <- tableby.control(
  test = TRUE,
  total = TRUE,
  numeric.test = "kwt",
  cat.test = "fe",
  numeric.stats = c("medianq1q3"),
  cat.stats = c("countpct"),
  cat.simplify = TRUE,
  stats.labels = list(
    medianq1q3 = "Median (25th, 75th percentile)"
  )
)

# ================================================================
# -------------------- CREATE TABLE ------------------------------
# ================================================================

tbl_df <- tableby(
  ~ enrollageinweeks_dev +
    efgestagewk +
    birthWgtGm_dev +
    age_sample +
    efsex +
    race_white_dev +
    efdelivery +
    any_breastfd_dev +
    mat_asthma_dev +
    young_sib_dev +
    early_life_anti_dev +
    early_life_smoke_dev +
    insurance_dev +
    petany_dev +
    efrural_dev +
    any_lrti_dev +
    mold_dev +
    season_collection,
  data = metadata,
  control = tbl_control
)

tbl_1 <- summary(
  tbl_df,
  title = "Table 1. Baseline characteristics of children included in the study (n=148).",
  text = TRUE,
  digits = 2,
  digits.pct = 2,
  digits.p = 2,
  cat.simplify = TRUE
)

tbl_1

# ================================================================
# -------------------- EXPORT ------------------------------------
# ================================================================

write2word(tbl_1, file.path(DIR_TBL, "Table_1.doc"))

