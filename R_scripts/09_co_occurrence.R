# ================================================================
# Probabilistic co-ocurrence
# Figure - Heatmap of co-ocurrence
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

load_pkg("microbiome")
load_pkg("CooccurrenceAffinity")

# ================================================================
# -------------------- LOAD PHYLOSEQ -----------------------------
# ================================================================

ps <- readRDS(file.path(DIR_OUT, "ps_counts.rds")) ########## Import non-normalized dataset

# ================================================================
# -------------------- COLLAPSE TO GENUS -------------------------
# ================================================================

# Collapse

genus <- tax_glom(ps, "Genus")

# Check

if (any(taxa_sums(genus) == 0)) stop("Zero-count taxa detected.")
if (any(sample_sums(genus) == 0)) stop("Zero-count samples detected.")

# ================================================================
# -------------------- EXAMINE CO-OCCURRENCE ---------------------
# ================================================================

# ---------------------------
# 1) Parameters
# ---------------------------

MIN_PREV_PROP <- 0.10
MAX_PREV_PROP <- 1.00

# ---------------------------
# 2) Presence/absence matrix
# ---------------------------

mat <- as(otu_table(genus), "matrix")
if (taxa_are_rows(genus)) mat <- t(mat)

pa <- (mat > 0) * 1L

tax <- as.data.frame(tax_table(genus))
gen <- tax[colnames(mat), "Genus", drop = TRUE]
colnames(pa) <- gen

# ---------------------------
# 3) Prevalence filter
# ---------------------------

prev_prop <- colSums(pa) / nrow(pa)
keep <- prev_prop >= MIN_PREV_PROP & prev_prop <= MAX_PREV_PROP

pa <- pa[, keep, drop = FALSE]

if (ncol(pa) < 2) stop("Not enough taxa after filtering.")

# ---------------------------
# 4) Run CooccurrenceAffinity
# ---------------------------

res <- affinity(
  data       = pa,
  row.or.col = "col",
  datatype   = "binary"
)

print(head(res$all))

# ---------------------------
# 5) Examine results
# ---------------------------

tbl_all <- res$all %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  select(entity_1, entity_2, alpha_mle, ci_blaker, p_adj) %>%
  as.data.frame()

tbl_filt <- tbl_all %>%
  dplyr::filter(p_adj < FDR) %>%
  dplyr::arrange(entity_1, entity_2) %>%
  dplyr::mutate(
    alpha_mle = formatC(alpha_mle, format = "f", digits = 2),
    p_adj     = formatC(p_adj,     format = "f", digits = 2),
    ci_lo = readr::parse_number(ci_blaker),
    ci_hi = readr::parse_number(sub(".*,\\s*", "", sub("\\]$", "", sub("^\\[", "", ci_blaker)))),
    ci_blaker = paste0(
      formatC(ci_lo, format = "f", digits = 2),
      ", ",
      formatC(ci_hi, format = "f", digits = 2)
    )
  ) %>%
  dplyr::select(-ci_lo, -ci_hi) %>%
  dplyr::rename(
    "Genus 1" = entity_1,
    "Genus 2" = entity_2,
    "α̂" = alpha_mle,
    "95% confidence interval" = ci_blaker,
    "q-value" = p_adj
  )

print(head(tbl_filt))

# ================================================================
# -------------------- CREATING DATASETS FOR FIGURES -------------
# ================================================================

saveRDS(
  res,
  file.path(DIR_OUT, "heatmap_cooccurrence_fig_df.rds")
)

# ================================================================
# -------------------- EXPORT ------------------------------------
# ================================================================

write.csv(tbl_filt, file.path(DIR_TBL, "Table_cooccurrence.csv"), row.names = FALSE)
