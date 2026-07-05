# ================================================================
# Figure - Krona chart
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

# ================================================================
# -------------------- LOAD PHYLOSEQ -----------------------------
# ================================================================

ps <- readRDS(file.path(DIR_OUT, "ps_counts.rds")) ########## Import non-normalized dataset

# ================================================================
# -------------------- COLLAPSE TO GENUS -------------------------
# ================================================================

# Collapse

genus <- tax_glom(ps, taxrank = "Genus")

# Check

if (any(taxa_sums(genus) == 0)) stop("Zero-count taxa detected.")
if (any(sample_sums(genus) == 0)) stop("Zero-count samples detected.")

# ================================================================
# -------------------- DATA WRANGLING ----------------------------
# ================================================================

# Extract counts

counts <- as(otu_table(genus), "matrix")
if (taxa_are_rows(genus)) counts <- t(counts)

# Total counts per genus group

genus_counts <- colSums(counts)

# Get taxonomy from the tax_glom output

tax <- as.data.frame(tax_table(genus))

# Add genus_counts to taxonomy table using taxon IDs

tax$taxon_id <- rownames(tax)

df <- tax %>%
  dplyr::mutate(count = genus_counts[taxon_id]) %>%
  dplyr::select(count, Kingdom, Phylum, Class, Order, Family, Genus)

krona_df <- df %>% dplyr::arrange(dplyr::desc(count))

# ================================================================
# -------------------- EXPORT ------------------------------------
# ================================================================

write.table(
  krona_df,
  file.path(DIR_OUT, "krona_fig_df.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# ================================================================
# -------------------- TERMINAL ----------------------------------
# ================================================================

# conda install -c bioconda krona
# ktImportText -o ~/Desktop/ITS/Figures/Fig_krona_chart.html \
# ~/Desktop/ITS/Output/krona_fig_df.txt
