# ================================================================
# Core microbiome
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
# -------------------- CORE MICROBIOME ---------------------------
# ================================================================

genus <- microbiome::transform(genus, "compositional")

unique(rowSums(as.data.frame(otu_table(genus)))) # This should be 1

# Definition 1

core25 <- core_members(genus, detection = 0, prevalence = 25/100) # >0 in ≥ 25% samples (only using prevalence
# based on data characteristics)

print(core25)

tax_table(genus)[c("ASV4", "ASV5"), ]

# Definition 2

core50 <- core_members(genus, detection = 0, prevalence = 50/100) # >0 in ≥ 50% samples (only using prevalence
# based on data characteristics)

print(core50) # No genera identified

# Examining prevalences

prev <- prevalence(genus)

# Convert to %

prev_percent <- prev*100

# View taxa prevalence

prev_percent[c("ASV4", "ASV5")]
