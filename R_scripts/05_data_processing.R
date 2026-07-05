# ================================================================
# Import all ITS seq data and build a phyloseq object
# Export an index table of ASV full sequences and their ASV ID numbers
# Fix taxonomy table
# Data filtering
# Create some new variables
# Create rarefied dataset
# ================================================================

# ================================================================
# -------------------- CLEAN ENVIRONMENT -------------------------
# ================================================================

rm(list = ls())
gc()

# ================================================================
# --------------------SET WORKING DIRECTORY ----------------------
# ================================================================

setwd("~/Desktop/ITS")

# ================================================================
# -------------------- LOAD SOURCES AND PACKAGES -----------------
# ================================================================

source("./R_scripts/02_global_configuration.R")
source("./R_scripts/03_global_packages.R") # phyloseq, tidyverse, vegan, Hmisc, svglite

load_pkg("qiime2R")
load_pkg("lubridate")

# ================================================================
# -------------------- LOAD AND CLEAN TAXONOMY -------------------
# ================================================================

taxonomy <- read.csv(file.path(DIR_DATA, "taxonomy.csv"),
                           na.strings = c("", "NA"),
                           stringsAsFactors = FALSE) %>%
  rename(Feature.ID = X.OTUID,
         Taxon = taxonomy)

taxonomy_parsed <- parse_taxonomy(taxonomy) %>%
  as.data.frame()

taxonomy_parsed %>%
  dplyr::summarize(across(everything(), ~ mean(is.na(.) | . == ""))) # Proportion of missing taxonomy by rank

fill_unclassified <- function(child, parent) {
  if_else(
    child == "",
    if_else(
      grepl("_unclassified$", parent),
      parent,
      paste0(parent, "_unclassified")
    ),
    child
  )
}

taxonomy_clean <- taxonomy_parsed %>%
  mutate(across(1:7, as.character)) %>%
  mutate(
    across(1:7, ~na_if(.x, "unidentified")),
    across(1:7, ~ifelse(.x == "Unassigned", "Other", .x))
  ) %>%
  replace(is.na(.), "") %>%
  mutate(
    Phylum  = fill_unclassified(Phylum,  Kingdom),
    Class   = fill_unclassified(Class,   Phylum),
    Order   = fill_unclassified(Order,   Class),
    Family  = fill_unclassified(Family,  Order),
    Genus   = fill_unclassified(Genus,   Family),
    Species = fill_unclassified(Species, Genus)
  )

# ================================================================
# -------------------- LOAD METADATA -----------------------------
# ================================================================

metadata_raw <- read.csv(file.path(DIR_DATA, "metadata.csv"),
                         row.names = 1,
                         na.strings = c("", "NA"),
                         stringsAsFactors = FALSE)

# ================================================================
# -------------------- LOAD AND CLEAN OTU TABLE ------------------
# ================================================================

otu <- read.csv(file.path(DIR_DATA, "asvtable.csv"),
                          row.names = 1,
                          check.names = FALSE,
                          na.strings = c("", "NA"),
                          stringsAsFactors = FALSE)

# Fix incorrect ID (confirmed with lab)

rownames(otu) <- ifelse(
  rownames(otu) == "TH10_338_r2",
  "TH10_338",
  rownames(otu)
)

# Ensure samples are rows

if (!all(rownames(otu) %in% rownames(metadata_raw))) {
  otu <- t(otu)
}

if (!all(rownames(otu) %in% rownames(metadata_raw))) {
  stop("Sample names do not match between OTU table and metadata")
}

# Reordering metadata

metadata_sub <- metadata_raw[rownames(otu), ]

# ================================================================
# -------------------- BUILD PHYLOSEQ OBJECT ---------------------
# ================================================================

if (!identical(rownames(otu), rownames(metadata_sub))) {
  stop("Sample order mismatch between OTU table and metadata")
}

OTU <- otu_table(as.matrix(otu), taxa_are_rows = FALSE)
TAX <- tax_table(as.matrix(taxonomy_clean))
SAM <- sample_data(metadata_sub)

ps <- phyloseq(OTU, TAX, SAM)

# ================================================================
# -------------------- ADD INITIAL TOTAL NUMBER OF READS ---------
# ================================================================

sample_data(ps)$total_reads <- sample_sums(ps)
sample_data(ps)$log10_total_reads <- log10(sample_data(ps)$total_reads)

# ================================================================
# -------------------- CHECKING FOR SINGLETONS -------------------
# ================================================================

# Extract otu table

otu_mat <- as(otu_table(ps), "matrix")

# Examine if singleton exists anywhere

sum(otu_mat == 1) # Singletons were filtered out during DADA2

# ================================================================
# -------------------- RENAME ASVs -------------------------------
# ================================================================

asv_table <- data.frame(QIIME2_Name = taxa_names(ps)) %>%
  mutate(ASV_Name = paste0("ASV", row_number())) # Rename to facilitate visualization

taxa_names(ps) <- asv_table$ASV_Name

# ================================================================
# -------------------- ADD ASVs NUMBER TO SPECIES ----------------
# ================================================================

tax_table(ps)[, "Species"] <- paste0(
  tax_table(ps)[, "Species"],
  "_",
  taxa_names(ps)
)

# ================================================================
# -------------------- FILTER ------------------------------------
# ================================================================

# Keep only sequences classified as fungi down to phylum level

tax_ranks <- c("Kingdom", "Phylum", "Class", "Order")

lapply(tax_ranks, function(rank) {
  cat("\n", rank, "\n")
  print(table(tax_table(ps)[, rank]))
})

ps <- subset_taxa(ps, Kingdom=="Fungi" & Phylum!="Fungi_unclassified")

lapply(tax_ranks, function(rank) {
  cat("\n", rank, "\n")
  print(table(tax_table(ps)[, rank]))
})

# Keep only samples with library sizes >500 

ps <- prune_samples(sample_sums(ps)>500, ps) 

# Keep only taxa that appear in at least one sample

ps <- prune_taxa(taxa_sums(ps)>0, ps)

# Check

if (any(taxa_sums(ps) == 0)) stop("Zero-count taxa detected.")
if (any(sample_sums(ps) == 0)) stop("Zero-count samples detected.")

# ================================================================
# -------------------- ADD SCALED AGE VARIABLE -------------------
# ================================================================

sample_data(ps)$scaled_age_sample <-
  as.numeric(scale(as.numeric(sample_data(ps)$age_sample)))

print(head(as.data.frame(sample_data(ps))))

# ================================================================
# -------------------- RAREFYING TO LOWEST LIBRARY ---------------
# ================================================================

# Extract OTU table

count <- as.data.frame(otu_table(ps))

# Rarefaction depth

depth_min <- min(rowSums(count))

# Repeat rarefaction 400 times

set.seed(123)
n_rarefy <- 400
rarefied_list_min <- replicate(n_rarefy, rrarefy(count, depth_min), simplify = FALSE)

# Average counts across rarefactions

rarefied_array_min <- simplify2array(rarefied_list_min) # Creates array
averaged_min <- apply(rarefied_array_min, c(1, 2), mean) # Averages across iterations

# Smart rounding function

smart.round <- function(x) { 
  y <- floor(x)
  idx <- tail(order(x - y), round(sum(x)) - sum(y))
  y[idx] <- y[idx] + 1
  y
}

# Apply smart rounding row-wise

averaged_r_min <- t(apply(averaged_min, 1, smart.round))

# Convert to data frame (optional)

averaged_r_min <- as.data.frame(averaged_r_min)

Hmisc::describe(rowSums(averaged_r_min)) # it works!

# Create new phyloseq

rar_min <- ps

otu_table(rar_min) <- otu_table(
  as.matrix(averaged_r_min),
  taxa_are_rows = FALSE
)

# Keep only taxa that appear in at least one sample

rar_min <- prune_taxa(taxa_sums(rar_min)>0, rar_min)

# Check

if (any(taxa_sums(rar_min) == 0)) stop("Zero-count taxa detected.")
if (any(sample_sums(rar_min) == 0)) stop("Zero-count samples detected.")

# ================================================================
# -------------------- RAREFYING TO SPECIFIC DEPTH ---------------
# ================================================================

# Rarefaction depth

depth_sp <- 1000 # Per reviewers' comments

# Filtering by specific depth

filt <- count[rowSums(count) >= depth_sp, ]

# Repeat rarefaction 400 times

set.seed(123)
rarefied_list_sp <- replicate(n_rarefy, rrarefy(filt, depth_sp), simplify = FALSE)

# Average counts across rarefactions

rarefied_array_sp <- simplify2array(rarefied_list_sp) # creates array
averaged_sp <- apply(rarefied_array_sp, c(1, 2), mean) # averages across iterations

# Apply smart rounding row-wise

averaged_r_sp <- t(apply(averaged_sp, 1, smart.round))

# Convert to data frame (optional)

averaged_r_sp <- as.data.frame(averaged_r_sp)

Hmisc::describe(rowSums(averaged_r_sp))

# Create new phyloseq

rar_sp <- ps

otu_table(rar_sp) <- otu_table(
  as.matrix(averaged_r_sp),
  taxa_are_rows = FALSE
)

# Keep only taxa that appear in at least one sample

rar_sp <- prune_taxa(taxa_sums(rar_sp)>0, rar_sp)

# Check

if (any(taxa_sums(rar_sp) == 0)) stop("Zero-count taxa detected.")
if (any(sample_sums(rar_sp) == 0)) stop("Zero-count samples detected.")

# ================================================================
# -------------------- EXPORT ------------------------------------
# ================================================================

write.csv(asv_table, file.path(DIR_OUT, "list_of_full_sequences_and_asv_numbers.csv"), row.names = FALSE)
saveRDS(ps, file.path(DIR_OUT, "ps_counts.rds"))
saveRDS(rar_min, file.path(DIR_OUT, "ps_rar_min.rds"))
saveRDS(rar_sp, file.path(DIR_OUT, "ps_rar_sp.rds"))
