# ================================================================
# Differential prevalence analyses
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
# -------------------- PREVALENCE FILTERING ----------------------
# ================================================================

# Keep onlu taxa above minimum prevalence

keep <- taxa_names(genus)[colMeans(otu_table(genus) > 0) >= MIN_PREV]

filt <- prune_taxa(keep, genus)

# Keep only taxa that appear in at least one sample

filt <- prune_taxa(taxa_sums(filt)>0, filt)

# Keep only samples with reads

filt <- prune_samples(sample_sums(filt) > 0, filt)

# Check

if (any(taxa_sums(filt) == 0)) stop("Zero-count taxa detected.")
if (any(sample_sums(filt) == 0)) stop("Zero-count samples detected.")

# ================================================================
# -------------------- EXTRACT DATASETS --------------------------
# ================================================================

otu  <- data.frame(otu_table(filt))

metadata <- data.frame(sample_data(filt))

tax  <- data.frame(tax_table(filt))

tax$ASV <- rownames(tax)

tax$Genus <- as.character(tax$Genus)

# ================================================================
# -------------------- DATA WRANGLING ----------------------------
# ================================================================

# Rename OTU columns to genus names

colnames(otu) <- tax$Genus[match(colnames(otu), tax$ASV)]

# Factor conversion and relevel

factor_vars <- c(
  "efsex","race_white_dev","efdelivery",
  "any_breastfd_dev","mat_asthma_dev","young_sib_dev",
  "early_life_anti_dev","early_life_smoke_dev",
  "insurance_dev","petany_dev","efrural_dev",
  "any_lrti_dev","season_collection","mold_dev"
)

metadata[factor_vars] <- lapply(metadata[factor_vars], factor)

metadata$season_collection <- relevel(metadata$season_collection, ref = "Winter")

metadata <- droplevels(metadata)

# ================================================================
# -------------------- MAASLIN3 GROUP-WISE DIFFERENCES -----------
# ================================================================

load_pkg("maaslin3")

# Set up contrasts

options(contrasts = c("contr.treatment", "contr.poly"))

print(getOption("contrasts"))

# Parameters

DIS_MAS <- file.path(DIR_OUT, "MaAslin3")

if (!dir.exists(DIS_MAS)) {
  dir.create(DIS_MAS, recursive = TRUE)
}

# Predictors of interest

preds <- c(
  "season_collection",
  "early_life_anti_dev",
  "age_sample",
  "efsex",
  "race_white_dev",
  "log10_total_reads"
)

# Check reference

stopifnot(levels(metadata$season_collection)[1] == "Winter")

# Mulitvariate model

set.seed(SEED)

m_group <- maaslin3(
  input_data     = otu,
  input_metadata = metadata,
  output         = file.path(DIS_MAS, "Main"),
  formula        = "~ group(season_collection) + early_life_anti_dev + age_sample + efsex + race_white_dev + log10_total_reads",
  max_significance = FDR,
  standardize      = TRUE,
  max_pngs         = 250
)

# Helper function to extract prevalence results and recalculate q-value as per MaAsLin3 instructions

extract_prev_full <- function(fit, tax_table) {
  fit$fit_data_prevalence$results %>%
    rename(ASV = feature) %>%
    group_by(metadata) %>%
    mutate(
      qval_correct = p.adjust(pval_individual, method = "BH")
    ) %>%
    ungroup() %>%
    left_join(tax_table, by = "ASV") %>%
    arrange(metadata, pval_individual)
}

# Extract all prevalence results

results_full_prev <- extract_prev_full(m_group, tax)

# Explore only significant results

results_full_prev %>%
  select("ASV", "metadata", "coef", "stderr", "N", "N_not_zero", "pval_individual", "qval_correct") %>%
  filter(qval_correct < FDR) %>%
  print()

# Create table of results for season of collection

results_season_prev <- results_full_prev %>%
  filter(metadata=="season_collection") %>%
  mutate(
    p_value = vapply(pval_individual, format_p, character(1)),
    q_value = vapply(qval_correct, format_p, character(1))
    ) %>%
  select(ASV, N, N_not_zero, p_value, q_value) %>%
  arrange(ASV) %>%
  rename(
    Genus = ASV,
    "Prevalence count" = N_not_zero,
    "p-value" = p_value,
    "q-value" = q_value
    )

# ================================================================
# -------------------- MAASLIN3 CONTRASTS ------------------------
# ================================================================

# Multvariate model

set.seed(SEED)

m_contrast <- maaslin3(
  input_data     = otu,
  input_metadata = metadata,
  output         = file.path(DIS_MAS, "Contrasts"),
  formula        = "~ season_collection + early_life_anti_dev + age_sample + efsex + race_white_dev + log10_total_reads",
  max_significance = FDR,
  standardize      = TRUE,
  max_pngs         = 250
)

# Identify the exact coefficient names for season levels in the prevalence results

season_coef_names <- m_contrast$fit_data_prevalence$results %>%
  filter(metadata == "season_collection") %>%
  distinct(name) %>%
  pull(name)

season_coef_names

levels(metadata$season_collection)

level_map <- c(
  Spring = grep("Spring", season_coef_names, value = TRUE),
  Summer = grep("Summer", season_coef_names, value = TRUE),
  Fall   = grep("Fall", season_coef_names, value = TRUE)
)

level_map

# Create a 3 x 3 contrast matrix for Summer–Fall, Summer–Spring, Fall–Spring

contrast_mat <- matrix(
  0, nrow = 3, 
  ncol = length(season_coef_names),
  dimnames = list(
    c("Summer_vs_Fall","Summer_vs_Spring","Fall_vs_Spring"),
    season_coef_names
  )
)

contrast_mat["Summer_vs_Fall", 
             c(level_map["Summer"], level_map["Fall"])] <- c(1, -1)

contrast_mat["Summer_vs_Spring", 
             c(level_map["Summer"], level_map["Spring"])] <- c(1, -1)

contrast_mat["Fall_vs_Spring", 
             c(level_map["Fall"], level_map["Spring"])] <- c(1, -1)

# Run contrasts

results_contrast  <- maaslin_contrast_test(
  maaslin3_fit = m_contrast,
  contrast_mat = contrast_mat,
  correction = "BH",
  max_significance = FDR
)

# Extract contrasts

prev_vs_winter <- m_contrast$fit_data_prevalence$results %>%
  filter(metadata == "season_collection") %>%
  transmute(
    feature,
    contrast = paste0(value, "_vs_Winter"),
    coef,
    stderr,
    pval_individual,
    N,
    N_not_zero
  )

results_contrast_list <- list(
  Summer_vs_Fall   = results_contrast$fit_data_prevalence$results %>% filter(test == "Summer_vs_Fall"),
  Summer_vs_Spring = results_contrast$fit_data_prevalence$results %>% filter(test == "Summer_vs_Spring"),
  Fall_vs_Spring   = results_contrast$fit_data_prevalence$results %>% filter(test == "Fall_vs_Spring"),
  Summer_vs_Winter = prev_vs_winter %>% filter(grepl("Summer", contrast)),
  Spring_vs_Winter = prev_vs_winter %>% filter(grepl("Spring", contrast)),
  Fall_vs_Winter   = prev_vs_winter %>% filter(grepl("Fall", contrast))
)

results_contrast_list <- results_contrast_list[order(names(results_contrast_list))]

# Explore results

data.frame(
  variable = names(results_contrast_list),
  n_rows   = sapply(results_contrast_list, nrow),
  has_data = sapply(results_contrast_list, function(x) nrow(x) > 0)
)

# Combine datasets

prev_all <- bind_rows(prev_vs_winter, results_contrast$fit_data_prevalence$results %>%
                        filter(grepl("vs", test)) %>%
                        transmute(
                          feature,
                          contrast = test,
                          coef,
                          stderr,
                          pval_individual,
                          N,
                          N_not_zero
                        ))

# Taxa significant in group-wise model

sig_taxa <- c("Ascomycota_unclassified", "Curvularia", "Phlebia", "Pleosporales_unclassified")

# Recompute q-val within each contrast

results_contrast_season <- prev_all %>%
  filter(feature %in% sig_taxa) %>%
  group_by(contrast) %>%
  mutate(qval_correct = p.adjust(pval_individual, method = "BH")) %>%
  ungroup()

# Explore only significant results

results_contrast_season %>%
  filter(qval_correct<FDR) %>%
  print()

# Create table of results for significant taxa

results_contrast_prev <- results_contrast_season %>%
  mutate(
    coef_r = round(coef, 2),
    ci_low = coef - 1.96 * stderr,
    ci_high = coef + 1.96 * stderr,
    ci_low_r = round(ci_low, 2),
    ci_high_r = round(ci_high, 2),
    p_value = vapply(pval_individual, format_p, character(1)),
    q_value = vapply(qval_correct, format_p, character(1)),
    contrast_label = case_when(
      contrast == "Summer_vs_Winter" ~ "Summer vs. winter",
      contrast == "Summer_vs_Fall" ~ "Summer vs. fall",
      contrast == "Summer_vs_Spring" ~ "Summer vs. spring",
      contrast == "Spring_vs_Winter" ~ "Spring vs. winter",
      contrast == "Fall_vs_Winter" ~ "Fall vs. winter",
      contrast == "Fall_vs_Spring" ~ "Fall vs. spring")
    ) %>%
  select(feature, contrast_label, coef_r, ci_low_r, ci_high_r, N, N_not_zero, p_value, q_value) %>%
  arrange(feature, contrast_label) %>%
  rename(
    Genus = feature,
    Contrast = contrast_label,
    "MaAsLin3 effect size" = coef_r,
    "Lower 95% confidence interval" = ci_low_r,
    "Upper 95% confidence interval" = ci_high_r,
    "Prevalence count" = N_not_zero,
    "p-value" = p_value,
    "q-value" = q_value
  )

# ================================================================
# -------------------- EXPORT ------------------------------------
# ================================================================

write.csv(results_season_prev, file.path(DIR_TBL, "Table_diff_prevalence.csv"), row.names = FALSE)
write.csv(results_contrast_prev, file.path(DIR_TBL, "Table_contrasts.csv"), row.names = FALSE)
saveRDS(results_contrast_season, file.path(DIR_OUT, "diff_prevalence_results_contrasts_tbl.rds"))
