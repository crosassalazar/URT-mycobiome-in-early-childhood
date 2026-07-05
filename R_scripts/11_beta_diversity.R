# ================================================================
# Beta-diversity analyses
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

load_pkg("broom")

# ================================================================
# -------------------- LOAD PHYLOSEQ -----------------------------
# ================================================================

ps <- readRDS(file.path(DIR_OUT, "ps_rar_min.rds")) ########## Import normalized dataset

# ================================================================
# -------------------- DATA WRANGLING ----------------------------
# ================================================================

# Predictors

factor_vars <- c(
  "efdelivery","any_breastfd_dev","mat_asthma_dev",
  "young_sib_dev","early_life_anti_dev","early_life_smoke_dev",
  "insurance_dev","petany_dev","efrural_dev",
  "any_lrti_dev","mold_dev","season_collection"
)

# Identify complete-cases (needed for downstream analyses)

Hmisc::describe(as.data.frame(sample_data(ps))) # Minimal missing in this dataset

keep <- complete.cases(as.data.frame(sample_data(ps))[, factor_vars, drop = FALSE])

# Prune samples to non-missing

ps <- prune_samples(keep, ps)

# Drop taxa that became zero after sample pruning

ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Check

if (any(taxa_sums(ps) == 0)) stop("Zero-count taxa detected.")
if (any(sample_sums(ps) == 0)) stop("Zero-count samples detected.")

# Factor conversion and relevel

sample_data(ps)[, factor_vars] <- lapply(sample_data(ps)[, factor_vars, drop = FALSE], factor)

sample_data(ps)$season_collection <- relevel(sample_data(ps)$season_collection, ref = "Winter")

# Extract metadata

metadata <- data.frame(sample_data(ps))

# ================================================================
# -------------------- BETA MEASURES AND NMDS ORDINATION ---------
# ================================================================

dist_jc <- phyloseq::distance(ps, method = "jaccard", binary = TRUE)

set.seed(124)

nmds <- ordinate(ps, 
                 method = "NMDS",
                 distance = "jaccard",
                 binary = TRUE,
                 k = 2, 
                 trymax = 999)

cat("NMDS stress:", nmds$stress, "\n")

stressplot(nmds)

plot_ordination(ps, nmds, color="season_collection", title="NMDS for Jaccard") + 
  theme(aspect.ratio=1) + stat_ellipse(aes(group = season_collection), linetype = 2)

# ================================================================
# -------------------- UNIVARIATE MODELS -------------------------
# ================================================================

preds <- c(
  "efdelivery","any_breastfd_dev","mat_asthma_dev",
  "young_sib_dev","early_life_anti_dev","early_life_smoke_dev",
  "insurance_dev","petany_dev","efrural_dev",
  "any_lrti_dev","mold_dev","season_collection"
)

set.seed(SEED)

res_univ <- lapply(preds, function(v) {
  fit <- adonis2(as.formula(paste("dist_jc ~", v)),
                 data = metadata,
                 permutations = 999, 
                 by = "margin")
  
  tibble(
    term = v,
    R2   = fit$R2[1],
    F    = fit$F[1],
    p    = fit$`Pr(>F)`[1]
  )
})

permanova_univ <- bind_rows(res_univ)

permanova_univ %>%
  filter(p<0.20) %>%
  print()

# ================================================================
# -------------------- MULTIVARIATE MODELS -----------------------
# ================================================================

set.seed(SEED)

res_multi <- adonis2(
  dist_jc ~ scaled_age_sample + efsex + race_white_dev +
    any_breastfd_dev + early_life_anti_dev +
    any_lrti_dev + season_collection,
  data = metadata,
  permutations = 999, 
  by = "margin"
)

permanova_multi <- res_multi %>%
  as.data.frame() %>%
  rownames_to_column("term")

permanova_multi %>%
  filter(`Pr(>F)`<SIG) %>%
  print()

# ================================================================
# -------------------- NMDS DATASETS FOR FIGURES -----------------
# ================================================================

nmds_coords <- as.data.frame(scores(nmds, display = "sites")) %>%
  rownames_to_column("SampleID")

nmds_metadata <- data.frame(sample_data(ps)) %>%
  rownames_to_column("SampleID")

nmds_fig_df <- left_join(nmds_coords, nmds_metadata, by = "SampleID")

# ================================================================
# -------------------- ESTIMATE CENTROIDS ------------------------
# ================================================================

get_centroids <- function(df, group_var) {
  df %>%
    group_by(.data[[group_var]]) %>%
    summarise(across(c(NMDS1, NMDS2), mean), .groups="drop")
}

centroids_season <- get_centroids(nmds_fig_df, "season_collection")

# ================================================================
# -------------------- EXPORT ------------------------------------
# ================================================================

saveRDS(nmds_fig_df, file.path(DIR_OUT, "beta_fig_df.rds"))
saveRDS(centroids_season, file.path(DIR_OUT, "beta_centroids_season_tbl.rds"))
saveRDS(permanova_univ, file.path(DIR_OUT, "beta_univariate_results_tbl.rds"))
saveRDS(permanova_multi, file.path(DIR_OUT, "beta_multivariate_results_tbl.rds"))
