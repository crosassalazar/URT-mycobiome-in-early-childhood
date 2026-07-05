# ================================================================
# Alpha-diversity analyses
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
load_pkg("car")

# ================================================================
# -------------------- LOAD PHYLOSEQ -----------------------------
# ================================================================

ps <- readRDS(file.path(DIR_OUT, "ps_rar_min.rds")) ########## Import normalized dataset

# ================================================================
# -------------------- ALPHA MEASURES AND DATA WRANGLING ---------
# ================================================================

# Extract metadata

metadata <- data.frame(sample_data(ps))

# Estimate richness

alpha <- estimate_richness(ps, measures = c("Observed"))

# Merge

df <- metadata %>%
  rownames_to_column("SampleID") %>%
  left_join(alpha %>% rownames_to_column("SampleID"), by = "SampleID")

# Factor conversion and relevel

factor_vars <- c(
  "efsex","race_white_dev","efdelivery",
  "any_breastfd_dev","mat_asthma_dev","young_sib_dev",
  "early_life_anti_dev","early_life_smoke_dev",
  "insurance_dev","petany_dev","efrural_dev",
  "any_lrti_dev","season_collection","mold_dev"
)

df[factor_vars] <- lapply(df[factor_vars], as.factor)

df$season_collection <- relevel(df$season_collection, ref = "Winter")

# ================================================================
# -------------------- UNIVARIATE MODELS -------------------------
# ================================================================

preds <- c(
  "efdelivery","any_breastfd_dev","mat_asthma_dev",
  "young_sib_dev","early_life_anti_dev","early_life_smoke_dev",
  "insurance_dev","petany_dev","efrural_dev",
  "any_lrti_dev","mold_dev","season_collection"
)

results_single <- map_df(preds, function(v) {
  lm(as.formula(paste("Observed ~", v)), data = df) %>%
    tidy(conf.int = TRUE) %>%
    mutate(predictor = v)
})

results_single %>%
  filter(term!="(Intercept)") %>%
  filter(p.value<0.20) %>%
  print()

# ================================================================
# -------------------- MULTIVARIATE MODELS -----------------------
# ================================================================

selected <- c(
  "scaled_age_sample","efsex","race_white_dev",
  "young_sib_dev","early_life_anti_dev","season_collection"
)

fit_multi <- lm(
  as.formula(paste("Observed ~", paste(selected, collapse = " + "))),
  data = df
)

results_multi <- tidy(fit_multi, conf.int = TRUE)

results_multi %>%
  filter(term!="(Intercept)") %>%
  filter(p.value<SIG) %>%
  print()

# ================================================================
# -------------------- UNIVARIATE TERM P-VALUE FOR SEASON --------
# ================================================================

# Term p-value for season of sample collection to be used in figure

options(contrasts = c("contr.sum", "contr.poly"))

season_p_uni <- lm(Observed ~ season_collection, data = df)
car::Anova(season_p_uni, type = 3)

season_p_multi <- lm(Observed ~ season_collection + early_life_anti_dev + efsex + race_white_dev + scaled_age_sample,
                     data = df, na.action = na.omit)

car::Anova(season_p_multi, type = 3)

# ================================================================
# -------------------- EXPORT ------------------------------------
# ================================================================

saveRDS(df, file.path(DIR_OUT, "alpha_fig_df.rds"))
saveRDS(results_single, file.path(DIR_OUT, "alpha_univariate_results_tbl.rds"))
saveRDS(results_multi, file.path(DIR_OUT, "alpha_multivariate_results_tbl.rds"))
