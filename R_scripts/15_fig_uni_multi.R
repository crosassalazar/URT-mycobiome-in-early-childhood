# ================================================================
# Figure - Heatmap of alpha- and beta-diversity
# Figure - Forest plot of alpha-diversity
# Figure - Bar plot of beta-diversity
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

load_pkg("ggpubr")
load_pkg("scales")
load_pkg("patchwork")

# ================================================================
# -------------------- LOAD DATASETS -----------------------------
# ================================================================

alpha_uni <- readRDS(file.path(DIR_OUT, "alpha_univariate_results_tbl.rds"))
alpha_multi <- readRDS(file.path(DIR_OUT, "alpha_multivariate_results_tbl.rds"))
beta_uni <- readRDS(file.path(DIR_OUT, "beta_univariate_results_tbl.rds"))
beta_multi <- readRDS(file.path(DIR_OUT, "beta_multivariate_results_tbl.rds"))

# ================================================================
# -------------------- HEATMAP OF ALPHA- AND BETA-DIVERSITY ------
# ================================================================

# ---------------------------
# 1) Label map
# ---------------------------

label_map_hm <- c(
  efdelivery                   = "Birth by C-section",
  any_breastfd_dev             = "Any breastfeeding",
  mat_asthma_dev               = "Maternal asthma",
  young_sib_dev                = "Presence of household siblings ages\n6 years or younger at birth",
  early_life_anti_dev          = "Exposure to antibiotics in utero or\nduring early infancy",
  early_life_smoke_dev         = "Exposure to tobacco smoke in utero or\nduring early infancy",
  insurance_dev                = "Private insurance during early infancy",
  petany_dev                   = "Pet ownership during early infancy",
  efrural_dev                  = "Residence in rural area during early infancy",
  any_lrti_dev                 = "History of lower respiratory tract infection\nin the first year of life",
  mold_dev                     = "Residential mold damage during early infancy",
  season_collection            = "Season of sample collection"
)

# ---------------------------
# 2) Format p-values
# ---------------------------

round_half_up_uni <- function(x, digits = 0) {
  pos <- 10^digits
  floor(x * pos + 0.5) / pos
}

format_p_uni <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  decimals <- if (p < 0.20) 2 else 1
  rounded  <- round_half_up_uni(p, digits = decimals)
  sprintf(paste0("%.", decimals, "f"), rounded)
}

# ---------------------------
# 3) Build table for figure
# ---------------------------

alpha_uni <- alpha_uni %>%
  filter(predictor != "season_collection", term != "(Intercept)") %>%
  select(predictor, p.value) %>%
  rename(p_value_alpha = p.value) %>%
  add_row(predictor = "season_collection", p_value_alpha = 0.08)

beta_uni <- beta_uni %>%
  select(term, p) %>%
  rename(predictor = term, p_value_beta = p)

tbl_hm <- left_join(alpha_uni, beta_uni, by = "predictor") %>%
  mutate(
    p_fmt_alpha = vapply(p_value_alpha, format_p_uni, character(1)),
    p_fmt_beta  = vapply(p_value_beta, format_p_uni, character(1))
  ) %>%
  pivot_longer(
    cols = c(p_value_alpha, p_value_beta),
    names_to = "metric",
    values_to = "p_value"
  ) %>%
  mutate(
    label = ifelse(metric == "p_value_alpha", p_fmt_alpha, p_fmt_beta),
    metric_pretty = factor(
      metric,
      levels = c("p_value_alpha", "p_value_beta"),
      labels = c("\u03B1-diversity", "\u03B2-diversity")
    ),
    predictor = factor(predictor, levels = unique(predictor))
  )

# ---------------------------
# 4) Figure
# ---------------------------

fig_hm <- ggplot(tbl_hm, aes(x = metric_pretty, y = predictor, fill = p_value)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 4.0) +
  scale_fill_gradientn(
    name = "p-value",
    colours = c("#08306B", "#2171B5", "#6BAED6", "#D9D9D9", "#F2F2F2"),
    values = scales::rescale(c(0.00, 0.10, 0.20, 0.60, 1.00)),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = sprintf("%.1f", seq(0, 1, by = 0.2))
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(
    limits = rev(levels(tbl_hm$predictor)),
    labels = label_map_hm,
    expand = c(0, 0),
    drop = FALSE
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right")

print(fig_hm)

# ================================================================
# -------------------- FOREST PLOT OF ALPHA-DIVERSITY ------------
# ================================================================

# ---------------------------
# 1) Label map
# ---------------------------

label_map_fp <- c(
  scaled_age_sample            = "Age at sample collection (years)",
  efsexmale                    = "Male sex",
  race_white_dev1              = "Non-Hispanic White",
  young_sib_dev1               = "Presence of household siblings ages\n6 years or younger at birth",
  early_life_anti_dev1         = "Exposure to antibiotics in utero or\nduring early infancy",
  season_collectionWinter      = "Season of sample collection: Winter",
  season_collectionSpring      = "Spring",
  season_collectionSummer      = "Summer",
  season_collectionFall        = "Fall"
)

# ---------------------------
# 2) Build table for figure
# ---------------------------

tbl_fp <- alpha_multi %>%
  filter(term != "(Intercept)") %>%
  mutate(
    sig = as.integer(p.value < SIG),
    p_fmt = vapply(p.value, format_p, character(1)),
    `Estimate (95% CI), p-value` = sprintf(
      "%.2f (%.2f, %.2f), p=%s",
      estimate, conf.low, conf.high, p_fmt
    )
  ) %>%
  add_row(
    term = "season_collectionWinter",
    estimate = 0, conf.low = 0, conf.high = 0,
    p.value = 1, sig = 0,
    `Estimate (95% CI), p-value` = "Reference"
  )

current_levels <- unique(as.character(tbl_fp$term))
no_winter <- setdiff(current_levels, "season_collectionWinter")
fall_pos <- match("season_collectionFall", no_winter)

desired_levels <- if (!is.na(fall_pos)) {
  append(no_winter, "season_collectionWinter", after = fall_pos - 1)
} else {
  c(no_winter, "season_collectionWinter")
}

tbl_fp <- tbl_fp %>%
  mutate(term = factor(term, levels = desired_levels)) %>%
  arrange(term)

# ---------------------------
# 3) Figure
# ---------------------------

x_min <- min(tbl_fp$conf.low, na.rm = TRUE)
x_max <- max(tbl_fp$conf.high, na.rm = TRUE)
x_range <- x_max - x_min

right_pad <- x_max + 0.08 * x_range

fig_fp <- ggplot(tbl_fp, aes(x = estimate, y = term, color = factor(sig))) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.25, color = "grey60") +
  geom_point(shape = 15, size = 5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey85") +
  geom_text(
    aes(x = right_pad, y = term, label = `Estimate (95% CI), p-value`),
    hjust = 0, size = 4.2, color = "black", inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = c(`0` = "grey60", `1` = "#d62728"),
    labels = c(`0` = "p>=0.05", `1` = "p<0.05"),
    name = NULL
  ) +
  coord_cartesian(xlim = c(-5, 7.5), clip = "off") +
  scale_x_continuous(breaks = c(-5, -2.5, 0, 2.5, 5, 7.5), expand = c(0, 0)) +
  labs(x = "β coefficient (95% confidence interval)", y = NULL) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.line.x = element_line(color = "black"),
    plot.margin = margin(t = 10, r = 180, b = 10, l = 60),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(3, "pt")
  ) +
  scale_y_discrete(
    limits = rev(levels(tbl_fp$term)),
    labels = label_map_fp
  )

print(fig_fp)

# ================================================================
# -------------------- BARPLOT OF BETA-DIVERSITY -----------------
# ================================================================

label_map_bp <- c(
  scaled_age_sample            = "Age at sample collection (years)",
  efsex                        = "Male sex",
  race_white_dev               = "Non-Hispanic White",
  any_breastfd_dev             = "Any breastfeeding",
  early_life_anti_dev          = "Exposure to antibiotics in utero or\nduring early infancy",
  any_lrti_dev                 = "History of lower respiratory tract infection\nin the first year of life",
  season_collection            = "Season of sample collection"
)

format_r2 <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = 3))
}

tbl_bp <- beta_multi %>%
  filter(!term %in% c("Residual", "Total")) %>%
  mutate(
    p_value = `Pr(>F)`,
    sig = factor(as.integer(p_value < SIG), levels = c(1, 0),
                 labels = c("Significant", "NS")),
    p_fmt = vapply(p_value, format_p, character(1)),
    stat = paste0("R\u00B2=", format_r2(R2), ", p=", p_fmt)
  )

tbl_bp$label_rev <- factor(tbl_bp$term, levels = rev(unique(tbl_bp$term)))

fig_bp <- ggplot(tbl_bp, aes(x = R2, y = label_rev)) +
  geom_col(aes(fill = sig)) +
  geom_text(aes(x = Inf, label = stat), hjust = -0.05, size = 4) +
  scale_x_continuous(
    labels = scales::label_number(accuracy = 0.001),
    expand = expansion(mult = c(0.02, 0.35))
  ) +
  scale_fill_manual(values = c("Significant" = "#d62728", "NS" = "gray60")) +
  labs(
    x = expression(R^2 ~ "(proportion of variance explained)"),
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 5, r = 120, b = 5, l = 60),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(3, "pt")
  ) +
  scale_y_discrete(labels = label_map_bp)

print(fig_bp)

# ================================================================
# -------------------- PANELING ----------------------------------
# ================================================================

fig_combined <- (fig_fp / fig_bp) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

print(fig_combined)

# ================================================================
# -------------------- EXPORT ------------------------------------
# ================================================================

# SVG

svglite(file.path(DIR_FIG, "Fig_heatmap_uni.svg"), width = 7.5, height = 3.6)
print(fig_hm)
dev.off()

svglite(file.path(DIR_FIG, "Fig_panel_multi.svg"), width = 9.7, height = 7.5)
print(fig_combined)
dev.off()

# PDF

pdf(file.path(DIR_FIG, "Fig_heatmap_uni.pdf"), width = 7.5, height = 3.6)
print(fig_hm)
dev.off() # Note that there is an issue with labels in PDF

pdf(file.path(DIR_FIG, "Fig_panel_multi.pdf"), width = 9.7, height = 7.5)
print(fig_combined)
dev.off() # Note that there is an issue with labels in PDF
