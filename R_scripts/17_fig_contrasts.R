# ================================================================
# Figure - Forest plot for contrasts
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

df_contrasts <- readRDS(file.path(DIR_OUT, "diff_prevalence_results_contrasts_tbl.rds"))

# ================================================================
# -------------------- CREATING FIGURE FOR CONTRAST --------------
# ================================================================

# Parameters

sig_taxa <- c("Phlebia", "Ascomycota_unclassified", "Pleosporales_unclassified", "Curvularia")

# Order for contrasts

contrast_order <- c("Summer_vs_Winter", "Summer_vs_Fall", "Summer_vs_Spring", "Spring_vs_Winter", "Fall_vs_Winter", "Fall_vs_Spring")

# Data wrangling

tbl_contrasts <- df_contrasts %>%
  filter(feature %in% sig_taxa) %>%
  mutate(
    ci_low = coef - 1.96 * stderr,
    ci_high = coef + 1.96 * stderr,
    sig = qval_correct < FDR,
    feature = factor(feature, levels = sig_taxa),
    contrast = factor(contrast, levels = contrast_order),
    contrast_label = case_when(
      contrast == "Summer_vs_Winter" ~ "Summer vs. winter",
      contrast == "Summer_vs_Fall" ~ "Summer vs. fall",
      contrast == "Summer_vs_Spring" ~ "Summer vs. spring",
      contrast == "Spring_vs_Winter" ~ "Spring vs. winter",
      contrast == "Fall_vs_Winter" ~ "Fall vs. winter",
      contrast == "Fall_vs_Spring" ~ "Fall vs. spring"
      ))

min(tbl_contrasts$ci_low)
max(tbl_contrasts$ci_high)

# Ascomycota

ascomycota <- tbl_contrasts %>%
  filter(feature == "Ascomycota_unclassified") %>%
  ggplot(aes(x = contrast_label, y = coef, color = sig)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.4, linewidth = 0.5) +
  geom_point(size = 5, shape = 15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(y = "MaAsLin3 effect size (95% confidence interval)") +
  coord_flip() +
  scale_color_manual(
    values = c(`FALSE` = "grey60", `TRUE` = "#d62728"),
    name = paste0("FDR < ", FDR),
    labels = c(`FALSE` = "No", `TRUE` = "Yes")
  ) +
  scale_y_continuous(
    limits = c(-4.5, 6.5),
    breaks = c(-2.5, 0, 2.5, 5)
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, face = "plain", color = "black"),
    axis.text.y = element_text(size = 12, face = "plain", color = "black")
  )

ascomycota

# Curvularia

curvularia <- tbl_contrasts %>%
  filter(feature == "Curvularia") %>%
  ggplot(aes(x = contrast_label, y = coef, color = sig)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.4, linewidth = 0.5) +
  geom_point(size = 5, shape = 15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(y = "MaAsLin3 effect size (95% confidence interval)") +
  coord_flip() +
  scale_color_manual(
    values = c(`FALSE` = "grey60", `TRUE` = "#d62728"),
    name = paste0("FDR < ", FDR),
    labels = c(`FALSE` = "No", `TRUE` = "Yes")
  ) +
  scale_y_continuous(
    limits = c(-4.5, 6.5),
    breaks = c(-2.5, 0, 2.5, 5)
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, face = "plain", color = "black"),
    axis.text.y = element_text(size = 12, face = "plain", color = "black")
  )

curvularia

# Phlebia

phlebia <- tbl_contrasts %>%
  filter(feature == "Phlebia") %>%
  ggplot(aes(x = contrast_label, y = coef, color = sig)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.4, linewidth = 0.5) +
  geom_point(size = 5, shape = 15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(y = "MaAsLin3 effect size (95% confidence interval)") +
  coord_flip() +
  scale_color_manual(
    values = c(`FALSE` = "grey60", `TRUE` = "#d62728"),
    name = paste0("FDR < ", FDR),
    labels = c(`FALSE` = "No", `TRUE` = "Yes")
  ) +
  scale_y_continuous(
    limits = c(-4.5, 6.5),
    breaks = c(-2.5, 0, 2.5, 5)
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, face = "plain", color = "black"),
    axis.text.y = element_text(size = 12, face = "plain", color = "black")
  )

phlebia

# Pleosporales

pleosporales <- tbl_contrasts %>%
  filter(feature == "Pleosporales_unclassified") %>%
  ggplot(aes(x = contrast_label, y = coef, color = sig)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.4, linewidth = 0.5) +
  geom_point(size = 5, shape = 15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(y = "MaAsLin3 effect size (95% confidence interval)") +
  coord_flip() +
  scale_color_manual(
    values = c(`FALSE` = "grey60", `TRUE` = "#d62728"),
    name = paste0("FDR < ", FDR),
    labels = c(`FALSE` = "No", `TRUE` = "Yes")
  ) +
  scale_y_continuous(
    limits = c(-4.5, 6.5),
    breaks = c(-2.5, 0, 2.5, 5)
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, face = "plain", color = "black"),
    axis.text.y = element_text(size = 12, face = "plain", color = "black")
  )

pleosporales

# Paneling

fig_combined <- (ascomycota / curvularia / phlebia / pleosporales) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

fig_combined

# ================================================================
# -------------------- EXPORT ------------------------------------
# ================================================================

# SVG

svglite(file.path(DIR_FIG, "Fig_panel_contrasts.svg"), width = 6, height = 9)
print(fig_combined)
dev.off()

# PDF

pdf(file.path(DIR_FIG, "Fig_panel_contrasts.pdf"), width = 6, height = 9)
print(fig_combined)
dev.off()
