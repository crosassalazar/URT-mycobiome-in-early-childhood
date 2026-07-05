# ================================================================
# Figure - Boxplot alpha season of sample collection
# Figure - NMDS beta season of sample collection
# Figure - Heatmap differential abundance season of sample collection
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

ps <- readRDS(file.path(DIR_OUT, "ps_counts.rds"))
alpha_df <- readRDS(file.path(DIR_OUT, "alpha_fig_df.rds"))
beta_df <- readRDS(file.path(DIR_OUT, "beta_fig_df.rds"))
beta_centroids <- readRDS(file.path(DIR_OUT, "beta_centroids_season_tbl.rds"))

# ================================================================
# -------------------- ALPHA_DIVERSITY BOX-AND-WHISKER PLOT ------
# ================================================================

fig_bp <- ggplot(
  subset(alpha_df, !is.na(season_collection) & !is.na(Observed)),
  aes(x = season_collection, y = Observed, fill = season_collection)
) +
  stat_boxplot(geom = "errorbar", width = 0.25) +
  geom_boxplot(color = "black", width = 0.7, outlier.shape = NA) +
  geom_point(
    shape = 1, size = 3, color = "black",
    position = position_jitter(width = 0.1),
    show.legend = FALSE
  ) +
  stat_summary(
    fun = mean, geom = "point",
    shape = 23, size = 3,
    fill = "white", color = "black",
    show.legend = FALSE
  ) +
  scale_fill_manual(
    name = "Season of sample collection",
    labels = c("Winter", "Fall", "Spring", "Summer"),
    values = c("#B3CDE3", "#CCEBC5", "#FED9A6", "#DECBE4")
  ) +
  labs(
    x = NULL,
    y = "Observed Richness index"
  ) +
  annotate(
    "text",
    x = 1.0,
    y = max(alpha_df$Observed, na.rm = TRUE) * 1.05,
    label = "p=0.08",
    size = 4
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

print(fig_bp)

# ================================================================
# -------------------- BETA-DIVERSITY NMDS PLOT ------------------
# ================================================================

fig_nmds <- ggplot(data = beta_df,
                   aes(x = NMDS1, 
                       y = NMDS2)) +
  stat_ellipse(level = 0.95,
               geom = "polygon",
               aes(color = season_collection,
                   fill = season_collection),
               alpha = 0.2,
               linewidth = 0.8,
               linetype= "solid") +
  geom_point(aes(fill = season_collection),
             shape = 21, 
             size = 4) +
  scale_color_manual(name = "Season of sample collection",
                     values = c("#B3CDE3", "#CCEBC5", "#FED9A6", "#DECBE4"),
                     labels = c("Winter","Fall","Spring","Summer")) + 
  scale_fill_manual(name = "Season of sample collection",
                    values = scales::alpha(c("#B3CDE3", "#CCEBC5", "#FED9A6", "#DECBE4"), 0.6),
                    labels = c("Winter","Fall","Spring","Summer")) + 
  geom_point(data = beta_centroids,
             aes(x = NMDS1, y = NMDS2),
             size = 9,
             shape = 21,
             colour = "black",
             fill = c("#7ba8ce", "#9ad78c", "#fdb95a", "#be99ca")) +
  labs(x = "NMDS1",
       y = "NMDS2") +
  annotate("text",
           x = -220, y = Inf,
           label = "p=0.09",
           hjust = 1.1, vjust = 1.5,
           size = 4) +
  theme(aspect.ratio = 1,        
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(size = 12, hjust = 0.0),
        legend.text = element_text(size = 12)
  )

print(fig_nmds)

# ================================================================
# -------------------- DIFFERENTIAL ABUNDANCE HEATMAP ------------
# ================================================================

# Collapse to genus

genus <- tax_glom(ps, taxrank = "Genus")

# Check

if (any(taxa_sums(genus) == 0)) stop("Zero-count taxa detected.")
if (any(sample_sums(genus) == 0)) stop("Zero-count samples detected.")

# Rename taxa

taxa <- as.data.frame(tax_table(genus))

taxa_names(genus) <- taxa$Genus

# Extract otu table

otu <- data.frame(otu_table(genus))

# Extract metadata

metadata <- data.frame(sample_data(genus))

# Significant taxa

sig_taxa <- c("Ascomycota_unclassified", "Curvularia", "Phlebia", "Pleosporales_unclassified")

# Build table for figure

tbl_diff <- otu %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(cols = any_of(sig_taxa),
               names_to = "taxon",
               values_to = "abundance") %>%
  mutate(present = as.integer(abundance > 0)) %>%
  left_join(
    tibble(
      sample_id = rownames(metadata),
      season_collection = metadata$season_collection
    ),
    by = "sample_id"
  )

# Compute prevalence per taxon × season

tbl_prev <- tbl_diff %>%
  group_by(taxon, season_collection) %>%
  dplyr::summarize(
    prevalence = 100 * mean(present, na.rm = TRUE),
    .groups = "drop"
  )

max(tbl_prev$prevalence)

# Factor conversion

tbl_prev$season_collection <- factor(tbl_prev$season_collection,
                                     levels = c("Winter", "Fall", "Spring", "Summer"))

# Heatmap

Hmisc::describe(as.numeric(tbl_prev$prevalence))

annot_df <- tbl_prev %>%
  distinct(taxon) %>%
  mutate(
    season_collection = NA,         # not used
    x_pos = length(levels(tbl_prev$season_collection)) + 0.5,
    label = "q=0.09"
  )

fig_hm <- ggplot(tbl_prev, aes(x = season_collection, y = taxon, fill = prevalence)) +
  geom_tile(color = "white", width = 0.99, height = 0.99) +
  geom_text(aes(label = sprintf("%.1f%%", prevalence)),
            size = 4, color = "white") +
  geom_text(
    data = annot_df,
    aes(x = x_pos, y = taxon, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4
  ) +
  scale_fill_gradientn(
    colors = c("#03051A", "#FEC980", "#9E0142"),
    limits = c(0, 40),
    name = "Prevalence (%)"
  ) +
  labs(x = "", y = "", title = "") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = c(1.25, 0.5),
    legend.justification = c("left", "center"),
    plot.margin = margin(5.5, 160, 5.5, 5.5),
    axis.text.x = element_text(hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black")
  ) +
  coord_fixed(ratio = 0.5, clip = "off")

print(fig_hm)

# ================================================================
# -------------------- PANELING ----------------------------------
# ================================================================

fig_combined_top <- (fig_bp | fig_nmds) +
  # plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

# ================================================================
# -------------------- EXPORT ------------------------------------
# ================================================================

# SVG

svglite(file.path(DIR_FIG, "Fig_panel_season_collection_top.svg"), width = 12, height = 5)
print(fig_combined_top)
dev.off()

svglite(file.path(DIR_FIG, "Fig_panel_season_collection_bottom.svg"), width = 11, height = 3)
print(fig_hm)
dev.off()

# PDF

pdf(file.path(DIR_FIG, "Fig_panel_season_collection_top.pdf"), width = 12, height = 5)
print(fig_combined_top)
dev.off()

pdf(file.path(DIR_FIG, "Fig_panel_season_collection_bottom.pdf"), width = 11, height = 3)
print(fig_hm)
dev.off()

# Further paneling in Inkscape
