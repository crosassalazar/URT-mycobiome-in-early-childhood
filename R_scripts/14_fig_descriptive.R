# ================================================================
# Figure - Histogram reads
# Figure - Histogram genera
# Figure - Stacked bar plot
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

load_pkg("ggthemes")

# ================================================================
# -------------------- LOAD DATASETS -----------------------------
# ================================================================

hist_read_df <- readRDS(file.path(DIR_OUT, "descriptive_histogram_reads_fig_df.rds"))
hist_gen_df <- readRDS(file.path(DIR_OUT, "descriptive_histogram_genus_fig_df.rds"))
sbp_gen_df <- readRDS(file.path(DIR_OUT, "descriptive_stacked_barplot_genus_fig_df.rds"))

# ================================================================
# -------------------- CREATING HISTOGRAM FOR READS --------------
# ================================================================

fig_hist_lib <- ggplot(hist_read_df, aes(x = readdepth)) +
  geom_histogram(
    bins = 100, 
    fill = "grey60", 
    color = "black") +
  scale_x_continuous(
    breaks = seq(0, max(hist_read_df$readdepth, na.rm = TRUE), by = 10000),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    expand = c(0, 0)) +
  theme_minimal() +
  labs(
    x = "Total number of reads per sample",
    y = "Frequency"
  ) +
  geom_vline(
    xintercept = median(hist_read_df$readdepth),
    linetype = "dotted",
    color = "black",
    linewidth = 0.5
  ) +
  annotate(
    "text",
    x = median(hist_read_df$readdepth),
    y = Inf,
    label = "tilde(x)~'('~25^th~' percentile, '~75^th~' percentile) = 8,869 (3,862, 13,625)'",
    parse = TRUE,
    hjust = -0.05,
    vjust = 1.5,
    size = 4
  ) +
  theme(
    legend.position = "none",
    panel.grid = element_line(color = "grey85"),
    axis.title = element_text(size = 12, color = "black"),
    axis.text  = element_text(size = 12, color = "black"),
    axis.line  = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.border = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

print(fig_hist_lib)

# ================================================================
# -------------------- CREATING HISTOGRAM FOR GENERA -------------
# ================================================================

fig_hist_gen <- ggplot(hist_gen_df, aes(x = observed)) +
  geom_histogram(
    aes(y = after_stat(count)),
    fill = "grey60",
    color = "white",
    alpha = 0.6,
    bins = 35
  ) +
  geom_density(
    aes(y = after_stat(scaled) * max(after_stat(count))),
    color = "black",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = median(hist_gen_df$observed),
    linetype = "dotted",
    color = "black",
    linewidth = 0.5
  ) +
  labs(x = "Total number of genera per sample", 
       y = "Frequency") +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = seq(0, 40, by = 5)
  ) +
  scale_y_continuous(
    expand = c(0, 0)) +
  theme_minimal() +
  annotate(
    "text",
    x = median(hist_gen_df$observed),
    y = Inf,
    label = "tilde(x)~'('~25^th~' percentile, '~75^th~' percentile) = 3 (1, 5)'",
    parse = TRUE,
    hjust = -0.05,
    vjust = 1.5,
    size = 4
  ) +
  theme(
    legend.position = "none",
    panel.grid = element_line(color = "grey85"),
    axis.title = element_text(size = 12, color = "black"),
    axis.text  = element_text(size = 12, color = "black"),
    axis.line  = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.border = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

print(fig_hist_gen)

# ================================================================
# -------------------- CREATING STACKED BARPLOT ------------------
# ================================================================

fig_sb <- ggplot(sbp_gen_df,
                 aes(fill = top_gen_plot,
                     y = total_pct,
                     x = sample)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = "Participant samples", y = "Relative abundance (%)") +
  scale_fill_manual(
    name   = "Genera",
    values = tableau_color_pal(palette = "Tableau 20", direction = -1)(length(levels(sbp_gen_df$top_gen_plot)))
  ) +
  scale_y_continuous(expand = c(0, 2)) +
  theme(
    axis.title.x   = element_text(size = 9),
    axis.title.y   = element_text(size = 9),
    axis.text.y    = element_text(size = 7),
    axis.text.x    = element_blank(),
    axis.ticks.x   = element_blank(),
    legend.title   = element_text(size = 8),
    legend.text    = element_text(size = 7, face = "italic"),
    strip.background = element_rect(fill = "#e2e2e2", colour = NA),
    legend.key.size = unit(0.7, "line")
  )

print(fig_sb)

# ================================================================
# -------------------- EXPORT ------------------------------------
# ================================================================

# SVG

svglite(file.path(DIR_FIG, "Fig_histogram_reads.svg"), width = 10, height = 5)
print(fig_hist_lib)
dev.off()

svglite(file.path(DIR_FIG, "Fig_histogram_genus.svg"), width = 10, height = 5)
print(fig_hist_gen)
dev.off()

svglite(file.path(DIR_FIG, "Fig_stacked_barplot_top_genus.svg"), width = 10, height = 4)
print(fig_sb)
dev.off()

# PDF

pdf(file.path(DIR_FIG, "Fig_histogram_reads.pdf"), width = 10, height = 5)
print(fig_hist_lib)
dev.off()

pdf(file.path(DIR_FIG, "Fig_histogram_genus.pdf"), width = 10, height = 5)
print(fig_hist_gen)
dev.off()

pdf(file.path(DIR_FIG, "Fig_stacked_barplot_top_genus.pdf"), width = 10, height = 5)
print(fig_sb)
dev.off()
