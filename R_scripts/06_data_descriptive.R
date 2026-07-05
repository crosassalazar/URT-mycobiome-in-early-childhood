# ================================================================
# Data descriptive statistics at ASV level
# Rarefaction and coverage
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

# ================================================================
# -------------------- LOAD PHYLOSEQ -----------------------------
# ================================================================

ps_asv <- readRDS(file.path(DIR_OUT, "ps_counts.rds"))

# ================================================================
# -------------------- ASV LEVEL ---------------------------------
# ================================================================

# --------------------------------
# 1) Basic descriptive statistics
# --------------------------------

# Core counts

n_asv      <- ntaxa(ps_asv)
n_samples  <- nsamples(ps_asv)
total_reads <- sum(sample_sums(ps_asv))

cat("Number of ASVs:", n_asv, "\n")
cat("Number of samples:", n_samples, "\n")
cat("Total reads:", total_reads, "\n")

# Library size per sample

lib_sizes <- sample_sums(ps_asv)

Hmisc::describe(lib_sizes)

# Total reads per ASV and per sample

read_sums_df <- bind_rows(
  tibble(
    total_reads = sort(taxa_sums(ps_asv), decreasing = TRUE),
    order_place = seq_len(ntaxa(ps_asv)),
    type        = "ASVs"
  ),
  tibble(
    total_reads = sort(sample_sums(ps_asv), decreasing = TRUE),
    order_place = seq_len(nsamples(ps_asv)),
    type        = "Samples"
  )
)

plot_qc <- ggplot(read_sums_df, aes(x = order_place, y = total_reads)) +
  geom_bar(stat = "identity") +
  ggtitle("Total number of reads per ASV and per sample") +
  facet_wrap(~type, nrow = 1, scales = "free") +
  labs(x = "", y = "Total number of reads")

print(plot_qc)

# Number of ASVs per sample

ntaxabysample <- as.data.frame(
  apply(as.data.frame(otu_table(ps_asv)) > 0, 1, sum)
)

names(ntaxabysample)[1] <- "observed"

Hmisc::describe(ntaxabysample)

# ---------------------
# 2) Rarefaction curve
# ---------------------

# Create matrix

otu_mat <- as(otu_table(ps_asv), "matrix")

# Figure

rarecurve(
  otu_mat,
  step = 1000,
  cex = 0.6,
  label = FALSE,
  xlab = "Sequencing depth",
  ylab = "Observed ASVs",
  xaxt = "n",             # remove default axis
  xlim = c(0, quantile(lib_sizes, 0.75))      # force x-axis max
)

ticks <- seq(0, quantile(lib_sizes, 0.75), by = 1000)

axis(
  side = 1,
  at = ticks,
  labels = ticks
)

abline(v = min(lib_sizes), col = "blue", lty = 2, lwd = 2)
abline(v = 1000, col = "red", lty = 2, lwd = 2)

# Define depth grid

depth_seq <- seq(1, min(rowSums(otu_mat)), by = 100)

# Build rarefaction data

rare_list <- lapply(seq_len(nrow(otu_mat)), function(i) {
  r <- rarefy(x = otu_mat[i, ], sample = depth_seq)
  
  data.frame(
    sample   = rownames(otu_mat)[i],
    depth    = depth_seq,
    richness = as.numeric(r)
  )
})

rare_df <- bind_rows(rare_list)

# Compute slope

slope_check <- rare_df %>%
  group_by(sample) %>%
  arrange(depth) %>%
  mutate(slope = (richness - lag(richness)) / (depth - lag(depth))) %>%
  filter(depth >= 500 & depth <= 700) %>%
  summarise(mean_slope = mean(slope, na.rm = TRUE))

mean(slope_check$mean_slope < 0.01) # This gives % plateaued samples

# Mean slope < 0.01 indicates that fewer than ~1 additional taxon would be expected per 100 additional reads

summary(slope_check$mean_slope) 

# ---------------------
# 3) Predominant taxa
# ---------------------

top_asv_all <- tapply(
  taxa_sums(ps_asv),
  tax_table(ps_asv)[, "Species"],
  sum,
  na.rm = TRUE
) %>%
  sort(decreasing = TRUE) %>%
  as.data.frame()

names(top_asv_all)[1] <- "total_ct_sample"
top_asv_all$total_pct_all <- 100 * top_asv_all$total_ct_sample / sum(top_asv_all$total_ct_sample)

head(top_asv_all, 15) # Top 15 ASVs

# Proportion of reads represented by top 15 ASVs

prop_top15_all <- sum(top_asv_all$total_pct_all[1:15])
cat("Top 15 ASVs (% of reads):", prop_top15_all, "\n")

# --------------------------------
# 4) Histogram
# --------------------------------

hist(lib_sizes,
     main  = "Histogram of read counts per sample",
     xlab  = "Total number of reads per sample",
     border = "black",
     col    = "lightblue",
     breaks = 100,
     las    = 1)

# ---------------------------------
# 5) Stacked bar plot
# ---------------------------------

# Create new dataset for plotting

ps_plot_ct <- ps_asv

# Add "top_asv_plot" variable to taxonomy and replace with top 15 ASVs

tax_asv_plot <- as.data.frame(tax_table(ps_plot_ct))

tax_asv_plot$top_asv_plot <- "Other"

top_asv_ids <- names(sort(taxa_sums(ps_plot_ct), decreasing = TRUE)[1:15]) # Identify top 15 ASVs

tax_asv_plot[top_asv_ids, "top_asv_plot"] <- as.character(tax_asv_plot[top_asv_ids, "Species"])

print(head(tax_asv_plot))

# Replace taxonomy table

tax_table(ps_plot_ct) <- as.matrix(tax_asv_plot)

# Convert to relative abundance

ps_plot_rel <- transform_sample_counts(ps_plot_ct, function(x) 100 * x / sum(x))

# Stacked barplot

p_all_asv <- plot_bar(ps_plot_rel, fill = "top_asv_plot") +
  geom_bar(aes(color = top_asv_plot, fill = top_asv_plot),
           stat = "identity", position = "stack") +
  labs(x = "Samples", y = "Relative abundance (%)") +
  ggtitle("Top 15 ASVs in all") +
  scale_fill_discrete(name = "Species") +
  scale_color_discrete(name = "Species") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

print(p_all_asv)

# ================================================================
# -------------------- GENUS LEVEL -------------------------------
# ================================================================

# --------------------------------
# 1) Collapse to genus
# --------------------------------

ps_genus <- tax_glom(ps_asv, "Genus")

# --------------------------------
# 2) Check
# --------------------------------

if (any(taxa_sums(ps_genus) == 0)) stop("Zero-count taxa detected.")
if (any(sample_sums(ps_genus) == 0)) stop("Zero-count samples detected.")

# --------------------------------
# 3) Basic descriptive statistics
# --------------------------------

# Core counts

n_genera <- ntaxa(ps_genus)

cat("Number of genera:", n_genera, "\n")

# Number of genera per sample

ngenerabysample <- as.data.frame(
  apply(as.data.frame(otu_table(ps_genus)) > 0, 1, sum)
)

names(ngenerabysample)[1] <- "observed"

Hmisc::describe(ngenerabysample)

# --------------------------------------------
# 4) Looking for common taxa in other studies
# --------------------------------------------

# Abundance

top_gen_all <- tapply(
  taxa_sums(ps_genus),
  tax_table(ps_genus)[, "Genus"],
  sum,
  na.rm = TRUE
) %>%
  sort(decreasing = TRUE) %>%
  as.data.frame()

names(top_gen_all)[1] <- "total_ct_sample"
top_gen_all$total_pct_all <- 100 * top_gen_all$total_ct_sample / sum(top_gen_all$total_ct_sample)

target_genera <- c("Candida", 
                   "Debaryomyces",
                   "Trametes", 
                   "Malassezia",  
                   "Rigidoporus", # Only taxa not present
                   "Starmerella",
                   "Cladosporium",
                   "Helotiales_unclassified") # Note that this is an order

top_gen_all %>%
  tibble::rownames_to_column(var = "Genus") %>%
  mutate(row_id = row_number()) %>%
  filter(Genus %in% target_genera) %>%
  print()

# Prevalence

otu_common <- as.data.frame(otu_table(ps_genus))
tax_common <- as.data.frame(tax_table(ps_genus))

otu_common <- otu_common > 0

tax_common$ASV <- rownames(tax_common)

prev_counts_common <- tax_common %>%
  filter(Genus %in% target_genera) %>%
  group_by(Genus) %>%
  summarise(
    prevalence_count = sum(rowSums(otu_common[, ASV, drop = FALSE]) > 0),
    .groups = "drop"
  )

n_total <- nsamples(ps_genus)

prev_counts_common <- prev_counts_common %>%
  mutate(prevalence_pct = (prevalence_count / n_total) * 100)

print(prev_counts_common)

# ---------------------
# 5) Predominant taxa
# ---------------------

head(top_gen_all, 15) # Top 15 genera

# Proportion of reads represented by top 15 genera

prop_top15_all <- sum(top_gen_all$total_pct_all[1:15])
cat("Top 15 genera (% of reads):", prop_top15_all, "\n")

# ---------------------
# 6) Histogram
# ---------------------

hist(
  ngenerabysample$observed,
  main = "Histogram of observed genera per sample",
  xlab = "Total number of genera per sample",
  ylab = "Frequency",
  col = "lightblue",
  border = "black",
  breaks = 35,
  las = 1
)

# ---------------------------------
# 7) Stacked bar plot
# ---------------------------------

# Create new dataset for plotting

ps_plot_ct <- ps_genus

# Add "top_gen_plot" variable to taxonomy and replace with top 15 genera

tax_gen_plot <- as.data.frame(tax_table(ps_plot_ct))

tax_gen_plot$top_gen_plot <- "Other"

top_gen_ids <- names(sort(taxa_sums(ps_plot_ct), decreasing = TRUE)[1:15]) # Identify top 15 genera

tax_gen_plot[top_gen_ids, "top_gen_plot"] <- as.character(tax_gen_plot[top_gen_ids, "Genus"])

# Replace taxonomy table

tax_table(ps_plot_ct) <- as.matrix(tax_gen_plot)

# Convert to relative abundance

ps_plot_rel <- transform_sample_counts(ps_plot_ct, function(x) 100 * x / sum(x))

# Stacked barplot

p_all_genus <- plot_bar(ps_plot_rel, fill = "top_gen_plot") +
  geom_bar(aes(color = top_gen_plot, fill = top_gen_plot),
           stat = "identity", position = "stack") +
  labs(x = "Samples", y = "Relative abundance (%)") +
  ggtitle("Top 15 genera in all") +
  scale_fill_discrete(name = "Species") +
  scale_color_discrete(name = "Species") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

print(p_all_genus)

# ================================================================
# -------------------- OTHER TAXONOMIC LEVELS --------------------
# ================================================================

count_unique_by_rank <- function(x) {
  x <- trimws(x)
  x <- x[!is.na(x) & x != ""]
  length(unique(x))
}

tax_oth <- as.data.frame(tax_table(ps_asv))

tax_present <- sapply(colnames(tax_oth), function(r) count_unique_by_rank(tax_oth[[r]]))

tax_present

# ================================================================
# -------------------- CREATING DATASETS FOR FIGURES -------------
# ================================================================

# Creating dataset for publication-ready histogram for reads

hist_read_df <- data.frame(
  sample = names(lib_sizes),
  readdepth = lib_sizes,
  row.names = NULL
)

# Creating dataset for publication-ready histogram for genera

hist_gen_df <- data.frame(
  sample = rownames(ngenerabysample),
  observed = ngenerabysample$observed,
  row.names = NULL
)

# Creating dataset for publication-ready stacked bar plot

sbp_gen_df <- ps_plot_ct %>%
  psmelt() %>%
  rename(sample = Sample) %>%
  group_by(sample, top_gen_plot) %>%
  dplyr::summarize(total_ct = sum(Abundance), .groups = "drop") %>%
  filter(total_ct > 0) %>%
  group_by(sample) %>%
  mutate(total_pct = 100 * total_ct / sum(total_ct)) %>%
  ungroup()

label_order <- sbp_gen_df %>%
  group_by(top_gen_plot) %>%
  dplyr::summarize(total_ct = sum(total_ct), .groups = "drop") %>%
  arrange(total_ct) # Reorder taxa by total abundance

sbp_gen_df$top_gen_plot <- factor(sbp_gen_df$top_gen_plot,
                             levels = label_order$top_gen_plot)

# ================================================================
# -------------------- EXPORT ------------------------------------
# ================================================================

saveRDS(hist_read_df, file.path(DIR_OUT, "descriptive_histogram_reads_fig_df.rds"))
saveRDS(hist_gen_df, file.path(DIR_OUT, "descriptive_histogram_genus_fig_df.rds"))
saveRDS(sbp_gen_df, file.path(DIR_OUT, "descriptive_stacked_barplot_genus_fig_df.rds"))
