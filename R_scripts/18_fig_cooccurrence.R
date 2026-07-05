# ================================================================
# Figure - Heatmap cooccurrence
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

load_pkg("CooccurrenceAffinity")

# ================================================================
# -------------------- LOAD DATASETS -----------------------------
# ================================================================

res <- readRDS(file.path(DIR_OUT, "heatmap_cooccurrence_fig_df.rds"))

# ================================================================
# -------------------- CREATING HEATMAP --------------------------
# ================================================================

fig <- plotgg(
  data = res,
  variable = "alpha_mle",
  legendlimit = "balanced"
) +
  labs(fill = "Co-occurrence affinity (α̂)") +
  theme(
    axis.text.x = element_text(size = 10, color = "black", angle = 35, hjust = 1),
    axis.text.y = element_text(size = 10, color = "black")
  )

fig

# ================================================================
# -------------------- EXPORT ------------------------------------
# ================================================================

# SVG

svglite(file.path(DIR_FIG, "Fig_heatmap_cooccurrence.svg"), width = 7, height = 5)
print(fig)
dev.off()

# PDF

pdf(file.path(DIR_FIG, "Fig_heatmap_cooccurrence.pdf"), width = 7, height = 5)
print(fig)
dev.off() # Note that there is an issue with labels in PDF

# Further annotation in Inkscape
