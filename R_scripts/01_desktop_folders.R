# ================================================================
# Create project folder structure in desktop
# ================================================================

# ================================================================
# -------------------- SET WORKING DIRECTORY ---------------------
# ================================================================

DIR_BASE <- "~/Desktop/ITS"

# ================================================================
# -------------------- CREATE FOLDERS ----------------------------
# ================================================================

# Create main folder

if (!dir.exists(DIR_BASE)) {
  dir.create(DIR_BASE, recursive = TRUE)
}

# Create subfolders

subfolders <- c(
  #"Datasets",
  "Figures",
  "Output",
  "Current_logs",
  "Tables"
)

for (folder in subfolders) {
  dir.create(file.path(DIR_BASE, folder), showWarnings = FALSE)
}

# Confirmation message

cat("Folder structure created successfully on desktop.\n")
