# ================================================================
# Run all scripts
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
# -------------------- RUN ENTIRE PIPELINE -----------------------
# ================================================================

# Create log folder

dir.create("Current_logs")

# Logging setup

DIR_LOG <- "~/Desktop/ITS/Current_logs"

log_file <- file.path(DIR_LOG, "log_info.txt")

log_con <- file(log_file, open = "wt")

sink(log_con, split = TRUE)
sink(log_con, type = "message")

on.exit({
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

# Identify scripts

scripts <- list.files(
  path = "./R_scripts/",
  pattern = "^\\d+_.*\\.R$",
  full.names = TRUE
)

scripts <- sort(scripts)

# Exclude this file

scripts <- scripts[!grepl("00_run_all.R", scripts)]

# Run scripts sequentially

start_time <- Sys.time()

for (s in scripts) {
  cat("\n====================================\n")
  cat("Running:", s, "\n")
  cat("====================================\n")
  source(s, local = new.env())
}

end_time <- Sys.time()

# Total runtime

cat("\nTotal runtime:\n")
print(end_time - start_time)

# ================================================================
# -------------------- SAVE SESSION INFORMATION ------------------
# ================================================================

cat("Saving session information...\n")
writeLines(capture.output(sessionInfo()), file.path(DIR_LOG, "session_info.txt"))
cat("All analyses completed successfully.\n")
