# ================================================================
# Include global variable labels
# ================================================================

# ================================================================
# -------------------- GLOBAL VARIABLE LABELS --------------------
# ================================================================

label_list <- list(
  enrollageinweeks_dev = "Age at enrollment (weeks)",
  efgestagewk = "Gestational age (weeks)",
  birthWgtGm_dev = "Birth weight (grams)",
  age_sample = "Age at sample collection (years)",
  efsex = "Male sex",
  race_white_dev = "Non-Hispanic White",
  efdelivery = "Birth by C-section",
  any_breastfd_dev = "Any breastfeeding",
  mat_asthma_dev = "Maternal asthma",
  young_sib_dev = "Presence of household siblings ages 6 years or younger at birth",
  early_life_anti_dev = "Exposure to antibiotics in utero or during early infancy",
  early_life_smoke_dev = "Exposure to tobacco smoke in utero or during early infancy",
  insurance_dev = "Private insurance during early infancy",
  petany_dev = "Pet ownership during early infancy",
  efrural_dev = "Rural residence during early infancy",
  any_lrti_dev = "History of lower respiratory tract infection in the first year of life",
  mold_dev = "Residential mold damage during early infancy",
  season_collection = "Season of sample collection"
)

# ================================================================
# -------------------- LOADING FUNCTION --------------------------
# ================================================================

apply_labels <- function(df) {
  for (var in names(label_list)) {
    if (var %in% names(df)) {
      Hmisc::label(df[[var]]) <- label_list[[var]]
    }
    }
  return(df)
  }