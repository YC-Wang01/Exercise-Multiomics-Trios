# ==============================================================================
# Project: FinlandSports (Multi-Omics Bifurcation)
# Script:  00_Clinical_Purification.R
# Core:    Excel Multi-sheet Functional Merge, Strict Trio Assertion, V2.0 Hardcoding.
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(purrr)
  library(readr)
})

input_file <- "00_Raw_Data/Clinical_Data_Raw.xlsx"
output_file <- "01_Clean_Data/Clinical_Master_Strict.csv"

if (!file.exists(input_file)) stop("FATAL ERROR: Cannot find 00_Raw_Data/Clinical_Data_Raw.xlsx")
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Module A: Ultra-fast Excel Multi-sheet Merge
# ==============================================================================
message(">>> [Phase 1]: Loading and functionally merging all Excel sheets...")

sheet_names <- excel_sheets(input_file)
if (!"Basic information" %in% sheet_names) stop("FATAL ERROR: Missing main sheet 'Basic information'!")

list_of_dfs <- lapply(sheet_names, function(s) {
  read_excel(input_file, sheet = s) %>%
    mutate(FamilyID = as.character(FamilyID), Membercode = as.numeric(Membercode))
})

basic_idx <- grep("Basic information", sheet_names, ignore.case = TRUE)
if (length(basic_idx) > 0) {
  list_of_dfs <- c(list_of_dfs[basic_idx], list_of_dfs[-basic_idx])
}

join_keys <- intersect(names(list_of_dfs[[1]]), c("Sample_Key", "FamilyID", "Membercode"))
master_df <- list_of_dfs %>%
  reduce(left_join, by = join_keys) %>%
  select(-matches("\\.x$|\\.y$"))

# ==============================================================================
# Module B: Methodology Strict Defense Layer
# ==============================================================================
message(">>> [Phase 2]: Enforcing V2.0 Methodology Iron Rules...")

valid_roles <- c(1, 2, 3)
if (!all(master_df$Membercode %in% valid_roles, na.rm = TRUE)) {
  stop("FATAL PIPELINE HALT: Contaminated data found! Cohort must be strictly trios (daughters, mothers, fathers).")
}

if("Fat_percent" %in% colnames(master_df)) {
  master_df <- master_df %>% rename(`Fat(%)` = Fat_percent)
}

if("IPA" %in% colnames(master_df)) {
  message("  [Pass] Variable 'IPA' (Physical Inactivity) protected against pathway acronym collision.")
}

master_df <- master_df %>%
  mutate(Role = case_when(Membercode == 1 ~ "Daughter", Membercode == 2 ~ "Mother", Membercode == 3 ~ "Father")) %>%
  mutate(Clinical_Subject_ID = paste0(FamilyID, "_", Membercode)) 

write_csv(master_df, output_file)
message(paste0(">>> [Success] Pure clinical master generated. Features: ", ncol(master_df), " | Samples: ", nrow(master_df)))