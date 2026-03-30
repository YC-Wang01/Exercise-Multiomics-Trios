# ==============================================================================
# Project: FinlandSports (Multi-Omics Bifurcation)
# Script:  01_Omics_QC_and_Normalization.R
# Core:    Omics-specific normalization (Perseus-style MS imputation vs Quantile)
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(stringr)
  library(limma)
})

input_dir <- "00_Raw_Data"
output_dir <- "01_Clean_Data"
qc_dir <- "01_Clean_Data/QC_Metrics"
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Core Processing Engine
# ==============================================================================
process_omics <- function(file_name, type, sheet = 1) {
  message("\n>>> Processing: ", file_name, " [", type, "]")
  full_path <- file.path(input_dir, file_name)
  if (!file.exists(full_path)) {
    warning("  [Skip] File not found: ", file_name)
    return(NULL)
  }
  
  # 1. Robust loading and Feature ID locking
  df <- tryCatch({ read_excel(full_path, sheet = sheet, na = c("", "NA", "NaN", "0")) },
                 error = function(e) { read_excel(full_path, sheet = 1, na = c("", "NA", "NaN", "0")) })
  
  colnames(df)[1] <- "FeatureID"
  df <- df %>% filter(!is.na(FeatureID) & str_trim(FeatureID) != "") %>% distinct(FeatureID, .keep_all = TRUE)
  
  gene_ids <- as.character(df[[1]])
  expr <- df[, -1] %>% mutate(across(everything(), as.numeric))
  colnames(expr) <- tolower(sub("\\.\\.\\.[0-9]+$", "", colnames(expr)))
  
  mat <- as.matrix(expr)
  rownames(mat) <- gene_ids
  
  if(any(duplicated(colnames(mat)))) {
    colnames(mat) <- make.unique(colnames(mat), sep = "_twin.")
  }
  
  # 2. Contaminant Removal (KRT for Proteomics)
  if (type %in% c("Tissue_Prot", "Serum_Prot")) {
    is_krt <- grepl("^KRT[0-9]+", rownames(mat), ignore.case = TRUE)
    if (sum(is_krt) > 0) {
      mat <- mat[!is_krt, ]
      message("  [-] Removed ", sum(is_krt), " Keratin contaminants.")
    }
  }
  
  # 3. NA filtering and algorithmic bifurcation
  if (type %in% c("Tissue_Prot", "Serum_Prot", "Metabolomics")) {
    # Mass Spec: Allow minor NAs, remove highly missing features (>20% NA)
    keep_na <- rowSums(!is.na(mat)) >= (ncol(mat) * 0.8)
    mat <- mat[keep_na, , drop = FALSE]
    
    # Log2 Transformation
    mat <- log2(mat)
    mat[is.infinite(mat)] <- NA
    
    # Perseus-style Imputation (down-shift to simulate below limit of detection)
    if(sum(is.na(mat)) > 0) {
      valid_vals <- mat[!is.na(mat)]
      fill_mu <- mean(valid_vals) - 1.8 * sd(valid_vals)
      fill_sigma <- 0.3 * sd(valid_vals)
      na_idx <- which(is.na(mat))
      set.seed(123) 
      mat[na_idx] <- rnorm(length(na_idx), mean = fill_mu, sd = fill_sigma)
      message("  [+] Executed Perseus-style down-shift imputation.")
    }
    
    # Median Sweep Normalization
    sm <- apply(mat, 2, median, na.rm=TRUE)
    mat <- sweep(mat, 2, median(sm) - sm, "+")
    
  } else {
    # Microarray/Methylation: Quantile Normalization
    rnames <- rownames(mat); cnames <- colnames(mat)
    mat <- limma::normalizeQuantiles(mat)
    rownames(mat) <- rnames; colnames(mat) <- cnames
    message("  [+] Executed limma::normalizeQuantiles.")
  }
  
  # 4. Export purified expression matrix
  clean_name <- tools::file_path_sans_ext(file_name)
  df_out <- as.data.frame(mat) %>% tibble::rownames_to_column("FeatureID")
  write_csv(df_out, file.path(output_dir, paste0("Cleaned_", clean_name, ".csv")))
  message("  [v] Exported cleaned matrix: ", ncol(mat), " samples, ", nrow(mat), " features.")
}

# ==============================================================================
# Batch Execution Array
# ==============================================================================
message("\n>>> INITIATING OMICS NORMALIZATION PIPELINE...")

# Adipose Tissue
process_omics("Adipose_Microarray.xlsx",  type = "Microarray")
process_omics("Adipose_Proteomics.xlsx",  type = "Tissue_Prot")
process_omics("Adipose_Methylation.xlsx", type = "Methylation")

# Muscle Tissue
process_omics("Muscle_Microarray.xlsx",   type = "Microarray")
process_omics("Muscle_Proteomics.xlsx",   type = "Tissue_Prot")
process_omics("Muscle_Methylation.xlsx",  type = "Methylation")

# Serum Omics
process_omics("Serum_Proteomics.xlsx",    type = "Serum_Prot")
process_omics("Serum_Metabonomics.xlsx",  type = "Metabolomics") 
# Auto fallback for Sheet1/2 is handled within the core engine
message("\n>>> ALL OMICS MATRICES PURIFIED AND NORMALIZED.")