# ==============================================================================
# Project: FinlandSports V2.0 
# Script:  02_Fig1D_PCA_Distance.R
# Panel:   Fig 1D (Middle) - PCA Quadrant & Familial Euclidean Distance
# Core:    1:1 Aspect Ratio, Clean Data Only, Auto-Parsing, Fat(%) Assertion.
# ==============================================================================

# [0. Lock working directory]
# setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

# [1. Environment setup and package loading]
set.seed(42)

if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, dplyr, stringr, ggplot2, ggpubr, scales, tibble, tidyr, grid, cowplot, openxlsx)

# Define paths (All point to Clean Data, completely deprecating Raw Data's Sample Sheet)
INPUT_DIR <- "01_Clean_Data"
OUT_DIR   <- "03_Results/Fig_1/Panel_D_PCA"
if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Classic layout theme
my_theme <- theme_classic() +
  theme(
    aspect.ratio = 1, 
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1.2),
    axis.line = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 12, color = "black", face="bold"), 
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none" 
  )

# [2. Automated sample parsing engine (Completely replaces old Sample Sheet)]
parse_samples <- function(sample_ids) {
  df <- data.frame(OriginalSampleID = sample_ids, stringsAsFactors = FALSE)
  df$TimePoint_Norm <- "unknown"
  df$TimePoint_Norm[grepl("fast", sample_ids, ignore.case = TRUE)] <- "fast"
  df$TimePoint_Norm[grepl("pre", sample_ids, ignore.case = TRUE)] <- "pre"
  
  df$FamilyID <- str_extract(sample_ids, "^[0-9]+")
  raw_role <- str_extract(sample_ids, "(?i)(father|mother|daughter)")
  df$Role <- tools::toTitleCase(tolower(raw_role))
  
  df$Membercode <- case_when(
    df$Role == "Daughter" ~ 1,
    df$Role == "Mother"   ~ 2,
    df$Role == "Father"   ~ 3,
    TRUE ~ NA_real_
  )
  df$Clinical_Subject_ID <- paste0(df$FamilyID, "_", df$Membercode)
  return(df)
}

# [3. Load pure clinical master table and extract Fat(%) grouping]
cat("Step 1: Loading clinical master...\n")
clin_path <- file.path(INPUT_DIR, "Clinical_Master_Strict.csv")
if(!file.exists(clin_path)) stop(paste("FATAL ERROR: Cannot find clinical table! Path:", clin_path))

clin_df <- read_csv(clin_path, show_col_types = FALSE) %>% 
  select(Clinical_Subject_ID, `Fat(%)`) %>%
  mutate(Status = ifelse(`Fat(%)` > 30, "Obese", "Lean"))

# [4. Core distance algorithm]
calc_pca_distance <- function(pc_matrix, labels, is_real = TRUE, n_samples = 2000) {
  if (is_real) {
    dists <- c()
    for (fam in unique(labels)) {
      idx <- which(labels == fam)
      if (length(idx) >= 2) {
        d_mat <- dist(pc_matrix[idx, , drop=FALSE], method = "euclidean")
        dists <- c(dists, as.numeric(d_mat))
      }
    }
    return(dists)
  } else {
    dists <- c()
    all_idx <- 1:nrow(pc_matrix)
    count <- 0
    while (count < n_samples) {
      pair <- sample(all_idx, 2)
      if (labels[pair[1]] != labels[pair[2]]) {
        d <- dist(pc_matrix[pair, , drop=FALSE])
        dists <- c(dists, as.numeric(d))
        count <- count + 1
      }
    }
    return(dists)
  }
}

# [5. Execute omics loop calculation]
target_datasets <- c("Serum_Proteomics", "Adipose_Proteomics", "Muscle_Proteomics", "Serum_Metabolomics", "Serum_Metabonomics") 
boxplot_targets <- c("Serum_Proteomics", "Adipose_Proteomics", "Muscle_Proteomics")

wb <- createWorkbook()

distinct_colors <- c("#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231", 
                     "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE", 
                     "#008080", "#E6BEFF", "#800000", "#AAFFC3", "#808000")
fam_color_map <- NULL 

for (ds_name in target_datasets) {
  f_path <- file.path(INPUT_DIR, paste0("Cleaned_", ds_name, ".csv"))
  if(!file.exists(f_path)) next
  
  cat(paste0("\n>>> Processing dataset: ", ds_name, "...\n"))
  
  df <- read_csv(f_path, show_col_types = FALSE)
  feat_ids <- df[[1]]; expr_mat <- as.matrix(df[, -1]); rownames(expr_mat) <- feat_ids
  
  # Use parsing engine to process column names and get mapping info directly
  clean_sample_names <- gsub("_twin\\.[0-9]+", "", colnames(expr_mat))
  col_meta <- parse_samples(clean_sample_names)
  col_meta$OriginalSampleID <- colnames(expr_mat)
  
  # Merge clinical data (directly dock with Clinical_Master_Strict.csv)
  col_meta <- left_join(col_meta, clin_df, by = "Clinical_Subject_ID")
  
  # Intercept Baseline data and remove samples with missing phenotypes or roles
  valid_tps <- c("pre", "fasting", "fast")
  meta_sub <- col_meta %>% 
    filter(TimePoint_Norm %in% valid_tps & !is.na(Status) & !is.na(Role))
  
  common_samples <- meta_sub$OriginalSampleID
  
  if(length(common_samples) < 10) next
  
  t_mat <- t(expr_mat[, common_samples])
  
  if(is.null(fam_color_map)) {
    unique_fams <- sort(unique(meta_sub$FamilyID))
    fam_color_map <- setNames(colorRampPalette(distinct_colors)(length(unique_fams)), unique_fams)
  }
  
  pca_res <- prcomp(t_mat, center = TRUE, scale. = TRUE)
  var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)
  n_pcs <- min(20, ncol(pca_res$x))
  pc_for_dist <- pca_res$x[, 1:n_pcs]
  
  plot_df <- meta_sub %>% mutate(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2])
  addWorksheet(wb, paste0("PCA_", substring(ds_name, 1, 20)))
  writeData(wb, paste0("PCA_", substring(ds_name, 1, 20)), plot_df)
  
  p_pca <- ggplot() +
    geom_point(data = plot_df %>% filter(Status == "Lean"),
               aes(x=PC1, y=PC2, color=FamilyID, shape=Role),
               fill="white", size=4.5, stroke=1.2) +
    geom_point(data = plot_df %>% filter(Status == "Obese"),
               aes(x=PC1, y=PC2, fill=FamilyID, shape=Role),
               color="black", size=4.5, stroke=0.6) +
    scale_shape_manual(values = c("Daughter"=24, "Mother"=21, "Father"=22)) +
    scale_color_manual(values = fam_color_map) +
    scale_fill_manual(values = fam_color_map) +
    labs(title = paste0(ds_name, "\nPCA Quadrant Plot"), x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)")) +
    my_theme
  
  pca_filename <- paste0("Fig-1D-PCA-", gsub("_", "-", ds_name), ".pdf")
  ggsave(file.path(OUT_DIR, pca_filename), p_pca, width=5.5, height=5.5)
  cat(paste0("    [Saved] ", pca_filename, "\n"))
  
  if (ds_name %in% boxplot_targets) {
    real_d <- calc_pca_distance(pc_for_dist, plot_df$FamilyID, is_real = TRUE)
    pseudo_d <- calc_pca_distance(pc_for_dist, plot_df$FamilyID, is_real = FALSE, n_samples = 2000)
    
    df_dist <- data.frame(
      Type = factor(c(rep("Real Families", length(real_d)), rep("Pseudo Families", length(pseudo_d))), 
                    levels = c("Real Families", "Pseudo Families")),
      Distance = c(real_d, pseudo_d)
    )
    
    addWorksheet(wb, paste0("Dist_", substring(ds_name, 1, 20)))
    writeData(wb, paste0("Dist_", substring(ds_name, 1, 20)), df_dist)
    
    p_box <- ggplot(df_dist, aes(x=Type, y=Distance, color=Type)) +
      stat_boxplot(geom="errorbar", width=0.1, linewidth=0.8) + 
      geom_boxplot(width=0.25, linewidth=0.8, fill=NA, outlier.shape=NA) + 
      geom_jitter(width=0.1, size=1.5, alpha=0.1) +
      stat_compare_means(comparisons = list(c("Real Families", "Pseudo Families")), 
                         method = "wilcox.test", method.args = list(alternative = "less"),
                         symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
                         size = 6, vjust = 0.5) +
      scale_color_manual(values = c("Real Families"="#E76F51", "Pseudo Families"="#2A9D8F")) +
      labs(title = paste0(ds_name, "\nFamilial Similarity"), y = "Pairwise Euclidean Distance", x = NULL) +
      my_theme
    
    box_filename <- paste0("Fig-1D-Dist-", gsub("_", "-", ds_name), ".pdf")
    ggsave(file.path(OUT_DIR, box_filename), p_box, width=4.5, height=5.5)
    cat(paste0("    [Saved] ", box_filename, "\n"))
  }
}

# [6. Generate global common legend separately]
cat("\nStep 4: Extracting standard legend...\n")
dummy_df <- data.frame(X=1:2, Y=1:2, Role=c("Daughter","Mother"), Status=c("Lean","Obese"))
p_legend <- ggplot(dummy_df, aes(x=X, y=Y)) +
  geom_point(aes(shape=Role, fill=Status, color=Status), size=5, stroke=1) +
  scale_shape_manual(name = "Family Member", values = c("Daughter"=24, "Mother"=21, "Father"=22)) +
  scale_fill_manual(name = "Fat(%) Status", values = c("Lean"="white", "Obese"="grey40")) +
  scale_color_manual(name = "Fat(%) Status", values = c("Lean"="grey40", "Obese"="black")) +
  theme_void() + 
  theme(legend.position = "right", 
        legend.title = element_text(face="bold", size=12), 
        legend.text = element_text(size=11),
        legend.box = "horizontal")

legend_grob <- get_legend(p_legend)
ggsave(file.path(OUT_DIR, "Fig-1D-Common-Legend.pdf"), as_ggplot(legend_grob), width=6, height=3)

saveWorkbook(wb, file.path(OUT_DIR, "Fig-1D-SourceData.xlsx"), overwrite = TRUE)

cat("\n>>> Analysis complete! Files delivered to 03_Results/Fig_1/Panel_D_PCA!\n")
