# ==============================================================================
# Project: FinlandSports V2.0 
# Script:  01_FigS1_PCA_Limma.R
# Panel:   Fig 1E (Trajectories) & Fig S1 (Changes Proportion & Loadings)
# Features: 离群点检测, 运动前后配对轨迹, Limma差异分析, NPG配色
# ==============================================================================

# [0. 锁定工作目录]
setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

# 1. 环境准备与包加载
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, readxl, dplyr, tidyr, stringr, tibble, limma, ggplot2, ggrepel, openxlsx)

# [铁律 1: 路径标准化]
INPUT_DIR <- "01_Clean_Data"
OUT_FIG   <- "03_Results/Fig_S1"      # 附件图存放地
OUT_STATS <- "03_Results/Summary_Stats" # 统计表存放地
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_STATS, recursive = TRUE, showWarnings = FALSE)

# 2. 严谨的样本解析函数 (对接运动时间点)
parse_samples <- function(sample_ids) {
  df <- data.frame(SampleID = sample_ids, stringsAsFactors = FALSE)
  df$RawTime <- "Unknown"
  df$RawTime[grepl("fast", sample_ids, ignore.case = TRUE)] <- "fast"
  df$RawTime[grepl("pre", sample_ids, ignore.case = TRUE)] <- "pre"
  df$RawTime[grepl("post1h", sample_ids, ignore.case = TRUE)] <- "post1h"
  df$RawTime[grepl("post3h", sample_ids, ignore.case = TRUE)] <- "post3h"
  df$RawTime[grepl("post", sample_ids, ignore.case = TRUE) & df$RawTime == "Unknown"] <- "post3h" 
  
  df$FamilyID <- str_extract(sample_ids, "^[0-9]+")
  raw_role <- str_extract(sample_ids, "(?i)(father|mother|daughter)")
  df$Role <- tools::toTitleCase(tolower(raw_role))
  df$Subject_ID <- paste0(df$FamilyID, "_", df$Role)
  return(df)
}

files <- list.files(INPUT_DIR, pattern = "^Cleaned_.*\\.csv", full.names = TRUE)
prop_list <- list()
loadings_list <- list() 

# ==============================================================================
# 3. 多组学循环处理引擎
# ==============================================================================
for (f in files) {
  dataset_name <- stringr::str_remove(basename(f), "^Cleaned_")
  dataset_name <- stringr::str_remove(dataset_name, "\\.csv$")
  if(grepl("Clinical|Sample", dataset_name, ignore.case=TRUE)) next
  
  message("\n>>> Analyzing Trajectories & Diff: ", dataset_name)
  
  # 读取数据
  df <- read_csv(f, show_col_types = FALSE)
  feat_ids <- df[[1]]; expr_mat <- as.matrix(df[, -1]); rownames(expr_mat) <- feat_ids
  
  clean_sample_names <- gsub("_twin\\.[0-9]+", "", colnames(expr_mat))
  col_meta <- parse_samples(clean_sample_names)
  col_meta$OriginalSampleID <- colnames(expr_mat)
  
  # 定义分析组别 (Baseline vs Exercise)
  base_tp <- if("pre" %in% col_meta$RawTime) "pre" else if("fast" %in% col_meta$RawTime) "fast" else "Unknown"
  exe_tp  <- if("post3h" %in% col_meta$RawTime) "post3h" else "Unknown"
  col_meta$AnalysisGroup <- "Unknown"
  col_meta$AnalysisGroup[col_meta$RawTime == base_tp] <- "Baseline"
  col_meta$AnalysisGroup[col_meta$RawTime == exe_tp] <- "Exercise"
  
  valid_idx <- which(col_meta$AnalysisGroup %in% c("Baseline", "Exercise") & col_meta$Role != "Unknown")
  if(length(valid_idx) < 3) next
  mat_valid <- expr_mat[, valid_idx]; meta_valid <- col_meta[valid_idx, ]
  
  # 4. 离群点检测与 PCA
  pca_res <- prcomp(t(mat_valid), scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x[, 1:2]) %>% mutate(OriginalSampleID = rownames(.)) %>%
    left_join(meta_valid, by="OriginalSampleID")
  
  # 5. 提取 PC Loadings (保存至统计目录)
  loadings <- as.data.frame(pca_res$rotation[, 1:2]) %>% rownames_to_column("Feature")
  top_pc1 <- loadings %>% arrange(desc(abs(PC1))) %>% head(30) %>% mutate(Component="PC1")
  top_pc2 <- loadings %>% arrange(desc(abs(PC2))) %>% head(30) %>% mutate(Component="PC2")
  loadings_list[[dataset_name]] <- bind_rows(top_pc1, top_pc2)
  
  # 6. 绘图：运动轨迹 (Fig 1E)
  role_colors <- c("Father" = "#4DBBD5", "Mother" = "#E64B35", "Daughter" = "#00A087")
  
  pca_paired <- pca_df %>% group_by(Subject_ID) %>% filter(n() == 2) %>% 
    mutate(AnalysisGroup = factor(AnalysisGroup, levels = c("Baseline", "Exercise"))) %>% arrange(Subject_ID, AnalysisGroup)
  
  p_base <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_path(data = pca_paired, aes(group = Subject_ID), color = "gray70", alpha = 0.5, 
              arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
    geom_point(aes(fill = Role, shape = AnalysisGroup), size = 3.5, color = "black", stroke = 0.5) +
    scale_shape_manual(values = c("Baseline" = 21, "Exercise" = 24)) + 
    scale_fill_manual(values = role_colors) +
    theme_bw() + theme(aspect.ratio = 1, plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(title = paste0(dataset_name, "\nIndividual Trajectories"),
         x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 1), "%)"),
         y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 1), "%)")) +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
  
  ggsave(file.path(OUT_FIG, paste0("Fig1E_Trajectory_", dataset_name, ".pdf")), p_base, width = 6, height = 5)
  
  # 7. Limma 差异分析 (Exercise vs Baseline)
  if(length(unique(meta_valid$AnalysisGroup)) == 2) {
    Subject <- factor(meta_valid$Subject_ID)
    Group   <- factor(meta_valid$AnalysisGroup, levels = c("Baseline", "Exercise"))
    design  <- model.matrix(~ Group)
    corfit  <- duplicateCorrelation(mat_valid, design, block = Subject)
    fit     <- lmFit(mat_valid, design, block = Subject, correlation = corfit$consensus)
    fit     <- eBayes(fit)
    res     <- topTable(fit, coef = 2, number = Inf) %>% rownames_to_column("Feature")
    
    write_csv(res, file.path(OUT_STATS, paste0("Limma_Exercise_vs_Baseline_", dataset_name, ".csv")))
    
    # 统计上调/下调比例 (FC > 1.2 即 logFC > 0.263)
    total <- nrow(res)
    n_up  <- sum(res$P.Value < 0.05 & res$logFC > 0.263)
    n_dn  <- sum(res$P.Value < 0.05 & res$logFC < -0.263)
    prop_list[[dataset_name]] <- data.frame(Dataset = dataset_name, Up = n_up, Down = n_dn, Total = total)
  }
}

# 8. 导出汇总结果
if(length(loadings_list) > 0) write.xlsx(loadings_list, file.path(OUT_STATS, "PCA_Top_Loadings.xlsx"))

# 9. 绘制 Proportion of Changes 柱状图 (Fig S1)
if(length(prop_list) > 0) {
  df_bar <- bind_rows(prop_list) %>%
    mutate(Pct_Up = (Up/Total)*100, Pct_Dn = -(Down/Total)*100) %>%
    pivot_longer(cols = starts_with("Pct"), names_to = "Direction", values_to = "Percentage")
  
  p_bar <- ggplot(df_bar, aes(x = reorder(Dataset, Percentage), y = Percentage, fill = Direction)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    coord_flip() + geom_hline(yintercept = 0) +
    scale_fill_manual(values = c("Pct_Dn" = "#4DBBD5", "Pct_Up" = "#E64B35"), labels = c("Down", "Up")) +
    theme_classic() + labs(title = "Proportion of Significant Changes", y = "Percentage of features (P < 0.05)", x = "")
  
  ggsave(file.path(OUT_FIG, "FigS1_Proportion_of_Changes.pdf"), p_bar, width = 8, height = 5)
}

message("\n>>> [Success] Fig 1E & S1 Pipeline Completed.")