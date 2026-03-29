# ==============================================================================
# Project: FinlandSports V2.0 
# Script:  01_Fig1D_CV_Analysis.R
# Description: 计算基线变异度，包含总体(All)、肥胖(Obese)与消瘦(Lean)三个独立分组
# Features: V2.0 相对路径体系, 对接严格临床表, 纯净 Fat(%) 变量驱动
# Core: 极速版栅格化 (Unified Rasterization) + 底部高级横线分组标注
# ==============================================================================

# [0. 核心纠正：锁定战区大本营坐标 (Project Root)]
# 这句话是全篇唯一的绝对路径锚点，执行后，后续所有代码都将完美使用相对路径！
setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

# 1. 环境准备与包加载
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, readxl, dplyr, tidyr, stringr, tibble, ggplot2, matrixStats, ggrastr)

# [铁律 1: 从大本营出发的纯粹相对路径]
CONF_CLEAN_DIR <- "01_Clean_Data"
CONF_CV_DIR    <- "03_Results/Fig_1/Panel_D_CV"
dir.create(CONF_CV_DIR, recursive = TRUE, showWarnings = FALSE)

# 2. 样本解析与临床数据合并 (带智能角色映射引擎)
parse_samples <- function(sample_ids) {
  df <- data.frame(SampleID = sample_ids, stringsAsFactors = FALSE)
  df$RawTime <- "Unknown"
  df$RawTime[grepl("fast", sample_ids, ignore.case = TRUE)] <- "fast"
  df$RawTime[grepl("pre", sample_ids, ignore.case = TRUE)] <- "pre"
  
  df$FamilyID <- str_extract(sample_ids, "^[0-9]+")
  raw_role <- str_extract(sample_ids, "(?i)(father|mother|daughter)")
  df$Role <- tools::toTitleCase(tolower(raw_role))
  
  # 智能映射组学文本到临床数字主键
  df$Membercode <- case_when(
    df$Role == "Daughter" ~ 1,
    df$Role == "Mother"   ~ 2,
    df$Role == "Father"   ~ 3,
    TRUE ~ NA_real_
  )
  df$Clinical_Subject_ID <- paste0(df$FamilyID, "_", df$Membercode)
  return(df)
}

# [铁律 2: 挂载 V2.0 纯净临床表与 Fat(%) 断言]
clin_file <- file.path(CONF_CLEAN_DIR, "Clinical_Master_Strict.csv")
if(!file.exists(clin_file)) {
  stop("FATAL ERROR: 找不到临床主表！请确保 Clinical_Master_Strict.csv 位于 01_Clean_Data 中！")
}

clin_df <- read_csv(clin_file, show_col_types = FALSE) %>%
  select(Clinical_Subject_ID, `Fat(%)`)

# 定义全局高级组学配色
omics_colors <- c(
  "Proteomics"   = "#9A7EAF", # 丁香紫 (Lilac Purple)
  "Metabonomics" = "#D96A70", # 蔷薇红 (Muted Rose)
  "Methylation"  = "#5C7699", # 钢青蓝 (Steel Blue)
  "Microarray"   = "#8B9D76"  # 鼠尾草绿 (Sage Green)
)

dataset_order <- c("Serum_Proteomics", "Serum_Metabonomics", 
                   "Adipose_Microarray", "Adipose_Methylation", "Adipose_Proteomics", 
                   "Muscle_Microarray", "Muscle_Methylation", "Muscle_Proteomics")

files <- list.files(CONF_CLEAN_DIR, pattern = "^Cleaned_.*\\.csv", full.names = TRUE)
cv_base_list <- list()

# ==============================================================================
# 3. 数据循环与 CV 计算 (带进度条)
# ==============================================================================
message("\n>>> [1/2] 正在加载多组学数据并计算 CV (共 ", length(files), " 个数据集)...")
pb_data <- txtProgressBar(min = 0, max = length(files), style = 3)

for (i in seq_along(files)) {
  f <- files[i]
  
  dataset_name <- stringr::str_remove(basename(f), "^Cleaned_")
  dataset_name <- stringr::str_remove(dataset_name, "\\.csv$")
  if(grepl("Clinical|Sample", dataset_name, ignore.case=TRUE)) next
  
  parts <- str_split(dataset_name, "_")[[1]]
  formatted_label <- paste0(parts[1], " (", tolower(parts[2]), ")")
  
  df <- read_csv(f, show_col_types = FALSE)
  feat_ids <- df[[1]]; expr_mat <- as.matrix(df[, -1]); rownames(expr_mat) <- feat_ids
  
  clean_sample_names <- gsub("_twin\\.[0-9]+", "", colnames(expr_mat))
  col_meta <- parse_samples(clean_sample_names)
  col_meta$OriginalSampleID <- colnames(expr_mat)
  
  # 锚定 Baseline
  base_tp <- if("pre" %in% col_meta$RawTime) "pre" else if("fast" %in% col_meta$RawTime) "fast" else "Unknown"
  col_meta$AnalysisGroup <- ifelse(col_meta$RawTime == base_tp, "Baseline", "Unknown")
  
  # 连接体脂数据并划定胖瘦组 (完全依赖 Fat(%) 变量)
  col_meta <- left_join(col_meta, clin_df, by = "Clinical_Subject_ID")
  col_meta <- col_meta %>% 
    mutate(Phenotype = case_when(
      `Fat(%)` > 30 ~ "Obese",
      `Fat(%)` <= 30 ~ "Lean",
      TRUE ~ NA_character_
    ))
  
  valid_idx <- which(col_meta$AnalysisGroup == "Baseline" & col_meta$Role != "Unknown" & !is.na(col_meta$Phenotype))
  if(length(valid_idx) < 3) next
  mat_valid <- expr_mat[, valid_idx]; meta_valid <- col_meta[valid_idx, ]
  
  # 剔除离群点并逆向对数还原
  pca_res <- prcomp(t(mat_valid), scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x[, 1:2]) %>% mutate(OriginalSampleID = rownames(.))
  pca_df <- left_join(pca_df, meta_valid, by="OriginalSampleID") %>% 
    mutate(Z_PC1 = scale(PC1), Z_PC2 = scale(PC2), is_far = abs(Z_PC1) > 3 | abs(Z_PC2) > 3)
  to_remove <- pca_df$OriginalSampleID[pca_df$is_far & ave(pca_df$OriginalSampleID, pca_df$FamilyID, FUN=length) <= 1]
  
  if(length(to_remove) > 0) {
    mat_valid <- mat_valid[, !colnames(mat_valid) %in% to_remove]
    meta_valid <- meta_valid[!meta_valid$OriginalSampleID %in% to_remove, ]
  }
  
  if(!grepl("Methylation", dataset_name, ignore.case = TRUE) && max(mat_valid, na.rm = TRUE) < 50) mat_valid <- 2^mat_valid
  non_zero_vals <- mat_valid[mat_valid > 1e-6 & !is.na(mat_valid)]
  if(length(non_zero_vals) > 0) mat_valid[mat_valid <= 1e-6 | mat_valid <= min(non_zero_vals) * 1.001] <- NA
  
  # 分别针对 All (总体), Lean 和 Obese 计算 CV
  for(pheno in c("All", "Lean", "Obese")) {
    
    if (pheno == "All") {
      base_cols <- meta_valid$OriginalSampleID[!is.na(meta_valid$Phenotype)]
    } else {
      base_cols <- meta_valid$OriginalSampleID[meta_valid$Phenotype == pheno]
    }
    
    if(length(base_cols) < 3) next
    
    mat_base <- mat_valid[, base_cols, drop=FALSE]
    valid_features <- rowSums(!is.na(mat_base)) >= 3
    base_mean <- rowMeans(mat_base, na.rm = TRUE)
    base_sd <- rowSds(mat_base, na.rm = TRUE)
    keep_mask <- valid_features & base_mean > 1e-4
    
    cv_base <- (base_sd[keep_mask] / base_mean[keep_mask]) * 100
    if(length(cv_base) > 0) {
      cv_base_list[[paste0(dataset_name, "_", pheno)]] <- data.frame(
        Dataset = dataset_name, Omics = parts[2], Label = formatted_label,
        Phenotype = pheno, Feature = names(cv_base), CV_Percent = cv_base
      )
    }
  }
  setTxtProgressBar(pb_data, i)
}
close(pb_data)

# ==============================================================================
# 4. 汇总统计、终极绘图与结果导出
# ==============================================================================

if(length(cv_base_list) > 0) {
  cv_df <- bind_rows(cv_base_list)
  
  # 1. 导出汇总统计表格
  cv_summary <- cv_df %>%
    group_by(Dataset, Omics, Label, Phenotype) %>%
    summarise(
      Feature_Count = n(),
      Mean_CV = mean(CV_Percent, na.rm = TRUE),
      Median_CV = median(CV_Percent, na.rm = TRUE),
      SD_CV = sd(CV_Percent, na.rm = TRUE),
      Min_CV = min(CV_Percent, na.rm = TRUE),
      Max_CV = max(CV_Percent, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(factor(Phenotype, levels = c("All", "Obese", "Lean")), Label)
  
  write_csv(cv_summary, file.path(CONF_CV_DIR, "CV_Summary_Stats_Baseline_All_Obese_Lean.csv"))
  write_csv(cv_df, file.path(CONF_CV_DIR, "CV_Data_Baseline_All_Obese_Lean_Raw.csv"))
  
  # 2. 数据平滑与因子排序
  cv_df <- cv_df %>% group_by(Label, Phenotype) %>% filter(CV_Percent <= quantile(CV_Percent, 0.99)) %>% ungroup()
  label_order <- sapply(dataset_order, function(x) paste0(str_split(x, "_")[[1]][1], " (", tolower(str_split(x, "_")[[1]][2]), ")"))
  cv_df$Label <- factor(cv_df$Label, levels = label_order[label_order %in% unique(cv_df$Label)])
  cv_df$Phenotype <- factor(cv_df$Phenotype, levels = c("All", "Obese", "Lean"))
  cv_df$Omics <- factor(cv_df$Omics, levels = names(omics_colors))
  
  # 3. 独立绘图并分别保存 (分层极速渲染版 + 严格坐标与顺序锁定)
  pheno_groups <- c("All", "Obese", "Lean")
  
  message("\n>>> [2/2] 开始渲染全量栅格化高清多组学散点图 (极速渲染引擎启动)...")
  pb_plot <- txtProgressBar(min = 0, max = length(pheno_groups), style = 3)
  
  for (i in seq_along(pheno_groups)) {
    p_type <- pheno_groups[i]
    
    plot_data <- cv_df %>% filter(Phenotype == p_type & CV_Percent <= 200)
    if(nrow(plot_data) == 0) next
    
    # 拆分数据
    plot_data_serum  <- plot_data %>% filter(grepl("Serum", Label))
    plot_data_tissue <- plot_data %>% filter(!grepl("Serum", Label))
    
    plot_title <- paste0("Baseline Inter-individual Variability (", p_type, ")")
    
    p <- ggplot(plot_data, aes(x = Label, y = CV_Percent, color = Omics)) +
      
      # 强制锁定 X 轴顺序
      scale_x_discrete(limits = levels(cv_df$Label)) +
      
      # 图层 1：非血清组织 
      geom_jitter_rast(
        data = plot_data_tissue,
        size = 0.5, alpha = 0.25, 
        raster.dpi = 300,
        width = 0.25, 
        stroke = 0
      ) +
      
      # 图层 2：血清组学 
      geom_jitter_rast(
        data = plot_data_serum,
        size = 1.2, alpha = 0.85, 
        raster.dpi = 300,
        width = 0.25, 
        stroke = 0
      ) +
      
      # 顶层：空心箱线图
      geom_boxplot(fill = NA, width = 0.6, linewidth = 0.8, outlier.shape = NA) +
      
      # 映射颜色方案
      scale_color_manual(values = omics_colors) +
      scale_y_continuous(breaks = seq(0, 200, by = 25)) +
      coord_cartesian(ylim = c(0, 200), clip = "off") +
      
      # 底部高级分组横线标注
      annotate("segment", x = 0.6, xend = 2.4, y = -12, yend = -12, color = "black", linewidth = 1) +
      annotate("text", x = 1.5, y = -13, label = "Serum", fontface = "bold", size = 5, vjust = 1) +
      annotate("segment", x = 2.6, xend = 5.4, y = -12, yend = -12, color = "black", linewidth = 1) +
      annotate("text", x = 4.0, y = -13, label = "Adipose", fontface = "bold", size = 5, vjust = 1) +
      annotate("segment", x = 5.6, xend = 8.4, y = -12, yend = -12, color = "black", linewidth = 1) +
      annotate("text", x = 7.0, y = -13, label = "Muscle", fontface = "bold", size = 5, vjust = 1) +
      
      theme_classic(base_size = 14) +
      labs(title = plot_title, y = "Coefficient of Variation (CV %)", x = "", color = "Omics Type") +
      guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) + 
      theme(
        aspect.ratio = 1, 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 11),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 35, l = 10), 
        axis.text.y = element_text(face = "bold", color="black"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
      )
    
    file_name <- paste0("FigS1_CV_Baseline_", p_type, "_ScatterBox.pdf")
    ggsave(file.path(CONF_CV_DIR, file_name), p, width = 9, height = 7.5)
    
    setTxtProgressBar(pb_plot, i)
  }
  close(pb_plot)
  message("\n>>> 🎉 所有分析与全量极速光栅化图表生成已完美结束！结果已送达 03_Results/Fig_1/Panel_D_CV")
}