# ==============================================================================
# Project: FinlandSports V2.0
# Script:  03_FigS2_Clinical_GSEA_Crosstalk.R
# Author:  Sorcier_W (Assisted by AI)
# Description: 
#   1. Limma Analysis (Blood Corrected + Time Filtered)
#   2. GSEA Dotplots (Classic Bubble Plot)
#   3. Crosstalk Plots (4 Key Combinations)
# Features: V2.0 I/O Standardization, Column-Parsing Engine, Pure White Theme.
# ==============================================================================

# [0] 环境初始化
setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  limma, readr, readxl, dplyr, stringr, 
  clusterProfiler, org.Hs.eg.db, ggplot2, 
  ggrepel, tibble, grid, openxlsx,
  enrichplot, cowplot
)

options(timeout = 600)

# ==============================================================================
# [1] 用户配置区域 (User Configuration)
# ==============================================================================
DIR_DATA  <- "01_Clean_Data"
DIR_OUT   <- "03_Results/Fig_S2/Clinical_GSEA"

TIME_FILTER_KEYWORD <- "pre" 
TARGET_VARS <- c("BMI", "HOMA_IR", "S_TotalOC", "CRP") 

KEYWORD_DICT <- list(
  "BMI" = c("Lipid", "Fatty acid", "Adipocytokine", "PPAR", "Insulin", "Thermogenesis", "Mitochondria", "Oxidative"),
  "HOMA_IR" = c("Insulin", "Glycolysis", "Gluconeogenesis", "Diabetes", "FoxO", "AMPK", "mTOR", "PI3K"),
  "S_TotalOC" = c("Osteoclast", "Wnt", "Insulin", "Adipocytokine", "Calcium", "Thyroid", "Parathyroid", "Signaling"),
  "CRP" = c("TNF", "NF-kappa", "Inflammation", "Cytokine", "Chemokine", "Toll-like", "IL-17", "NOD-like", "Immune"),
  "DEFAULT" = c("Metabolic", "Signaling", "Immune")
)

# ==============================================================================
# [2] 核心函数定义
# ==============================================================================

# --- 2.1 血液评分计算 (Blood Score) ---
calculate_blood_score <- function(expr_matrix) {
  blood_genes <- c("HBB", "HBA1", "HBA2") 
  found_genes <- intersect(rownames(expr_matrix), blood_genes)
  if (length(found_genes) == 0) {
    return(rep(0, ncol(expr_matrix)))
  } else {
    blood_expr <- expr_matrix[found_genes, , drop = FALSE]
    raw_score <- if (nrow(blood_expr) > 1) colSums(blood_expr, na.rm = T) else as.numeric(blood_expr)
    return(as.numeric(scale(raw_score)))
  }
}

# --- 2.2 数据加载 (废弃 Sample_Sheet，直连 V2.0 临床表) ---
prepare_metadata <- function() {
  message(">>> [Init] Loading Clinical Metadata...")
  clinical <- read_csv(file.path(DIR_DATA, "Clinical_Master_Strict.csv"), show_col_types = F) %>% 
    dplyr::mutate(
      Role_Orig = dplyr::case_when(Membercode == 1 ~ "Daughter", Membercode == 2 ~ "Mother", Membercode == 3 ~ "Father", TRUE ~ "Unknown"),
      Subject_ID = tolower(paste0(FamilyID, "_", Role_Orig))
    )
  return(clinical)
}

# --- 2.3 分析流: Limma -> GSEA -> Dotplot ---
run_gsea_analysis <- function(filename, prefix, clin_master, target_vars) {
  full_path <- file.path(DIR_DATA, filename)
  if(!file.exists(full_path)) { warning("File not found: ", full_path); return(NULL) }
  
  message(paste0("\n>>> [Analysis] Processing: ", prefix))
  
  # 1. 读取 & 清洗
  data_raw <- read_csv(full_path, show_col_types = F)
  colnames(data_raw)[1] <- "FeatureID"
  data_raw <- data_raw %>% dplyr::filter(!is.na(FeatureID) & str_trim(FeatureID) != "") %>% dplyr::distinct(FeatureID, .keep_all = T)
  expr_mat <- as.matrix(data_raw[, -1]); rownames(expr_mat) <- data_raw$FeatureID
  
  # 2. 列名直读解析引擎 (完美对接 V2.0 矩阵)
  meta_sub <- data.frame(Matrix_Col = colnames(expr_mat), stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      Clean_Col = tolower(gsub("_twin\\.[0-9]+", "", Matrix_Col)),
      TimeType = dplyr::case_when(grepl("pre", Clean_Col) ~ "pre", grepl("fast", Clean_Col) ~ "fast", TRUE ~ "Other"),
      Subject_ID = gsub("_fast|_pre|_post1h|_post3h", "", Clean_Col)
    ) %>%
    dplyr::inner_join(clin_master, by = "Subject_ID")
  
  # 过滤时间点
  if (!is.null(TIME_FILTER_KEYWORD) && TIME_FILTER_KEYWORD != "") {
    keep_idx <- grepl(TIME_FILTER_KEYWORD, meta_sub$TimeType, ignore.case = TRUE)
    if (sum(keep_idx) < 5) { message("    [Skip] Too few samples after time filter."); return(NULL) }
    meta_sub <- meta_sub[keep_idx, ]
    # 防止同一个家系抽出两条 pre
    meta_sub <- meta_sub %>% dplyr::group_by(Subject_ID) %>% dplyr::slice(1) %>% dplyr::ungroup()
  }
  
  mat_sub <- expr_mat[, meta_sub$Matrix_Col]
  
  # 3. 预处理 & 血液校正
  mat_sub[is.na(mat_sub)] <- min(mat_sub, na.rm = T)/2
  if(max(mat_sub, na.rm = T) > 100) mat_sub <- log2(mat_sub + 1)
  blood_score_vec <- calculate_blood_score(mat_sub)
  
  # 4. 循环临床指标
  for (trait in target_vars) {
    if (!trait %in% colnames(meta_sub)) next
    valid_idx <- which(!is.na(meta_sub[[trait]]))
    if (length(valid_idx) < 5) next
    
    curr_mat <- mat_sub[, valid_idx]; curr_blood <- blood_score_vec[valid_idx]
    curr_trait <- as.numeric(meta_sub[[trait]][valid_idx])
    
    # Limma
    design <- model.matrix(~ curr_trait + curr_blood)
    fit <- eBayes(lmFit(curr_mat, design))
    res_limma <- topTable(fit, coef = 2, number = Inf) %>% tibble::rownames_to_column("Gene")
    
    # GSEA
    gene_map <- tryCatch(bitr(res_limma$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"), error=function(e) NULL)
    if(!is.null(gene_map)) {
      gsea_input <- res_limma %>% dplyr::inner_join(gene_map, by=c("Gene"="SYMBOL")) %>% 
        dplyr::arrange(desc(t)) %>% dplyr::distinct(ENTREZID, .keep_all=T)
      gene_list <- gsea_input$t; names(gene_list) <- gsea_input$ENTREZID
      
      gsea_res <- gseKEGG(geneList=gene_list, organism='hsa', pvalueCutoff=1, verbose=F, seed=123)
      
      if(!is.null(gsea_res) && nrow(gsea_res) > 0) {
        gsea_res <- setReadable(gsea_res, OrgDb=org.Hs.eg.db, keyType="ENTREZID")
        
        # 结果目录
        out_sub_dir <- file.path(DIR_OUT, trait)
        if(!dir.exists(out_sub_dir)) dir.create(out_sub_dir, recursive = T)
        
        # 保存 CSV
        res_df <- as.data.frame(gsea_res)
        res_df$Count <- sapply(strsplit(res_df$core_enrichment, "/"), length)
        write.csv(res_df, file.path(out_sub_dir, paste0("FigS2_", prefix, "_", trait, "_GSEA.csv")), row.names = F)
        
        # === 绘制经典气泡图 (Dotplot) ===
        try({
          p_dot <- dotplot(gsea_res, showCategory=15, split=".sign") + 
            facet_grid(.~.sign) +
            scale_color_gradient(low = "#E64B35", high = "#3C5488") +
            ggtitle(paste0(prefix, " - ", trait)) +
            theme_bw(base_size = 12) +
            theme(
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
              axis.ticks = element_line(color = "black"), axis.text = element_text(color = "black", face="bold"),
              strip.background = element_blank(), strip.text = element_text(color = "black", face="bold")
            )
          
          ggsave(file.path(out_sub_dir, paste0("FigS2_Dotplot_", prefix, "_", trait, ".pdf")), p_dot, width = 9, height = 7)
        })
      }
    }
  }
}

# --- 2.4 绘制 Crosstalk 图 (最终修正版) ---
plot_crosstalk_final <- function(prefix_x, prefix_y, trait, label_x, label_y) {
  
  file_x <- file.path(DIR_OUT, trait, paste0("FigS2_", prefix_x, "_", trait, "_GSEA.csv"))
  file_y <- file.path(DIR_OUT, trait, paste0("FigS2_", prefix_y, "_", trait, "_GSEA.csv"))
  
  if(!file.exists(file_x) | !file.exists(file_y)) return(NULL)
  message(paste0("    -> Crosstalk: ", label_y, " vs ", label_x))
  
  df_x <- read_csv(file_x, show_col_types=F) %>% dplyr::select(ID, Description, NES, p.adjust, Count) %>% dplyr::rename(NES_X=NES, FDR_X=p.adjust, Count_X=Count)
  df_y <- read_csv(file_y, show_col_types=F) %>% dplyr::select(ID, NES, p.adjust, Count) %>% dplyr::rename(NES_Y=NES, FDR_Y=p.adjust, Count_Y=Count)
  df_merge <- dplyr::inner_join(df_x, df_y, by = "ID")
  
  if(nrow(df_merge) == 0) return(NULL)
  
  keywords <- KEYWORD_DICT[[trait]]
  if(is.null(keywords)) keywords <- KEYWORD_DICT[["DEFAULT"]]
  pattern_red <- paste(keywords, collapse = "|")
  
  df_merge <- df_merge %>%
    dplyr::mutate(
      Is_Significant = (FDR_X < 0.05 | FDR_Y < 0.05),
      Avg_Count = (Count_X + Count_Y) / 2,
      Quadrant = dplyr::case_when(
        NES_X > 0 & NES_Y > 0 ~ "Q1: Co-Up",
        NES_X < 0 & NES_Y > 0 ~ "Q2: Y-Up/X-Down",
        NES_X < 0 & NES_Y < 0 ~ "Q3: Co-Down",
        NES_X > 0 & NES_Y < 0 ~ "Q4: X-Up/Y-Down",
        TRUE ~ "Center"
      ),
      Color_Group = ifelse(Is_Significant, Quadrant, "NS"),
      Is_Target = grepl(pattern_red, Description, ignore.case = TRUE),
      Rank_Score = (abs(NES_X) + abs(NES_Y)) * (-log10(pmin(FDR_X, FDR_Y) + 1e-10))
    )
  
  labels_red_df <- df_merge %>% dplyr::filter(Is_Significant & Is_Target) %>% dplyr::arrange(desc(Rank_Score)) %>% dplyr::group_by(Quadrant) %>% dplyr::slice_head(n = 3) %>% dplyr::ungroup()
  red_ids <- labels_red_df$ID
  
  labels_black_df <- df_merge %>% dplyr::filter(Is_Significant, !ID %in% red_ids) %>% dplyr::group_by(Quadrant) %>% dplyr::arrange(desc(Rank_Score)) %>% dplyr::slice_head(n = 3) %>% dplyr::ungroup()
  label_ids <- c(red_ids, labels_black_df$ID)
  
  df_merge <- df_merge %>%
    dplyr::mutate(
      Label_Text = ifelse(ID %in% label_ids, str_wrap(Description, 15), NA),
      Label_Color = dplyr::case_when(ID %in% red_ids ~ "Target", ID %in% labels_black_df$ID ~ "Top", TRUE ~ "None")
    )
  
  limit_val <- max(abs(c(df_merge$NES_X, df_merge$NES_Y)), na.rm = T) * 1.1
  if(limit_val < 1) limit_val <- 1
  
  p <- ggplot(df_merge, aes(x = NES_X, y = NES_Y)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
    geom_point(aes(size = Avg_Count, color = Color_Group, alpha = Is_Significant)) +
    geom_text_repel(aes(label = Label_Text, color = Label_Color), 
                    size = 3, max.overlaps = 50, seed = 42, 
                    bg.color = "white", bg.r = 0.1, box.padding = 0.3, min.segment.length = 0) +
    scale_color_manual(values = c(
      "Q1: Co-Up"="#d73027", "Q2: Y-Up/X-Down"="#fdae61", 
      "Q3: Co-Down"="#4575b4", "Q4: X-Up/Y-Down"="#abd9e9", 
      "NS" = "grey90", "Target" = "#D00000", "Top" = "black"    
    ), breaks = c("Q1: Co-Up", "Q2: Y-Up/X-Down", "Q3: Co-Down", "Q4: X-Up/Y-Down")) + 
    scale_alpha_manual(values = c(`FALSE`=0.4, `TRUE`=0.9), guide="none") +
    scale_size_continuous(range = c(2, 8), name = "Gene Count") +
    labs(title = paste0("Crosstalk: ", label_y, " vs ", label_x), 
         subtitle = paste0("Trait: ", trait, " (", TIME_FILTER_KEYWORD, ")"),
         x = paste0("NES (", label_x, ")"), y = paste0("NES (", label_y, ")")) +
    theme_bw(base_size = 12) + 
    theme(
      panel.grid = element_blank(), legend.position = "right", aspect.ratio = 1,
      panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
      axis.ticks = element_line(color = "black"), axis.text = element_text(color = "black", face="bold")
    ) +
    coord_fixed(xlim=c(-limit_val, limit_val), ylim=c(-limit_val, limit_val))
  
  fname <- paste0("FigS2_Crosstalk_", trait, "_", str_remove(label_y, " "), "_vs_", str_remove(label_x, " "), ".pdf")
  ggsave(file.path(DIR_OUT, trait, fname), p, width = 7.5, height = 7)
}

# ==============================================================================
# [3] 执行主程序
# ==============================================================================

if(!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = T)
meta_data <- prepare_metadata()

message("\n========== Step 1: Running Analysis & Dotplots ==========")
run_gsea_analysis("Cleaned_Adipose_Microarray.csv", "Adi_RNA", meta_data, TARGET_VARS)
run_gsea_analysis("Cleaned_Muscle_Microarray.csv",  "Mus_RNA", meta_data, TARGET_VARS)
run_gsea_analysis("Cleaned_Adipose_Proteomics.csv", "Adi_Prot", meta_data, TARGET_VARS)
run_gsea_analysis("Cleaned_Muscle_Proteomics.csv",  "Mus_Prot", meta_data, TARGET_VARS)

message("\n========== Step 2: Generating 4 Key Crosstalk Plots ==========")
for(trait in TARGET_VARS) {
  plot_crosstalk_final("Mus_RNA", "Adi_RNA", trait, "Muscle RNA", "Adipose RNA")
  plot_crosstalk_final("Mus_Prot", "Adi_Prot", trait, "Muscle Prot", "Adipose Prot")
  plot_crosstalk_final("Adi_Prot", "Adi_RNA", trait, "Adi Prot", "Adi RNA")
  plot_crosstalk_final("Mus_Prot", "Mus_RNA", trait, "Mus Prot", "Mus RNA")
}

message("\n>>> All Done! Check folder: ", DIR_OUT)