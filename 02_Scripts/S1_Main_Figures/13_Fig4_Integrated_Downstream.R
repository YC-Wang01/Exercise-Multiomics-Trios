# ==============================================================================
# Project: FinlandSports V2.0 - Figure 4 Integrated Downstream Pipeline
# Description: 
#   1. 3-Axis Core Screening (Two-way ANOVA + limma with Blood Correction)
#   2. 3-Pillars GO-BP ORA & Full-Spectrum Interaction GSEA
#   3. Rasterized Quadrant Plots (Tissue & Serum Separated)
#   4. Horizontal Paired Heatmap (Pathway Core)
#   5. HOMA-IR Centric Separate Networks (Obese vs Lean)
# Features: V2.0 In-Memory Delta Engine, 100% English, Strict Namespacing
# ==============================================================================

# --- 0. Environment Setup ---
setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readr, dplyr, tidyr, stringr, tibble, ggplot2, openxlsx,
  limma, car, clusterProfiler, org.Hs.eg.db, enrichplot,
  ComplexHeatmap, circlize, grid, cowplot, igraph, ggraph, ggforce,
  ggrastr, Hmisc, reshape2, ggrepel
)

# --- 1. Unified Path Configuration ---
DIR_DATA <- "01_Clean_Data"
DIR_OUT  <- "03_Results/Fig_4/Fig4_Integrated_Downstream"
if(!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)

message(">>> Directories Ready. Output set to: ", DIR_OUT)

# --- 2. V2.0 Core Data Loaders ---
message("\n=== STEP 1: Loading V2.0 Clinical Metadata & Engines ===")

clinical <- read_csv(file.path(DIR_DATA, "Clinical_Master_Strict.csv"), show_col_types = FALSE) %>% 
  dplyr::mutate(
    Role_Orig = dplyr::case_when(Membercode == 1 ~ "Daughter", Membercode == 2 ~ "Mother", Membercode == 3 ~ "Father", TRUE ~ "Unknown"),
    Subject_ID = tolower(paste0(FamilyID, "_", Role_Orig)),
    Fat_percent = as.numeric(`Fat(%)`),
    Fat_Group = ifelse(Fat_percent >= 30, "Obese", "Lean")
  )

# In-Memory Delta Matrix Generator
get_delta_matrix <- function(filename) {
  fpath <- file.path(DIR_DATA, filename)
  if(!file.exists(fpath)) return(NULL)
  
  data_raw <- read_csv(fpath, show_col_types = FALSE)
  expr_mat <- as.matrix(data_raw[, -1]); rownames(expr_mat) <- data_raw[[1]]
  
  clean_cols <- tolower(str_remove(colnames(expr_mat), "_twin\\.\\d+$")) %>% str_trim()
  col_meta <- data.frame(Matrix_Col = colnames(expr_mat), Clean_Col = clean_cols, stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      TimePoint = dplyr::case_when(grepl("pre", Clean_Col) ~ "pre", grepl("fast", Clean_Col) ~ "fast", grepl("post3h", Clean_Col) ~ "post3h", TRUE ~ "Other"),
      Subject_ID = gsub("_fast|_pre|_post1h|_post3h", "", Clean_Col)
    )
  
  subjects <- unique(clinical$Subject_ID)
  delta_list <- list()
  
  for(s in subjects) {
    pre_ids <- col_meta %>% dplyr::filter(Subject_ID == s, TimePoint %in% c("pre", "fast")) %>% dplyr::pull(Matrix_Col)
    post_ids <- col_meta %>% dplyr::filter(Subject_ID == s, TimePoint == "post3h") %>% dplyr::pull(Matrix_Col)
    
    if(length(pre_ids) > 0 && length(post_ids) > 0) {
      v_pre <- if(length(pre_ids) > 1) rowMeans(expr_mat[, pre_ids, drop=F], na.rm=T) else expr_mat[, pre_ids]
      v_post <- if(length(post_ids) > 1) rowMeans(expr_mat[, post_ids, drop=F], na.rm=T) else expr_mat[, post_ids]
      
      if(max(v_pre, na.rm=T) > 50) v_pre <- log2(v_pre + 1)
      if(max(v_post, na.rm=T) > 50) v_post <- log2(v_post + 1)
      
      delta_list[[s]] <- v_post - v_pre
    }
  }
  if(length(delta_list) < 3) return(NULL)
  return(do.call(cbind, delta_list))
}

# In-Memory Baseline Matrix Generator (For Microarrays)
get_baseline_matrix <- function(filename) {
  fpath <- file.path(DIR_DATA, filename)
  if(!file.exists(fpath)) return(NULL)
  data_raw <- read_csv(fpath, show_col_types = FALSE)
  expr_mat <- as.matrix(data_raw[, -1]); rownames(expr_mat) <- data_raw[[1]]
  
  clean_cols <- tolower(str_remove(colnames(expr_mat), "_twin\\.\\d+$")) %>% str_trim()
  col_meta <- data.frame(Matrix_Col = colnames(expr_mat), Clean_Col = clean_cols, stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      TimePoint = dplyr::case_when(grepl("pre", Clean_Col) ~ "pre", grepl("fast", Clean_Col) ~ "fast", TRUE ~ "Other"),
      Subject_ID = gsub("_fast|_pre|_post1h|_post3h", "", Clean_Col)
    ) %>% 
    dplyr::filter(TimePoint %in% c("pre", "fast")) %>%
    dplyr::group_by(Subject_ID) %>% dplyr::slice(1) %>% dplyr::ungroup()
  
  mat <- expr_mat[, col_meta$Matrix_Col, drop=FALSE]
  if(max(mat, na.rm=T) > 50) mat <- log2(mat + 1)
  colnames(mat) <- col_meta$Subject_ID
  return(mat)
}

# ==============================================================================
# MODULE 1: 3-Axis Core Screening Pipeline (Math & Stats Only)
# ==============================================================================
message("\n=== MODULE 1: 3-Axis Core Math Engine (ANOVA & Limma) ===")

datasets <- c(
  "Cleaned_Adipose_Proteomics", "Cleaned_Muscle_Proteomics", 
  "Cleaned_Serum_Proteomics", "Cleaned_Serum_Metabonomics",
  "Cleaned_Adipose_Methylation", "Cleaned_Muscle_Methylation",
  "Cleaned_Adipose_Microarray", "Cleaned_Muscle_Microarray"
)

# We focus ONLY on HOMA per user request, but pipeline supports expansibility
targets <- list(
  list(col="HOMA_IR", label="HOMA", out_file="Table_Fig4_HOMA_MathOnly.xlsx")
)

MATH_RESULTS <- list()

for (t_info in targets) {
  t_col <- t_info$col; t_lbl <- t_info$label; t_out <- t_info$out_file
  message(sprintf("  > Scanning Axis: %s", t_lbl))
  
  clin_df <- clinical[!is.na(clinical[[t_col]]), ]
  t_median <- median(clin_df[[t_col]], na.rm = TRUE)
  clin_df$Target_Group <- ifelse(clin_df[[t_col]] >= t_median, paste0("High", t_lbl), paste0("Low", t_lbl))
  
  wb <- createWorkbook()
  
  # --- Clinical Baseline (ANOVA) ---
  clin_results <- data.frame()
  num_cols <- names(clin_df)[sapply(clin_df, is.numeric)]
  exclude_cols <- c("FamilyID", "Membercode", "Unique_Key", "Fat_percent", "HOMA_IR", "Age", "Gender")
  num_cols <- setdiff(num_cols, exclude_cols)
  
  for (tgt in num_cols) {
    plot_data <- clin_df[!is.na(clin_df[[tgt]]), ]
    if (nrow(plot_data) < 20) next
    formula_str <- paste0("`", tgt, "` ~ Fat_Group * Target_Group")
    aov_res <- tryCatch(car::Anova(aov(as.formula(formula_str), data = plot_data), type = 2), error = function(e) NULL)
    
    if (!is.null(aov_res)) {
      if (any(c(aov_res["Fat_Group", "Pr(>F)"], aov_res["Target_Group", "Pr(>F)"], aov_res["Fat_Group:Target_Group", "Pr(>F)"]) < 0.05, na.rm = TRUE)) {
        clin_results <- rbind(clin_results, data.frame(Feature = tgt, P_Fat = aov_res["Fat_Group", "Pr(>F)"], P_Target = aov_res["Target_Group", "Pr(>F)"], P_Interaction = aov_res["Fat_Group:Target_Group", "Pr(>F)"]))
      }
    }
  }
  if(nrow(clin_results) > 0) {
    clin_results <- clin_results[order(clin_results$P_Interaction), ]
    addWorksheet(wb, "Clinical_Baseline")
    writeData(wb, "Clinical_Baseline", clin_results)
  }
  
  # --- Omics Data (Limma) ---
  for (ds in datasets) {
    is_microarray <- grepl("Microarray", ds)
    is_serum <- grepl("Serum", ds)
    
    mat <- if(is_microarray) get_baseline_matrix(paste0(ds, ".csv")) else get_delta_matrix(paste0(ds, ".csv"))
    if(is.null(mat)) next
    
    common_ids <- intersect(colnames(mat), clin_df$Subject_ID)
    if(length(common_ids) < 5) next
    
    mat_sub <- mat[, common_ids, drop=FALSE]
    pheno <- clin_df %>% dplyr::filter(Subject_ID %in% common_ids) %>% dplyr::arrange(match(Subject_ID, common_ids))
    
    pheno$Fat_Group <- factor(pheno$Fat_Group, levels = c("Lean", "Obese"))
    pheno$Target_Group <- factor(pheno$Target_Group, levels = c(paste0("Low", t_lbl), paste0("High", t_lbl)))
    
    has_HBB <- "HBB" %in% rownames(mat_sub)
    if(has_HBB && !is_serum) {
      pheno$Blood_Cov <- mat_sub["HBB", ]
      formula_limma <- "~ Fat_Group * Target_Group + Blood_Cov"
    } else {
      formula_limma <- "~ Fat_Group * Target_Group"
    }
    
    design <- model.matrix(as.formula(formula_limma), data = pheno)
    fit <- eBayes(lmFit(mat_sub, design))
    
    p_vals <- fit$p.value
    coef_fat <- "Fat_GroupObese"; coef_tar <- paste0("Target_GroupHigh", t_lbl); coef_int <- paste0("Fat_GroupObese:Target_GroupHigh", t_lbl)
    
    res <- data.frame(
      Feature = rownames(p_vals), P_Fat = p_vals[, coef_fat], P_Target = p_vals[, coef_tar], P_Interaction = p_vals[, coef_int],
      check.names = FALSE, stringsAsFactors = FALSE
    )
    
    sig_res <- res %>% dplyr::filter(P_Fat < 0.05 | P_Target < 0.05 | P_Interaction < 0.05) %>% dplyr::arrange(P_Interaction)
    
    sheet_name <- substr(gsub("Cleaned_", "", ds), 1, 31)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, sig_res)
    
    MATH_RESULTS[[paste0(t_lbl, "_", sheet_name)]] <- res # Save full result for GSEA
  }
  saveWorkbook(wb, file.path(DIR_OUT, t_out), overwrite = TRUE)
}

# ==============================================================================
# MODULE 2: GO ORA & Interaction GSEA
# ==============================================================================
message("\n=== MODULE 2: 3-Pillars GO-BP ORA & Full-Spectrum GSEA ===")

clean_genes <- function(g_vec) unique(na.omit(sapply(strsplit(as.character(g_vec), ";"), `[`, 1)))
TOP_N_ORA <- 7

for (t_info in targets) {
  t_lbl <- t_info$label
  
  # --- GO ORA ---
  sheet_key <- paste0(t_lbl, "_Muscle_Proteomics")
  if(sheet_key %in% names(MATH_RESULTS)) {
    df <- MATH_RESULTS[[sheet_key]]
    df$Feature <- clean_genes(df$Feature)
    
    gene_clusters <- list()
    g_fat <- df$Feature[df$P_Fat < 0.05]; g_tar <- df$Feature[df$P_Target < 0.05]; g_int <- df$Feature[df$P_Interaction < 0.05]
    
    if(length(g_fat)>0) gene_clusters[["Obesity"]] <- bitr(g_fat, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)$ENTREZID
    if(length(g_tar)>0) gene_clusters[[paste0(t_lbl)]] <- bitr(g_tar, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)$ENTREZID
    if(length(g_int)>0) gene_clusters[["Interaction"]] <- bitr(g_int, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)$ENTREZID
    
    if(length(gene_clusters) > 0) {
      ck_go <- suppressMessages(compareCluster(geneCluster = gene_clusters, fun = "enrichGO", OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05))
      if(!is.null(ck_go) && nrow(as.data.frame(ck_go)) > 0) {
        ck_go <- setReadable(ck_go, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
        df_go <- as.data.frame(ck_go) %>% dplyr::group_by(Cluster) %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = TOP_N_ORA) %>% dplyr::ungroup()
        df_go$GeneRatio_Num <- sapply(df_go$GeneRatio, function(x) { p <- as.numeric(strsplit(x, "/")[[1]]); p[1]/p[2] })
        df_go$Description <- paste0(toupper(substr(df_go$Description, 1, 1)), substring(df_go$Description, 2))
        
        pillar_levels <- c("Obesity", t_lbl, "Interaction")
        df_go$Cluster <- factor(df_go$Cluster, levels = pillar_levels)
        df_go$Description <- str_wrap(df_go$Description, width = 45)
        df_go$Description <- factor(df_go$Description, levels = unique(df_go$Description[order(df_go$p.adjust, decreasing = TRUE)]))
        
        p_go <- ggplot(df_go, aes(x = Cluster, y = Description)) +
          geom_point(aes(size = GeneRatio_Num, fill = p.adjust), shape = 21, color = "black", stroke = 0.7) +
          scale_fill_gradientn(colors = c("#E74C3C", "#F39C12", "#F1C40F", "#3498DB"), name = "p.adjust") +
          scale_size_continuous(range = c(4, 10)) + 
          scale_x_discrete(drop = FALSE) + scale_y_discrete(position = "right") +
          theme_classic(base_size = 14) + 
          labs(title = paste("GO-BP 3-Pillars:", t_lbl, "Axis"), x=NULL, y=NULL) +
          theme(panel.grid.major.x = element_line(color="gray80", linetype="dotted"), 
                panel.grid.major.y = element_line(color="gray90", linetype="dashed"),
                panel.border = element_rect(color="black", fill=NA, linewidth=1.5), axis.line=element_blank(),
                plot.title=element_text(face="bold", hjust=0.5, size=16), 
                axis.text.x=element_text(face="bold", color="black"),
                axis.text.y.right=element_text(face="bold", color="black", size=10, hjust=0),
                aspect.ratio = max(1.0, nrow(df_go) * 0.15)) 
        ggsave(file.path(DIR_OUT, paste0("Fig4_3Pillars_GO_", t_lbl, ".pdf")), p_go, width=8.5, height=max(3.5, nrow(df_go)*0.2 + 1.5), limitsize=FALSE)
      }
    }
  }
  
  # --- GSEA Interaction ---
  # Reconstruct t-values using Limma
  mat <- get_delta_matrix("Cleaned_Muscle_Proteomics.csv")
  clin_df <- clinical[!is.na(clinical[[t_info$col]]), ]
  clin_df$Target_Group <- ifelse(clin_df[[t_info$col]] >= median(clin_df[[t_info$col]], na.rm=T), "High", "Low")
  common_ids <- intersect(colnames(mat), clin_df$Subject_ID)
  
  mat_sub <- mat[, common_ids, drop=FALSE]
  pheno <- clin_df %>% dplyr::filter(Subject_ID %in% common_ids) %>% dplyr::arrange(match(Subject_ID, common_ids))
  pheno$Fat_Group <- factor(pheno$Fat_Group, levels=c("Lean", "Obese"))
  pheno$Target_Group <- factor(pheno$Target_Group, levels=c("Low", "High"))
  
  if("HBB" %in% rownames(mat_sub)) {
    pheno$Blood_Cov <- mat_sub["HBB", ]
    design <- model.matrix(~ Fat_Group * Target_Group + Blood_Cov, data=pheno)
  } else { design <- model.matrix(~ Fat_Group * Target_Group, data=pheno) }
  
  fit <- eBayes(lmFit(mat_sub, design))
  t_scores <- fit$t[, "Fat_GroupObese:Target_GroupHigh"]
  
  gene_list <- suppressMessages(bitr(names(t_scores), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=FALSE)) %>% na.omit()
  gene_list$t_value <- t_scores[gene_list$SYMBOL]
  ranked_list <- sort(setNames(gene_list$t_value, gene_list$ENTREZID), decreasing = TRUE)
  
  set.seed(123)
  gsea_go <- suppressMessages(gseGO(geneList = ranked_list, OrgDb = org.Hs.eg.db, ont = "BP", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE))
  
  if(!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
    gsea_go <- setReadable(gsea_go, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    write.xlsx(as.data.frame(gsea_go), file.path(DIR_OUT, paste0("Table_GSEA_Interaction_", t_lbl, ".xlsx")))
    
    top_n <- min(4, nrow(as.data.frame(gsea_go)))
    for(i in 1:top_n) {
      top_pathway_name <- as.data.frame(gsea_go)$Description[i]
      p_gsea <- gseaplot2(gsea_go, geneSetID = as.data.frame(gsea_go)$ID[i], title = paste(t_lbl, "Interaction:", top_pathway_name), pvalue_table = TRUE)
      safe_name <- gsub("[^A-Za-z0-9_.-]", "_", top_pathway_name)
      ggsave(file.path(DIR_OUT, paste0("Fig4_GSEA_Plot_", t_lbl, "_Rank", i, "_", substr(safe_name, 1, 30), ".pdf")), p_gsea, width=8, height=6.5)
    }
  }
}

# ==============================================================================
# MODULE 3: Rasterized Quadrant Plots
# ==============================================================================
message("\n=== MODULE 3: Generating Rasterized Quadrant Plots ===")

layout_config <- list(
  "Tissue" = list(datasets = c("Adipose" = "Cleaned_Adipose_Proteomics.csv", "Muscle" = "Cleaned_Muscle_Proteomics.csv"), colors = c("NS" = "gray85", "Adipose_Sig" = "#F79647", "Muscle_Sig" = "#3DA6AE"), borders = c("NS" = "gray70", "Adipose_Sig" = "black", "Muscle_Sig" = "black"), sizes = c("NS" = 1.5, "Adipose_Sig" = 3.5, "Muscle_Sig" = 3.5), subtitle = "Adipose Prot (Orange) vs. Muscle Prot (Cyan)"),
  "Serum" = list(datasets = c("SerumProt" = "Cleaned_Serum_Proteomics.csv", "SerumMetab" = "Cleaned_Serum_Metabonomics.csv"), colors = c("NS" = "gray85", "SerumProt_Sig" = "#9B59B6", "SerumMetab_Sig" = "#E74C3C"), borders = c("NS" = "gray70", "SerumProt_Sig" = "black", "SerumMetab_Sig" = "black"), sizes = c("NS" = 1.5, "SerumProt_Sig" = 3.5, "SerumMetab_Sig" = 3.5), subtitle = "Serum Proteomics (Purple) & Metabonomics (Red)")
)

for (t_info in targets) {
  t_col <- t_info$col; t_lbl <- t_info$label
  clin_df <- clinical[!is.na(clinical[[t_col]]), ]
  clin_df$Target_Group <- ifelse(clin_df[[t_col]] >= median(clin_df[[t_col]], na.rm=T), "High", "Low")
  
  for (plot_type in names(layout_config)) {
    conf <- layout_config[[plot_type]]
    plot_data <- data.frame()
    
    for (ds_name in names(conf$datasets)) {
      mat <- get_delta_matrix(conf$datasets[[ds_name]])
      if(is.null(mat)) next
      
      id_lean_low <- clin_df$Subject_ID[clin_df$Fat_Group == "Lean" & clin_df$Target_Group == "Low"]
      id_obese_high <- clin_df$Subject_ID[clin_df$Fat_Group == "Obese" & clin_df$Target_Group == "High"]
      
      common_lean <- intersect(id_lean_low, colnames(mat))
      common_obese <- intersect(id_obese_high, colnames(mat))
      if(length(common_lean)<3 || length(common_obese)<3) next
      
      mean_lean_low <- rowMeans(mat[, common_lean, drop=FALSE], na.rm=TRUE)
      mean_obese_high <- rowMeans(mat[, common_obese, drop=FALSE], na.rm=TRUE)
      
      # Use pre-calculated P-values from Math Engine
      p_vals <- rep(1, nrow(mat)); names(p_vals) <- rownames(mat)
      math_key <- paste0(t_lbl, "_", substr(gsub("Cleaned_", "", conf$datasets[[ds_name]]), 1, 31))
      if(math_key %in% names(MATH_RESULTS)) {
        res_df <- MATH_RESULTS[[math_key]]
        match_idx <- match(rownames(mat), res_df$Feature)
        p_vals[!is.na(match_idx)] <- res_df$P_Interaction[na.omit(match_idx)]
      }
      
      plot_data <- rbind(plot_data, data.frame(Feature = rownames(mat), DatasetType = ds_name, Mean_Lean_Low = mean_lean_low, Mean_Obese_High = mean_obese_high, P_Interaction = p_vals, stringsAsFactors = FALSE))
    }
    
    if(nrow(plot_data) == 0) next
    
    plot_data <- plot_data[!is.na(plot_data$P_Interaction), ]
    plot_data$Category <- "NS"
    sig_classes <- names(conf$colors)[names(conf$colors) != "NS"]
    for (cls in sig_classes) plot_data$Category[plot_data$P_Interaction < 0.05 & plot_data$DatasetType == gsub("_Sig", "", cls)] <- cls
    plot_data$Category <- factor(plot_data$Category, levels = c("NS", sig_classes))
    plot_data <- plot_data[order(plot_data$Category), ]
    
    top_labels <- plot_data %>% dplyr::filter(Category != "NS") %>% dplyr::group_by(DatasetType) %>% dplyr::slice_min(order_by = P_Interaction, n = 8)
    axis_limit <- max(abs(c(plot_data$Mean_Lean_Low, plot_data$Mean_Obese_High)), na.rm = TRUE) * 1.1
    
    p_quad <- ggplot(plot_data, aes(x = Mean_Lean_Low, y = Mean_Obese_High)) +
      geom_hline(yintercept = 0, color = "gray40", linetype = "dotted", size = 0.8) +
      geom_vline(xintercept = 0, color = "gray40", linetype = "dotted", size = 0.8) +
      geom_abline(intercept = 0, slope = 1, color = "gray40", linetype = "dashed", size = 0.8) +
      rasterise(geom_point(aes(fill = Category, color = Category, size = Category), shape = 21, stroke = 0.4, alpha = 0.85), dpi = 300) +
      scale_fill_manual(values = conf$colors) + scale_color_manual(values = conf$borders) + scale_size_manual(values = conf$sizes) +
      geom_text_repel(data = top_labels, aes(label = Feature), size = 4.5, fontface = "bold", color = "black", box.padding = 0.5, point.padding = 0.5, max.overlaps = 30) +
      scale_x_continuous(limits = c(-axis_limit, axis_limit)) + scale_y_continuous(limits = c(-axis_limit, axis_limit)) +
      theme_classic(base_size = 15) +
      labs(title = paste("Systemic Crosstalk (Delta):", t_lbl, "Axis"), subtitle = conf$subtitle, x = paste("Change in Lean & Low", t_lbl), y = paste("Change in Obese & High", t_lbl)) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16), plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic", color = "darkred"), legend.position = "none", aspect.ratio = 1, panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2), axis.line = element_blank())
    
    ggsave(file.path(DIR_OUT, paste0("Fig4_QuadrantPlot_", plot_type, "_", t_lbl, ".pdf")), p_quad, width = 7.5, height = 7.5)
  }
}

# ==============================================================================
# MODULE 4: Horizontal Paired Heatmap (Pathway Core)
# ==============================================================================
message("\n=== MODULE 4: Generating Horizontal Paired Heatmap ===")

mat_adi <- get_delta_matrix("Cleaned_Adipose_Proteomics.csv")
mat_mus <- get_delta_matrix("Cleaned_Muscle_Proteomics.csv")
common_subjects <- intersect(clinical$Subject_ID[!is.na(clinical$HOMA_IR)], intersect(colnames(mat_adi), colnames(mat_mus)))

anno_df <- clinical %>% dplyr::filter(Subject_ID %in% common_subjects) %>% dplyr::arrange(factor(Fat_Group, levels = c("Lean", "Obese")), HOMA_IR)
mat_adi <- mat_adi[, anno_df$Subject_ID, drop=FALSE]; mat_mus <- mat_mus[, anno_df$Subject_ID, drop=FALSE]

pathway_genes <- list(
  "Mitochondrial Translation" = c("MRPL40", "MRPL12", "MRPL55", "MRPL38", "TFAM", "SHMT2"),
  "Mitophagy & Proteolysis"   = c("MUL1", "RPS27L", "SMAD3", "PYCARD", "VCP", "DDX3X", "IDH3G"),
  "NAD & Nucleotide Synth"    = c("NMNAT3", "CTPS1", "PRPS1", "COX11", "PDK4", "NDUFS5")
)

col_list <- list(); anno_list <- list()
for(p_name in names(pathway_genes)) {
  for(g in pathway_genes[[p_name]]) {
    if(g %in% rownames(mat_adi) || g %in% rownames(mat_mus)) {
      # Adipose Side
      if(g %in% rownames(mat_adi)) {
        v_adi <- mat_adi[g, anno_df$Subject_ID]
        z_adi <- if(sd(v_adi, na.rm=T) > 0) (v_adi - mean(v_adi, na.rm=T))/sd(v_adi, na.rm=T) else rep(0, length(v_adi))
        ct_adi <- suppressWarnings(cor.test(z_adi, anno_df$HOMA_IR, method="spearman", exact=FALSE))
        r_adi <- unname(ct_adi$estimate); p_adi <- unname(ct_adi$p.value)
      } else { z_adi <- rep(NA, nrow(anno_df)); r_adi <- NA; p_adi <- NA }
      
      # Muscle Side
      if(g %in% rownames(mat_mus)) {
        v_mus <- mat_mus[g, anno_df$Subject_ID]
        z_mus <- if(sd(v_mus, na.rm=T) > 0) (v_mus - mean(v_mus, na.rm=T))/sd(v_mus, na.rm=T) else rep(0, length(v_mus))
        ct_mus <- suppressWarnings(cor.test(z_mus, anno_df$HOMA_IR, method="spearman", exact=FALSE))
        r_mus <- unname(ct_mus$estimate); p_mus <- unname(ct_mus$p.value)
      } else { z_mus <- rep(NA, nrow(anno_df)); r_mus <- NA; p_mus <- NA }
      
      col_list[[paste0(g, "_Adipose")]] <- z_adi; col_list[[paste0(g, "_Muscle")]] <- z_mus
      anno_list[[length(anno_list) + 1]] <- data.frame(Feature = paste0(g, "_Adipose"), Gene = g, Tissue = "Adipose", Pathway = p_name, R_value = r_adi, P_value = p_adi, stringsAsFactors=FALSE)
      anno_list[[length(anno_list) + 1]] <- data.frame(Feature = paste0(g, "_Muscle"), Gene = g, Tissue = "Muscle", Pathway = p_name, R_value = r_mus, P_value = p_mus, stringsAsFactors=FALSE)
    }
  }
}

z_mat_all <- do.call(rbind, col_list)
colnames(z_mat_all) <- anno_df$Subject_ID
row_anno_df <- do.call(rbind, anno_list)

ht_mat <- t(z_mat_all)
gene_labels <- ifelse(row_anno_df$Tissue == "Adipose", row_anno_df$Gene, "")
col_fun_heatmap <- colorRamp2(c(-1.2, 0, 1.2), c("#00468B", "white", "#ED0000"))
col_homa <- colorRamp2(c(min(anno_df$HOMA_IR), median(anno_df$HOMA_IR), max(anno_df$HOMA_IR)), c("white", "#FFC107", "#DC0000"))

left_anno <- rowAnnotation(
  `Fat(%)` = anno_df$Fat_Group, `HOMA-IR` = anno_df$HOMA_IR,
  col = list(`Fat(%)` = c("Obese" = "#F79647", "Lean" = "#3DA6AE"), `HOMA-IR` = col_homa),
  show_annotation_name = TRUE, annotation_name_side = "bottom", annotation_name_rot = 45, annotation_name_gp = gpar(fontsize = 12, fontface = "bold", col = "black")
)

safe_r_value <- ifelse(is.na(row_anno_df$R_value), 0, row_anno_df$R_value)
top_anno <- HeatmapAnnotation(
  Tissue = row_anno_df$Tissue, Cor_R = anno_barplot(safe_r_value, bar_width = 1, gp = gpar(fill = ifelse(safe_r_value > 0, "#ED0000", "#00468B"), col=NA)),
  col = list(Tissue = c("Muscle" = "#3C5488", "Adipose" = "#E64B35")), show_annotation_name = TRUE, annotation_name_side = "left", annotation_name_rot = c(Tissue = 90, Cor_R = 0), simple_anno_size = unit(4, "mm")
)

ht <- Heatmap(
  ht_mat, name = "Z-Score\n(Delta)", col = col_fun_heatmap, na_col = "#E0E0E0", 
  width = ncol(ht_mat) * unit(6, "mm"), height = nrow(ht_mat) * unit(6, "mm"),
  left_annotation = left_anno, top_annotation = top_anno,
  row_split = factor(anno_df$Fat_Group, levels = c("Lean", "Obese")), column_split = factor(row_anno_df$Pathway, levels = unique(row_anno_df$Pathway)),
  row_gap = unit(3, "mm"), column_gap = unit(4, "mm"), cluster_rows = FALSE, cluster_columns = FALSE, 
  show_row_names = FALSE, show_column_names = TRUE, column_labels = gene_labels, column_names_side = "bottom", column_names_rot = 45, column_names_gp = gpar(fontsize = 12, fontface = "bold.italic", col = "black"),
  row_title_gp = gpar(fontsize = 12, fontface = "bold", rot = 90), column_title_gp = gpar(fontsize = 12, fontface = "bold"), rect_gp = gpar(col = "white", lwd = 1.5), border = TRUE
)

pdf(file.path(DIR_OUT, "Fig4_Waterfall_Horizontal_Paired.pdf"), width = max(12, ncol(ht_mat) * 0.35 + 5), height = max(7, nrow(ht_mat) * 0.35 + 4))
draw(ht, merge_legend = TRUE, column_title = "Cross-Tissue Core Pathway Responses vs HOMA-IR Gradient", column_title_gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()

# ==============================================================================
# MODULE 5: HOMA-IR Centric Separate Networks (Fixed Duplicate & Variables)
# ==============================================================================
message("\n=== MODULE 5: Generating HOMA-IR Centric Networks ===")

# 明确提取四个组学的 Delta 矩阵，统一命名规则，彻底杜绝找不到对象的问题
mat_mus  <- get_delta_matrix("Cleaned_Muscle_Proteomics.csv")
mat_adi  <- get_delta_matrix("Cleaned_Adipose_Proteomics.csv")
mat_met  <- get_delta_matrix("Cleaned_Serum_Metabonomics.csv")
mat_prot <- get_delta_matrix("Cleaned_Serum_Proteomics.csv")                                   

run_homa_network <- function(group_name) {
  # 统一使用上面刚声明的干净变量名
  mats_list <- list(mat_mus, mat_adi, mat_met, mat_prot)
  valid_mats <- mats_list[!sapply(mats_list, is.null)]
  
  all_subs <- Reduce(union, lapply(valid_mats, colnames))
  target_subs <- clinical %>% dplyr::filter(Subject_ID %in% all_subs, Fat_Group == group_name) %>% dplyr::pull(Subject_ID) %>% unique()
  
  if(length(target_subs) < 4) return(NULL)
  
  safe_subset <- function(mat, subs) {
    if (is.null(mat)) return(NULL)
    res <- matrix(NA, nrow = nrow(mat), ncol = length(subs)); rownames(res) <- rownames(mat); colnames(res) <- subs
    avail <- intersect(colnames(mat), subs); res[, avail] <- mat[, avail]; return(res)
  }
  
  mat_mus_sub  <- safe_subset(mat_mus, target_subs)
  mat_adi_sub  <- safe_subset(mat_adi, target_subs)
  mat_met_sub  <- safe_subset(mat_met, target_subs)
  mat_prot_sub <- safe_subset(mat_prot, target_subs)
  
  clin_mat <- matrix(NA, nrow=1, ncol=length(target_subs), dimnames=list("HOMA_IR", target_subs))
  clin_homa_vec <- clinical$HOMA_IR[match(target_subs, clinical$Subject_ID)]
  clin_mat[1, ] <- clin_homa_vec
  
  select_top_by_cor <- function(mat_sub, n) {
    if(is.null(mat_sub)) return(NULL)
    cor_vals <- apply(mat_sub, 1, function(x) {
      if(sum(complete.cases(x, clin_homa_vec)) < 4 || sd(x, na.rm=TRUE) == 0) return(0)
      abs(cor(x, clin_homa_vec, use="pairwise.complete.obs", method="spearman"))
    })
    names(sort(cor_vals, decreasing = TRUE, na.last = TRUE)[1:min(n, length(cor_vals))])
  }
  
  t_mus <- select_top_by_cor(mat_mus_sub, 15); t_adi <- select_top_by_cor(mat_adi_sub, 15)
  t_met <- select_top_by_cor(mat_met_sub, 15); t_prot <- select_top_by_cor(mat_prot_sub, 15)
  
  comb_list <- list()
  if(!is.null(t_mus)) comb_list[["Muscle"]] <- mat_mus_sub[t_mus, , drop=FALSE]
  if(!is.null(t_adi)) comb_list[["Adipose"]] <- mat_adi_sub[t_adi, , drop=FALSE]
  if(!is.null(t_met)) comb_list[["Metab"]] <- mat_met_sub[t_met, , drop=FALSE]
  if(!is.null(t_prot)) comb_list[["Prot"]] <- mat_prot_sub[t_prot, , drop=FALSE]
  comb_list[["Clinical"]] <- clin_mat
  
  comb_mat <- do.call(rbind, comb_list)
  
  edges <- list()
  for(i in 1:(nrow(comb_mat)-1)) {
    for(j in (i+1):nrow(comb_mat)) {
      x <- comb_mat[i, ]; y <- comb_mat[j, ]; valid <- complete.cases(x, y)
      if(sum(valid) >= 4) {
        ct <- suppressWarnings(cor.test(x[valid], y[valid], method="spearman", exact=FALSE))
        if(ct$p.value < 0.01 && abs(ct$estimate) >= 0.65) {
          edges[[length(edges)+1]] <- data.frame(Node1=rownames(comb_mat)[i], Node2=rownames(comb_mat)[j], R=ct$estimate, P=ct$p.value)
        }
      }
    }
  }
  if(length(edges) == 0) return(NULL)
  edges_df <- do.call(rbind, edges)
  
  # 去重修复逻辑保持不动
  node_info <- data.frame(name = rownames(comb_mat)) %>%
    dplyr::mutate(
      Category = dplyr::case_when(name %in% t_adi ~ "Adipose", name %in% t_mus ~ "Muscle", name %in% t_met ~ "Metabolomics", name %in% t_prot ~ "Proteomics", TRUE ~ "Clinical"),
      Shape = ifelse(Category == "Clinical", "Clinical", "Molecule")
    ) %>% 
    dplyr::filter(name %in% c(edges_df$Node1, edges_df$Node2)) %>%
    dplyr::distinct(name, .keep_all = TRUE) 
  
  g_full <- graph_from_data_frame(edges_df, directed = F, vertices = node_info)
  if (!"HOMA_IR" %in% V(g_full)$name) return(NULL)
  
  g <- make_ego_graph(g_full, order = 2, nodes = "HOMA_IR")[[1]]
  g <- delete_vertices(g, V(g)[degree(g) < 1]) 
  if (vcount(g) == 0) return(NULL)
  
  V(g)$Degree <- degree(g) 
  layer_colors <- c("Adipose" = "#F39B7F", "Muscle" = "#4DBBD5", "Metabolomics" = "#8491B4", "Proteomics" = "#91D1C2", "Clinical" = "#E64B35")
  
  p_base <- ggraph(g, layout = if(vcount(g) <= 6) "circle" else "kk") +  
    ggforce::geom_mark_hull(data = ~ subset(., Category != "Clinical"), aes(x, y, group = Category, fill = Category, color = Category, label = Category), concavity = 2, expand = unit(4, "mm"), radius = unit(3, "mm"), alpha = 0.08, linetype = "dashed", linewidth = 0.6, label.fontsize = 10, label.fontface = "bold", label.fill = "white", label.margin = margin(5, 5, 5, 5), label.buffer = unit(25, "mm"), con.colour = "grey50", con.cap = unit(1, "mm")) +
    geom_edge_link(aes(edge_width = abs(R), color = R > 0, linetype = R > 0), alpha = 0.6) +
    scale_edge_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#3C5488"), labels = c("TRUE" = "Positive (+)", "FALSE" = "Negative (-)"), name = "Interaction Trend") +
    scale_edge_width(range = c(0.4, 1.8), name = "Correlation |R|") + scale_edge_linetype_manual(values = c("TRUE"="solid", "FALSE"="dashed"), guide = "none") +
    geom_node_point(aes(fill = Category, shape = Shape, size = Degree), color = "white", stroke = 1.2) +
    scale_fill_manual(values = layer_colors, name = "Biological Layer") + scale_color_manual(values = layer_colors, guide = "none") +
    scale_shape_manual(values = c("Molecule" = 21, "Clinical" = 23), guide = "none") + scale_size_continuous(range = c(5, 14), name = "Node Degree") +
    geom_node_text(aes(label = name), repel = TRUE, size = 3.5, family = "sans", fontface = "bold", bg.color = "white", bg.r = 0.15, segment.color = "grey60", segment.size = 0.5) +
    coord_cartesian(clip = "off") + theme_void() + 
    theme(text = element_text(family = "sans"), plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b=10)), plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40", margin = margin(b=20)), plot.margin = margin(30, 30, 30, 30), panel.background = element_rect(fill = "white", color = NA), legend.position = "right", legend.box = "vertical", legend.title = element_text(face = "bold", size = 11), legend.text = element_text(size = 10)) +
    labs(title = paste0("HOMA-IR Centric Core Crosstalk: ", group_name), subtitle = "Exercise-Induced Tissue-Serum Cascades Linked to Insulin Resistance")
  
  ggsave(file.path(DIR_OUT, paste0("Plot_Network_HOMAcentric_", group_name, "_Main.pdf")), p_base + theme(legend.position = "none"), width = 12, height = 10, dpi = 300)
  
  legend_obj <- cowplot::get_legend(p_base)
  if(!is.null(legend_obj)) ggsave(file.path(DIR_OUT, paste0("Plot_Network_HOMAcentric_", group_name, "_LegendOnly.pdf")), cowplot::ggdraw(legend_obj), width = 4, height = 7, dpi = 300)
  
  out_edges <- igraph::as_data_frame(g, what = "edges") %>%
    dplyr::mutate(
      Group = group_name,
      R = round(R, 4),       
      P = signif(P, 4),
      FDR = signif(p.adjust(P, method = "BH"), 4),
      Trend = ifelse(R > 0, "Positive", "Negative") 
    ) %>%
    dplyr::select(Group, Node1 = from, Node2 = to, Correlation_R = R, P_value = P, FDR, Trend) %>%
    dplyr::arrange(P_value)
  
  write_csv(out_edges, file.path(DIR_OUT, paste0("Table_S_HOMAcentric_Edges_", group_name, ".csv")))
  
  out_nodes <- data.frame(
    Molecule_Name = V(g)$name,
    Biological_Layer = V(g)$Category,
    Degree = V(g)$Degree
  ) %>%
    dplyr::arrange(Biological_Layer, desc(Degree)) %>% 
    dplyr::mutate(Group = group_name) %>%
    dplyr::select(Group, Molecule_Name, Biological_Layer, Degree)
  
  write_csv(out_nodes, file.path(DIR_OUT, paste0("Table_S_HOMAcentric_Nodes_", group_name, ".csv")))
  message("  > Generated network and saved tables for: ", group_name)
}

run_homa_network("Obese")
run_homa_network("Lean")

message("\n=======================================================")
message(">>> MODULE 5 COMPLETED! Check output directory.")
message("=======================================================")