# ==============================================================================
# Project: FinlandSports V2.0 - Fig 5 D&E: CRP Systemic Integration
# Description:
#   1. Strict Intersection: Only subjects with ALL THREE tissues are plotted.
#   2. Row Clustering: Enabled within-tissue clustering to show red/blue blocks.
#   3. Bubble Plot: Reduced max size, moved legend to right, loosened Y-spacing.
#   4. Gene Names: Cleaned up combined names (e.g., "ZNF76;ZNF143" -> "ZNF76").
#   5. Features: V2.0 Path Standard, 100% English, Strict Namespacing.
# ==============================================================================

# --- 0. Environment Setup ---
# setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readr, dplyr, stringr, limma, tibble, tidyr,
  ggplot2, ggrepel, ggsci, ggpubr, grid, 
  ComplexHeatmap, circlize, clusterProfiler, 
  org.Hs.eg.db, enrichplot, openxlsx, scales
)

# --- 1. PATH CONFIGURATION (V2.0) ---
DIR_IN     <- "01_Clean_Data"
DIR_OUT    <- "03_Results/Fig_5/Fig5DE_CRP_Downstream"
TABLE_DIR  <- file.path(DIR_OUT, "Results_Tables")

dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

TARGETS <- list(
  Adipose = "Cleaned_Adipose_Proteomics.csv", 
  Muscle  = "Cleaned_Muscle_Proteomics.csv", 
  Serum   = "Cleaned_Serum_Proteomics.csv"
)

# ==============================================================================
# --- 2. PLOT DIMENSIONS & AESTHETICS CONFIGURATION ---
# ==============================================================================
# Bubble Plot Dimensions (inches)
BUBBLE_PDF_WIDTH  <- 6.5
BUBBLE_PDF_HEIGHT <- 5.5

# Heatmap Cell Dimensions (mm)
HM_CELL_WIDTH  <- 4  
HM_CELL_HEIGHT <- 4  

# Heatmap Base Dimensions (inches)
HM_PDF_BASE_WIDTH  <- 4.5 
HM_PDF_BASE_HEIGHT <- 2.5 

# Top Significant Proteins per Tissue
HM_TOP_N <- 8
# ==============================================================================

# --- 3. METADATA LOADING (V2.0 Strict Engine) ---
message(">>> Loading Metadata...")
clinical <- read_csv(file.path(DIR_IN, "Clinical_Master_Strict.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(
    Role_Orig = dplyr::case_when(Membercode == 1 ~ "Daughter", Membercode == 2 ~ "Mother", Membercode == 3 ~ "Father", TRUE ~ "Unknown"),
    Sample_Key = tolower(paste0(FamilyID, "_", Role_Orig))
  )

meta_merged <- clinical %>%
  dplyr::rename(CRP = CRP) %>% # Explicitly map CRP if needed based on strict columns
  dplyr::filter(!is.na(CRP)) %>%
  dplyr::mutate(
    Log_CRP = log(CRP), 
    Sex = factor(Gender), 
    Age = as.numeric(Age),
    Fat_percent = as.numeric(`Fat(%)`),
    Fat_Group = ifelse(!is.na(Fat_percent) & Fat_percent >= 30, "Obese", "Lean")
  )

parse_sample_v2 <- function(col_name) {
  x_clean <- tolower(str_remove(col_name, "_twin\\.\\d+$")) %>% str_trim()
  if(str_detect(x_clean, "post")) return(NULL) # Baseline only
  s_key <- gsub("_fast|_pre|_post1h|_post3h", "", x_clean)
  
  idx <- which(meta_merged$Sample_Key == s_key)
  if(length(idx) != 1) return(NULL)
  
  list(OriginalCol=col_name, UniqueID=s_key)
}

# --- 4. MASTER PIPELINE WITH CACHE SYSTEM ---
results_store <- list()
for(tissue in names(TARGETS)) {
  message(paste0("\n>>> Processing Tissue: ", tissue))
  file_limma <- file.path(TABLE_DIR, paste0("Stats_Limma_", tissue, ".xlsx"))
  
  fpath <- file.path(DIR_IN, TARGETS[[tissue]])
  if(!file.exists(fpath)) { message("  [Skip] File not found: ", fpath); next }
  
  raw <- read_csv(fpath, show_col_types = FALSE)
  feat_ids <- raw[[1]]; expr_raw <- as.matrix(raw[, -1])
  
  col_infos <- lapply(colnames(expr_raw), parse_sample_v2)
  valid_idx <- which(!sapply(col_infos, is.null))
  if(length(valid_idx) == 0) next
  
  expr_valid <- expr_raw[, valid_idx, drop=FALSE]
  col_meta <- dplyr::bind_rows(col_infos[valid_idx]) %>% 
    dplyr::left_join(meta_merged, by = c("UniqueID" = "Sample_Key")) %>%
    dplyr::distinct(UniqueID, .keep_all = TRUE)
  
  mat <- expr_valid[, col_meta$OriginalCol, drop=FALSE]
  rownames(mat) <- feat_ids
  mat <- mat[rowSums(is.na(mat)) < (ncol(mat)*0.5), ] 
  mat[is.na(mat)] <- min(mat, na.rm=T)/2
  if(max(mat, na.rm=T) > 50) mat <- log2(mat + 1)
  colnames(mat) <- col_meta$UniqueID
  
  if(file.exists(file_limma)) {
    message("    [CACHE HIT] Loaded existing Limma results.")
    res <- read.xlsx(file_limma)
  } else {
    message("    [CALCULATING] Running Limma model...")
    design <- model.matrix(~ Log_CRP + Age + Sex, data = col_meta)
    fit <- eBayes(lmFit(mat, design))
    res <- topTable(fit, coef="Log_CRP", number=Inf) %>% tibble::rownames_to_column("Protein") %>% dplyr::mutate(Tissue = tissue)
    write.xlsx(res, file_limma)
  }
  
  results_store[[tissue]] <- list(res = res, mat = mat, meta = col_meta)
}

# ==============================================================================
# PART 1: THEMED GSEA LANDSCAPE (Panel D)
# ==============================================================================
message("\n>>> Generating Panel D: Themed GSEA Landscape...")

gsea_list <- list()
for(t in names(results_store)) {
  file_gsea <- file.path(TABLE_DIR, paste0("Stats_GSEA_", t, ".xlsx"))
  if(file.exists(file_gsea)) {
    sig_paths <- read.xlsx(file_gsea)
    gsea_list[[t]] <- sig_paths
  } else {
    res <- results_store[[t]]$res
    ids <- tryCatch(bitr(res$Protein, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"), error=function(e) NULL)
    if(is.null(ids)) next
    
    res_mapped <- dplyr::inner_join(res, ids, by=c("Protein"="SYMBOL")) %>% dplyr::arrange(desc(t)) %>% dplyr::filter(!duplicated(ENTREZID))
    gene_list <- res_mapped$t; names(gene_list) <- res_mapped$ENTREZID
    try({
      gse <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, ont = "BP", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, verbose = FALSE, seed = 123)
      if(!is.null(gse) && nrow(as.data.frame(gse)) > 0) {
        gse_res <- gse@result %>% dplyr::mutate(Tissue = t)
        sig_paths <- gse_res %>% dplyr::filter(p.adjust < 0.25)
        if(nrow(sig_paths) < 3) sig_paths <- gse_res %>% dplyr::filter(pvalue < 0.05)
        if(nrow(sig_paths) < 3) sig_paths <- gse_res %>% dplyr::arrange(pvalue) %>% dplyr::slice_head(n = 20)
        write.xlsx(sig_paths, file_gsea)
        gsea_list[[t]] <- sig_paths
      }
    })
  }
}

if(length(gsea_list) > 0) {
  all_gsea <- dplyr::bind_rows(gsea_list)
  theme_dict <- list(
    "Energy &\nMitochondria" = c("mitochondrial", "respiration", "ATP synthesis", "oxidative", "electron transport"),
    "Immune &\nInflammation" = c("immune", "inflammatory", "complement", "cytokine", "leukocyte", "chemotaxis", "proteolysis"),
    "Structure &\nTransport" = c("adhesion", "integrin", "lipid transport", "organic hydroxy", "cell-substrate")
  )
  
  categorize_pathway <- function(desc) {
    for(thm in names(theme_dict)) { if(any(str_detect(tolower(desc), str_replace_all(theme_dict[[thm]], "\n", "")))) return(thm) }
    return("Other")
  }
  all_gsea$Theme <- sapply(all_gsea$Description, categorize_pathway)
  
  plot_data <- all_gsea %>% dplyr::filter(Theme != "Other") %>% dplyr::group_by(Tissue, Theme) %>% dplyr::arrange(pvalue) %>% dplyr::slice_head(n = 3) %>% dplyr::ungroup() 
  plot_data$Description <- str_wrap(str_to_sentence(str_remove(plot_data$Description, "regulation of ")), width = 35)
  plot_data$Tissue <- factor(plot_data$Tissue, levels = c("Serum", "Adipose", "Muscle"))
  plot_data$Theme <- factor(plot_data$Theme, levels = names(theme_dict))
  plot_data <- plot_data %>% dplyr::arrange(Theme, NES) %>% dplyr::mutate(Description = factor(Description, levels = unique(Description)))
  
  write.xlsx(plot_data, file.path(TABLE_DIR, "Table_Panel_D_BubbleData.xlsx"))
  
  p_bubble <- ggplot(plot_data, aes(x = Tissue, y = Description)) +
    geom_point(aes(size = -log10(pvalue), fill = NES), shape = 21, color = "black", stroke = 0.6) +
    scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", midpoint = 0, name = "NES") +
    scale_size_continuous(range = c(2.5, 6), name = "-Log P") +
    scale_x_discrete(drop = FALSE, expand = expansion(mult = 0.15)) + 
    scale_y_discrete(position = "right") + 
    facet_grid(Theme ~ ., scales = "free_y", switch = "y") +
    theme_classic(base_size = 12) +
    labs(title = "Systemic CRP Associations", x = NULL, y = NULL) +
    theme(
      panel.grid.major = element_blank(), panel.border = element_rect(color="black", fill=NA, linewidth=1), 
      axis.line = element_blank(), strip.background = element_blank(),
      strip.text.y.left = element_text(face = "bold", angle = 90, size = 10, hjust = 0.5), 
      strip.placement = "outside",
      plot.title = element_text(face="bold", hjust=0.5, size=13, margin = margin(b = 5)),
      axis.text.x = element_text(face="bold", color="black", angle = 45, hjust = 1),
      axis.text.y.right = element_text(color="black", size=9, hjust=0),
      panel.spacing = unit(0.3, "lines"), 
      plot.margin = margin(t=2, r=2, b=0, l=2, unit="mm"),
      legend.position = "right", 
      legend.box = "vertical",
      legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold"),
      aspect.ratio = 1
    )
  
  ggsave(file.path(DIR_OUT, "Fig5D_GSEA_Landscape.pdf"), p_bubble, width = BUBBLE_PDF_WIDTH, height = BUBBLE_PDF_HEIGHT)
  message("  -> Panel D (GSEA Bubble) Saved.")
}

# ==============================================================================
# PART 2: SQUARE-CELL HEATMAP (Panel E: Strict Intersection)
# ==============================================================================
message("\n>>> Generating Panel E: Clustered Strict-Intersection Heatmap...")

# 1. Extract VIP Subjects (Must have data across all 3 tissues)
subs_adi <- if(!is.null(results_store[["Adipose"]])) colnames(results_store[["Adipose"]]$mat) else character(0)
subs_mus <- if(!is.null(results_store[["Muscle"]])) colnames(results_store[["Muscle"]]$mat) else character(0)
subs_ser <- if(!is.null(results_store[["Serum"]])) colnames(results_store[["Serum"]]$mat) else character(0)

valid_tissue_subs <- intersect(intersect(subs_adi, subs_mus), subs_ser)

# 2. Strict Subject Ordering from metadata
all_subjects <- meta_merged %>% 
  dplyr::filter(Sample_Key %in% valid_tissue_subs) %>% 
  dplyr::arrange(factor(Fat_Group, levels=c("Lean", "Obese")), CRP) %>% 
  dplyr::pull(Sample_Key)

message(sprintf("   -> Strict Intersection (Adi + Mus + Ser): %d subjects retained.", length(all_subjects)))

mat_list <- list(); anno_list <- list()

for(t in c("Muscle", "Serum", "Adipose")) {
  if(!(t %in% names(results_store))) next
  dat <- results_store[[t]]
  
  sig_genes <- dat$res %>% dplyr::filter(P.Value < 0.05) %>% dplyr::arrange(P.Value) %>% head(HM_TOP_N)
  if(nrow(sig_genes) == 0) next
  
  mat_sub <- dat$mat[sig_genes$Protein, all_subjects, drop=FALSE]
  
  # --- Robust Z-Score & Winsorization ---
  z_mat <- t(apply(mat_sub, 1, function(x) {
    if(sum(!is.na(x)) > 1) {
      if (mad(x, na.rm=T) > 0) { z <- (x - median(x, na.rm=T)) / mad(x, na.rm=T) }
      else if (sd(x, na.rm=T) > 0) { z <- (x - mean(x, na.rm=T)) / sd(x, na.rm=T) }
      else { z <- rep(0, length(x)) }
    } else { z <- rep(0, length(x)) }
    return(z)
  }))
  
  z_mat[z_mat > 2.5] <- 2.5
  z_mat[z_mat < -2.5] <- -2.5
  # --------------------------------------
  
  mat_list[[t]] <- z_mat
  
  for(i in 1:nrow(sig_genes)) {
    clean_name <- str_replace(sig_genes$Protein[i], ";.*", "")
    anno_list[[length(anno_list)+1]] <- data.frame(
      Feature = sig_genes$Protein[i], 
      Tissue = t, 
      logFC = sig_genes$logFC[i],
      Display_Name = clean_name, 
      stringsAsFactors = FALSE
    )
  }
}

if(length(mat_list) > 0) {
  unified_mat <- do.call(rbind, mat_list)
  row_anno_df <- do.call(rbind, anno_list)
  rownames(unified_mat) <- row_anno_df$Display_Name
  
  valid_cols <- apply(unified_mat, 2, function(x) sum(!is.na(x)) > 0)
  unified_mat <- unified_mat[, valid_cols]
  final_subs <- colnames(unified_mat)
  final_meta <- meta_merged[match(final_subs, meta_merged$Sample_Key), ]
  final_meta$Fat_Group <- factor(final_meta$Fat_Group, levels = c("Lean", "Obese"))
  
  # --- Annotations ---
  col_crp <- colorRamp2(c(min(final_meta$CRP, na.rm=T), median(final_meta$CRP, na.rm=T), max(final_meta$CRP, na.rm=T)), c("white", "#FFC107", "#DC0000"))
  top_anno <- HeatmapAnnotation(
    Fat_Group = final_meta$Fat_Group,
    CRP = final_meta$CRP,
    col = list(Fat_Group = c("Obese" = "#F79647", "Lean" = "#3DA6AE"), CRP = col_crp),
    annotation_label = c(Fat_Group = "Fat(%)", CRP = "CRP"),
    show_annotation_name = TRUE, annotation_name_side = "left", simple_anno_size = unit(3, "mm")
  )
  
  row_anno_df$Tissue <- factor(row_anno_df$Tissue, levels = c("Serum", "Adipose", "Muscle"))
  tissue_colors <- c("Adipose"="#F39B7F", "Muscle"="#4DBBD5", "Serum"="#00A087")
  
  left_anno <- rowAnnotation(
    Tissue = row_anno_df$Tissue,
    Beta = anno_barplot(row_anno_df$logFC, bar_width = 1, gp = gpar(fill = ifelse(row_anno_df$logFC > 0, "#E64B35", "#3C5488"), col=NA)),
    col = list(Tissue = tissue_colors),
    show_annotation_name = TRUE, annotation_name_rot = 45
  )
  
  # --- Heatmap Plotting ---
  ht <- Heatmap(
    unified_mat, 
    name = "Z-Score",
    col = colorRamp2(c(-1.5, 0, 1.5), c("#3C5488", "white", "#E64B35")),
    na_col = "#E8E8E8", 
    top_annotation = top_anno,
    left_annotation = left_anno,
    column_split = final_meta$Fat_Group,
    row_split = row_anno_df$Tissue,
    cluster_columns = FALSE, 
    cluster_rows = TRUE, 
    cluster_row_slices = FALSE, 
    column_gap = unit(2, "mm"), row_gap = unit(2, "mm"),
    width = ncol(unified_mat) * unit(HM_CELL_WIDTH, "mm"),
    height = nrow(unified_mat) * unit(HM_CELL_HEIGHT, "mm"),
    show_column_names = FALSE,
    row_names_side = "right",
    row_names_gp = gpar(fontsize = 10, fontface="bold.italic", col="black"),
    rect_gp = gpar(col = "white", lwd = 1), 
    border = TRUE,
    column_title = paste0("Strict Intersection Response (N=", ncol(unified_mat), ")"),
    column_title_gp = gpar(fontsize=12, fontface="bold")
  )
  
  final_pdf_w <- HM_PDF_BASE_WIDTH + (ncol(unified_mat) * HM_CELL_WIDTH / 25.4)
  final_pdf_h <- HM_PDF_BASE_HEIGHT + (nrow(unified_mat) * HM_CELL_HEIGHT / 25.4)
  
  pdf(file.path(DIR_OUT, "Fig5E_StrictIntersection_Heatmap.pdf"), width = final_pdf_w, height = final_pdf_h)
  draw(ht, merge_legend = TRUE)
  dev.off()
  message(sprintf("   -> Panel E (Heatmap) Saved: %.2f x %.2f inches", final_pdf_w, final_pdf_h))
}

# ==============================================================================
# PART 3: GENERATING PUBLICATION-READY SUPPLEMENTARY TABLE (Data_S9)
# ==============================================================================
message("\n>>> Generating Publication-Ready Supplementary Table (Data_S9)...")

wb_sup <- createWorkbook()
header_style <- createStyle(fontName = "Arial", fontSize = 11, fontColour = "white", fgFill = "#4F81BD", textDecoration = "bold")
sig_red      <- createStyle(fontColour = "#B2182B", textDecoration = "bold")
sig_blue     <- createStyle(fontColour = "#2166AC", textDecoration = "bold")

# --- Sheet 1-3: Limma Associations ---
for(t_name in names(results_store)) {
  df_limma <- results_store[[t_name]]$res %>%
    dplyr::mutate(
      Significance = dplyr::case_when(P.Value < 0.001 ~ "***", P.Value < 0.01 ~ "**", P.Value < 0.05 ~ "*", TRUE ~ "ns"),
      Clean_Protein = str_replace(Protein, ";.*", "") 
    ) %>%
    dplyr::arrange(P.Value) %>%
    dplyr::select(Tissue, Feature = Protein, Clean_Feature = Clean_Protein, 
                  Association_Coefficient = logFC, t_statistic = t, 
                  P_Value = P.Value, adj.P_Val_FDR = adj.P.Val, Significance)
  
  sheet_name <- paste0("Limma_", t_name)
  addWorksheet(wb_sup, sheet_name)
  writeData(wb_sup, sheet_name, df_limma)
  addStyle(wb_sup, sheet_name, style = header_style, rows = 1, cols = 1:ncol(df_limma), gridExpand = TRUE)
  conditionalFormatting(wb_sup, sheet_name, cols = 4, rows = 2:(nrow(df_limma)+1), rule = ">0", style = sig_red)
  conditionalFormatting(wb_sup, sheet_name, cols = 4, rows = 2:(nrow(df_limma)+1), rule = "<0", style = sig_blue)
  conditionalFormatting(wb_sup, sheet_name, cols = 8, rows = 2:(nrow(df_limma)+1), rule = "*", type = "contains", style = sig_red)
}

# --- Sheet 4-6: GSEA Enrichment Results ---
for(t_name in names(gsea_list)) {
  df_gsea <- gsea_list[[t_name]] %>%
    dplyr::mutate(
      Tissue = t_name, 
      Theme = sapply(Description, categorize_pathway) 
    ) %>%
    dplyr::arrange(pvalue) %>%
    dplyr::select(Tissue, Theme, ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, core_enrichment)
  
  sheet_name <- paste0("GSEA_", t_name)
  addWorksheet(wb_sup, sheet_name)
  writeData(wb_sup, sheet_name, df_gsea)
  addStyle(wb_sup, sheet_name, style = header_style, rows = 1, cols = 1:ncol(df_gsea), gridExpand = TRUE)
}

# --- Sheet 7: Heatmap Z-Score Tidy Matrix ---
if(exists("unified_mat") && nrow(unified_mat) > 0) {
  tidy_hm <- as.data.frame(unified_mat) %>%
    tibble::rownames_to_column("Feature") %>%
    tidyr::pivot_longer(cols = -Feature, names_to = "Subject_ID", values_to = "Z_Score") %>%
    dplyr::left_join(row_anno_df %>% dplyr::select(Display_Name, Tissue, Association_Coefficient = logFC), by = c("Feature" = "Display_Name")) %>%
    dplyr::left_join(final_meta %>% dplyr::select(Sample_Key, Fat_Group, CRP, Log_CRP), by = c("Subject_ID" = "Sample_Key")) %>%
    dplyr::select(Tissue, Feature, Subject_ID, Fat_Group, CRP, Log_CRP, Association_Coefficient, Z_Score) %>%
    dplyr::arrange(factor(Tissue, levels = c("Serum", "Adipose", "Muscle")), Feature, factor(Fat_Group, levels = c("Lean", "Obese")), CRP)
  
  addWorksheet(wb_sup, "Heatmap_Z_Scores")
  writeData(wb_sup, "Heatmap_Z_Scores", tidy_hm)
  addStyle(wb_sup, "Heatmap_Z_Scores", style = header_style, rows = 1, cols = 1:ncol(tidy_hm), gridExpand = TRUE)
}

# --- Save Final Table ---
sup_out_name <- file.path(DIR_OUT, "Data_S9_CRP_Integrated_Downstream.xlsx")
saveWorkbook(wb_sup, sup_out_name, overwrite = TRUE)

message("\n=================================================================")
message("🏁 FIG 5 D&E AND DATA_S9 EXPORTED SUCCESSFULLY!")
message(paste("Check your files at:", DIR_OUT))
message("=================================================================")
