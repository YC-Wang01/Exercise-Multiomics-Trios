# ==============================================================================
# Project: FinlandSports V2.0 - Fig S4 Continuous GSEA Comprehensive
# Description: 
#   1. Continuous GSEA: Protein vs HOMA-IR (Baseline & Response).
#   2. Quadrant Assembly: Merges Muscle & Adipose NES to create Quadrant Data.
#   3. Quadrant Plot: Generates Anti-Crowding quadrant scatter plots.
#   4. Smart Selection: Ranks high-value pathways and plots GSEA line charts.
#   5. Cnetplot: Generates circular Gene-Concept network for target pathway.
# Features: 100% English, V2.0 Paths, Strict dplyr:: Namespacing, No Missing Data.
# ==============================================================================

# --- 0. Environment Setup ---
setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")
options(timeout = 600) # Prevent KEGG API timeout

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readr, readxl, dplyr, stringr, tibble, ggplot2, openxlsx, 
  clusterProfiler, org.Hs.eg.db, enrichplot, cowplot, grid, ggrepel
)

DIR_DATA <- "01_Clean_Data"
DIR_OUT  <- "03_Results/Fig_S4/Continuous_GSEA"
if(!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)

message(">>> Directories Ready. Output set to: ", DIR_OUT)

# --- 1. Load Clinical Data ---
message("\n=== STEP 1: Loading Clinical Data & Calculating HOMA-IR ===")

clinical <- read_csv(file.path(DIR_DATA, "Clinical_Master_Strict.csv"), show_col_types = FALSE) %>% 
  dplyr::mutate(
    Role_Orig = dplyr::case_when(Membercode == 1 ~ "Daughter", Membercode == 2 ~ "Mother", Membercode == 3 ~ "Father", TRUE ~ "Unknown"),
    Subject_ID = tolower(paste0(FamilyID, "_", Role_Orig))
  ) %>% 
  dplyr::filter(!is.na(HOMA_IR)) %>%
  dplyr::distinct(Subject_ID, .keep_all = TRUE)

# ==============================================================================
# === STEP 2: Continuous Correlation & GSEA Execution ===
# ==============================================================================
message("\n=== STEP 2: Running Continuous GSEA (Protein vs HOMA-IR) ===")

GSEA_OBJS <- list()
GSEA_GENELISTS <- list()
GSEA_RESULTS_DF <- data.frame()

run_tissue_gsea <- function(filename, tissue_name) {
  fpath <- file.path(DIR_DATA, filename)
  if(!file.exists(fpath)) return(NULL)
  
  raw <- read_csv(fpath, show_col_types = FALSE)
  mat <- as.matrix(raw[, -1]); rownames(mat) <- raw[[1]]
  
  clean_cols <- tolower(str_remove(colnames(mat), "_twin\\.\\d+$")) %>% str_trim()
  meta <- data.frame(Matrix_Col = colnames(mat), Clean_Col = clean_cols, stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      TimePoint = dplyr::case_when(grepl("pre", Clean_Col) ~ "pre", grepl("fast", Clean_Col) ~ "fast", grepl("post3h", Clean_Col) ~ "post3h", TRUE ~ "Other"),
      Subject_ID = gsub("_fast|_pre|_post1h|_post3h", "", Clean_Col)
    ) %>%
    dplyr::inner_join(clinical[, c("Subject_ID", "HOMA_IR")], by = "Subject_ID")
  
  mat_sub <- mat[, meta$Matrix_Col, drop=FALSE]
  if(max(mat_sub, na.rm=T) > 50) mat_sub <- log2(mat_sub + 1)
  
  run_gsea_core <- function(gene_cor_vec, analysis_name) {
    if(length(gene_cor_vec) < 10) return(NULL)
    gene_cor_vec <- sort(gene_cor_vec, decreasing = TRUE)
    
    ids <- tryCatch(bitr(names(gene_cor_vec), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db), error=function(e) NULL)
    if(is.null(ids)) return(NULL)
    
    gene_list_entrez <- gene_cor_vec[names(gene_cor_vec) %in% ids$SYMBOL]
    names(gene_list_entrez) <- ids$ENTREZID[match(names(gene_list_entrez), ids$SYMBOL)]
    
    gsea_res <- tryCatch(gseKEGG(geneList=gene_list_entrez, organism='hsa', pvalueCutoff=1, verbose=FALSE, seed=123), error=function(e) NULL)
    if(is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) return(NULL)
    
    gsea_res <- setReadable(gsea_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    res_df <- as.data.frame(gsea_res)
    res_df$Analysis <- analysis_name
    
    # Store globally for later plots
    GSEA_OBJS[[analysis_name]] <<- gsea_res
    GSEA_GENELISTS[[analysis_name]] <<- gene_cor_vec # Keep original Symbol version for cnetplot
    GSEA_RESULTS_DF <<- rbind(GSEA_RESULTS_DF, res_df)
  }
  
  # A. Baseline
  meta_base <- meta %>% dplyr::filter(TimePoint %in% c("pre", "fast")) %>% dplyr::distinct(Subject_ID, .keep_all = TRUE)
  if(nrow(meta_base) >= 6) {
    mat_base <- mat_sub[, meta_base$Matrix_Col, drop=FALSE]
    cor_vals <- apply(mat_base, 1, function(y) {
      if(sum(!is.na(y)) < 3 || sd(y,na.rm=T)==0) return(NA)
      suppressWarnings(cor(y, meta_base$HOMA_IR, use="pairwise.complete.obs"))
    })
    cor_vals <- cor_vals[!is.na(cor_vals)]
    message("   -> Running GSEA: ", tissue_name, " Baseline...")
    run_gsea_core(cor_vals, paste0(tissue_name, "_Baseline"))
  }
  
  # B. Response (Delta)
  subjects_paired <- intersect(meta$Subject_ID[meta$TimePoint %in% c("pre", "fast")], meta$Subject_ID[meta$TimePoint == "post3h"])
  if(length(subjects_paired) >= 6) {
    delta_list <- list(); homa_list <- c()
    for(sub in subjects_paired) {
      col_pre  <- meta$Matrix_Col[meta$Subject_ID == sub & meta$TimePoint %in% c("pre", "fast")][1]
      col_post <- meta$Matrix_Col[meta$Subject_ID == sub & meta$TimePoint == "post3h"][1]
      delta_list[[sub]] <- mat_sub[, col_post] - mat_sub[, col_pre]
      homa_list <- c(homa_list, meta$HOMA_IR[meta$Subject_ID == sub][1])
    }
    mat_delta <- do.call(cbind, delta_list)
    cor_vals_d <- apply(mat_delta, 1, function(y) {
      if(sum(!is.na(y)) < 3 || sd(y,na.rm=T)==0) return(NA)
      suppressWarnings(cor(y, homa_list, use="pairwise.complete.obs"))
    })
    cor_vals_d <- cor_vals_d[!is.na(cor_vals_d)]
    message("   -> Running GSEA: ", tissue_name, " Response...")
    run_gsea_core(cor_vals_d, paste0(tissue_name, "_Response"))
  }
}

run_tissue_gsea("Cleaned_Muscle_Proteomics.csv", "Muscle")
run_tissue_gsea("Cleaned_Adipose_Proteomics.csv", "Adipose")

write.xlsx(GSEA_RESULTS_DF, file.path(DIR_OUT, "Table_S_Continuous_GSEA_All_Results.xlsx"))
message("  > Global GSEA results saved.")

# ==============================================================================
# === STEP 3: Quadrant Assembly & Anti-Crowding Plot ===
# ==============================================================================
message("\n=== STEP 3: Assembling Quadrant Data & Generating Scatter Plots ===")

build_and_plot_quadrant <- function(mode_label) {
  df_mode <- GSEA_RESULTS_DF %>% dplyr::filter(grepl(mode_label, Analysis))
  if(nrow(df_mode) == 0) return(NULL)
  
  mus <- df_mode %>% dplyr::filter(grepl("Muscle", Analysis)) %>% dplyr::select(ID, Description, NES_Mus=NES, FDR_Mus=p.adjust)
  adi <- df_mode %>% dplyr::filter(grepl("Adipose", Analysis)) %>% dplyr::select(ID, Description, NES_Adi=NES, FDR_Adi=p.adjust)
  
  merged <- dplyr::full_join(mus, adi, by = c("ID", "Description")) %>%
    dplyr::mutate(
      NES_Mus = ifelse(is.na(NES_Mus), 0, NES_Mus),
      NES_Adi = ifelse(is.na(NES_Adi), 0, NES_Adi),
      FDR_Mus = ifelse(is.na(FDR_Mus), 1, FDR_Mus),
      FDR_Adi = ifelse(is.na(FDR_Adi), 1, FDR_Adi),
      Is_Significant = (FDR_Mus < 0.05 | FDR_Adi < 0.05),
      Quadrant = dplyr::case_when(
        !Is_Significant ~ "Not Significant",
        NES_Mus > 0 & NES_Adi > 0 ~ "Q1: Co-Activation",
        NES_Mus < 0 & NES_Adi > 0 ~ "Q2: Adi Up / Mus Down",
        NES_Mus < 0 & NES_Adi < 0 ~ "Q3: Co-Suppression",
        NES_Mus > 0 & NES_Adi < 0 ~ "Q4: Mus Up / Adi Down",
        TRUE ~ "Not Significant"
      )
    )
  
  out_name <- paste0("Table_S_FigS4_Quadrant_", mode_label, "_Refined_Data.xlsx")
  write.xlsx(merged, file.path(DIR_OUT, out_name))
  message("  > Generated Quadrant Data: ", out_name)
  
  # --- Plotting Engine ---
  KEYWORDS_HOMA <- c("Insulin", "Diabetes", "Diabetic", "Adipocytokine", "Fatty acid", "Lipid", "PPAR", "AMPK", "Gluconeogenesis", "Glycolysis", "Glucagon", "Oxidative phosphorylation", "Citrate cycle", "TCA", "Mitochondria", "Thermogenesis")
  KEYWORDS_IGNORE <- c("Carbon metabolism", "Metabolic pathways", "Biosynthesis of amino acids", "Propanoate metabolism", "Butanoate metabolism", "Pyruvate metabolism", "Valine", "2-Oxocarboxylic acid")
  
  pattern_homa <- paste(KEYWORDS_HOMA, collapse = "|")
  pattern_ignore <- paste(KEYWORDS_IGNORE, collapse = "|")
  
  df_plot <- merged %>%
    dplyr::mutate(
      Is_HOMA = grepl(pattern_homa, Description, ignore.case = TRUE),
      Is_Ignored = grepl(pattern_ignore, Description, ignore.case = TRUE),
      Vector_Strength = sqrt(NES_Mus^2 + NES_Adi^2)
    )
  
  labels_homa <- df_plot %>% dplyr::filter(Is_Significant & Is_HOMA & !Is_Ignored) %>% dplyr::group_by(Quadrant) %>% dplyr::arrange(desc(Vector_Strength)) %>% dplyr::slice_head(n = 6) %>% dplyr::pull(Description)
  labels_top  <- df_plot %>% dplyr::filter(Is_Significant & !Is_HOMA & !Is_Ignored) %>% dplyr::group_by(Quadrant) %>% dplyr::arrange(desc(Vector_Strength)) %>% dplyr::slice_head(n = 4) %>% dplyr::pull(Description)
  
  df_plot <- df_plot %>%
    dplyr::mutate(
      Label_Type = dplyr::case_when(Description %in% labels_homa ~ "HOMA", Description %in% labels_top ~ "Top5", TRUE ~ "None"),
      Label_Formatted = ifelse(Label_Type != "None", str_wrap(Description, 20), NA)
    )
  
  max_val <- max(abs(c(df_plot$NES_Mus, df_plot$NES_Adi)), na.rm = TRUE) * 1.15
  limit_range <- c(-max_val, max_val)
  my_colors <- c("Q1: Co-Activation" = "#d73027", "Q2: Adi Up / Mus Down" = "#fdae61", "Q3: Co-Suppression" = "#4575b4", "Q4: Mus Up / Adi Down" = "#abd9e9", "Not Significant" = "grey95")
  
  p <- ggplot(df_plot, aes(x = NES_Mus, y = NES_Adi)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") + geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(data = subset(df_plot, !Is_Significant), color = "grey92", size = 3.5, alpha = 0.6) +
    geom_point(data = subset(df_plot, Is_Significant), aes(color = Quadrant, size = Vector_Strength), alpha = 0.85, stroke = 0.2) +
    geom_text_repel(data = subset(df_plot, Label_Type == "Top5"), aes(label = Label_Formatted), color = "black", size = 3.8, fontface = "italic", bg.color = "white", bg.r = 0.15, box.padding = 0.6, force = 5, seed = 1) +
    geom_text_repel(data = subset(df_plot, Label_Type == "HOMA"), aes(label = Label_Formatted), color = "#D00000", size = 4.2, fontface = "bold", bg.color = "white", bg.r = 0.15, box.padding = 0.8, force = 10, seed = 2, min.segment.length = 0) +
    scale_color_manual(values = my_colors) + scale_size_continuous(range = c(5, 12), guide = "none") + scale_x_continuous(limits = limit_range) + scale_y_continuous(limits = limit_range) +
    labs(title = paste(mode_label, "Crosstalk: Functional Shift"), subtitle = "Red = Metabolic/HOMA (Top 6 per Quad) | Black = Top 4 Others", x = "Muscle NES (Correlation with HOMA)", y = "Adipose NES (Correlation with HOMA)", color = NULL) +
    theme_bw(base_size = 15) + theme(panel.grid = element_blank(), legend.position = "bottom", plot.title = element_text(face = "bold", size = 18), aspect.ratio = 1)
  
  ggsave(file.path(DIR_OUT, paste0("FigS4_Quadrant_", mode_label, "_Final.pdf")), p, width = 10, height = 10)
  message("  > Generated Quadrant Plot for: ", mode_label)
}

build_and_plot_quadrant("Baseline")
build_and_plot_quadrant("Response")

# ==============================================================================
# === STEP 4: Smart Selection & Supplementary GSEA Plots ===
# ==============================================================================
message("\n=== STEP 4: Smart Selection & GSEA Line Plots ===")

KEYWORDS_HIGH_VALUE <- c("Oxidative phosphorylation", "TCA cycle", "Citrate cycle", "Insulin", "Diabetes", "Diabetic", "Cardiomyopathy", "Fatty acid metabolism", "Thermogenesis", "PPAR", "AMPK", "AGE-RAGE", "Alzheimer", "Inflammation")

rank_pathways <- function(mode_select) {
  df_mode <- GSEA_RESULTS_DF %>% dplyr::filter(grepl(mode_select, Analysis))
  if(nrow(df_mode)==0) return(data.frame())
  
  mus <- df_mode %>% dplyr::filter(grepl("Muscle", Analysis)) %>% dplyr::select(ID, Description, NES_Mus=NES, FDR_Mus=p.adjust)
  adi <- df_mode %>% dplyr::filter(grepl("Adipose", Analysis)) %>% dplyr::select(ID, Description, NES_Adi=NES, FDR_Adi=p.adjust)
  
  merged <- dplyr::inner_join(mus, adi, by = c("ID", "Description")) %>%
    dplyr::mutate(
      Is_Key = grepl(paste(KEYWORDS_HIGH_VALUE, collapse="|"), Description, ignore.case=TRUE),
      Is_Opposite = sign(NES_Mus) != sign(NES_Adi),
      Sig_Score = -log10(FDR_Mus) - log10(FDR_Adi),
      Total_Score = (as.numeric(Is_Key) * 100) + (as.numeric(Is_Opposite) * 50) + Sig_Score
    ) %>% dplyr::arrange(desc(Total_Score)) %>% head(6)
  return(merged)
}

df_base_top <- rank_pathways("Baseline") %>% dplyr::mutate(Select_Group = "Baseline Top")
df_resp_top <- rank_pathways("Response") %>% dplyr::mutate(Select_Group = "Response Top")
final_selection <- rbind(df_base_top, df_resp_top)

if(nrow(final_selection) > 0) {
  write.xlsx(final_selection, file.path(DIR_OUT, "Table_S_Smart_Selected_Pathways.xlsx"))
  target_ids <- unique(final_selection$ID)
  
  make_grob <- function(obj_name, tid, title) {
    obj <- GSEA_OBJS[[obj_name]]
    if(is.null(obj) || !tid %in% obj$ID) return(grid::grid.grabExpr(print(ggplot()+theme_void()+labs(title=paste(title, "\n(Not Significant)")))))
    p <- gseaplot2(obj, geneSetID=tid, title=title, color="#4daf4a", base_size=10, rel_heights=c(1.5, 0.5, 1), subplots=1:2)
    return(grid::grid.grabExpr(print(p)))
  }
  
  for(tid in target_ids) {
    desc <- final_selection %>% dplyr::filter(ID == tid) %>% dplyr::pull(Description) %>% head(1)
    safe_name <- gsub("[^A-Za-z0-9]", "_", str_trunc(desc, 30))
    
    g1 <- make_grob("Muscle_Baseline", tid, paste0("Mus Base\n", desc))
    g2 <- make_grob("Adipose_Baseline", tid, paste0("Adi Base\n", desc))
    g3 <- make_grob("Muscle_Response", tid, paste0("Mus Resp\n", desc))
    g4 <- make_grob("Adipose_Response", tid, paste0("Adi Resp\n", desc))
    
    final_plot <- plot_grid(g1, g2, g3, g4, ncol=2, labels="AUTO")
    ggsave(file.path(DIR_OUT, paste0("FigS4_GSEA_LinePlot_", safe_name, ".pdf")), final_plot, width=12, height=10)
  }
  message("  > Saved High-Value GSEA Line Plots.")
}

# ==============================================================================
# === STEP 5: Circular Gene-Concept Network (Cnetplot) [VERSION FIXED] ===
# ==============================================================================
message("\n=== STEP 5: Generating Circular Cnetplot for Oxidative Phosphorylation ===")

TARGET_ID <- "hsa00190"  # Oxidative phosphorylation

make_cnet_plot <- function(analysis_name, title_text) {
  obj <- GSEA_OBJS[[analysis_name]]
  fc <- GSEA_GENELISTS[[analysis_name]] # Original correlations
  
  if(is.null(obj) || !TARGET_ID %in% obj$ID) {
    return(ggplot() + theme_void() + labs(title=paste(title_text, "\n(Pathway Not Found)")))
  }
  target_desc <- obj@result$Description[obj@result$ID == TARGET_ID]
  
  # 【核心修正】：使用安全参数 layout = "circle" 完美替代被废弃的 circular = TRUE
  p <- cnetplot(obj, showCategory = target_desc, foldChange = fc, layout = "circle") +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", name = "Corr") +
    labs(title = title_text) + 
    theme(legend.position = "right")
  
  return(p)
}

p1 <- make_cnet_plot("Muscle_Baseline", "A. Muscle Baseline (Pre)")
p2 <- make_cnet_plot("Adipose_Baseline", "B. Adipose Baseline (Pre)")
p3 <- make_cnet_plot("Muscle_Response", "C. Muscle Response (Delta)")
p4 <- make_cnet_plot("Adipose_Response", "D. Adipose Response (Delta)")

g1 <- grid::grid.grabExpr(print(p1))
g2 <- grid::grid.grabExpr(print(p2))
g3 <- grid::grid.grabExpr(print(p3))
g4 <- grid::grid.grabExpr(print(p4))

final_cnet <- cowplot::plot_grid(g1, g2, g3, g4, ncol = 2)
ggsave(file.path(DIR_OUT, paste0("FigS4_Cnetplot_", TARGET_ID, "_Circular.pdf")), final_cnet, width = 16, height = 14)

message("  > Saved Circular Cnetplot.")
message("\n=======================================================")
message(">>> FIG S4 PIPELINE COMPLETED! All outputs saved successfully.")
message("=======================================================")