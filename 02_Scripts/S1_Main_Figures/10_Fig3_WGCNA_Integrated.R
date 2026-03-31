# ==============================================================================
# Project: FinlandSports V2.0 
# Script:  10_Fig3_WGCNA_Integrated.R
# Panels:  WGCNA Bubble Plots & Systemic Inflammation (CRP) Correlation Barplots
# Core:    Unified R Pipeline (Delta Calc -> WGCNA -> GO Enrich -> Plotting)
# Features: 100% English, Strict dplyr:: Namespacing, Unified Excel Export
# ==============================================================================

# ------------------------------------------------------------------------------
# [0] Environment & Path Initialization
# ------------------------------------------------------------------------------
# setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readr, dplyr, stringr, tibble, ggplot2, WGCNA, 
  clusterProfiler, org.Hs.eg.db, reshape2, openxlsx
)

enableWGCNAThreads()
cor <- WGCNA::cor # Force WGCNA correlation function to prevent conflicts

DIR_DATA <- "01_Clean_Data"
DIR_OUT  <- "03_Results/Fig_3/WGCNA_Analysis"
if(!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)

SUPP_TABLE_LIST <<- list()

message(">>> Output directory set to: ", DIR_OUT)

# ------------------------------------------------------------------------------
# [1] Data Preparation (In-Memory Delta Matrix Calculation)
# ------------------------------------------------------------------------------
message("\n=== STEP 1: Preparing Clinical Data & Delta Matrices ===")

# Load clinical master
clinical <- read_csv(file.path(DIR_DATA, "Clinical_Master_Strict.csv"), show_col_types = FALSE) %>% 
  dplyr::mutate(
    Role_Orig = dplyr::case_when(Membercode == 1 ~ "Daughter", Membercode == 2 ~ "Mother", Membercode == 3 ~ "Father", TRUE ~ "Unknown"),
    Subject_ID = tolower(paste0(FamilyID, "_", Role_Orig)),
    Fat_percent = as.numeric(`Fat(%)`)
  )

# Calculate Delta Matrix (Post3h - Pre/Fast) directly in R
get_delta_matrix <- function(filename) {
  fpath <- file.path(DIR_DATA, filename)
  if(!file.exists(fpath)) return(NULL)
  
  data_raw <- read_csv(fpath, show_col_types = FALSE)
  feat_ids <- data_raw[[1]]; expr_mat <- as.matrix(data_raw[, -1]); rownames(expr_mat) <- feat_ids
  
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
    
    idx_pre <- which(colnames(expr_mat) %in% pre_ids)
    idx_post <- which(colnames(expr_mat) %in% post_ids)
    
    if(length(idx_pre) > 0 && length(idx_post) > 0) {
      v_pre <- if(length(idx_pre) > 1) rowMeans(expr_mat[, idx_pre, drop=F], na.rm=T) else expr_mat[, idx_pre]
      v_post <- if(length(idx_post) > 1) rowMeans(expr_mat[, idx_post, drop=F], na.rm=T) else expr_mat[, idx_post]
      
      if(max(v_pre, na.rm=T) > 50) v_pre <- log2(v_pre + 1)
      if(max(v_post, na.rm=T) > 50) v_post <- log2(v_post + 1)
      
      delta_list[[s]] <- v_post - v_pre
    }
  }
  if(length(delta_list) < 5) return(NULL)
  return(do.call(cbind, delta_list))
}

mat_mus_delta <- get_delta_matrix("Cleaned_Muscle_Proteomics.csv")
mat_adi_delta <- get_delta_matrix("Cleaned_Adipose_Proteomics.csv")

# Extract clinical traits matching exactly with delta matrix subjects
trait_cols <- c("Age", "Fat_percent", "HOMA_IR", "VO2max", "TRIGLY", "HDL", "LDL", "CRP")
datTraits <- clinical %>% 
  dplyr::select(Subject_ID, dplyr::any_of(trait_cols)) %>% 
  dplyr::distinct(Subject_ID, .keep_all = TRUE) %>% 
  tibble::column_to_rownames("Subject_ID")

# ------------------------------------------------------------------------------
# [2] Shared Functions: Plotting (WGCNA Bubble + CRP Barplot)
# ------------------------------------------------------------------------------
category_colors <- c(
  "Protein Homeostasis" = "#E64B35",               
  "Energy Metabolism" = "#4DBBD5",                 
  "Immune & Inflammatory Response" = "#00A087",    
  "Structural & Vesicular Transport" = "#3C5488",  
  "Other" = "#B0B0B0"                              
)

categorize_go <- function(go_term, tissue) {
  go_term <- str_trim(str_replace_all(go_term, "[\r\n]+", " "))
  if (tissue == "Adipose") {
    if (go_term %in% c("Regulation of protein stability", "Establishment of protein localization to vacuole", "Proteasome-mediated ubiquitin-dependent protein catabolic process")) return("Protein Homeostasis")
    if (go_term %in% c("Generation of precursor metabolites and energy", "Cellular respiration", "Nucleotide metabolic process")) return("Energy Metabolism")
    if (go_term %in% c("Complement activation, classical pathway", "Complement activation", "B cell receptor signaling pathway", "Neutrophil mediated immunity")) return("Immune & Inflammatory Response")
    if (go_term %in% c("Membrane fusion", "Golgi vesicle transport", "Cytoplasmic translation", "Actin cytoskeleton organization")) return("Structural & Vesicular Transport")
  } else if (tissue == "Muscle") {
    if (go_term %in% c("Proteasome-mediated ubiquitin-dependent protein catabolic process", "Regulation of protein stability")) return("Protein Homeostasis")
    if (go_term %in% c("Cellular respiration", "Mitochondrial gene expression", "Nucleotide metabolic process")) return("Energy Metabolism")
    if (go_term %in% c("Complement activation, classical pathway", "Neutrophil mediated immunity")) return("Immune & Inflammatory Response")
    if (go_term %in% c("Cytoplasmic translation", "Endoplasmic reticulum membrane organization", "Muscle contraction")) return("Structural & Vesicular Transport")
  }
  return("Other")
}

plot_crp_correlation <- function(plot_data, tissue_name, title, output_prefix) {
  df_filtered <- plot_data %>%
    dplyr::filter(Trait == "CRP") %>%
    dplyr::mutate(
      GO_Clean = str_trim(str_replace_all(Module_Label, "[\r\n]+", " ")),
      Category = sapply(GO_Clean, categorize_go, tissue = tissue_name),
      Abs_Cor = abs(Correlation)
    ) %>%
    dplyr::arrange(Abs_Cor) %>%
    dplyr::mutate(GO_Clean = factor(GO_Clean, levels = GO_Clean))
  
  if (nrow(df_filtered) == 0) return(NULL)
  
  p <- ggplot(df_filtered, aes(x = Correlation, y = GO_Clean, fill = Category)) +
    geom_col(width = 0.7, color = "black", size = 0.3) +
    scale_fill_manual(values = category_colors) +
    geom_text(
      aes(
        x = ifelse(Correlation > 0, Correlation + 0.05, Correlation - 0.05),
        label = sprintf("italic(P) == '%.3f'", Pvalue)
      ),
      parse = TRUE, size = 3.5, hjust = ifelse(df_filtered$Correlation > 0, 0, 1)
    ) +
    coord_cartesian(xlim = c(min(df_filtered$Correlation) * 1.3, max(df_filtered$Correlation) * 1.3)) +
    theme_classic(base_size = 12) +
    labs(title = title, x = "Correlation with CRP", y = NULL) +
    theme(
      axis.text.y = element_text(size = 11, color = "black"),
      axis.text.x = element_text(size = 11, color = "black"),
      axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      panel.grid.major.x = element_line(color = "grey90", linetype = "dashed")
    )
  
  ggsave(file.path(DIR_OUT, paste0(output_prefix, ".pdf")), plot = p, width = 10, height = 6, dpi = 300)
}

# ------------------------------------------------------------------------------
# [3] Unified WGCNA Engine
# ------------------------------------------------------------------------------
run_wgcna_pipeline <- function(delta_mat, tissue_name, pwr, min_mod_size) {
  if(is.null(delta_mat)) return()
  message(sprintf("\n=== STEP 2: Running %s WGCNA (Power=%d) ===", tissue_name, pwr))
  
  # Transpose and filter
  datExpr_all <- t(delta_mat)
  keepSamples <- (rowSums(!is.na(datExpr_all)) > (ncol(datExpr_all) * 0.5))
  datExpr <- datExpr_all[keepSamples, ]
  
  gsg = goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
  
  common_subs <- intersect(rownames(datExpr), rownames(datTraits))
  datExpr <- datExpr[common_subs, ]; traits <- datTraits[common_subs, ]
  
  # Network Construction
  net = blockwiseModules(datExpr, power = pwr, TOMType = "unsigned", minModuleSize = min_mod_size, 
                         reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, 
                         pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 0, corType = "pearson")
  
  moduleColors = labels2colors(net$colors)
  MEs = orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
  geneInfo = data.frame(Protein = colnames(datExpr), ModuleColor = moduleColors)
  
  # Module-Trait Correlation
  moduleTraitCor = WGCNA::cor(MEs, traits, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))
  
  df_cor_long <- as.data.frame(as.table(moduleTraitCor)) %>% dplyr::rename(Module=Var1, Trait=Var2, Correlation=Freq)
  df_p_long <- as.data.frame(as.table(moduleTraitPvalue)) %>% dplyr::rename(Module=Var1, Trait=Var2, Pvalue=Freq)
  df_long <- merge(df_cor_long, df_p_long, by=c("Module", "Trait"))
  
  # GO Enrichment
  message("  > Performing GO Enrichment...")
  module_stats <- df_long %>% dplyr::group_by(Module) %>% dplyr::summarise(Min_Pvalue = min(Pvalue, na.rm = TRUE))
  modules <- unique(geneInfo$ModuleColor); modules <- modules[modules != "grey"]
  module_label_df <- data.frame(Module = character(), Label = character(), stringsAsFactors = FALSE)
  
  for (mod in modules) {
    proteins <- geneInfo$Protein[geneInfo$ModuleColor == mod]
    entrez_ids <- tryCatch(bitr(proteins, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"), error = function(e) return(NULL))
    
    final_term <- "No Annotation"
    if (!is.null(entrez_ids) && nrow(entrez_ids) > 0) {
      ego <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
      if ((is.null(ego) || nrow(ego) == 0)) {
        ego <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "none", pvalueCutoff = 0.5, readable = TRUE)
      }
      if (!is.null(ego) && nrow(ego) > 0) {
        raw_term <- ego@result$Description[1]
        raw_term_cap <- paste0(toupper(substr(raw_term, 1, 1)), substr(raw_term, 2, nchar(raw_term)))
        final_term <- str_wrap(raw_term_cap, width = 35)
      }
    }
    module_label_df <- rbind(module_label_df, data.frame(Module = paste0("ME", mod), Label = final_term))
  }
  
  # Deduplication
  merged_info <- merge(module_label_df, module_stats, by = "Module") %>% dplyr::filter(Label != "No Annotation")
  
  if(tissue_name == "Adipose") { merged_info <- merged_info %>% dplyr::filter(Min_Pvalue < 0.05) }
  
  best_modules_df <- merged_info %>% dplyr::group_by(Label) %>% dplyr::arrange(Min_Pvalue) %>% dplyr::slice(1) %>% dplyr::ungroup()
  if (nrow(best_modules_df) > 10) best_modules_df <- best_modules_df %>% dplyr::arrange(Min_Pvalue) %>% head(10)
  
  if (nrow(best_modules_df) == 0) { message("  ! No valid modules left after filtering."); return(NULL) }
  
  # Prepare Plot Data
  plot_data <- df_long %>% dplyr::filter(Module %in% best_modules_df$Module)
  plot_data$Module_Label <- best_modules_df$Label[match(plot_data$Module, best_modules_df$Module)]
  
  SUPP_TABLE_LIST[[tissue_name]] <<- plot_data %>% dplyr::rename(Module_Color = Module, GO_Annotation = Module_Label)
  
  # Axis Ordering
  trait_order <- trait_cols
  current_traits <- unique(plot_data$Trait)
  plot_data$Trait <- factor(plot_data$Trait, levels = c(intersect(trait_order, current_traits), setdiff(current_traits, trait_order)))
  
  mat_cor <- dcast(plot_data, Module_Label ~ Trait, value.var = "Correlation")
  rownames(mat_cor) <- mat_cor$Module_Label; mat_cor <- mat_cor[, -1]; mat_cor[is.na(mat_cor)] <- 0
  
  final_y_levels <- if(nrow(mat_cor) > 2) hclust(dist(mat_cor, method = "euclidean"), method = "ward.D2")$labels[hclust(dist(mat_cor, method = "euclidean"), method = "ward.D2")$order] else rownames(mat_cor)
  
  plot_data$Module_Label <- factor(plot_data$Module_Label, levels = final_y_levels)
  plot_data$LogP <- pmin(-log10(plot_data$Pvalue), 5)
  
  # Plot WGCNA Bubble
  p_bubble <- ggplot(plot_data, aes(x = Trait, y = Module_Label)) +
    geom_point(aes(size = LogP, fill = Correlation), shape = 21, color = "black", stroke = 0.5) +
    scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C", midpoint = 0, limit = c(-0.6, 0.6), name = "Correlation", oob = scales::squish) +
    scale_size_continuous(range = c(3, 9), name = "-log10(P-value)") + 
    coord_fixed(ratio = 1) + 
    theme_minimal(base_size = 14, base_family = "sans") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color="black", face="bold"),
      axis.text.y = element_text(color="black", size = 11, lineheight = 0.8),
      panel.grid.major = element_line(color = "grey95", linewidth = 0.2),
      panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
      legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(x = "", y = "", title = paste0(tissue_name, " Module-Trait Associations"))
  
  calc_width <- 3 + (length(unique(plot_data$Trait)) * 0.6); calc_height <- 2 + (length(unique(plot_data$Module_Label)) * 0.6)
  ggsave(file.path(DIR_OUT, paste0("Plot_Bubble_WGCNA_", tissue_name, ".pdf")), p_bubble, width = calc_width, height = calc_height)
  
  # Plot CRP Barplot
  plot_crp_correlation(plot_data, tissue_name, paste0(tissue_name, " Proteomic Networks & Systemic Inflammation (CRP)"), paste0("Plot_Bar_CRP_", tissue_name))
  
  message("  > Finished pipeline for ", tissue_name)
}

# ------------------------------------------------------------------------------
# [4] Execute All Pipelines
# ------------------------------------------------------------------------------
run_wgcna_pipeline(mat_adi_delta, "Adipose", pwr = 5, min_mod_size = 25)
run_wgcna_pipeline(mat_mus_delta, "Muscle", pwr = 3, min_mod_size = 30)

# Export Supplementary Table
write.xlsx(SUPP_TABLE_LIST, file.path(DIR_OUT, "Data_WGCNA_Module_Traits.xlsx"), overwrite = TRUE)

message("\n=======================================================")
message(">>> WGCNA Integrated Pipeline Completed Successfully! <<<")
message(">>> All outputs securely stored in: ", DIR_OUT)
message("=======================================================")
