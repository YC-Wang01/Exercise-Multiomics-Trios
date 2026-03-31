# ==============================================================================
# Project: FinlandSports V2.0
# Script:  05_Fig2CD_S2_GSEA_Network.R
# Panels:  Fig 2C (GSEA Matrix), Fig 2D (Network), Fig S2 (Crosstalk)
# Style:   Original Theme with Black Axis/Borders + No Grids.
# ==============================================================================

# [0] Set working directory and initialize environment
# setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  limma, readr, readxl, dplyr, stringr, 
  clusterProfiler, org.Hs.eg.db, ggplot2, 
  ggrepel, tibble, grid, openxlsx, enrichplot, cowplot,
  igraph, ggraph, Hmisc, tidyr, ggforce
)

options(timeout = 600)

DIR_DATA    <- "01_Clean_Data"
DIR_OUT_F2  <- "03_Results/Fig_2"
DIR_OUT_S2  <- "03_Results/Fig_S2"
DIR_OUT_CSV <- "03_Results/Fig_2/Serum-meta" 
DIR_TABLES  <- file.path(DIR_OUT_F2, "Tables")

invisible(lapply(c(DIR_OUT_F2, DIR_OUT_S2, DIR_OUT_CSV, DIR_TABLES), dir.create, recursive = TRUE, showWarnings = FALSE))

my_clean_theme <- theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_blank(), 
    panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA), 
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black", face = "bold"),
    axis.title = element_text(color = "black", face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", face = "bold", size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# ------------------------------------------------------------------------------
# [1] Clinical Data Integration (Restore variable names for consistency)
# ------------------------------------------------------------------------------
message(">>> Loading...")
clin_master <- read_csv(file.path(DIR_DATA, "Clinical_Master_Strict.csv"), show_col_types = FALSE) %>%
  dplyr::filter(!is.na(`Fat(%)`)) %>%
  dplyr::mutate(
    Role_Orig = dplyr::case_when(Membercode == 1 ~ "Daughter", Membercode == 2 ~ "Mother", Membercode == 3 ~ "Father", TRUE ~ "Unknown"),
    Subject_ID = tolower(paste0(FamilyID, "_", Role_Orig)),
    Fat_percent = `Fat(%)`, 
    Group = ifelse(Fat_percent >= 30, "Obese", "Lean") 
  )

# ------------------------------------------------------------------------------
# [2] Module A: Limma-GSEA Analysis Engine
# ------------------------------------------------------------------------------
gsea_full_list <- list() 

generate_full_gsea <- function(filename, prefix) {
  full_path <- file.path(DIR_DATA, filename)
  out_csv   <- file.path(DIR_OUT_CSV, paste0(prefix, "_GSEA_Full.csv"))
  if(file.exists(out_csv)) {
    message("  [Check] GSEA results already exist, skipping: ", basename(out_csv))
    return(TRUE)
  }
  if(!file.exists(full_path)) return(FALSE)
  
  message(paste("\n[Analysis] Running GSEA for:", prefix, "..."))
  data_raw <- read_csv(full_path, show_col_types = F)
  colnames(data_raw)[1] <- "FeatureID"
  expr_mat <- as.matrix(data_raw %>% dplyr::distinct(FeatureID, .keep_all = T) %>% tibble::column_to_rownames("FeatureID"))
  
  matched_meta <- data.frame(Omics_ID = colnames(expr_mat), stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      Clean_Name = tolower(gsub("_twin\\.[0-9]+", "", Omics_ID)),
      TimeType = dplyr::case_when(grepl("pre", Clean_Name) ~ "pre", grepl("fast", Clean_Name) ~ "fast", TRUE ~ "Other"),
      Subject_ID = gsub("_fast|_pre|_post1h|_post3h", "", Clean_Name)
    ) %>%
    dplyr::inner_join(clin_master, by = "Subject_ID") %>%
    dplyr::filter(TimeType %in% c("pre", "fast")) %>%
    dplyr::group_by(Subject_ID) %>% dplyr::slice(1) %>% dplyr::ungroup()
  
  if(nrow(matched_meta) < 6) return(FALSE)
  
  mat_sub <- expr_mat[, matched_meta$Omics_ID]
  mat_sub[is.na(mat_sub)] <- min(mat_sub, na.rm = T)/2
  if(max(mat_sub, na.rm = T) > 100) mat_sub <- log2(mat_sub + 1)
  
  found_genes <- intersect(rownames(mat_sub), c("HBB", "HBA1", "HBA2"))
  blood_score_vec <- if(length(found_genes) == 0) rep(0, ncol(mat_sub)) else as.numeric(scale(colSums(mat_sub[found_genes, , drop=F], na.rm=T)))
  
  groups <- factor(matched_meta$Group, levels = c("Lean", "Obese"))
  design <- if(sd(blood_score_vec, na.rm=T) == 0) model.matrix(~0 + groups) else model.matrix(~0 + groups + blood_score_vec)
  colnames(design)[1:2] <- levels(groups)
  
  fit <- eBayes(contrasts.fit(lmFit(mat_sub, design), makeContrasts(Obese - Lean, levels = design)))
  res_limma <- topTable(fit, number = Inf) %>% tibble::rownames_to_column("Gene")
  
  gene_map <- tryCatch(bitr(res_limma$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"), error=function(e) NULL)
  if(is.null(gene_map)) return(FALSE)
  
  gsea_input <- res_limma %>% dplyr::inner_join(gene_map, by = c("Gene" = "SYMBOL")) %>% dplyr::arrange(desc(t)) %>% dplyr::distinct(ENTREZID, .keep_all = T)
  gene_list <- gsea_input$t; names(gene_list) <- gsea_input$ENTREZID
  gsea_res <- gseKEGG(geneList = gene_list, organism = 'hsa', pvalueCutoff = 1, minGSSize = 3, maxGSSize = 2000, verbose = F, seed = 123)
  if(is.null(gsea_res)) return(FALSE)
  
  write_csv(as.data.frame(setReadable(gsea_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")), out_csv)
  return(TRUE)
}

generate_full_gsea("Cleaned_Adipose_Microarray.csv", "Fig2A_Adipose_Transcriptome")
generate_full_gsea("Cleaned_Adipose_Proteomics.csv", "Fig2B_Adipose_Proteomics")
generate_full_gsea("Cleaned_Muscle_Microarray.csv", "Fig2C_Muscle_Transcriptome")
generate_full_gsea("Cleaned_Muscle_Proteomics.csv", "Fig2D_Muscle_Proteomics")

# ------------------------------------------------------------------------------
# [3] Module B: Crosstalk Quadrant Plots (Fig_S2)
# ------------------------------------------------------------------------------
plot_compact_crosstalk <- function(prefix_adi, prefix_mus, title_text) {
  file_adi <- file.path(DIR_OUT_CSV, paste0(prefix_adi, "_GSEA_Full.csv"))
  file_mus <- file.path(DIR_OUT_CSV, paste0(prefix_mus, "_GSEA_Full.csv"))
  if(!file.exists(file_adi) | !file.exists(file_mus)) return(NULL)
  
  df_merge <- dplyr::inner_join(
    read_csv(file_adi, show_col_types = F) %>% dplyr::select(ID, Description, NES, p.adjust, setSize) %>% dplyr::rename(NES_Adi = NES, FDR_Adi = p.adjust),
    read_csv(file_mus, show_col_types = F) %>% dplyr::select(ID, NES, p.adjust) %>% dplyr::rename(NES_Mus = NES, FDR_Mus = p.adjust),
    by = "ID"
  ) %>% dplyr::mutate(
    Is_Significant = (FDR_Adi < 0.05 | FDR_Mus < 0.05),
    Quadrant = dplyr::case_when(NES_Mus > 0 & NES_Adi > 0 ~ "Q1: Co-Activation", NES_Mus < 0 & NES_Adi > 0 ~ "Q2: Adi Up / Mus Down", 
                                NES_Mus < 0 & NES_Adi < 0 ~ "Q3: Co-Suppression", NES_Mus > 0 & NES_Adi < 0 ~ "Q4: Mus Up / Adi Down", TRUE ~ "Center"),
    Vector_Strength = sqrt(NES_Mus^2 + NES_Adi^2),
    Is_Target = grepl("Insulin|Lipid|TCA|Metabolic|Inflammation|Stress|Thermogenesis|Oxidative", Description, ignore.case = T)
  )
  
  labels_red <- df_merge %>% dplyr::filter(Is_Significant & Is_Target) %>% dplyr::group_by(Quadrant) %>% dplyr::slice_max(order_by = Vector_Strength, n = 6) %>% dplyr::pull(Description)
  df_merge <- df_merge %>% dplyr::mutate(Label_Text = ifelse(Description %in% labels_red, str_wrap(Description, width = 15), NA))
  limit_val <- max(abs(c(df_merge$NES_Mus, df_merge$NES_Adi)), na.rm = T) * 1.05  
  
  p <- ggplot(df_merge, aes(x = NES_Mus, y = NES_Adi)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_point(data = subset(df_merge, !Is_Significant), aes(size = setSize), shape = 16, color = "grey92", alpha = 0.5) +
    geom_point(data = subset(df_merge, Is_Significant), aes(color = Quadrant, size = setSize), shape = 16, alpha = 0.85) +
    geom_text_repel(aes(label = Label_Text), color = "#D00000", size = 3.5, fontface = "bold", box.padding = 0.4, max.overlaps = 50) +
    scale_color_manual(values = c("Q1: Co-Activation"="#d73027", "Q2: Adi Up / Mus Down"="#fdae61", "Q3: Co-Suppression"="#4575b4", "Q4: Mus Up / Adi Down"="#abd9e9", "Center"="grey90")) +
    my_clean_theme + coord_fixed(ratio = 1, xlim = c(-limit_val, limit_val), ylim = c(-limit_val, limit_val)) +
    labs(title = paste0(title_text, " Crosstalk"), x = "Muscle NES", y = "Adipose NES") +
    theme(legend.position = "none", aspect.ratio = 1) 
  
  ggsave(file.path(DIR_OUT_S2, paste0("FigS2_Crosstalk_", title_text, "_Final.pdf")), p, width = 8, height = 8)
  return(df_merge)
}

crosstalk_list <- list()
crosstalk_list[["Transcriptome"]] <- plot_compact_crosstalk("Fig2A_Adipose_Transcriptome", "Fig2C_Muscle_Transcriptome", "Transcriptome")
crosstalk_list[["Proteomics"]]    <- plot_compact_crosstalk("Fig2B_Adipose_Proteomics", "Fig2D_Muscle_Proteomics", "Proteomics")

# ------------------------------------------------------------------------------
# [4] Module C: Ultra-Compact GSEA Bubble Matrix (Fig. 2C)
# ------------------------------------------------------------------------------
message(">>> Generating GSEA comparison matrix...")
files_to_load <- list("Adipose RNA" = "Fig2A_Adipose_Transcriptome_GSEA_Full.csv", "Adipose Protein" = "Fig2B_Adipose_Proteomics_GSEA_Full.csv", "Muscle RNA" = "Fig2C_Muscle_Transcriptome_GSEA_Full.csv", "Muscle Protein" = "Fig2D_Muscle_Proteomics_GSEA_Full.csv")
df_all <- dplyr::bind_rows(lapply(names(files_to_load), function(n) { p <- file.path(DIR_OUT_CSV, files_to_load[[n]]); if(file.exists(p)) read_csv(p, show_col_types=F) %>% dplyr::mutate(Dataset=n) else NULL }))

if(nrow(df_all) > 0) {
  curated_dict <- data.frame(
    Description = c("Phagosome", "Hematopoietic cell lineage", "Leukocyte transendothelial migration", "Natural killer cell mediated cytotoxicity", "Platelet activation", "Type I diabetes mellitus", "Oxidative phosphorylation", "Citrate cycle (TCA cycle)", "Thermogenesis", "Carbon metabolism", "Valine, leucine and isoleucine degradation", "2-Oxocarboxylic acid metabolism", "Glyoxylate and dicarboxylate metabolism"),
    Pathway_Group = factor(c(rep("Immune", 6), rep("Metabolism", 7)), levels = c("Immune", "Metabolism"))
  )
  dataset_labels <- c("Adipose RNA" = "Adipose\n(RNA)", "Adipose Protein" = "Adipose\n(Protein)", "Muscle RNA" = "Muscle\n(RNA)", "Muscle Protein" = "Muscle\n(Protein)")
  
  plot_data <- df_all %>% dplyr::inner_join(curated_dict, by = "Description") %>%
    dplyr::mutate(minus_log10_padj = -log10(p.adjust + 1e-10), Is_Significant = p.adjust < 0.05, Dataset_Label = factor(dataset_labels[Dataset], levels = dataset_labels)) %>%
    dplyr::group_by(Description) %>% dplyr::mutate(Importance_Score = sum(abs(NES) * minus_log10_padj, na.rm = T)) %>% dplyr::ungroup() %>%
    dplyr::mutate(Description = reorder(Description, Importance_Score))
  
  p_matrix <- ggplot(plot_data, aes(x = Dataset_Label, y = Description)) +
    geom_point(aes(size = minus_log10_padj, fill = NES, color = Is_Significant, alpha = Is_Significant), shape = 21, stroke = 0.5) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-3, 3)) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "transparent"), guide = "none") +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4), guide = "none") +
    scale_size_continuous(range = c(2, 7)) +
    facet_grid(Pathway_Group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_discrete(position = "right") + my_clean_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", color = "black", size = 11), 
      axis.text.y = element_text(color = "black", size = 11),
      axis.title = element_blank(),
      strip.placement = "outside", 
      strip.text.y.left = element_text(face = "bold", angle = 90, color = "black", size = 13), 
      plot.margin = margin(t = 10, r = 10, b = 20, l = 15, unit = "pt")
    )
  ggsave(file.path(DIR_OUT_F2, "Fig2C_GSEA_Comparison_Matrix_Ultra_Compact.pdf"), p_matrix, width = 6.8, height = 4.2)
}

# ------------------------------------------------------------------------------
# [5] Module D: Global Spearman Network (Fig. 2D)
# ------------------------------------------------------------------------------
message(">>> Building multi-omics network...")
get_expression_safe <- function(filename, hub_genes=NULL) {
  raw_mat <- read_csv(file.path(DIR_DATA, filename), show_col_types = F); colnames(raw_mat)[1] <- "Gene"
  if(!is.null(hub_genes)) raw_mat <- raw_mat %>% dplyr::filter(Gene %in% hub_genes)
  
  expr_data <- raw_mat %>% dplyr::distinct(Gene, .keep_all = T) %>% tibble::column_to_rownames("Gene") %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("Omics_ID")
  
  dict <- data.frame(Omics_ID = expr_data$Omics_ID, stringsAsFactors = FALSE) %>%
    dplyr::mutate(Clean_Name = tolower(gsub("_twin\\.[0-9]+", "", Omics_ID)), Subject_ID = gsub("_fast|_pre|_post1h|_post3h", "", Clean_Name)) %>%
    dplyr::inner_join(clin_master %>% dplyr::select(Subject_ID), by = "Subject_ID")
  
  expr_data %>% dplyr::inner_join(dict, by = "Omics_ID") %>% dplyr::select(-Omics_ID, -Clean_Name) %>% 
    dplyr::mutate(dplyr::across(-Subject_ID, as.numeric)) %>% dplyr::group_by(Subject_ID) %>% dplyr::summarise(dplyr::across(dplyr::everything(), \(x) mean(x, na.rm = T)), .groups = 'drop')
}

adi_hubs <- c("ITGAM", "CD14", "CYBB", "ITGB2", "FCGR3A", "HLA-DRA", "CD36", "PLCG2", "FCER1G", "ITGA5")
mus_hubs <- c("ATP5F1A", "NDUFA9", "SDHA", "MDH1", "ACO2", "IDH2", "CS", "DLAT", "ATP5F1B", "ATP5F1C")

tissue_df <- get_expression_safe("Cleaned_Adipose_Microarray.csv", adi_hubs) %>% dplyr::inner_join(get_expression_safe("Cleaned_Muscle_Proteomics.csv", mus_hubs), by="Subject_ID")
s_prot <- get_expression_safe(basename(list.files(DIR_DATA, pattern="Serum.*prot", ignore.case=T, full.names=T)[1]))
s_metab <- get_expression_safe(basename(list.files(DIR_DATA, pattern="Serum.*metab", ignore.case=T, full.names=T)[1]))

master_df <- clin_master %>% dplyr::select(Subject_ID, Fat_percent, HOMA_IR, CRP, Leptin, Adiponectin, TRIGLY, NEFA, VO2max) %>% dplyr::distinct() %>% dplyr::inner_join(tissue_df, by="Subject_ID")
if(nrow(s_prot) > 0) master_df <- master_df %>% dplyr::inner_join(s_prot, by="Subject_ID")
if(nrow(s_metab) > 0) master_df <- master_df %>% dplyr::inner_join(s_metab, by="Subject_ID")

rcorr_res <- rcorr(as.matrix(master_df %>% tibble::column_to_rownames("Subject_ID") %>% dplyr::mutate_all(as.numeric)), type = "spearman")
edges <- rcorr_res$r %>% as.data.frame() %>% tibble::rownames_to_column("Node1") %>% tidyr::pivot_longer(cols = -Node1, names_to = "Node2", values_to = "R") %>% dplyr::filter(Node1 < Node2)
p_edges <- rcorr_res$P %>% as.data.frame() %>% tibble::rownames_to_column("Node1") %>% tidyr::pivot_longer(cols = -Node1, names_to = "Node2", values_to = "Pval") %>% dplyr::filter(Node1 < Node2)

node_cats <- data.frame(Name = colnames(master_df)[-1]) %>% dplyr::mutate(Category = dplyr::case_when(Name %in% adi_hubs ~ "Adipose Inflammation", Name %in% mus_hubs ~ "Muscle Metabolism", Name %in% c("Fat_percent", "HOMA_IR", "CRP", "Leptin", "Adiponectin", "TRIGLY", "NEFA", "VO2max") ~ "Clinical Phenotype", TRUE ~ "Serum Factors"))

# ================= Core Network Construction =================
net_data <- edges %>% dplyr::inner_join(p_edges, by = c("Node1", "Node2")) %>% 
  dplyr::left_join(node_cats, by = c("Node1" = "Name")) %>% dplyr::rename(Cat1 = Category) %>% 
  dplyr::left_join(node_cats, by = c("Node2" = "Name")) %>% dplyr::rename(Cat2 = Category) %>% 
  dplyr::filter(Pval < 0.05, abs(R) > 0.35, Cat1 != Cat2) %>% 
  dplyr::mutate(Cor_Direction = ifelse(R > 0, "Positive", "Negative")) 

g <- graph_from_data_frame(d = net_data, vertices = node_cats %>% dplyr::filter(Name %in% unique(c(net_data$Node1, net_data$Node2))), directed = F)
g <- delete_vertices(g, V(g)[degree(g) < 2])
comps <- components(g)
g <- induced_subgraph(g, V(g)[comps$membership == which.max(comps$csize)])
V(g)$Degree <- degree(g)
E(g)$Weight <- abs(E(g)$R)
V(g)$Lipid_Cluster <- dplyr::case_when(grepl("VLDL", V(g)$name) ~ "VLDL Subfractions", grepl("HDL", V(g)$name) ~ "HDL Subfractions", grepl("LDL", V(g)$name) & !grepl("VLDL", V(g)$name) ~ "LDL Subfractions", TRUE ~ NA_character_ )

volcano_file <- file.path(DIR_OUT_CSV, "Crosscheck_Network_Nodes_in_Volcano.csv")
sig_features <- if(file.exists(volcano_file)) { read_csv(volcano_file, show_col_types = FALSE) %>% dplyr::filter(grepl("✅", Status_in_Volcano)) %>% dplyr::pull(Feature) } else { c() }

temp_label <- ifelse(V(g)$Degree >= 5 | V(g)$name %in% c("Fat_percent", "HOMA_IR", "CRP", "Leptin", "VO2max") | V(g)$name %in% sig_features, V(g)$name, "")
V(g)$Show_Label <- ifelse(!is.na(V(g)$Lipid_Cluster), "", ifelse(temp_label == "Fat_percent", "Fat%", temp_label))

set.seed(121) 
p_net <- ggraph(g, layout = 'stress') + 
  geom_edge_link(aes(edge_width = Weight, alpha = Weight, color = Cor_Direction)) +
  scale_edge_color_manual(values = c("Positive" = "#E64B35", "Negative" = "#4DBBD5")) +
  scale_edge_alpha_continuous(range = c(0.1, 0.45), guide = "none") + scale_edge_width_continuous(range = c(0.3, 1.0), guide = "none") +
  geom_mark_hull(aes(x = x, y = y, filter = !is.na(Lipid_Cluster), fill = Lipid_Cluster, label = Lipid_Cluster), concavity = 1.5, expand = unit(3, "mm"), radius = unit(3, "mm"), alpha = 0.15, color = "gray50", linetype = "dashed", label.fontsize = 11, label.fontface = "bold", label.fill = "white", label.buffer = unit(5, "mm"), con.colour = "gray50", con.cap = unit(0, "mm")) +
  geom_node_point(aes(size = Degree, fill = Category), shape = 21, color = "white", stroke = 1.0) +
  scale_size_continuous(range = c(3, 8), guide = "none") + 
  geom_node_text(aes(label = Show_Label), size = 4.0, fontface = "bold", color = "black", repel = TRUE, force = 1.0, point.padding = 0.1, box.padding = 0.3, min.segment.length = unit(0.5, "lines"), segment.color = "grey60", segment.size = 0.3, bg.color = "white", bg.r = 0.15, max.overlaps = Inf) +
  scale_fill_manual(values = c("Adipose Inflammation" = "#E64B35", "Muscle Metabolism" = "#4DBBD5", "Clinical Phenotype" = "#00A087", "Serum Factors" = "#F39C12", "VLDL Subfractions" = "#FBC497", "HDL Subfractions"  = "#88C5CA", "LDL Subfractions"  = "#A1D99B")) + 
  theme_void(base_size = 14) + theme(legend.position = "right", plot.margin = margin(15, 15, 15, 15))

ggsave(file.path(DIR_OUT_F2, "Fig2D_Ultimate_Crosstalk_Network.pdf"), p_net, width = 11, height = 8.5)

# ------------------------------------------------------------------------------
# [6] Export Supplementary Data (Without Table_S Prefix)
# ------------------------------------------------------------------------------
message(">>> Exporting data tables...")

wb_gsea <- createWorkbook()
prefixes <- c("Fig2A_Adipose_Transcriptome", "Fig2B_Adipose_Proteomics", "Fig2C_Muscle_Transcriptome", "Fig2D_Muscle_Proteomics")
for(p in prefixes) {
  csv_file <- file.path(DIR_OUT_CSV, paste0(p, "_GSEA_Full.csv"))
  if(file.exists(csv_file)) {
    sheet_name <- substr(p, 7, 30)
    addWorksheet(wb_gsea, sheet_name); writeData(wb_gsea, sheet_name, read_csv(csv_file, show_col_types=F))
  }
}
saveWorkbook(wb_gsea, file.path(DIR_TABLES, "GSEA_Full_Results.xlsx"), overwrite = TRUE)

wb_net <- createWorkbook()
if(exists("plot_data")) { addWorksheet(wb_net, "GSEA_Matrix_Data"); writeData(wb_net, "GSEA_Matrix_Data", plot_data) }
if(exists("net_data")) { addWorksheet(wb_net, "Network_Edges"); writeData(wb_net, "Network_Edges", net_data) }
for(n in names(crosstalk_list)) {
  if(!is.null(crosstalk_list[[n]])) { addWorksheet(wb_net, paste0("Crosstalk_", n)); writeData(wb_net, paste0("Crosstalk_", n), crosstalk_list[[n]]) }
}
saveWorkbook(wb_net, file.path(DIR_TABLES, "Network_and_Crosstalk_Details.xlsx"), overwrite = TRUE)

message("\n>>> Task Completed Successfully!")
