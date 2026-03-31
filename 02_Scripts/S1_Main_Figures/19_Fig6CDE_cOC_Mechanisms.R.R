# ==============================================================================
# Project: FinlandSports V2.0 - Fig 6 C, D, E: cOC Mechanisms & Networks
# Description: 
#   1. Fig 6C Heatmap: Strict 18 features (12 Adipose, 6 Muscle) vs cOC gradient.
#   2. Fig 6D/E Networks: cOC vs tOC pure multi-omics network (Lean & Obese).
#   3. Aesthetics: Anti-crowding, centered labels, bright NPG palette.
#   4. V2.0 Standard: 100% English, Pure CSV loading, strict namespacing.
# ==============================================================================

# --- 0. ENVIRONMENT SETUP ---
# setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readr, dplyr, tidyr, stringr, openxlsx, tibble,
  ComplexHeatmap, circlize, grid,
  igraph, ggraph, ggforce, concaveman, Hmisc, cowplot, ggplot2
)

# --- 1. GLOBAL CONFIGURATION ---
DIR_IN  <- "01_Clean_Data"
DIR_OUT <- "03_Results/Fig_6/Fig6CDE_cOC_Mechanisms"
if(!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)

message(">>> Directories Ready. Output set to: ", DIR_OUT)

# ==============================================================================
# --- 2. V2.0 SHARED DATA ENGINE (Clinical Baseline + Omics Delta) ---
# ==============================================================================
message("\n=== STEP 1: Loading Clinical Metadata & Calculating Omics Deltas ===")

# 1. Clinical Data (FIXED: Using exact column names from Clinical_Master_Strict.csv)
df_clin <- read_csv(file.path(DIR_IN, "Clinical_Master_Strict.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(
    Role_Orig = dplyr::case_when(Membercode == 1 ~ "Daughter", Membercode == 2 ~ "Mother", Membercode == 3 ~ "Father", TRUE ~ "Unknown"),
    Sample_Key = tolower(paste0(FamilyID, "_", Role_Orig)),
    Fat_percent = as.numeric(`Fat(%)`),
    Fat_Group = ifelse(!is.na(Fat_percent) & Fat_percent >= 30, "Obese", "Lean"),
    S_cOC = as.numeric(S_cOC),             # Fixed here
    Clin_S_TotalOC = as.numeric(S_TotalOC) # Fixed here
  ) %>%
  dplyr::filter(!is.na(S_cOC) & !is.na(Clin_S_TotalOC) & !is.na(Fat_Group)) %>%
  dplyr::distinct(Sample_Key, .keep_all = TRUE)

# 2. Omics Delta Engine
get_omics_delta_v2 <- function(filename, suffix) {
  fpath <- file.path(DIR_IN, filename)
  if(!file.exists(fpath)) return(NULL)
  
  raw <- read_csv(fpath, show_col_types = FALSE)
  feat_ids <- raw[[1]]; expr_raw <- as.matrix(raw[, -1])
  rownames(expr_raw) <- feat_ids
  
  clean_cols <- tolower(str_remove(colnames(expr_raw), "_twin\\.\\d+$")) %>% str_trim()
  col_meta <- data.frame(Matrix_Col = colnames(expr_raw), Clean_Col = clean_cols, stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      TimePoint = dplyr::case_when(grepl("pre", Clean_Col) ~ "Pre", grepl("fast", Clean_Col) ~ "Pre", grepl("post3h", Clean_Col) ~ "Post3H", TRUE ~ "Other"),
      Subject_ID = gsub("_fast|_pre|_post1h|_post3h", "", Clean_Col)
    )
  
  df_wide <- col_meta %>% dplyr::filter(TimePoint %in% c("Pre", "Post3H")) %>%
    dplyr::group_by(Subject_ID, TimePoint) %>% dplyr::slice(1) %>% dplyr::ungroup()
  
  res_delta <- data.frame(Subject_ID = unique(df_wide$Subject_ID))
  
  for(m in rownames(expr_raw)) {
    v_pre <- expr_raw[m, df_wide$Matrix_Col[df_wide$TimePoint == "Pre"]]
    names(v_pre) <- df_wide$Subject_ID[df_wide$TimePoint == "Pre"]
    v_post <- expr_raw[m, df_wide$Matrix_Col[df_wide$TimePoint == "Post3H"]]
    names(v_post) <- df_wide$Subject_ID[df_wide$TimePoint == "Post3H"]
    
    delta_vec <- c()
    for(s in res_delta$Subject_ID) {
      if(!is.na(v_pre[s]) && !is.na(v_post[s])) {
        if(v_pre[s] > 50) delta_vec <- c(delta_vec, log2(v_post[s]+1) - log2(v_pre[s]+1))
        else delta_vec <- c(delta_vec, v_post[s] - v_pre[s])
      } else { delta_vec <- c(delta_vec, NA) }
    }
    res_delta[[paste0(m, suffix)]] <- delta_vec
  }
  return(res_delta)
}

df_mus_delta <- get_omics_delta_v2("Cleaned_Muscle_Proteomics.csv", "_Mus")
df_adi_delta <- get_omics_delta_v2("Cleaned_Adipose_Proteomics.csv", "_Adi")

# Merge into Global Master
df_master_global <- df_clin %>% dplyr::select(Sample_Key, Fat_Group, S_cOC, Clin_S_TotalOC) %>%
  dplyr::inner_join(df_mus_delta, by=c("Sample_Key"="Subject_ID")) %>%
  dplyr::inner_join(df_adi_delta, by=c("Sample_Key"="Subject_ID"))

# ==============================================================================
# === MODULE 1: FIG 6C cOC DELTA WATERFALL HEATMAP ===
# ==============================================================================
message("\n=== STEP 2: Rendering Fig 6C (Waterfall Heatmap) ===")

core_mus <- c("PABPN1", "FAM3C")
top_mus  <- c("MYLK", "TNNT1", "ATP2A2", "MYOZ1")
core_adi <- c("CRACDL", "TGFB1")
top_adi  <- c("ADIPOQ", "FABP4", "LEP", "LPL", "PLIN1", "FASN", "SCD", "ACACA", "ACLY", "SREBF1") 

mus_targets <- paste0(c(core_mus, top_mus), "_Mus")
adi_targets <- paste0(c(core_adi, top_adi), "_Adi")
all_targets <- c(mus_targets, adi_targets)

# 1. Structure Matrix
anno_df_hm <- df_master_global %>% dplyr::select(Sample_Key, Fat_Group, S_cOC) %>% na.omit() %>%
  dplyr::arrange(factor(Fat_Group, levels = c("Lean", "Obese")), S_cOC)

valid_targets <- intersect(all_targets, colnames(df_master_global))
expr_mat <- t(df_master_global[, valid_targets])
colnames(expr_mat) <- df_master_global$Sample_Key
expr_mat <- expr_mat[, anno_df_hm$Sample_Key]

# 2. Z-score & Correlation
z_mat <- t(apply(expr_mat, 1, function(x) {
  if (sd(x, na.rm=TRUE) > 0) return((x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)) else return(rep(0, length(x)))
}))

cor_res <- apply(z_mat, 1, function(x) {
  if(sum(!is.na(x)) > 5 && sd(x, na.rm=TRUE) > 0) {
    ct <- suppressWarnings(cor.test(x, anno_df_hm$S_cOC, method="spearman", exact=FALSE))
    c(R_value = unname(ct$estimate), P_value = unname(ct$p.value))
  } else { c(R_value = 0, P_value = 1) }
})
cor_df <- as.data.frame(t(cor_res))

row_anno_df <- data.frame(Feature = rownames(z_mat), Tissue = ifelse(grepl("_Mus", rownames(z_mat)), "Muscle", "Adipose"), R_value = cor_df$R_value, P_value = cor_df$P_value, stringsAsFactors = FALSE) %>% 
  dplyr::mutate(Significance = dplyr::case_when(P_value < 0.001 ~ "***", P_value < 0.01 ~ "**", P_value < 0.05 ~ "*", TRUE ~ "ns")) %>%
  dplyr::arrange(factor(Tissue, levels = c("Adipose", "Muscle")), desc(R_value))

z_mat <- z_mat[row_anno_df$Feature, ]
rownames(z_mat) <- gsub("_Mus|_Adi", "", rownames(z_mat))

# 3. Plotting
col_fun_heatmap <- colorRamp2(c(-1.5, 0, 1.5), c("#3C5488", "white", "#DC0000"))
col_coc <- colorRamp2(c(min(anno_df_hm$S_cOC), median(anno_df_hm$S_cOC), max(anno_df_hm$S_cOC)), c("white", "#FFC107", "#DC0000"))

top_anno <- HeatmapAnnotation(
  Fat_Group = anno_df_hm$Fat_Group, cOC_Level = anno_df_hm$S_cOC,
  col = list(Fat_Group = c("Obese" = "#F79647", "Lean" = "#3DA6AE"), cOC_Level = col_coc),
  annotation_label = c(Fat_Group = "Fat(%)", cOC_Level = "cOC"),
  annotation_name_side = "right", show_annotation_name = TRUE,
  annotation_legend_param = list(cOC_Level = list(title = "cOC Baseline", at = round(seq(min(anno_df_hm$S_cOC), max(anno_df_hm$S_cOC), length.out=3), 2)))
)

left_anno <- rowAnnotation(
  Tissue = row_anno_df$Tissue,
  Cor_R = anno_barplot(row_anno_df$R_value, bar_width = 1, gp = gpar(fill = ifelse(row_anno_df$R_value > 0, "#E64B35", "#3C5488"), col=NA), axis_param = list(side = "top")),
  col = list(Tissue = c("Muscle" = "#3C5488", "Adipose" = "#E64B35")), show_annotation_name = TRUE, annotation_name_rot = 45
)

ht <- Heatmap(
  z_mat, name = "Z-Score (Delta)", col = col_fun_heatmap,
  width = ncol(z_mat) * unit(6, "mm"), height = nrow(z_mat) * unit(6, "mm"),
  top_annotation = top_anno, left_annotation = left_anno,
  column_split = anno_df_hm$Fat_Group, row_split = row_anno_df$Tissue,
  column_gap = unit(2, "mm"), row_gap = unit(2, "mm"),
  cluster_columns = FALSE, cluster_rows = FALSE, show_column_names = FALSE, row_names_side = "right",
  row_names_gp = gpar(col = "black", fontsize = 10, fontface = "italic"),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"), row_title_gp = gpar(fontsize = 14, fontface = "bold"), border = TRUE
)

pdf(file.path(DIR_OUT, "Fig6C_cOC_Delta_Waterfall_Heatmap.pdf"), width = max(12, ncol(z_mat)*0.3+5), height = max(8, nrow(z_mat)*0.3+3))
draw(ht, merge_legend = TRUE, column_title = "Exercise-Induced Dynamic Responses (Post3H - Pre) vs cOC Gradient", column_title_gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()
message("  -> Fig 6C Heatmap Saved.")

# ==============================================================================
# === MODULE 2: FIG 6D & 6E GALAXY NETWORKS ===
# ==============================================================================
message("\n=== STEP 3: Rendering Fig 6D (Lean) & Fig 6E (Obese) Networks ===")

# Parameters
NODE_SIZE_cOC  <- 16; NODE_SIZE_SEED <- 11; NODE_SIZE_CLIN <- 10; NODE_SIZE_BASE <- 6         
LBL_SIZE_CORE <- 5.0; LBL_SIZE_NET  <- 3.5        
LBL_FORCE <- 0.5; LBL_POINT_PADDING <- 0.0; LBL_BOX_PADDING <- 0.2; LBL_SEGMENT_SIZE <- 0.3   

cat_colors <- c(
  "Center: Osteocalcin (cOC)" = "#DC0000", "Clinical Phenotype (tOC)" = "#4DBBD5",
  "Muscle Core Seed" = "#00A087", "Muscle Network" = "#71D0C4",   
  "Adipose Core Seed" = "#F39B7F", "Adipose Network" = "#FFC107",   
  "Adipose System (Fibrosis & Stress)" = "#F39B7F", "Muscle System (Heme & Oxygenation)" = "#00A087"
)

base_core_edges <- list(
  c("S_cOC", "Clin_S_TotalOC"), c("S_cOC", "PABPN1_Mus"), c("S_cOC", "FAM3C_Mus"),
  c("S_cOC", "CRACDL_Adi"), c("S_cOC", "TGFB1_Adi"), c("PABPN1_Mus", "FAM3C_Mus"), c("CRACDL_Adi", "TGFB1_Adi")
)
core_strings <- sapply(base_core_edges, function(x) paste(sort(x), collapse="--"))

target_cohorts <- c("Obese", "Lean")
excel_edges_list <- list(); excel_nodes_list <- list()

for(cohort in target_cohorts) {
  df_master <- df_master_global %>% dplyr::filter(Fat_Group == cohort) %>% dplyr::select(-Fat_Group)
  core_seeds <- intersect(c("S_cOC", "Clin_S_TotalOC", "PABPN1_Mus", "FAM3C_Mus", "CRACDL_Adi", "TGFB1_Adi"), colnames(df_master))
  
  recruited_nodes <- c(core_seeds)
  all_candidates <- setdiff(colnames(df_master), c("Sample_Key", core_seeds))
  valid_candidates <- all_candidates[sapply(df_master[all_candidates], function(x) sum(!is.na(x)) >= 10)]
  
  for(seed in core_seeds) {
    cors <- apply(df_master[valid_candidates], 2, function(x) {
      tmp <- na.omit(data.frame(A=df_master[[seed]], B=x))
      if(nrow(tmp) > 8) { ct <- cor.test(tmp$A, tmp$B, method="spearman", exact=FALSE); return(c(Rho = unname(ct$estimate), P = unname(ct$p.value))) } 
      else { return(c(Rho=0, P=1)) }
    })
    cor_df <- as.data.frame(t(cors)) %>% tibble::rownames_to_column("Target")
    top_targets <- cor_df %>% dplyr::filter(P < 0.05) %>% dplyr::arrange(desc(abs(Rho))) %>% head(12) %>% dplyr::pull(Target)
    recruited_nodes <- unique(c(recruited_nodes, top_targets))
  }
  
  mat_for_net <- df_master %>% dplyr::select(dplyr::all_of(recruited_nodes)) %>% as.matrix()
  rcorr_res <- Hmisc::rcorr(mat_for_net, type = "spearman")
  
  edges <- rcorr_res$r %>% as.data.frame() %>% tibble::rownames_to_column("Node1") %>% tidyr::pivot_longer(cols = -Node1, names_to = "Node2", values_to = "R") %>% dplyr::filter(Node1 < Node2)
  p_edges <- rcorr_res$P %>% as.data.frame() %>% tibble::rownames_to_column("Node1") %>% tidyr::pivot_longer(cols = -Node1, names_to = "Node2", values_to = "Pval") %>% dplyr::filter(Node1 < Node2)
  
  net_data <- edges %>% dplyr::inner_join(p_edges, by = c("Node1", "Node2")) %>% 
    dplyr::rowwise() %>% dplyr::mutate(Edge_ID = paste(sort(c(Node1, Node2)), collapse="--")) %>% dplyr::ungroup() %>%
    dplyr::filter((Pval < 0.05 & abs(R) > 0.45) | Edge_ID %in% core_strings) %>% 
    dplyr::mutate(Direction = ifelse(R > 0, "Positive", "Negative"), Plot_Weight = abs(R) * 2, Stress_Weight = abs(R))
  
  if(nrow(net_data) < 5) next
  
  g_raw <- igraph::graph_from_data_frame(d = net_data, directed = FALSE)
  comp <- igraph::components(g_raw)
  active_nodes <- names(comp$membership[comp$membership == which.max(comp$csize)])
  
  nodes_df <- data.frame(Name = active_nodes, stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      Clean_Label = gsub("Clin_S_TotalOC", "tOC", gsub("S_cOC", "cOC", gsub("_Mus|_Adi", "", Name))),
      Category = dplyr::case_when(Name == "S_cOC" ~ "Center: Osteocalcin (cOC)", Name == "Clin_S_TotalOC" ~ "Clinical Phenotype (tOC)", Name %in% c("PABPN1_Mus", "FAM3C_Mus") ~ "Muscle Core Seed", Name %in% c("CRACDL_Adi", "TGFB1_Adi") ~ "Adipose Core Seed", grepl("_Adi", Name) ~ "Adipose Network", grepl("_Mus", Name) ~ "Muscle Network", TRUE ~ "Other"),
      Hull_Group = dplyr::case_when(grepl("_Adi", Name) ~ "Adipose System (Fibrosis & Stress)", grepl("_Mus", Name) ~ "Muscle System (Heme & Oxygenation)", TRUE ~ NA_character_),
      Shape_Type = ifelse(Category == "Clinical Phenotype (tOC)", "Square", "Circle"),
      Node_Size = dplyr::case_when(Name == "S_cOC" ~ NODE_SIZE_cOC, Category == "Clinical Phenotype (tOC)" ~ NODE_SIZE_CLIN, Name %in% core_seeds ~ NODE_SIZE_SEED, TRUE ~ NODE_SIZE_BASE)
    )
  
  g <- igraph::graph_from_data_frame(d = net_data %>% dplyr::filter(Node1 %in% active_nodes & Node2 %in% active_nodes), vertices = nodes_df, directed = FALSE)
  igraph::V(g)$Degree <- igraph::degree(g)
  igraph::V(g)$Final_Size <- pmax(igraph::V(g)$Node_Size, igraph::V(g)$Degree * 0.4 + 3)
  
  set.seed(142) 
  p_net_full <- ggraph::ggraph(g, layout = 'stress', weights = Stress_Weight, bbox = 20) + 
    ggforce::geom_mark_hull(aes(x = x, y = y, filter = !is.na(Hull_Group), fill = Hull_Group, color = Hull_Group, label = Hull_Group), concavity = 1.5, expand = unit(6, "mm"), radius = unit(3, "mm"), alpha = 0.08, linetype = "dashed", linewidth = 0.8, label.fontsize = 13, label.fontface = "bold", label.fill = "white", label.buffer = unit(5, "mm"), con.colour = "gray50", con.cap = unit(0, "mm")) +
    ggraph::geom_edge_link(aes(edge_width = Plot_Weight, alpha = Plot_Weight, color = Direction)) +
    ggraph::scale_edge_color_manual(values = c("Positive" = "#E64B35", "Negative" = "#3C5488")) +
    ggraph::scale_edge_alpha_continuous(range = c(0.2, 0.8), guide = "none") + 
    ggraph::scale_edge_width_continuous(range = c(0.25, 1.25), guide = "none") +
    ggraph::geom_node_point(aes(size = Final_Size, fill = Category, shape = Shape_Type), color = "white", stroke = 1.2) +
    ggplot2::scale_shape_manual(values = c("Square" = 22, "Circle" = 21), guide = "none") + ggplot2::scale_size_identity() + 
    ggplot2::scale_fill_manual(values = cat_colors, name = "Node Categories") + ggplot2::scale_color_manual(values = cat_colors, guide = "none") + 
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 21, size = 6, color = "black"), ncol = 1)) +
    ggraph::geom_node_text(aes(label = Clean_Label, size = ifelse(Category %in% c("Center: Osteocalcin (cOC)", "Muscle Core Seed", "Adipose Core Seed", "Clinical Phenotype (tOC)"), LBL_SIZE_CORE, LBL_SIZE_NET), fontface = ifelse(Category %in% c("Center: Osteocalcin (cOC)", "Muscle Core Seed", "Adipose Core Seed", "Clinical Phenotype (tOC)"), "bold", "plain")), color = "black", repel = TRUE, force = LBL_FORCE, point.padding = LBL_POINT_PADDING, box.padding = LBL_BOX_PADDING, segment.color = "grey60", segment.size = LBL_SEGMENT_SIZE, min.segment.length = 0.5, bg.color = "white", bg.r = 0.15, max.overlaps = Inf) +
    ggplot2::labs(title = paste("Osteocalcin Multi-Organ Crosstalk (", cohort, " Cohort)", sep=""), subtitle = "Bright NPG Palette | Centered Labels | Narrow Edges") +
    ggplot2::theme_void(base_size = 14) + 
    ggplot2::theme(legend.position = "right", legend.title = element_text(face="bold", size=12), legend.text = element_text(size=11), plot.title = element_text(hjust = 0.5, face = "bold", size = 18, color = ifelse(cohort=="Obese", "#F39B7F", "#00A087")), plot.subtitle = element_text(hjust = 0.5, size = 13, color="gray40"), plot.margin = margin(10, 10, 10, 10))
  
  leg_plot <- cowplot::get_legend(p_net_full)
  if(!is.null(leg_plot)) ggsave(file.path(DIR_OUT, paste0("Fig6_Legend_Only_", cohort, ".pdf")), plot = leg_plot, width = 4, height = 6)
  
  fig_label <- ifelse(cohort == "Lean", "6D", "6E")
  ggsave(file.path(DIR_OUT, paste0("Fig", fig_label, "_Polished_Network_", cohort, ".pdf")), p_net_full + theme(legend.position = "none"), width = 12, height = 10)
  
  excel_edges_list[[cohort]] <- net_data %>% dplyr::arrange(desc(abs(R)))
  excel_nodes_list[[cohort]] <- nodes_df
  message("  -> Fig ", fig_label, " (", cohort, ") Network Saved.")
}

# ==============================================================================
# === MODULE 3: PUBLICATION-READY EXPORT (Data_S14 & Data_S15) ===
# ==============================================================================
message("\n=== STEP 4: Generating Supplementary Excel Tables ===")

if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)
header_style <- createStyle(fontName = "Arial", fontSize = 11, fontColour = "white", fgFill = "#4F81BD", textDecoration = "bold")

# Export Data_S14 (Heatmap)
wb_14 <- createWorkbook()
tidy_table <- as.data.frame(z_mat) %>% tibble::rownames_to_column("Feature_Clean") %>% tidyr::pivot_longer(cols = -Feature_Clean, names_to = "Subject_ID", values_to = "Z_Score_Delta") %>% dplyr::left_join(row_anno_df %>% dplyr::mutate(Feature_Clean = gsub("_Mus|_Adi", "", Feature)), by = "Feature_Clean") %>% dplyr::left_join(anno_df_hm, by = c("Subject_ID"="Sample_Key")) %>% dplyr::select(Tissue, Feature = Feature_Clean, Subject_ID, Fat_Group, cOC_Baseline = S_cOC, Z_Score_Delta, R_value, P_value, Significance) %>% dplyr::arrange(factor(Tissue, levels = c("Adipose", "Muscle")), desc(R_value), factor(Fat_Group, levels = c("Lean", "Obese")), cOC_Baseline)
summary_df <- row_anno_df %>% dplyr::mutate(Feature = gsub("_Mus|_Adi", "", Feature)) %>% dplyr::select(Tissue, Feature, R_value, P_value, Significance)

addWorksheet(wb_14, "Heatmap_Tidy_Data"); writeData(wb_14, "Heatmap_Tidy_Data", tidy_table)
addWorksheet(wb_14, "Gene_Stats_Summary"); writeData(wb_14, "Gene_Stats_Summary", summary_df)
for(sheet in c("Heatmap_Tidy_Data", "Gene_Stats_Summary")) { addStyle(wb_14, sheet, style = header_style, rows = 1, cols = 1:10, gridExpand = TRUE) }
saveWorkbook(wb_14, file.path(DIR_OUT, "Data_S14_cOC_Heatmap_Stats.xlsx"), overwrite = TRUE)

# Export Data_S15 (Network)
wb_15 <- createWorkbook()
for(g in target_cohorts) {
  addWorksheet(wb_15, paste0(g, "_Nodes")); writeData(wb_15, paste0(g, "_Nodes"), excel_nodes_list[[g]])
  addWorksheet(wb_15, paste0(g, "_Edges")); writeData(wb_15, paste0(g, "_Edges"), excel_edges_list[[g]])
  addStyle(wb_15, paste0(g, "_Nodes"), style = header_style, rows = 1, cols = 1:10, gridExpand = TRUE)
  addStyle(wb_15, paste0(g, "_Edges"), style = header_style, rows = 1, cols = 1:10, gridExpand = TRUE)
}
saveWorkbook(wb_15, file.path(DIR_OUT, "Data_S15_cOC_Network_Data.xlsx"), overwrite = TRUE)

message("\n=================================================================")
message("🏁 FIG 6C, 6D, 6E & DATA_S14, DATA_S15 COMPLETED PERFECTLY!")
message(">>> Output Directory: ", DIR_OUT)
message("=================================================================")
