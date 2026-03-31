# ==============================================================================
# Project: FinlandSports V2.0 - Fig 5F & 5G: CRP Galaxy Networks
# Description:
#   1. Strict Filter: |R| >= 0.55 and P < 0.05 (Spearman).
#   2. Core Highlight: CRP, ZNF76, TRAPPC13, ACTB, JUP, C1S, PLTP.
#   3. Assignment: Lean mapped to Fig 5F, Obese mapped to Fig 5G.
#   4. V2.0 Integration: Automatically bridges with Limma results from Fig 5D/E.
# ==============================================================================

# --- 0. Environment Setup ---
# setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readr, dplyr, stringr, tidyr, tibble, openxlsx,
  igraph, ggraph, tidygraph, Hmisc, ggforce, concaveman,
  ggplot2, ggsci
)

# --- 1. PATH CONFIGURATION (V2.0) ---
DIR_IN   <- "01_Clean_Data"
DIR_PREV <- "03_Results/Fig_5/Fig5DE_CRP_Downstream/Results_Tables"
DIR_OUT  <- "03_Results/Fig_5/Fig5FG_CRP_Networks"

if(!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)

# --- 2. LOAD METADATA & CRP ---
message("\n>>> STEP 1: Loading Clinical Metadata...")
clinical <- read_csv(file.path(DIR_IN, "Clinical_Master_Strict.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(
    Role_Orig = dplyr::case_when(Membercode == 1 ~ "Daughter", Membercode == 2 ~ "Mother", Membercode == 3 ~ "Father", TRUE ~ "Unknown"),
    Sample_Key = tolower(paste0(FamilyID, "_", Role_Orig))
  )

meta_merged <- clinical %>%
  dplyr::rename(CRP = CRP) %>% # Replace with specific CRP column name if needed
  dplyr::filter(!is.na(CRP), !is.na(`Fat(%)`)) %>%
  dplyr::mutate(
    Log_CRP = log(CRP),
    Fat_Group = ifelse(`Fat(%)` >= 30, "Obese", "Lean")
  )

# --- 3. DEFINE CORE HIGHLIGHT GENES (The 7 Pillars) ---
core_highlights <- c("CRP", "ZNF76", "TRAPPC13", "ACTB", "JUP", "C1S", "PLTP")

# --- 4. LOAD TOP 50 GENES FROM PREVIOUS LIMMA RUN ---
message(">>> STEP 2: Loading Top 50 Significant Genes from Limma...")
TARGETS <- c("Adipose", "Muscle", "Serum")
target_genes <- list()
clean_dict <- c()
tissue_dict <- c()

for(t in TARGETS) {
  res_file <- file.path(DIR_PREV, paste0("Stats_Limma_", t, ".xlsx"))
  if(!file.exists(res_file)) {
    message("  [Warning] Limma results not found for ", t, ". Ensure Fig5D/E script was run.")
    next
  }
  
  res <- openxlsx::read.xlsx(res_file)
  top50 <- res %>% dplyr::filter(P.Value < 0.05) %>% dplyr::arrange(P.Value) %>% head(50)
  
  for(i in 1:nrow(top50)) {
    orig_name <- top50$Protein[i]
    clean_name <- stringr::str_replace(orig_name, ";.*", "")
    target_genes[[t]] <- c(target_genes[[t]], orig_name)
    clean_dict[orig_name] <- clean_name
    tissue_dict[clean_name] <- t
  }
}

# --- 5. EXTRACT EXPRESSION MATRIX ---
message(">>> STEP 3: Extracting Omics Baseline Expression Matrices...")

parse_sample_v2 <- function(col_name) {
  x_clean <- tolower(stringr::str_remove(col_name, "_twin\\.\\d+$")) %>% stringr::str_trim()
  if(stringr::str_detect(x_clean, "post")) return(NULL) 
  s_key <- gsub("_fast|_pre|_post1h|_post3h", "", x_clean)
  
  idx <- which(meta_merged$Sample_Key == s_key)
  if(length(idx) != 1) return(NULL)
  
  list(OriginalCol=col_name, UniqueID=s_key)
}

mat_list <- list()
file_targets <- list(Adipose = "Cleaned_Adipose_Proteomics.csv", Muscle = "Cleaned_Muscle_Proteomics.csv", Serum = "Cleaned_Serum_Proteomics.csv")

for(t in names(file_targets)) {
  if(is.null(target_genes[[t]])) next
  fpath <- file.path(DIR_IN, file_targets[[t]])
  if(!file.exists(fpath)) next
  
  raw <- readr::read_csv(fpath, show_col_types = FALSE)
  feat_ids <- raw[[1]]
  expr_raw <- as.matrix(raw[, -1])
  
  col_infos <- lapply(colnames(expr_raw), parse_sample_v2)
  valid_idx <- which(!sapply(col_infos, is.null))
  
  col_meta <- dplyr::bind_rows(col_infos[valid_idx])
  expr_valid <- expr_raw[, valid_idx]
  mat <- suppressWarnings(as.data.frame(lapply(as.data.frame(expr_valid), function(x) as.numeric(as.character(x)))))
  rownames(mat) <- feat_ids
  
  genes_to_keep <- intersect(target_genes[[t]], rownames(mat))
  mat_sub <- mat[genes_to_keep, , drop=FALSE]
  colnames(mat_sub) <- col_meta$UniqueID
  
  rownames(mat_sub) <- clean_dict[rownames(mat_sub)]
  mat_list[[t]] <- mat_sub
}

subs_adi <- if(!is.null(mat_list[["Adipose"]])) colnames(mat_list[["Adipose"]]) else character(0)
subs_mus <- if(!is.null(mat_list[["Muscle"]])) colnames(mat_list[["Muscle"]]) else character(0)
subs_ser <- if(!is.null(mat_list[["Serum"]])) colnames(mat_list[["Serum"]]) else character(0)
valid_tissue_subs <- intersect(intersect(subs_adi, subs_mus), subs_ser)

combined_mat <- matrix(NA, nrow=0, ncol=length(valid_tissue_subs), dimnames=list(NULL, valid_tissue_subs))
for(m in mat_list) { combined_mat <- rbind(combined_mat, m[, valid_tissue_subs]) }

crp_row <- matrix(meta_merged$CRP[match(valid_tissue_subs, meta_merged$Sample_Key)], nrow=1, dimnames=list("CRP", valid_tissue_subs))
combined_mat <- rbind(combined_mat, crp_row)
tissue_dict["CRP"] <- "Clinical"

# --- 6. NETWORK CONSTRUCTION FUNCTION ---
build_network <- function(mat, group_name) {
  message(paste0("   -> Calculating Spearman Correlation for: ", group_name))
  mat_f <- mat[rowSums(!is.na(mat)) > 5, ] 
  
  cor_res <- Hmisc::rcorr(t(mat_f), type = "spearman")
  r_mat <- cor_res$r; p_mat <- cor_res$P
  
  edges <- which(upper.tri(r_mat), arr.ind = TRUE)
  edge_df <- data.frame(
    from = rownames(r_mat)[edges[,1]], to = rownames(r_mat)[edges[,2]],
    R = r_mat[edges], P = p_mat[edges], stringsAsFactors = FALSE
  )
  
  # Strict filter: |R| >= 0.55 & P < 0.05
  edge_df <- edge_df %>% dplyr::filter(abs(R) >= 0.55 & P < 0.05)
  if(nrow(edge_df) == 0) return(list(graph = NULL, edges = NULL))
  
  edge_df$Weight <- abs(edge_df$R)
  edge_df$Direction <- ifelse(edge_df$R > 0, "Positive", "Negative")
  
  graph <- igraph::graph_from_data_frame(edge_df, directed = FALSE)
  igraph::V(graph)$Tissue <- tissue_dict[igraph::V(graph)$name]
  igraph::V(graph)$Degree <- igraph::degree(graph)
  igraph::V(graph)$Is_Core <- ifelse(igraph::V(graph)$name %in% core_highlights, TRUE, FALSE)
  igraph::V(graph)$Shape <- ifelse(igraph::V(graph)$name == "CRP", "square", "circle")
  
  return(list(graph = graph, edges = edge_df %>% dplyr::arrange(dplyr::desc(Weight))))
}

# --- 7. PLOT NETWORKS & EXPORT COMBINED EXCEL ---
message("\n>>> STEP 4: Generating Fig 5F (Lean) & Fig 5G (Obese) Networks...")
groups <- c("Lean", "Obese")
excel_export_list <- list() 

for(g in groups) {
  # Panel Assignment logic
  panel_label <- ifelse(g == "Lean", "F", "G")
  
  subs_g <- meta_merged %>% dplyr::filter(Sample_Key %in% valid_tissue_subs & Fat_Group == g) %>% dplyr::pull(Sample_Key)
  if(length(subs_g) < 5) {
    message("  [Warning] Not enough samples for ", g, " group.")
    next
  }
  
  mat_g <- combined_mat[, subs_g, drop=FALSE]
  net_res <- build_network(mat_g, g)
  
  if(is.null(net_res$graph)) {
    message("  [Warning] No significant edges found for ", g, " network.")
    next
  }
  
  g_graph <- net_res$graph
  excel_export_list[[g]] <- net_res$edges 
  
  node_colors <- c("Adipose"="#F39B7F", "Muscle"="#4DBBD5", "Serum"="#00A087", "Clinical"="#E64B35")
  edge_colors <- c("Positive"="#E64B35", "Negative"="#3C5488")
  
  set.seed(42) 
  layout_pos <- ggraph::create_layout(g_graph, layout = "fr")
  
  p_net <- ggraph::ggraph(layout_pos) + 
    ggforce::geom_mark_hull(aes(x = x, y = y, fill = Tissue, color = Tissue, label = Tissue, filter = Tissue != "Clinical"),
                            concavity = 2, expand = unit(4, "mm"), alpha = 0.15,
                            linetype = "dashed", linewidth = 0.6, 
                            label.fontsize = 12, label.fontface = "bold", label.fill = alpha("white", 0.8),
                            con.cap = 0, show.legend = FALSE) +
    
    ggraph::geom_edge_link(aes(edge_alpha = Weight, edge_width = Weight, color = Direction), show.legend = c(edge_width=FALSE, edge_alpha=FALSE)) +
    ggraph::scale_edge_width(range = c(0.1, 0.6)) +
    ggraph::scale_edge_alpha(range = c(0.4, 0.8)) +
    ggraph::scale_edge_color_manual(values = edge_colors) +
    
    ggraph::geom_node_point(aes(size = Is_Core, fill = Tissue, shape = Shape), color = "white", stroke = 1.0) +
    ggplot2::scale_shape_manual(values = c("circle" = 21, "square" = 22), guide="none") +
    ggplot2::scale_fill_manual(values = node_colors) +
    ggplot2::scale_color_manual(values = node_colors) + 
    ggplot2::scale_size_manual(values = c("TRUE" = 8, "FALSE" = 3.5), guide="none") + 
    
    # Generic Labels (Black)
    ggraph::geom_node_text(aes(label = ifelse(!Is_Core, name, "")), 
                           repel = TRUE, point.padding = 0.1, force = 0.5,
                           size = 2.8, fontface = "bold", color = "black",
                           bg.color = "white", bg.r = 0.1, max.overlaps = 150) +
    
    # Core Highlights (Black, larger)
    ggraph::geom_node_text(aes(label = ifelse(Is_Core, name, "")), 
                           repel = TRUE, point.padding = 0.3, force = 1.5,
                           size = 4.8, fontface = "bold", color = "black",
                           bg.color = "white", bg.r = 0.15, max.overlaps = 50) +
    
    ggplot2::theme_void(base_size = 14) +
    ggplot2::labs(
      title = paste0("Fig 5", panel_label, ": ", g, " Systemic Crosstalk"),
      subtitle = "Nodes: Bright NPG | Hulls: Dashed Polygons | Edges: |R| \u2265 0.55, P < 0.05"
    ) +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(face="bold", size=12),
      plot.title = ggplot2::element_text(face="bold", hjust=0.5, size=18, color="black"),
      plot.subtitle = ggplot2::element_text(hjust=0.5, size=11, color="grey50"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(15, 15, 15, 15)
    )
  
  ggplot2::ggsave(file.path(DIR_OUT, paste0("Fig5", panel_label, "_Network_", g, ".pdf")), p_net, width = 9, height = 7)
  message(paste0("   -> Fig 5", panel_label, " (", g, ") Network saved."))
}

# --- 8. EXPORT EDGES DATA ---
if(length(excel_export_list) > 0) {
  openxlsx::write.xlsx(excel_export_list, file.path(DIR_OUT, "Data_S10_CRP_Network_Edges.xlsx"))
  message("\n   -> Edge tables combined and saved to Data_S10_CRP_Network_Edges.xlsx")
}

message("\n=======================================================")
message(">>> FIG 5 F & G NETWORK ANALYSIS COMPLETE!")
message(">>> Check Folder: ", DIR_OUT)
message("=======================================================")
