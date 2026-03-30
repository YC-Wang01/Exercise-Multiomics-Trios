# ==============================================================================
# Project: FinlandSports (ATM) - Figure 3 Separate Networks
# Target: Tissue (Delta) -> Serum (Delta) -> Clinical Phenotype
# Features: Straight lines, Separate Legend, Distributed Hull Labels, Export CSVs
# ==============================================================================

# ------------------------------------------------------------------------------
# [0] Environment and Path Initialization
# ------------------------------------------------------------------------------
setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  limma, readr, dplyr, stringr, tibble, ggplot2, 
  igraph, ggraph, Hmisc, tidyr, cowplot, ggsci, ggforce
)

# V2.0 Standard Relative Paths
DIR_DATA <- "01_Clean_Data"
DIR_OUT  <- "03_Results/Fig_3/Fig3_Network"
if(!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)

# ------------------------------------------------------------------------------
# [1] Data Integration and Grouping
# ------------------------------------------------------------------------------
# In V2.0, use the unified strict clinical table instead of the raw sample sheet
clinical <- read_csv(file.path(DIR_DATA, "Clinical_Master_Strict.csv"), show_col_types = F)

meta_merged <- clinical %>% 
  dplyr::mutate(
    Role_Orig = dplyr::case_when(Membercode == 1 ~ "Daughter", Membercode == 2 ~ "Mother", Membercode == 3 ~ "Father", TRUE ~ "Unknown"),
    Subject_ID = tolower(paste0(FamilyID, "_", Role_Orig)),
    Fat_percent = as.numeric(`Fat(%)`),
    Group = dplyr::case_when(
      !is.na(Fat_percent) & Fat_percent >= 30 ~ "Obese",
      !is.na(Fat_percent) & Fat_percent < 30 ~ "Lean",
      TRUE ~ NA_character_
    )
  )

# ------------------------------------------------------------------------------
# [2] Delta Matrix Calculation (Post3h - Pre)
# ------------------------------------------------------------------------------
get_delta_matrix <- function(filename, label) {
  fpath <- file.path(DIR_DATA, filename)
  if(!file.exists(fpath)) return(NULL)
  
  data_raw <- read_csv(fpath, show_col_types = F)
  colnames(data_raw)[1] <- "FeatureID"
  
  expr_mat <- as.matrix(data_raw[, -1])
  rownames(expr_mat) <- data_raw$FeatureID
  clean_cols <- tolower(str_remove(colnames(expr_mat), "_twin\\.\\d+$")) %>% str_trim()
  
  # Parse column names into metadata dataframe
  col_meta <- data.frame(Matrix_Col = colnames(expr_mat), Clean_Col = clean_cols, stringsAsFactors = FALSE) %>%
    mutate(
      TimePoint = case_when(grepl("pre", Clean_Col) ~ "pre", grepl("fast", Clean_Col) ~ "fast", grepl("post3h", Clean_Col) ~ "post3h", TRUE ~ "Other"),
      Subject_ID = gsub("_fast|_pre|_post1h|_post3h", "", Clean_Col)
    )
  
  subjects <- unique(meta_merged$Subject_ID)
  delta_list <- list()
  
  for(s in subjects) {
    pre_ids <- col_meta %>% filter(Subject_ID == s, TimePoint %in% c("pre", "fast")) %>% pull(Matrix_Col)
    post_ids <- col_meta %>% filter(Subject_ID == s, TimePoint == "post3h") %>% pull(Matrix_Col)
    
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

mat_mus  <- get_delta_matrix("Cleaned_Muscle_Proteomics.csv", "Muscle")
mat_adi  <- get_delta_matrix("Cleaned_Adipose_Proteomics.csv", "Adipose") 
mat_met  <- get_delta_matrix("Cleaned_Serum_Metabonomics.csv", "Metabonomics")
mat_prot <- get_delta_matrix("Cleaned_Serum_Proteomics.csv", "Serum_Prot")

# ------------------------------------------------------------------------------
# [3] Core Plotting Function and Data Export
# ------------------------------------------------------------------------------
run_separate_network <- function(group_name) {
  message("\n========================================")
  message(">>> Generating network for: ", group_name)
  
  mats_list <- list(mat_mus, mat_adi, mat_met, mat_prot)
  valid_mats <- mats_list[!sapply(mats_list, is.null)]
  if(length(valid_mats) == 0) stop("FATAL ERROR: No omics data available!")
  
  common_subs <- Reduce(intersect, lapply(valid_mats, colnames))
  target_subs <- meta_merged %>% filter(Subject_ID %in% common_subs, Group == group_name) %>% pull(Subject_ID) %>% unique()
  
  if(length(target_subs) < 5) {
    message(" -> [Skip] Insufficient subjects (n < 5) for group: ", group_name)
    return(NULL)
  }
  
  select_top <- function(mat, n) {
    if(is.null(mat)) return(NULL)
    names(sort(apply(mat[, target_subs, drop=FALSE], 1, sd, na.rm=T), decreasing = T)[1:min(n, nrow(mat))])
  }
  
  t_mus  <- select_top(mat_mus, 15)
  t_adi  <- select_top(mat_adi, 15)
  t_met  <- select_top(mat_met, 20)
  t_prot <- select_top(mat_prot, 20)
  
  clin_traits <- c("HOMA_IR", "VO2max", "CRP", "Fat_percent", "TRIGLY")
  clin_mat <- meta_merged %>% 
    filter(Subject_ID %in% target_subs) %>% 
    select(Subject_ID, any_of(clin_traits)) %>% 
    distinct(Subject_ID, .keep_all = TRUE) %>% 
    column_to_rownames("Subject_ID") %>% t()
  
  comb_list <- list()
  if(!is.null(t_mus)) comb_list[["Muscle"]] <- mat_mus[t_mus, target_subs, drop=FALSE]
  if(!is.null(t_adi)) comb_list[["Adipose"]] <- mat_adi[t_adi, target_subs, drop=FALSE]
  if(!is.null(t_met)) comb_list[["Metab"]] <- mat_met[t_met, target_subs, drop=FALSE]
  if(!is.null(t_prot)) comb_list[["Prot"]] <- mat_prot[t_prot, target_subs, drop=FALSE]
  comb_list[["Clinical"]] <- clin_mat[, target_subs, drop=FALSE]
  
  comb_mat <- do.call(rbind, comb_list)
  rownames(comb_mat)[rownames(comb_mat) == "Fat_percent"] <- "Fat%"
  
  cor_res <- rcorr(t(comb_mat), type = "spearman")
  edges <- cor_res$r %>% as.data.frame() %>% rownames_to_column("Node1") %>%
    pivot_longer(-Node1, names_to = "Node2", values_to = "R") %>%
    inner_join(cor_res$P %>% as.data.frame() %>% rownames_to_column("Node1") %>%
                 pivot_longer(-Node1, names_to = "Node2", values_to = "P"), by=c("Node1","Node2")) %>%
    filter(Node1 < Node2, P < 0.05, abs(R) > 0.5) 
  
  node_info <- data.frame(name = rownames(comb_mat)) %>%
    mutate(
      Category = case_when(
        name %in% t_adi ~ "Adipose Tissue (Source)", 
        name %in% t_mus ~ "Muscle Tissue (Source)", 
        name %in% t_met ~ "Serum Metabolite (Mediator)", 
        name %in% t_prot ~ "Serum Protein (Mediator)",
        TRUE ~ "Clinical Phenotype (Target)" 
      ),
      Shape = ifelse(Category == "Clinical Phenotype (Target)", "Clinical", "Molecule")
    )
  
  g <- graph_from_data_frame(edges, directed = F, vertices = node_info)
  g <- delete_vertices(g, V(g)[degree(g) < 1]) 
  V(g)$Degree <- degree(g) # Save Degree directly into graph object
  
  layer_colors <- c(
    "Adipose Tissue (Source)" = "#F39B7F", 
    "Muscle Tissue (Source)"  = "#4DBBD5", 
    "Serum Metabolite (Mediator)" = "#8491B4", 
    "Serum Protein (Mediator)"    = "#91D1C2", 
    "Clinical Phenotype (Target)" = "#E64B35" 
  )
  
  # ================= Build Base Plot Including Legend =================
  p_base <- ggraph(g, layout = 'kk') +  
    
    # Dashed hull annotations: using label.margin and label.buffer
    ggforce::geom_mark_hull(
      aes(x, y, group = Category, fill = Category, color = Category, label = Category), 
      concavity = 2, expand = unit(4, "mm"), radius = unit(3, "mm"), 
      alpha = 0.08, linetype = "dashed", linewidth = 0.6,
      label.fontsize = 9, label.fontface = "bold", label.fill = "white",
      label.margin = margin(5, 5, 5, 5), label.buffer = unit(10, "mm"), 
      con.colour = "grey50", con.cap = unit(1, "mm")
    ) +
    
    # Straight edges
    geom_edge_link(
      aes(edge_width = abs(R), color = R > 0, linetype = R > 0), 
      alpha = 0.6
    ) +
    scale_edge_color_manual(
      values = c("TRUE" = "#E64B35", "FALSE" = "#3C5488"), 
      labels = c("TRUE" = "Positive (+)", "FALSE" = "Negative (-)"),
      name = "Interaction Trend"
    ) +
    scale_edge_width(range = c(0.4, 1.8), name = "Correlation |R|") + 
    scale_edge_linetype_manual(values = c("TRUE"="solid", "FALSE"="dashed"), guide = "none") +
    
    # Node rendering
    geom_node_point(
      aes(fill = Category, shape = Shape, size = Degree), 
      color = "white", stroke = 1.2
    ) +
    scale_fill_manual(values = layer_colors, name = "Biological Layer") + 
    scale_color_manual(values = layer_colors, guide = "none") +
    scale_shape_manual(values = c("Molecule" = 21, "Clinical" = 23), guide = "none") +
    scale_size_continuous(range = c(5, 14), name = "Node Degree") +
    
    # Node text
    geom_node_text(
      aes(label = name), repel = TRUE, size = 3.5, 
      family = "sans", fontface = "bold", 
      bg.color = "white", bg.r = 0.15, segment.color = "grey60", segment.size = 0.5
    ) +
    
    coord_cartesian(clip = "off") + 
    theme_void() + 
    theme(
      text = element_text(family = "sans"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b=10)),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40", margin = margin(b=20)),
      plot.margin = margin(30, 30, 30, 30), 
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "right",          
      legend.box = "vertical",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10)
    ) +
    labs(
      title = paste0("Tissue-Serum-Clinical Crosstalk: ", group_name),
      subtitle = "Multi-omics integration mapped by functional layers"
    )
  
  # ================= Extract and Save Legend and Main Plot =================
  legend_obj <- cowplot::get_legend(p_base)
  p_main <- p_base + theme(legend.position = "none")
  
  file_main <- file.path(DIR_OUT, paste0("Fig3_Network_", group_name, "_MainPlot.pdf"))
  ggsave(file_main, p_main, width = 12, height = 10, dpi = 300)
  
  file_legend <- file.path(DIR_OUT, paste0("Fig3_Network_", group_name, "_LegendOnly.pdf"))
  if(!is.null(legend_obj)) {
    ggsave(file_legend, cowplot::ggdraw(legend_obj), width = 4, height = 7, dpi = 300)
  }
  
  message(" -> Main plot and standalone legend saved successfully.")
  
  # ================= Export Network Data as Supplementary Tables =================
  out_edges <- edges %>%
    arrange(P) %>% 
    mutate(
      Group = group_name,
      R = round(R, 4),       
      P = signif(P, 4),      
      Trend = ifelse(R > 0, "Positive", "Negative") 
    ) %>%
    select(Group, Node1, Node2, Correlation_R = R, P_value = P, Trend)
  
  file_edges <- file.path(DIR_OUT, paste0("Table_Network_Edges_", group_name, ".csv"))
  write_csv(out_edges, file_edges)
  
  # Extract node information directly from the rendered graph object
  out_nodes <- data.frame(
    Molecule_Name = V(g)$name,
    Biological_Layer = V(g)$Category,
    Degree = V(g)$Degree
  ) %>%
    arrange(Biological_Layer, desc(Degree)) %>% 
    mutate(Group = group_name) %>%
    select(Group, Molecule_Name, Biological_Layer, Degree)
  
  file_nodes <- file.path(DIR_OUT, paste0("Table_Network_Nodes_", group_name, ".csv"))
  write_csv(out_nodes, file_nodes)
  
  message(" -> Network nodes and edges exported as CSV supplementary tables.")
}

# ------------------------------------------------------------------------------
# [4] Execute Analysis
# ------------------------------------------------------------------------------
run_separate_network("Obese")
run_separate_network("Lean")

message("\n>>> Figure 3 components generation completed successfully! All results exported.")