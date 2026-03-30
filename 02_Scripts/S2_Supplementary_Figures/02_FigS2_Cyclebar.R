# ==============================================================================
# Project: FinlandSports V2.0
# Script:  02_FigS2_Cyclebar.R
# Author:  Gemini & Sorcier_W
# Description: 
#   Merges circular omics baseline differential analysis for Serum, Adipose, 
#   and Muscle tissues (Fig. S2). Includes automated differential analysis (Limma), 
#   Nature-style circular plotting, and legend generation.
# Features: V2.0 I/O Standardization, 100% Original Logic Replication.
# ==============================================================================

# ==============================================================================
# [0] PATHS & CONFIGURATION
# ==============================================================================
setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")
CONF_OUT_ROOT   <- "03_Results/Fig_S2/Cyclebar"
CONF_OMICS_DIR  <- "01_Clean_Data"
CONF_META_FILE  <- file.path(CONF_OMICS_DIR, "Clinical_Master_Strict.csv")
message("========================================================")
message(">>> [INPUT]  Meta File: ", CONF_META_FILE)
message(">>> [INPUT]  Omics Dir: ", CONF_OMICS_DIR)
message(">>> [OUTPUT] Results:   ", CONF_OUT_ROOT)
message("========================================================")

PLOT_CONFIG <- list(
  COLOR_UP_MAX        = "#CC0000", 
  COLOR_DN_MAX        = "#00008B", 
  COLOR_MIN_INTENSITY = 0.2,       
  
  SORT_BY             = "PValue", 
  R_IN                = 0.5,   
  R_WIDTH_TISSUE      = 4.5,   
  GAP_RING_BAR        = 0.15,  
  CENTER_HOLE_SIZE    = 1.5, 
  HEIGHT_SCALE        = 3.5,   
  BAR_WIDTH           = 0.9,   
  MAX_FC_CAP          = 2.5,   
  
  LABEL_SIZE          = 3.0,   
  LABEL_FONT_FACE     = "plain",
  SHOW_GRID           = FALSE, 
  SHOW_BAR_BORDER     = FALSE  
)

# Task definition: Tissue -> Filename -> Sector colors -> Order
TISSUE_TASKS <- list(
  "Serum" = list(
    files = c(Metabolomics = "Cleaned_Serum_Metabonomics.csv", Proteomics = "Cleaned_Serum_Proteomics.csv"),
    colors = c(Metabolomics = "#F39B7F", Proteomics = "#4DBBD5"),
    order = c("Metabolomics", "Proteomics"),
    total_slots = 110
  ),
  "Adipose" = list(
    files = c(Methylation = "Cleaned_Adipose_Methylation.csv", Transcriptome = "Cleaned_Adipose_Microarray.csv", Proteomics = "Cleaned_Adipose_Proteomics.csv"),
    colors = c(Methylation = "#E64B35", Transcriptome = "#00A087", Proteomics = "#4DBBD5"),
    order = c("Methylation", "Transcriptome", "Proteomics"),
    total_slots = 110
  ),
  "Muscle" = list(
    files = c(Methylation = "Cleaned_Muscle_Methylation.csv", Transcriptome = "Cleaned_Muscle_Microarray.csv", Proteomics = "Cleaned_Muscle_Proteomics.csv"),
    colors = c(Methylation = "#E64B35", Transcriptome = "#00A087", Proteomics = "#4DBBD5"),
    order = c("Methylation", "Transcriptome", "Proteomics"),
    total_slots = 110
  )
)

# ==============================================================================
# [1] Environment Initialization
# ==============================================================================
if (!require("pacman")) install.packages("pacman")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("ComplexHeatmap")) BiocManager::install("ComplexHeatmap")

pacman::p_load(readr, dplyr, stringr, limma, tibble, ggplot2, ggrepel, ggsci, 
               grid, cowplot, scales, openxlsx, ComplexHeatmap, circlize)

tryCatch(library(shadowtext), error=function(e) NULL) 

CONF_FONT_FAMILY <- "sans"
dir.create(CONF_OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# [2] Data Loading and General Analysis Functions (V2.0 logic)
# ==============================================================================
if(!file.exists(CONF_META_FILE)) stop(paste("Meta file not found:", CONF_META_FILE))
meta_basic <- read_csv(CONF_META_FILE, show_col_types = FALSE) %>%
  dplyr::filter(!is.na(`Fat(%)`)) %>%
  dplyr::mutate(
    Role_Orig = dplyr::case_when(Membercode == 1 ~ "Daughter", Membercode == 2 ~ "Mother", Membercode == 3 ~ "Father", TRUE ~ "Unknown"),
    Sample_Key = tolower(paste0(FamilyID, "_", Role_Orig)),
    Group = ifelse(`Fat(%)` >= 30, "Obese", "Lean"),
    Sex = factor(Gender),
    Age = as.numeric(Age)
  )

parse_sample <- function(col_name) {
  x_clean <- tolower(str_remove(col_name, "\\.\\.\\.\\d+$") %>% str_trim())
  if(grepl("post", x_clean)) return(NULL) 
  
  s_key <- gsub("_fast|_pre|_post1h|_post3h", "", x_clean)
  idx <- which(meta_basic$Sample_Key == s_key)
  if(length(idx) != 1) return(NULL)
  return(meta_basic[idx, ])
}

run_limma <- function(omics_type, filename) {
  fpath <- file.path(CONF_OMICS_DIR, filename)
  if(!file.exists(fpath)) { warning(paste("File missing:", filename)); return(NULL) }
  
  raw <- read_csv(fpath, show_col_types = FALSE)
  feat_ids <- raw[[1]]; expr <- raw[,-1]
  
  col_infos <- lapply(colnames(expr), parse_sample)
  valid_idx <- which(!sapply(col_infos, is.null))
  if(length(valid_idx) < 10) return(NULL)
  
  mat <- suppressWarnings(as.data.frame(lapply(expr[,valid_idx], function(x) as.numeric(as.character(x)))))
  rownames(mat) <- feat_ids
  
  mat <- mat[rowSums(is.na(mat)) < (ncol(mat)*0.5), ]
  mat[is.na(mat)] <- min(mat, na.rm=T)/2
  
  if(omics_type != "Methylation" && max(mat, na.rm=T) > 50) { mat <- log2(mat + 1) }
  
  col_meta <- bind_rows(col_infos[valid_idx])
  Group <- factor(col_meta$Group, levels=c("Lean","Obese")); Age <- col_meta$Age
  design <- if(length(unique(col_meta$Sex))>1) model.matrix(~0+Group+Age+col_meta$Sex) else model.matrix(~0+Group+Age)
  colnames(design)[1:2] <- c("Lean", "Obese"); colnames(design) <- make.names(colnames(design))
  
  fit <- lmFit(mat, design)
  fit <- eBayes(contrasts.fit(fit, makeContrasts(Diff=Obese-Lean, levels=design)))
  res <- topTable(fit, coef="Diff", number=Inf) %>% 
    rownames_to_column("Feature") %>% 
    mutate(Omics = omics_type)
  
  res$Feature <- gsub(";.*", "", res$Feature)
  return(list(res = res))
}

# ==============================================================================
# [3] Core Plotting Function (1:1 Replication of Visual Mapping)
# ==============================================================================
draw_custom_cycle <- function(data_store, cfg, tissue_cfg) {
  total_slots <- tissue_cfg$total_slots
  r_in <- cfg$R_IN; r_mid <- r_in + cfg$R_WIDTH_TISSUE; r_out <- r_mid + cfg$GAP_RING_BAR
  bar_border_col <- if(cfg$SHOW_BAR_BORDER) "black" else NA
  
  mix_white <- function(col, intensity) { colorRampPalette(c("white", col))(100)[max(1, round(intensity * 100))] }
  col_up_min <- mix_white(cfg$COLOR_UP_MAX, cfg$COLOR_MIN_INTENSITY)
  col_dn_min <- mix_white(cfg$COLOR_DN_MAX, cfg$COLOR_MIN_INTENSITY)
  
  stats_list <- lapply(data_store, function(x) if(!is.null(x)) nrow(filter(x$res, P.Value < 0.05)) else 0)
  total_sig <- sum(unlist(stats_list))
  if(total_sig == 0) return(NULL)
  
  processed_list <- list()
  for(t in tissue_cfg$order) {
    if(!t %in% names(data_store)) next
    df <- data_store[[t]]$res %>% filter(P.Value < 0.05)
    if(nrow(df)==0) next
    
    ratio <- nrow(df) / total_sig
    n_slots <- max(15, round(total_slots * ratio))
    
    n_up <- sum(df$logFC > 0); n_dn <- sum(df$logFC < 0)
    ratio_up <- n_up / (n_up + n_dn)
    n_up_slots <- round((n_slots-3) * ratio_up)
    n_dn_slots <- (n_slots-3) - n_up_slots
    if(n_up > 0 && n_up_slots < 2) n_up_slots <- 2
    if(n_dn > 0 && n_dn_slots < 2) n_dn_slots <- 2
    
    if (cfg$SORT_BY == "PValue") {
      top_up <- df %>% filter(logFC > 0) %>% arrange(P.Value) %>% head(n_up_slots) %>% arrange(P.Value)
      top_dn <- df %>% filter(logFC < 0) %>% arrange(P.Value) %>% head(n_dn_slots) %>% arrange(desc(P.Value))
    } else {
      top_up <- df %>% filter(logFC > 0) %>% arrange(desc(logFC)) %>% head(n_up_slots) %>% arrange(desc(logFC))
      top_dn <- df %>% filter(logFC < 0) %>% arrange(logFC) %>% head(n_dn_slots) %>% arrange(desc(logFC))
    }
    
    if(nrow(top_up) > 0) top_up$Type <- "Up"
    if(nrow(top_dn) > 0) top_dn$Type <- "Down"
    
    processed_list[[t]] <- bind_rows(top_up, data.frame(Feature=c(".","."), logFC=0, P.Value=1, Omics=t, Type="Gap"), top_dn)
  }
  
  plot_df <- bind_rows(processed_list) %>% mutate(id = 1:n())
  
  plot_df$HeightRaw <- abs(plot_df$logFC)
  plot_df$HeightRaw[plot_df$Type == "Gap"] <- 0
  plot_df$HeightRaw <- ifelse(plot_df$HeightRaw > cfg$MAX_FC_CAP, cfg$MAX_FC_CAP, plot_df$HeightRaw)
  plot_df$HeightRaw[plot_df$Type != "Gap" & plot_df$HeightRaw < 0.1] <- 0.1
  plot_df$HeightScaled <- plot_df$HeightRaw * cfg$HEIGHT_SCALE
  plot_df$NegLogP <- -log10(plot_df$P.Value)
  
  real_data <- plot_df %>% filter(Type != "Gap")
  min_p <- min(real_data$NegLogP, na.rm=T); max_p <- max(real_data$NegLogP, na.rm=T)
  
  plot_df$Fill_Color <- NA
  if(any(plot_df$Type == "Up", na.rm=T)) {
    idx <- which(plot_df$Type == "Up")
    plot_df$Fill_Color[idx] <- scales::col_numeric(c(col_up_min, cfg$COLOR_UP_MAX), domain=c(min_p, max_p))(plot_df$NegLogP[idx])
  }
  if(any(plot_df$Type == "Down", na.rm=T)) {
    idx <- which(plot_df$Type == "Down")
    plot_df$Fill_Color[idx] <- scales::col_numeric(c(col_dn_min, cfg$COLOR_DN_MAX), domain=c(min_p, max_p))(plot_df$NegLogP[idx])
  }
  
  n_bars <- nrow(plot_df)
  angle <- 90 - 360 * (plot_df$id - 0.5) / n_bars
  plot_df$hjust <- ifelse(angle < -90, 1, 0)
  plot_df$angle <- ifelse(angle < -90, angle + 180, angle)
  
  tissue_segs <- plot_df %>% group_by(Omics) %>% 
    summarize(start_id = min(id), end_id = max(id), mid_id = mean(c(min(id), max(id)))) %>% 
    mutate(col = tissue_cfg$colors[Omics])
  
  p <- ggplot(plot_df) +
    geom_rect(data=tissue_segs, aes(xmin=start_id-0.5, xmax=end_id+0.5, ymin=r_in, ymax=r_mid, fill=col), color=NA) +
    geom_rect(data=filter(plot_df, Type!="Gap"), aes(xmin=id-cfg$BAR_WIDTH/2, xmax=id+cfg$BAR_WIDTH/2, ymin=r_out, ymax=r_out+HeightScaled, fill=Fill_Color), color=bar_border_col, size=0.05) +
    geom_point(data=filter(plot_df, Type=="Gap"), aes(x=id, y=r_out + 0.5*cfg$HEIGHT_SCALE), shape=16, size=1, color="grey40") +
    geom_text(data=filter(plot_df, Type!="Gap"), aes(x=id, y=r_out+HeightScaled+0.5, label=Feature, angle=angle, hjust=hjust), size=cfg$LABEL_SIZE, fontface=cfg$LABEL_FONT_FACE, family=CONF_FONT_FAMILY) +
    geom_text(data=tissue_segs, aes(x=mid_id, y=(r_in+r_mid)/2, label=Omics), size=4, fontface="bold", color="white", family=CONF_FONT_FAMILY) +
    scale_fill_identity() + coord_polar() + theme_void() +
    scale_y_continuous(limits = c(-1 * cfg$CENTER_HOLE_SIZE, r_out + (cfg$MAX_FC_CAP * cfg$HEIGHT_SCALE) + 4))
  
  return(list(plot=p, min_p=min_p, max_p=max_p, data=plot_df))
}

# ==============================================================================
# [4] Legend Generation (Separate Legends)
# ==============================================================================
draw_separate_legends <- function(min_p, max_p, cfg, out_dir, t_name) {
  col_dn_min <- colorRampPalette(c("white", cfg$COLOR_DN_MAX))(100)[20]
  col_up_min <- colorRampPalette(c("white", cfg$COLOR_UP_MAX))(100)[20]
  
  # 1. Colorbar Legend
  breaks_vec <- c(-max_p, -min_p, 0, min_p, max_p)
  colors_vec <- c(cfg$COLOR_DN_MAX, col_dn_min, "white", col_up_min, cfg$COLOR_UP_MAX)
  col_fun_combined <- colorRamp2(breaks_vec, colors_vec)
  at_display <- c(-max_p, -min_p, min_p, max_p)
  labels_expr <- parse(text = paste0("10^-", round(abs(at_display), 1)))
  
  lgd_color <- Legend(col_fun = col_fun_combined, title = expression(bold("Significance") ~ (-log[10]*P)), 
                      at = at_display, labels = labels_expr, direction = "horizontal", 
                      legend_width = unit(6, "cm"), border = "black")
  
  pdf(file.path(out_dir, paste0("FigS2_Legend_Colorbar_", t_name, ".pdf")), width = 4, height = 2)
  draw(lgd_color)
  dev.off()
  
  # 2. Ruler Legend
  df_ruler <- data.frame(y = 0:3)
  p_ruler <- ggplot(df_ruler, aes(x = 0, y = y)) +
    geom_segment(aes(x = 0, xend = 0, y = 0, yend = 3), linewidth = 1) +
    geom_segment(aes(x = 0, xend = 0.2, y = y, yend = y), linewidth = 0.8) +
    geom_text(aes(x = 0.3, y = y, label = y), hjust = 0, size = 5, family = CONF_FONT_FAMILY) +
    annotate("text", x = 0, y = 3.3, label = "Log2 FC", fontface = "bold", size = 5, hjust=0, family = CONF_FONT_FAMILY) +
    scale_x_continuous(limits = c(0, 1.5), expand = c(0,0)) + scale_y_continuous(limits = c(0, 3.5), expand = c(0,0)) + theme_void()
  
  ggsave(file.path(out_dir, paste0("FigS2_Legend_Ruler_", t_name, ".pdf")), p_ruler, width = 2, height = 4)
}

# ==============================================================================
# [5] Main Loop Execution
# ==============================================================================
for(t_name in names(TISSUE_TASKS)) {
  message("\n########################################################")
  message(">>> Processing Tissue: ", t_name)
  task <- TISSUE_TASKS[[t_name]]
  out_dir <- file.path(CONF_OUT_ROOT, paste0("ColorCycle_", t_name))
  dir.create(file.path(out_dir, "Results_Tables"), recursive=T, showWarnings=F)
  
  results <- list()
  for(omics_n in names(task$files)) {
    message("   - Analyzing: ", omics_n)
    res_obj <- run_limma(omics_n, task$files[[omics_n]])
    if(!is.null(res_obj)) {
      results[[omics_n]] <- res_obj
      write.xlsx(res_obj$res, file.path(out_dir, "Results_Tables", paste0("Limma_", omics_n, "_ALL.xlsx")))
      write.xlsx(filter(res_obj$res, P.Value < 0.05), file.path(out_dir, "Results_Tables", paste0("Limma_", omics_n, "_SIG.xlsx")))
    }
  }
  
  res_plot <- draw_custom_cycle(results, PLOT_CONFIG, task)
  if(!is.null(res_plot)) {
    ggsave(file.path(out_dir, paste0("FigS2_Cycle_", t_name, ".pdf")), res_plot$plot, width=11, height=11, device=cairo_pdf)
    write.xlsx(res_plot$data %>% filter(Type!="Gap"), file.path(out_dir, paste0("SourceData_", t_name, ".xlsx")))
    draw_separate_legends(res_plot$min_p, res_plot$max_p, PLOT_CONFIG, out_dir, t_name)
    
    message(">>> Success! Outputs saved in: ", out_dir)
  } else {
    message(">>> Warning: No significant features found for ", t_name)
  }
}

message("\n########################################################")
message("Done. All tissue cycles generated at: ", CONF_OUT_ROOT)