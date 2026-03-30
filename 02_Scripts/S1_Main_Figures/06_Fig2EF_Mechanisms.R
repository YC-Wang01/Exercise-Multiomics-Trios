# ==============================================================================
# Project: FinlandSports V2.0 - Full Pipeline
# Target: Fig 2E & 2F Complete (Mechanisms & Heritability)
# Description: 
#   PART 0: Transcriptome-Proteome Correlation (Generates Driver Files!)
#   PART 1: Heritability Analysis (Uses Driver Files from Part 0)
#   PART 2: Hub Genes Visualization (Square Family Plots)
#   PART 3: Methylation Mechanism (Quadrant Plots)
# Features: V2.0 I/O Standardization, 100% Original Logic Replication.
# ==============================================================================

# [0] Global Environment Setup
setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readxl, dplyr, stringr, ggplot2, tidyr, readr, 
  pheatmap, grid, scales, openxlsx, ggpubr, 
  patchwork, ggrepel, cowplot, gridExtra, limma, tibble
)

# --- GLOBAL PATH CONFIGURATION (V2.0 Standard) ---

DIR_DATA_OMICS <- "01_Clean_Data"
DIR_ROOT_OUT   <- "03_Results/Fig_2/Fig2EF_Mechanisms"
DIR_TABLES     <- "03_Results/Fig_2/Tables"

invisible(lapply(c(DIR_ROOT_OUT, DIR_TABLES), dir.create, recursive = TRUE, showWarnings = FALSE))

message(">>> Path Configuration:")
message("   Data Dir:   ", DIR_DATA_OMICS)
message("   Output Dir: ", DIR_ROOT_OUT)
message("   Tables Dir: ", DIR_TABLES)

# ==============================================================================
# PART 0: Transcriptome-Proteome Correlation (Generating Drivers)
# ==============================================================================
message("\n>>> Starting PART 0: Generating Drivers (Transcriptome-Proteome Correlation)...")
SUPP_TABLE_LIST <<- list() 

# 0.1 Load Metadata (Docking with V2.0 Clinical Master)
FILE_META <- file.path(DIR_DATA_OMICS, "Clinical_Master_Strict.csv")
if(!file.exists(FILE_META)) stop("Metadata file missing at: ", FILE_META)

meta_data <- read_csv(FILE_META, show_col_types = FALSE) %>%
  mutate(
    Sample_Key = paste(FamilyID, Membercode, sep="_"),
    Fat_percent = `Fat(%)`, 
    Group = case_when(
      Fat_percent >= 30 ~ "Obese",
      Fat_percent < 30 ~ "Lean",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group))

# 0.2 Helper Functions for Part 0
parse_and_filter_columns <- function(col_names) {
  valid_indices <- c(); sample_groups <- c()
  for (i in seq_along(col_names)) {
    col <- col_names[i]
    clean_col <- str_remove(col, "\\.\\.\\.\\d+$") %>% str_trim()
    if (grepl("pre|fast", clean_col, ignore.case = TRUE) && !grepl("post", clean_col, ignore.case = TRUE)) {
      parts <- str_split(clean_col, "[_\\s]")[[1]]
      if(length(parts) >= 2) {
        fid <- parts[1]; role <- tolower(parts[2])
        mcode <- case_when(grepl("father", role) ~ "2", grepl("mother", role) ~ "3", TRUE ~ "1")
        s_key <- paste0(fid, "_", mcode)
        match_idx <- which(meta_data$Sample_Key == s_key)
        if (length(match_idx) == 1) {
          valid_indices <- c(valid_indices, i)
          sample_groups <- c(sample_groups, meta_data$Group[match_idx])
        }
      }
    }
  }
  return(list(indices = valid_indices, groups = sample_groups))
}

run_limma <- function(filename, omics_type) {
  full_path <- file.path(DIR_DATA_OMICS, filename) 
  if (!file.exists(full_path)) { message("   [Missing] ", filename, " in ", DIR_DATA_OMICS); return(NULL) }
  
  message(paste("   Running Limma on:", filename))
  if (grepl(".csv$", filename)) raw <- read.csv(full_path, check.names = FALSE) else raw <- read_excel(full_path, .name_repair = "unique")
  
  feat_ids <- raw[[1]]; expr_mat <- raw[, -1]
  parsed <- parse_and_filter_columns(colnames(expr_mat))
  
  if (length(parsed$indices) < 6) return(NULL)
  
  mat <- expr_mat[, parsed$indices]
  mat <- suppressWarnings(as.data.frame(lapply(mat, function(x) as.numeric(as.character(x)))))
  rownames(mat) <- feat_ids
  mat <- mat[rowSums(is.na(mat)) < (ncol(mat) * 0.5), ]; mat[is.na(mat)] <- min(mat, na.rm = TRUE) / 2
  if (max(mat, na.rm = TRUE) > 50) mat <- log2(mat + 1)
  
  groups <- factor(parsed$groups, levels = c("Lean", "Obese"))
  design <- model.matrix(~0 + groups); colnames(design) <- levels(groups)
  
  fit <- lmFit(mat, design)
  fit2 <- eBayes(contrasts.fit(fit, makeContrasts(Diff = Obese - Lean, levels = design)))
  
  res <- topTable(fit2, coef = "Diff", number = Inf) %>%
    rownames_to_column("Gene") %>% mutate(Omics = omics_type)
  res$Gene <- gsub(";.*", "", res$Gene)
  return(res)
}

plot_omics_correlation <- function(limma_trans, limma_prot, tissue_name) {
  df_trans <- limma_trans %>% dplyr::select(Gene, logFC, P.Value) %>% rename(logFC_Trans = logFC, P_Trans = P.Value)
  df_prot <- limma_prot %>% dplyr::select(Gene, logFC, P.Value) %>% rename(logFC_Prot = logFC, P_Prot = P.Value)
  merged_df <- inner_join(df_trans, df_prot, by = "Gene")
  
  if(nrow(merged_df) < 10) return(NULL)
  
  P_CUTOFF <- 0.05
  merged_df$Category <- "Not Sig"
  merged_df$Category[merged_df$P_Trans < P_CUTOFF & merged_df$logFC_Trans > 0 & merged_df$P_Prot < P_CUTOFF & merged_df$logFC_Prot > 0] <- "Both Up"
  merged_df$Category[merged_df$P_Trans < P_CUTOFF & merged_df$logFC_Trans < 0 & merged_df$P_Prot < P_CUTOFF & merged_df$logFC_Prot < 0] <- "Both Down"
  merged_df$Category[merged_df$P_Trans < P_CUTOFF & merged_df$P_Prot < P_CUTOFF & sign(merged_df$logFC_Trans) != sign(merged_df$logFC_Prot)] <- "Discordant"
  merged_df$Category[merged_df$Category == "Not Sig" & merged_df$P_Trans < P_CUTOFF] <- "Trans Only"
  merged_df$Category[merged_df$Category == "Not Sig" & merged_df$P_Prot < P_CUTOFF] <- "Prot Only"
  
  # Save Plot
  my_colors <- c("Not Sig"="grey90", "Trans Only"="#F39B7F", "Prot Only"="#8491B4", "Discordant"="#00A087", "Both Up"="#DC0000", "Both Down"="#3C5488")
  cor_res <- cor.test(merged_df$logFC_Trans, merged_df$logFC_Prot, method = "pearson")
  
  p <- ggplot(merged_df, aes(x = logFC_Trans, y = logFC_Prot)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") + geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = Category), alpha = 0.8) +
    annotate("text", x = min(merged_df$logFC_Trans), y = max(merged_df$logFC_Prot), label = paste0("R = ", round(cor_res$estimate, 3)), hjust = 0, vjust = 1) +
    scale_color_manual(values = my_colors) + theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank()) +
    labs(title = paste0(tissue_name, ": Transcriptome vs Proteomics"))
  
  ggsave(file.path(DIR_ROOT_OUT, paste0("Fig2E_Correlation_", tissue_name, ".pdf")), p, width = 8, height = 7)
  
  # SAVE DRIVER FILE
  consistent_genes <- merged_df %>% 
    filter(Category %in% c("Both Up", "Both Down", "Discordant")) %>%
    arrange(desc(Category), desc(abs(logFC_Prot)))
  
  driver_filename <- paste0("Fig2E_Correlation_", tissue_name, "_Drivers.xlsx")
  write.xlsx(consistent_genes, file.path(DIR_ROOT_OUT, driver_filename))
  SUPP_TABLE_LIST[[paste0("Drivers_", tissue_name)]] <<- consistent_genes
  message(paste("   >>> Driver File Generated:", driver_filename))
}

# 0.3 Execute Part
adi_trans <- run_limma("Cleaned_Adipose_Microarray.csv", "Transcriptome")
adi_prot  <- run_limma("Cleaned_Adipose_Proteomics.csv", "Proteomics")
if(!is.null(adi_trans) & !is.null(adi_prot)) plot_omics_correlation(adi_trans, adi_prot, "Adipose")

mus_trans <- run_limma("Cleaned_Muscle_Microarray.csv", "Transcriptome")
mus_prot  <- run_limma("Cleaned_Muscle_Proteomics.csv", "Proteomics")
if(!is.null(mus_trans) & !is.null(mus_prot)) plot_omics_correlation(mus_trans, mus_prot, "Muscle")


# ==============================================================================
# PART 1: Heritability Analysis (Microarray & Proteomics)
# ==============================================================================
message("\n>>> Starting PART 1: Heritability Analysis (Stats, Heatmaps, Violins)...")
USER_CONFIG_P1 <- list(
  list(Tissue = "Adipose", Type = "Microarray", File = "Cleaned_Adipose_Microarray.csv",  Driver = "Fig2E_Correlation_Adipose_Drivers.xlsx"),
  list(Tissue = "Adipose", Type = "Proteomics", File = "Cleaned_Adipose_Proteomics.csv",  Driver = "Fig2E_Correlation_Adipose_Drivers.xlsx"),
  list(Tissue = "Muscle",  Type = "Microarray", File = "Cleaned_Muscle_Microarray.csv",   Driver = "Fig2E_Correlation_Muscle_Drivers.xlsx"),
  list(Tissue = "Muscle",  Type = "Proteomics", File = "Cleaned_Muscle_Proteomics.csv",   Driver = "Fig2E_Correlation_Muscle_Drivers.xlsx")
)

OBESITY_SIG_GENES <- c(
  "SPARC", "MAP3K5", "MALSU1", "CLIP1", 
  "HTRA1", "POLR2I", "CRYAB", "FTH1", "HSPB7", "LRIG1", "SLC14A1", "SEC14L4",
  "ACLY", "BCAT2", "LDHB", "ITGAV", "EGFL6", "UCHL1", "CES1", "CA3",
  "GSTO1", "RRN3", "MYBPH", "NDEL1", "RING1", "GPX3", "EZR", "CALM1",
  "MYH1", "MSTN", "NNMT", "ACBD5", "DHRS7", "DYNLL1", "HABP2", "APRT"
)

load_data_part1 <- function(filename) {
  full_path <- file.path(DIR_DATA_OMICS, filename)
  if(!file.exists(full_path)) { message("   [Missing Data] ", filename); return(NULL) }
  
  if (grepl(".csv$", filename)) raw <- read.csv(full_path, check.names = FALSE) else raw <- read_excel(full_path, .name_repair = "unique")
  raw <- as.data.frame(raw, check.names = FALSE) 
  genes <- gsub(";.*", "", as.character(raw[[1]]))
  expr <- suppressWarnings(as.data.frame(lapply(raw[, -1], as.numeric), check.names = FALSE))
  if(any(duplicated(genes))) { keep <- !duplicated(genes); genes <- genes[keep]; expr <- expr[keep, ] }
  rownames(expr) <- genes
  if(max(expr, na.rm = TRUE) > 100) expr <- log2(expr + 1)
  
  cols <- colnames(expr)
  pre_cols <- cols[grepl("pre|fast", cols, ignore.case = TRUE) & !grepl("post", cols, ignore.case = TRUE)]
  sample_map <- data.frame(ColName = pre_cols, stringsAsFactors = FALSE)
  sample_map$FamilyID <- str_extract(sample_map$ColName, "\\d+") 
  sample_map$Role[grepl("daughter|son|child", sample_map$ColName, ignore.case = TRUE)] <- "Offspring"
  sample_map$Role[grepl("mother", sample_map$ColName, ignore.case = TRUE)] <- "Mother"
  sample_map$Role[grepl("father", sample_map$ColName, ignore.case = TRUE)] <- "Father"
  
  trio_map <- sample_map %>% filter(!is.na(Role)) %>% 
    pivot_wider(names_from = Role, values_from = ColName, values_fn = function(x) x[1]) %>%
    filter(!is.na(Offspring)) %>% as.data.frame()
  return(list(expr = expr, map = trio_map))
}

calc_correlations_p1 <- function(gene_list, data_obj) {
  expr_mat <- data_obj$expr; trio_map <- data_obj$map; res_list <- list()
  valid_genes <- intersect(gene_list, rownames(expr_mat))
  
  for(g in valid_genes) {
    vals <- expr_mat[g, , drop = FALSE]
    o_vals <- unname(as.numeric(vals[1, as.character(trio_map$Offspring)]))
    m_vals <- unname(sapply(as.character(trio_map$Mother), function(x) if(is.na(x)) NA else as.numeric(vals[1, x])))
    f_vals <- unname(sapply(as.character(trio_map$Father), function(x) if(is.na(x)) NA else as.numeric(vals[1, x])))
    
    df_m <- data.frame(O = o_vals, P = m_vals) %>% na.omit()
    df_f <- data.frame(O = o_vals, P = f_vals) %>% na.omit()
    
    cor_m <- if(nrow(df_m) >= 3) cor.test(df_m$O, df_m$P, method = "pearson")$estimate else 0
    p_m <- if(nrow(df_m) >= 3) cor.test(df_m$O, df_m$P, method = "pearson")$p.value else NA
    cor_f <- if(nrow(df_f) >= 3) cor.test(df_f$O, df_f$P, method = "pearson")$estimate else 0
    p_f <- if(nrow(df_f) >= 3) cor.test(df_f$O, df_f$P, method = "pearson")$p.value else NA
    
    res_list[[g]] <- data.frame(Gene = g, Cor_Mother = as.numeric(cor_m), Pval_Mother = as.numeric(p_m),
                                Cor_Father = as.numeric(cor_f), Pval_Father = as.numeric(p_f))
  }
  if(length(res_list) == 0) return(NULL)
  return(bind_rows(res_list))
}

process_task_p1 <- function(task) {
  message(paste("   Processing:", task$Tissue, task$Type))
  prefix <- paste0(task$Tissue, "_", task$Type)
  data_obj <- load_data_part1(task$File)
  if(is.null(data_obj)) return()
  
  driver_path <- file.path(DIR_ROOT_OUT, task$Driver)
  if(!file.exists(driver_path)) { message("      [Warning] Driver file not found: ", driver_path); return() }
  df_drivers <- read.xlsx(driver_path)
  
  res_stats <- calc_correlations_p1(df_drivers$Gene, data_obj)
  if(is.null(res_stats)) return()
  
  res_stats$Bias <- res_stats$Cor_Mother - res_stats$Cor_Father
  SUPP_TABLE_LIST[[paste0("Herit_", prefix)]] <<- res_stats
  
  sig_in_data <- res_stats %>% filter(Gene %in% OBESITY_SIG_GENES)
  if(nrow(sig_in_data) > 0) {
    # Heatmap
    plot_mat <- sig_in_data %>% arrange(desc(Bias)) %>% dplyr::select(Gene, Cor_Mother, Cor_Father)
    rownames(plot_mat) <- plot_mat$Gene; plot_mat$Gene <- NULL
    pdf(file.path(DIR_ROOT_OUT, paste0("Fig2F_Heatmap_ObesitySigGenes_", prefix, ".pdf")), width = 3, height = max(3, nrow(plot_mat)*0.25))
    pheatmap(plot_mat, cluster_rows=F, cluster_cols=F, display_numbers=T, number_color="black",
             color=colorRampPalette(c("#3C5488", "white", "#DC0000"))(100), breaks=seq(-1, 1, length.out=101),
             cellwidth=30, cellheight=12, angle_col="45", legend=F)
    dev.off()
    
    # Scatter
    p_scatter <- ggplot(res_stats, aes(x = Cor_Father, y = Cor_Mother)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
      geom_point(aes(color = Bias), alpha = 0.4, size = 1.5) +
      geom_point(data = sig_in_data, aes(color = Bias), size = 3, shape = 1, stroke = 1, color = "black") +
      geom_point(data = sig_in_data, aes(color = Bias), size = 3, alpha = 1) +
      scale_color_gradient2(low = "#3C5488", mid = "grey95", high = "#DC0000", midpoint = 0) +
      geom_text_repel(data = sig_in_data, aes(label = Gene), size = 4, fontface = "bold") +
      theme_bw() + theme(aspect.ratio = 1, legend.position = "none", panel.grid = element_blank()) +
      coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1))
    ggsave(file.path(DIR_ROOT_OUT, paste0("Fig2F_Scatter_ObesitySigGenes_", prefix, ".pdf")), p_scatter, width = 5, height = 5)
  }
}

# Generate Legend for Heatmap
grad_df <- data.frame(y = seq(-1, 1, length.out = 200))
p_heat_legend <- ggplot(grad_df, aes(x = 1, y = y, fill = y)) +
  geom_tile() + scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#DC0000", midpoint = 0) +
  scale_y_continuous(breaks = c(-1, 0, 1), labels = c("-1", "0", "1"), expand = c(0,0)) +
  theme_classic() + theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank())
ggsave(file.path(DIR_ROOT_OUT, "Fig2F_Legend_Heatmap.pdf"), p_heat_legend, width = 1.2, height = 3)

for(task in USER_CONFIG_P1) { tryCatch({ process_task_p1(task) }, error = function(e) { message("Error in Part 1 task: ", e) }) }

# ==============================================================================
# PART 2: Hub Genes Visualization (Square Family Plots)
# ==============================================================================
message("\n>>> Starting PART 2: Hub Genes Visualization (Square Layout)...")

FILE_RAW_HUB <- file.path(DIR_DATA_OMICS, "Cleaned_Adipose_Proteomics.csv")

VIS_CONFIG <- list(
  colors = c("Father" = "#4DBBD5", "Mother" = "#E64B35", "Daughter" = "#00A087"),
  shapes = c("Father" = 15, "Mother" = 16, "Daughter" = 17),
  point_size = 4.5, line_width = 1.2, text_color = "black" 
)
TARGET_GENES_HUB <- c("SPARC", "MAP3K5", "MALSU1", "CLIP1")

load_data_hub <- function(filepath) {
  if(!file.exists(filepath)) { message("   [Missing] ", filepath); return(NULL) }
  if (grepl(".csv$", filepath)) raw <- read.csv(filepath, check.names = FALSE) else raw <- read_excel(filepath, .name_repair = "unique")
  raw <- as.data.frame(raw, check.names = FALSE) 
  genes <- gsub(";.*", "", as.character(raw[[1]]))
  expr <- suppressWarnings(as.data.frame(lapply(raw[, -1], as.numeric), check.names = FALSE))
  if(any(duplicated(genes))) { keep <- !duplicated(genes); genes <- genes[keep]; expr <- expr[keep, ] }
  rownames(expr) <- genes
  if(max(expr, na.rm = TRUE) > 100) expr <- log2(expr + 1)
  
  cols <- colnames(expr)
  pre_cols <- cols[grepl("pre|fast", cols, ignore.case = TRUE) & !grepl("post", cols, ignore.case = TRUE)]
  sample_map <- data.frame(ColName = pre_cols, stringsAsFactors = FALSE)
  sample_map$FamilyID <- str_extract(sample_map$ColName, "\\d+") 
  sample_map$Role[grepl("daughter|son|child", sample_map$ColName, ignore.case = TRUE)] <- "Daughter"
  sample_map$Role[grepl("mother", sample_map$ColName, ignore.case = TRUE)] <- "Mother"
  sample_map$Role[grepl("father", sample_map$ColName, ignore.case = TRUE)] <- "Father"
  
  trio_map <- sample_map %>% filter(!is.na(Role)) %>% 
    pivot_wider(names_from = Role, values_from = ColName, values_fn = function(x) x[1]) %>%
    as.data.frame()
  return(list(expr = expr, map = trio_map))
}

run_hub_vis <- function() {
  data_obj <- load_data_hub(FILE_RAW_HUB)
  if(is.null(data_obj)) return()
  
  expr_mat <- data_obj$expr; trio_map <- data_obj$map
  all_fams <- sort(as.numeric(unique(trio_map$FamilyID)))
  family_map_df <- data.frame(OriginalID = as.character(all_fams), NewID = as.character(seq_along(all_fams)), stringsAsFactors = FALSE)
  
  for(g in TARGET_GENES_HUB) {
    if(!g %in% rownames(expr_mat)) { message("   [Warning] Hub Gene ", g, " not found."); next }
    vals <- expr_mat[g, , drop = FALSE]
    
    df_wide <- data.frame(FamilyID = character(), Father = numeric(), Mother = numeric(), Daughter = numeric(), stringsAsFactors = F)
    for(i in 1:nrow(trio_map)) {
      fid <- trio_map$FamilyID[i]
      v_m <- if(!is.null(trio_map$Mother) && !is.na(trio_map$Mother[i])) vals[1, trio_map$Mother[i]] else NA
      v_f <- if(!is.null(trio_map$Father) && !is.na(trio_map$Father[i])) vals[1, trio_map$Father[i]] else NA
      v_d <- if(!is.null(trio_map$Daughter) && !is.na(trio_map$Daughter[i])) vals[1, trio_map$Daughter[i]] else NA
      if(all(is.na(c(v_m, v_f, v_d)))) next 
      df_wide <- rbind(df_wide, data.frame(FamilyID = as.character(fid), Father = v_f, Mother = v_m, Daughter = v_d))
    }
    if(nrow(df_wide) == 0) next
    
    df_wide$MemberCount <- rowSums(!is.na(df_wide[, c("Father", "Mother", "Daughter")]))
    df_wide <- df_wide %>% filter(MemberCount >= 2)
    if(nrow(df_wide) == 0) next
    
    df_wide <- df_wide %>% left_join(family_map_df, by = c("FamilyID" = "OriginalID")) %>% mutate(DisplayID = NewID) %>% dplyr::select(-NewID)
    
    # Stats
    df_md <- df_wide %>% filter(!is.na(Mother) & !is.na(Daughter))
    if(nrow(df_md) > 2) {
      ct_m <- cor.test(df_md$Mother, df_md$Daughter)
      p_str <- ifelse(ct_m$p.value < 0.001, "< 0.001", paste0("= ", signif(ct_m$p.value, 2)))
      stats_m <- paste0("Mother-Daughter: r = ", round(ct_m$estimate, 2), ", P ", p_str)
    } else { stats_m <- "Mother-Daughter: N too small" }
    
    df_fd <- df_wide %>% filter(!is.na(Father) & !is.na(Daughter))
    if(nrow(df_fd) > 2) {
      ct_f <- cor.test(df_fd$Father, df_fd$Daughter)
      p_str <- ifelse(ct_f$p.value < 0.001, "< 0.001", paste0("= ", signif(ct_f$p.value, 2)))
      stats_f <- paste0("Father-Daughter: r = ", round(ct_f$estimate, 2), ", P ", p_str)
    } else { stats_f <- "Father-Daughter: N too small" }
    
    df_wide <- df_wide %>% 
      mutate(GroupWeight = case_when(!is.na(Father) & !is.na(Mother) & !is.na(Daughter) ~ 3, is.na(Father) & !is.na(Mother) & !is.na(Daughter) ~ 2, TRUE ~ 1),
             SortKey = ifelse(!is.na(Mother), Mother, Father)) %>% 
      arrange(GroupWeight, SortKey) %>% mutate(DisplayID = factor(DisplayID, levels = unique(DisplayID)))
    
    df_lines_md <- df_wide %>% filter(GroupWeight == 2)
    df_long <- df_wide %>% pivot_longer(cols = c("Father", "Mother", "Daughter"), names_to = "Role", values_to = "Expression") %>% na.omit()
    
    p <- ggplot() +
      geom_segment(data = subset(df_wide, !is.na(Father) & !is.na(Mother)), aes(x = Father, xend = Mother, y = DisplayID, yend = DisplayID), color = "grey85", size = VIS_CONFIG$line_width) +
      geom_segment(data = df_lines_md, aes(x = Mother, xend = Daughter, y = DisplayID, yend = DisplayID), color = VIS_CONFIG$colors["Mother"], alpha = 0.3, size = 0.5, linetype = "dotted") +
      geom_point(data = df_long, aes(x = Expression, y = DisplayID, color = Role, shape = Role), size = VIS_CONFIG$point_size, alpha = 0.9) +
      scale_color_manual(values = VIS_CONFIG$colors) + scale_shape_manual(values = VIS_CONFIG$shapes) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
      labs(title = paste0(g, " (Adipose Proteomics)"), subtitle = paste(stats_m, "\n", stats_f), x = "Log2(Protein expression)", y = "Family") +
      theme_bw() +
      theme(aspect.ratio = 1, plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = VIS_CONFIG$text_color),
            plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic", color = VIS_CONFIG$text_color),
            axis.title = element_text(face = "bold", color = VIS_CONFIG$text_color, size = 12),
            axis.text = element_text(color = VIS_CONFIG$text_color, size = 10),
            axis.ticks = element_line(color = VIS_CONFIG$text_color),
            panel.grid = element_blank(), 
            panel.border = element_rect(color = VIS_CONFIG$text_color, fill = NA, size = 1),
            legend.position = "top", legend.title = element_blank(), legend.text = element_text(color = VIS_CONFIG$text_color))
    
    ggsave(file.path(DIR_ROOT_OUT, paste0("Fig2F_Square_", g, ".pdf")), p, width = 6, height = 6)
  }
}
tryCatch({ run_hub_vis() }, error = function(e) { message("Error in Part 2: ", e) })


# ==============================================================================
# PART 3: Methylation Mechanism (Quadrant Plots)
# ==============================================================================
message("\n>>> Starting PART 3: Methylation Mechanism Analysis...")
TARGETS_MET_ADI <- c("SPARC", "MAP3K5", "MALSU1", "CLIP1", "HTRA1", "POLR2I", "CRYAB", "FTH1", "HSPB7", "LRIG1", "SLC14A1", "SEC14L4", "ACLY", "BCAT2", "LDHB", "ITGAV", "EGFL6", "UCHL1", "CES1", "CA3")
TARGETS_MET_MUS <- c("SPARC", "MAP3K5", "MALSU1", "CLIP1", "GSTO1", "RRN3", "MYBPH", "NDEL1", "RING1", "GPX3", "EZR", "CALM1", "MYH1", "MSTN", "NNMT", "ACBD5", "DHRS7", "DYNLL1", "HABP2", "APRT")

load_data_met <- function(filename) {
  full_path <- file.path(DIR_DATA_OMICS, filename)
  if(!file.exists(full_path)) return(NULL)
  if (grepl(".csv$", filename)) raw <- read.csv(full_path, check.names = FALSE) else raw <- read_excel(full_path, .name_repair = "unique")
  raw <- as.data.frame(raw, check.names = FALSE)
  ids <- gsub(";.*", "", as.character(raw[[1]]))
  expr <- suppressWarnings(as.data.frame(lapply(raw[, -1], as.numeric), check.names = FALSE))
  if(any(duplicated(ids))) { expr <- expr[!duplicated(ids), ]; ids <- ids[!duplicated(ids)] }
  rownames(expr) <- ids
  return(expr)
}

get_trio_map_met <- function(colnames_vec) {
  pre_cols <- colnames_vec[grepl("pre|fast", colnames_vec, ignore.case = TRUE) & !grepl("post", colnames_vec, ignore.case = TRUE)]
  map <- data.frame(ColName = pre_cols, stringsAsFactors = FALSE)
  map$FamilyID <- str_extract(map$ColName, "\\d+")
  map$Role <- NA
  map$Role[grepl("daughter|son|child", map$ColName, ignore.case = TRUE)] <- "Offspring"
  map$Role[grepl("mother", map$ColName, ignore.case = TRUE)] <- "Mother"
  map$Role[grepl("father", map$ColName, ignore.case = TRUE)] <- "Father"
  map %>% filter(!is.na(Role)) %>% pivot_wider(names_from = Role, values_from = ColName, values_fn = function(x) x[1]) %>% filter(!is.na(Offspring)) %>% as.data.frame()
}

analyze_methylation <- function(tissue, gene_list, meth_file, expr_file) {
  message(paste("   Analyzing Methylation for:", tissue))
  dat_meth <- load_data_met(meth_file); dat_expr <- load_data_met(expr_file)
  if(is.null(dat_meth) || is.null(dat_expr)) return(NULL)
  
  map_meth <- get_trio_map_met(colnames(dat_meth)); map_expr <- get_trio_map_met(colnames(dat_expr))
  results_list <- list()
  
  for(g in gene_list) {
    probe_idx <- which(rownames(dat_meth) == g)
    if(length(probe_idx) == 0) probe_idx <- grep(paste0("^", g), rownames(dat_meth))
    if(length(probe_idx) == 0) next
    
    for(i in probe_idx) {
      probe_id <- rownames(dat_meth)[i]
      common_fams <- intersect(map_meth$FamilyID, map_expr$FamilyID)
      if(length(common_fams) < 5) next
      
      vec_m_d <- c(); vec_m_m <- c(); vec_m_f <- c(); vec_e_d <- c()
      for(fid in common_fams) {
        col_m_d <- map_meth$Offspring[map_meth$FamilyID == fid]
        col_m_m <- map_meth$Mother[map_meth$FamilyID == fid]
        col_m_f <- map_meth$Father[map_meth$FamilyID == fid]
        col_e_d <- map_expr$Offspring[map_expr$FamilyID == fid]
        
        if(length(col_m_d)>0 && length(col_e_d)>0) {
          val_m_d <- as.numeric(dat_meth[probe_id, col_m_d])
          val_m_m <- if(length(col_m_m)>0) as.numeric(dat_meth[probe_id, col_m_m]) else NA
          val_m_f <- if(length(col_m_f)>0) as.numeric(dat_meth[probe_id, col_m_f]) else NA
          val_e_d <- if(g %in% rownames(dat_expr)) as.numeric(dat_expr[g, col_e_d]) else NA
          vec_m_d <- c(vec_m_d, val_m_d); vec_m_m <- c(vec_m_m, val_m_m)
          vec_m_f <- c(vec_m_f, val_m_f); vec_e_d <- c(vec_e_d, val_e_d)
        }
      }
      cor_h_m <- cor(vec_m_d, vec_m_m, use = "pairwise.complete.obs", method = "pearson")
      cor_h_f <- cor(vec_m_d, vec_m_f, use = "pairwise.complete.obs", method = "pearson")
      cor_mech <- cor(vec_m_d, vec_e_d, use = "pairwise.complete.obs", method = "pearson")
      if(!is.na(cor_h_m) && !is.na(cor_h_f)) {
        results_list[[paste0(probe_id)]] <- data.frame(Gene=g, Probe=probe_id, Meth_Herit_Mother=cor_h_m, Meth_Herit_Father=cor_h_f, Bias=cor_h_m - cor_h_f, Meth_Expr_Cor=cor_mech)
      }
    }
  }
  if(length(results_list) > 0) return(bind_rows(results_list)) else return(NULL)
}

plot_mechanism_v2 <- function(df, tissue) {
  if(is.null(df)) return()
  df$Label <- ifelse((df$Meth_Expr_Cor < -0.15 & abs(df$Bias) > 0.15) | abs(df$Bias) > 0.4, df$Gene, NA)
  
  p <- ggplot(df, aes(x = Meth_Expr_Cor, y = Bias)) +
    annotate("rect", xmin = -1, xmax = 0, ymin = 0, ymax = 1, fill = "#DC0000", alpha = 0.08) + 
    annotate("rect", xmin = -1, xmax = 0, ymin = -1, ymax = 0, fill = "#3C5488", alpha = 0.08) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(aes(size = abs(Bias), color = Bias), alpha = 0.85, stroke=0.2) +
    scale_color_gradient2(low = "#3C5488", mid = "grey80", high = "#DC0000", limits=c(-1,1)) +
    scale_size_continuous(range = c(2, 6), limits=c(0, 1)) +
    geom_text_repel(aes(label = Label), size = 4, fontface = "bold", box.padding = 0.6, point.padding = 0.3, max.overlaps = 50, segment.color = "grey30", segment.size = 0.3) +
    theme_bw() + theme(aspect.ratio = 1, panel.grid = element_blank(), legend.position = "none", plot.title = element_blank(), axis.title = element_text(face = "bold", size = 11), axis.text = element_text(color = "black")) +
    coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
    labs(x = "Mechanism: Correlation (Methylation vs Expression)", y = "Origin: Heritability Bias (Mother - Father)")
  ggsave(file.path(DIR_ROOT_OUT, paste0("Fig2F_Mechanism_Plot_", tissue, ".pdf")), p, width = 5, height = 5)
}

# Generate Legend
dummy <- data.frame(x=1:5, y=1:5, Bias=seq(-1, 1, length.out=5), Size=seq(0, 1, length.out=5))
p_dummy <- ggplot(dummy, aes(x, y, color=Bias, size=Size)) +
  geom_point() +
  scale_color_gradient2(low = "#3C5488", mid = "grey80", high = "#DC0000", name = "Heritability Bias\n(Mother - Father)") +
  scale_size_continuous(name = "|Bias Magnitude|") + theme_bw() + theme(legend.position = "right")
leg <- get_legend(p_dummy)
pdf(file.path(DIR_ROOT_OUT, "Fig2F_Legend_Mechanism.pdf"), width = 3, height = 4); grid.newpage(); grid.draw(leg); dev.off()

res_adi <- analyze_methylation("Adipose", TARGETS_MET_ADI, "Cleaned_Adipose_Methylation.csv", "Cleaned_Adipose_Microarray.csv")
if(!is.null(res_adi)) { 
  SUPP_TABLE_LIST[["Meth_Mech_Adipose"]] <- res_adi 
  plot_mechanism_v2(res_adi, "Adipose") 
}

res_mus <- analyze_methylation("Muscle", TARGETS_MET_MUS, "Cleaned_Muscle_Methylation.csv", "Cleaned_Muscle_Microarray.csv")
if(!is.null(res_mus)) { 
  SUPP_TABLE_LIST[["Meth_Mech_Muscle"]] <- res_mus
  plot_mechanism_v2(res_mus, "Muscle") 
}

# ==============================================================================
# [Final Output] Save supplementary data tables
# ==============================================================================
write.xlsx(SUPP_TABLE_LIST, file.path(DIR_TABLES, "Data_Molecular_Mechanisms.xlsx"), overwrite = TRUE)
message("  -> Data (Molecular Mechanisms) summary table saved successfully!")
message("\n========================================================")
message(">>> ALL TASKS COMPLETED SUCCESSFULLY.")
message(">>> Results saved in: ", DIR_ROOT_OUT)