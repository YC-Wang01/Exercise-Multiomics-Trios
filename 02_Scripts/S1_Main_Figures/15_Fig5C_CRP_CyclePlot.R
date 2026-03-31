# ==============================================================================
# Project: FinlandSports V2.0 - Fig 5C Standardized Circular Heatmap
# Target: Generate Fig5 Style Circular Plot for ANY Clinical Variable (Default CRP)
# Changes: Fixed Target Column Name (CRP), Strict Namespacing, English UI.
# ==============================================================================

# --- 0. Environment Initialization ---
# setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readr, dplyr, stringr, limma, tibble, 
  ggplot2, ggrepel, ggsci, ggpubr, grid, cowplot,
  scales, shadowtext, grDevices, openxlsx
)

# Set font (Register Arial on Windows)
if(Sys.info()[['sysname']] == "Windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

# ==============================================================================
# ### CONFIGURATION ###
# ==============================================================================

# --- A. Paths (V2.0 Standard) ---
DIR_IN  <- "01_Clean_Data"
DIR_OUT <- "03_Results/Fig_5/Fig5C_CRP_CyclePlot"

# --- B. Target Configuration ---
CONF_TARGET_COL   <- "CRP"           # [FIXED] Exact column name in Clinical_Master_Strict.csv
CONF_CENTER_LABEL <- "CRP"           # Label to display in the center of the plot
CONF_DO_LOG_TRANS <- TRUE            # Log transform CRP values

# --- C. Omics File Mapping (V2.0 CSVs) ---
CONF_OMICS_FILES <- list(
  Adipose = "Cleaned_Adipose_Proteomics.csv",
  Muscle  = "Cleaned_Muscle_Proteomics.csv",
  Serum   = "Cleaned_Serum_Proteomics.csv"
)

# --- D. Visual Colors (Nature NPG) ---
CONF_COLORS_TISSUE <- c(
  "Adipose" = "#E64B35", 
  "Muscle"  = "#00A087", 
  "Serum"   = "#4DBBD5"
)
# LogFC Gradient (Blue-White-Red)
CONF_COLORS_GRADIENT <- c("#3B4992", "#FFFFFF", "#EE0000") 
CONF_FONT_FAMILY     <- "Arial"

# Create output directory
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# PART 1: DATA LOADING & PREPARATION (V2.0 Standardized)
# ==============================================================================
message(">>> STEP 1: Loading Metadata & Target Variable [", CONF_TARGET_COL, "]...")

# 1. Load basic info and target from Strict Master Table
meta_merged <- read_csv(file.path(DIR_IN, "Clinical_Master_Strict.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(
    Role_Orig = dplyr::case_when(Membercode == 1 ~ "Daughter", Membercode == 2 ~ "Mother", Membercode == 3 ~ "Father", TRUE ~ "Unknown"),
    Subject_ID = tolower(paste0(FamilyID, "_", Role_Orig))
  ) %>%
  dplyr::rename(Target_Raw = dplyr::all_of(CONF_TARGET_COL)) %>%
  dplyr::filter(!is.na(Target_Raw)) %>%
  dplyr::mutate(
    Sex = factor(Gender), 
    Age = as.numeric(Age)
  )

# Dynamic Log Transformation
if(CONF_DO_LOG_TRANS) {
  meta_merged$Analysis_Var <- log(meta_merged$Target_Raw)
} else {
  meta_merged$Analysis_Var <- meta_merged$Target_Raw
}

# ==============================================================================
# PART 2: LIMMA ANALYSIS (Dynamic Model)
# ==============================================================================

run_limma_standardized <- function(tissue, filename) {
  fpath <- file.path(DIR_IN, filename)
  if(!file.exists(fpath)) return(NULL)
  
  raw <- read_csv(fpath, show_col_types = FALSE)
  feat_ids <- raw[[1]]; expr_raw <- as.matrix(raw[, -1])
  
  # Parse columns to match Subject_ID (Baseline only)
  clean_cols <- tolower(str_remove(colnames(expr_raw), "_twin\\.\\d+$")) %>% str_trim()
  col_meta <- data.frame(Matrix_Col = colnames(expr_raw), Clean_Col = clean_cols, stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      TimePoint = dplyr::case_when(grepl("pre", Clean_Col) ~ "pre", grepl("fast", Clean_Col) ~ "fast", TRUE ~ "Other"),
      Subject_ID = gsub("_fast|_pre|_post1h|_post3h", "", Clean_Col)
    ) %>% 
    dplyr::filter(TimePoint %in% c("pre", "fast")) %>%
    dplyr::inner_join(meta_merged, by = "Subject_ID") %>%
    dplyr::distinct(Subject_ID, .keep_all = TRUE)
  
  if(nrow(col_meta) < 10) return(NULL)
  
  mat <- expr_raw[, col_meta$Matrix_Col, drop=FALSE]
  rownames(mat) <- feat_ids
  
  # Filter & Impute
  mat <- mat[rowSums(is.na(mat)) < (ncol(mat)*0.5), ] 
  mat[is.na(mat)] <- min(mat, na.rm=T)/2
  if(max(mat, na.rm=T) > 50) mat <- log2(mat + 1)
  
  # --- Dynamic Formula ---
  design <- model.matrix(~ Analysis_Var + Age + Sex, data = col_meta)
  fit <- eBayes(lmFit(mat, design))
  
  # Extract Analysis_Var coefficients
  res <- topTable(fit, coef="Analysis_Var", number=Inf) %>% 
    tibble::rownames_to_column("Protein") %>% dplyr::mutate(Tissue = tissue)
  
  return(list(res = res))
}

message(">>> STEP 2: Running Limma on [", CONF_TARGET_COL, "]...")
results_store <- list()
for(t in names(CONF_OMICS_FILES)) { 
  results_store[[t]] <- run_limma_standardized(t, CONF_OMICS_FILES[[t]]) 
}

# ==============================================================================
# PART 3: PLOTTING FUNCTION (Standardized Visuals)
# ==============================================================================

draw_legend_only <- function(out_path) {
  # Gradient Legend
  df_dummy <- data.frame(x=1:100, y=1, z=seq(-2, 2, length.out=100))
  p1 <- ggplot(df_dummy, aes(x=x, y=y, fill=z)) +
    geom_tile() +
    scale_fill_gradient2(low=CONF_COLORS_GRADIENT[1], mid=CONF_COLORS_GRADIENT[2], high=CONF_COLORS_GRADIENT[3], 
                         midpoint=0, name="Log2 FC") +
    theme_void() + 
    theme(legend.position="right", 
          legend.title=element_text(face="bold", family=CONF_FONT_FAMILY),
          legend.text=element_text(family=CONF_FONT_FAMILY))
  leg1 <- cowplot::get_legend(p1)
  
  # Tissue Legend
  df_cat <- data.frame(Tissue=names(CONF_COLORS_TISSUE), v=1)
  p2 <- ggplot(df_cat, aes(x=Tissue, y=v, color=Tissue)) + 
    geom_line(size=2) +
    scale_color_manual(values=CONF_COLORS_TISSUE) +
    theme_void() + 
    theme(legend.position="right", 
          legend.title=element_text(face="bold", family=CONF_FONT_FAMILY),
          legend.text=element_text(family=CONF_FONT_FAMILY))
  leg2 <- cowplot::get_legend(p2)
  
  fin <- cowplot::plot_grid(leg2, leg1, ncol=1, align="v")
  ggsave(out_path, fin, width=2.5, height=4, device=cairo_pdf) 
  message("   -> Legend saved: ", out_path)
}

draw_standardized_cycle <- function(data_store, total_bars = 120) {
  
  message(">>> STEP 3: Designing Cycle Plot...")
  
  # --- A. Allocation ---
  get_cnt <- function(t) nrow(dplyr::filter(data_store[[t]]$res, P.Value < 0.05))
  cnts <- c(Adipose = get_cnt("Adipose"), Muscle = get_cnt("Muscle"), Serum = get_cnt("Serum"))
  if(sum(cnts)==0) cnts[] <- c(40,30,30)
  
  alloc <- round(total_bars * (cnts / sum(cnts)))
  if(alloc["Serum"] < 8) {
    alloc["Serum"] <- 8
    rem <- total_bars - 8
    alloc["Adipose"] <- round(rem * (cnts["Adipose"]/(cnts["Adipose"]+cnts["Muscle"])))
    alloc["Muscle"] <- rem - alloc["Adipose"]
  }
  
  t_order <- c("Adipose", "Muscle", "Serum")
  
  # --- B. Data Processing ---
  processed_list <- list()
  for(tissue_name in t_order) {
    if(tissue_name %in% names(data_store)) {
      df <- data_store[[tissue_name]]$res
      limit <- alloc[[tissue_name]]
      top_df <- df %>% 
        dplyr::arrange(P.Value) %>% 
        head(limit) %>%
        dplyr::mutate(
          Tissue = tissue_name,
          Height = -log10(P.Value),
          LogFC_Value = logFC
        )
      
      if(nrow(top_df) > 0) {
        top_df$Height <- ifelse(top_df$Height > 10, 10, top_df$Height)
        top_df <- top_df %>% dplyr::arrange(LogFC_Value) 
        processed_list[[tissue_name]] <- top_df
      }
    }
  }
  
  plot_df <- dplyr::bind_rows(processed_list)
  plot_df$id <- seq(1, nrow(plot_df))
  
  # --- C. Colors ---
  plot_df$Tissue_Color <- CONF_COLORS_TISSUE[plot_df$Tissue]
  
  max_fc <- max(abs(plot_df$LogFC_Value), na.rm=T)
  pal_func <- scales::col_numeric(
    palette = CONF_COLORS_GRADIENT, 
    domain = c(-max_fc, max_fc)
  )
  plot_df$fill_color <- pal_func(plot_df$LogFC_Value)
  
  # --- D. Radius & Labels ---
  r_inner_start <- 10
  r_inner_end   <- 13
  r_outer_start <- 14
  
  tissue_segments <- plot_df %>%
    dplyr::group_by(Tissue) %>%
    dplyr::summarize(start = min(id), end = max(id), mid = mean(c(min(id), max(id)))) %>%
    dplyr::ungroup()
  tissue_segments$col <- CONF_COLORS_TISSUE[as.character(tissue_segments$Tissue)]
  
  n_bars <- nrow(plot_df)
  angle <- 90 - 360 * (plot_df$id - 0.5) / n_bars
  plot_df$hjust <- ifelse(angle < -90, 1, 0)
  plot_df$angle <- ifelse(angle < -90, angle + 180, angle)
  
  # --- E. Plotting ---
  p <- ggplot() +
    
    # 1. Inner Ring
    geom_rect(data = plot_df,
              aes(xmin = id - 0.5, xmax = id + 0.5,
                  ymin = r_inner_start, ymax = r_inner_end,
                  fill = Tissue_Color, color = Tissue_Color), 
              size = 0.5) + 
    
    # 2. Outer Ring
    geom_rect(data = plot_df,
              aes(xmin = id - 0.5, xmax = id + 0.5,
                  ymin = r_outer_start, ymax = r_outer_start + Height,
                  fill = fill_color),
              color = "black", size = 0.05) + 
    
    # 3. Gene Labels
    geom_text(data = plot_df %>% dplyr::filter(!is.na(Protein)),
              aes(x = id, y = r_outer_start + Height + 0.5, 
                  label = Protein, angle = angle, hjust = hjust),
              color = "black", size = 2.2, fontface = "plain", family = CONF_FONT_FAMILY) +
    
    # 4. Tissue Labels
    geom_text(data = tissue_segments,
              aes(x = mid, y = 7, label = Tissue, color = col),
              size = 5, fontface = "bold", family = CONF_FONT_FAMILY) +
    
    # 5. Center Label
    annotate("text", x = 1, y = 0, label = CONF_CENTER_LABEL, 
             size = 7, fontface = "bold", family = CONF_FONT_FAMILY) +
    
    scale_fill_identity() +
    scale_color_identity() +
    scale_x_continuous(limits = c(0.5, n_bars + 0.5), expand = c(0, 0)) +
    ylim(0, r_outer_start + 12) + 
    coord_polar(theta = "x") +
    theme_void() +
    theme(plot.margin = unit(c(-1,-1,-1,-1), "cm"))
  
  return(p)
}

# --- 4. Execution ---
if(exists("results_store")) {
  p_final <- draw_standardized_cycle(results_store, total_bars = 120)
  
  out_name <- paste0("Fig5C_CyclePlot_", CONF_TARGET_COL, "_Nature.pdf")
  out_file <- file.path(DIR_OUT, out_name)
  ggsave(out_file, p_final, width = 10, height = 10, device = cairo_pdf)
  
  draw_legend_only(file.path(DIR_OUT, paste0("Fig5C_CyclePlot_Legend_", CONF_TARGET_COL, ".pdf")))
  
  message("\n>>> DONE! Cycle Plot Generation Complete.")
  message("    Output saved to: ", out_file)
} else {
  message("Error: results_store missing.")
}

# ==============================================================================
# 5. Generate Publication-Ready Supplementary Table (Data_S8)
# ==============================================================================
message("\n>>> STEP 4: Generating Publication-Ready Supplementary Table...")

wb_sup <- createWorkbook()
header_style <- createStyle(fontName = "Arial", fontSize = 11, fontColour = "white", fgFill = "#4F81BD", textDecoration = "bold")
sig_red      <- createStyle(fontColour = "#B2182B", textDecoration = "bold")
sig_blue     <- createStyle(fontColour = "#2166AC", textDecoration = "bold")

for(t_name in names(results_store)) {
  
  df_full <- results_store[[t_name]]$res %>%
    dplyr::mutate(
      Significance = dplyr::case_when(P.Value < 0.001 ~ "***", P.Value < 0.01 ~ "**", P.Value < 0.05 ~ "*", TRUE ~ "ns"),
      Protein = gsub("\\..*", "", Protein) 
    ) %>%
    dplyr::arrange(P.Value) %>% 
    dplyr::select(Tissue, Protein, Association_Coefficient = logFC, t_statistic = t, P_Value = P.Value, adj.P_Val_FDR = adj.P.Val, Significance)
  
  if(nrow(df_full) > 0) {
    sheet_n <- paste0("Limma_", t_name)
    addWorksheet(wb_sup, sheet_n)
    writeData(wb_sup, sheet_n, df_full)
    
    addStyle(wb_sup, sheet_n, style = header_style, rows = 1, cols = 1:ncol(df_full), gridExpand = TRUE)
    
    conditionalFormatting(wb_sup, sheet_n, cols = 3, rows = 2:(nrow(df_full)+1), rule = ">0", style = sig_red)
    conditionalFormatting(wb_sup, sheet_n, cols = 3, rows = 2:(nrow(df_full)+1), rule = "<0", style = sig_blue)
    conditionalFormatting(wb_sup, sheet_n, cols = 7, rows = 2:(nrow(df_full)+1), rule = "*", type = "contains", style = sig_red)
    
    setColWidths(wb_sup, sheet_n, cols = 1:ncol(df_full), widths = "auto")
  }
}

sup_out_name <- file.path(DIR_OUT, paste0("Data_S8_", CONF_TARGET_COL, "_Associations.xlsx"))
saveWorkbook(wb_sup, sup_out_name, overwrite = TRUE)
message("   -> Saved Full Supplementary Table to: ", sup_out_name)
message("========================================================================")
