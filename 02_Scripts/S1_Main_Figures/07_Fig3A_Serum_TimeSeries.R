# ==============================================================================
# Project: FinlandSports V2.0 - Full Pipeline
# Script:  07_Fig3ABC_Serum_TimeSeries.R
# Panels:  Fig 3A (Clusters), Fig 3B (Targets), Fig 3C (Double Ring Heatmap)
# Core:    Mfuzz Time-Series Clustering, Target Extraction, V2.0 Path Config.
# ==============================================================================

# [0] Global Environment Setup
setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

if (!require("pacman")) install.packages("pacman")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("Mfuzz")) BiocManager::install("Mfuzz")
if (!require("ComplexHeatmap")) BiocManager::install("ComplexHeatmap")

pacman::p_load(
  tidyverse, Mfuzz, ComplexHeatmap, circlize, 
  grid, scales, openxlsx
)

# ------------------------------------------------------------------------------
# 1. Unified Path Configuration (V2.0 Standard)
# ------------------------------------------------------------------------------
DIR_DATA <- "01_Clean_Data"
DIR_OUT  <- "03_Results/Fig_3/Fig3_Serum_TimeSeries"
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. Load Clinical and Sample Metadata
# ------------------------------------------------------------------------------
clin_df <- read_csv(file.path(DIR_DATA, "Clinical_Master_Strict.csv"), show_col_types = FALSE) %>%
  mutate(
    Role_Orig = case_when(Membercode == 1 ~ "Daughter", Membercode == 2 ~ "Mother", Membercode == 3 ~ "Father", TRUE ~ "Unknown"),
    Clinical_ID = tolower(paste0(FamilyID, "_", Role_Orig)),
    Fat_percent = as.numeric(`Fat(%)`),
    Obesity_Group = ifelse(Fat_percent >= 30, "Obese", "Lean")
  ) %>%
  drop_na(Obesity_Group) %>%
  select(Clinical_ID, Fat_percent, Obesity_Group)

time_levels <- c("fast", "pre", "post1h", "post3h")
nature_colors <- c("1" = "#E64B35FF", "2" = "#4DBBD5FF", "3" = "#00A087FF", "4" = "#3C5488FF")
group_colors <- c("Obese" = "#E64B35FF", "Lean" = "#4DBBD5FF")

theme_publication <- function() {
  theme_classic(base_size = 14) +
    theme(
      axis.line = element_line(linewidth = 0.8, color = "black"),
      axis.ticks = element_line(linewidth = 0.8, color = "black"),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", face = "bold", size = 14),
      strip.background = element_blank(), 
      strip.text = element_text(face = "bold", size = 14, hjust = 0.5), 
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      panel.spacing = unit(1.5, "lines"),
      panel.grid = element_blank() # Strictly remove background grid lines
    )
}

# Initialize Global List for Supplementary Table
SUPP_TABLE_LIST <<- list()

# ------------------------------------------------------------------------------
# 3. Core Grand Unified Execution Function
# ------------------------------------------------------------------------------
analyze_full_pipeline <- function(omics_name) {
  message(paste0("\n>>> Initiating Pipeline for: ", omics_name, " ..."))
  
  file_path <- file.path(DIR_DATA, paste0(omics_name, ".csv"))
  if(!file.exists(file_path)) {
    message("    [Skip] File not found: ", file_path)
    return(NULL)
  }
  
  # a. Data Loading & Baseline Alignment (V2.0 Column Parsing)
  omics_df <- read_csv(file_path, show_col_types = FALSE)
  colnames(omics_df)[1] <- "symbol"
  
  expr_mat <- as.matrix(omics_df[,-1])
  rownames(expr_mat) <- omics_df$symbol
  
  parsed_meta <- data.frame(Omics_ID = colnames(expr_mat), stringsAsFactors = FALSE) %>%
    mutate(
      Clean_Col = tolower(gsub("_twin\\.[0-9]+", "", Omics_ID)),
      TimePoint = case_when(
        grepl("pre", Clean_Col) ~ "pre", 
        grepl("fast", Clean_Col) ~ "fast",
        grepl("post1h", Clean_Col) ~ "post1h",
        grepl("post3h", Clean_Col) ~ "post3h",
        TRUE ~ "Other"
      ),
      Clinical_ID = gsub("_fast|_pre|_post1h|_post3h", "", Clean_Col)
    ) %>%
    inner_join(clin_df, by = "Clinical_ID") %>%
    filter(TimePoint %in% time_levels)
  
  long_df <- omics_df %>%
    pivot_longer(cols = -symbol, names_to = "Omics_ID", values_to = "Value") %>%
    inner_join(parsed_meta, by = "Omics_ID") %>%
    drop_na(Value)
  
  baseline_df <- long_df %>%
    filter(TimePoint == "fast") %>%
    group_by(symbol, Clinical_ID) %>%
    summarise(Base_Value = mean(Value, na.rm = TRUE), .groups = 'drop')
  
  fc_df <- long_df %>%
    inner_join(baseline_df, by = c("symbol", "Clinical_ID")) %>%
    mutate(Log2FC = Value - Base_Value) %>%
    filter(TimePoint != "fast")
  
  fast_fc <- long_df %>% filter(TimePoint == "fast") %>% mutate(Log2FC = 0, Base_Value = Value)
  
  full_fc_df <- bind_rows(fast_fc, fc_df) %>%
    mutate(TimePoint = factor(TimePoint, levels = time_levels))
  
  # b. Screen Top 25% Core Molecules
  sig_symbols <- full_fc_df %>%
    group_by(symbol) %>%
    summarise(variance = var(Log2FC, na.rm = TRUE), .groups = 'drop') %>%
    slice_max(order_by = variance, prop = 0.25, with_ties = FALSE) %>%
    pull(symbol)
  
  sig_fc_df <- full_fc_df %>% filter(symbol %in% sig_symbols)
  
  # c. Mfuzz Clustering
  time_mean_df <- sig_fc_df %>%
    group_by(symbol, TimePoint) %>%
    summarise(Mean_FC = mean(Log2FC, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = TimePoint, values_from = Mean_FC) %>%
    drop_na()
  
  mat_mfuzz <- as.matrix(time_mean_df[,-1])
  rownames(mat_mfuzz) <- time_mean_df$symbol
  
  eset <- new("ExpressionSet", exprs = mat_mfuzz)
  eset <- Mfuzz::standardise(eset)
  set.seed(123)
  cl <- Mfuzz::mfuzz(eset, c = 4, m = 1.25)
  
  cluster_info <- data.frame(
    symbol = rownames(mat_mfuzz),
    Cluster = as.character(cl$cluster), 
    Membership = apply(cl$membership, 1, max)
  )
  sig_fc_df <- sig_fc_df %>% left_join(cluster_info, by = "symbol") %>% drop_na(Cluster)
  
  # d. [Export] Save Statistical Tables
  stats_export <- sig_fc_df %>%
    group_by(symbol, Cluster, Membership, Obesity_Group, TimePoint) %>%
    summarise(Mean_Log2FC = mean(Log2FC, na.rm=TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = c(Obesity_Group, TimePoint), values_from = Mean_Log2FC) %>%
    arrange(Cluster, desc(Membership)) # Prioritize high-confidence core molecules
  
  # Clean sheet names for publication quality
  sheet_n <- gsub("_", " ", sub("^Cleaned_", "", omics_name))
  if(nchar(sheet_n) > 31) sheet_n <- substr(sheet_n, 1, 31) 
  
  # Add to global list
  SUPP_TABLE_LIST[[sheet_n]] <<- stats_export
  
  # =========================================================================
  # e. [Plot] Fig 3A: Gradient Core Trajectory (2x2 Layout)
  # =========================================================================
  plot_A_data <- sig_fc_df %>%
    group_by(Cluster, symbol, TimePoint, Membership) %>%
    summarise(Mean_FC = mean(Log2FC, na.rm=TRUE), .groups='drop') %>%
    mutate(Cluster_Label = paste("Cluster", Cluster))
  
  pA <- ggplot(plot_A_data, aes(x = as.numeric(TimePoint), y = Mean_FC)) +
    annotate("rect", xmin = 1, xmax = 2, ymin = -Inf, ymax = Inf, fill = "gray95", alpha = 0.8) +
    annotate("rect", xmin = 2, xmax = 3, ymin = -Inf, ymax = Inf, fill = "gray85", alpha = 0.6) +
    annotate("rect", xmin = 3, xmax = 4, ymin = -Inf, ymax = Inf, fill = "gray95", alpha = 0.8) +
    geom_line(aes(group = symbol, color = Cluster, alpha = Membership), linewidth = 0.6) +
    stat_summary(aes(group = Cluster, color = Cluster), fun = mean, geom = "line", linewidth = 2, alpha = 1) +
    scale_color_manual(values = nature_colors) +
    scale_alpha_continuous(range = c(0.05, 0.5), guide = "none") +
    facet_wrap(~ Cluster_Label, scales = "free_y", nrow = 2, ncol = 2) +
    scale_x_continuous(breaks = 1:4, labels = time_levels, expand = expansion(add = c(0.05, 0.2))) +
    scale_y_continuous(n.breaks = 5) + 
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none", 
      strip.background = element_blank(), 
      strip.text = element_text(face="bold"),
      aspect.ratio = 1, # Force square
      panel.grid = element_blank() # Strictly remove background grid lines
    ) +
    labs(x = "Time Point", y = "Log2 Fold Change (vs Fast)")
  
  ggsave(file.path(DIR_OUT, paste0(omics_name, "_Fig3A_Clusters_Gradient.pdf")), pA, width = 8, height = 8)
  
  # =========================================================================
  # f. [Plot] Fig 3B: New Curated Target Molecules (Deep Biological Significance)
  # =========================================================================
  if(grepl("Metabonomics", omics_name)) {
    target_mols <- c("bOHBut", "LA", "SM", "FALen")
  } else {
    target_mols <- c("PRG4", "LRRC8C", "PARP14", "CFD")
  }
  
  actual_targets <- intersect(target_mols, unique(sig_fc_df$symbol))
  
  if(length(actual_targets) > 0) {
    plot_B_data <- sig_fc_df %>%
      filter(symbol %in% actual_targets) %>%
      mutate(symbol = factor(symbol, levels = target_mols)) %>% 
      group_by(symbol, TimePoint, Obesity_Group) %>%
      summarise(
        Mean = mean(Log2FC, na.rm = TRUE),
        SEM = sd(Log2FC, na.rm = TRUE) / sqrt(n()),
        .groups = 'drop'
      )
    
    pB <- ggplot(plot_B_data, aes(x = TimePoint, y = Mean, color = Obesity_Group, fill = Obesity_Group, group = Obesity_Group)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) + 
      geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.15, linewidth = 0.8) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3, shape = 21, fill = "white", stroke = 1.2) + 
      scale_color_manual(values = group_colors) +
      scale_fill_manual(values = group_colors) +
      facet_wrap(~ symbol, scales = "free_y", nrow = 2, ncol = 4) + 
      theme_publication() +
      theme(aspect.ratio = 0.8) + 
      labs(x = "Time (Post-intervention)", y = "Mean Log2 Fold Change")
    
    ggsave(file.path(DIR_OUT, paste0(omics_name, "_Fig3B_Representative_Targets.pdf")), pB, width = 12, height = 7)
  }
  
  # =========================================================================
  # g. [Plot] Fig 3C: Ultimate Seamless Double Ring Heatmap
  # =========================================================================
  elite_mols <- cluster_info %>%
    group_by(Cluster) %>%
    slice_max(order_by = Membership, n = 12, with_ties = FALSE) %>%
    pull(symbol)
  
  heat_df <- sig_fc_df %>%
    filter(symbol %in% elite_mols) %>%
    group_by(symbol, Obesity_Group, TimePoint) %>%
    summarise(Mean_FC = mean(Log2FC, na.rm = TRUE), .groups = 'drop') %>%
    unite("Group_Time", Obesity_Group, TimePoint, sep = "_") %>%
    pivot_wider(names_from = Group_Time, values_from = Mean_FC) %>%
    left_join(cluster_info, by = "symbol")
  
  mat_heat <- as.matrix(heat_df %>% select(-symbol, -Cluster, -Membership))
  rownames(mat_heat) <- heat_df$symbol
  
  col_lean_rev <- intersect(paste0("Lean_", rev(time_levels)), colnames(mat_heat))
  col_obese_rev <- intersect(paste0("Obese_", rev(time_levels)), colnames(mat_heat))
  
  mat_combined <- cbind(mat_heat[, col_lean_rev], mat_heat[, col_obese_rev])
  
  col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#3C5488FF", "white", "#E64B35FF")) 
  row_split <- factor(heat_df$Cluster)
  
  pdf(file.path(DIR_OUT, paste0(omics_name, "_Fig3C_DoubleRing_Integrated.pdf")), width = 8, height = 8)
  
  circos.clear()
  circos.par(start.degree = 90, gap.degree = c(rep(2, nlevels(row_split)-1), 20))
  
  circos.heatmap(mat_combined, 
                 split = row_split, 
                 col = col_fun, 
                 cluster = TRUE,             
                 dend.side = "inside",       
                 dend.track.height = 0.15,   
                 track.height = 0.28,        
                 bg.border = "black", bg.lwd = 1.2,
                 rownames.side = "outside",  
                 rownames.cex = 0.8,         
                 rownames.font = 3)          
  
  circos.clear()
  
  lgd = Legend(title = "Log2FC\n(vs Fast)", col_fun = col_fun, direction = "vertical")
  draw(lgd, x = unit(10, "mm"), y = unit(10, "mm"), just = c("left", "bottom"))
  
  grid.text(
    "► Outer Ring (4 columns): Lean\n► Inner Ring (4 columns): Obese\n(fast → pre → post1h → post3h)", 
    x = unit(10, "mm"), y = unit(0.95, "npc"), 
    just = c("left", "top"),
    gp = gpar(fontsize = 12, fontface = "bold", col = "black")
  )
  
  dev.off()
  
  message(paste0(">>> Process completed for ", omics_name, " !"))
}

# ------------------------------------------------------------------------------
# 4. Execute Main Pipeline
# ------------------------------------------------------------------------------
# Run pipeline (Data will be automatically appended to the list and PDFs generated)
analyze_full_pipeline("Cleaned_Serum_Metabonomics")
analyze_full_pipeline("Cleaned_Serum_Proteomics")

# Compress the list into a clean multi-sheet Excel file
write.xlsx(SUPP_TABLE_LIST, file.path(DIR_OUT, "Data_TimeSeries_Clusters.xlsx"), overwrite = TRUE)
message("  -> Data (TimeSeries Clusters) ultimate supplementary table saved successfully!")

message("\n=======================================================")
message("ALL DONE! Check your Fig3-Serum folder for the final outputs.")
message("=======================================================")