# ==============================================================================
# Project: FinlandSports V2.0 
# Script:  03_Fig1H_Clinical_Phenotypes.R
# Description: Clinical phenotypes dual-table merge and boxplot formatting (Fig 1H)
# Features: Enforced 1:1 square ratio, smart decimal scale engine, Fat(%) assertion, strict type alignment
# ==============================================================================

# [0. Lock global working directory]
# setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, readxl, dplyr, stringr, ggplot2, ggpubr, gridExtra, scales)

# ==============================================================================
# 1. Environment and relative path setup
# ==============================================================================
IN_CLINICAL <- "01_Clean_Data/Clinical_Master_Strict.csv"
VOL_FILE    <- "00_Raw_Data/VolunteerData.xlsx"  
OUT_DIR     <- "03_Results/Fig_1/Panel_H_Clinical"

if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

pheno_colors <- c("Obese" = "#F79647", "Lean" = "#3DA6AE")

theme_clean <- theme_bw(base_size = 15) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color = "black", face = "bold", size = 16),
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", color = "black", size = 13), 
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title.y = element_text(face = "bold", color = "black", margin = margin(r = 15)),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(), 
    legend.position = "none"
  )

# ==============================================================================
# 2. Smart join and Fat(%) grouping (with strict type alignment engine)
# ==============================================================================
cat("\n>>> [Step 1] Loading and Joining Data...\n")

if(!file.exists(IN_CLINICAL)) stop(paste("FATAL ERROR: Cannot find clinical master table, please check:", IN_CLINICAL))
if(!file.exists(VOL_FILE)) stop(paste("FATAL ERROR: Cannot find volunteer table, please check:", VOL_FILE))

# [Core Fix]: Force primary keys to Character and Numeric during clinical master load
clin_master <- read_csv(IN_CLINICAL, show_col_types = FALSE) %>%
  mutate(
    FamilyID = as.character(FamilyID), 
    Membercode = as.numeric(Membercode)
  ) %>%
  filter(!is.na(`Fat(%)`)) %>%
  mutate(Phenotype = ifelse(`Fat(%)` >= 30, "Obese", "Lean"))

# Load specified Excel file (and fix potential typos like FaimlyID)
vol_df <- read_excel(VOL_FILE)
colnames(vol_df)[1:2] <- c("FamilyID", "Membercode")
vol_df <- vol_df %>% 
  mutate(
    FamilyID = as.character(FamilyID), 
    Membercode = as.numeric(Membercode)
  )

# Remove potential overlapping columns (keep new data from vol_df)
overlap_cols <- setdiff(intersect(colnames(clin_master), colnames(vol_df)), c("FamilyID", "Membercode"))
if(length(overlap_cols) > 0) clin_master <- clin_master %>% dplyr::select(-all_of(overlap_cols))

# Ultimate merge and factorization (types are now strictly aligned to prevent errors)
merged_df <- clin_master %>%
  inner_join(vol_df, by = c("FamilyID", "Membercode")) %>%
  mutate(
    MemberLabel = factor(Membercode, levels = c(1, 2, 3), labels = c("Daughter", "Mother", "Father")),
    Combo = factor(paste0(Phenotype, "_", MemberLabel),
                   levels = c("Lean_Daughter", "Obese_Daughter", 
                              "Lean_Mother", "Obese_Mother", 
                              "Lean_Father", "Obese_Father"))
  )

# ==============================================================================
# 3. Variable dictionary and fully automated analysis engine
# ==============================================================================
vars_to_plot <- c(
  "Fat mass (kg)", "Lean mass (kg)", "Bone mass (kg)", "BMI", 
  "Visceral fat (kg)", "Subcutaneous fat (kg)", "Liverfat (%)", "Waist circumference (cm)", 
  "TEE (kcal/day)", "REE (kcal/day)", "AEE (kcal/day)", "VO2max (ml/kg/min)", 
  "IEE (kcal/day)", "Fat (E%)", "Protein (E%)", "Carbohydrate (E%)"
)

custom_label_fmt <- function(breaks) {
  valid_breaks <- breaks[!is.na(breaks)]
  is_all_integer <- all(abs(valid_breaks - round(valid_breaks)) < 1e-9)
  
  sapply(breaks, function(x) {
    if (is.na(x)) return("")
    if (is_all_integer) {
      return(as.character(as.integer(round(x)))) 
    } else {
      return(sprintf("%.1f", x)) 
    }
  })
}

create_boxplot <- function(var_name_key) {
  all_cols <- colnames(merged_df)
  clean_target <- tolower(gsub("[^a-zA-Z0-9]", "", var_name_key))
  clean_cols   <- tolower(gsub("[^a-zA-Z0-9]", "", all_cols))
  
  match_idx <- which(clean_cols == clean_target)
  if (length(match_idx) == 0) return(NULL)
  actual_col <- all_cols[match_idx[1]]
  
  df_sub <- merged_df %>% 
    dplyr::select(MemberLabel, Phenotype, Combo, Value = all_of(actual_col)) %>%
    filter(!is.na(Phenotype), !is.na(Value))
  
  df_sub$Value <- as.numeric(as.character(df_sub$Value))
  
  df_sub <- df_sub %>%
    group_by(Combo) %>%
    mutate(
      Q1 = quantile(Value, 0.25, na.rm = TRUE),
      Q3 = quantile(Value, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      IsOutlier = Value < (Q1 - 1.5 * IQR) | Value > (Q3 + 1.5 * IQR)
    ) %>%
    filter(!IsOutlier) %>% 
    ungroup()
  
  if(nrow(df_sub) < 6) return(NULL)
  
  clean_min <- min(df_sub$Value, na.rm = TRUE)
  clean_max <- max(df_sub$Value, na.rm = TRUE)
  
  n_intervals <- 4  
  
  R <- clean_max - clean_min
  if (R == 0) R <- 1
  raw_step <- R / n_intervals
  mag <- 10^floor(log10(raw_step))
  norm_step <- raw_step / mag
  
  if (norm_step <= 1) { nice_mult <- 1 } else
    if (norm_step <= 1.5) { nice_mult <- 1.5 } else
      if (norm_step <= 2) { nice_mult <- 2 } else
        if (norm_step <= 2.5) { nice_mult <- 2.5 } else
          if (norm_step <= 5) { nice_mult <- 5 } else { nice_mult <- 10 }
  
  step <- nice_mult * mag
  if (step < 0.5) step <- 0.5
  
  start_val <- floor(clean_min / step) * step
  
  while (start_val + n_intervals * step < clean_max) {
    if (nice_mult == 1) nice_mult <- 1.5
    else if (nice_mult == 1.5) nice_mult <- 2
    else if (nice_mult == 2) nice_mult <- 2.5
    else if (nice_mult == 2.5) nice_mult <- 5
    else if (nice_mult == 5) nice_mult <- 10
    else nice_mult <- nice_mult * 2 
    
    step <- nice_mult * mag
    if (step < 0.5) step <- 0.5
    start_val <- floor(clean_min / step) * step
  }
  
  my_breaks <- seq(start_val, by = step, length.out = n_intervals + 1)
  y_min_plot <- min(my_breaks)
  y_max_plot <- max(my_breaks)
  
  sig_comps <- list()
  roles <- c("Daughter", "Mother", "Father")
  for(role in roles) {
    grp_obese <- paste0("Obese_", role)
    grp_lean  <- paste0("Lean_", role)
    
    vals_obese <- df_sub$Value[df_sub$Combo == grp_obese]
    vals_lean  <- df_sub$Value[df_sub$Combo == grp_lean]
    
    if(length(vals_obese) >= 3 && length(vals_lean) >= 3) {
      pval <- wilcox.test(vals_obese, vals_lean, exact = FALSE)$p.value
      if(!is.na(pval) && pval < 0.05) {
        sig_comps <- append(sig_comps, list(c(grp_lean, grp_obese)))
      }
    }
  }
  
  y_bracket <- clean_max + (clean_max - clean_min) * 0.05
  
  suppressWarnings({
    p <- ggplot(df_sub, aes(x = Combo, y = Value, fill = Phenotype)) +
      stat_boxplot(geom = "errorbar", width = 0.25, color = "black", linewidth = 0.6) +
      geom_boxplot(width = 0.6, alpha = 0.6, outlier.shape = NA, color = "black", linewidth = 0.6) +
      geom_jitter(aes(color = Phenotype), width = 0.15, size = 2.5, alpha = 0.7, shape = 16) +
      
      scale_fill_manual(values = pheno_colors) +
      scale_color_manual(values = pheno_colors) +
      scale_x_discrete(labels = c("DL", "DO", "ML", "MO", "FL", "FO"), drop = FALSE) +
      
      scale_y_continuous(
        limits = c(y_min_plot, y_max_plot), 
        breaks = my_breaks, labels = custom_label_fmt(my_breaks),
        expand = expansion(mult = c(0.05, 0.05)) 
      ) +
      
      theme_clean +
      theme(aspect.ratio = 1) + 
      labs(title = var_name_key, y = var_name_key) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    
    if(length(sig_comps) > 0) {
      p <- p + stat_compare_means(
        comparisons = sig_comps, method = "wilcox.test", 
        method.args = list(exact = FALSE),
        symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
        label = "p.signif", bracket.size = 0.6, tip.length = 0.02, 
        size = 6, color = "black", vjust = 0.4,
        label.y = rep(y_bracket, length(sig_comps)),
        step.increase = 0 
      )
    }
  })
  
  return(p)
}

# ==============================================================================
# 4. Batch generation and super-patchwork
# ==============================================================================
cat("\n>>> [Step 2] Generating Plots with Smart Decimal Logic...\n")
plot_list <- list()
for(var in vars_to_plot) {
  p <- create_boxplot(var)
  if(!is.null(p)) plot_list[[var]] <- p
}

if(length(plot_list) > 0) {
  for(var in names(plot_list)) {
    safe_name <- gsub("[^A-Za-z0-9]", "_", var)
    ggsave(file.path(OUT_DIR, paste0("Box_Phenotype_", safe_name, ".pdf")), plot_list[[var]], width=5, height=5)
  }
  
  pdf(file.path(OUT_DIR, "Combined_Fig1H_Clinical_Comparison.pdf"), width = 14, height = 10)
  n_plots <- length(plot_list); batch_size <- 6 
  for(i in seq(1, n_plots, by=batch_size)) {
    end_idx <- min(i + batch_size - 1, n_plots)
    subset_plots <- plot_list[i:end_idx]
    do.call(grid.arrange, c(subset_plots, ncol=3))
  }
  dev.off()
  cat("\n>>> Task completed successfully!\nSaved to:", OUT_DIR, "\n")
} else {
  cat("\nError: No plots generated, please check the joined table or variable name matching.\n")
}
