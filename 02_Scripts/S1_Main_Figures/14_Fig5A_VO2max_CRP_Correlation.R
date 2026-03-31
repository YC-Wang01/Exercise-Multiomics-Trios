# ==============================================================================
# Project: FinlandSports V2.0 - Fig 5A VO2max vs CRP Correlations
# Description: 
#   Generates 6 Plots in total (Categorical & Continuous Analysis):
#   [Categorical Analysis - Phenotype]
#     1. Violin_Grouped: Standard Median Split (High vs Low).
#     2. Violin_Boxplot: Enhanced with Boxplot overlay.
#     3. Violin_Q1vsQ4: Extreme Contrast (Lowest 25% vs Highest 25%).
#   [Continuous Analysis - Mechanism]
#     4. Scatter_Basic: Simple Log(CRP) vs VO2max.
#     5. Scatter_Marginal: Scatter with Density plots on axes.
#     6. Scatter_Binned: Scatter with Binned Means Trend.
# Features: V2.0 Path Standard, Single Master Clinical Table, Strict Namespacing
# ==============================================================================

# --- 0. Environment Setup ---
# setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, dplyr, tidyr, ggplot2, ggpubr, cowplot, scales, grid, stringr, openxlsx)

# --- 1. CONFIGURATION & THEMES ---
DIR_DATA <- "01_Clean_Data"
DIR_OUT  <- "03_Results/Fig_5/Fig5A_VO2max_CRP"

if(!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)

# Colors
COL_LOW  <- "#3C5488" # Blue (Healthy / Low CRP)
COL_HIGH <- "#E64B35" # Red (Inflammation / High CRP)
GRADIENT_COLS <- c("#4DBBD5", "#E64B35")

# Theme for Violin Plots
THEME_VIOLIN <- theme_classic(base_family = "sans", base_size = 14) +
  theme(
    axis.line = element_line(linewidth = 1),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "none"
  )

# Theme for Scatter Plots
THEME_SCATTER <- theme_bw(base_family = "sans", base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.title = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

# --- 2. DATA LOADING & PROCESSING (V2.0 Engine) ---
message(">>> Loading Clinical Data...")

# Direct extraction from the unified Strict Master Table
df_merged <- read_csv(file.path(DIR_DATA, "Clinical_Master_Strict.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(
    Role_Orig = dplyr::case_when(Membercode == 1 ~ "Daughter", Membercode == 2 ~ "Mother", Membercode == 3 ~ "Father", TRUE ~ "Unknown"),
    Subject_ID = tolower(paste0(FamilyID, "_", Role_Orig))
  ) %>%
  dplyr::select(Subject_ID, CRP, VO2max) %>%
  tidyr::drop_na(CRP, VO2max) %>%
  dplyr::mutate(Log_CRP = log(CRP))

message(paste(">>> Data Ready. Total Valid Subjects:", nrow(df_merged)))

# --- 3. CATEGORICAL PLOTS (Violins) ---
message(">>> Generating Categorical Plots...")

# Prepare Median Split Data
median_val <- median(df_merged$CRP)
df_median <- df_merged %>%
  dplyr::mutate(
    Group = ifelse(CRP > median_val, "High CRP", "Low CRP"),
    Group = factor(Group, levels = c("Low CRP", "High CRP"))
  )

# ==============================================================================
# Calculate Statistics & Effect Size (Cohen's d) and Export
# ==============================================================================
vo2_low <- df_median$VO2max[df_median$Group == "Low CRP"]
vo2_high <- df_median$VO2max[df_median$Group == "High CRP"]

mean_low <- mean(vo2_low, na.rm = TRUE)
mean_high <- mean(vo2_high, na.rm = TRUE)
mean_diff <- mean_low - mean_high

n_low <- length(na.omit(vo2_low))
n_high <- length(na.omit(vo2_high))
var_low <- var(vo2_low, na.rm = TRUE)
var_high <- var(vo2_high, na.rm = TRUE)
pooled_sd <- sqrt(((n_low - 1) * var_low + (n_high - 1) * var_high) / (n_low + n_high - 2))
cohens_d <- abs(mean_diff) / pooled_sd

wilcox_res <- wilcox.test(VO2max ~ Group, data = df_median)

stats_summary <- data.frame(
  Grouping_Method = paste0("Median Split (CRP cutoff = ", round(median_val, 2), " mg/L)"),
  Mean_Low_CRP = round(mean_low, 2),
  Mean_High_CRP = round(mean_high, 2),
  Mean_Difference = round(mean_diff, 2),
  Cohens_d = round(cohens_d, 3),
  P_Value_Wilcoxon = signif(wilcox_res$p.value, 4)
)

openxlsx::write.xlsx(
  list(Statistics = stats_summary, SourceData = df_median),
  file.path(DIR_OUT, "Table_S_Fig5A_VO2max_Stats.xlsx")
)
message("  > Excel Table Saved: Table_S_Fig5A_VO2max_Stats.xlsx")

# --- Plot A1: Standard Grouped Violin (Median Split) ---
p_a1 <- ggplot(df_median, aes(x = Group, y = VO2max, color = Group, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0, linewidth = 1, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.2, size = 3, shape = 16, alpha = 0.8) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center", vjust = -1, size = 5) +
  scale_color_manual(values = c("Low CRP" = COL_LOW, "High CRP" = COL_HIGH)) +
  scale_fill_manual(values = c("Low CRP" = COL_LOW, "High CRP" = COL_HIGH)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  labs(title = "VO2max: Low vs High CRP\n(Standard)", x = NULL, y = "VO2max (ml/kg/min)") +
  THEME_VIOLIN

ggsave(file.path(DIR_OUT, "Fig5A_Plot_A1_Violin_Grouped.pdf"), p_a1, width = 4, height = 5)

# --- Plot A2: Violin + Boxplot Enhanced (Median Split) ---
p_a2 <- ggplot(df_median, aes(x = Group, y = VO2max, color = Group, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.1, linewidth = 0.8) +
  geom_boxplot(width = 0.2, color = "black", alpha = 0.8, fill = "white", outlier.shape = NA, linewidth = 0.8) +
  geom_jitter(width = 0.15, size = 2, shape = 16, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center", vjust = -1, size = 5) +
  scale_color_manual(values = c("Low CRP" = COL_LOW, "High CRP" = COL_HIGH)) +
  scale_fill_manual(values = c("Low CRP" = COL_LOW, "High CRP" = COL_HIGH)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  labs(title = "VO2max: Low vs High CRP\n(Boxplot Enhanced)", x = NULL, y = "VO2max (ml/kg/min)") +
  THEME_VIOLIN

ggsave(file.path(DIR_OUT, "Fig5A_Plot_A2_Violin_Boxplot.pdf"), p_a2, width = 4, height = 5)

# --- Plot A3: Q1 vs Q4 Extreme Contrast ---
qs <- quantile(df_merged$CRP, probs = c(0.25, 0.75))
df_quartile <- df_merged %>%
  dplyr::mutate(Group = dplyr::case_when(
    CRP <= qs[1] ~ "Q1 (Lowest)",
    CRP >= qs[2] ~ "Q4 (Highest)",
    TRUE ~ "Middle"
  )) %>%
  dplyr::filter(Group != "Middle") %>%
  dplyr::mutate(Group = factor(Group, levels = c("Q1 (Lowest)", "Q4 (Highest)")))

p_a3 <- ggplot(df_quartile, aes(x = Group, y = VO2max, color = Group, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.1, linewidth = 0.8) +
  geom_boxplot(width = 0.2, color = "black", alpha = 0.8, fill = "white", outlier.shape = NA, linewidth = 0.8) +
  geom_jitter(width = 0.15, size = 2, shape = 16, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center", vjust = -1, size = 5) +
  scale_color_manual(values = c("Q1 (Lowest)" = COL_LOW, "Q4 (Highest)" = COL_HIGH)) +
  scale_fill_manual(values = c("Q1 (Lowest)" = COL_LOW, "Q4 (Highest)" = COL_HIGH)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  labs(title = "Extreme Contrast\n(Lowest 25% vs Highest 25%)", x = NULL, y = "VO2max (ml/kg/min)") +
  THEME_VIOLIN

ggsave(file.path(DIR_OUT, "Fig5A_Plot_A3_Violin_Q1vsQ4.pdf"), p_a3, width = 4, height = 5)

# --- 4. CONTINUOUS PLOTS (Scatters) ---
message(">>> Generating Continuous Plots...")

# Common Statistics
stats_fit <- cor.test(df_merged$Log_CRP, df_merged$VO2max)
stats_label <- paste0("P = ", signif(stats_fit$p.value, 2), "\nR = ", signif(stats_fit$estimate, 2))

# --- Plot B1: Basic Continuous Scatter ---
p_b1 <- ggplot(df_merged, aes(x = Log_CRP, y = VO2max)) +
  geom_point(aes(fill = CRP), shape = 21, size = 4, color = "black", stroke = 0.5, alpha = 0.9) +
  scale_fill_gradientn(colors = GRADIENT_COLS, name = "CRP (mg/L)") +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", fill = "grey", alpha = 0.2, linewidth = 1) +
  annotate("text", x = Inf, y = Inf, label = stats_label, hjust = 1.1, vjust = 1.5, size = 4.5, fontface = "italic") +
  labs(title = "Association between\nCRP and VO2max", x = "Log(CRP)", y = "VO2max (ml/kg/min)") +
  THEME_SCATTER

ggsave(file.path(DIR_OUT, "Fig5A_Plot_B1_Scatter_Basic.pdf"), p_b1, width = 5.5, height = 5)

# --- Plot B2: Marginal Densities ---
p_main <- p_b1 + theme(legend.position = "none")
p_top <- ggplot(df_merged, aes(x = Log_CRP)) +
  geom_density(fill = COL_HIGH, alpha = 0.3, color = NA) + theme_void() + theme(plot.margin = margin(0, 0, -5, 0))
p_right <- ggplot(df_merged, aes(x = VO2max)) +
  geom_density(fill = COL_LOW, alpha = 0.3, color = NA) + coord_flip() + theme_void() + theme(plot.margin = margin(0, 0, 0, -5))

p_b2 <- cowplot::plot_grid(
  p_top, NULL, 
  p_main, p_right, 
  ncol = 2, nrow = 2, align = "hv", 
  rel_widths = c(3, 1), rel_heights = c(1, 3), axis = "lb"
)

ggsave(file.path(DIR_OUT, "Fig5A_Plot_B2_Scatter_Marginal.pdf"), p_b2, width = 6, height = 6)

# --- Plot B3: Binned Trends (The "Scatter Rescue") ---
df_binned <- df_merged %>%
  dplyr::mutate(Bin = ggplot2::cut_number(Log_CRP, n = 5)) %>%
  dplyr::group_by(Bin) %>%
  dplyr::summarise(
    Mean_X = mean(Log_CRP),
    Mean_Y = mean(VO2max),
    SE_Y = sd(VO2max) / sqrt(n()),
    Count = n(),
    .groups = 'drop'
  )

p_b3 <- ggplot(df_merged, aes(x = Log_CRP, y = VO2max)) +
  geom_point(aes(fill = CRP), shape = 21, size = 3, color = "grey80", stroke = 0, alpha = 0.3) +
  geom_errorbar(data = df_binned, aes(x = Mean_X, ymin = Mean_Y - SE_Y, ymax = Mean_Y + SE_Y), 
                width = 0.05, linewidth = 0.8, color = "black", inherit.aes = FALSE) +
  geom_line(data = df_binned, aes(x = Mean_X, y = Mean_Y), 
            color = "black", linewidth = 1, inherit.aes = FALSE) +
  geom_point(data = df_binned, aes(x = Mean_X, y = Mean_Y), 
             size = 4, shape = 21, fill = "white", color = "black", stroke = 1.5, inherit.aes = FALSE) +
  annotate("text", x = Inf, y = Inf, label = stats_label, hjust = 1.1, vjust = 1.5, fontface = "italic") +
  scale_fill_gradientn(colors = GRADIENT_COLS, name = "CRP") +
  labs(title = "Average Trend: Inflammation vs VO2max", subtitle = "Points: Individual Values | Line: Binned Average ± SEM",
       x = "Log(CRP)", y = "VO2max (ml/kg/min)") +
  THEME_SCATTER

ggsave(file.path(DIR_OUT, "Fig5A_Plot_B3_Scatter_BinnedTrend.pdf"), p_b3, width = 6, height = 5)

# --- 5. SUMMARY ---
message("\n=======================================================")
message(">>> FIG 5A PLOTS GENERATED SUCCESSFULLY!")
message(paste(">>> Output Folder:", DIR_OUT))
message("=======================================================")
