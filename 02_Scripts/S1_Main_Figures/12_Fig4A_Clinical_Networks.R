# ==============================================================================
# Project: FinlandSports - Fig 4A Clinical Excavation (Fat Network & Family)
# Output:  03_Results/Fig_4/Fig4A_Networks
# Upgraded: Partial Spearman ρ, Age/Gender Covariates, FDR, Smart NA Handling
# Features: 100% English, V2.0 Relative Paths, Strict Namespacing, read_csv fix
# ==============================================================================

# --- 0. Environment Setup ---
# setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readr, dplyr, tidyr, ggplot2, reshape2, patchwork, ggsci, 
  FactoMineR, factoextra, umap, stringr, ppcor, ggtext, openxlsx
)

# --- 1. Path & Data Loading (V2.0 Standard) ---
DIR_IN  <- "01_Clean_Data"
DIR_OUT <- "03_Results/Fig_4/Fig4A_Networks"

if(!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)

INPUT_FILE <- file.path(DIR_IN, "Clinical_Master_Strict.csv")
# FIXED: Used read_csv instead of read.csv to prevent mangling "Fat(%)" into "Fat..."
df <- read_csv(INPUT_FILE, show_col_types = FALSE) %>% as.data.frame()

# Basic cleaning for column names
colnames(df) <- gsub("/", ".", colnames(df))
colnames(df) <- gsub(" ", "_", colnames(df))
colnames(df)[colnames(df) == "Android_FM"] <- "Android(FM)"

# ================= 🌟 Ultimate Standard Cleaning Engine (V2.0) 🌟 =================
# Scope: General clinical data processing for Fig 4, Fig 5, Fig 6
# Function: 1. Remove subjective scores 2. Map abbreviations 3. Strip suffixes 4. Grouping
# ==============================================================================

# 1. Protect essential columns from unintended modification
protect_cols <- c("FamilyID", "Membercode", "Fat(%)", "Fat_percent", "HOMA_IR", "Unique_Key", "Age", "Gender")

# 2. Blacklist subjective noise and force removal
noise_blacklist <- c("SelfRatedHealthScore", "Health_History_Score", "Family_History_Score")
df <- df %>% dplyr::select(-dplyr::any_of(noise_blacklist))

# Extract columns allowed to be renamed
rename_idx <- !colnames(df) %in% protect_cols
new_names <- colnames(df)[rename_idx]

# 3. Base cleaning: unify columns representing energy percentages
new_names <- gsub("\\.E\\.\\.", "(E%)", new_names) 
new_names <- gsub("\\.E\\.$", "(E%)", new_names)   

# 4. Ultimate academic abbreviation dictionary (Aligned with AbbreviationList)
rename_dict <- c(
  # Core body composition and phenotype
  "^WAISTLINE$" = "WC",
  "^Liverfat$" = "HFC",
  "^liverfat$" = "HFC",
  "^A\\.L\\.Ratio$" = "A/L Ratio",
  "^Energy_kcal$" = "IEE",
  "^IPA$" = "IPA", 
  
  # Inflammation, Hepatic and Endocrine
  "^CRP$" = "hsCRP",
  "^ALT$" = "S-ALAT",
  "^Testosterone$" = "Te",
  "^Estradiol$" = "E2",
  "^IGF1$" = "IGF-1",
  
  # Osteocalcin family
  "^S_cOC$" = "cOC",
  "^S_TotalOC$" = "tOC",
  "^cOC_totalOC_ratio$" = "cOC/tOC Ratio",
  
  # Muscle strength and fitness
  "^HandGrapeStrength$" = "HGS",
  "^ElbowFlex$" = "EFS",
  "^KneeExtention$" = "NES",
  "^JumpHeight$" = "VJH",
  "^PA_Hour\\.Week$" = "PAD",
  "^PA_Times\\.Week$" = "PAF",
  
  # Fat subpopulations
  "^VISCERAL$" = "Fat Vi",
  "^Visceral$" = "Fat Vi",
  "^Visceralfat$" = "Fat Vi",
  "^Vertebralfat$" = "Fat Ve",
  "^SubcutaneousFat$" = "Fat Su",
  "^IntraperitonealFat$" = "Fat Ip",
  "^RetroperitonealFat$" = "Fat Rp",
  "^Muscle_totalfat$" = "Fat Mu",
  "^Android_FM$" = "FM An",
  "^Android\\(FM\\)$" = "FM An",
  "^Gynoid_FM$" = "FM Gy",
  "^Gynoid\\(FM\\)$" = "FM Gy"
)

# Execute batch replacement
for(pattern in names(rename_dict)) {
  new_names <- str_replace_all(new_names, pattern, rename_dict[pattern])
}
colnames(df)[rename_idx] <- new_names

# 5. New Mechanism Grouping Dictionary
var_map <- data.frame(
  Variable = c(
    "cOC", "IGF-1", "NES",          
    "Fat Vi", "Fat Ve", "A/L Ratio", 
    "Leptin", "AST", "CHOL",                    
    "REE", "AEE", "IPA"                         
  ),
  Group = c(
    rep("Bone-Muscle Axis", 3),               
    rep("Ectopic Adiposity", 3),              
    rep("Hepatic & Systemic Metabolism", 3),  
    rep("Energy & Lifestyle", 3)              
  )
)
# ==============================================================================

message(">>> Generating Combined Network vs HOMA-IR (Obese vs Lean)...")

# --- 2. Advanced Data Preprocessing ---
# Ensure Age and Gender are numeric to be used as covariates for partial correlation
if("Age" %in% colnames(df)) df$Age <- as.numeric(df$Age)
if("Gender" %in% colnames(df)) df$Gender <- as.numeric(as.factor(df$Gender))

# Dynamically locate the fat column
fat_col <- "Fat(%)"
if(!"Fat(%)" %in% colnames(df) && "Fat_percent" %in% colnames(df)) {
  fat_col <- "Fat_percent"
}

if(!fat_col %in% colnames(df)) {
  stop("FATAL ERROR: Could not find Fat percentage column. Please check Clinical_Master_Strict.csv headers.")
}

df_obese <- df %>% dplyr::filter(.data[[fat_col]] >= 30)
df_lean  <- df %>% dplyr::filter(.data[[fat_col]] < 30)

# Grab numeric variables (exclude IDs, grouping criteria, and covariates)
ignore_cols <- c("FamilyID", "Membercode", fat_col, "HOMA_IR", "Unique_Key", "Age", "Gender")
num_cols <- names(df)[sapply(df, is.numeric)]
test_vars <- setdiff(num_cols, ignore_cols)

# --- 3. Core Calculation Module (Partial Spearman + FDR + Smart NA Handling) ---
calc_partial_with_cov <- function(df_sub) {
  res_list <- list()
  for(col in test_vars) {
    tmp <- df_sub %>% dplyr::select(HOMA_IR, dplyr::all_of(col), Age, Gender) %>% na.omit()
    n_samples <- nrow(tmp)
    
    if (n_samples > 10 && sd(tmp[[col]]) > 0 && sd(tmp$HOMA_IR) > 0) {
      try({
        if(length(unique(tmp$Gender)) > 1) {
          ct <- pcor.test(tmp[[col]], tmp$HOMA_IR, tmp[, c("Age", "Gender")], method = "spearman")
        } else {
          ct <- pcor.test(tmp[[col]], tmp$HOMA_IR, tmp[, c("Age")], method = "spearman")
        }
        res_list[[col]] <- data.frame(
          Variable = col, 
          Cor = ct$estimate, 
          P_Value = ct$p.value,
          N = n_samples
        )
      }, silent=TRUE)
    }
  }
  
  res_df <- dplyr::bind_rows(res_list)
  if(nrow(res_df) > 0) {
    res_df$FDR <- p.adjust(res_df$P_Value, method = "BH")
  }
  return(res_df)
}

message("  > Calculating Obese Cohort...")
stats_obese_full <- calc_partial_with_cov(df_obese)
message("  > Calculating Lean Cohort...")
stats_lean_full  <- calc_partial_with_cov(df_lean)

# --- 4 & 5 & 6 & 7. Variable Extraction and Dual Plotting (Functionalized Refactoring) ---

valid_both <- intersect(stats_obese_full$Variable, stats_lean_full$Variable)
sig_obese_full <- stats_obese_full %>% dplyr::filter(P_Value < 0.05) %>% dplyr::pull(Variable)
sig_lean_full  <- stats_lean_full %>% dplyr::filter(P_Value < 0.05) %>% dplyr::pull(Variable)
sig_vars_full_union <- intersect(union(sig_obese_full, sig_lean_full), valid_both)

target_top12 <- c(
  "cOC", "IGF-1", "NES", 
  "Fat Vi", "Fat Ve", "A/L Ratio", 
  "Leptin", "AST", "CHOL", 
  "REE", "AEE", "IPA"
)
sig_vars_golden12 <- intersect(target_top12, valid_both)

new_vars <- setdiff(sig_vars_full_union, var_map$Variable)
if(length(new_vars) > 0) {
  var_map <- dplyr::bind_rows(var_map, data.frame(Variable = new_vars, Group = "New_Excavated_Traits"))
}

generate_split_network <- function(selected_vars, file_suffix) {
  
  add_stars_and_map <- function(df_stats) {
    df_stats %>% 
      dplyr::filter(Variable %in% selected_vars) %>%
      dplyr::mutate(Sig = dplyr::case_when(P_Value < 0.001 ~ "***", P_Value < 0.01 ~ "**", P_Value < 0.05 ~ "*", TRUE ~ "")) %>%
      dplyr::left_join(var_map, by = "Variable")
  }
  
  biological_module_order <- c(
    "cOC", "IGF-1", "NES",
    "Fat Vi", "Fat Ve", "A/L Ratio",
    "Leptin", "AST", "CHOL",
    "REE", "AEE", "IPA"
  )
  
  if (file_suffix == "Golden12") {
    g_order <- intersect(biological_module_order, selected_vars)
  } else {
    tmp_lean <- add_stars_and_map(stats_lean_full) %>% dplyr::arrange(Cor)
    g_order <- tmp_lean$Variable
  }
  
  s_lean <- add_stars_and_map(stats_lean_full)
  s_lean$Variable <- factor(s_lean$Variable, levels = g_order)
  s_obese <- add_stars_and_map(stats_obese_full)
  s_obese$Variable <- factor(s_obese$Variable, levels = g_order)
  
  m_cor <- max(abs(c(s_obese$Cor, s_lean$Cor)), na.rm = TRUE)
  b_style <- element_rect(color = "black", fill = NA, linewidth = 1)
  
  # Plotting Module 1: Lean Heatmap
  c_L <- cor(df_lean[, as.character(g_order)], use = "pairwise.complete.obs", method = "spearman")
  m_L <- melt(c_L); m_L$Var1 <- factor(m_L$Var1, levels = g_order); m_L$Var2 <- factor(m_L$Var2, levels = g_order)
  p_h_L <- ggplot(m_L, aes(Var1, Var2, fill = value)) + geom_tile(color = "black", linewidth = 0.5) +
    scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#DC0000", midpoint=0, limits=c(-1,1), name=expression("Partial "*rho)) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme_minimal(base_size = 14) +
    labs(title = "Collinearity Matrix (Lean)") +
    theme(aspect.ratio = 1, panel.border = b_style, 
          axis.text.x = element_markdown(angle = 45, hjust = 1, color="black", size=15), 
          axis.text.y = element_markdown(color="black", size=14),                       
          axis.title = element_blank(),
          plot.title = element_text(face="bold", hjust=0.5, size=15), legend.position = "none", panel.grid = element_blank())
  
  # Plotting Module 2: Lean Lollipop
  p_b_L <- ggplot(s_lean, aes(x = Variable, y = Cor)) +
    geom_segment(aes(x=Variable, xend=Variable, y=0, yend=Cor), color="grey60", linewidth=1.2) +
    geom_point(aes(fill=Group, size=abs(Cor)), shape=21, color="white", stroke=0.8) +
    geom_text(aes(label=Sig, y=Cor, hjust=ifelse(Cor>=0, -2, 2)), size=11, color="black", vjust=0.75) +
    scale_fill_npg() + scale_size(range = c(5, 11), guide = "none") + geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_y_continuous(limits = c(-m_cor, m_cor), expand = expansion(mult = c(0.45, 0.45))) + coord_flip() + theme_classic(base_size = 14) +
    labs(title = "Associations with HOMA-IR (Lean)", y = expression("Partial "*rho), x = "") + 
    theme(aspect.ratio = 1, panel.border = b_style, axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(),
          axis.text.x = element_markdown(color="black", size=15),  
          plot.title = element_text(face="bold", hjust=0.5, size=15),
          legend.position = "bottom", legend.title = element_text(face="bold"))
  
  # Plotting Module 3: Obese Heatmap
  c_O <- cor(df_obese[, as.character(g_order)], use = "pairwise.complete.obs", method = "spearman")
  m_O <- melt(c_O); m_O$Var1 <- factor(m_O$Var1, levels = g_order); m_O$Var2 <- factor(m_O$Var2, levels = g_order)
  p_h_O <- ggplot(m_O, aes(Var1, Var2, fill = value)) + geom_tile(color = "black", linewidth = 0.5) +
    scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#DC0000", midpoint=0, limits=c(-1,1), name=expression("Partial "*rho)) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme_minimal(base_size = 14) +
    labs(title = "Collinearity Matrix (Obese)") +
    theme(aspect.ratio = 1, panel.border = b_style, 
          axis.text.x = element_markdown(angle = 45, hjust = 1, color="black", size=15),
          axis.text.y = element_markdown(color="black", size=14), 
          axis.title = element_blank(),
          plot.title = element_text(face="bold", hjust=0.5, size=15), legend.position = "none", panel.grid = element_blank())
  
  # Plotting Module 4: Obese Lollipop
  p_b_O <- ggplot(s_obese, aes(x = Variable, y = Cor)) +
    geom_segment(aes(x=Variable, xend=Variable, y=0, yend=Cor), color="grey60", linewidth=1.2) +
    geom_point(aes(fill=Group, size=abs(Cor)), shape=21, color="white", stroke=0.8) +
    geom_text(aes(label=Sig, y=Cor, hjust=ifelse(Cor>=0, -2, 2)), size=11, color="black", vjust=0.75) +
    scale_fill_npg() + scale_size(range = c(5, 11), guide = "none") + geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_y_continuous(limits = c(-m_cor, m_cor), expand = expansion(mult = c(0.45, 0.45))) + coord_flip() + theme_classic(base_size = 14) +
    labs(title = "Associations with HOMA-IR (Obese)", y = expression("Partial "*rho), x = "") +
    theme(aspect.ratio = 1, panel.border = b_style, axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(),
          axis.text.x = element_markdown(color="black", size=15), 
          plot.title = element_text(face="bold", hjust=0.5, size=15), legend.position = "none")
  
  # Ultimate Assembly and Output
  p_final <- (p_h_L | p_b_L | p_h_O | p_b_O) + plot_layout(guides = 'collect') &
    theme(legend.position = 'bottom', legend.box = "horizontal", legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 13)) &
    guides(fill = guide_legend(override.aes = list(size = 8)))
  
  out_name <- paste0("Fig4A_Metabolic_SplitNetwork_", file_suffix, ".pdf")
  ggsave(file.path(DIR_OUT, out_name), p_final, width = 22, height = 8)
  message(">>> successfully saved: ", out_name)
}

# ====================================================================
# Execute Dual Task Publication
# ====================================================================
message("\n>>> Generating Comparison Version (All Significant Traits)...")
generate_split_network(sig_vars_full_union, "FULL_Comparison")

message("\n>>> Generating Main Figure Version (Golden 12 Core)...")
generate_split_network(sig_vars_golden12, "Golden12")

# ==============================================================================
# 8. Special Analysis: Lifestyle Independent Network and Lollipop Chart
# ==============================================================================
message("\n>>> Generating Lifestyle-Specific Network (Fig4A_Lifestyle)...")

# Replace "Energy(kcal)" with "IEE" and define lifestyle vars
base_life_vars <- c("PAF", "PAD", "IPA", 
                    "IEE", "Smoking", "Alcohol(%)", "Alcohol(E%)",  
                    "Protein(%)", "Fat(%)", "Safa(%)", "Mufa(%)", "Pufa(%)", "Carbohydrate(%)", "Sucrose(%)",
                    "Protein(E%)", "Fat(E%)", "Safa(E%)", "Mufa(E%)", "Pufa(E%)", "Carbohydrate(E%)", "Sucrose(E%)")

lifestyle_vars <- intersect(base_life_vars, colnames(df))

# Retain all Lifestyle variables present
sig_life_union <- intersect(lifestyle_vars, valid_both)

add_stars_life <- function(df_stats, cohort_name) {
  df_stats %>% 
    dplyr::filter(Variable %in% sig_life_union) %>%
    dplyr::mutate(
      Sig = dplyr::case_when(
        P_Value < 0.001 ~ "***",
        P_Value < 0.01  ~ "**",
        P_Value < 0.05  ~ "*",
        TRUE ~ ""
      ),
      Group = cohort_name 
    )
}

stats_lean_life  <- add_stars_life(stats_lean_full, "Lean") %>% dplyr::arrange(Cor)
golden_order_life <- stats_lean_life$Variable
stats_lean_life$Variable <- factor(stats_lean_life$Variable, levels = golden_order_life)

stats_obese_life <- add_stars_life(stats_obese_full, "Obese")
stats_obese_life$Variable <- factor(stats_obese_life$Variable, levels = golden_order_life)

# Plotting parameters
max_cor_life <- max(abs(c(stats_obese_life$Cor, stats_lean_life$Cor)), na.rm = TRUE) * 1.1
border_style <- element_rect(color = "black", fill = NA, linewidth = 1)

# ----------------- [Obese Plots] -----------------
cor_mat_O_life <- cor(df_obese[, as.character(golden_order_life)], use = "pairwise.complete.obs", method = "spearman")
melt_O_life <- melt(cor_mat_O_life)
melt_O_life$Var1 <- factor(melt_O_life$Var1, levels = golden_order_life)
melt_O_life$Var2 <- factor(melt_O_life$Var2, levels = golden_order_life)

p_heat_O_life <- ggplot(melt_O_life, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#DC0000", midpoint=0, limits=c(-1,1), name=expression("Partial "*rho)) +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal(base_size = 14) + labs(title = "Lifestyle Collinearity (Obese)") +
  theme(aspect.ratio = 1, panel.border = border_style, 
        axis.text.x = element_markdown(angle = 45, hjust = 1, color="black", size=15),
        axis.text.y = element_markdown(color="black", size=14), axis.title = element_blank(),
        plot.title = element_text(face="bold", hjust=0.5, size=15), legend.position = "none", panel.grid = element_blank())

p_bar_O_life <- ggplot(stats_obese_life, aes(x = Variable, y = Cor)) +
  geom_segment(aes(x=Variable, xend=Variable, y=0, yend=Cor), color="grey60", linewidth=1.2) +
  geom_point(aes(fill=Group, size=abs(Cor)), shape=21, color="white", stroke=0.8) +
  geom_text(aes(label=Sig, y=Cor, hjust=ifelse(Cor>=0, -1.5, 1.5)), size=11, color="black", vjust=0.75) +
  scale_fill_manual(name = "Cohort", values = c("Obese" = "#F79647", "Lean" = "#3DA6AE")) + scale_size(range = c(5, 11), guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(-max_cor_life, max_cor_life), expand = expansion(mult = c(0.4, 0.4))) +
  coord_flip() + theme_classic(base_size = 14) + labs(title = "HOMA-IR Associations (Obese)", y = expression("Partial "*rho), x = "") +
  theme(aspect.ratio = 1, panel.border = border_style, axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.line = element_blank(), axis.text.x = element_markdown(color="black", size=15),
        plot.title = element_text(face="bold", hjust=0.5, size=15), legend.position = "none")

# ----------------- [Lean Plots] -----------------
cor_mat_L_life <- cor(df_lean[, as.character(golden_order_life)], use = "pairwise.complete.obs", method = "spearman")
melt_L_life <- melt(cor_mat_L_life)
melt_L_life$Var1 <- factor(melt_L_life$Var1, levels = golden_order_life)
melt_L_life$Var2 <- factor(melt_L_life$Var2, levels = golden_order_life)

p_heat_L_life <- ggplot(melt_L_life, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#DC0000", midpoint=0, limits=c(-1,1), name=expression("Partial "*rho)) +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal(base_size = 14) + labs(title = "Lifestyle Collinearity (Lean)") +
  theme(aspect.ratio = 1, panel.border = border_style, 
        axis.text.x = element_markdown(angle = 45, hjust = 1, color="black", size=15),
        axis.text.y = element_markdown(color="black", size=14), axis.title = element_blank(),
        plot.title = element_text(face="bold", hjust=0.5, size=15), legend.position = "none", panel.grid = element_blank())

p_bar_L_life <- ggplot(stats_lean_life, aes(x = Variable, y = Cor)) +
  geom_segment(aes(x=Variable, xend=Variable, y=0, yend=Cor), color="grey60", linewidth=1.2) +
  geom_point(aes(fill=Group, size=abs(Cor)), shape=21, color="white", stroke=0.8) +
  geom_text(aes(label=Sig, y=Cor, hjust=ifelse(Cor>=0, -1.5, 1.5)), size=11, color="black", vjust=0.75) +
  scale_fill_manual(name = "Cohort", values = c("Obese" = "#F79647", "Lean" = "#3DA6AE")) + scale_size(range = c(5, 11), guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(-max_cor_life, max_cor_life), expand = expansion(mult = c(0.4, 0.4))) +
  coord_flip() + theme_classic(base_size = 14) + labs(title = "HOMA-IR Associations (Lean)", y = expression("Partial "*rho), x = "") +
  theme(aspect.ratio = 1, panel.border = border_style, axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.line = element_blank(), axis.text.x = element_markdown(color="black", size=15),
        plot.title = element_text(face="bold", hjust=0.5, size=15), legend.position = "bottom", legend.title = element_text(face="bold"))

# 5. Assembly and Output
p_life_final <- (p_heat_L_life | p_bar_L_life | p_heat_O_life | p_bar_O_life) + 
  plot_layout(guides = 'collect') &
  theme(
    legend.position = 'bottom', 
    legend.box = "horizontal",
    legend.title = element_text(size = 14, face = "bold"), 
    legend.text = element_text(size = 13)                 
  ) &
  guides(fill = guide_legend(override.aes = list(size = 8)))

ggsave(file.path(DIR_OUT, "Fig4A_Lifestyle_SplitNetwork.pdf"), p_life_final, width = 22, height = 8)
message(">>> Successfully generated Lifestyle Analysis (Fig4A).")

# ==============================================================================
# 9. Generate Supplementary Table for Fig4A (Comprehensive Clinical Associations)
# ==============================================================================
message("\n>>> Generating Publication-Ready Supplementary Table: Data_S7_Clinical_Associations...")

# ------------------------------------------------------------------------------
# 1. Data arrangement and sorting
# ------------------------------------------------------------------------------
df_sup_lean <- stats_lean_full %>%
  dplyr::arrange(P_Value) %>%
  dplyr::mutate(Significance = dplyr::case_when(
    P_Value < 0.001 ~ "***",
    P_Value < 0.01  ~ "**",
    P_Value < 0.05  ~ "*",
    TRUE ~ "ns"
  )) %>%
  dplyr::select(Variable, Partial_Rho = Cor, P_Value, adj.P_Val_FDR = FDR, Significance, N_Samples = N)

df_sup_obese <- stats_obese_full %>%
  dplyr::arrange(P_Value) %>%
  dplyr::mutate(Significance = dplyr::case_when(
    P_Value < 0.001 ~ "***",
    P_Value < 0.01  ~ "**",
    P_Value < 0.05  ~ "*",
    TRUE ~ "ns"
  )) %>%
  dplyr::select(Variable, Partial_Rho = Cor, P_Value, adj.P_Val_FDR = FDR, Significance, N_Samples = N)

# ------------------------------------------------------------------------------
# 2. Create Excel Workbook and Save
# ------------------------------------------------------------------------------
wb_fig4a <- createWorkbook()

header_style <- createStyle(
  fontName = "Arial", fontSize = 11, fontColour = "white", 
  fgFill = "#4F81BD", textDecoration = "bold", halign = "center"
)

addWorksheet(wb_fig4a, "Lean_Cohort")
writeData(wb_fig4a, "Lean_Cohort", df_sup_lean)
addStyle(wb_fig4a, "Lean_Cohort", header_style, rows = 1, cols = 1:ncol(df_sup_lean), gridExpand = TRUE)

addWorksheet(wb_fig4a, "Obese_Cohort")
writeData(wb_fig4a, "Obese_Cohort", df_sup_obese)
addStyle(wb_fig4a, "Obese_Cohort", header_style, rows = 1, cols = 1:ncol(df_sup_obese), gridExpand = TRUE)

setColWidths(wb_fig4a, "Lean_Cohort", cols = 1:ncol(df_sup_lean), widths = "auto")
setColWidths(wb_fig4a, "Obese_Cohort", cols = 1:ncol(df_sup_obese), widths = "auto")

saveWorkbook(wb_fig4a, file.path(DIR_OUT, "Data_S7_Clinical_Associations.xlsx"), overwrite = TRUE)

message("  -> [SUCCESS] Data_S7_Clinical_Associations.xlsx successfully saved to: ", DIR_OUT)
message("========================================================================")
