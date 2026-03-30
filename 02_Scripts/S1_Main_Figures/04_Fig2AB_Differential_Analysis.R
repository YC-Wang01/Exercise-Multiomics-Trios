# ==============================================================================
# Project: FinlandSports V2.0 
# Script:  04_Fig2AB_Full_Differential_Analysis.R
# Panel:   Fig. 2A, 2B & Baseline Landscape/Volcanoes
# Style:   Pure White Background, Black Axis, 100% Original Logic Replication
# ==============================================================================

# [0. Lock working directory]
setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

# 1. Environment setup and package loading
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, dplyr, tidyr, stringr, ggplot2, ggsci, ggpubr, ggrepel, openxlsx, limma, ggrastr)

# [Rule 1: Standardized path configuration]
INPUT_DIR  <- "01_Clean_Data"
CLIN_FILE  <- file.path(INPUT_DIR, "Clinical_Master_Strict.csv")
OUT_DIR    <- "03_Results/Fig_2"
if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# [Core style customization]: Pure black and white academic style
my_clean_theme <- theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),       
    panel.grid.minor = element_blank(),       
    axis.line = element_line(color = "black", linewidth = 0.8), 
    axis.ticks = element_line(color = "black"),                 
    axis.text = element_text(color = "black", face = "bold"),   
    axis.title = element_text(color = "black", face = "bold"),
    strip.background = element_blank(),                         
    strip.text = element_text(color = "black", face = "bold", size = 16),
    panel.border = element_blank(),                             
    legend.position = "none",
    aspect.ratio = 1                                            
  )

# ==============================================================================
# 2. Clinical data loading (1:1 replicated ID construction logic)
# ==============================================================================
cat("\n>>> [Step 1] Loading Clinical Data...")
clin <- read_csv(CLIN_FILE, show_col_types = FALSE) %>%
  filter(!is.na(`Fat(%)`)) %>%
  mutate(
    Role_Label = case_when(Membercode == 1 ~ "Daughter", Membercode == 2 ~ "Mother", Membercode == 3 ~ "Father", TRUE ~ "Unknown"),
    # Strictly replicate original ID: FamilyID_Role (lowercase)
    Subject_ID = tolower(paste0(FamilyID, "_", Role_Label)),
    Ind_Phenotype = ifelse(`Fat(%)` >= 30, "Obese", "Lean")
  )

# ==============================================================================
# 3. PART I: Limma paired differential analysis sweep (Pre vs Post exercise)
# ==============================================================================
omics_files <- list.files(INPUT_DIR, pattern = "^Cleaned_.*\\.csv", full.names = TRUE)
all_diff_results <- list()
summary_counts <- list()

for (f in omics_files) {
  ds_name <- str_remove(basename(f), "^Cleaned_") %>% str_remove("\\.csv$")
  if(grepl("Clinical|Sample", ds_name, ignore.case = TRUE)) next
  
  df <- read_csv(f, show_col_types = FALSE)
  mat <- as.matrix(df[,-1]); rownames(mat) <- as.character(df[[1]])
  if(max(mat, na.rm = TRUE) > 100) mat <- log2(mat + 1)
  
  # Original sample parsing logic
  col_meta <- data.frame(OriginalCol = colnames(mat), stringsAsFactors = FALSE) %>%
    mutate(
      Clean_Name = tolower(gsub("_twin\\.[0-9]+", "", OriginalCol)),
      TimeType = case_when(grepl("pre", Clean_Name) ~ "pre", grepl("fast", Clean_Name) ~ "fast",
                           grepl("post3h", Clean_Name) ~ "post3h", TRUE ~ "Other"),
      Subject_ID = gsub("_fast|_pre|_post1h|_post3h", "", Clean_Name)
    ) %>% inner_join(clin, by = "Subject_ID")
  
  # Strict Baseline & Post3h selection
  base_cands <- col_meta %>% filter(TimeType %in% c("pre", "fast")) %>% arrange(Subject_ID, TimeType) %>% 
    group_by(Subject_ID) %>% dplyr::slice(1) %>% mutate(Time = "Baseline")
  post_cands <- col_meta %>% filter(TimeType == "post3h") %>% mutate(Time = "Post3h")
  
  col_meta_clean <- bind_rows(base_cands, post_cands) %>% ungroup()
  paired_subjects <- col_meta_clean %>% group_by(Subject_ID) %>% filter(n() == 2) %>% pull(Subject_ID) %>% unique()
  
  if(length(paired_subjects) < 3) next
  
  for(pheno in c("Lean", "Obese")) {
    meta_sub <- col_meta_clean %>% filter(Subject_ID %in% paired_subjects, Ind_Phenotype == pheno)
    if(length(unique(meta_sub$Subject_ID)) >= 3) {
      mat_sub <- mat[, meta_sub$OriginalCol]
      design <- model.matrix(~ factor(Subject_ID) + factor(Time, levels = c("Baseline", "Post3h")), data = meta_sub)
      fit <- eBayes(lmFit(mat_sub, design))
      res <- topTable(fit, coef = ncol(design), number = Inf, sort.by = "none")
      res$Feature <- rownames(res); res$Dataset <- ds_name; res$Phenotype <- pheno
      all_diff_results[[paste(ds_name, pheno, sep="_")]] <- res
      
      n_up <- sum(res$P.Value < 0.05 & res$logFC > 0.263, na.rm = TRUE)
      n_down <- sum(res$P.Value < 0.05 & res$logFC < -0.263, na.rm = TRUE)
      summary_counts[[length(summary_counts)+1]] <- data.frame(
        Dataset = ds_name, Phenotype = pheno, Up = n_up, Down = n_down, Total = nrow(res),
        Pct_Up = (n_up/nrow(res))*100, Pct_Down = (n_down/nrow(res))*100
      )
    }
  }
}

df_res_all <- bind_rows(all_diff_results)
df_summary <- bind_rows(summary_counts)

# ==============================================================================
# 4. PART II: Fig. 2A Global landscape bar plot
# ==============================================================================
cat("\n>>> Plotting Fig. 2A Landscape...")
df_bar <- df_summary %>%
  mutate(Label = paste0(str_split(Dataset, "_", simplify = T)[, 1], "\n(", str_split(Dataset, "_", simplify = T)[, 2], ")"),
         Phenotype = factor(Phenotype, levels = c("Obese", "Lean"))) %>%
  pivot_longer(cols = c(Pct_Up, Pct_Down), names_to = "Direction", values_to = "Percentage") %>%
  mutate(Plot_Value = ifelse(Direction == "Pct_Down", -Percentage, Percentage))

p_landscape <- ggplot(df_bar, aes(x = Label, y = Plot_Value, fill = paste(Phenotype, Direction, sep="_"))) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  facet_wrap(~ Phenotype) + coord_flip() + my_clean_theme +
  scale_fill_manual(values = c("Lean_Pct_Up"="#3DA6AE", "Lean_Pct_Down"="#88C5CA", "Obese_Pct_Up"="#F79647", "Obese_Pct_Down"="#FBC497")) +
  labs(title = "Global Exercise Response Landscape", y = "Percentage of Altered Features (%)")

ggsave(file.path(OUT_DIR, "Fig2A_Response_Landscape.pdf"), p_landscape, width = 11, height = 7)

# ==============================================================================
# 5. PART III: Fig. 2B Serum metabolomics volcano plot
# ==============================================================================
cat("\n>>> Plotting Fig. 2B Serum Volcano...")
df_vol_meta <- df_res_all %>% filter(Dataset == "Serum_Metabonomics") %>%
  mutate(Sig = case_when(P.Value < 0.05 & logFC > 0.263 ~ "Up", P.Value < 0.05 & logFC < -0.263 ~ "Down", TRUE ~ "NS"),
         Color_G = ifelse(Sig == "NS", "NS", paste(Phenotype, Sig, sep="_")), logP = -log10(P.Value),
         Phenotype = factor(Phenotype, levels = c("Obese", "Lean")))

top_labels_meta <- df_vol_meta %>% filter(Sig != "NS") %>% group_by(Phenotype) %>% slice_max(order_by = logP, n = 10)

p_vol_meta <- ggplot(df_vol_meta, aes(x = logFC, y = logP)) +
  geom_point(aes(color = Color_G), alpha = 0.6, size = 1.8) +
  geom_text_repel(data = top_labels_meta, aes(label = Feature), size = 4, fontface = "bold", max.overlaps = 30) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.263, 0.263), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Lean_Up"="#3DA6AE", "Lean_Down"="#88C5CA", "Obese_Up"="#F79647", "Obese_Down"="#FBC497", "NS"="#D3D3D3")) +
  facet_wrap(~ Phenotype) + my_clean_theme +
  labs(title = "Serum Metabonomics: Exercise Response", x = "Log2 Fold Change", y = "-Log10(P-value)")

ggsave(file.path(OUT_DIR, "Fig2B_Serum_Volcano.pdf"), p_vol_meta, width = 12, height = 6)

# ==============================================================================
# 6. PART IV: Baseline comparison (Obese vs Lean)
# ==============================================================================
cat("\n>>> Running Baseline Analysis...")
baseline_results <- list()
for (f in omics_files) {
  ds_name <- str_remove(basename(f), "^Cleaned_") %>% str_remove("\\.csv$")
  if(grepl("Clinical|Sample", ds_name, ignore.case = TRUE)) next
  df <- read_csv(f, show_col_types = FALSE); mat <- as.matrix(df[,-1]); rownames(mat) <- as.character(df[[1]])
  if(max(mat, na.rm = TRUE) > 100) mat <- log2(mat + 1)
  
  col_meta <- data.frame(OriginalCol = colnames(mat), stringsAsFactors = FALSE) %>%
    mutate(Clean_Name = tolower(gsub("_twin\\.[0-9]+", "", OriginalCol)),
           TimeType = case_when(grepl("pre", Clean_Name) ~ "pre", grepl("fast", Clean_Name) ~ "fast", TRUE ~ "Other"),
           Subject_ID = gsub("_fast|_pre|_post1h|_post3h", "", Clean_Name)) %>% inner_join(clin, by = "Subject_ID")
  
  meta_base <- col_meta %>% filter(TimeType %in% c("pre", "fast")) %>% arrange(Subject_ID, TimeType) %>% group_by(Subject_ID) %>% slice(1) %>% ungroup()
  if(length(unique(meta_base$Ind_Phenotype)) == 2 && nrow(meta_base) >= 6) {
    meta_base$Ind_Phenotype <- factor(meta_base$Ind_Phenotype, levels = c("Lean", "Obese"))
    design_base <- model.matrix(~ Ind_Phenotype, data = meta_base)
    fit_base <- eBayes(lmFit(mat[, meta_base$OriginalCol], design_base))
    res_base <- topTable(fit_base, coef = "Ind_PhenotypeObese", number = Inf, sort.by = "none")
    res_base$Feature <- rownames(res_base); res_base$Dataset <- ds_name
    baseline_results[[ds_name]] <- res_base
  }
}
df_base_all <- bind_rows(baseline_results)

# Plot baseline global landscape
df_base_sum <- df_base_all %>% group_by(Dataset) %>%
  summarise(Total = n(), Up_Obese = sum(P.Value < 0.05 & logFC > 0.263, na.rm=T), Up_Lean = sum(P.Value < 0.05 & logFC < -0.263, na.rm=T), .groups="drop") %>%
  mutate(Label = paste0(str_split(Dataset, "_", simplify = T)[,1], "\n(", str_split(Dataset, "_", simplify = T)[,2], ")"),
         Pct_Obese = (Up_Obese/Total)*100, Pct_Lean = (Up_Lean/Total)*100) %>%
  pivot_longer(cols = c(Pct_Obese, Pct_Lean), names_to = "Dir", values_to = "Pct") %>%
  mutate(PlotV = ifelse(Dir == "Pct_Lean", -Pct, Pct))

p_base_landscape <- ggplot(df_base_sum, aes(x = Label, y = PlotV, fill = Dir)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  scale_fill_manual(values = c("Pct_Obese"="#F79647", "Pct_Lean"="#3DA6AE")) +
  coord_flip() + my_clean_theme + labs(title = "Global Baseline Differences (Obese vs. Lean)", x = "", y = "Percentage (%)")

ggsave(file.path(OUT_DIR, "Fig2_Baseline_Global_Landscape.pdf"), p_base_landscape, width = 8, height = 7)

# ==============================================================================
# 7. PART V: Baseline volcano plots by tissue (1:1 replicated advanced rendering)
# ==============================================================================
cat("\n>>> Generating Baseline Tissue Volcanos (ggrastr)...")
df_base_plot <- df_base_all %>%
  mutate(Sig = case_when(P.Value < 0.05 & logFC > 0.263 ~ "Up_in_Obese", P.Value < 0.05 & logFC < -0.263 ~ "Up_in_Lean", TRUE ~ "NS"),
         logP = -log10(P.Value),
         Facet_Name = paste0(str_split(Dataset, "_", simplify = T)[,1], "\n", str_split(Dataset, "_", simplify = T)[,2]),
         Plot_Group = ifelse(Dataset %in% c("Adipose_Microarray", "Serum_Metabonomics", "Muscle_Proteomics"), "Main_Triad", "Supp"))

df_base_plot$Facet_Name <- factor(df_base_plot$Facet_Name, levels = c("Serum\nMetabonomics", "Adipose\nMicroarray","Muscle\nProteomics",
                                                                      "Adipose\nProteomics", "Adipose\nMethylation", "Muscle\nMicroarray", "Muscle\nMethylation", "Serum\nProteomics"))

for(g in c("Main_Triad", "Supp")) {
  df_g <- df_base_plot %>% filter(Plot_Group == g)
  top_g <- df_g %>% filter(Sig != "NS") %>% group_by(Facet_Name, Sig) %>% slice_max(order_by = logP, n = 4)
  
  p_g <- ggplot(df_g, aes(x = logFC, y = logP)) +
    geom_point_rast(aes(color = Sig), alpha = 0.7, size = 1.5, raster.dpi = 300) +
    geom_text_repel(data = top_g, aes(label = Feature), size = 3.8, fontface = "bold", max.overlaps = 20) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-0.263, 0.263), linetype = "dashed", color = "black") +
    scale_color_manual(values = c("Up_in_Obese"="#F79647", "Up_in_Lean"="#3DA6AE", "NS"="#D3D3D3")) +
    facet_wrap(~ Facet_Name, ncol = 3, scales = "free_x") + my_clean_theme +
    labs(title = paste0("Baseline: Obese vs Lean (", g, ")"), x = "Log2 Fold Change", y = "-Log10(P-value)")
  
  ggsave(file.path(OUT_DIR, paste0("Fig2_Baseline_Volcano_", g, ".pdf")), p_g, width = 14, height = 9)
}

# ==============================================================================
# 8. PART VI: Save supplementary data tables separately
# ==============================================================================
cat("\n>>> Exporting Supplementary Tables Separately...")
write.xlsx(list("Summary" = df_summary, "Detailed_Results" = df_res_all), 
           file.path(OUT_DIR, "Data_Exercise_Response_Full.xlsx"), overwrite = TRUE)

write.xlsx(list("Summary" = df_base_sum, "Detailed_Results" = df_base_all), 
           file.path(OUT_DIR, "Data_Baseline_Differences_Full.xlsx"), overwrite = TRUE)

message("\n>>> [Success] Fig. 2 Full Analysis Completed. Results in: ", OUT_DIR)