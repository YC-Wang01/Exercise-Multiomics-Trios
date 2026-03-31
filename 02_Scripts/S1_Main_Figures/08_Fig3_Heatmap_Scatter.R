# ==============================================================================
# Project: FinlandSports (ATM) - Figure 3 Components
# Author: Gemini (Integrated)
# Description: 
#   Strict reproduction of Limma & Heatmap Crosstalk logic.
#   Changes:
#     - Upgraded to V2.0 I/O standards (Cleaned_*.csv & Clinical_Master_Strict)
#     - Removed hardcoded ABC panel labels for flexible patchwork
#     - Integrated Limma & Heatmap Crosstalk into comprehensive data table
#     - Restored readxl for Temp Excel reading and updated deprecated ggplot2 syntax
# ==============================================================================

# ------------------------------------------------------------------------------
# STEP 0: SETUP & LIBRARIES
# ------------------------------------------------------------------------------
# setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readxl, readr, dplyr, stringr, tibble, limma, 
  ggplot2, ggpubr, ggrepel, RColorBrewer, gridExtra, 
  matrixStats, openxlsx, ComplexHeatmap, circlize, showtext,
  reshape2, igraph, ggraph, corrplot, ggsci, cowplot, Hmisc
)

# Font Setup
tryCatch({
  font_add("Arial", "arial.ttf") 
  showtext_auto()
}, error = function(e) message("Font loading failed, falling back to default."))

options(stringsAsFactors = FALSE)

# --- DIRECTORY CONFIGURATION (V2.0 Standard) ---
INPUT_DIR    <- "01_Clean_Data"
META_FILE    <- file.path(INPUT_DIR, "Clinical_Master_Strict.csv")
MAIN_OUT_DIR <- "03_Results/Fig_3/Heatmap_and_Scatter"

if(!dir.exists(MAIN_OUT_DIR)) dir.create(MAIN_OUT_DIR, recursive = TRUE)

# Initialize global list for subsequent Limma table merging
LIMMA_RESULTS_LIST <<- list()

message(">>> Directories Ready. Output set to: ", MAIN_OUT_DIR)

# --- GLOBAL HELPERS ---

# 1. Metadata Loader (V2.0 adaptation)
if(!file.exists(META_FILE)) stop(paste("Metadata file not found at:", META_FILE))
meta <- read_csv(META_FILE, show_col_types = FALSE) %>%
  mutate(
    Group = case_when(`Fat(%)` >= 30 ~ "Obese", !is.na(`Fat(%)`) ~ "Lean", TRUE ~ NA_character_),
    FamilyID = str_trim(as.character(FamilyID)), 
    Membercode = str_trim(as.character(Membercode))
  ) %>% filter(!is.na(Group))

# 2. STRICT PARSER 
parse_sample_info_strict <- function(col_name) {
  x_clean <- str_remove(col_name, "\\.\\.\\.\\d+$") %>% str_trim()
  parts <- str_split(x_clean, "[_\\s]")[[1]]
  if(length(parts) < 2) return(NULL)
  
  fid <- str_trim(parts[1]); p_lower <- tolower(parts)
  mcode <- NA
  if(any(str_detect(p_lower, "father"))) mcode <- "3" # V2.0 Role: Father=3
  else if(any(str_detect(p_lower, "mother"))) mcode <- "2" # Mother=2
  else if(any(str_detect(p_lower, "daughter|son|child"))) mcode <- "1" # Daughter=1
  
  timepoint <- NA
  if(any(str_detect(p_lower, "post3h|post 3h|post_3h"))) timepoint <- "Post3h"
  else if(any(str_detect(p_lower, "pre|fast"))) timepoint <- "Pre"
  
  if(is.na(mcode) || is.na(timepoint)) return(NULL)
  idx <- which(meta$FamilyID == fid & meta$Membercode == mcode)
  if(length(idx) != 1) return(NULL)
  
  return(data.frame(
    OriginalCol = col_name, 
    SubjectID = paste0("S", fid, "_", mcode), 
    Timepoint = timepoint, 
    Group = meta$Group[idx], 
    stringsAsFactors = FALSE
  ))
}

# 3. Metabolite Classifier
classify_metabolite <- function(met_names) {
  sapply(met_names, function(x) {
    lc <- tolower(x)
    if(str_detect(lc, "hdl|ldl|vldl|mufa|pufa|sfa|fa|lipid|chol|tg|mobch")) return("Lipid")
    if(str_detect(lc, "val|leu|ile|ala|gly|tyr|phe|his|gln|glu|aa|amino")) return("Amino Acid")
    if(str_detect(lc, "pyr|lac|cit|bohbut|acace|glc|glol|ketone")) return("Energy")
    if(str_detect(lc, "gp|glyca|c-reactive")) return("Inflammation")
    return("Other")
  })
}

# ------------------------------------------------------------------------------
# STEP 1: LIMMA ANALYSIS
# ------------------------------------------------------------------------------
message("\n=== STEP 1: Differential Analysis (Limma) ===")

# Adapted to V2.0 Cleaned CSVs
targets <- list(
  Mus_Prot = "Cleaned_Muscle_Proteomics.csv",
  Adi_Prot = "Cleaned_Adipose_Proteomics.csv",
  Ser_Meta = "Cleaned_Serum_Metabonomics.csv",
  Ser_Prot = "Cleaned_Serum_Proteomics.csv" 
)

process_limma_dataset <- function(label, filename) {
  fpath <- file.path(INPUT_DIR, filename)
  if(!file.exists(fpath)) { message("  [Skip] File not found: ", filename); return(NULL) }
  
  # Load & Parse
  raw <- read_csv(fpath, show_col_types = FALSE)
  # Safely replace underscores with hyphens in feature names here only
  feat_ids <- str_replace_all(raw[[1]], "_", "-")
  expr_raw <- raw[, -1]
  
  col_infos <- list()
  for(col in colnames(expr_raw)) {
    info <- parse_sample_info_strict(col)
    if(!is.null(info)) col_infos[[col]] <- info
  }
  if(length(col_infos) == 0) return(NULL)
  
  col_info <- bind_rows(col_infos)
  expr_valid <- expr_raw[, col_info$OriginalCol]
  
  expr_num <- suppressWarnings(as.data.frame(lapply(expr_valid, function(x) as.numeric(as.character(x)))))
  rownames(expr_num) <- feat_ids
  sample_meta <- col_info
  
  # Handle Duplicates
  sample_meta$UniqueID <- paste0(sample_meta$SubjectID, "_", sample_meta$Timepoint)
  if(any(duplicated(sample_meta$UniqueID))) {
    t_expr <- t(expr_num)
    agg_res <- aggregate(t_expr, by=list(ID=sample_meta$UniqueID), FUN=mean, na.rm=TRUE)
    mat_final <- t(agg_res[,-1]); colnames(mat_final) <- agg_res$ID
    sample_meta_final <- sample_meta %>% distinct(UniqueID, .keep_all = TRUE)
    rownames(sample_meta_final) <- sample_meta_final$UniqueID
    mat_final <- mat_final[, sample_meta_final$UniqueID]
  } else {
    mat_final <- as.matrix(expr_num); colnames(mat_final) <- sample_meta$UniqueID
    sample_meta_final <- sample_meta
  }
  
  # Normalize
  mat_final <- mat_final[rowMeans(is.na(mat_final)) < 0.5, ]
  min_val <- min(mat_final, na.rm=TRUE); if(!is.finite(min_val)) min_val <- 0
  mat_final[is.na(mat_final)] <- min_val / 2
  if(max(mat_final, na.rm=TRUE) > 50) mat_final <- log2(mat_final + 1)
  
  # Pairing
  target_pre <- "Pre"; target_post <- "Post3h"
  subs <- unique(sample_meta_final$SubjectID)
  paired_subs <- c()
  for(s in subs) {
    times <- sample_meta_final$Timepoint[sample_meta_final$SubjectID == s]
    if(target_pre %in% times && target_post %in% times) paired_subs <- c(paired_subs, s)
  }
  if(length(paired_subs) < 3) return(NULL)
  
  use_ids <- sample_meta_final$UniqueID[sample_meta_final$SubjectID %in% paired_subs & sample_meta_final$Timepoint %in% c(target_pre, target_post)]
  pair_mat <- mat_final[, use_ids]
  pair_meta <- sample_meta_final %>% filter(UniqueID %in% use_ids)
  
  # Limma engine
  run_limma_engine <- function(sub_meta, sub_mat, comp_name) {
    Time <- factor(sub_meta$Timepoint, levels=c(target_pre, target_post))
    Subject <- factor(sub_meta$SubjectID)
    if(nlevels(Time) < 2) return(NULL)
    design <- model.matrix(~0 + Time + Subject)
    colnames(design) <- str_replace_all(colnames(design), "Time", "")
    cm <- makeContrasts(contrasts = paste0(target_post, "-", target_pre), levels=design)
    fit <- tryCatch({ eBayes(contrasts.fit(lmFit(sub_mat, design), cm)) }, error=function(e) NULL)
    if(is.null(fit)) return(NULL)
    res <- topTable(fit, number=Inf) %>% rownames_to_column("FeatureID")
    res$Comparison <- comp_name
    return(res)
  }
  
  res_all <- run_limma_engine(pair_meta, pair_mat, "Overall")
  res_ob  <- run_limma_engine(pair_meta %>% filter(Group=="Obese"), pair_mat[,pair_meta$Group=="Obese"], "Obese")
  res_ln  <- run_limma_engine(pair_meta %>% filter(Group=="Lean"), pair_mat[,pair_meta$Group=="Lean"], "Lean")
  
  full_res <- bind_rows(res_all, res_ob, res_ln)
  
  # Save original differential results and append to comprehensive list
  write.csv(full_res, file.path(MAIN_OUT_DIR, paste0("Limma_Paired_", label, ".csv")), row.names=F)
  LIMMA_RESULTS_LIST[[label]] <<- full_res 
  
  message(paste0("  [Done] Processed: ", label))
}

for(n in names(targets)) process_limma_dataset(n, targets[[n]])

# ------------------------------------------------------------------------------
# STEP 2: DELTA MATRIX CALCULATION
# ------------------------------------------------------------------------------
message("\n=== STEP 2: Preparing Delta Matrices ===")

get_delta_matrix <- function(filename) {
  fpath <- file.path(INPUT_DIR, filename)
  if(!file.exists(fpath)) return(NULL)
  raw <- read_csv(fpath, show_col_types = FALSE)
  feat_ids <- str_replace_all(raw[[1]], "_", "-") 
  expr_raw <- raw[, -1]
  
  col_infos <- list()
  for(col in colnames(expr_raw)) {
    info <- parse_sample_info_strict(col)
    if(!is.null(info)) col_infos[[col]] <- info
  }
  if(length(col_infos) == 0) return(NULL)
  
  sample_meta <- bind_rows(col_infos)
  mat <- suppressWarnings(as.matrix(sapply(expr_raw[, sample_meta$OriginalCol], as.numeric)))
  rownames(mat) <- feat_ids
  if(max(mat, na.rm=T) > 50) mat <- log2(mat + 1)
  
  subjects <- unique(sample_meta$SubjectID)
  delta_list <- list()
  for(s in subjects) {
    pre <- sample_meta %>% dplyr::filter(SubjectID==s, Timepoint=="Pre") %>% pull(OriginalCol)
    post <- sample_meta %>% dplyr::filter(SubjectID==s, Timepoint=="Post3h") %>% pull(OriginalCol)
    if(length(pre)>0 && length(post)>0) {
      val_pre <- if(length(pre)>1) rowMeans(mat[, pre], na.rm=T) else mat[, pre]
      val_post <- if(length(post)>1) rowMeans(mat[, post], na.rm=T) else mat[, post]
      delta_list[[s]] <- val_post - val_pre
    }
  }
  if(length(delta_list)<5) return(NULL)
  do.call(cbind, delta_list)
}

mat_mus <- get_delta_matrix(targets$Mus_Prot)
mat_adi <- get_delta_matrix(targets$Adi_Prot)
mat_ser <- get_delta_matrix(targets$Ser_Meta) 

# ------------------------------------------------------------------------------
# STEP 3: EXTRACT TOP FEATURES
# ------------------------------------------------------------------------------
message("\n=== STEP 3: Extracting Top Features ===")

get_top_features <- function(limma_file, n_top=30) {
  fpath <- file.path(MAIN_OUT_DIR, limma_file) 
  if(!file.exists(fpath)) return(NULL)
  res <- read.csv(fpath)
  res %>% dplyr::filter(Comparison == "Overall") %>% dplyr::arrange(P.Value) %>% head(n_top) %>% pull(FeatureID)
}

# ------------------------------------------------------------------------------
# STEP 4: BLOCKED & LOCKED HEATMAPS
# ------------------------------------------------------------------------------
message("\n=== STEP 4: Blocked Heatmaps ===")

run_heatmap_blocked <- function() {
  if(is.null(mat_mus) || is.null(mat_adi) || is.null(mat_ser)) return()
  
  t_m <- get_top_features("Limma_Paired_Mus_Prot.csv", 18)
  t_a <- get_top_features("Limma_Paired_Adi_Prot.csv", 18)
  t_s <- get_top_features("Limma_Paired_Ser_Meta.csv", 20)
  
  M <- mat_mus[intersect(rownames(mat_mus), t_m), , drop=FALSE]
  A <- mat_adi[intersect(rownames(mat_adi), t_a), , drop=FALSE]
  S <- mat_ser[intersect(rownames(mat_ser), t_s), , drop=FALSE]
  
  common_subs <- intersect(intersect(colnames(M), colnames(A)), colnames(S))
  
  # Helper to calculate R and P
  calc_stats <- function(mT, mS) {
    res <- rcorr(as.matrix(cbind(t(mT), t(mS))), type="spearman")
    nT <- nrow(mT); nS <- nrow(mS)
    r_out <- res$r[1:nT, (nT+1):(nT+nS), drop=FALSE]
    p_out <- res$P[1:nT, (nT+1):(nT+nS), drop=FALSE]
    r_out[is.na(r_out)] <- 0 
    p_out[is.na(p_out)] <- 1 
    list(r=r_out, p=p_out)
  }
  
  # --- OBESE (Reference) ---
  subs_ob <- meta %>% mutate(SID=paste0("S",FamilyID,"_",Membercode)) %>% 
    filter(Group=="Obese", SID %in% common_subs) %>% pull(SID)
  
  stat_Mob <- calc_stats(M[,subs_ob], S[,subs_ob])
  stat_Aob <- calc_stats(A[,subs_ob], S[,subs_ob])
  
  # Ordering
  serum_cls <- classify_metabolite(colnames(stat_Mob$r))
  serum_cls <- factor(serum_cls, levels = c("Energy", "Amino Acid", "Lipid", "Inflammation", "Other"))
  ord_A <- hclust(dist(stat_Aob$r))$order
  ord_M <- hclust(dist(stat_Mob$r))$order
  
  # Sort
  sort_mat <- function(m, ord_r) m[ord_r, ]
  r_Mob <- sort_mat(stat_Mob$r, ord_M); p_Mob <- sort_mat(stat_Mob$p, ord_M)
  r_Aob <- sort_mat(stat_Aob$r, ord_A); p_Aob <- sort_mat(stat_Aob$p, ord_A)
  
  # --- LEAN (Locked) ---
  subs_ln <- meta %>% mutate(SID=paste0("S",FamilyID,"_",Membercode)) %>% 
    filter(Group=="Lean", SID %in% common_subs) %>% pull(SID)
  stat_Mln <- calc_stats(M[,subs_ln], S[,subs_ln])
  stat_Aln <- calc_stats(A[,subs_ln], S[,subs_ln])
  
  # Lock Order
  r_Mln <- stat_Mln$r[rownames(r_Mob), ]; p_Mln <- stat_Mln$p[rownames(r_Mob), ]
  r_Aln <- stat_Aln$r[rownames(r_Aob), ]; p_Aln <- stat_Aln$p[rownames(r_Aob), ]
  
  # Draw
  draw_it <- function(rM, rA, title, fname) {
    col_fun = colorRamp2(c(-0.6, 0, 0.6), c("#2166AC", "#FFFFFF", "#B2182B"))
    ha <- HeatmapAnnotation(Class=serum_cls, 
                            col=list(Class=c("Energy"="#F39B7F", "Amino Acid"="#00A087", "Lipid"="#3C5488", "Inflammation"="#DC0000", "Other"="grey")),
                            show_legend = TRUE, show_annotation_name = FALSE, simple_anno_size = unit(3, "mm"))
    
    ht1 <- Heatmap(rA, name="Corr (Adi)", column_title=title, col=col_fun, top_annotation=ha,
                   column_split=serum_cls, cluster_columns=T, cluster_column_slices=F, cluster_rows=F,
                   show_column_names = FALSE,
                   row_names_gp = gpar(fontsize = 9, fontfamily="Arial"),
                   column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                   rect_gp = gpar(col = "white", lwd = 0.5),
                   row_title = "Adipose", row_title_gp = gpar(col = "#B2182B", fontface="bold"))
    
    ht2 <- Heatmap(rM, name="Corr (Mus)", col=col_fun, column_split=serum_cls, 
                   cluster_columns=T, cluster_column_slices=F, cluster_rows=F,
                   show_column_names = TRUE, 
                   column_names_gp = gpar(fontsize = 9, fontfamily="Arial", fontface="bold"), 
                   column_names_rot = 45, 
                   row_names_gp = gpar(fontsize = 9, fontfamily="Arial"),
                   rect_gp = gpar(col = "white", lwd = 0.5),
                   row_title = "Muscle", row_title_gp = gpar(col = "#2166AC", fontface="bold"))
    
    pdf(file.path(MAIN_OUT_DIR, fname), width=8.5, height=10)
    draw(ht1 %v% ht2, ht_gap = unit(2, "mm"))
    dev.off()
  }
  
  draw_it(r_Mob, r_Aob, "Obese: Tissue-Serum Crosstalk (Reference)", "Plot_Heatmap_Crosstalk_Obese.pdf")
  draw_it(r_Mln, r_Aln, "Lean: Tissue-Serum Crosstalk (Locked)", "Plot_Heatmap_Crosstalk_Lean.pdf")
  
  # SAVE EXCEL (Temp files for aggregation)
  save_wb <- function(rM, pM, rA, pA, fname) {
    wb <- createWorkbook()
    addWorksheet(wb, "Muscle_R"); writeData(wb, "Muscle_R", as.data.frame(rM), rowNames=T)
    addWorksheet(wb, "Muscle_P"); writeData(wb, "Muscle_P", as.data.frame(pM), rowNames=T)
    addWorksheet(wb, "Adipose_R"); writeData(wb, "Adipose_R", as.data.frame(rA), rowNames=T)
    addWorksheet(wb, "Adipose_P"); writeData(wb, "Adipose_P", as.data.frame(pA), rowNames=T)
    saveWorkbook(wb, file.path(MAIN_OUT_DIR, fname), overwrite=T)
  }
  save_wb(r_Mob, p_Mob, r_Aob, p_Aob, "Temp_Heatmap_Obese.xlsx")
  save_wb(r_Mln, p_Mln, r_Aln, p_Aln, "Temp_Heatmap_Lean.xlsx")
}
run_heatmap_blocked()

# ------------------------------------------------------------------------------
# STEP 5: CORE GENE SCATTER PLOTS
# ------------------------------------------------------------------------------
message("\n=== STEP 5: Core Gene Scatter Plots ===")

generate_core_plot <- function(met, adi_prot, mus_prot, suffix) {
  if(!met %in% rownames(mat_ser)) return()
  common_subs <- intersect(intersect(colnames(mat_adi), colnames(mat_mus)), colnames(mat_ser))
  
  df_plot <- data.frame(
    SID = common_subs,
    Serum = as.numeric(mat_ser[met, common_subs]),
    Adipose = as.numeric(mat_adi[adi_prot, common_subs]),
    Muscle = as.numeric(mat_mus[mus_prot, common_subs])
  ) %>%
    left_join(meta %>% mutate(SID=paste0("S",FamilyID,"_",Membercode)), by="SID") %>%
    filter(!is.na(Group))
  
  draw_p <- function(y_col, y_id) {
    ggplot(df_plot, aes(x = Serum, y = .data[[y_col]], color = Group, fill = Group)) +
      geom_smooth(method="lm", se=F, fullrange=T) +
      geom_point(shape=21, color="white", size=3) +
      stat_cor(method="pearson", label.y.npc="top") +
      scale_color_manual(values=c("Obese"="#E64B35", "Lean"="#4DBBD5")) +
      scale_fill_manual(values=c("Obese"="#E64B35", "Lean"="#4DBBD5")) +
      theme_pubr() + 
      theme(panel.grid = element_blank()) + # Ensure no grid lines
      labs(x=paste0("Serum ", met, " (\u0394)"), y=paste0(y_id, " (\u0394)"))
  }
  
  p <- ggarrange(draw_p("Adipose", adi_prot), draw_p("Muscle", mus_prot), ncol=2, common.legend=T)
  ggsave(file.path(MAIN_OUT_DIR, paste0("Plot_Scatter_", suffix, "_", met, ".pdf")), p, width=10, height=5)
}

generate_core_plot("Phe", "IRS1", "SMIM1", "1_IR")
generate_core_plot("AcAce", "RENBP", "LTF", "2_Stress")
generate_core_plot("Gly", "LSM12", "S100A9", "3_AA")

# ------------------------------------------------------------------------------
# STEP 6: COMPREHENSIVE SUPPLEMENTARY TABLE
# ------------------------------------------------------------------------------
message("\n=== STEP 6: Generating Comprehensive Supplementary Table ===")

generate_sup_table <- function() {
  f_ob <- file.path(MAIN_OUT_DIR, "Temp_Heatmap_Obese.xlsx")
  f_ln <- file.path(MAIN_OUT_DIR, "Temp_Heatmap_Lean.xlsx")
  
  if(!file.exists(f_ob) || !file.exists(f_ln)) { message("  [Skip] Source Data missing."); return() }
  
  read_melt <- function(f, sheet_r, sheet_p, grp) {
    mat_r <- read_excel(f, sheet=sheet_r); p_ord <- mat_r[[1]]; colnames(mat_r)[1]<-"Protein_ID"
    mat_p <- read_excel(f, sheet=sheet_p); colnames(mat_p)[1]<-"Protein_ID"
    df_r <- melt(mat_r, id="Protein_ID", variable.name="Met", value.name="R")
    df_p <- melt(mat_p, id="Protein_ID", variable.name="Met", value.name="P")
    df_r %>% mutate(P=df_p$P, Protein_ID=factor(Protein_ID, levels=p_ord)) %>% rename_with(~paste0(grp, "_", .), c("R","P"))
  }
  
  proc_tissue <- function(lbl) {
    df <- merge(read_melt(f_ob, paste0(lbl,"_R"), paste0(lbl,"_P"), "Obese"),
                read_melt(f_ln, paste0(lbl,"_R"), paste0(lbl,"_P"), "Lean"), by=c("Protein_ID","Met"))
    df$Tissue <- lbl
    df$Class <- classify_metabolite(as.character(df$Met))
    df %>% filter(Obese_P < 0.05 | Lean_P < 0.05) %>%
      arrange(Protein_ID, factor(Class, levels=c("Energy","Amino Acid","Lipid","Inflammation","Other")))
  }
  
  final_tab <- bind_rows(proc_tissue("Adipose"), proc_tissue("Muscle")) %>%
    mutate(Obese_Sig = ifelse(Obese_P<0.05, "*", ""), Lean_Sig = ifelse(Lean_P<0.05, "*", "")) %>%
    select(Tissue, Class, Protein_ID, Met, Obese_R, Obese_P, Obese_Sig, Lean_R, Lean_P, Lean_Sig)
  
  wb <- createWorkbook()
  
  # 1. Write heatmap mechanism data with conditional formatting
  addWorksheet(wb, "Heatmap_Crosstalk")
  writeData(wb, "Heatmap_Crosstalk", final_tab)
  
  sig_red <- createStyle(fontColour="#B2182B", textDecoration="bold")
  sig_blue <- createStyle(fontColour="#2166AC", textDecoration="bold")
  conditionalFormatting(wb, "Heatmap_Crosstalk", cols=5, rows=2:(nrow(final_tab)+1), rule=">0", style=sig_red)
  conditionalFormatting(wb, "Heatmap_Crosstalk", cols=5, rows=2:(nrow(final_tab)+1), rule="<0", style=sig_blue)
  conditionalFormatting(wb, "Heatmap_Crosstalk", cols=8, rows=2:(nrow(final_tab)+1), rule=">0", style=sig_red)
  conditionalFormatting(wb, "Heatmap_Crosstalk", cols=8, rows=2:(nrow(final_tab)+1), rule="<0", style=sig_blue)
  
  # 2. Append Limma differential data
  for(lbl in names(LIMMA_RESULTS_LIST)) {
    sheet_name <- paste0("Limma_", lbl)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, LIMMA_RESULTS_LIST[[lbl]])
  }
  
  saveWorkbook(wb, file.path(MAIN_OUT_DIR, "Data_Comprehensive_Limma_Crosstalk.xlsx"), overwrite=TRUE)
  
  # Clean up: Delete intermediate Excel files generated for heatmap merging to keep the folder clean
  unlink(f_ob)
  unlink(f_ln)
}

generate_sup_table()

message("\n>>> PIPELINE COMPLETE. All results saved cleanly in: ", MAIN_OUT_DIR)
