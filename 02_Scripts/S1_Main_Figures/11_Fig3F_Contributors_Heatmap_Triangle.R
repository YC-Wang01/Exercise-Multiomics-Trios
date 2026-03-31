# ==============================================================================
# Project: FinlandSports (ATM) - Fig 3 Contributors (Heatmap + Triangle)
# Description: 
#   1. Aligned output directory to V2.0 standard: 03_Results/Fig_3/Fig3_Contributors
#   2. Kept Square Heatmaps (with white borders) and Triangle Plots original.
#   3. Merged Sup Table logic: Generates Data_S6 with 8 independent Sheets
#      (4 for Heatmap every-point data, 4 for Triangle Full Set data).
#   4. Label formatting: Cleaned gene names (.) and metadata names (E%, PA).
#   5. Enforced dplyr:: prefix to prevent namespace masking.
# ==============================================================================

# --- 0. Set Working Directory ---
# setwd("C:/Users/Sorcier_W/Desktop/ATM/Exercise-Multiomics-Trios")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readxl, dplyr, stringr, tibble, tidyr, openxlsx, 
  limma, lme4, glmnet, pheatmap, ggplot2, ggpubr, grid, scales, ggrepel,
  ggsci, RColorBrewer, Hmisc, cowplot, gridExtra
)

# --- 1. Parameter Settings ---
TOP_N <- 5       # 5 features * 3 categories = 15 rows
CELL_SIZE <- 12  # Fixed square cell size

# --- 2. Unified Path Configuration ---
ROOT_DIR  <- "C:/Users/Sorcier_W/Desktop/ATM"
DIR_IN    <- file.path(ROOT_DIR, "Preprocessed_Data/Origindata_FinalEdition")
DIR_PHENO <- file.path(ROOT_DIR, "Preprocessed_Data/Processed_Volunteer_Data_HistoryMerge.xlsx")

# Updated to V2.0 relative output path
DIR_OUT   <- "03_Results/Fig_3/Fig3_Contributors"

if(!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive=TRUE)

# Global lists to collect per-point data for heatmaps and full-set data for ternary plots
HEATMAP_RESULTS_LIST  <- list()
TRIANGLE_RESULTS_LIST <- list()
HEATMAP_ORDER_LIST    <- list()

message(">>> Directories Ready. Output set to: ", DIR_OUT)

# --- 3. Core Functions ---

# Dedicated beautification functions for cleaning row and column names
clean_rownames_fn <- function(x) {
  # Remove dots and subsequent characters from gene names (e.g., RBMB4.RBM4 -> RBMB4)
  gsub("\\..*", "", x)
}

clean_colnames_fn <- function(x) {
  # 1. Physical Activity
  x <- gsub("PA_Hour\\.Week", "PA(Duration)", x)
  x <- gsub("PA_Times\\.Week", "PA(Frequence)", x)
  x <- gsub("IPA", "IPA", x) 
  
  # 2. Clinical and Adiposity Metrics
  x <- gsub("Fat_percent", "Fat%", x)
  x <- gsub("WAISTLINE", "WC", x)
  x <- gsub("HOMA_IR", "HOMA-IR", x)
  x <- gsub("Liverfat", "HFC", x, ignore.case = TRUE)
  
  # 3. Dietary Characteristics (Strip E and rename)
  x <- gsub("Energy_kcal", "IEE", x)
  x <- gsub("Fat\\.E\\.\\.", "Fats", x)
  x <- gsub("Protein\\.E\\.\\.", "Protein", x)
  x <- gsub("Carbohydrate\\.E\\.\\.", "Carbohydrate", x)
  x <- gsub("Sucrose\\.E\\.\\.", "Sucrose", x)
  x <- gsub("Alcohol\\.E\\.\\.", "Alcohol", x)
  
  return(x)
}

coord_ternary <- function(a,b,c) { 
  s <- a + b + c; s[s==0] <- 1
  data.frame(x = 0.5 * (2 * b/s + c/s), y = (sqrt(3)/2) * c/s) 
}

calc_lm_stats <- function(df, y_col, x_cols) {
  valid_x <- intersect(x_cols, colnames(df))
  if(length(valid_x) == 0) return(c(0, 1))
  
  dat_sub <- df[, c(y_col, valid_x)]
  dat_sub <- na.omit(dat_sub)
  
  if(nrow(dat_sub) < 5) return(c(0, 1))
  y_vec <- dat_sub[[y_col]]
  x_mat <- as.matrix(dat_sub[, valid_x, drop=FALSE])
  
  if(sd(y_vec) == 0 || any(apply(x_mat, 2, sd) == 0)) return(c(0, 1))
  
  tryCatch({
    fit <- lm(y_vec ~ x_mat)
    s <- summary(fit)
    f <- s$fstatistic
    r2 <- s$r.squared
    p_val <- if (is.null(f)) 1 else pf(f[1], f[2], f[3], lower.tail = FALSE)
    return(c(as.numeric(r2), as.numeric(p_val)))
  }, error = function(e) c(0, 1))
}

get_responders_pool <- function(raw_data, label) {
  orig_cols <- colnames(raw_data); meta_list <- list()
  for(i in 2:length(orig_cols)) {
    match <- str_match(orig_cols[i], "(?i)(\\d+)_([A-Za-z]+)_([A-Za-z0-9]+)")
    if(!is.na(match[1,1])) {
      time_std <- if(tolower(match[1,4]) %in% c("pre","fast")) "Pre" else if(str_detect(tolower(match[1,4]), "post")) "Post" else NA
      if(!is.na(time_std)) meta_list[[length(meta_list)+1]] <- data.frame(SampleID=orig_cols[i], SubjectID=paste0(match[1,2],"_",str_to_title(match[1,3])), Timepoint=time_std, stringsAsFactors=F)
    }
  }
  if(length(meta_list)==0) return(NULL)
  meta <- bind_rows(meta_list)
  paired <- meta %>% dplyr::group_by(SubjectID) %>% dplyr::filter(n()>=2 & all(c("Pre","Post")%in%Timepoint)) %>% dplyr::pull(SubjectID) %>% unique()
  if(length(paired)<3) return(NULL)
  meta_sub <- meta %>% dplyr::filter(SubjectID %in% paired) %>% dplyr::group_by(SubjectID, Timepoint) %>% dplyr::slice(1) %>% dplyr::ungroup() %>% dplyr::arrange(SubjectID, Timepoint)
  
  row_ids <- make.names(as.character(raw_data[[1]]), unique=TRUE)
  mat <- raw_data[, meta_sub$SampleID]; mat <- suppressWarnings(apply(as.matrix(mat), 2, as.numeric)); rownames(mat) <- row_ids
  mat <- mat[rowSums(is.na(mat))<(ncol(mat)*0.8) & apply(mat,1,sd,na.rm=T)>0, ] 
  if(nrow(mat)<5) return(NULL)
  if(max(mat, na.rm=T)>20) mat <- log2(mat+1)
  
  design <- model.matrix(~ 0 + Timepoint + SubjectID, data=meta_sub)
  fit <- lmFit(mat, design); cm <- makeContrasts(diff=TimepointPost-TimepointPre, levels=design)
  res <- topTable(eBayes(contrasts.fit(fit, cm)), coef="diff", number=Inf) %>% tibble::rownames_to_column("Feature")
  
  sig_feats <- res %>% dplyr::filter(P.Value < 0.05) %>% dplyr::pull(Feature)
  candidate_feats <- if(length(sig_feats) < (TOP_N*3)) head(res$Feature[order(res$P.Value)], 50) else sig_feats
  delta <- mat[candidate_feats, meta_sub$Timepoint=="Post"] - mat[candidate_feats, meta_sub$Timepoint=="Pre"]
  colnames(delta) <- meta_sub$SubjectID[meta_sub$Timepoint=="Post"]
  return(list(delta=as.data.frame(t(delta), check.names=F), sig_res=res, feats=candidate_feats, source=label))
}

# --- 4. Clinical Metadata Loading ---
message("\n>>> Loading Phenotypes...")
map_id <- function(df) { df$SubjectID <- paste0(str_split(df$Sample_Key,"_",simplify=T)[,1],"_",dplyr::case_when(str_split(df$Sample_Key,"_",simplify=T)[,2]=="1"~"Daughter",str_split(df$Sample_Key,"_",simplify=T)[,2]=="2"~"Mother",T~"Father")); df$SubjectID<-str_trim(df$SubjectID); return(df) }

lifestyle <- map_id(read_excel(DIR_PHENO, sheet="Lifestyle"))
colnames(lifestyle) <- make.names(colnames(lifestyle), unique=TRUE)
basic <- map_id(read_excel(DIR_PHENO, sheet="Basic information"))
colnames(basic) <- make.names(colnames(basic), unique=TRUE)
endocrine <- map_id(read_excel(DIR_PHENO, sheet="Blood_Endocrine"))
colnames(endocrine) <- make.names(colnames(endocrine), unique=TRUE)
adiposity <- map_id(read_excel(DIR_PHENO, sheet="Ectopic adiposity"))
colnames(adiposity) <- make.names(colnames(adiposity), unique=TRUE)

diet_raw <- read_excel(list.files(DIR_IN, pattern="Diet_Checked", full.names=T)[1])
diet_cat <- as.data.frame(t(diet_raw[,-1]), check.names=F)
colnames(diet_cat) <- make.names(diet_raw[[1]], unique=TRUE)
diet_cat$SubjectID <- str_trim(str_remove(rownames(diet_cat), "_Pre$"))

vars_act <- c("PA_Hour.Week", "PA_Times.Week", "IPA")
vars_mac <- c("Energy_kcal", "Protein.E..", "Fat.E..", "Carbohydrate.E..", "Sucrose.E..", "Alcohol.E..")
vars_cli <- c("BMI", "Fat_percent", "WAISTLINE", "HOMA_IR", "Liverfat") 
vars_fod <- intersect(colnames(diet_cat), c("Beef", "Pork", "Sausages", "Fish", "Eggs", "Milk", "Sour_milk", "Cheese", "Fruits", "Sugar", "Coffee", "Tea"))

pheno_master <- lifestyle %>% dplyr::select(SubjectID, dplyr::any_of(vars_act), dplyr::any_of(vars_mac)) %>%
  dplyr::left_join(basic %>% dplyr::select(SubjectID, dplyr::any_of(vars_cli)), by="SubjectID") %>%
  dplyr::left_join(endocrine %>% dplyr::select(SubjectID, HOMA_IR), by="SubjectID") %>%
  dplyr::left_join(adiposity %>% dplyr::select(SubjectID, Liverfat), by="SubjectID") %>%
  dplyr::left_join(diet_cat %>% dplyr::select(SubjectID, dplyr::any_of(vars_fod)), by="SubjectID") %>% dplyr::distinct(SubjectID, .keep_all=T)

# --- 5. Main Calculation Loop ---
groups <- list(Serum_Metabonomics=c("Serum_Metabonomics"), Serum_Proteomics=c("Serum_Proteomics"), Muscle=c("Muscle_Proteomics"), Adipose=c("Adipose_Proteomics"))

for(g_name in names(groups)) {
  message(">>> Processing Group: ", g_name)
  
  # A. Prep Data
  group_delta <- data.frame(); group_meta <- data.frame()
  for(pat in groups[[g_name]]) {
    f_path <- list.files(DIR_IN, pattern=pat, full.names=T, ignore.case=T)[1]; if(is.na(f_path)) next
    sub_label <- if(str_detect(pat,"Metabo")) "Serum_Metabo" else if(str_detect(pat,"Serum_Prot")) "Serum_Prot" else g_name
    res_obj <- get_responders_pool(suppressMessages(read_excel(f_path)), sub_label); if(is.null(res_obj)) next
    group_meta <- rbind(group_meta, data.frame(Feature=colnames(res_obj$delta), Source=sub_label))
    if(nrow(group_delta)==0) group_delta <- res_obj$delta else group_delta <- merge(group_delta, res_obj$delta, by="row.names", all=T) %>% tibble::column_to_rownames("Row.names")
  }
  if(nrow(group_delta)==0) { message("   [Skip] No delta data."); next }
  
  master_t <- group_delta %>% tibble::rownames_to_column("SubjectID") %>% dplyr::inner_join(pheno_master, by="SubjectID")
  master_t$FamilyID <- str_split(master_t$SubjectID, "_", simplify=T)[,1]
  
  # B. Feature Selection
  calc_max_cor <- function(feat, var_list) {
    if(!feat %in% colnames(master_t)) return(0)
    valid_vars <- intersect(var_list, colnames(master_t)); if(length(valid_vars) == 0) return(0)
    sub_dat <- master_t %>% dplyr::select(dplyr::all_of(feat), dplyr::all_of(valid_vars)) %>% na.omit()
    if(nrow(sub_dat) < 5 || ncol(sub_dat) < 2) return(0)
    if(sd(sub_dat[[1]]) == 0) return(0)
    tryCatch({ max(abs(cor(sub_dat[[1]], sub_dat[,-1], method="spearman")), na.rm=T) }, error=function(e) 0)
  }
  
  all_feats <- intersect(colnames(group_delta), colnames(master_t))
  scores <- data.frame(Feature=all_feats, Score_Act=0, Score_Diet=0, Score_Clin=0)
  for(i in 1:nrow(scores)) {
    f <- scores$Feature[i]
    scores$Score_Act[i] <- calc_max_cor(f, vars_act)
    scores$Score_Diet[i] <- calc_max_cor(f, vars_mac)
    scores$Score_Clin[i] <- calc_max_cor(f, vars_cli)
  }
  
  top_act <- scores %>% dplyr::arrange(desc(Score_Act)) %>% head(TOP_N) %>% dplyr::pull(Feature)
  top_diet <- scores %>% dplyr::filter(!Feature %in% top_act) %>% dplyr::arrange(desc(Score_Diet)) %>% head(TOP_N) %>% dplyr::pull(Feature)
  top_clin <- scores %>% dplyr::filter(!Feature %in% c(top_act, top_diet)) %>% dplyr::arrange(desc(Score_Clin)) %>% head(TOP_N) %>% dplyr::pull(Feature)
  main_feats <- unique(c(top_act, top_diet, top_clin))
  if(length(main_feats) < 3) next
  
  # --- 1. Plot HEATMAP and extract data for each point ---
  cols_main <- unique(c(vars_act, vars_cli, vars_mac))
  valid_cols_main <- intersect(cols_main, colnames(master_t))
  
  if(length(valid_cols_main) > 0) {
    cor_res <- rcorr(as.matrix(master_t[, c(main_feats, valid_cols_main)]), type="spearman")
    cor_mat <- cor_res$r[main_feats, valid_cols_main, drop=F]; cor_mat[is.na(cor_mat)] <- 0
    p_mat <- cor_res$P[main_feats, valid_cols_main, drop=F]; p_mat[is.na(p_mat)] <- 1
    
    if(ncol(cor_mat) > 0) {
      colnames(cor_mat) <- make.names(colnames(cor_mat), unique=T); colnames(p_mat) <- colnames(cor_mat)
      
      col_anno <- data.frame(Category = rep("Other", ncol(cor_mat)), stringsAsFactors=FALSE)
      rownames(col_anno) <- colnames(cor_mat)
      for(nm in rownames(col_anno)) {
        if(any(sapply(make.names(vars_act), function(x) x == nm))) col_anno[nm, "Category"] <- "PA"
        if(any(sapply(make.names(vars_mac), function(x) x == nm))) col_anno[nm, "Category"] <- "Dietary intake"
        if(any(sapply(make.names(vars_cli), function(x) x == nm))) col_anno[nm, "Category"] <- "Adiposity"
      }
      # Map commander-specified color scheme
      anno_colors <- list(Category=c(`PA`="#FFD700", `Adiposity`="#4DBBD5", `Dietary intake`="#00A087"))
      star_mat <- matrix("", nrow=nrow(p_mat), ncol=ncol(p_mat)); star_mat[p_mat<0.05]<-"*"; star_mat[p_mat<0.01]<-"**"
      
      # [Core Cleaning]: Apply beautification functions to clean row and column names
      rownames(cor_mat) <- clean_rownames_fn(rownames(cor_mat))
      rownames(p_mat)   <- clean_rownames_fn(rownames(p_mat))
      
      colnames(cor_mat)  <- clean_colnames_fn(colnames(cor_mat))
      colnames(p_mat)    <- clean_colnames_fn(colnames(p_mat))
      rownames(col_anno) <- clean_colnames_fn(rownames(col_anno))
      
      # [Extract Full Heatmap Data]: Convert to long format for supplementary table
      df_r <- as.data.frame(as.table(cor_mat))
      colnames(df_r) <- c("Feature", "Metadata", "Spearman_R")
      df_p <- as.data.frame(as.table(p_mat))
      colnames(df_p) <- c("Feature", "Metadata", "P_Value")
      
      # =========================================================
      # Core Modification: Draw heatmap first, generate object p
      # =========================================================
      pdf_w <- (ncol(cor_mat) * CELL_SIZE / 72) + 3.5
      pdf_h <- (nrow(cor_mat) * CELL_SIZE / 72) + 2.5
      
      pdf(file.path(DIR_OUT, paste0(g_name, "_Heatmap.pdf")), width=pdf_w, height=pdf_h)
      p <- pheatmap(cor_mat, main=paste(g_name), 
                    annotation_col=col_anno, annotation_colors=anno_colors, annotation_row=NULL, 
                    color=colorRampPalette(c("#3C5488","white","#E64B35"))(100),
                    cluster_rows=TRUE, cluster_cols=FALSE, treeheight_row = 0, 
                    display_numbers=star_mat, number_color="black", fontsize_number=10, 
                    border_color="white", 
                    cellwidth = CELL_SIZE, cellheight = CELL_SIZE, 
                    gaps_col=cumsum(table(col_anno$Category)[unique(col_anno$Category)]), 
                    legend=FALSE,               # <--- Set to FALSE
                    annotation_legend=FALSE,    # <--- Set to FALSE
                    fontsize=8, angle_col=45)
      dev.off() # Close device
      
      # =========================================================
      # Then use object p to extract order and generate curated df_hm table
      # =========================================================
      hm_order <- rownames(cor_mat)[p$tree_row$order]
      HEATMAP_ORDER_LIST[[g_name]] <- hm_order
      
      df_hm <- df_r %>%
        dplyr::left_join(df_p, by = c("Feature", "Metadata")) %>%
        dplyr::mutate(
          Feature = factor(Feature, levels = hm_order), # Convert to factor to fix order
          Metadata = as.character(Metadata),
          Significance = dplyr::case_when(P_Value < 0.001 ~ "***", P_Value < 0.01 ~ "**", P_Value < 0.05 ~ "*", TRUE ~ "ns")
        ) %>%
        dplyr::arrange(Feature) # This ensures the table lists genes from the top of the heatmap downwards
      
      # Load into heatmap-specific global list
      HEATMAP_RESULTS_LIST[[g_name]] <- df_hm
      
      # Draw heatmap again to save order internally
      pdf_w <- (ncol(cor_mat) * CELL_SIZE / 72) + 3.5
      pdf_h <- (nrow(cor_mat) * CELL_SIZE / 72) + 2.5
      
      pdf(file.path(DIR_OUT, paste0(g_name, "_Heatmap.pdf")), width=pdf_w, height=pdf_h)
      p <- pheatmap(cor_mat, main=paste(g_name), 
                    annotation_col=col_anno, annotation_colors=anno_colors, annotation_row=NULL, 
                    color=colorRampPalette(c("#3C5488","white","#E64B35"))(100),
                    cluster_rows=TRUE, cluster_cols=FALSE, treeheight_row = 0, 
                    display_numbers=star_mat, number_color="black", fontsize_number=10, 
                    border_color="white", 
                    cellwidth = CELL_SIZE, cellheight = CELL_SIZE, 
                    gaps_col=cumsum(table(col_anno$Category)[unique(col_anno$Category)]), 
                    legend=FALSE,               # <--- Set to FALSE
                    annotation_legend=FALSE,    # <--- Set to FALSE
                    fontsize=8, angle_col=45)
      HEATMAP_ORDER_LIST[[g_name]] <- rownames(cor_mat)[p$tree_row$order]
      
      dev.off()
    }
  }
  
  # --- 2. Calculate, Save, and Load TRIANGLE Data ---
  g_tri <- data.frame()
  for(feat in all_feats) {
    act_s <- calc_lm_stats(master_t, feat, vars_act)
    diet_s <- calc_lm_stats(master_t, feat, vars_mac)
    clin_s <- calc_lm_stats(master_t, feat, vars_cli)
    
    if(act_s[1] + diet_s[1] + clin_s[1] > 0) {
      g_tri <- rbind(g_tri, data.frame(
        Feature=feat, 
        Activity=act_s[1], Macros=diet_s[1], Clinical=clin_s[1],
        Activity_P=act_s[2], Macros_P=diet_s[2], Clinical_P=clin_s[2]
      ))
    }
  }
  
  if(nrow(g_tri)>0) {
    TRIANGLE_RESULTS_LIST[[g_name]] <- g_tri
    
    coords <- coord_ternary(g_tri$Activity, g_tri$Macros, g_tri$Clinical); coords[is.na(coords)]<-0
    g_tri$Dominant <- apply(g_tri[,c("Activity","Macros","Clinical")], 1, function(x) if(sum(x)==0) "None" else names(x)[which.max(x)])
    g_tri$TotalR2 <- rowSums(g_tri[,c("Activity","Macros","Clinical")]); 
    plot_data <- cbind(g_tri, coords) %>% dplyr::filter(Dominant!="None")
    
    plot_data$dist_act <- sqrt((plot_data$x - 0)^2 + (plot_data$y - 0)^2)              
    plot_data$dist_diet <- sqrt((plot_data$x - 1)^2 + (plot_data$y - 0)^2)            
    plot_data$dist_clin <- sqrt((plot_data$x - 0.5)^2 + (plot_data$y - sqrt(3)/2)^2) 
    
    ids_act <- plot_data %>% dplyr::arrange(dist_act) %>% head(3) %>% dplyr::pull(Feature)
    ids_diet <- plot_data %>% dplyr::arrange(dist_diet) %>% head(3) %>% dplyr::pull(Feature)
    ids_clin <- plot_data %>% dplyr::arrange(dist_clin) %>% head(3) %>% dplyr::pull(Feature)
    all_highlight_ids <- unique(c(main_feats, ids_act, ids_diet, ids_clin))
    
    # Label cleaning: Use clean labels on ternary plots as well
    plot_data$CleanFeature <- clean_rownames_fn(plot_data$Feature)
    
    plot_data$IsHighlight <- plot_data$Feature %in% all_highlight_ids
    plot_data$Dominant <- dplyr::recode(plot_data$Dominant, "Activity"="PA", "Macros"="Dietary intake", "Clinical"="Adiposity")
    
    tri_bg <- data.frame(x=c(0, 1, 0.5, 0), y=c(0, 0, sqrt(3)/2, 0))
    # Modify vertex labels of the ternary plot
    tri_labels <- data.frame(x = c(0, 1, 0.5), y = c(-0.05, -0.05, sqrt(3)/2 + 0.05), 
                             label = c("PA", "Dietary intake", "Adiposity"), 
                             hjust = c(0.5, 0.5, 0.5), vjust = c(1, 1, 0))
    
    p_tri <- ggplot() +
      geom_polygon(data=tri_bg, aes(x, y), fill="white", color="black", linewidth=0.8) +
      geom_text(data=tri_labels, aes(x=x, y=y, label=label, hjust=hjust, vjust=vjust), size=5, fontface="bold") +
      geom_point(data=subset(plot_data, !IsHighlight), aes(x=x, y=y, size=TotalR2), color="grey85", alpha=0.5) +
      geom_point(data=subset(plot_data, IsHighlight), aes(x=x, y=y, fill=Dominant, size=TotalR2), shape=21, color="white", stroke=0.5, alpha=0.9) +
      geom_text_repel(data=subset(plot_data, IsHighlight), aes(x=x, y=y, label=CleanFeature), size=3, max.overlaps=50, bg.color="white", bg.r=0.1, box.padding=0.3) +
      scale_fill_manual(values=c("PA"="#FFD700", "Adiposity"="#4DBBD5", "Dietary intake"="#00A087")) +
      scale_size(range=c(1, 6)) + theme_void() + coord_fixed(clip="off") + labs(title=NULL) + theme(legend.position="right", plot.margin=margin(2,2,2,2,"cm"))
    
    ggsave(file.path(DIR_OUT, paste0(g_name, "_Triangle_Plot.pdf")), p_tri, width=10, height=9)
  }
}

# --- 6. Generate Consolidated Supplementary Table (8 Sheets: 4 Heatmap + 4 Triangle) ---
message("\n=== Generating Consolidated Supplementary Table (8 Sheets: 4 Heatmap + 4 Triangle) ===")

wb <- createWorkbook()
header_style <- createStyle(fontName = "Arial", fontSize = 10, fontColour = "white", fgFill = "#4F81BD", textDecoration = "bold")
sig_style    <- createStyle(fontColour = "#B2182B", textDecoration = "bold")
sig_red      <- createStyle(fontColour = "#B2182B", textDecoration = "bold")
sig_blue     <- createStyle(fontColour = "#2166AC", textDecoration = "bold")
add_stars <- function(p) dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")

# [Write Part 1]: Write the first 4 Heatmap data sheets
for(label in names(HEATMAP_RESULTS_LIST)) {
  # Ensure order is verified again
  df_hm <- HEATMAP_RESULTS_LIST[[label]] %>%
    dplyr::arrange(factor(Feature, levels = HEATMAP_ORDER_LIST[[label]]))
  
  sheet_n <- paste0("HM ", gsub("_", " ", label))
  if(nchar(sheet_n) > 31) sheet_n <- substr(sheet_n, 1, 31)
  
  addWorksheet(wb, sheet_n)
  writeData(wb, sheet_n, df_hm)
  # Apply general styling
  addStyle(wb, sheet_n, header_style, rows = 1, cols = 1:ncol(df_hm), gridExpand = TRUE)
  
  # Positive correlation in red, negative in blue, significant stars in red
  conditionalFormatting(wb, sheet_n, cols = 3, rows = 2:(nrow(df_hm)+1), rule = ">0", style = sig_red)
  conditionalFormatting(wb, sheet_n, cols = 3, rows = 2:(nrow(df_hm)+1), rule = "<0", style = sig_blue)
  conditionalFormatting(wb, sheet_n, cols = 5, rows = 2:(nrow(df_hm)+1), rule = "*", type = "contains", style = sig_style)
}

# [Write Part 2]: Write the last 4 Triangle full-set data sheets
for(label in names(TRIANGLE_RESULTS_LIST)) {
  df <- TRIANGLE_RESULTS_LIST[[label]]
  
  coords <- coord_ternary(df$Activity, df$Macros, df$Clinical)
  df$x <- coords$x; df$y <- coords$y
  df$Dist_Activity <- sqrt((df$x - 0)^2 + (df$y - 0)^2)
  df$Dist_Diet     <- sqrt((df$x - 1)^2 + (df$y - 0)^2)
  df$Dist_Clinical <- sqrt((df$x - 0.5)^2 + (df$y - sqrt(3)/2)^2)
  
  df$Driver_Class <- apply(df[, c("Dist_Activity", "Dist_Diet", "Dist_Clinical")], 1, function(x) {
    c("Physical Activity", "Dietary Composition", "Metabolic Phenotypes")[which.min(x)]
  })
  
  df$Dominant_R2 <- mapply(function(cls, a, d, c) { if(cls == "Physical Activity") return(a) else if(cls == "Dietary Composition") return(d) else return(c) }, df$Driver_Class, df$Activity, df$Macros, df$Clinical)
  df$Dominant_P <- mapply(function(cls, a_p, d_p, c_p) { if(cls == "Physical Activity") return(a_p) else if(cls == "Dietary Composition") return(d_p) else return(c_p) }, df$Driver_Class, df$Activity_P, df$Macros_P, df$Clinical_P)
  
  # No truncation, no filtering (keep full set!), and apply uniform name cleaning
  final_table <- df %>% 
    dplyr::mutate(Act_Sig = add_stars(Activity_P), Diet_Sig = add_stars(Macros_P), Clin_Sig = add_stars(Clinical_P)) %>%
    dplyr::mutate(Feature = clean_rownames_fn(Feature)) %>% 
    dplyr::arrange(factor(Feature, levels = HEATMAP_ORDER_LIST[[label]])) %>%
    dplyr::select(Driver_Class, Feature, Activity_R2 = Activity, Activity_P, Act_Sig, Diet_R2 = Macros, Diet_P = Macros_P, Diet_Sig, Clinical_R2 = Clinical, Clinical_P = Clinical_P, Clin_Sig)
  
  if(nrow(final_table) > 0) {
    sheet_n <- paste0("Tri ", gsub("_", " ", label))
    if(nchar(sheet_n) > 31) sheet_n <- substr(sheet_n, 1, 31)
    
    addWorksheet(wb, sheet_n)
    writeData(wb, sheet_n, final_table)
    addStyle(wb, sheet_n, header_style, rows = 1, cols = 1:ncol(final_table), gridExpand = TRUE)
    
    for(col_idx in which(colnames(final_table) %in% c("Act_Sig", "Diet_Sig", "Clin_Sig"))) {
      conditionalFormatting(wb, sheet_n, cols = col_idx, rows = 2:(nrow(final_table)+1), rule = "*", type = "contains", style = sig_style)
    }
  }
}

out_file <- file.path(DIR_OUT, "Data_S6_Triangle_Drivers.xlsx")
saveWorkbook(wb, out_file, overwrite = TRUE)

# --- 7. Standalone Legend Export ---
message("\n=== Generating Standalone Legend ===")

# Force ensure packages are loaded
if (!require("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
if (!require("circlize", quietly = TRUE)) install.packages("circlize")
library(ComplexHeatmap)
library(circlize)

pdf(file.path(DIR_OUT, "Fig3F_Standalone_Legend.pdf"), width = 3, height = 4)

# 1. Build correlation color bar legend
col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("#3C5488", "white", "#E64B35"))
heatmap_legend <- Legend(col_fun = col_fun, title = "Correlation", direction = "horizontal")

# 2. Build category color block legend
anno_legend <- Legend(labels = c("PA", "Adiposity", "Dietary intake"), 
                      title = "Category", 
                      legend_gp = gpar(fill = c("#FFD700", "#4DBBD5", "#00A087")))

# 3. Pack, draw, and save
draw(packLegend(heatmap_legend, anno_legend, direction = "vertical"))
dev.off()

message("=======================================================")
message("PIPELINE COMPLETED: Figures & Data_S6 saved cleanly in V2.0 Directory!")
message("=======================================================")
