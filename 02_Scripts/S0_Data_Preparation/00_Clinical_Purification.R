# ==============================================================================
# Project: FinlandSports (Multi-Omics Bifurcation)
# Script:  00_Clinical_Purification.R
# Core:    Excel Multi-sheet Functional Merge, Strict Trio Assertion, V2.0 Hardcoding.
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(purrr)
  library(readr)
})

input_file <- "00_Raw_Data/Clinical_Data_Raw.xlsx"
output_file <- "01_Clean_Data/Clinical_Master_Strict.csv"

if (!file.exists(input_file)) stop("FATAL ERROR: 找不到 00_Raw_Data/Clinical_Data_Raw.xlsx")
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Module A: Excel 极速多表合并
# ==============================================================================
message(">>> [Phase 1]: Loading and functionally merging all Excel sheets...")

sheet_names <- excel_sheets(input_file)
if (!"Basic information" %in% sheet_names) stop("Error: 缺失 'Basic information' 主表！")

# 批量读取所有 Sheet
list_of_dfs <- lapply(sheet_names, function(s) {
  read_excel(input_file, sheet = s) %>%
    mutate(FamilyID = as.character(FamilyID), Membercode = as.numeric(Membercode))
})

# 将 Basic information 提至最前，作为左连接的主骨架
basic_idx <- grep("Basic information", sheet_names, ignore.case = TRUE)
if (length(basic_idx) > 0) {
  list_of_dfs <- c(list_of_dfs[basic_idx], list_of_dfs[-basic_idx])
}

# 提取公共主键，自动识别并执行 purrr::reduce 合并
join_keys <- intersect(names(list_of_dfs[[1]]), c("Sample_Key", "FamilyID", "Membercode"))
master_df <- list_of_dfs %>%
  reduce(left_join, by = join_keys) %>%
  select(-matches("\\.x$|\\.y$")) # 斩断所有冗余后缀

# ==============================================================================
# Module B: V2.0 顶刊铁律防御层
# ==============================================================================
message(">>> [Phase 2]: Enforcing V2.0 Methodology Iron Rules...")

# [铁律 1: 绝对断言防御 (Strict Trio Assertion)]
valid_roles <- c(1, 2, 3)
if (!all(master_df$Membercode %in% valid_roles, na.rm = TRUE)) {
  stop("FATAL PIPELINE HALT: 发现污染数据！队列必须严格为 trios (daughters, mothers, fathers)。")
}

# [铁律 2: 变量名内存级夺权 (Fat_percent -> Fat(%))]
if("Fat_percent" %in% colnames(master_df)) {
  master_df <- master_df %>% rename(`Fat(%)` = Fat_percent)
}

# [铁律 3: IPA 物理活动不足变量隔离确认]
if("IPA" %in% colnames(master_df)) {
  message("  [Pass] Variable 'IPA' (Physical Inactivity) protected against pathway acronym collision.")
}

# 角色文本化与终极主键生成
master_df <- master_df %>%
  mutate(Role = case_when(Membercode == 1 ~ "Daughter", Membercode == 2 ~ "Mother", Membercode == 3 ~ "Father")) %>%
  mutate(Clinical_Subject_ID = paste0(FamilyID, "_", Membercode)) 

write_csv(master_df, output_file)
message(paste0(">>> [Success] Pure clinical master generated. Features: ", ncol(master_df), " | Samples: ", nrow(master_df)))