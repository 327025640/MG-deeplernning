# 02b_build_metadata_injection_v3.R
# Step 02b：构建逐样本 metadata 注入查找表（453 行），供 Step 03 merge 时注入。
# 规则（用户 2026-08-28 裁决）：
#   - Merge_Include 全部 TRUE（GSE141383 对象内 FALSE 在 merge 时被本表覆盖）。
#   - Treatment_Status 归一化：No_reported_default -> Not_reported；TAGORE 描述串 -> No；
#     Unknown 保持；对象无列 -> Not_reported。原始值存 Treatment_Status_Original。
#   - Disease_Abbrev_Analysis：14 套取权威表；GSE103224 按 series matrix（GBM grade IV/recurrent -> rGBM，
#     Anaplastic Astrocytoma -> Anaplastic_astrocytoma）；GSE141383 全部 rGBM；权威表外样本 -> Not_reported。
#   - Sex/Age/脑区参考列：14 套取权威表；GSE103224/141383 取 series matrix（缺则 NA）。
# 输出: outputV3/merge_seurat/00_manifest_and_preflight/sample_metadata_injection_v3.csv

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

root <- "D:/data/microglia_new"
out_dir <- file.path(root, "outputV3/merge_seurat/00_manifest_and_preflight")
audit_dir <- file.path(root, "outputV3/merge_seurat/00_dataset_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 权威表（含中文目录与非 UTF-8 字节：ASCII 缓存复制 + fread，同 00a 方案）
acands <- Sys.glob(file.path(root, "output", "01_metadata_overview", "*", "final_atlas_sample_clinical_metadata.csv"))
if (length(acands) != 1L) stop("Authority CSV must resolve to exactly one file")
auth_cache <- file.path(audit_dir, "tmp_authority_cache.csv")
if (!file.copy(acands[1], auth_cache, overwrite = TRUE)) stop("Cannot copy authority CSV")
auth <- fread(auth_cache, encoding = "UTF-8")

samp <- fread(file.path(audit_dir, "dataset16_sample_level_all.csv"))
audit <- fread(file.path(audit_dir, "dataset16_metadata_audit_summary.csv"))
clin_dir <- file.path(root, "output", "01_metadata_overview")
# series matrix 临床表（中文父目录 fread 打不开：复制为 ASCII 缓存后读取，源文件只读）
dspath <- Sys.glob(file.path(root, "output", "01_metadata_overview", "*", "dataset"))
if (length(dspath) != 1L) stop("dataset evidence dir must resolve to exactly one path")
clin_103_cache <- file.path(audit_dir, "tmp_clinical_GSE103224_cache.csv")
clin_141_cache <- file.path(audit_dir, "tmp_clinical_GSE141383_cache.csv")
if (!file.copy(file.path(dspath[1], "GSE103224", "GSE103224_sample_clinical_from_series_matrix.csv"), clin_103_cache, overwrite = TRUE)) stop("Cannot copy GSE103224 clinical CSV")
if (!file.copy(file.path(dspath[1], "GSE141383", "GSE141383_sample_clinical_from_series_matrix.csv"), clin_141_cache, overwrite = TRUE)) stop("Cannot copy GSE141383 clinical CSV")
clin_103 <- fread(clin_103_cache, encoding = "UTF-8")
clin_141 <- fread(clin_141_cache, encoding = "UTF-8")

# --- GSE_ID ---
gse_map <- audit[, .(DatasetID, GSE_ID)]

# --- 逐样本基础表 ---
inj <- samp[, .(DatasetID, SampleID, DonorID)]
inj[, SampleID_Global := paste(DatasetID, SampleID, sep = "__")]
inj[, Merge_Include := TRUE]
inj <- merge(inj, gse_map, by = "DatasetID", all.x = TRUE, sort = FALSE)

# --- Disease_Abbrev_Analysis ---
disease_col <- function(dt, ds) {
  if (ds == "GSE103224") {
    m <- match(dt[DatasetID == ds, SampleID], clin_103$SampleID)
    dx <- clin_103$diagnosis[m]
    out <- ifelse(grepl("Anaplastic Astrocytoma", dx, ignore.case = TRUE), "Anaplastic_astrocytoma",
           ifelse(grepl("Glioblastoma", dx, ignore.case = TRUE), "rGBM", NA_character_))
    return(data.table(DatasetID = ds, SampleID = dt[DatasetID == ds, SampleID],
                      Disease_Abbrev_Analysis = out, Disease_Source = "series_matrix_rule"))
  }
  if (ds == "GSE141383") {
    return(data.table(DatasetID = ds, SampleID = dt[DatasetID == ds, SampleID],
                      Disease_Abbrev_Analysis = "rGBM", Disease_Source = "series_matrix_rule"))
  }
  a <- auth[DatasetID == ds, .(SampleID, Disease_Abbrev_Analysis)]
  x <- merge(dt[DatasetID == ds, .(DatasetID, SampleID)], a, by = "SampleID", all.x = TRUE, sort = FALSE)
  x[is.na(Disease_Abbrev_Analysis), Disease_Abbrev_Analysis := "Not_reported"]
  x[, Disease_Source := ifelse(Disease_Abbrev_Analysis == "Not_reported", "not_in_authority", "authority_table")]
  x
}
dlist <- lapply(unique(inj$DatasetID), function(ds) disease_col(inj, ds))
dall <- rbindlist(dlist, fill = TRUE)
inj <- merge(inj, dall, by = c("DatasetID", "SampleID"), all.x = TRUE, sort = FALSE)

# --- Treatment_Status 归一化 ---
norm_treat <- function(raw, has_col) {
  if (is.na(raw) || !has_col) return("Not_reported")
  if (grepl("No_reported_default", raw, fixed = TRUE)) return("Not_reported")
  if (grepl("^No ", raw) || raw == "No") return("No")
  if (raw == "Unknown") return("Unknown")
  if (raw %in% c("No", "Not_reported", "Unknown", "")) return(raw)
  "Not_reported"
}
tmap <- audit[, .(DatasetID, Treat_Col = Treatment_Status_Column, Treat_Val = Treatment_Status_Values)]
inj <- merge(inj, tmap, by = "DatasetID", all.x = TRUE, sort = FALSE)
inj[, Treatment_Status_Original := fifelse(is.na(Treat_Col) | Treat_Col == "", "no_column_in_object", Treat_Val)]
inj[, Treatment_Status := mapply(norm_treat, Treat_Val, !is.na(Treat_Col) & Treat_Col != "")]
inj[, c("Treat_Col", "Treat_Val") := NULL]

# --- Sex/Age/脑区参考列 ---
inj[, Sex_Authority := NA_character_]
inj[, Age_Authority := NA_character_]
inj[, Brain_Region_Authority := NA_character_]
inj[, Brain_Region_Group_Authority := NA_character_]
for (ds in unique(inj$DatasetID)) {
  idx <- inj$DatasetID == ds
  if (ds == "GSE103224") {
    m <- match(inj$SampleID[idx], clin_103$SampleID)
    inj$Sex_Authority[idx] <- clin_103$gender[m]
    inj$Age_Authority[idx] <- clin_103$age[m]
    inj$Brain_Region_Authority[idx] <- clin_103$location[m]
  } else if (ds == "GSE141383") {
    m <- match(inj$SampleID[idx], clin_141$SampleID)
    inj$Sex_Authority[idx] <- clin_141$sex[m]
  } else {
    a <- auth[DatasetID == ds]
    m <- match(inj$SampleID[idx], a$SampleID)
    inj$Sex_Authority[idx] <- a$Sex[m]
    inj$Age_Authority[idx] <- a$Age[m]
    inj$Brain_Region_Authority[idx] <- a$Brain_Region[m]
    inj$Brain_Region_Group_Authority[idx] <- a$Brain_Region_Group[m]
  }
}

# --- 权威表外 20 样本的疾病覆盖（用户 2026-08-29 裁决：按对象内疾病标签归入对应疾病）---
# 标签映射沿用图谱词表：AD_LOAD->AD, Control->Ctrl, Lung_cancer_brain_metastasis->LCBM,
# Medulloblastoma->MB, 其余（ALS/ALS_FTD/MS/FTD_MND/DLB）直接使用；FTD_MND 与 DLB 为独立标签。
override <- data.table(
  DatasetID = c(rep("TUDDENHAM_2024", 10), rep("LUNG_BRAIN_MET_2022", 5),
                rep("WANG_C9_DLPFC_PNAS_2025", 2), "MEDULLOBLASTOMA_2021", "MEDULLOBLASTOMA_2021", "KISS_2023"),
  SampleID = c("GSM6186021", "GSM6186050", "GSM6186051", "GSM6186052", "GSM6186053",
               "GSM6186066", "GSM6186067", "GSM6186077", "GSM6186080", "GSM6186087",
               "LC02", "LC03", "LC07", "LC08", "LC10",
               "GSM8735066", "GSM8735075", "943", "966-recurrence", "GSM7108093"),
  Dis_Override = c("AD", "ALS", "ALS", "ALS", "ALS",
                   "AD", "AD", "Ctrl", "ALS_FTD", "AD",
                   "LCBM", "LCBM", "LCBM", "LCBM", "LCBM",
                   "FTD_MND", "DLB", "MB", "MB", "MS")
)
inj[override, Disease_Abbrev_Analysis := i.Dis_Override, on = .(DatasetID, SampleID)]
inj[override, on = .(DatasetID, SampleID), Disease_Source := "object_label_user_approved_2026-08-29"]
if (any(is.na(inj$Disease_Abbrev_Analysis))) stop("Disease override left NA values")

setcolorder(inj, c("DatasetID", "SampleID", "SampleID_Global", "DonorID", "GSE_ID",
                   "Merge_Include", "Disease_Abbrev_Analysis", "Disease_Source",
                   "Treatment_Status", "Treatment_Status_Original",
                   "Sex_Authority", "Age_Authority", "Brain_Region_Authority", "Brain_Region_Group_Authority"))
setorder(inj, DatasetID, SampleID)
if (anyDuplicated(inj$SampleID_Global)) stop("Duplicate SampleID_Global in injection table")
fwrite(inj, file.path(out_dir, "sample_metadata_injection_v3.csv"))

cat("rows:", nrow(inj), " datasets:", uniqueN(inj$DatasetID), "\n")
cat("disease:", paste(names(sort(table(inj$Disease_Abbrev_Analysis), decreasing = TRUE)), collapse = ", "), "\n")
cat("treatment:", paste(names(sort(table(inj$Treatment_Status), decreasing = TRUE)), collapse = ", "), "\n")
