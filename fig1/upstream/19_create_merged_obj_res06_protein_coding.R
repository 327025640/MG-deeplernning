# =============================================================================
# 19_create_merged_obj_res06_protein_coding.R
# 目的：
#   1. 复用正式对象中已保存的 RNA_snn 图，冻结 resolution = 0.6；
#   2. 将 active.ident 更新为全局 resolution 0.6 聚类；
#   3. 仅保留 GENCODE v50 protein-coding genes；
#   4. 建立主图疾病大类 Disease，并保留细分疾病 disease_analyse；
#   5. 将 meta.data 精简为下游分析所需的核心字段；
#   6. 保存为 E:/micro_output/v3_dims16/objects/merged_obj.rds。
#
# 说明：本脚本不会重新运行 FindNeighbors、PCA、Harmony 或 UMAP。
# =============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

# -----------------------------------------------------------------------------
# 1. 冻结路径与参数
# -----------------------------------------------------------------------------
INPUT_PATH <- "E:/micro_output/v3_dims16/objects/microglia_atlas_harmony_dims16_v3.rds"
PROTEIN_CODING_REF <- "D:/data/microglia_new/authority/03_reference_configs/gencode_v50_protein_coding_genes.csv"
AUTHORITY_INJECTION <- "D:/data/microglia_new/outputV3/merge_seurat/00_manifest_and_preflight/sample_metadata_injection_v3.csv"
RES06_SOURCE <- "E:/micro_output/v3_dims16/resolution_determine/source_data/umap_source_res0.6.csv"
OUTPUT_DIR <- "E:/micro_output/v3_dims16/objects"
OUTPUT_PATH <- file.path(OUTPUT_DIR, "merged_obj.rds")
AUDIT_PATH <- file.path(OUTPUT_DIR, "merged_obj_audit.csv")
DISEASE_MAPPING_PATH <- file.path(OUTPUT_DIR, "merged_obj_disease_mapping.csv")
SAMPLE_DISEASE_AUDIT_PATH <- file.path(OUTPUT_DIR, "merged_obj_sample_disease_authority_audit.csv")
SESSION_PATH <- file.path(OUTPUT_DIR, "merged_obj_sessionInfo.txt")
INTERPRETATION_PATH <- "E:/micro_output/v3_dims16/解读.md"

SNN_GRAPH <- "RNA_snn"
RESOLUTION <- 0.6
RESOLUTION_COL <- "RNA_snn_res.0.6"
RANDOM_SEED <- 42L
EXPECTED_CELLS <- 424222L
EXPECTED_PROTEIN_CODING_GENES <- 18687L

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

required_paths <- c(INPUT_PATH, PROTEIN_CODING_REF, AUTHORITY_INJECTION, RES06_SOURCE)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("缺少必需输入文件：", paste(missing_paths, collapse = "; "))
}
if (file.exists(OUTPUT_PATH)) {
  stop("目标文件已存在，为避免覆盖而停止：", OUTPUT_PATH)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("缺少 digest 包，无法生成正式对象 SHA256 审计。")
}

log_msg <- function(...) {
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
             paste0(..., collapse = ""), "\n"))
}

sort_cluster_levels <- function(x) {
  levels_chr <- unique(as.character(x))
  levels_num <- suppressWarnings(as.numeric(levels_chr))
  if (all(is.finite(levels_num))) {
    levels_chr[order(levels_num)]
  } else {
    sort(levels_chr)
  }
}

# -----------------------------------------------------------------------------
# 2. 读取并核验正式 Seurat 对象
# -----------------------------------------------------------------------------
log_msg("读取正式对象：", INPUT_PATH)
obj <- readRDS(INPUT_PATH)

if (!inherits(obj, "Seurat")) stop("输入不是 Seurat 对象。")
if (ncol(obj) != EXPECTED_CELLS) {
  stop("细胞数不符合正式对象契约：expected=", EXPECTED_CELLS,
       ", observed=", ncol(obj))
}
if (!"RNA" %in% Assays(obj)) stop("输入对象缺少 RNA assay。")
if (!all(c("counts", "data") %in% Layers(obj[["RNA"]]))) {
  stop("RNA assay 必须同时包含 counts 和 data layer。")
}
if (!SNN_GRAPH %in% Graphs(obj)) stop("输入对象缺少已保存的 RNA_snn 图。")
if (anyDuplicated(colnames(obj))) stop("Seurat 对象存在重复 CellID。")

required_meta <- c(
  "CellID", "Original_Barcode", "DatasetID", "GSM", "SampleID",
  "SampleID_Global", "DonorID", "Disease_Abbrev_Analysis",
  "Disease_original", "Disease_Class", "Disease_stage", "Grade",
  "Age", "Age_Status", "Sex", "Brain_Region", "sc_or_sn",
  "Tissue_Status", "Treatment_Status", "Experimental_Perturbation_Status",
  "Merge_Include", "User_Approved_Include", "nCount_RNA", "nFeature_RNA",
  "percent.mt", "percent.ribo", "scDblFinder.class", "scDblFinder.score",
  "BroadLineage", "MicrogliaGate"
)
missing_meta <- setdiff(required_meta, colnames(obj[[]]))
if (length(missing_meta) > 0L) {
  stop("输入对象缺少核心 metadata 列：", paste(missing_meta, collapse = ", "))
}
if (!identical(as.character(obj$CellID), colnames(obj))) {
  stop("metadata CellID 与 Seurat 细胞列名不完全一致。")
}

# -----------------------------------------------------------------------------
# 3. 复用 RNA_snn 图计算并冻结全局 resolution 0.6
# -----------------------------------------------------------------------------
log_msg("复用 ", SNN_GRAPH, " 运行 FindClusters(resolution=0.6)；不重算 neighbors。")
obj <- FindClusters(
  object = obj,
  graph.name = SNN_GRAPH,
  algorithm = 1,
  resolution = RESOLUTION,
  random.seed = RANDOM_SEED,
  verbose = FALSE
)

if (!RESOLUTION_COL %in% colnames(obj[[]])) {
  stop("FindClusters 后未生成预期列：", RESOLUTION_COL)
}

cluster_res06 <- as.character(obj[[RESOLUTION_COL, drop = TRUE]])
if (length(cluster_res06) != ncol(obj) || anyNA(cluster_res06)) {
  stop("resolution 0.6 聚类标签长度错误或包含 NA。")
}
cluster_levels <- sort_cluster_levels(cluster_res06)
obj[[RESOLUTION_COL]] <- factor(cluster_res06, levels = cluster_levels)
obj$seurat_clusters <- obj[[RESOLUTION_COL, drop = TRUE]]
Idents(obj) <- "seurat_clusters"

# 与用户用于选择分辨率的正式 res0.6 UMAP source data 逐细胞核对。
res06_source <- fread(
  RES06_SOURCE,
  colClasses = list(character = c("CellID", "Cluster")),
  select = c("CellID", "Cluster")
)
if (nrow(res06_source) != ncol(obj) || anyDuplicated(res06_source$CellID)) {
  stop("res0.6 source data 的细胞数或 CellID 唯一性异常。")
}
source_index <- match(colnames(obj), res06_source$CellID)
if (anyNA(source_index)) stop("res0.6 source data 未覆盖全部 Seurat 细胞。")
expected_cluster <- res06_source$Cluster[source_index]
if (!identical(as.character(Idents(obj)), expected_cluster)) {
  mismatch_n <- sum(as.character(Idents(obj)) != expected_cluster)
  stop("重算的 resolution 0.6 与已审阅图形标签不一致；不一致细胞数：", mismatch_n)
}
rm(res06_source, source_index, expected_cluster)
invisible(gc())

# -----------------------------------------------------------------------------
# 4. 建立主图疾病大类与细分疾病列
# -----------------------------------------------------------------------------
# disease_analyse 保留当前权威疾病/亚型标签；Disease 用于后续主图。
obj$disease_analyse <- as.character(obj$Disease_Abbrev_Analysis)
obj$Disease <- obj$disease_analyse
obj$Disease[obj$disease_analyse %in% c("rGBM", "pGBM")] <- "GBM"
obj$Disease[obj$disease_analyse == "Anaplastic_astrocytoma"] <- "LGG"

if (anyNA(obj$disease_analyse) || any(trimws(obj$disease_analyse) == "")) {
  stop("disease_analyse 包含 NA 或空值。")
}
if (anyNA(obj$Disease) || any(trimws(obj$Disease) == "")) {
  stop("Disease 包含 NA 或空值。")
}
if (any(obj$Disease %in% c("rGBM", "pGBM", "Anaplastic_astrocytoma"))) {
  stop("Disease 主图列仍包含应被合并的疾病亚型。")
}
if (!all(obj$Disease[obj$disease_analyse %in% c("rGBM", "pGBM")] == "GBM")) {
  stop("rGBM/pGBM -> GBM 映射失败。")
}
if (!all(obj$Disease[obj$disease_analyse == "Anaplastic_astrocytoma"] == "LGG")) {
  stop("Anaplastic_astrocytoma -> LGG 映射失败。")
}

disease_mapping_audit <- data.table(
  disease_analyse = obj$disease_analyse,
  Disease = obj$Disease
)[, .(n_cells = .N), by = .(disease_analyse, Disease)]
setorder(disease_mapping_audit, Disease, disease_analyse)

# 按 SampleID_Global 与冻结的453样本权威注入表逐样本核对疾病。
authority <- fread(AUTHORITY_INJECTION)
required_authority_cols <- c(
  "SampleID_Global", "DatasetID", "SampleID", "Disease_Abbrev_Analysis"
)
missing_authority_cols <- setdiff(required_authority_cols, names(authority))
if (length(missing_authority_cols) > 0L) {
  stop("权威注入表缺少列：", paste(missing_authority_cols, collapse = ", "))
}
authority <- authority[, ..required_authority_cols]
if (anyDuplicated(authority$SampleID_Global)) {
  stop("权威注入表的 SampleID_Global 不唯一。")
}

object_sample_disease <- unique(data.table(
  SampleID_Global = as.character(obj$SampleID_Global),
  DatasetID = as.character(obj$DatasetID),
  SampleID = as.character(obj$SampleID),
  disease_analyse = as.character(obj$disease_analyse)
))
if (anyDuplicated(object_sample_disease$SampleID_Global)) {
  stop("对象内至少一个 SampleID_Global 对应多个 DatasetID、SampleID 或疾病标签。")
}

authority_index <- match(object_sample_disease$SampleID_Global, authority$SampleID_Global)
sample_disease_audit <- copy(object_sample_disease)
sample_disease_audit[, `:=`(
  Authority_DatasetID = authority$DatasetID[authority_index],
  Authority_SampleID = authority$SampleID[authority_index],
  Authority_Disease = authority$Disease_Abbrev_Analysis[authority_index]
)]
sample_disease_audit[, `:=`(
  SampleID_Global_found = !is.na(authority_index),
  DatasetID_match = !is.na(authority_index) & DatasetID == Authority_DatasetID,
  SampleID_match = !is.na(authority_index) & SampleID == Authority_SampleID,
  Disease_match = !is.na(authority_index) & disease_analyse == Authority_Disease
)]
sample_disease_audit[, Status := fifelse(
  SampleID_Global_found & DatasetID_match & SampleID_match & Disease_match,
  "PASS", "FAIL"
)]
setorder(sample_disease_audit, DatasetID, SampleID_Global)

if (nrow(sample_disease_audit) != 453L) {
  stop("对象中的 SampleID_Global 数量不是冻结的453：observed=", nrow(sample_disease_audit))
}
if (any(sample_disease_audit$Status != "PASS")) {
  fwrite(sample_disease_audit, SAMPLE_DISEASE_AUDIT_PATH)
  stop("逐样本疾病权威核验失败；详见：", SAMPLE_DISEASE_AUDIT_PATH)
}
rm(authority, authority_index, object_sample_disease)
invisible(gc())

# -----------------------------------------------------------------------------
# 5. 仅保留 GENCODE v50 protein-coding genes
# -----------------------------------------------------------------------------
log_msg("读取 GENCODE v50 protein-coding gene reference。")
pc_ref <- fread(PROTEIN_CODING_REF)
required_pc_cols <- c("Gene_Symbol", "Gene_Type")
missing_pc_cols <- setdiff(required_pc_cols, names(pc_ref))
if (length(missing_pc_cols) > 0L) {
  stop("protein-coding reference 缺少列：", paste(missing_pc_cols, collapse = ", "))
}

pc_symbols <- unique(pc_ref[Gene_Type == "protein_coding" & !is.na(Gene_Symbol), Gene_Symbol])
pc_use <- intersect(rownames(obj), pc_symbols)
if (length(pc_use) != EXPECTED_PROTEIN_CODING_GENES) {
  stop("对象中的 protein-coding gene 数与冻结审计不一致：expected=",
       EXPECTED_PROTEIN_CODING_GENES, ", observed=", length(pc_use))
}

original_nCount_RNA <- obj$nCount_RNA
original_nFeature_RNA <- obj$nFeature_RNA
old_variable_features <- VariableFeatures(obj)
obj <- subset(obj, features = pc_use)
VariableFeatures(obj) <- intersect(old_variable_features, rownames(obj))
# 用户决定：QC列必须原样继承全基因输入对象，不随protein-coding子集重算。
obj$nCount_RNA <- original_nCount_RNA
obj$nFeature_RNA <- original_nFeature_RNA

if (nrow(obj) != EXPECTED_PROTEIN_CODING_GENES ||
    !all(rownames(obj) %in% pc_symbols)) {
  stop("protein-coding gene 子集核验失败。")
}
rm(pc_ref, pc_symbols, pc_use, old_variable_features,
   original_nCount_RNA, original_nFeature_RNA)
invisible(gc())

# -----------------------------------------------------------------------------
# 6. 精简 metadata 与 assay layers
# -----------------------------------------------------------------------------
keep_meta <- c(
  # 细胞、数据集、样本和供体追溯
  "CellID", "Original_Barcode", "DatasetID", "GSM", "SampleID",
  "SampleID_Global", "DonorID",
  # 正式疾病与临床/组织信息
  "Disease", "disease_analyse", "Disease_Abbrev_Analysis",
  "Disease_original", "Disease_Class",
  "Disease_stage", "Grade", "Age", "Age_Status", "Sex",
  "Brain_Region", "sc_or_sn", "Tissue_Status",
  # 纳入、治疗和实验扰动审计
  "Treatment_Status", "Experimental_Perturbation_Status",
  "Merge_Include", "User_Approved_Include",
  # 细胞级 QC 与 doublet 信息
  "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo",
  "scDblFinder.class", "scDblFinder.score",
  # 原 microglia 纳入证据
  "BroadLineage", "MicrogliaGate",
  # 最终全局 resolution 0.6
  RESOLUTION_COL, "seurat_clusters"
)

obj@meta.data <- obj@meta.data[, keep_meta, drop = FALSE]

# counts/data供表达与下游统计使用；删除可重建的scale.data以减小对象体积。
obj <- DietSeurat(
  object = obj,
  assays = "RNA",
  layers = c("counts", "data"),
  dimreducs = Reductions(obj),
  graphs = Graphs(obj),
  misc = TRUE
)
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "seurat_clusters"

# -----------------------------------------------------------------------------
# 7. 保存前核验
# -----------------------------------------------------------------------------
if (!identical(colnames(obj[[]]), keep_meta)) stop("最终 metadata 列不符合冻结清单。")
if (!identical(rownames(obj[[]]), colnames(obj))) stop("最终 metadata 行名与 CellID 不一致。")
if (!identical(as.character(Idents(obj)), as.character(obj$seurat_clusters))) {
  stop("最终 active.ident 与 seurat_clusters 不一致。")
}
if (anyNA(obj$seurat_clusters)) stop("最终 seurat_clusters 包含 NA。")
if (!setequal(Layers(obj[["RNA"]]), c("counts", "data"))) {
  stop("最终 RNA assay 应只包含 counts 和 data layer。")
}
if (!isTRUE(validObject(obj))) stop("最终 Seurat 对象 validObject() 未通过。")

audit <- data.table(
  item = c(
    "input_path", "output_path", "n_cells", "n_genes_protein_coding",
    "n_clusters_res0.6", "active_ident", "metadata_columns_n",
    "metadata_columns", "rna_layers", "assays", "reductions", "graphs",
    "neighbors_recomputed", "resolution", "random_seed",
    "main_disease_levels", "disease_analyse_levels",
    "disease_mapping_rule_1", "disease_mapping_rule_2",
    "sample_disease_authority_rows", "sample_disease_authority_status",
    "n_datasets", "n_samples", "qc_columns_definition"
  ),
  value = c(
    INPUT_PATH, OUTPUT_PATH, ncol(obj), nrow(obj), length(levels(Idents(obj))),
    "seurat_clusters (= RNA_snn_res.0.6)", ncol(obj[[]]),
    paste(colnames(obj[[]]), collapse = ";"),
    paste(Layers(obj[["RNA"]]), collapse = ";"),
    paste(Assays(obj), collapse = ";"), paste(Reductions(obj), collapse = ";"),
    paste(Graphs(obj), collapse = ";"), FALSE, RESOLUTION, RANDOM_SEED,
    length(unique(obj$Disease)), length(unique(obj$disease_analyse)),
    "rGBM+pGBM -> GBM", "Anaplastic_astrocytoma -> LGG",
    nrow(sample_disease_audit),
    if (all(sample_disease_audit$Status == "PASS")) "PASS" else "FAIL",
    uniqueN(obj$DatasetID), uniqueN(obj$SampleID_Global),
    "nCount_RNA/nFeature_RNA retained unchanged from input full-gene object metadata"
  )
)

# -----------------------------------------------------------------------------
# 8. 原子保存、重读核验与 SHA256
# -----------------------------------------------------------------------------
tmp_path <- tempfile(pattern = "merged_obj_", tmpdir = OUTPUT_DIR, fileext = ".tmp.rds")
on.exit({
  if (file.exists(tmp_path)) unlink(tmp_path)
}, add = TRUE)

log_msg("保存临时对象：", tmp_path)
saveRDS(obj, tmp_path, compress = "gzip")
expected_cells <- colnames(obj)
expected_meta <- colnames(obj[[]])
expected_nCount_RNA <- obj$nCount_RNA
expected_nFeature_RNA <- obj$nFeature_RNA
rm(obj)
invisible(gc())

log_msg("重读临时对象并执行最终核验。")
check_obj <- readRDS(tmp_path)
if (!inherits(check_obj, "Seurat") || !isTRUE(validObject(check_obj))) {
  stop("保存后的对象类型或 validObject() 核验失败。")
}
if (nrow(check_obj) != EXPECTED_PROTEIN_CODING_GENES ||
    ncol(check_obj) != EXPECTED_CELLS ||
    !identical(colnames(check_obj), expected_cells) ||
    !identical(colnames(check_obj[[]]), expected_meta) ||
    !identical(check_obj$nCount_RNA, expected_nCount_RNA) ||
    !identical(check_obj$nFeature_RNA, expected_nFeature_RNA) ||
    anyNA(check_obj$seurat_clusters) ||
    anyNA(check_obj$Disease) || anyNA(check_obj$disease_analyse) ||
    any(check_obj$Disease %in% c("rGBM", "pGBM", "Anaplastic_astrocytoma")) ||
    !identical(as.character(Idents(check_obj)), as.character(check_obj$seurat_clusters))) {
  stop("保存后重读核验失败。")
}

rm(check_obj, expected_cells, expected_meta,
   expected_nCount_RNA, expected_nFeature_RNA)
invisible(gc())

sha256 <- digest::digest(file = tmp_path, algo = "sha256", serialize = FALSE)
audit <- rbind(
  audit,
  data.table(
    item = c("output_size_bytes", "sha256"),
    value = c(file.info(tmp_path)$size, sha256)
  )
)
# 先完成审计sidecar，再将已核验RDS原子重命名为正式对象。
fwrite(audit, AUDIT_PATH)
fwrite(disease_mapping_audit, DISEASE_MAPPING_PATH)
fwrite(sample_disease_audit, SAMPLE_DISEASE_AUDIT_PATH)
capture.output(sessionInfo(), file = SESSION_PATH)
if (!file.rename(tmp_path, OUTPUT_PATH)) stop("无法将临时文件重命名为正式输出。")

interpretation_lines <- c(
  "",
  "## resolution 0.6 protein-coding 精简对象",
  "",
  "| 项目 | 内容 |",
  "|---|---|",
  "| 分析目的 | 冻结用户选定的全局 resolution 0.6，并生成仅含蛋白编码基因和核心 metadata 的下游对象。 |",
  paste0("| 输入数据 | `", INPUT_PATH, "`；", EXPECTED_CELLS, " 个细胞。 |"),
  paste0("| 方法 | 复用已保存的 `RNA_snn` 图运行 `FindClusters(resolution=0.6, seed=42)`；不重算 neighbors；按 GENCODE v50 筛选 protein-coding genes。 |"),
  paste0("| 结果 | 输出 ", EXPECTED_PROTEIN_CODING_GENES, " 个 protein-coding genes、", EXPECTED_CELLS, " 个细胞、30个resolution 0.6簇、453个样本、16个数据集和 ", audit[item == "metadata_columns_n", value], " 个核心 metadata 字段；正式对象为 `", OUTPUT_PATH, "`。 |"),
  "| 疾病字段 | 主图使用 `Disease`：rGBM+pGBM 合并为 GBM，Anaplastic astrocytoma 并入 LGG；亚型分析使用保留原权威标签的 `disease_analyse`。 |",
  "| QC字段 | `nCount_RNA/nFeature_RNA`原样继承输入对象metadata，代表原始全基因对象QC，不按protein-coding子集重算。 |",
  "| 生物学意义 | 后续 cluster marker、NMF、pseudobulk 和疾病比较统一使用 resolution 0.6 身份及未整合 RNA counts/data；主图与亚型分析使用不同疾病粒度。 |",
  "| 限制 | 聚类来源于既有 Harmony dims16 SNN 图；protein-coding 子集不会重新计算 PCA、Harmony 或 UMAP。 |"
)
cat(interpretation_lines, sep = "\n", file = INTERPRETATION_PATH, append = file.exists(INTERPRETATION_PATH))

log_msg("完成：", OUTPUT_PATH)
log_msg("细胞数=", EXPECTED_CELLS,
        "；protein-coding genes=", EXPECTED_PROTEIN_CODING_GENES,
        "；clusters=", audit[item == "n_clusters_res0.6", value],
        "；SHA256=", sha256)
