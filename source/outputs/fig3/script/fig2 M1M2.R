###############################################################################
# ============================================================================#
# Fig2F / Fig2G
# Microglia-adapted M1-like / M2-like AddModuleScore UMAP
#
# 目标：
#   Fig2F = M1-like score UMAP
#   Fig2G = M2-like score UMAP
#
# 主参考：
# Single-cell transcriptome analysis reveals the immune heterogeneity and
# the repopulation of microglia by Hif1α in mice after spinal cord injury
# Cell Death & Disease, 2022
# DOI: 10.1038/s41419-022-04864-z
#
# 文献明确使用：
#
# M1 markers:
# TNF, Il1b, Icam1, Fas, Ccl3, Cd24a
#
# M2 markers:
# Itga6, Tgfb1, Arg1, Mrc1, Igf1, Ccl5
#
# 人类数据 ortholog symbol：
#
# M1-like:
# TNF, IL1B, ICAM1, FAS, CCL3, CD24
#
# M2-like:
# ITGA6, TGFB1, ARG1, MRC1, IGF1, CCL5
#
# 注意：
# M1/M2 terminology here represents transcriptional programs,
# NOT binary classification of microglial cells.
#
#
# 分析原则：
#
# 1. 先建立文献 gene -> human gene mapping
# 2. 再和正式 RNA assay rownames 取交集
# 3. 所有 missing genes 保存到 audit，但不静默替换
# 4. 只使用实际存在于 RNA assay 的 genes 打分
# 5. 若某 gene set 保留不足，停止
# 6. 全部 424,222 个细胞参与 AddModuleScore
# 7. 使用现成的 umap_harmony16
# 8. 不重新运行 Harmony / UMAP / clustering
# 9. 不修改 merged_obj.rds
#
#
# 输出：
#
# E:/micro_output/fig3/
#
# ├── figures/
# ├── source_data/
# ├── tables/
# └── data/
#
# ============================================================================#
###############################################################################


###############################################################################
# 0. 路径
###############################################################################

INPUT_RDS <- "E:/micro_output/v3_dims16/objects/merged_obj.rds"

OUTPUT_ROOT <- "E:/micro_output/fig3"


FIG_DIR <- file.path(
  OUTPUT_ROOT,
  "figures"
)


SOURCE_DIR <- file.path(
  OUTPUT_ROOT,
  "source_data"
)


TABLE_DIR <- file.path(
  OUTPUT_ROOT,
  "tables"
)


DATA_DIR <- file.path(
  OUTPUT_ROOT,
  "data"
)


###############################################################################
# 1. 固定参数
###############################################################################

ASSAY_USE <- "RNA"

REDUCTION_USE <- "umap_harmony16"

CLUSTER_COL <- "RNA_snn_res.0.6"

SEED <- 42


set.seed(
  SEED
)


###############################################################################
# AddModuleScore 参数
#
# 与 Seurat 默认值保持一致，但显式写出方便复现
###############################################################################

MODULE_NBIN <- 24L

MODULE_CTRL <- 100L


###############################################################################
# Gene set 最低保留标准
#
# 当前每套 signature 原始为6个基因。
#
# 允许少量基因因为对象基因集合原因缺失，
# 但至少需要：
#
# >=4 genes
# 且 retention >=50%
###############################################################################

MIN_GENES_REQUIRED <- 4L

MIN_RETENTION_REQUIRED <- 0.50


###############################################################################
# UMAP绘图 percentile clipping
#
# 原始 AddModuleScore 完全保存。
#
# 这里只对绘图色阶进行 1%-99% clipping，
# 防止少量极端值压缩整体视觉动态范围。
###############################################################################

PLOT_LOW_QUANTILE <- 0.01

PLOT_HIGH_QUANTILE <- 0.99


###############################################################################
# Fig2F M1-like 颜色
#
# 高 M1-like = 蓝色
###############################################################################

M1_LOW_COLOR <- "#F2F2F2"

M1_HIGH_COLOR <- "#2166AC"


###############################################################################
# Fig2G M2-like 颜色
#
# 高 M2-like = 经典 SCI 红色
###############################################################################

M2_LOW_COLOR <- "#F2F2F2"

M2_HIGH_COLOR <- "#B2182B"


###############################################################################
# 字体
###############################################################################

FONT <- "Arial"


###############################################################################
###############################################################################
#                   PART A
#             MICROGLIA SIGNATURE DEFINITION
###############################################################################
###############################################################################


###############################################################################
# 2. 原始文献中的 mouse symbols
#
# 保存原始写法，确保来源可追踪
###############################################################################

M1_REFERENCE_MOUSE <- c(
  "Tnf",
  "Il1b",
  "Icam1",
  "Fas",
  "Ccl3",
  "Cd24a"
)


M2_REFERENCE_MOUSE <- c(
  "Itga6",
  "Tgfb1",
  "Arg1",
  "Mrc1",
  "Igf1",
  "Ccl5"
)


###############################################################################
# 3. Human ortholog symbols
#
# 对于你的正式 human RNA 数据：
#
# Cd24a -> CD24
#
# 其他都是直接的人类 ortholog symbol。
###############################################################################

M1_HUMAN_CANDIDATES <- c(
  "TNF",
  "IL1B",
  "ICAM1",
  "FAS",
  "CCL3",
  "CD24"
)


M2_HUMAN_CANDIDATES <- c(
  "ITGA6",
  "TGFB1",
  "ARG1",
  "MRC1",
  "IGF1",
  "CCL5"
)


###############################################################################
# 4. 建立来源 -> human ortholog 映射表
###############################################################################

M1_MAPPING_REFERENCE <- data.table::data.table(

  signature = "M1-like",

  reference_gene = M1_REFERENCE_MOUSE,

  human_gene = M1_HUMAN_CANDIDATES,

  mapping_type = c(
    "mouse-to-human ortholog",
    "mouse-to-human ortholog",
    "mouse-to-human ortholog",
    "mouse-to-human ortholog",
    "mouse-to-human ortholog",
    "mouse-to-human ortholog: Cd24a -> CD24"
  )
)


M2_MAPPING_REFERENCE <- data.table::data.table(

  signature = "M2-like",

  reference_gene = M2_REFERENCE_MOUSE,

  human_gene = M2_HUMAN_CANDIDATES,

  mapping_type = rep(
    "mouse-to-human ortholog",
    length(M2_REFERENCE_MOUSE)
  )
)


SIGNATURE_REFERENCE <- data.table::rbindlist(
  list(
    M1_MAPPING_REFERENCE,
    M2_MAPPING_REFERENCE
  )
)


SIGNATURE_REFERENCE[
  ,
  reference :=
    "Cell Death & Disease 2022"
]


SIGNATURE_REFERENCE[
  ,
  DOI :=
    "10.1038/s41419-022-04864-z"
]


###############################################################################
# 5. 必要 R 包
###############################################################################

required_packages <- c(
  "Seurat",
  "data.table",
  "ggplot2"
)


missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]


if (length(
  missing_packages
) > 0L) {

  stop(
    paste0(
      "缺少必要 R 包：",
      paste(
        missing_packages,
        collapse = ", "
      ),
      "\n请手动安装后重新运行。",
      "\n本脚本不会自动安装 R 包。"
    )
  )
}


suppressPackageStartupMessages({

  library(Seurat)

  library(data.table)

  library(ggplot2)

})


###############################################################################
# 6. 输入检查
###############################################################################

if (!file.exists(
  INPUT_RDS
)) {

  stop(
    paste0(
      "输入 Seurat 对象不存在：",
      INPUT_RDS
    )
  )
}


###############################################################################
# 7. 创建输出目录
###############################################################################

for (dir_i in c(
  FIG_DIR,
  SOURCE_DIR,
  TABLE_DIR,
  DATA_DIR
)) {

  dir.create(
    dir_i,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


###############################################################################
###############################################################################
#                       PART B
#                  SEURAT OBJECT QC
###############################################################################
###############################################################################


###############################################################################
# 8. 读取正式对象
###############################################################################

message("")
message("============================================================")
message("读取正式 Seurat 对象")
message("============================================================")

message(
  INPUT_RDS
)


obj <- readRDS(
  INPUT_RDS
)


###############################################################################
# 9. 正式对象尺寸
###############################################################################

if (ncol(
  obj
) != 424222L) {

  stop(
    paste0(
      "对象细胞数异常。",
      "\n预期：424222",
      "\n实际：",
      ncol(obj)
    )
  )
}


if (nrow(
  obj
) != 18687L) {

  stop(
    paste0(
      "对象基因数异常。",
      "\n预期：18687",
      "\n实际：",
      nrow(obj)
    )
  )
}


###############################################################################
# 10. RNA assay
###############################################################################

if (!ASSAY_USE %in%
    Assays(
      obj
    )) {

  stop(
    paste0(
      "对象不存在 assay：",
      ASSAY_USE
    )
  )
}


DefaultAssay(
  obj
) <- ASSAY_USE


###############################################################################
# 11. RNA layer
###############################################################################

required_layers <- c(
  "counts",
  "data"
)


available_layers <- Layers(
  obj[[ASSAY_USE]]
)


missing_layers <- setdiff(
  required_layers,
  available_layers
)


if (length(
  missing_layers
) > 0L) {

  stop(
    paste0(
      "RNA assay 缺少 layer：",
      paste(
        missing_layers,
        collapse = ", "
      )
    )
  )
}


###############################################################################
# 12. Reduction
###############################################################################

if (!REDUCTION_USE %in%
    Reductions(
      obj
    )) {

  stop(
    paste0(
      "对象不存在 reduction：",
      REDUCTION_USE
    )
  )
}


###############################################################################
# 13. Metadata
###############################################################################

required_metadata <- c(
  "CellID",
  "DatasetID",
  "SampleID_Global",
  "Disease",
  CLUSTER_COL,
  "seurat_clusters"
)


missing_metadata <- setdiff(
  required_metadata,
  colnames(
    obj[[]]
  )
)


if (length(
  missing_metadata
) > 0L) {

  stop(
    paste0(
      "缺少必要 metadata：",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  )
}


###############################################################################
# 14. 正式 cluster 一致性
###############################################################################

if (!identical(

  as.character(
    obj$seurat_clusters
  ),

  as.character(
    obj[[CLUSTER_COL]][, 1]
  )

)) {

  stop(
    "seurat_clusters 与 RNA_snn_res.0.6 不一致。"
  )
}


Idents(
  obj
) <- CLUSTER_COL


###############################################################################
# 15. CellID
###############################################################################

md <- data.table::as.data.table(
  obj[[]]
)


if (anyDuplicated(
  md$CellID
) > 0L) {

  stop(
    "CellID 不是全局唯一。"
  )
}


if (!identical(

  as.character(
    md$CellID
  ),

  colnames(
    obj
  )

)) {

  stop(
    "metadata$CellID 与 colnames(obj) 未严格一致。"
  )
}


###############################################################################
# 16. UMAP embeddings
###############################################################################

umap_matrix <- Seurat::Embeddings(
  object = obj,
  reduction = REDUCTION_USE
)


if (ncol(
  umap_matrix
) < 2L) {

  stop(
    "umap_harmony16 少于2个维度。"
  )
}


if (!identical(

  rownames(
    umap_matrix
  ),

  colnames(
    obj
  )

)) {

  stop(
    "UMAP embeddings 与 cell 顺序不一致。"
  )
}


###############################################################################
# 17. RNA gene universe
###############################################################################

RNA_GENES <- rownames(
  obj[[ASSAY_USE]]
)


if (length(
  RNA_GENES
) != 18687L) {

  stop(
    "RNA gene universe 数量异常。"
  )
}


###############################################################################
###############################################################################
#                       PART C
#                  GENE MAPPING AUDIT
###############################################################################
###############################################################################


###############################################################################
# 18. 与正式 RNA assay 取交集
#
# 不要求全部存在。
# missing gene 会明确保存。
###############################################################################

GENE_AUDIT <- data.table::copy(
  SIGNATURE_REFERENCE
)


GENE_AUDIT[
  ,
  present_in_RNA :=
    human_gene %in%
    RNA_GENES
]


GENE_AUDIT[
  ,
  used_in_AddModuleScore :=
    present_in_RNA
]


GENE_AUDIT[
  ,
  status :=
    ifelse(
      present_in_RNA,
      "present",
      "absent"
    )
]


###############################################################################
# 19. 保存 gene audit
###############################################################################

GENE_AUDIT_FILE <- file.path(
  DATA_DIR,
  "22_Fig2fg_microglia_M1M2_gene_audit.csv"
)


data.table::fwrite(
  GENE_AUDIT,
  GENE_AUDIT_FILE
)


data.table::fwrite(
  GENE_AUDIT,
  file.path(
    TABLE_DIR,
    "Fig2fg_microglia_M1M2_gene_audit.csv"
  )
)


###############################################################################
# 20. 打印 mapping
###############################################################################

message("")
message("============================================================")
message("Microglia M1-like / M2-like gene mapping")
message("============================================================")


print(
  GENE_AUDIT[
    ,
    .(
      signature,
      reference_gene,
      human_gene,
      present_in_RNA,
      used_in_AddModuleScore,
      mapping_type
    )
  ]
)


###############################################################################
# 21. 正式用于打分的 genes
###############################################################################

M1_GENES_FINAL <- GENE_AUDIT[
  signature == "M1-like" &
    used_in_AddModuleScore == TRUE,
  human_gene
]


M2_GENES_FINAL <- GENE_AUDIT[
  signature == "M2-like" &
    used_in_AddModuleScore == TRUE,
  human_gene
]


###############################################################################
# Missing
###############################################################################

M1_GENES_MISSING <- GENE_AUDIT[
  signature == "M1-like" &
    present_in_RNA == FALSE,
  human_gene
]


M2_GENES_MISSING <- GENE_AUDIT[
  signature == "M2-like" &
    present_in_RNA == FALSE,
  human_gene
]


###############################################################################
# 22. Retention
###############################################################################

M1_RETENTION <-
  length(
    M1_GENES_FINAL
  ) /
  length(
    M1_HUMAN_CANDIDATES
  )


M2_RETENTION <-
  length(
    M2_GENES_FINAL
  ) /
  length(
    M2_HUMAN_CANDIDATES
  )


###############################################################################
# Summary
###############################################################################

GENESET_SUMMARY <- data.table::data.table(

  signature = c(
    "M1-like",
    "M2-like"
  ),

  source_genes = c(
    length(
      M1_HUMAN_CANDIDATES
    ),
    length(
      M2_HUMAN_CANDIDATES
    )
  ),

  genes_present = c(
    length(
      M1_GENES_FINAL
    ),
    length(
      M2_GENES_FINAL
    )
  ),

  genes_missing = c(
    length(
      M1_GENES_MISSING
    ),
    length(
      M2_GENES_MISSING
    )
  ),

  retention_fraction = c(
    M1_RETENTION,
    M2_RETENTION
  ),

  retention_percent = c(
    M1_RETENTION * 100,
    M2_RETENTION * 100
  ),

  genes_used = c(
    paste(
      M1_GENES_FINAL,
      collapse = ";"
    ),
    paste(
      M2_GENES_FINAL,
      collapse = ";"
    )
  ),

  missing_genes = c(
    paste(
      M1_GENES_MISSING,
      collapse = ";"
    ),
    paste(
      M2_GENES_MISSING,
      collapse = ";"
    )
  )
)


data.table::fwrite(
  GENESET_SUMMARY,
  file.path(
    TABLE_DIR,
    "Fig2fg_microglia_M1M2_geneset_summary.csv"
  )
)


message("")
print(
  GENESET_SUMMARY
)


###############################################################################
# 23. 最低 gene-set QC
###############################################################################

if (
  length(
    M1_GENES_FINAL
  ) <
  MIN_GENES_REQUIRED
) {

  stop(
    paste0(
      "M1-like gene set 保留基因过少。",
      "\n要求至少：",
      MIN_GENES_REQUIRED,
      "\n实际：",
      length(
        M1_GENES_FINAL
      ),
      "\nMissing：",
      paste(
        M1_GENES_MISSING,
        collapse = ", "
      )
    )
  )
}


if (
  length(
    M2_GENES_FINAL
  ) <
  MIN_GENES_REQUIRED
) {

  stop(
    paste0(
      "M2-like gene set 保留基因过少。",
      "\n要求至少：",
      MIN_GENES_REQUIRED,
      "\n实际：",
      length(
        M2_GENES_FINAL
      ),
      "\nMissing：",
      paste(
        M2_GENES_MISSING,
        collapse = ", "
      )
    )
  )
}


if (
  M1_RETENTION <
  MIN_RETENTION_REQUIRED
) {

  stop(
    paste0(
      "M1-like retention不足：",
      round(
        M1_RETENTION * 100,
        1
      ),
      "%"
    )
  )
}


if (
  M2_RETENTION <
  MIN_RETENTION_REQUIRED
) {

  stop(
    paste0(
      "M2-like retention不足：",
      round(
        M2_RETENTION * 100,
        1
      ),
      "%"
    )
  )
}


###############################################################################
# 24. 去重检查
###############################################################################

if (anyDuplicated(
  M1_GENES_FINAL
) > 0L) {

  stop(
    "M1-like final genes 存在重复。"
  )
}


if (anyDuplicated(
  M2_GENES_FINAL
) > 0L) {

  stop(
    "M2-like final genes 存在重复。"
  )
}


###############################################################################
# 25. 检查两个 gene set 是否存在交集
###############################################################################

M1_M2_OVERLAP <- intersect(
  M1_GENES_FINAL,
  M2_GENES_FINAL
)


if (length(
  M1_M2_OVERLAP
) > 0L) {

  warning(
    paste0(
      "M1-like 与 M2-like gene set 存在重叠基因：",
      paste(
        M1_M2_OVERLAP,
        collapse = ", "
      )
    )
  )
}


###############################################################################
# 26. 保存最终 gene set
###############################################################################

FINAL_SIGNATURE_TABLE <- data.table::rbindlist(
  list(

    data.table::data.table(

      signature =
        "M1-like",

      gene =
        M1_GENES_FINAL
    ),

    data.table::data.table(

      signature =
        "M2-like",

      gene =
        M2_GENES_FINAL
    )

  )
)


FINAL_SIGNATURE_TABLE[
  ,
  reference :=
    "Cell Death & Disease 2022"
]


FINAL_SIGNATURE_TABLE[
  ,
  DOI :=
    "10.1038/s41419-022-04864-z"
]


data.table::fwrite(
  FINAL_SIGNATURE_TABLE,
  file.path(
    DATA_DIR,
    "23_Fig2fg_microglia_M1M2_final_signatures.csv"
  )
)


###############################################################################
# Console
###############################################################################

message("")
message("最终用于 AddModuleScore：")

message(
  "M1-like: ",
  paste(
    M1_GENES_FINAL,
    collapse = ", "
  )
)


message(
  "M2-like: ",
  paste(
    M2_GENES_FINAL,
    collapse = ", "
  )
)


if (
  length(
    M1_GENES_MISSING
  ) > 0L
) {

  message(
    "M1-like missing: ",
    paste(
      M1_GENES_MISSING,
      collapse = ", "
    )
  )
}


if (
  length(
    M2_GENES_MISSING
  ) > 0L
) {

  message(
    "M2-like missing: ",
    paste(
      M2_GENES_MISSING,
      collapse = ", "
    )
  )
}


###############################################################################
###############################################################################
#                       PART D
#                   ADD MODULE SCORE
###############################################################################
###############################################################################


###############################################################################
# 27. 清除当前内存中可能存在的同名临时列
#
# 不修改硬盘上的 merged_obj.rds
###############################################################################

temporary_columns <- c(
  "M1Microglia1",
  "M2Microglia1",
  "M1_like_score",
  "M2_like_score"
)


for (col_i in temporary_columns) {

  if (
    col_i %in%
    colnames(
      obj[[]]
    )
  ) {

    obj[[col_i]] <- NULL
  }
}


###############################################################################
# 28. M1-like AddModuleScore
###############################################################################

message("")
message("============================================================")
message("计算 M1-like AddModuleScore")
message("============================================================")


set.seed(
  SEED
)


obj <- Seurat::AddModuleScore(

  object =
    obj,

  features =
    list(
      M1_GENES_FINAL
    ),

  pool =
    RNA_GENES,

  nbin =
    MODULE_NBIN,

  ctrl =
    MODULE_CTRL,

  assay =
    ASSAY_USE,

  name =
    "M1Microglia",

  seed =
    SEED,

  search =
    FALSE
)


###############################################################################
# 应生成 M1Microglia1
###############################################################################

if (!"M1Microglia1" %in%
    colnames(
      obj[[]]
    )) {

  stop(
    "AddModuleScore 未生成 M1Microglia1。"
  )
}


obj$M1_like_score <-
  obj$M1Microglia1


###############################################################################
# 29. M2-like AddModuleScore
###############################################################################

message("")
message("============================================================")
message("计算 M2-like AddModuleScore")
message("============================================================")


set.seed(
  SEED
)


obj <- Seurat::AddModuleScore(

  object =
    obj,

  features =
    list(
      M2_GENES_FINAL
    ),

  pool =
    RNA_GENES,

  nbin =
    MODULE_NBIN,

  ctrl =
    MODULE_CTRL,

  assay =
    ASSAY_USE,

  name =
    "M2Microglia",

  seed =
    SEED,

  search =
    FALSE
)


###############################################################################
# 应生成 M2Microglia1
###############################################################################

if (!"M2Microglia1" %in%
    colnames(
      obj[[]]
    )) {

  stop(
    "AddModuleScore 未生成 M2Microglia1。"
  )
}


obj$M2_like_score <-
  obj$M2Microglia1


###############################################################################
# 30. Score QC
###############################################################################

if (
  anyNA(
    obj$M1_like_score
  )
) {

  stop(
    "M1_like_score 存在 NA。"
  )
}


if (
  anyNA(
    obj$M2_like_score
  )
) {

  stop(
    "M2_like_score 存在 NA。"
  )
}


if (
  !all(
    is.finite(
      obj$M1_like_score
    )
  )
) {

  stop(
    "M1_like_score 存在非有限值。"
  )
}


if (
  !all(
    is.finite(
      obj$M2_like_score
    )
  )
) {

  stop(
    "M2_like_score 存在非有限值。"
  )
}


###############################################################################
###############################################################################
#                       PART E
#                  CELL-LEVEL DATA
###############################################################################
###############################################################################


###############################################################################
# 31. Cluster label
###############################################################################

cluster_numeric <- as.integer(
  as.character(
    obj[[CLUSTER_COL]][, 1]
  )
)


if (
  anyNA(
    cluster_numeric
  )
) {

  stop(
    "无法将 cluster 转换为整数。"
  )
}


cluster_label <- sprintf(
  "C%02d",
  cluster_numeric
)


###############################################################################
# 32. Cell-level score table
###############################################################################

score_data <- data.table::data.table(

  CellID =
    colnames(
      obj
    ),

  UMAP_1 =
    umap_matrix[
      ,
      1
    ],

  UMAP_2 =
    umap_matrix[
      ,
      2
    ],

  M1_like_score =
    as.numeric(
      obj$M1_like_score
    ),

  M2_like_score =
    as.numeric(
      obj$M2_like_score
    ),

  Cluster =
    cluster_label,

  Disease =
    as.character(
      obj$Disease
    ),

  DatasetID =
    as.character(
      obj$DatasetID
    ),

  SampleID_Global =
    as.character(
      obj$SampleID_Global
    )
)


###############################################################################
# Cell number
###############################################################################

if (
  nrow(
    score_data
  ) != 424222L
) {

  stop(
    paste0(
      "score_data 行数异常：",
      nrow(
        score_data
      )
    )
  )
}


###############################################################################
# ID顺序
###############################################################################

if (!identical(

  score_data$CellID,

  colnames(
    obj
  )

)) {

  stop(
    "score_data CellID 顺序异常。"
  )
}


###############################################################################
# 33. 计算 plotting limits
#
# M1和M2分别使用自身1%-99%范围。
#
# 原始score不修改。
###############################################################################

M1_LIMITS <- as.numeric(
  stats::quantile(
    score_data$M1_like_score,
    probs = c(
      PLOT_LOW_QUANTILE,
      PLOT_HIGH_QUANTILE
    ),
    na.rm = TRUE,
    names = FALSE
  )
)


M2_LIMITS <- as.numeric(
  stats::quantile(
    score_data$M2_like_score,
    probs = c(
      PLOT_LOW_QUANTILE,
      PLOT_HIGH_QUANTILE
    ),
    na.rm = TRUE,
    names = FALSE
  )
)


###############################################################################
# 防止异常
###############################################################################

if (
  !all(
    is.finite(
      M1_LIMITS
    )
  ) ||
  M1_LIMITS[1] >=
  M1_LIMITS[2]
) {

  stop(
    "M1-like UMAP color limits 计算异常。"
  )
}


if (
  !all(
    is.finite(
      M2_LIMITS
    )
  ) ||
  M2_LIMITS[1] >=
  M2_LIMITS[2]
) {

  stop(
    "M2-like UMAP color limits 计算异常。"
  )
}


###############################################################################
# 34. Plot clipping
#
# 只影响画图。
###############################################################################

score_data[
  ,
  M1_like_score_plot :=
    pmax(
      M1_LIMITS[1],
      pmin(
        M1_LIMITS[2],
        M1_like_score
      )
    )
]


score_data[
  ,
  M2_like_score_plot :=
    pmax(
      M2_LIMITS[1],
      pmin(
        M2_LIMITS[2],
        M2_like_score
      )
    )
]


###############################################################################
# 35. 保存 cell-level score
#
# 后续 Fig2H 可以直接使用，不用重新 AddModuleScore
###############################################################################

data.table::fwrite(
  score_data,
  file.path(
    DATA_DIR,
    "24_Fig2fg_microglia_cell_M1M2_scores.csv"
  )
)


saveRDS(
  score_data,
  file.path(
    DATA_DIR,
    "25_Fig2fg_microglia_cell_M1M2_scores.rds"
  )
)


###############################################################################
# 36. Score distribution audit
###############################################################################

score_distribution <- data.table::rbindlist(
  list(

    data.table::data.table(

      signature =
        "M1-like",

      n_cells =
        length(
          score_data$M1_like_score
        ),

      min =
        min(
          score_data$M1_like_score
        ),

      q01 =
        as.numeric(
          stats::quantile(
            score_data$M1_like_score,
            0.01
          )
        ),

      q05 =
        as.numeric(
          stats::quantile(
            score_data$M1_like_score,
            0.05
          )
        ),

      median =
        stats::median(
          score_data$M1_like_score
        ),

      mean =
        mean(
          score_data$M1_like_score
        ),

      q95 =
        as.numeric(
          stats::quantile(
            score_data$M1_like_score,
            0.95
          )
        ),

      q99 =
        as.numeric(
          stats::quantile(
            score_data$M1_like_score,
            0.99
          )
        ),

      max =
        max(
          score_data$M1_like_score
        ),

      plot_min =
        M1_LIMITS[1],

      plot_max =
        M1_LIMITS[2]
    ),


    data.table::data.table(

      signature =
        "M2-like",

      n_cells =
        length(
          score_data$M2_like_score
        ),

      min =
        min(
          score_data$M2_like_score
        ),

      q01 =
        as.numeric(
          stats::quantile(
            score_data$M2_like_score,
            0.01
          )
        ),

      q05 =
        as.numeric(
          stats::quantile(
            score_data$M2_like_score,
            0.05
          )
        ),

      median =
        stats::median(
          score_data$M2_like_score
        ),

      mean =
        mean(
          score_data$M2_like_score
        ),

      q95 =
        as.numeric(
          stats::quantile(
            score_data$M2_like_score,
            0.95
          )
        ),

      q99 =
        as.numeric(
          stats::quantile(
            score_data$M2_like_score,
            0.99
          )
        ),

      max =
        max(
          score_data$M2_like_score
        ),

      plot_min =
        M2_LIMITS[1],

      plot_max =
        M2_LIMITS[2]
    )
  )
)


data.table::fwrite(
  score_distribution,
  file.path(
    TABLE_DIR,
    "Fig2fg_microglia_score_distribution.csv"
  )
)


###############################################################################
# 37. Cluster-level descriptive summary
#
# 为后续 Fig2H 提前保存
###############################################################################

cluster_score_summary <- score_data[
  ,
  .(

    n_cells =
      .N,

    M1_mean =
      mean(
        M1_like_score
      ),

    M1_median =
      stats::median(
        M1_like_score
      ),

    M1_sd =
      stats::sd(
        M1_like_score
      ),

    M1_q25 =
      as.numeric(
        stats::quantile(
          M1_like_score,
          0.25
        )
      ),

    M1_q75 =
      as.numeric(
        stats::quantile(
          M1_like_score,
          0.75
        )
      ),


    M2_mean =
      mean(
        M2_like_score
      ),

    M2_median =
      stats::median(
        M2_like_score
      ),

    M2_sd =
      stats::sd(
        M2_like_score
      ),

    M2_q25 =
      as.numeric(
        stats::quantile(
          M2_like_score,
          0.25
        )
      ),

    M2_q75 =
      as.numeric(
        stats::quantile(
          M2_like_score,
          0.75
        )
      )

  ),
  by = Cluster
]


cluster_score_summary[
  ,
  cluster_number :=
    as.integer(
      sub(
        "^C",
        "",
        Cluster
      )
    )
]


data.table::setorder(
  cluster_score_summary,
  cluster_number
)


cluster_score_summary[
  ,
  cluster_number := NULL
]


data.table::fwrite(
  cluster_score_summary,
  file.path(
    DATA_DIR,
    "26_Fig2fg_microglia_cluster_M1M2_summary.csv"
  )
)


###############################################################################
###############################################################################
#                         PART F
#                       UMAP PLOTS
###############################################################################
###############################################################################


###############################################################################
# 38. 为避免高score细胞被低score覆盖：
#
# 从低score -> 高score排序
###############################################################################

M1_PLOT_DATA <- data.table::copy(
  score_data
)


data.table::setorder(
  M1_PLOT_DATA,
  M1_like_score_plot
)


M2_PLOT_DATA <- data.table::copy(
  score_data
)


data.table::setorder(
  M2_PLOT_DATA,
  M2_like_score_plot
)


###############################################################################
# 39. SCI UMAP theme
###############################################################################

theme_score_umap <- ggplot2::theme_classic(
  base_family = FONT,
  base_size = 9
) +

  ggplot2::theme(

    axis.text =
      ggplot2::element_blank(),

    axis.ticks =
      ggplot2::element_blank(),

    axis.line =
      ggplot2::element_line(
        colour = "black",
        linewidth = 0.35
      ),

    axis.title =
      ggplot2::element_text(
        size = 8,
        colour = "black"
      ),

    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 10,
        hjust = 0,
        colour = "black"
      ),

    plot.subtitle =
      ggplot2::element_text(
        size = 7,
        colour = "grey30"
      ),

    legend.title =
      ggplot2::element_text(
        size = 8,
        colour = "black"
      ),

    legend.text =
      ggplot2::element_text(
        size = 7,
        colour = "black"
      ),

    legend.key.height =
      grid::unit(
        0.65,
        "cm"
      ),

    plot.background =
      ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),

    panel.background =
      ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),

    aspect.ratio =
      1
  )


###############################################################################
# 40. Fig2F
#
# M1-like score
###############################################################################

p_M1 <-

  ggplot2::ggplot(

    M1_PLOT_DATA,

    ggplot2::aes(
      x =
        UMAP_1,

      y =
        UMAP_2,

      colour =
        M1_like_score_plot
    )
  ) +

  ggplot2::geom_point(
    size = 0.10,
    alpha = 0.85,
    shape = 16,
    stroke = 0
  ) +

  ggplot2::scale_colour_gradient(

    name =
      "M1-like\nscore",

    low =
      M1_LOW_COLOR,

    high =
      M1_HIGH_COLOR,

    limits =
      M1_LIMITS,

    breaks = c(
      M1_LIMITS[1],
      mean(
        M1_LIMITS
      ),
      M1_LIMITS[2]
    ),

    labels =
      round(
        c(
          M1_LIMITS[1],
          mean(
            M1_LIMITS
          ),
          M1_LIMITS[2]
        ),
        2
      )
  ) +

  ggplot2::coord_equal() +

  ggplot2::labs(

    title =
      "M1-like score",

    x =
      "UMAP 1",

    y =
      "UMAP 2"
  ) +

  theme_score_umap


###############################################################################
# 41. Fig2G
#
# M2-like score
###############################################################################

p_M2 <-

  ggplot2::ggplot(

    M2_PLOT_DATA,

    ggplot2::aes(
      x =
        UMAP_1,

      y =
        UMAP_2,

      colour =
        M2_like_score_plot
    )
  ) +

  ggplot2::geom_point(
    size = 0.10,
    alpha = 0.85,
    shape = 16,
    stroke = 0
  ) +

  ggplot2::scale_colour_gradient(

    name =
      "M2-like\nscore",

    low =
      M2_LOW_COLOR,

    high =
      M2_HIGH_COLOR,

    limits =
      M2_LIMITS,

    breaks = c(
      M2_LIMITS[1],
      mean(
        M2_LIMITS
      ),
      M2_LIMITS[2]
    ),

    labels =
      round(
        c(
          M2_LIMITS[1],
          mean(
            M2_LIMITS
          ),
          M2_LIMITS[2]
        ),
        2
      )
  ) +

  ggplot2::coord_equal() +

  ggplot2::labs(

    title =
      "M2-like score",

    x =
      "UMAP 1",

    y =
      "UMAP 2"
  ) +

  theme_score_umap


###############################################################################
# 42. Source CSV
###############################################################################

data.table::fwrite(

  M1_PLOT_DATA[
    ,
    .(
      CellID,
      UMAP_1,
      UMAP_2,
      M1_like_score,
      M1_like_score_plot,
      Cluster,
      Disease,
      DatasetID,
      SampleID_Global
    )
  ],

  file.path(
    SOURCE_DIR,
    "Fig2f_microglia_M1like_UMAP_source.csv"
  )
)


data.table::fwrite(

  M2_PLOT_DATA[
    ,
    .(
      CellID,
      UMAP_1,
      UMAP_2,
      M2_like_score,
      M2_like_score_plot,
      Cluster,
      Disease,
      DatasetID,
      SampleID_Global
    )
  ],

  file.path(
    SOURCE_DIR,
    "Fig2g_microglia_M2like_UMAP_source.csv"
  )
)


###############################################################################
# 43. Plot保存函数
###############################################################################

save_formal_plot <- function(
    plot_object,
    basename,
    width = 5.5,
    height = 5.2
) {

  png_file <- file.path(
    FIG_DIR,
    paste0(
      basename,
      ".png"
    )
  )


  pdf_file <- file.path(
    FIG_DIR,
    paste0(
      basename,
      ".pdf"
    )
  )


  ###########################################################################
  # PNG
  ###########################################################################

  ggplot2::ggsave(
    filename =
      png_file,

    plot =
      plot_object,

    width =
      width,

    height =
      height,

    units =
      "in",

    dpi =
      600,

    bg =
      "white"
  )


  ###########################################################################
  # PDF
  ###########################################################################

  pdf_success <- FALSE


  if (
    capabilities(
      "cairo"
    )
  ) {

    pdf_success <- tryCatch(

      {

        ggplot2::ggsave(

          filename =
            pdf_file,

          plot =
            plot_object,

          width =
            width,

          height =
            height,

          units =
            "in",

          device =
            grDevices::cairo_pdf,

          family =
            FONT,

          bg =
            "white"
        )

        TRUE
      },

      error = function(e) {

        FALSE
      }
    )
  }


  if (!pdf_success) {

    warning(
      paste0(
        "Cairo PDF输出失败，使用sans fallback：",
        basename
      )
    )


    plot_fallback <-
      plot_object +

      ggplot2::theme(
        text =
          ggplot2::element_text(
            family = "sans"
          )
      )


    grDevices::pdf(
      file =
        pdf_file,

      width =
        width,

      height =
        height,

      family =
        "sans",

      useDingbats =
        FALSE
    )


    print(
      plot_fallback
    )


    grDevices::dev.off()
  }
}


###############################################################################
# 44. 保存 Fig2F
###############################################################################

save_formal_plot(

  plot_object =
    p_M1,

  basename =
    "Fig2f_microglia_M1like_AddModuleScore_UMAP",

  width =
    5.5,

  height =
    5.2
)


###############################################################################
# 45. 保存 Fig2G
###############################################################################

save_formal_plot(

  plot_object =
    p_M2,

  basename =
    "Fig2g_microglia_M2like_AddModuleScore_UMAP",

  width =
    5.5,

  height =
    5.2
)


###############################################################################
###############################################################################
#                       PART G
#                   METHOD / AUDIT
###############################################################################
###############################################################################


###############################################################################
# 46. 方法说明
###############################################################################

METHOD_NOTE <- c(

  "Fig2F / Fig2G microglia-adapted M1-like / M2-like scoring",

  "",

  "PRIMARY SIGNATURE REFERENCE",

  "Single-cell transcriptome analysis reveals the immune heterogeneity and the repopulation of microglia by Hif1a in mice after spinal cord injury.",

  "Cell Death & Disease, 2022.",

  "DOI: 10.1038/s41419-022-04864-z",

  "",

  "SOURCE MARKERS",

  "M1 markers reported in the reference: TNF, Il1b, Icam1, Fas, Ccl3, Cd24a.",

  "M2 markers reported in the reference: Itga6, Tgfb1, Arg1, Mrc1, Igf1, Ccl5.",

  "",

  "HUMAN ORTHOLOG SYMBOLS",

  paste0(
    "M1-like candidates: ",
    paste(
      M1_HUMAN_CANDIDATES,
      collapse = ", "
    )
  ),

  paste0(
    "M2-like candidates: ",
    paste(
      M2_HUMAN_CANDIDATES,
      collapse = ", "
    )
  ),

  "Mouse Cd24a was mapped explicitly to human CD24.",

  "",

  "RNA INTERSECTION",

  "Human ortholog gene sets were intersected with rownames of the formal RNA assay.",

  "Missing genes were recorded in the audit table and were not replaced by guessed aliases.",

  paste0(
    "M1-like genes retained: ",
    length(
      M1_GENES_FINAL
    ),
    "/",
    length(
      M1_HUMAN_CANDIDATES
    ),
    " (",
    round(
      M1_RETENTION * 100,
      1
    ),
    "%)"
  ),

  paste0(
    "M2-like genes retained: ",
    length(
      M2_GENES_FINAL
    ),
    "/",
    length(
      M2_HUMAN_CANDIDATES
    ),
    " (",
    round(
      M2_RETENTION * 100,
      1
    ),
    "%)"
  ),

  paste0(
    "Final M1-like genes: ",
    paste(
      M1_GENES_FINAL,
      collapse = ", "
    )
  ),

  paste0(
    "Final M2-like genes: ",
    paste(
      M2_GENES_FINAL,
      collapse = ", "
    )
  ),

  "",

  "SCORING",

  "Seurat AddModuleScore was used.",

  paste0(
    "Assay: ",
    ASSAY_USE
  ),

  paste0(
    "nbin: ",
    MODULE_NBIN
  ),

  paste0(
    "ctrl: ",
    MODULE_CTRL
  ),

  paste0(
    "seed: ",
    SEED
  ),

  "All 424,222 cells participated in scoring.",

  "",

  "UMAP",

  paste0(
    "Existing reduction used: ",
    REDUCTION_USE
  ),

  "No PCA, Harmony, UMAP, nearest-neighbor graph or clustering was recomputed.",

  "",

  "VISUALIZATION",

  paste0(
    "Plotting color range uses q",
    PLOT_LOW_QUANTILE * 100,
    " to q",
    PLOT_HIGH_QUANTILE * 100,
    " clipping."
  ),

  "The raw AddModuleScore values were not clipped and are preserved in the source and intermediate data.",

  "M1-like high score is shown in blue.",

  "M2-like high score is shown in red (#B2182B).",

  "",

  "INTERPRETATION",

  "M1-like and M2-like scores represent relative transcriptional programs.",

  "They should not be interpreted as binary M1/M2 cell identities.",

  "Microglial activation exists along heterogeneous and overlapping transcriptional states."
)


writeLines(
  METHOD_NOTE,
  file.path(
    TABLE_DIR,
    "Fig2fg_microglia_method_note.txt"
  )
)


###############################################################################
# 47. 数据字典
###############################################################################

DATA_DICTIONARY <- data.table::data.table(

  file = c(
    "22_Fig2fg_microglia_M1M2_gene_audit.csv",
    "23_Fig2fg_microglia_M1M2_final_signatures.csv",
    "24_Fig2fg_microglia_cell_M1M2_scores.csv",
    "25_Fig2fg_microglia_cell_M1M2_scores.rds",
    "26_Fig2fg_microglia_cluster_M1M2_summary.csv"
  ),

  purpose = c(

    "文献mouse marker到human ortholog以及是否存在于RNA assay的完整审计",

    "最终实际进入AddModuleScore的M1-like和M2-like genes",

    "424222个细胞的M1-like/M2-like score和UMAP坐标，供Agent及后续Fig2H使用",

    "与CSV相同的cell-level结果RDS，R中读取更快",

    "30个cluster的M1-like/M2-like描述性summary"
  )
)


data.table::fwrite(
  DATA_DICTIONARY,
  file.path(
    DATA_DIR,
    "27_Fig2fg_microglia_data_dictionary.csv"
  )
)


###############################################################################
# 48. sessionInfo
###############################################################################

capture.output(
  sessionInfo(),
  file = file.path(
    TABLE_DIR,
    "Fig2fg_microglia_sessionInfo.txt"
  )
)


###############################################################################
# 49. 最终控制台输出
###############################################################################

message("")
message("============================================================")
message("Fig2F / Fig2G 完成")
message("============================================================")


message(
  "总细胞数：",
  nrow(
    score_data
  )
)


message("")
message("M1-like")

message(
  "候选 genes：",
  length(
    M1_HUMAN_CANDIDATES
  )
)

message(
  "实际 genes：",
  length(
    M1_GENES_FINAL
  )
)

message(
  "Retention：",
  round(
    M1_RETENTION * 100,
    1
  ),
  "%"
)

message(
  "Used：",
  paste(
    M1_GENES_FINAL,
    collapse = ", "
  )
)


if (
  length(
    M1_GENES_MISSING
  ) > 0L
) {

  message(
    "Missing：",
    paste(
      M1_GENES_MISSING,
      collapse = ", "
    )
  )
}


message("")
message("M2-like")

message(
  "候选 genes：",
  length(
    M2_HUMAN_CANDIDATES
  )
)

message(
  "实际 genes：",
  length(
    M2_GENES_FINAL
  )
)

message(
  "Retention：",
  round(
    M2_RETENTION * 100,
    1
  ),
  "%"
)

message(
  "Used：",
  paste(
    M2_GENES_FINAL,
    collapse = ", "
  )
)


if (
  length(
    M2_GENES_MISSING
  ) > 0L
) {

  message(
    "Missing：",
    paste(
      M2_GENES_MISSING,
      collapse = ", "
    )
  )
}


message("")
message("Gene audit：")

message(
  GENE_AUDIT_FILE
)


message("")
message("Cell-level scores：")

message(
  file.path(
    DATA_DIR,
    "24_Fig2fg_microglia_cell_M1M2_scores.csv"
  )
)


message("")
message("Fig2F：")

message(
  file.path(
    FIG_DIR,
    "Fig2f_microglia_M1like_AddModuleScore_UMAP.png"
  )
)


message("")
message("Fig2G：")

message(
  file.path(
    FIG_DIR,
    "Fig2g_microglia_M2like_AddModuleScore_UMAP.png"
  )
)


message("")
message("============================================================")














###############################################################################
###############################################################################
#                       PART H
#                        FIG2H
#          M1-like / M2-like score by subcluster
###############################################################################
###############################################################################
#
# 说明：
#
# 1. 本段直接续接 Fig2F/G 脚本。
# 2. 不重新读取 Seurat object。
# 3. 不重新计算 AddModuleScore。
# 4. 直接使用已经存在于内存中的：
#
#       score_data
#
#    其中必须包含：
#
#       CellID
#       M1_like_score
#       M2_like_score
#       Cluster
#       Disease
#       DatasetID
#       SampleID_Global
#
# 5. Fig2H 主图按照原论文 Fig2h 的视觉逻辑：
#
#       每个 cluster 两个并排 boxplot
#       M1-like = blue
#       M2-like = red
#
# 6. 主图是 cell-level descriptive visualization。
#    不在这里对42万细胞直接进行显著性检验。
#
# 7. 另外保存 SampleID_Global × Cluster 的样本级平均score，
#    后续正式统计 / Agent解释优先使用该表。
#
###############################################################################


###############################################################################
# 50. 检查前序变量是否存在
###############################################################################

required_objects_fig2h <- c(
  "score_data",
  "FIG_DIR",
  "SOURCE_DIR",
  "TABLE_DIR",
  "DATA_DIR",
  "FONT",
  "M1_HIGH_COLOR",
  "M2_HIGH_COLOR"
)


missing_objects_fig2h <- required_objects_fig2h[
  !vapply(
    required_objects_fig2h,
    exists,
    logical(1),
    inherits = TRUE
  )
]


if (length(
  missing_objects_fig2h
) > 0L) {

  stop(
    paste0(
      "Fig2H缺少前序脚本对象：",
      paste(
        missing_objects_fig2h,
        collapse = ", "
      ),
      "\n请先完整运行Fig2F/G脚本，再续接本段代码。"
    )
  )
}


###############################################################################
# 51. score_data 必要列检查
###############################################################################

required_score_columns_fig2h <- c(
  "CellID",
  "M1_like_score",
  "M2_like_score",
  "Cluster",
  "Disease",
  "DatasetID",
  "SampleID_Global"
)


missing_score_columns_fig2h <- setdiff(
  required_score_columns_fig2h,
  colnames(
    score_data
  )
)


if (length(
  missing_score_columns_fig2h
) > 0L) {

  stop(
    paste0(
      "score_data缺少Fig2H必要列：",
      paste(
        missing_score_columns_fig2h,
        collapse = ", "
      )
    )
  )
}


###############################################################################
# 52. Cell-level score基本检查
###############################################################################

if (nrow(
  score_data
) != 424222L) {

  stop(
    paste0(
      "score_data细胞数异常。",
      "\n预期：424222",
      "\n实际：",
      nrow(
        score_data
      )
    )
  )
}


if (
  anyNA(
    score_data$M1_like_score
  ) ||
  anyNA(
    score_data$M2_like_score
  )
) {

  stop(
    "M1-like或M2-like score存在NA。"
  )
}


if (
  !all(
    is.finite(
      score_data$M1_like_score
    )
  ) ||
  !all(
    is.finite(
      score_data$M2_like_score
    )
  )
) {

  stop(
    "M1-like或M2-like score存在非有限值。"
  )
}


###############################################################################
# 53. Fig2H Cluster顺序
#
# 优先使用前面Fig2C/D/E冻结好的cluster行序：
#
# E:/micro_output/fig3/data/05_cluster_row_order.csv
#
# Fig2H从左到右使用相同顺序。
###############################################################################

FIG2H_CLUSTER_ORDER_FILE <- file.path(
  DATA_DIR,
  "05_cluster_row_order.csv"
)


if (!file.exists(
  FIG2H_CLUSTER_ORDER_FILE
)) {

  stop(
    paste0(
      "Fig2H要求与Fig2C/D/E保持cluster顺序一致，",
      "但未找到：\n",
      FIG2H_CLUSTER_ORDER_FILE
    )
  )
}


fig2h_cluster_order_table <- data.table::fread(
  FIG2H_CLUSTER_ORDER_FILE
)


###############################################################################
# 兼容 Cluster / cluster 两种列名
###############################################################################

if (
  "Cluster" %in%
  colnames(
    fig2h_cluster_order_table
  )
) {

  cluster_order_raw <-
    as.character(
      fig2h_cluster_order_table$Cluster
    )

} else if (
  "cluster" %in%
  colnames(
    fig2h_cluster_order_table
  )
) {

  cluster_order_raw <-
    as.character(
      fig2h_cluster_order_table$cluster
    )

} else {

  stop(
    "05_cluster_row_order.csv 中没有 Cluster 或 cluster 列。"
  )
}


###############################################################################
# 标准化为C00-C29
###############################################################################

normalize_cluster_label_fig2h <- function(x) {

  x <- as.character(
    x
  )


  result <- rep(
    NA_character_,
    length(
      x
    )
  )


  is_c_style <- grepl(
    "^C[0-9]+$",
    x
  )


  is_numeric_style <- grepl(
    "^[0-9]+$",
    x
  )


  result[
    is_c_style
  ] <- sprintf(
    "C%02d",
    as.integer(
      sub(
        "^C",
        "",
        x[
          is_c_style
        ]
      )
    )
  )


  result[
    is_numeric_style
  ] <- sprintf(
    "C%02d",
    as.integer(
      x[
        is_numeric_style
      ]
    )
  )


  if (
    anyNA(
      result
    )
  ) {

    stop(
      "Fig2H cluster order文件中存在无法识别的cluster标签。"
    )
  }


  return(
    result
  )
}


FIG2H_CLUSTER_ORDER <- normalize_cluster_label_fig2h(
  cluster_order_raw
)


###############################################################################
# 必须正好30个cluster
###############################################################################

if (length(
  FIG2H_CLUSTER_ORDER
) != 30L) {

  stop(
    paste0(
      "Fig2H cluster顺序不是30个cluster。",
      "\n实际：",
      length(
        FIG2H_CLUSTER_ORDER
      )
    )
  )
}


if (anyDuplicated(
  FIG2H_CLUSTER_ORDER
) > 0L) {

  stop(
    "Fig2H cluster顺序存在重复。"
  )
}


if (!setequal(

  unique(
    as.character(
      score_data$Cluster
    )
  ),

  FIG2H_CLUSTER_ORDER

)) {

  stop(
    "score_data中的cluster与Fig2C/D/E冻结顺序不一致。"
  )
}


###############################################################################
# 54. 构建cell-level long table
#
# 一个细胞变成两行：
#
# CellID | Cluster | Signature | Score
#
# M1-like
# M2-like
###############################################################################

fig2h_cell_long <- data.table::melt(

  data.table::copy(
    score_data
  ),

  id.vars = c(
    "CellID",
    "Cluster",
    "Disease",
    "DatasetID",
    "SampleID_Global"
  ),

  measure.vars = c(
    "M1_like_score",
    "M2_like_score"
  ),

  variable.name =
    "Signature_raw",

  value.name =
    "Score",

  variable.factor =
    FALSE
)


###############################################################################
# 55. Signature显示名称
###############################################################################

fig2h_cell_long[
  ,
  Signature :=
    data.table::fcase(

      Signature_raw ==
        "M1_like_score",
      "M1-like",

      Signature_raw ==
        "M2_like_score",
      "M2-like",

      default =
        NA_character_
    )
]


if (
  anyNA(
    fig2h_cell_long$Signature
  )
) {

  stop(
    "Fig2H Signature转换失败。"
  )
}


fig2h_cell_long[
  ,
  Signature_raw := NULL
]


###############################################################################
# 因子顺序
###############################################################################

fig2h_cell_long[
  ,
  Cluster :=
    factor(
      Cluster,
      levels =
        FIG2H_CLUSTER_ORDER
    )
]


fig2h_cell_long[
  ,
  Signature :=
    factor(
      Signature,
      levels = c(
        "M1-like",
        "M2-like"
      )
    )
]


###############################################################################
# 56. 保存完整cell-level long table
#
# 848444行。
#
# CSV会比较大，因此data目录保存RDS；
# 已有的24_Fig2fg_microglia_cell_M1M2_scores.csv
# 本身就是完整raw cell-level source，不需要再重复写一个超大CSV。
###############################################################################

saveRDS(
  fig2h_cell_long,
  file.path(
    DATA_DIR,
    "28_Fig2h_cell_scores_long.rds"
  )
)


###############################################################################
# 57. 计算Fig2H正式boxplot summary
#
# 这个CSV可以作为正式source_data。
###############################################################################

fig2h_boxplot_summary <- fig2h_cell_long[
  ,
  .(

    n_cells =
      .N,

    mean_score =
      mean(
        Score
      ),

    sd_score =
      stats::sd(
        Score
      ),

    min_score =
      min(
        Score
      ),

    q05 =
      as.numeric(
        stats::quantile(
          Score,
          0.05
        )
      ),

    q25 =
      as.numeric(
        stats::quantile(
          Score,
          0.25
        )
      ),

    median_score =
      stats::median(
        Score
      ),

    q75 =
      as.numeric(
        stats::quantile(
          Score,
          0.75
        )
      ),

    q95 =
      as.numeric(
        stats::quantile(
          Score,
          0.95
        )
      ),

    max_score =
      max(
        Score
      )

  ),
  by = .(
    Cluster,
    Signature
  )
]


###############################################################################
# IQR
###############################################################################

fig2h_boxplot_summary[
  ,
  IQR :=
    q75 -
    q25
]


###############################################################################
# ggplot boxplot whiskers
#
# lower whisker:
# max(min, Q1 - 1.5*IQR)
#
# upper whisker:
# min(max, Q3 + 1.5*IQR)
###############################################################################

fig2h_boxplot_summary[
  ,
  lower_whisker :=
    pmax(
      min_score,
      q25 -
        1.5 *
        IQR
    )
]


fig2h_boxplot_summary[
  ,
  upper_whisker :=
    pmin(
      max_score,
      q75 +
        1.5 *
        IQR
    )
]


###############################################################################
# 正式顺序
###############################################################################

fig2h_boxplot_summary[
  ,
  cluster_order :=
    match(
      as.character(
        Cluster
      ),
      FIG2H_CLUSTER_ORDER
    )
]


fig2h_boxplot_summary[
  ,
  signature_order :=
    match(
      as.character(
        Signature
      ),
      c(
        "M1-like",
        "M2-like"
      )
    )
]


data.table::setorder(
  fig2h_boxplot_summary,
  cluster_order,
  signature_order
)


###############################################################################
# Source CSV
###############################################################################

data.table::fwrite(
  fig2h_boxplot_summary,
  file.path(
    SOURCE_DIR,
    "Fig2h_M1like_M2like_boxplot_source.csv"
  )
)


###############################################################################
# 同时保存到data
###############################################################################

data.table::fwrite(
  fig2h_boxplot_summary,
  file.path(
    DATA_DIR,
    "29_Fig2h_boxplot_summary.csv"
  )
)


###############################################################################
###############################################################################
#             PART H2
#     SAMPLE-LEVEL DATA FOR FORMAL INTERPRETATION
###############################################################################
###############################################################################
#
# 原论文Fig2h是cell-level描述性boxplot。
#
# 但是如果后续需要做正式疾病/cluster统计，
# 不应该把42万细胞视为独立生物学重复。
#
# 因此这里额外生成：
#
# SampleID_Global × Cluster
#
# 的平均M1-like / M2-like score。
###############################################################################


###############################################################################
# 58. 每个样本总细胞数
###############################################################################

sample_cell_count_fig2h <- score_data[
  ,
  .(
    sample_total_cells =
      .N
  ),
  by =
    SampleID_Global
]


###############################################################################
# >=20 cells的正式样本
###############################################################################

sample_cell_count_fig2h[
  ,
  sample_eligible :=
    sample_total_cells >=
    20L
]


###############################################################################
# 59. Sample × Cluster summary
###############################################################################

fig2h_sample_cluster <- score_data[
  ,
  .(

    n_cells_cluster =
      .N,

    M1_like_mean =
      mean(
        M1_like_score
      ),

    M1_like_median =
      stats::median(
        M1_like_score
      ),

    M2_like_mean =
      mean(
        M2_like_score
      ),

    M2_like_median =
      stats::median(
        M2_like_score
      )

  ),
  by = .(
    SampleID_Global,
    DatasetID,
    Disease,
    Cluster
  )
]


###############################################################################
# 合并样本总细胞QC
###############################################################################

fig2h_sample_cluster <- merge(

  fig2h_sample_cluster,

  sample_cell_count_fig2h,

  by =
    "SampleID_Global",

  all.x =
    TRUE,

  sort =
    FALSE
)


###############################################################################
# Cluster order
###############################################################################

fig2h_sample_cluster[
  ,
  cluster_order :=
    match(
      Cluster,
      FIG2H_CLUSTER_ORDER
    )
]


data.table::setorder(
  fig2h_sample_cluster,
  cluster_order,
  Disease,
  SampleID_Global
)


###############################################################################
# 保存
###############################################################################

data.table::fwrite(
  fig2h_sample_cluster,
  file.path(
    DATA_DIR,
    "30_Fig2h_sample_cluster_M1M2_scores.csv"
  )
)


###############################################################################
# 60. 仅>=20 total cells样本版本
###############################################################################

fig2h_sample_cluster_eligible <- fig2h_sample_cluster[
  sample_eligible ==
    TRUE
]


data.table::fwrite(
  fig2h_sample_cluster_eligible,
  file.path(
    DATA_DIR,
    "31_Fig2h_sample_cluster_M1M2_scores_eligible_samples.csv"
  )
)


###############################################################################
# 61. Fig2H颜色
#
# 与原论文视觉逻辑保持一致：
#
# M1 = Blue
# M2 = Red
#
# 使用前序Fig2F/G已经定义好的SCI颜色。
###############################################################################

FIG2H_COLORS <- c(
  "M1-like" = "#2166AC",
  "M2-like" = "#B2182B"
)
###############################################################################
# 62. Fig2H主题
###############################################################################

theme_fig2h <- ggplot2::theme_classic(
  base_family =
    FONT,

  base_size =
    9
) +

  ggplot2::theme(

    ###########################################################################
    # X轴
    ###########################################################################

    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 7,
        colour = "black"
      ),


    ###########################################################################
    # Y轴
    ###########################################################################

    axis.text.y =
      ggplot2::element_text(
        size = 7,
        colour = "black"
      ),


    axis.title.x =
      ggplot2::element_blank(),


    axis.title.y =
      ggplot2::element_text(
        size = 8.5,
        colour = "black"
      ),


    ###########################################################################
    # 坐标轴线
    ###########################################################################

    axis.line =
      ggplot2::element_line(
        colour = "black",
        linewidth = 0.35
      ),


    axis.ticks =
      ggplot2::element_line(
        colour = "black",
        linewidth = 0.30
      ),


    ###########################################################################
    # 标题
    ###########################################################################

    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 10,
        hjust = 0,
        colour = "black"
      ),


    plot.subtitle =
      ggplot2::element_text(
        size = 7.5,
        colour = "grey30"
      ),


    ###########################################################################
    # Legend
    ###########################################################################

    legend.title =
      ggplot2::element_blank(),


    legend.text =
      ggplot2::element_text(
        size = 8,
        colour = "black"
      ),


    legend.position =
      "top",


    legend.justification =
      "left",


    legend.box.just =
      "left",


    ###########################################################################
    # 背景
    ###########################################################################

    panel.grid =
      ggplot2::element_blank(),


    panel.background =
      ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),


    plot.background =
      ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),


    plot.margin =
      ggplot2::margin(
        6,
        8,
        5,
        5
      )
  )


###############################################################################
# 63. Fig2H cell-level descriptive boxplot
#
# outlier.shape = NA：
# 42万细胞时不绘制海量离群点，
# 避免黑色outlier覆盖主体。
#
# 数据本身没有删除。
###############################################################################

BOX_POSITION <- ggplot2::position_dodge(
  width = 0.72
)


p_fig2h <-

  ggplot2::ggplot(

    fig2h_cell_long,

    ggplot2::aes(
      x =
        Cluster,

      y =
        Score,

      fill =
        Signature
    )
  ) +

  ggplot2::geom_boxplot(

    position =
      BOX_POSITION,

    width =
      0.62,

    linewidth =
      0.28,

    colour =
      "black",

    alpha =
      0.90,

    outlier.shape =
      NA,

    staplewidth =
      0.45
  ) +


  ###############################################################################
# M1 blue / M2 red
###############################################################################

ggplot2::scale_fill_manual(

  values =
    FIG2H_COLORS,

  breaks = c(
    "M1-like",
    "M2-like"
  ),

  drop =
    FALSE
) +


  ###############################################################################
# 保留全部30个cluster
###############################################################################

ggplot2::scale_x_discrete(
  drop =
    FALSE,

  expand =
    ggplot2::expansion(
      mult = c(
        0.01,
        0.015
      )
    )
) +


  ###############################################################################
# 标题
###############################################################################

ggplot2::labs(

  title =
    "Microglial polarization scores",

  subtitle =
    "M1-like and M2-like transcriptional programs across subclusters",

  y =
    "AddModuleScore"
) +


  theme_fig2h


###############################################################################
# 64. 保存Fig2H
#
# 30个cluster，因此使用横向宽图。
###############################################################################

save_formal_plot(

  plot_object =
    p_fig2h,

  basename =
    "Fig2h_microglia_M1like_M2like_score_boxplot",

  width =
    11.5,

  height =
    4.8
)


###############################################################################
# 65. 保存Fig2H cluster-level agent summary
#
# 这是最适合Agent快速解读每个cluster M1/M2倾向的文件。
###############################################################################

fig2h_cluster_agent_summary <- fig2h_boxplot_summary[
  ,
  .(
    Cluster =
      as.character(
        Cluster
      ),

    Signature =
      as.character(
        Signature
      ),

    n_cells,
    mean_score,
    median_score,
    q25,
    q75,
    cluster_order
  )
]


###############################################################################
# 转wide方便直接比较M1 vs M2
###############################################################################

fig2h_cluster_agent_wide <- data.table::dcast(

  fig2h_cluster_agent_summary,

  Cluster +
    cluster_order ~

    Signature,

  value.var = c(
    "n_cells",
    "mean_score",
    "median_score",
    "q25",
    "q75"
  )
)


###############################################################################
# 添加M1-M2差值
#
# 注意：
# AddModuleScore来自不同gene set，
# 差值只能作为描述性辅助指标，
# 不用于把细胞硬分类成M1/M2。
###############################################################################

fig2h_cluster_agent_wide[
  ,
  mean_M1_minus_M2 :=
    `mean_score_M1-like` -
    `mean_score_M2-like`
]


fig2h_cluster_agent_wide[
  ,
  median_M1_minus_M2 :=
    `median_score_M1-like` -
    `median_score_M2-like`
]


###############################################################################
# 描述性倾向
###############################################################################

fig2h_cluster_agent_wide[
  ,
  descriptive_tendency :=
    data.table::fcase(

      mean_M1_minus_M2 >
        0,
      "M1-like higher",

      mean_M1_minus_M2 <
        0,
      "M2-like higher",

      default =
        "Similar"
    )
]


data.table::setorder(
  fig2h_cluster_agent_wide,
  cluster_order
)


data.table::fwrite(
  fig2h_cluster_agent_wide,
  file.path(
    DATA_DIR,
    "32_Fig2h_cluster_M1M2_agent_summary.csv"
  )
)


###############################################################################
# 66. Fig2H方法说明
###############################################################################

FIG2H_METHOD_NOTE <- c(

  "Fig2H M1-like / M2-like score boxplot",

  "",

  "PLOT DESIGN",

  "The panel reproduces the visual logic of the reference Fig.2h.",

  "For each microglial subcluster, M1-like and M2-like AddModuleScore distributions are displayed side by side.",

  "M1-like scores are shown in blue and M2-like scores in red.",

  "",

  "DATA",

  "The plot directly uses the cell-level AddModuleScore values calculated in the preceding Fig2F/Fig2G workflow.",

  "No Seurat object was reloaded and AddModuleScore was not recalculated.",

  paste0(
    "Number of cells: ",
    nrow(
      score_data
    )
  ),

  "",

  "CLUSTER ORDER",

  paste0(
    "Cluster order is inherited from: ",
    FIG2H_CLUSTER_ORDER_FILE
  ),

  "This keeps Fig2C, Fig2D, Fig2E and Fig2H visually aligned.",

  "",

  "INTERPRETATION",

  "The cell-level boxplot is descriptive and reproduces the presentation logic of the reference figure.",

  "Cells are not treated as independent biological replicates for formal statistical inference.",

  "For downstream formal interpretation, SampleID_Global-level cluster summaries were additionally generated.",

  "M1-like and M2-like scores represent overlapping transcriptional programs rather than binary microglial identities.",

  "",

  "FORMAL SAMPLE-LEVEL DATA",

  paste0(
    "All sample-cluster summaries: ",
    file.path(
      DATA_DIR,
      "30_Fig2h_sample_cluster_M1M2_scores.csv"
    )
  ),

  paste0(
    "Eligible samples (>=20 total cells): ",
    file.path(
      DATA_DIR,
      "31_Fig2h_sample_cluster_M1M2_scores_eligible_samples.csv"
    )
  )
)


writeLines(
  FIG2H_METHOD_NOTE,
  file.path(
    TABLE_DIR,
    "Fig2h_method_note.txt"
  )
)


###############################################################################
# 67. 更新sessionInfo
###############################################################################

capture.output(
  sessionInfo(),
  file = file.path(
    TABLE_DIR,
    "Fig2fgh_microglia_sessionInfo.txt"
  )
)


###############################################################################
# 68. 最终输出
###############################################################################

message("")
message("============================================================")
message("Fig2H完成")
message("============================================================")


message(
  "Cluster数：",
  data.table::uniqueN(
    score_data$Cluster
  )
)


message(
  "Cell-level observations：",
  nrow(
    fig2h_cell_long
  )
)


message("")
message("Figure:")

message(
  file.path(
    FIG_DIR,
    "Fig2h_microglia_M1like_M2like_score_boxplot.png"
  )
)


message("")
message("Formal source data:")

message(
  file.path(
    SOURCE_DIR,
    "Fig2h_M1like_M2like_boxplot_source.csv"
  )
)


message("")
message("Agent推荐读取:")

message(
  file.path(
    DATA_DIR,
    "32_Fig2h_cluster_M1M2_agent_summary.csv"
  )
)


message(
  file.path(
    DATA_DIR,
    "31_Fig2h_sample_cluster_M1M2_scores_eligible_samples.csv"
  )
)


message("")
message("============================================================")



###############################################################################
# 保存 M1 / M2 得分 + 最终实际使用的 marker genes
###############################################################################

SAVE_DIR <- "E:/micro_output/fig3/data"

dir.create(
  SAVE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


###############################################################################
# 1. 保存每个细胞的 M1 / M2 score
###############################################################################

m1m2_score_export <- data.table::copy(
  score_data
)

data.table::fwrite(
  m1m2_score_export,
  file.path(
    SAVE_DIR,
    "M1_M2_cell_scores.csv"
  )
)


###############################################################################
# 同时保存 RDS，后续 R 中读取更快
###############################################################################

saveRDS(
  m1m2_score_export,
  file.path(
    SAVE_DIR,
    "M1_M2_cell_scores.rds"
  )
)


###############################################################################
# 2. 保存最终实际参与 AddModuleScore 的 marker genes
###############################################################################

m1m2_marker_export <- data.table::rbindlist(
  list(

    data.table::data.table(
      Signature = "M1-like",
      Gene = M1_GENES_FINAL
    ),

    data.table::data.table(
      Signature = "M2-like",
      Gene = M2_GENES_FINAL
    )

  )
)


data.table::fwrite(
  m1m2_marker_export,
  file.path(
    SAVE_DIR,
    "M1_M2_marker_genes.csv"
  )
)


###############################################################################
# 3. 如果还想同时保留候选marker、缺失marker和最终使用状态
###############################################################################

m1m2_marker_audit_export <- data.table::copy(
  GENE_AUDIT
)


data.table::fwrite(
  m1m2_marker_audit_export,
  file.path(
    SAVE_DIR,
    "M1_M2_marker_gene_audit.csv"
  )
)


###############################################################################
# 4. 控制台确认
###############################################################################

cat("\n保存完成：\n")

cat(
  file.path(
    SAVE_DIR,
    "M1_M2_cell_scores.csv"
  ),
  "\n"
)

cat(
  file.path(
    SAVE_DIR,
    "M1_M2_cell_scores.rds"
  ),
  "\n"
)

cat(
  file.path(
    SAVE_DIR,
    "M1_M2_marker_genes.csv"
  ),
  "\n"
)

cat(
  file.path(
    SAVE_DIR,
    "M1_M2_marker_gene_audit.csv"
  ),
  "\n"
)
