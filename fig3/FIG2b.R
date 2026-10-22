###############################################################################
# ============================================================================#
# 人脑多疾病小胶质细胞单细胞图谱
# Fig2B - Subcluster signature gene expression heatmap
#
# 输入：
# 1. 正式 Seurat 对象
# 2. 已完成的 FindAllMarkers 结果
# 3. Fig2C 已冻结的 cluster 行顺序
#
# 核心流程：
# 1. 不重新运行 FindAllMarkers
# 2. 每个 cluster 选择2个代表性 positive marker
# 3. 使用 RNA/data 计算代表基因在30个 cluster 中的平均表达
# 4. log1p 后按每个基因跨 cluster 做 Z-score
# 5. 画 signature heatmap
#
# 重要：
# - 这是描述性 marker heatmap，不是新的显著性检验
# - 不修改或覆盖 merged_obj.rds
# - 不重新计算 neighbors / clusters / UMAP / NMF
# - cluster 顺序与已经完成的 Fig2C 保持一致
#
# 输出根目录：
# E:/micro_output/fig3/
#
# figures/
# source_data/
# tables/
# data/
# ============================================================================#
###############################################################################


###############################################################################
# 0. 路径
###############################################################################

INPUT_RDS <- "E:/micro_output/v3_dims16/objects/merged_obj.rds"

MARKER_FILE <-
  "E:/micro_output/v3_dims16/findallmarkers/.staging_res0.6_20260831_174739/all_markers_res0.6.csv"

OUTPUT_ROOT <- "E:/micro_output/fig3"

FIG_DIR <- file.path(
  OUTPUT_ROOT,
  "figures"
)

TABLE_DIR <- file.path(
  OUTPUT_ROOT,
  "tables"
)

SOURCE_DIR <- file.path(
  OUTPUT_ROOT,
  "source_data"
)

DATA_DIR <- file.path(
  OUTPUT_ROOT,
  "data"
)


###############################################################################
# 已完成 Fig2C 的 cluster 排序文件
###############################################################################

PANEL_C_TIER_FILE <- file.path(
  OUTPUT_ROOT,
  "tables",
  "cluster_tier_map.csv"
)


###############################################################################
# 权威审计文件
###############################################################################

OBJECT_AUDIT <-
  "E:/micro_output/v3_dims16/objects/merged_obj_audit.csv"

DISEASE_MAPPING_AUDIT <-
  "E:/micro_output/v3_dims16/objects/merged_obj_disease_mapping.csv"

SAMPLE_DISEASE_AUDIT <-
  "E:/micro_output/v3_dims16/objects/merged_obj_sample_disease_authority_audit.csv"


###############################################################################
# 1. 固定分析参数
###############################################################################

CLUSTER_COL <- "RNA_snn_res.0.6"

SEED <- 42
set.seed(SEED)


###############################################################################
# 每个 cluster 展示几个 signature genes
#
# 文献 Fig2B 是紧凑的 signature heatmap，
# 默认每个 cluster 选2个代表 marker。
###############################################################################

MARKERS_PER_CLUSTER <- 3L


###############################################################################
# Marker 严格筛选阈值
###############################################################################

MARKER_PADJ_MAX <- 0.05

MARKER_LOG2FC_STRICT <- 0.50

MARKER_PCT1_STRICT <- 0.25

MARKER_DELTA_PCT_STRICT <- 0.10


###############################################################################
# Marker 展示时排除高度非特异性的线粒体 / 核糖体基因
#
# 注意：
# 这里只影响“Fig2B代表marker选择”，
# 不会从Seurat对象中删除任何基因。
#
# 如后续想全部保留，可改为：
#
# EXCLUDE_GENE_REGEX <- NULL
###############################################################################

EXCLUDE_GENE_REGEX <- "^(MT-|RPL|RPS)"


###############################################################################
# Heatmap Z-score 显示截断
###############################################################################

Z_CLIP <- 2.5


###############################################################################
# 正式 cluster
###############################################################################

CLUSTER_NUMBERS <- as.character(
  0:29
)

CLUSTER_LABELS <- sprintf(
  "C%02d",
  0:29
)


###############################################################################
# 2. 必要 R 包
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


if (length(missing_packages) > 0L) {

  stop(
    paste0(
      "缺少必要 R 包：",
      paste(
        missing_packages,
        collapse = ", "
      ),
      "\n请手动安装后重新运行。",
      "\n本脚本不会自动安装软件包。"
    )
  )
}


suppressPackageStartupMessages({

  library(Seurat)
  library(data.table)
  library(ggplot2)

})


###############################################################################
# 3. 输入文件检查
###############################################################################

required_files <- c(
  INPUT_RDS,
  MARKER_FILE,
  PANEL_C_TIER_FILE,
  OBJECT_AUDIT,
  DISEASE_MAPPING_AUDIT,
  SAMPLE_DISEASE_AUDIT
)


missing_files <- required_files[
  !file.exists(required_files)
]


if (length(missing_files) > 0L) {

  stop(
    paste0(
      "以下必要输入文件不存在：\n",
      paste0(
        "  - ",
        missing_files,
        collapse = "\n"
      )
    )
  )
}


###############################################################################
# 4. 创建输出目录
###############################################################################

for (dir_i in c(
  FIG_DIR,
  TABLE_DIR,
  SOURCE_DIR,
  DATA_DIR
)) {

  dir.create(
    dir_i,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


###############################################################################
# 5. 读取正式 Seurat 对象
###############################################################################

message("")
message("正在读取正式 Seurat 对象：")
message(INPUT_RDS)


obj <- readRDS(
  INPUT_RDS
)


###############################################################################
# 6. Assay 检查
###############################################################################

if (!"RNA" %in% Assays(obj)) {

  stop(
    "Seurat 对象中不存在正式 RNA assay。"
  )
}


DefaultAssay(obj) <- "RNA"


###############################################################################
# 7. RNA layer 检查
###############################################################################

required_layers <- c(
  "counts",
  "data"
)


available_layers <- Layers(
  obj[["RNA"]]
)


missing_layers <- setdiff(
  required_layers,
  available_layers
)


if (length(missing_layers) > 0L) {

  stop(
    paste0(
      "RNA assay 缺少必要 layer：",
      paste(
        missing_layers,
        collapse = ", "
      )
    )
  )
}


###############################################################################
# 8. Reduction 检查
###############################################################################

required_reductions <- c(
  "pca",
  "harmony_dims16",
  "umap_pca16",
  "umap_harmony16"
)


missing_reductions <- setdiff(
  required_reductions,
  Reductions(obj)
)


if (length(missing_reductions) > 0L) {

  stop(
    paste0(
      "Seurat 对象缺少必要 reduction：",
      paste(
        missing_reductions,
        collapse = ", "
      )
    )
  )
}


###############################################################################
# 9. Graph 检查
###############################################################################

required_graphs <- c(
  "RNA_nn",
  "RNA_snn"
)


missing_graphs <- setdiff(
  required_graphs,
  names(obj@graphs)
)


if (length(missing_graphs) > 0L) {

  stop(
    paste0(
      "Seurat 对象缺少必要 graph：",
      paste(
        missing_graphs,
        collapse = ", "
      )
    )
  )
}


###############################################################################
# 10. Metadata 检查
###############################################################################

required_metadata <- c(
  "CellID",
  "DatasetID",
  "SampleID_Global",
  "Disease",
  "disease_analyse",
  "RNA_snn_res.0.6",
  "seurat_clusters"
)


missing_metadata <- setdiff(
  required_metadata,
  colnames(obj[[]])
)


if (length(missing_metadata) > 0L) {

  stop(
    paste0(
      "Seurat 对象缺少必要 metadata 列：",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  )
}


###############################################################################
# 11. 正式对象基础检查
###############################################################################

stopifnot(
  ncol(obj) == 424222,
  nrow(obj) == 18687,
  "RNA_snn_res.0.6" %in% colnames(obj[[]]),
  "seurat_clusters" %in% colnames(obj[[]]),
  "Disease" %in% colnames(obj[[]]),
  "SampleID_Global" %in% colnames(obj[[]])
)


stopifnot(
  identical(
    as.character(
      obj$seurat_clusters
    ),
    as.character(
      obj$RNA_snn_res.0.6
    )
  )
)


Idents(obj) <- "RNA_snn_res.0.6"


###############################################################################
# 12. Metadata 转换
###############################################################################

md <- data.table::as.data.table(
  obj[[]]
)


###############################################################################
# 13. CellID 检查
###############################################################################

if (anyDuplicated(
  md$CellID
) > 0L) {

  stop(
    "metadata 中 CellID 不是全局唯一。"
  )
}


if (!identical(
  as.character(md$CellID),
  colnames(obj)
)) {

  stop(
    "metadata$CellID 与 colnames(obj) 未严格逐一对应。"
  )
}


###############################################################################
# 14. cluster 检查
###############################################################################

observed_clusters <- unique(
  as.character(
    md[[CLUSTER_COL]]
  )
)


if (!setequal(
  observed_clusters,
  CLUSTER_NUMBERS
)) {

  stop(
    paste0(
      "正式 cluster 不符合 0-29。",
      "\n实际：",
      paste(
        sort(observed_clusters),
        collapse = ", "
      )
    )
  )
}


###############################################################################
# 建立 C00-C29 显示标签
###############################################################################

md[
  ,
  ClusterNum :=
    as.integer(
      as.character(
        get(CLUSTER_COL)
      )
    )
]


if (anyNA(
  md$ClusterNum
)) {

  stop(
    "无法把正式聚类转换为整数 0-29。"
  )
}


md[
  ,
  Cluster :=
    sprintf(
      "C%02d",
      ClusterNum
    )
]


###############################################################################
# 15. 读取 Fig2C cluster 正式排序
###############################################################################

panelc_order <- data.table::fread(
  PANEL_C_TIER_FILE
)


required_panelc_columns <- c(
  "cluster",
  "tier",
  "n_dis_present",
  "total_cells"
)


missing_panelc_columns <- setdiff(
  required_panelc_columns,
  colnames(panelc_order)
)


if (length(missing_panelc_columns) > 0L) {

  stop(
    paste0(
      "cluster_tier_map.csv 缺少必要列：",
      paste(
        missing_panelc_columns,
        collapse = ", "
      )
    )
  )
}


###############################################################################
# cluster 名称标准化函数
###############################################################################

normalize_cluster_label <- function(x) {

  x <- as.character(x)

  output <- rep(
    NA_character_,
    length(x)
  )


  c_style <- grepl(
    "^C[0-9]+$",
    x
  )


  numeric_style <- grepl(
    "^[0-9]+$",
    x
  )


  output[c_style] <- sprintf(
    "C%02d",
    as.integer(
      sub(
        "^C",
        "",
        x[c_style]
      )
    )
  )


  output[numeric_style] <- sprintf(
    "C%02d",
    as.integer(
      x[numeric_style]
    )
  )


  if (anyNA(
    output
  )) {

    stop(
      "Fig2C cluster_tier_map.csv 中存在无法识别的 cluster 标签。"
    )
  }


  return(
    output
  )
}


panelc_order[
  ,
  Cluster :=
    normalize_cluster_label(
      cluster
    )
]


###############################################################################
# Fig2C 应为30个cluster
###############################################################################

if (nrow(
  panelc_order
) != 30L) {

  stop(
    paste0(
      "Fig2C cluster_tier_map.csv 应包含30个cluster，实际为 ",
      nrow(panelc_order)
    )
  )
}


if (anyDuplicated(
  panelc_order$Cluster
) > 0L) {

  stop(
    "Fig2C cluster_tier_map.csv 存在重复 cluster。"
  )
}


if (!setequal(
  panelc_order$Cluster,
  CLUSTER_LABELS
)) {

  stop(
    "Fig2C cluster 集合与当前正式 C00-C29 不一致。"
  )
}


###############################################################################
# Fig2B cluster 列顺序：
# 与 Fig2C 自上而下 cluster 顺序一致
###############################################################################

FIG2B_CLUSTER_ORDER <- panelc_order$Cluster


###############################################################################
# 保存 cluster order
###############################################################################

fig2b_cluster_order <- panelc_order[
  ,
  .(
    Cluster,
    tier =
      as.character(
        tier
      ),
    n_dis_present,
    total_cells
  )
]


fig2b_cluster_order[
  ,
  column_order :=
    seq_len(.N)
]


data.table::fwrite(
  fig2b_cluster_order,
  file.path(
    DATA_DIR,
    "10_Fig2b_cluster_column_order.csv"
  )
)


###############################################################################
# 16. 读取 FindAllMarkers 结果
###############################################################################

message("")
message("读取 FindAllMarkers：")
message(MARKER_FILE)


markers <- data.table::fread(
  MARKER_FILE
)


###############################################################################
# 17. Marker 文件列检查
###############################################################################

required_marker_columns <- c(
  "p_val",
  "avg_log2FC",
  "pct.1",
  "pct.2",
  "p_val_adj",
  "cluster",
  "gene"
)


missing_marker_columns <- setdiff(
  required_marker_columns,
  colnames(markers)
)


if (length(missing_marker_columns) > 0L) {

  stop(
    paste0(
      "FindAllMarkers 文件缺少必要列：",
      paste(
        missing_marker_columns,
        collapse = ", "
      )
    )
  )
}


###############################################################################
# 18. Marker 基础清理
###############################################################################

markers[
  ,
  gene :=
    trimws(
      as.character(
        gene
      )
    )
]


markers[
  ,
  cluster :=
    as.character(
      cluster
    )
]


###############################################################################
# 空 gene
###############################################################################

if (any(
  is.na(
    markers$gene
  ) |
  markers$gene == ""
)) {

  stop(
    "FindAllMarkers 文件中存在 NA 或空 gene。"
  )
}


###############################################################################
# cluster 必须0-29
###############################################################################

if (!setequal(
  unique(
    markers$cluster
  ),
  CLUSTER_NUMBERS
)) {

  stop(
    paste0(
      "FindAllMarkers 的 cluster 与正式0-29不一致。",
      "\nMarker文件实际 cluster：",
      paste(
        sort(
          unique(
            markers$cluster
          )
        ),
        collapse = ", "
      )
    )
  )
}


###############################################################################
# 数值列检查
###############################################################################

numeric_marker_columns <- c(
  "p_val",
  "avg_log2FC",
  "pct.1",
  "pct.2",
  "p_val_adj"
)


for (col_i in numeric_marker_columns) {

  if (!is.numeric(
    markers[[col_i]]
  )) {

    stop(
      paste0(
        "Marker列 ",
        col_i,
        " 不是 numeric。"
      )
    )
  }
}


###############################################################################
# 建立 C00-C29
###############################################################################

markers[
  ,
  MarkerCluster :=
    sprintf(
      "C%02d",
      as.integer(
        cluster
      )
    )
]


###############################################################################
# pct差值
###############################################################################

markers[
  ,
  delta_pct :=
    pct.1 -
    pct.2
]


###############################################################################
# 19. 检查 marker gene 是否属于当前 RNA assay
###############################################################################

rna_genes <- rownames(
  obj[["RNA"]]
)


markers[
  ,
  gene_in_RNA :=
    gene %chin%
    rna_genes
]


missing_marker_gene_table <- markers[
  gene_in_RNA == FALSE,
  .(
    n_marker_rows = .N
  ),
  by = gene
]


if (nrow(
  missing_marker_gene_table
) > 0L) {

  data.table::fwrite(
    missing_marker_gene_table,
    file.path(
      TABLE_DIR,
      "Fig2b_marker_genes_missing_from_RNA.csv"
    )
  )


  warning(
    paste0(
      nrow(
        missing_marker_gene_table
      ),
      " 个 marker gene 不存在于当前 RNA assay，",
      "这些基因不会参与 Fig2B marker 选择。"
    )
  )
}


markers <- markers[
  gene_in_RNA == TRUE
]


###############################################################################
# 20. 去掉 Cluster 内重复 gene
#
# 若同一 cluster/gene 出现多次，
# 保留 avg_log2FC 更大的记录。
###############################################################################

data.table::setorder(
  markers,
  MarkerCluster,
  gene,
  -avg_log2FC,
  p_val_adj
)


markers <- unique(
  markers,
  by = c(
    "MarkerCluster",
    "gene"
  )
)


###############################################################################
# 21. 仅用于 signature gene 展示的基因过滤
###############################################################################

markers[
  ,
  excluded_by_pattern :=
    FALSE
]


if (!is.null(
  EXCLUDE_GENE_REGEX
)) {

  markers[
    ,
    excluded_by_pattern :=
      grepl(
        EXCLUDE_GENE_REGEX,
        gene
      )
  ]
}


###############################################################################
# 22. 定义 marker selection level
#
# strict:
#   padj <= 0.05
#   log2FC >= 0.5
#   pct1 >= 0.25
#   pct1-pct2 >= 0.10
#
# relaxed:
#   padj <= 0.05
#   log2FC > 0
#   pct1 > pct2
#
# fallback:
#   log2FC > 0
#
# exclude:
#   其他情况
###############################################################################

markers[
  ,
  selection_level :=
    data.table::fcase(

      !excluded_by_pattern &
        !is.na(p_val_adj) &
        p_val_adj <=
        MARKER_PADJ_MAX &
        avg_log2FC >=
        MARKER_LOG2FC_STRICT &
        pct.1 >=
        MARKER_PCT1_STRICT &
        delta_pct >=
        MARKER_DELTA_PCT_STRICT,

      "strict",


      !excluded_by_pattern &
        !is.na(p_val_adj) &
        p_val_adj <=
        MARKER_PADJ_MAX &
        avg_log2FC >
        0 &
        pct.1 >
        pct.2,

      "relaxed",


      !excluded_by_pattern &
        avg_log2FC >
        0,

      "fallback",


      default =
        "exclude"
    )
]


markers[
  ,
  selection_priority :=
    match(
      selection_level,
      c(
        "strict",
        "relaxed",
        "fallback",
        "exclude"
      )
    )
]


###############################################################################
# 23. Marker candidate audit
###############################################################################

marker_candidate_audit <- markers[
  ,
  .(
    n_all =
      .N,

    n_strict =
      sum(
        selection_level ==
          "strict"
      ),

    n_relaxed =
      sum(
        selection_level ==
          "relaxed"
      ),

    n_fallback =
      sum(
        selection_level ==
          "fallback"
      ),

    n_excluded =
      sum(
        selection_level ==
          "exclude"
      )
  ),
  by = MarkerCluster
]


marker_candidate_audit[
  ,
  cluster_order :=
    match(
      MarkerCluster,
      FIG2B_CLUSTER_ORDER
    )
]


data.table::setorder(
  marker_candidate_audit,
  cluster_order
)


marker_candidate_audit[
  ,
  cluster_order := NULL
]


data.table::fwrite(
  marker_candidate_audit,
  file.path(
    TABLE_DIR,
    "Fig2b_marker_candidate_audit.csv"
  )
)


###############################################################################
# 24. 每个 cluster 选择代表 marker
#
# 注意：
# 不禁止同一个 gene 成为不同 cluster 的 marker。
#
# 因为如果一个基因确实是多个亚群的高表达 marker，
# 强行全局去重会人为改变 FindAllMarkers 结果。
#
# 若同一 gene 被多个 cluster 选中，
# heatmap 行标签会显示：
#
# GENE [Cxx]
#
# 从而保持透明。
###############################################################################

selected_marker_list <- vector(
  mode = "list",
  length = length(
    FIG2B_CLUSTER_ORDER
  )
)


for (i in seq_along(
  FIG2B_CLUSTER_ORDER
)) {

  cluster_i <- FIG2B_CLUSTER_ORDER[i]


  candidate_i <- markers[
    MarkerCluster ==
      cluster_i &
      selection_level !=
      "exclude"
  ]


  ###########################################################################
  # 排序优先级：
  # strict > relaxed > fallback
  # 然后：
  # avg_log2FC 高
  # delta_pct 高
  # pct.1 高
  # p_val_adj 小
  ###########################################################################

  data.table::setorder(
    candidate_i,
    selection_priority,
    -avg_log2FC,
    -delta_pct,
    -pct.1,
    p_val_adj,
    gene
  )


  if (nrow(
    candidate_i
  ) < MARKERS_PER_CLUSTER) {

    stop(
      paste0(
        "Cluster ",
        cluster_i,
        " 在当前 marker 筛选规则下不足 ",
        MARKERS_PER_CLUSTER,
        " 个可用 positive markers。",
        "\n实际只有 ",
        nrow(candidate_i),
        " 个。",
        "\n请检查 FindAllMarkers 结果或放宽顶部参数。"
      )
    )
  }


  selected_i <- candidate_i[
    seq_len(
      MARKERS_PER_CLUSTER
    )
  ]


  selected_i[
    ,
    marker_rank :=
      seq_len(.N)
  ]


  selected_marker_list[[i]] <- selected_i
}


###############################################################################
# 25. 合并 selected markers
###############################################################################

selected_markers <- data.table::rbindlist(
  selected_marker_list,
  use.names = TRUE,
  fill = TRUE
)


###############################################################################
# 应为 30 × MARKERS_PER_CLUSTER
###############################################################################

expected_marker_count <-
  30L *
  MARKERS_PER_CLUSTER


if (nrow(
  selected_markers
) != expected_marker_count) {

  stop(
    paste0(
      "Selected marker 数异常。",
      "\n预期：",
      expected_marker_count,
      "\n实际：",
      nrow(
        selected_markers
      )
    )
  )
}


###############################################################################
# 26. 标记重复 selected genes
###############################################################################

selected_markers[
  ,
  selected_occurrence_n :=
    .N,
  by = gene
]


selected_markers[
  ,
  GeneRow :=
    ifelse(
      selected_occurrence_n >
        1L,

      paste0(
        gene,
        " [",
        MarkerCluster,
        "]"
      ),

      gene
    )
]


###############################################################################
# GeneRow 必须唯一
###############################################################################

if (anyDuplicated(
  selected_markers$GeneRow
) > 0L) {

  stop(
    "生成的 GeneRow 标签不唯一。"
  )
}


###############################################################################
# 27. Selected marker 排序
###############################################################################

selected_markers[
  ,
  cluster_column_order :=
    match(
      MarkerCluster,
      FIG2B_CLUSTER_ORDER
    )
]


data.table::setorder(
  selected_markers,
  cluster_column_order,
  marker_rank
)


###############################################################################
# row_order = heatmap 自上而下
###############################################################################

selected_markers[
  ,
  gene_row_order :=
    seq_len(.N)
]


###############################################################################
# 28. 保存 selected markers
###############################################################################

selected_marker_export <- selected_markers[
  ,
  .(
    gene,
    GeneRow,
    MarkerCluster,
    marker_rank,
    gene_row_order,
    cluster_column_order,
    p_val,
    avg_log2FC,
    pct.1,
    pct.2,
    delta_pct,
    p_val_adj,
    selection_level,
    selected_occurrence_n
  )
]


data.table::fwrite(
  selected_marker_export,
  file.path(
    DATA_DIR,
    "11_Fig2b_selected_markers.csv"
  )
)


data.table::fwrite(
  selected_marker_export,
  file.path(
    SOURCE_DIR,
    "Fig2b_selected_signature_genes.csv"
  )
)


###############################################################################
# 29. selected marker 汇总 audit
###############################################################################

selected_marker_audit <- selected_markers[
  ,
  .(
    n_selected =
      .N,

    n_strict =
      sum(
        selection_level ==
          "strict"
      ),

    n_relaxed =
      sum(
        selection_level ==
          "relaxed"
      ),

    n_fallback =
      sum(
        selection_level ==
          "fallback"
      ),

    marker_genes =
      paste(
        gene,
        collapse = ";"
      )
  ),
  by = MarkerCluster
]


selected_marker_audit[
  ,
  cluster_order :=
    match(
      MarkerCluster,
      FIG2B_CLUSTER_ORDER
    )
]


data.table::setorder(
  selected_marker_audit,
  cluster_order
)


selected_marker_audit[
  ,
  cluster_order := NULL
]


data.table::fwrite(
  selected_marker_audit,
  file.path(
    TABLE_DIR,
    "Fig2b_selected_marker_audit.csv"
  )
)


###############################################################################
# 30. Heatmap 所需要的 unique genes
###############################################################################

selected_unique_genes <- unique(
  selected_markers$gene
)


if (length(
  selected_unique_genes
) == 0L) {

  stop(
    "没有 selected signature genes。"
  )
}


###############################################################################
# 所有selected genes必须存在于RNA assay
###############################################################################

missing_selected_genes <- setdiff(
  selected_unique_genes,
  rownames(
    obj[["RNA"]]
  )
)


if (length(
  missing_selected_genes
) > 0L) {

  stop(
    paste0(
      "以下 selected genes 不存在于 RNA assay：",
      paste(
        missing_selected_genes,
        collapse = ", "
      )
    )
  )
}


###############################################################################
# 31. 给 Seurat 对象增加临时绘图分组
#
# 仅修改当前R内存中的 obj，
# 不保存回 merged_obj.rds。
###############################################################################

obj$Fig2B_Cluster <- factor(
  md$Cluster,
  levels =
    FIG2B_CLUSTER_ORDER
)


###############################################################################
# 32. 每个 Cluster 的细胞数 audit
###############################################################################

cluster_cell_audit <- md[
  ,
  .(
    n_cells = .N
  ),
  by = Cluster
]


cluster_cell_audit[
  ,
  column_order :=
    match(
      Cluster,
      FIG2B_CLUSTER_ORDER
    )
]


data.table::setorder(
  cluster_cell_audit,
  column_order
)


data.table::fwrite(
  cluster_cell_audit,
  file.path(
    TABLE_DIR,
    "Fig2b_cluster_cell_count_audit.csv"
  )
)


###############################################################################
# 33. 使用 Seurat 官方 AverageExpression
#
# RNA/data：
# 当前正式 log-normalized RNA expression。
#
# AverageExpression(layer="data")：
# Seurat 会根据 data layer 计算每个 cluster 的平均表达。
#
# 这里只是描述性 heatmap，
# 不作为差异表达显著性统计。
###############################################################################

message("")
message(
  "正在计算 ",
  length(
    selected_unique_genes
  ),
  " 个 unique signature genes × 30 clusters 的平均表达..."
)


avg_expression_list <- Seurat::AverageExpression(
  object = obj,
  assays = "RNA",
  features = selected_unique_genes,
  group.by = "Fig2B_Cluster",
  layer = "data",
  return.seurat = FALSE,
  verbose = FALSE
)


###############################################################################
# 34. 提取 RNA average expression matrix
###############################################################################

if (!"RNA" %in% names(
  avg_expression_list
)) {

  stop(
    "AverageExpression 结果中不存在 RNA matrix。"
  )
}


avg_expression_matrix <- avg_expression_list$RNA


###############################################################################
# 矩阵基础检查
###############################################################################

if (!all(
  selected_unique_genes %in%
  rownames(
    avg_expression_matrix
  )
)) {

  stop(
    "AverageExpression 结果缺少部分 selected genes。"
  )
}


if (!setequal(
  colnames(
    avg_expression_matrix
  ),
  FIG2B_CLUSTER_ORDER
)) {

  stop(
    paste0(
      "AverageExpression cluster 列名与 Fig2C 顺序文件不一致。",
      "\n实际列：",
      paste(
        colnames(
          avg_expression_matrix
        ),
        collapse = ", "
      )
    )
  )
}


###############################################################################
# 35. 正式按照 Fig2C cluster 顺序重排
###############################################################################

avg_expression_matrix <- avg_expression_matrix[
  ,
  FIG2B_CLUSTER_ORDER,
  drop = FALSE
]


###############################################################################
# unique gene 顺序按 marker heatmap 首次出现顺序
###############################################################################

unique_gene_order <- unique(
  selected_markers$gene
)


avg_expression_matrix <- avg_expression_matrix[
  unique_gene_order,
  ,
  drop = FALSE
]


###############################################################################
# 36. 保存 average expression wide table
###############################################################################

avg_expression_wide <- data.table::data.table(
  gene =
    rownames(
      avg_expression_matrix
    )
)


avg_expression_wide <- cbind(
  avg_expression_wide,
  data.table::as.data.table(
    avg_expression_matrix
  )
)


data.table::fwrite(
  avg_expression_wide,
  file.path(
    DATA_DIR,
    "12_Fig2b_average_expression_wide.csv"
  )
)


###############################################################################
# 37. 转成长表
###############################################################################

avg_expression_long <- data.table::melt(
  avg_expression_wide,
  id.vars = "gene",
  variable.name = "ExpressionCluster",
  value.name = "avg_expression"
)


avg_expression_long[
  ,
  ExpressionCluster :=
    as.character(
      ExpressionCluster
    )
]


###############################################################################
# 38. log1p average expression
#
# AverageExpression 是平均表达，
# 为避免少数极高表达基因支配视觉范围，
# Z-score 前使用 log1p。
###############################################################################

avg_expression_long[
  ,
  avg_expression_log1p :=
    log1p(
      avg_expression
    )
]


###############################################################################
# 39. 将 selected marker 信息复制到每个 expression cluster
#
# 如果同一个gene被两个 MarkerCluster 选中，
# 它会得到两个 GeneRow，
# 但实际表达值完全来自同一真实gene。
###############################################################################

marker_for_join <- selected_markers[
  ,
  .(
    gene,
    GeneRow,
    MarkerCluster,
    marker_rank,
    gene_row_order,
    cluster_column_order,
    marker_avg_log2FC =
      avg_log2FC,
    marker_pct1 =
      pct.1,
    marker_pct2 =
      pct.2,
    marker_delta_pct =
      delta_pct,
    marker_p_val_adj =
      p_val_adj,
    selection_level
  )
]


fig2b_plot_data <- merge(
  marker_for_join,
  avg_expression_long,
  by = "gene",
  allow.cartesian = TRUE,
  sort = FALSE
)


###############################################################################
# 应为：
#
# selected marker rows × 30 expression clusters
###############################################################################

expected_heatmap_rows <-
  nrow(
    selected_markers
  ) *
  30L


if (nrow(
  fig2b_plot_data
) != expected_heatmap_rows) {

  stop(
    paste0(
      "Fig2B heatmap long table 行数异常。",
      "\n预期：",
      expected_heatmap_rows,
      "\n实际：",
      nrow(
        fig2b_plot_data
      )
    )
  )
}


###############################################################################
# 40. 每个 GeneRow 跨30个 cluster 做 Z-score
###############################################################################

fig2b_plot_data[
  ,
  expression_z :=
    {

      expression_vector <-
        avg_expression_log1p


      expression_mean <-
        mean(
          expression_vector,
          na.rm = TRUE
        )


      expression_sd <-
        stats::sd(
          expression_vector,
          na.rm = TRUE
        )


      if (
        !is.finite(
          expression_sd
        ) ||
        expression_sd ==
        0
      ) {

        rep(
          0,
          .N
        )

      } else {

        (
          expression_vector -
            expression_mean
        ) /
          expression_sd
      }

    },
  by = GeneRow
]


###############################################################################
# 41. Z-score 截断
###############################################################################

fig2b_plot_data[
  ,
  expression_z_plot :=
    pmax(
      -Z_CLIP,
      pmin(
        Z_CLIP,
        expression_z
      )
    )
]


###############################################################################
# 42. 表达 cluster 列顺序
###############################################################################

fig2b_plot_data[
  ,
  expression_cluster_order :=
    match(
      ExpressionCluster,
      FIG2B_CLUSTER_ORDER
    )
]


###############################################################################
# 排序：gene row -> cluster column
###############################################################################

data.table::setorder(
  fig2b_plot_data,
  gene_row_order,
  expression_cluster_order
)


###############################################################################
# 43. 保存最终绘图中间数据
#
# 这是以后 Agent 解读 Fig2B 的首选文件。
###############################################################################

data.table::fwrite(
  fig2b_plot_data,
  file.path(
    DATA_DIR,
    "13_Fig2b_plot_data_long.csv"
  )
)


###############################################################################
# 同时作为正式 source data
###############################################################################

data.table::fwrite(
  fig2b_plot_data,
  file.path(
    SOURCE_DIR,
    "Fig2b_signature_expression_heatmap_source.csv"
  )
)


###############################################################################
# 44. Heatmap gene row order
###############################################################################

gene_row_order_table <- selected_markers[
  ,
  .(
    GeneRow,
    gene,
    MarkerCluster,
    marker_rank,
    gene_row_order,
    cluster_column_order,
    avg_log2FC,
    pct.1,
    pct.2,
    delta_pct,
    p_val_adj,
    selection_level
  )
]


data.table::setorder(
  gene_row_order_table,
  gene_row_order
)


data.table::fwrite(
  gene_row_order_table,
  file.path(
    DATA_DIR,
    "14_Fig2b_gene_row_order.csv"
  )
)


###############################################################################
# 45. 准备 ggplot factor
#
# y轴：
# gene_row_order=1 要显示在最上方。
#
# ggplot 第一个 factor level 位于底部，
# 因此需要 reverse。
###############################################################################

gene_order_top_to_bottom <- gene_row_order_table$GeneRow


fig2b_plot_data[
  ,
  GeneRow :=
    factor(
      GeneRow,
      levels =
        rev(
          gene_order_top_to_bottom
        )
    )
]


fig2b_plot_data[
  ,
  ExpressionCluster :=
    factor(
      ExpressionCluster,
      levels =
        FIG2B_CLUSTER_ORDER
    )
]


###############################################################################
# 46. 计算不同 MarkerCluster 之间的横向分割线
###############################################################################

gene_group_sizes <- selected_markers[
  ,
  .(
    n_gene =
      .N
  ),
  by = .(
    MarkerCluster,
    cluster_column_order
  )
]


data.table::setorder(
  gene_group_sizes,
  cluster_column_order
)


gene_group_sizes[
  ,
  cumulative_gene :=
    cumsum(
      n_gene
    )
]


###############################################################################
# 最后一组之后不需要分割线
###############################################################################

if (nrow(
  gene_group_sizes
) > 1L) {

  cumulative_before_last <-
    gene_group_sizes[
      1:(.N - 1L),
      cumulative_gene
    ]


  horizontal_boundaries <-
    nrow(
      selected_markers
    ) -
    cumulative_before_last +
    0.5

} else {

  horizontal_boundaries <-
    numeric(0)
}


###############################################################################
###############################################################################
# 47. Fig2B SCI发表级热图主题
#
# 配色：
#   Blue  = 相对低表达
#   White = 接近平均表达
#   Red   = 相对高表达
#
# 每个小格：
#   强制正方形
#   带细白边
###############################################################################

FONT <- "Arial"


###############################################################################
# SCI级别 Blue-White-Red 配色
#
# 使用较低饱和度的经典 diverging colors，
# 避免纯 #FF0000 / #0000FF 过于刺眼。
###############################################################################

HEATMAP_LOW <- "#2166AC"
HEATMAP_MID <- "#F7F7F7"
HEATMAP_HIGH <- "#B2182B"


###############################################################################
# Heatmap主题
###############################################################################

theme_fig2b <- ggplot2::theme_bw(
  base_family = FONT,
  base_size = 8
) +

  ggplot2::theme(

    ###########################################################################
    # 不使用默认panel grid
    ###########################################################################

    panel.grid =
      ggplot2::element_blank(),


    ###########################################################################
    # 去掉外层黑框
    ###########################################################################

    panel.border =
      ggplot2::element_blank(),


    ###########################################################################
    # 不显示坐标轴线
    ###########################################################################

    axis.line =
      ggplot2::element_blank(),


    ###########################################################################
    # 不显示ticks
    ###########################################################################

    axis.ticks =
      ggplot2::element_blank(),


    ###########################################################################
    # 不需要axis title
    ###########################################################################

    axis.title =
      ggplot2::element_blank(),


    ###########################################################################
    # Cluster标签
    ###########################################################################

    axis.text.x =
      ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 6.8,
        colour = "black"
      ),


    ###########################################################################
    # Gene标签
    ###########################################################################

    axis.text.y =
      ggplot2::element_text(
        size = 5.8,
        colour = "black"
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
        0.55,
        "cm"
      ),


    ###########################################################################
    # 白色背景
    ###########################################################################

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


    ###########################################################################
    # 边距
    ###########################################################################

    plot.margin =
      ggplot2::margin(
        t = 6,
        r = 10,
        b = 4,
        l = 4
      )
  )


###############################################################################
# 48. Fig2B heatmap
#
# 关键变化：
#
# 1. geom_tile：
#    color="white"
#    给每个小格增加独立白边
#
# 2. coord_fixed(ratio=1)：
#    强制每个heatmap tile为正方形
#
# 3. scale_fill_gradient2：
#    对称Blue-White-Red
###############################################################################

p_fig2b <-

  ggplot2::ggplot(
    fig2b_plot_data,
    ggplot2::aes(
      x = ExpressionCluster,
      y = GeneRow,
      fill = expression_z_plot
    )
  ) +


  ###########################################################################
# 每个格子
#
# linewidth可以根据最终PDF效果微调：
# 0.15 更细
# 0.20 推荐
# 0.25 更明显
###########################################################################

ggplot2::geom_tile(
  colour = "white",
  linewidth = 0.20,
  width = 1,
  height = 1
) +


  ###########################################################################
# 不同marker来源cluster之间增加稍粗一点的白色横向分隔
#
# 这样可以看到每2个marker属于一个subcluster，
# 但不会破坏整体heatmap。
###########################################################################

ggplot2::geom_hline(
  yintercept = horizontal_boundaries,
  colour = "white",
  linewidth = 0.65
) +


  ###########################################################################
# SCI Blue-White-Red
#
# Z-score < 0 → blue
# Z-score = 0 → white
# Z-score > 0 → red
#
# limits严格对称
###########################################################################

ggplot2::scale_fill_gradient2(

  name = "Row Z-score",

  low = HEATMAP_LOW,

  mid = HEATMAP_MID,

  high = HEATMAP_HIGH,

  midpoint = 0,

  limits = c(
    -Z_CLIP,
    Z_CLIP
  ),

  breaks = c(
    -2,
    -1,
    0,
    1,
    2
  ),

  labels = c(
    "-2",
    "-1",
    "0",
    "1",
    "2"
  ),

  oob = scales::squish
) +


  ###########################################################################
# 不留tile间额外空间
###########################################################################

ggplot2::scale_x_discrete(
  expand = c(
    0,
    0
  )
) +


  ggplot2::scale_y_discrete(
    expand = c(
      0,
      0
    )
  ) +


  ###########################################################################
# 强制正方形格子
#
# x和y各一个离散单位拥有相同物理长度
###########################################################################

ggplot2::coord_fixed(
  ratio = 1,
  clip = "off"
) +


  ###########################################################################
# 标题
###########################################################################

ggplot2::labs(

  title =
    "Subcluster signature gene expression",

  subtitle =
    paste0(
      "Top ",
      MARKERS_PER_CLUSTER,
      " positive markers per cluster; ",
      "row-scaled average RNA expression"
    )
) +


  theme_fig2b


###############################################################################
# 49. 正式图片保存函数
###############################################################################

save_formal_plot <- function(
    plot_object,
    basename,
    width,
    height
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
  # PNG 600 dpi
  ###########################################################################

  png_success <- tryCatch(

    {

      ggplot2::ggsave(
        filename = png_file,
        plot = plot_object,
        width = width,
        height = height,
        units = "in",
        dpi = 600,
        bg = "white"
      )

      TRUE
    },

    error = function(e) {

      FALSE
    }
  )


  if (!png_success) {

    warning(
      paste0(
        "Arial PNG 输出失败，使用 sans 字体回退：",
        basename
      )
    )


    fallback_plot <-
      plot_object +

      ggplot2::theme(
        text =
          ggplot2::element_text(
            family = "sans"
          )
      )


    ggplot2::ggsave(
      filename = png_file,
      plot = fallback_plot,
      width = width,
      height = height,
      units = "in",
      dpi = 600,
      bg = "white"
    )
  }


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
          filename = pdf_file,
          plot = plot_object,
          width = width,
          height = height,
          units = "in",
          device = grDevices::cairo_pdf,
          family = FONT,
          bg = "white"
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
        "Arial/Cairo PDF 输出失败，使用 sans 字体回退：",
        basename
      )
    )


    fallback_plot <-
      plot_object +

      ggplot2::theme(
        text =
          ggplot2::element_text(
            family = "sans"
          )
      )


    grDevices::pdf(
      file = pdf_file,
      width = width,
      height = height,
      family = "sans",
      useDingbats = FALSE
    )


    print(
      fallback_plot
    )


    grDevices::dev.off()
  }
}


###############################################################################
# 50. 保存新版 Fig2B
#
# 30列 × 60行
#
# 因为每个格子要求正方形，
# heatmap主体的理论高宽比约为60/30 = 2。
#
# 考虑左侧gene名、底部cluster名、右侧legend，
# 推荐最终尺寸：
#
# width  = 8.5 inch
# height = 15.5 inch
###############################################################################

save_formal_plot(

  plot_object =
    p_fig2b,

  basename =
    "Fig2b_subcluster_signature_heatmap_SCI_redblue",

  width =
    8.5,

  height =
    15.5
)


###############################################################################
# 51. 控制台提示
###############################################################################

message("")
message("============================================================")
message("Fig2B SCI red-blue heatmap 已输出")
message("============================================================")

message(
  "PNG: ",
  file.path(
    FIG_DIR,
    "Fig2b_subcluster_signature_heatmap_SCI_redblue.png"
  )
)

message(
  "PDF: ",
  file.path(
    FIG_DIR,
    "Fig2b_subcluster_signature_heatmap_SCI_redblue.pdf"
  )
)

message("")
message(
  "Color scale: ",
  HEATMAP_LOW,
  " -> ",
  HEATMAP_MID,
  " -> ",
  HEATMAP_HIGH
)

message(
  "Z-score range: -",
  Z_CLIP,
  " to +",
  Z_CLIP
)

message(
  "Tile aspect ratio: 1:1"
)

message("============================================================")
