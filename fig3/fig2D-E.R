###############################################################################
# ============================================================================#
# 人脑多疾病小胶质细胞图谱
# Fig2D + Fig2E
#
# Fig2D:
#   Disease17 × Microglia subcluster
#   样本内计算 cluster proportion
#   Disease 内 SampleID_Global 等权平均
#   气泡大小 = mean proportion
#   气泡颜色 = <1%, 1-5%, 5-10%, >10%
#
# Fig2E:
#   Tumor / Neurological / Ctrl × Microglia subcluster
#   样本级 preference OR
#
# 行顺序：
#   严格读取已经完成的 Fig2C：
#   E:/micro_output/fig3/tables/cluster_tier_map.csv
#
# 输出：
#   figures/     正式 PNG + PDF
#   source_data/ 正式 source CSV
#   tables/      audit / QC / sessionInfo
#   data/        分析和绘图中间数据，供 Agent 解读
#
# 统计原则：
#   - SampleID_Global 为统计单位
#   - 每个样本 >=20 cells
#   - Disease 聚合时样本等权
#   - 不重新计算 neighbors / clusters / UMAP / NMF
#   - 不修改 merged_obj.rds
# ============================================================================#
###############################################################################


###############################################################################
# 0. 固定路径
###############################################################################

INPUT_RDS <- "E:/micro_output/v3_dims16/objects/merged_obj.rds"

# 项目根目录
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
# Fig2C 已冻结结果
#
# 该文件决定 Fig2C / Fig2D / Fig2E 的最终 cluster 行顺序
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
# 权威颜色文件
###############################################################################

DISEASE_PALETTE_FILE <-
  "D:/data/microglia_new/authority/03_reference_configs/disease_palette_contract.csv"

DATASET_PALETTE_FILE <-
  "D:/data/microglia_new/authority/03_reference_configs/dataset_palette_contract.csv"


###############################################################################
# 其他参数
###############################################################################

CLUSTER_COL <- "RNA_snn_res.0.6"

GROUP_COL <- "Disease_Group3"

MIN_CELLS_PER_SAMPLE <- 20L

# 与 Fig2C 完全一致
MIN_CELLS_PRESENT <- 20L
PAN_MIN <- 8L
SPEC_MAX <- 2L

SEED <- 42

set.seed(SEED)


###############################################################################
# 1. Disease17 正式顺序
###############################################################################

DISEASE_ORDER_CONTRACT <- c(
  "Ctrl",
  "LGG",
  "GBM",
  "MB",
  "LCBM",
  "AD",
  "EOAD",
  "PD",
  "DLB",
  "DLBD_PD",
  "ALS",
  "ALS_FTD",
  "FTD",
  "FTD_MND",
  "PSP",
  "HD",
  "MS"
)


###############################################################################
# Disease17 正式细胞数
###############################################################################

EXPECTED_DISEASE_COUNTS <- c(
  Ctrl     = 49321L,
  LGG      = 10116L,
  GBM      = 103114L,
  MB       = 3385L,
  LCBM     = 32213L,
  AD       = 83652L,
  EOAD     = 8266L,
  PD       = 14612L,
  DLB      = 85L,
  DLBD_PD  = 14226L,
  ALS      = 56958L,
  ALS_FTD  = 12332L,
  FTD      = 5681L,
  FTD_MND  = 336L,
  PSP      = 6265L,
  HD       = 4529L,
  MS       = 19131L
)


###############################################################################
# 2. 三大疾病类别
#
# Tumor:
#   LGG / GBM / MB / LCBM
#
# Neurological:
#   其余神经系统疾病
#
# Ctrl:
#   Ctrl
###############################################################################

TUMOR_DISEASES <- c(
  "LGG",
  "GBM",
  "MB",
  "LCBM"
)


NEUROLOGICAL_DISEASES <- c(
  "AD",
  "EOAD",
  "PD",
  "DLB",
  "DLBD_PD",
  "ALS",
  "ALS_FTD",
  "FTD",
  "FTD_MND",
  "PSP",
  "HD",
  "MS"
)


GROUP_ORDER <- c(
  "Tumor",
  "Neurological",
  "Ctrl"
)


###############################################################################
# 3. Cluster 定义
###############################################################################

CLUSTER_NUMBERS <- as.character(
  0:29
)

CLUSTER_LABELS <- sprintf(
  "C%02d",
  0:29
)


###############################################################################
# 4. Fig2D 固定配色
#
# 与你已经完成的 Fig2C/PanelD 脚本保持一致
###############################################################################

COL_BLUES <- c(
  "<1%"   = "#D9F0D3",
  "1-5%"  = "#D9EAF2",
  "5-10%" = "#4DB6AC",
  ">10%"  = "#2166AC"
)


###############################################################################
# 与文献 / 已完成脚本保持一致的气泡图例
###############################################################################

PCT_BREAKS <- c(
  1,
  10,
  25,
  50
)


###############################################################################
# 5. Fig2E OR 配色
###############################################################################

OR_COLORS <- c(
  "#FFF7EC",
  "#FEE8C8",
  "#FDBB84",
  "#E34A33"
)

OR_EPS <- 1e-6
OR_PLOT_CAP <- 3


###############################################################################
# 6. 必要 R 包检查
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
# 7. 输入文件检查
###############################################################################

required_files <- c(
  INPUT_RDS,
  OBJECT_AUDIT,
  DISEASE_MAPPING_AUDIT,
  SAMPLE_DISEASE_AUDIT,
  DISEASE_PALETTE_FILE,
  DATASET_PALETTE_FILE,
  PANEL_C_TIER_FILE
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
# 8. 创建输出目录
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
# 9. 读取正式 Seurat 对象
###############################################################################

message("")
message("读取正式 Seurat 对象：")
message(INPUT_RDS)


obj <- readRDS(
  INPUT_RDS
)


###############################################################################
# 10. Assay 检查
###############################################################################

if (!"RNA" %in% Assays(obj)) {

  stop(
    "Seurat 对象中不存在正式 RNA assay。"
  )
}


DefaultAssay(obj) <- "RNA"


###############################################################################
# 11. RNA layer 检查
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
# 12. Reduction 检查
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
# 13. Graph 检查
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
# 14. Metadata 检查
###############################################################################

required_metadata <- c(
  "CellID",
  "Original_Barcode",
  "DatasetID",
  "GSM",
  "SampleID",
  "SampleID_Global",
  "DonorID",
  "Disease",
  "disease_analyse",
  "Disease_Abbrev_Analysis",
  "Disease_original",
  "Disease_Class",
  "Disease_stage",
  "Grade",
  "Age",
  "Age_Status",
  "Sex",
  "Brain_Region",
  "sc_or_sn",
  "Tissue_Status",
  "Treatment_Status",
  "Experimental_Perturbation_Status",
  "Merge_Include",
  "User_Approved_Include",
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt",
  "percent.ribo",
  "scDblFinder.class",
  "scDblFinder.score",
  "BroadLineage",
  "MicrogliaGate",
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
# 15. 正式对象基础检查
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
# 16. Metadata 转为 data.table
###############################################################################

md <- data.table::as.data.table(
  obj[[]]
)


###############################################################################
# 17. CellID 唯一性
###############################################################################

if (anyDuplicated(
  md$CellID
) > 0L) {

  stop(
    "metadata 中 CellID 不是全局唯一。"
  )
}


###############################################################################
# 18. CellID 与 colnames(obj) 对应检查
###############################################################################

if (!identical(
  as.character(md$CellID),
  colnames(obj)
)) {

  stop(
    paste0(
      "metadata$CellID 与 colnames(obj) 未严格逐一对应。",
      "\n停止分析，避免细胞身份错位。"
    )
  )
}


###############################################################################
# 19. 核心字段 NA / 空值检查
###############################################################################

core_columns <- c(
  "CellID",
  "DatasetID",
  "SampleID_Global",
  "Disease",
  CLUSTER_COL
)


for (col_i in core_columns) {

  value_i <- md[[col_i]]

  bad_i <-
    is.na(value_i) |
    trimws(
      as.character(value_i)
    ) == ""

  if (any(bad_i)) {

    stop(
      paste0(
        "核心字段 ",
        col_i,
        " 存在 NA 或空字符串，共 ",
        sum(bad_i),
        " 个细胞。"
      )
    )
  }
}


###############################################################################
# 20. Disease17 类别检查
###############################################################################

observed_diseases <- unique(
  as.character(
    md$Disease
  )
)


if (!setequal(
  observed_diseases,
  DISEASE_ORDER_CONTRACT
)) {

  stop(
    paste0(
      "当前 Disease 类别与正式 Disease17 不一致。\n",
      "实际：",
      paste(
        sort(observed_diseases),
        collapse = ", "
      )
    )
  )
}


###############################################################################
# 21. Disease17 正式细胞数审计
###############################################################################

observed_disease_counts <- md[
  ,
  .(
    observed_cells = .N
  ),
  by = Disease
]


expected_disease_counts <- data.table::data.table(
  Disease = names(
    EXPECTED_DISEASE_COUNTS
  ),
  expected_cells = as.integer(
    EXPECTED_DISEASE_COUNTS
  )
)


disease_count_audit <- merge(
  expected_disease_counts,
  observed_disease_counts,
  by = "Disease",
  all = TRUE,
  sort = FALSE
)


disease_count_audit[
  ,
  difference :=
    observed_cells -
    expected_cells
]


disease_count_audit[
  ,
  disease_order :=
    match(
      Disease,
      DISEASE_ORDER_CONTRACT
    )
]


data.table::setorder(
  disease_count_audit,
  disease_order
)


disease_count_audit[
  ,
  disease_order := NULL
]


data.table::fwrite(
  disease_count_audit,
  file.path(
    TABLE_DIR,
    "Fig2DE_disease17_cell_count_audit.csv"
  )
)


if (
  any(
    is.na(
      disease_count_audit$expected_cells
    )
  ) ||
  any(
    is.na(
      disease_count_audit$observed_cells
    )
  ) ||
  any(
    disease_count_audit$difference != 0L
  )
) {

  stop(
    "Disease17 细胞数与正式数据契约不一致。"
  )
}


###############################################################################
# 22. 正式 cluster 检查
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
# 23. 建立 C00-C29 显示标签
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
    "无法把正式 cluster 转换成整数 0-29。"
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
# 24. 三大疾病类别完整性检查
###############################################################################

defined_diseases <- c(
  "Ctrl",
  TUMOR_DISEASES,
  NEUROLOGICAL_DISEASES
)


if (!setequal(
  defined_diseases,
  DISEASE_ORDER_CONTRACT
)) {

  stop(
    "Tumor / Neurological / Ctrl 分类没有完整覆盖 Disease17。"
  )
}


if (anyDuplicated(
  defined_diseases
) > 0L) {

  stop(
    "三大疾病类别定义中存在重复 Disease。"
  )
}


###############################################################################
# 25. 构建 Tumor / Neurological / Ctrl
###############################################################################

md[
  ,
  Disease_Group3 :=
    data.table::fcase(

      Disease %chin%
        TUMOR_DISEASES,
      "Tumor",

      Disease %chin%
        NEUROLOGICAL_DISEASES,
      "Neurological",

      Disease == "Ctrl",
      "Ctrl",

      default =
        NA_character_
    )
]


if (anyNA(
  md$Disease_Group3
)) {

  stop(
    "存在无法映射至 Tumor / Neurological / Ctrl 的 Disease。"
  )
}


md[
  ,
  Disease_Group3 :=
    factor(
      Disease_Group3,
      levels = GROUP_ORDER
    )
]


###############################################################################
# 26. Disease17 -> 三大类映射表
###############################################################################

group_mapping <- unique(
  md[
    ,
    .(
      Disease,
      Disease_Group3 =
        as.character(
          Disease_Group3
        )
    )
  ]
)


group_mapping[
  ,
  disease_order :=
    match(
      Disease,
      DISEASE_ORDER_CONTRACT
    )
]


data.table::setorder(
  group_mapping,
  disease_order
)


group_mapping[
  ,
  disease_order := NULL
]


###############################################################################
# 保存到 data
###############################################################################

data.table::fwrite(
  group_mapping,
  file.path(
    DATA_DIR,
    "01_Disease17_group_mapping.csv"
  )
)


###############################################################################
# 27. 读取 Fig2C 的正式 cluster 行顺序
###############################################################################

panelc_tier <- data.table::fread(
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
  colnames(panelc_tier)
)


if (length(missing_panelc_columns) > 0L) {

  stop(
    paste0(
      "Fig2C cluster_tier_map.csv 缺少列：",
      paste(
        missing_panelc_columns,
        collapse = ", "
      )
    )
  )
}


###############################################################################
# 28. 统一 Fig2C cluster 标签格式
###############################################################################

normalize_cluster_label <- function(x) {

  x <- as.character(x)

  out <- rep(
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


  out[c_style] <- sprintf(
    "C%02d",
    as.integer(
      sub(
        "^C",
        "",
        x[c_style]
      )
    )
  )


  out[numeric_style] <- sprintf(
    "C%02d",
    as.integer(
      x[numeric_style]
    )
  )


  if (anyNA(out)) {

    stop(
      "Fig2C cluster_tier_map.csv 中存在无法识别的 cluster 标签。"
    )
  }


  return(out)
}


panelc_tier[
  ,
  Cluster :=
    normalize_cluster_label(
      cluster
    )
]


###############################################################################
# Fig2C 行顺序必须正好30个 cluster
###############################################################################

if (nrow(panelc_tier) != 30L) {

  stop(
    paste0(
      "Fig2C cluster_tier_map.csv 应有30行，实际为 ",
      nrow(panelc_tier),
      " 行。"
    )
  )
}


if (anyDuplicated(
  panelc_tier$Cluster
) > 0L) {

  stop(
    "Fig2C cluster_tier_map.csv 存在重复 cluster。"
  )
}


if (!setequal(
  panelc_tier$Cluster,
  CLUSTER_LABELS
)) {

  stop(
    "Fig2C cluster 顺序文件与当前 C00-C29 不完全一致。"
  )
}


###############################################################################
# Fig2C 文件当前行序就是正式展示行序
###############################################################################

PANEL_C_CLUSTER_ORDER <- panelc_tier$Cluster


###############################################################################
# 29. Disease 列顺序
#
# 与你现有 Fig2C 脚本完全一致：
# Ctrl 第一列，其余 Disease 字母序
###############################################################################

all_disease_sorted <- sort(
  unique(
    as.character(
      md$Disease
    )
  )
)


DISEASE_ORDER_FIG2C <- c(
  "Ctrl",
  setdiff(
    all_disease_sorted,
    "Ctrl"
  )
)


###############################################################################
# 30. 重算 cluster × Disease presence
#
# 用于：
# 1. Fig2D 哪些位置需要绘制气泡
# 2. 独立验证当前对象与 Fig2C 冻结排序完全一致
###############################################################################

cluster_disease_grid <- data.table::CJ(
  Cluster = CLUSTER_LABELS,
  Disease = DISEASE_ORDER_FIG2C,
  unique = TRUE
)


cluster_disease_counts <- md[
  ,
  .(
    n_cells = .N
  ),
  by = .(
    Cluster,
    Disease
  )
]


cluster_disease_presence <- merge(
  cluster_disease_grid,
  cluster_disease_counts,
  by = c(
    "Cluster",
    "Disease"
  ),
  all.x = TRUE,
  sort = FALSE
)


cluster_disease_presence[
  is.na(n_cells),
  n_cells := 0L
]


cluster_disease_totals <- md[
  ,
  .(
    disease_total_cells = .N
  ),
  by = Disease
]


cluster_disease_presence <- merge(
  cluster_disease_presence,
  cluster_disease_totals,
  by = "Disease",
  all.x = TRUE,
  sort = FALSE
)


cluster_disease_presence[
  ,
  cell_pooled_pct :=
    100 *
    n_cells /
    disease_total_cells
]


cluster_disease_presence[
  ,
  present :=
    as.integer(
      n_cells >=
        MIN_CELLS_PRESENT
    )
]


###############################################################################
# 31. 独立重算 Fig2C tier
###############################################################################

cluster_coverage_check <- cluster_disease_presence[
  ,
  .(
    n_dis_present =
      sum(
        present
      ),

    total_cells =
      sum(
        n_cells
      )
  ),
  by = Cluster
]


cluster_coverage_check[
  ,
  tier :=
    data.table::fcase(

      n_dis_present >= PAN_MIN,
      "Pan-disease type",

      n_dis_present <= SPEC_MAX,
      "Disease-specific type",

      default =
        "Shared type"
    )
]


###############################################################################
# 32. 与 Fig2C 冻结结果交叉验证
###############################################################################

panelc_check <- merge(

  panelc_tier[
    ,
    .(
      Cluster,
      tier_panelc =
        as.character(
          tier
        ),
      n_dis_present_panelc =
        as.integer(
          n_dis_present
        ),
      total_cells_panelc =
        as.integer(
          total_cells
        )
    )
  ],

  cluster_coverage_check[
    ,
    .(
      Cluster,
      tier_current =
        as.character(
          tier
        ),
      n_dis_present_current =
        as.integer(
          n_dis_present
        ),
      total_cells_current =
        as.integer(
          total_cells
        )
    )
  ],

  by = "Cluster",
  all = TRUE,
  sort = FALSE
)


panelc_check[
  ,
  identical_tier :=
    tier_panelc ==
    tier_current
]


panelc_check[
  ,
  identical_n_dis_present :=
    n_dis_present_panelc ==
    n_dis_present_current
]


panelc_check[
  ,
  identical_total_cells :=
    total_cells_panelc ==
    total_cells_current
]


data.table::fwrite(
  panelc_check,
  file.path(
    TABLE_DIR,
    "Fig2DE_PanelC_order_validation.csv"
  )
)


if (
  anyNA(
    panelc_check$identical_tier
  ) ||
  anyNA(
    panelc_check$identical_n_dis_present
  ) ||
  anyNA(
    panelc_check$identical_total_cells
  ) ||
  !all(
    panelc_check$identical_tier
  ) ||
  !all(
    panelc_check$identical_n_dis_present
  ) ||
  !all(
    panelc_check$identical_total_cells
  )
) {

  stop(
    paste0(
      "当前对象重新计算得到的 cluster tier 与已经完成的 Fig2C 不一致。",
      "\n请检查：",
      file.path(
        TABLE_DIR,
        "Fig2DE_PanelC_order_validation.csv"
      )
    )
  )
}


###############################################################################
# 33. 构建 Fig2C / D / E 正式行序表
###############################################################################

cluster_order_table <- panelc_tier[
  ,
  .(
    Cluster,
    tier =
      as.character(
        tier
      ),
    n_dis_present =
      as.integer(
        n_dis_present
      ),
    total_cells =
      as.integer(
        total_cells
      )
  )
]


cluster_order_table[
  ,
  row_order :=
    seq_len(.N)
]


###############################################################################
# row_order = 1 表示图的最上方
###############################################################################

data.table::fwrite(
  cluster_order_table,
  file.path(
    DATA_DIR,
    "05_cluster_row_order.csv"
  )
)


###############################################################################
# 同时保存到 tables 便于审计
###############################################################################

data.table::fwrite(
  cluster_order_table,
  file.path(
    TABLE_DIR,
    "Fig2CDE_cluster_row_order.csv"
  )
)


###############################################################################
# 34. ggplot y level
#
# ggplot discrete y 的第一个 level 位于底部，
# 所以需要反转。
###############################################################################

y_levels <- rev(
  PANEL_C_CLUSTER_ORDER
)


###############################################################################
# 35. 构建 tier 分界线
###############################################################################

final_order <- data.table::copy(
  cluster_order_table
)


K <- nrow(
  final_order
)


final_order[
  ,
  y_pos :=
    K -
    .I +
    1L
]


previous_tier <- data.table::shift(
  final_order$tier
)


boundary_index <- which(
  !is.na(
    previous_tier
  ) &
    final_order$tier !=
    previous_tier
)


bounds <- data.table::data.table(
  y =
    final_order$y_pos[
      boundary_index
    ] +
    0.5
)


###############################################################################
# 36. tier 标签中心位置
###############################################################################

tier_mid <- final_order[
  ,
  .(
    y_mid =
      mean(
        y_pos
      ),

    y_start =
      min(
        y_pos
      )
  ),
  by = tier
]


data.table::setorder(
  tier_mid,
  y_start
)


tier_mid[
  ,
  y_start := NULL
]


###############################################################################
# 37. 保存 Cluster × Disease presence 中间表
###############################################################################

cluster_disease_presence_export <- data.table::copy(
  cluster_disease_presence
)


cluster_disease_presence_export[
  ,
  cluster_row_order :=
    match(
      Cluster,
      PANEL_C_CLUSTER_ORDER
    )
]


cluster_disease_presence_export[
  ,
  disease_column_order :=
    match(
      Disease,
      DISEASE_ORDER_FIG2C
    )
]


data.table::setorder(
  cluster_disease_presence_export,
  cluster_row_order,
  disease_column_order
)


data.table::fwrite(
  cluster_disease_presence_export,
  file.path(
    DATA_DIR,
    "04_cluster_disease_presence.csv"
  )
)


###############################################################################
# 38. SampleID_Global 映射唯一性检查
###############################################################################

sample_mapping_qc <- md[
  ,
  .(
    n_cells = .N,

    n_DatasetID =
      data.table::uniqueN(
        DatasetID
      ),

    n_Disease =
      data.table::uniqueN(
        Disease
      ),

    n_Disease_Group3 =
      data.table::uniqueN(
        Disease_Group3
      ),

    n_disease_analyse =
      data.table::uniqueN(
        disease_analyse
      )
  ),
  by = SampleID_Global
]


if (
  sample_mapping_qc[
    ,
    any(
      n_DatasetID != 1L
    )
  ]
) {

  stop(
    "存在同一 SampleID_Global 对应多个 DatasetID。"
  )
}


if (
  sample_mapping_qc[
    ,
    any(
      n_Disease != 1L
    )
  ]
) {

  stop(
    "存在同一 SampleID_Global 对应多个 Disease。"
  )
}


if (
  sample_mapping_qc[
    ,
    any(
      n_Disease_Group3 != 1L
    )
  ]
) {

  stop(
    "存在同一 SampleID_Global 对应多个疾病大类。"
  )
}


###############################################################################
# 39. 构建样本 metadata
###############################################################################

sample_metadata <- unique(
  md[
    ,
    .(
      SampleID_Global,
      DatasetID,
      Disease,
      Disease_Group3 =
        as.character(
          Disease_Group3
        )
    )
  ]
)


if (
  nrow(sample_metadata) !=
  data.table::uniqueN(
    md$SampleID_Global
  )
) {

  stop(
    "sample_metadata 构建异常。"
  )
}


###############################################################################
# 40. 样本细胞数
###############################################################################

sample_qc <- md[
  ,
  .(
    n_cells = .N
  ),
  by = SampleID_Global
]


sample_qc <- merge(
  sample_qc,
  sample_metadata,
  by = "SampleID_Global",
  all.x = TRUE,
  sort = FALSE
)


###############################################################################
# 41. >=20 cells 样本门槛
###############################################################################

sample_qc[
  ,
  eligible :=
    n_cells >=
    MIN_CELLS_PER_SAMPLE
]


eligible_samples <- sample_qc[
  eligible == TRUE,
  SampleID_Global
]


excluded_samples <- sample_qc[
  eligible == FALSE,
  SampleID_Global
]


if (length(
  eligible_samples
) == 0L) {

  stop(
    "没有任何 SampleID_Global 达到 >=20 cells 门槛。"
  )
}


###############################################################################
# 42. 保存样本 QC 中间数据
###############################################################################

sample_qc_export <- data.table::copy(
  sample_qc
)


sample_qc_export[
  ,
  Disease_Group3 :=
    as.character(
      Disease_Group3
    )
]


data.table::fwrite(
  sample_qc_export,
  file.path(
    DATA_DIR,
    "02_sample_qc.csv"
  )
)


###############################################################################
# 43. Disease 样本审计
###############################################################################

disease_sample_audit <- sample_qc[
  ,
  .(
    total_cells =
      sum(
        n_cells
      ),

    total_samples =
      .N,

    eligible_samples =
      sum(
        eligible
      ),

    excluded_samples =
      sum(
        !eligible
      ),

    eligible_cells =
      sum(
        n_cells[
          eligible
        ]
      )
  ),
  by = Disease
]


disease_sample_audit[
  ,
  disease_order :=
    match(
      Disease,
      DISEASE_ORDER_CONTRACT
    )
]


data.table::setorder(
  disease_sample_audit,
  disease_order
)


disease_sample_audit[
  ,
  disease_order := NULL
]


data.table::fwrite(
  disease_sample_audit,
  file.path(
    TABLE_DIR,
    "Fig2DE_samples_by_Disease_audit.csv"
  )
)


###############################################################################
# 44. 三大类样本审计
###############################################################################

three_group_sample_audit <- sample_qc[
  ,
  .(
    total_cells =
      sum(
        n_cells
      ),

    total_samples =
      .N,

    eligible_samples =
      sum(
        eligible
      ),

    excluded_samples =
      sum(
        !eligible
      ),

    eligible_cells =
      sum(
        n_cells[
          eligible
        ]
      )
  ),
  by = Disease_Group3
]


three_group_sample_audit[
  ,
  group_order :=
    match(
      Disease_Group3,
      GROUP_ORDER
    )
]


data.table::setorder(
  three_group_sample_audit,
  group_order
)


three_group_sample_audit[
  ,
  group_order := NULL
]


data.table::fwrite(
  three_group_sample_audit,
  file.path(
    TABLE_DIR,
    "Fig2DE_samples_by_three_groups_audit.csv"
  )
)


###############################################################################
# 45. 总体 audit
###############################################################################

analysis_audit_summary <- data.table::data.table(

  item = c(
    "total_cells",
    "total_samples",
    "eligible_samples_min20",
    "excluded_samples_lt20",
    "n_diseases",
    "n_three_groups",
    "n_clusters"
  ),

  value = c(
    nrow(md),

    data.table::uniqueN(
      md$SampleID_Global
    ),

    length(
      eligible_samples
    ),

    length(
      excluded_samples
    ),

    data.table::uniqueN(
      md$Disease
    ),

    data.table::uniqueN(
      md$Disease_Group3
    ),

    data.table::uniqueN(
      md$Cluster
    )
  )
)


data.table::fwrite(
  analysis_audit_summary,
  file.path(
    TABLE_DIR,
    "Fig2DE_analysis_audit_summary.csv"
  )
)


###############################################################################
# 46. 仅保留合格样本
###############################################################################

eligible_md <- md[
  SampleID_Global %chin%
    eligible_samples
]


###############################################################################
# 47. 每个样本 × Cluster 细胞数
###############################################################################

sample_cluster_counts <- eligible_md[
  ,
  .(
    cluster_cells = .N
  ),
  by = .(
    SampleID_Global,
    Cluster
  )
]


###############################################################################
# 48. 建立完整 Sample × 30 Cluster 网格
###############################################################################

sample_cluster_grid <- data.table::CJ(
  SampleID_Global =
    eligible_samples,

  Cluster =
    CLUSTER_LABELS,

  unique = TRUE
)


###############################################################################
# 49. 合并 count
###############################################################################

sample_cluster_usage <- merge(
  sample_cluster_grid,
  sample_cluster_counts,
  by = c(
    "SampleID_Global",
    "Cluster"
  ),
  all.x = TRUE,
  sort = FALSE
)


sample_cluster_usage[
  is.na(
    cluster_cells
  ),
  cluster_cells := 0L
]


###############################################################################
# 50. 合并样本 metadata
###############################################################################

sample_cluster_usage <- merge(

  sample_cluster_usage,

  sample_qc[
    eligible == TRUE,
    .(
      SampleID_Global,
      DatasetID,
      Disease,
      Disease_Group3,
      n_cells
    )
  ],

  by = "SampleID_Global",
  all.x = TRUE,
  sort = FALSE
)


###############################################################################
# 51. 样本内 Cluster usage
###############################################################################

sample_cluster_usage[
  ,
  usage_frac :=
    cluster_cells /
    n_cells
]


sample_cluster_usage[
  ,
  usage_pct :=
    usage_frac *
    100
]


###############################################################################
# 52. 检查每个样本30个cluster比例之和
###############################################################################

usage_sum_check <- sample_cluster_usage[
  ,
  .(
    n_clusters = .N,

    usage_sum =
      sum(
        usage_frac
      )
  ),
  by = SampleID_Global
]


if (
  any(
    usage_sum_check$n_clusters != 30L
  ) ||
  any(
    abs(
      usage_sum_check$usage_sum -
      1
    ) > 1e-10
  )
) {

  stop(
    "样本级 cluster usage 完整性检查失败。"
  )
}


###############################################################################
# 53. 保存最核心中间数据：
#     SampleID_Global × Cluster usage
###############################################################################

sample_cluster_usage_export <- data.table::copy(
  sample_cluster_usage
)


sample_cluster_usage_export[
  ,
  Disease_Group3 :=
    as.character(
      Disease_Group3
    )
]


sample_cluster_usage_export[
  ,
  cluster_row_order :=
    match(
      Cluster,
      PANEL_C_CLUSTER_ORDER
    )
]


data.table::setorder(
  sample_cluster_usage_export,
  SampleID_Global,
  cluster_row_order
)


data.table::fwrite(
  sample_cluster_usage_export,
  file.path(
    DATA_DIR,
    "03_sample_cluster_usage.csv"
  )
)


###############################################################################
# 54. 字体设置
###############################################################################

FONT <- "Arial"


###############################################################################
# 55. 与 Fig2C 一致的绘图主题
###############################################################################

theme_landscape <- ggplot2::theme_bw(
  base_family = FONT,
  base_size = 9
) +

  ggplot2::theme(

    panel.grid =
      ggplot2::element_blank(),

    panel.border =
      ggplot2::element_blank(),

    axis.line =
      ggplot2::element_blank(),

    axis.ticks =
      ggplot2::element_blank(),

    axis.title =
      ggplot2::element_blank(),

    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = 7.5,
        colour = "black"
      ),

    axis.text.y =
      ggplot2::element_text(
        size = 7,
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
        size = 8,
        colour = "grey30"
      ),

    legend.title =
      ggplot2::element_text(
        size = 8,
        colour = "black"
      ),

    legend.text =
      ggplot2::element_text(
        size = 7.5,
        colour = "black"
      ),

    legend.key.size =
      grid::unit(
        0.28,
        "cm"
      ),

    legend.box =
      "vertical",

    plot.margin =
      ggplot2::margin(
        6,
        12,
        4,
        4
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
      )
  )


###############################################################################
# 56. 正式图片保存函数
#
# PNG = 600 dpi
# PDF = vector
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
  # PNG
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
        "Arial PNG 输出失败，使用 sans 回退：",
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
          device =
            grDevices::cairo_pdf,
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
        "Arial/Cairo PDF 输出失败，使用 sans 回退：",
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
###############################################################################
#                              FIGURE 2D
###############################################################################
###############################################################################


###############################################################################
# 57. 每个 Disease 内对样本等权汇总 Cluster usage
###############################################################################

fig2d_sample_summary <- sample_cluster_usage[
  ,
  .(
    n_samples =
      data.table::uniqueN(
        SampleID_Global
      ),

    mean_usage_frac =
      mean(
        usage_frac
      ),

    mean_usage_pct =
      mean(
        usage_pct
      ),

    median_usage_pct =
      stats::median(
        usage_pct
      ),

    sd_usage_pct =
      stats::sd(
        usage_pct
      ),

    min_usage_pct =
      min(
        usage_pct
      ),

    max_usage_pct =
      max(
        usage_pct
      )
  ),
  by = .(
    Disease,
    Cluster
  )
]


###############################################################################
# 58. 合并 Fig2C 的 cell-count presence
#
# 同时保留：
#   cell_pooled_pct
#   sample_equal_mean_pct
#
# agent 后面可直接比较两种口径
###############################################################################

fig2d_source <- merge(

  cluster_disease_presence[
    ,
    .(
      Cluster,
      Disease,
      n_cells,
      disease_total_cells,
      cell_pooled_pct,
      present
    )
  ],

  fig2d_sample_summary,

  by = c(
    "Cluster",
    "Disease"
  ),

  all.x = TRUE,
  sort = FALSE
)


###############################################################################
# 59. 比例颜色分档
###############################################################################

fig2d_source[
  ,
  pct_bin :=
    cut(
      mean_usage_pct,
      breaks = c(
        -Inf,
        1,
        5,
        10,
        Inf
      ),
      labels = c(
        "<1%",
        "1-5%",
        "5-10%",
        ">10%"
      ),
      right = FALSE
    )
]


###############################################################################
# 60. 写入正式绘图顺序
###############################################################################

fig2d_source[
  ,
  cluster_row_order :=
    match(
      Cluster,
      PANEL_C_CLUSTER_ORDER
    )
]


fig2d_source[
  ,
  disease_column_order :=
    match(
      Disease,
      DISEASE_ORDER_FIG2C
    )
]


data.table::setorder(
  fig2d_source,
  cluster_row_order,
  disease_column_order
)


###############################################################################
# 61. 正式 Fig2D source CSV
###############################################################################

data.table::fwrite(
  fig2d_source,
  file.path(
    SOURCE_DIR,
    "Fig2d_Disease17_subcluster_percentage_source.csv"
  )
)


###############################################################################
# 62. 实际送入 ggplot 的数据
#
# 与 Fig2C 一致：
# 只有 >=20 cells 的 cluster × Disease 才绘制气泡
###############################################################################

fig2d_plot_data <- data.table::copy(
  fig2d_source[
    present == 1L
  ]
)


###############################################################################
# 因子顺序
###############################################################################

fig2d_plot_data[
  ,
  Disease :=
    factor(
      Disease,
      levels =
        DISEASE_ORDER_FIG2C
    )
]


fig2d_plot_data[
  ,
  Cluster :=
    factor(
      Cluster,
      levels =
        y_levels
    )
]


fig2d_plot_data[
  ,
  pct_bin :=
    factor(
      as.character(
        pct_bin
      ),
      levels = c(
        "<1%",
        "1-5%",
        "5-10%",
        ">10%"
      )
    )
]


###############################################################################
# 63. 保存 Fig2D 真正绘图数据到 data/
###############################################################################

fig2d_plot_export <- data.table::copy(
  fig2d_plot_data
)


fig2d_plot_export[
  ,
  Disease :=
    as.character(
      Disease
    )
]


fig2d_plot_export[
  ,
  Cluster :=
    as.character(
      Cluster
    )
]


fig2d_plot_export[
  ,
  pct_bin :=
    as.character(
      pct_bin
    )
]


data.table::setorder(
  fig2d_plot_export,
  cluster_row_order,
  disease_column_order
)


data.table::fwrite(
  fig2d_plot_export,
  file.path(
    DATA_DIR,
    "06_Fig2d_plot_data.csv"
  )
)


###############################################################################
# 64. Fig2D size scale 上限
#
# 至少保留 50% 图例范围，与 Fig2C/文献图例一致
###############################################################################

max_fig2d_pct <- max(
  fig2d_plot_data$mean_usage_pct,
  na.rm = TRUE
)


size_upper_limit <- max(
  50,
  ceiling(
    max_fig2d_pct
  )
)


###############################################################################
# 65. Fig2D 布局参数
###############################################################################

x_expand_disease <- ggplot2::expansion(
  add = c(
    2.2,
    0.3
  )
)


x_gutter <- -0.75

x_dash_l <- 0.6

x_dash_r_disease <-
  length(
    DISEASE_ORDER_FIG2C
  ) +
  0.4


###############################################################################
# 66. Fig2D
#
# color 与你提供的脚本完全一致：
#
# <1%    #C6DBEF
# 1-5%   #6BAED6
# 5-10%  #2171B5
# >10%   #08306B
###############################################################################

p_fig2d <-

  ggplot2::ggplot(
    fig2d_plot_data,
    ggplot2::aes(
      x = Disease,
      y = Cluster
    )
  ) +

  ggplot2::geom_point(
    ggplot2::aes(
      size =
        mean_usage_pct,

      colour =
        pct_bin
    )
  ) +

  ggplot2::geom_segment(

    data =
      bounds,

    ggplot2::aes(
      x =
        x_dash_l,

      xend =
        x_dash_r_disease,

      y =
        y,

      yend =
        y
    ),

    inherit.aes = FALSE,

    linetype =
      "dashed",

    colour =
      "grey30",

    linewidth =
      0.35
  ) +

  ggplot2::annotate(
    "text",

    x =
      x_gutter,

    y =
      tier_mid$y_mid,

    label =
      tier_mid$tier,

    angle =
      90,

    size =
      8 /
      .pt,

    colour =
      "black"
  ) +

  ggplot2::scale_size_continuous(
    name =
      "Percentage (%)",

    breaks =
      PCT_BREAKS,

    limits = c(
      0,
      size_upper_limit
    ),

    range = c(
      0.3,
      7
    )
  ) +

  ggplot2::scale_colour_manual(
    name =
      "Percentage",

    values =
      COL_BLUES,

    drop =
      FALSE
  ) +

  ggplot2::scale_x_discrete(
    expand =
      x_expand_disease
  ) +

  ggplot2::coord_cartesian(
    clip =
      "off"
  ) +

  ggplot2::labs(
    title =
      "Subcluster percentage",

    subtitle =
      "Bubble size and color encode sample-equal-weighted cluster proportion within each disease"
  ) +

  theme_landscape


###############################################################################
# 67. 保存 Fig2D
###############################################################################

save_formal_plot(
  plot_object =
    p_fig2d,

  basename =
    "Fig2d_subcluster_percentage",

  width =
    8.2,

  height =
    8.3
)


###############################################################################
###############################################################################
#                              FIGURE 2E
###############################################################################
###############################################################################


###############################################################################
# 68. safe odds
###############################################################################

safe_odds <- function(
    p,
    eps = OR_EPS
) {

  p_safe <- pmin(
    pmax(
      p,
      eps
    ),
    1 -
      eps
  )


  return(
    p_safe /
      (
        1 -
          p_safe
      )
  )
}


###############################################################################
# 69. safe Wilcoxon
###############################################################################

safe_wilcox <- function(
    x,
    y
) {

  x <- x[
    is.finite(
      x
    )
  ]

  y <- y[
    is.finite(
      y
    )
  ]


  if (
    length(x) < 2L ||
    length(y) < 2L
  ) {

    return(
      NA_real_
    )
  }


  if (
    length(
      unique(
        c(
          x,
          y
        )
      )
    ) <= 1L
  ) {

    return(
      1
    )
  }


  p_value <- tryCatch(

    {

      suppressWarnings(

        stats::wilcox.test(
          x = x,
          y = y,
          alternative = "two.sided",
          exact = FALSE
        )$p.value

      )
    },

    error = function(e) {

      NA_real_
    }
  )


  return(
    p_value
  )
}


###############################################################################
# 70. 逐疾病大类 × Cluster 计算 preference OR
###############################################################################

preference_result_list <- vector(
  mode = "list",
  length =
    length(
      GROUP_ORDER
    ) *
    length(
      CLUSTER_LABELS
    )
)


result_index <- 1L


for (group_i in GROUP_ORDER) {

  message(
    "Fig2E 正在计算：",
    group_i
  )


  for (cluster_i in CLUSTER_LABELS) {


    ###########################################################################
    # 当前 cluster
    ###########################################################################

    current_dt <- sample_cluster_usage[
      Cluster ==
        cluster_i
    ]


    ###########################################################################
    # 当前疾病大类
    ###########################################################################

    usage_in <- current_dt[
      Disease_Group3 ==
        group_i,
      usage_frac
    ]


    ###########################################################################
    # 其他两个疾病大类
    ###########################################################################

    usage_out <- current_dt[
      Disease_Group3 !=
        group_i,
      usage_frac
    ]


    ###########################################################################
    # 基本检查
    ###########################################################################

    if (
      length(
        usage_in
      ) == 0L ||
      length(
        usage_out
      ) == 0L
    ) {

      preference_result_list[[result_index]] <-
        data.table::data.table(

          Disease_Group3 =
            group_i,

          Cluster =
            cluster_i,

          n_samples_in =
            length(
              usage_in
            ),

          n_samples_out =
            length(
              usage_out
            ),

          mean_usage_in =
            NA_real_,

          mean_usage_out =
            NA_real_,

          mean_usage_pct_in =
            NA_real_,

          mean_usage_pct_out =
            NA_real_,

          odds_in =
            NA_real_,

          odds_out =
            NA_real_,

          OR =
            NA_real_,

          log2_OR =
            NA_real_,

          wilcox_p =
            NA_real_
        )


      result_index <-
        result_index +
        1L


      next
    }


    ###########################################################################
    # 样本等权平均
    ###########################################################################

    mean_in <- mean(
      usage_in,
      na.rm = TRUE
    )


    mean_out <- mean(
      usage_out,
      na.rm = TRUE
    )


    ###########################################################################
    # Odds
    ###########################################################################

    odds_in <- safe_odds(
      mean_in
    )


    odds_out <- safe_odds(
      mean_out
    )


    ###########################################################################
    # Preference OR
    ###########################################################################

    preference_or <-
      odds_in /
      odds_out


    ###########################################################################
    # 样本级 Wilcoxon
    ###########################################################################

    wilcox_p <- safe_wilcox(
      usage_in,
      usage_out
    )


    ###########################################################################
    # 保存结果
    ###########################################################################

    preference_result_list[[result_index]] <-
      data.table::data.table(

        Disease_Group3 =
          group_i,

        Cluster =
          cluster_i,

        n_samples_in =
          length(
            usage_in
          ),

        n_samples_out =
          length(
            usage_out
          ),

        mean_usage_in =
          mean_in,

        mean_usage_out =
          mean_out,

        mean_usage_pct_in =
          mean_in *
          100,

        mean_usage_pct_out =
          mean_out *
          100,

        odds_in =
          odds_in,

        odds_out =
          odds_out,

        OR =
          preference_or,

        log2_OR =
          log2(
            preference_or
          ),

        wilcox_p =
          wilcox_p
      )


    result_index <-
      result_index +
      1L
  }
}


###############################################################################
# 71. 合并 Fig2E
###############################################################################

fig2e_source <- data.table::rbindlist(
  preference_result_list
)


###############################################################################
# 应为 3 × 30 = 90 行
###############################################################################

if (
  nrow(
    fig2e_source
  ) != 90L
) {

  stop(
    paste0(
      "Fig2E 结果行数异常。",
      "\n预期 90，实际 ",
      nrow(
        fig2e_source
      )
    )
  )
}


###############################################################################
# 72. BH FDR
###############################################################################

fig2e_source[
  ,
  wilcox_FDR :=
    stats::p.adjust(
      wilcox_p,
      method = "BH"
    )
]


###############################################################################
# 73. FDR 显著性标签
#
# 仅保存供 Agent 解读，
# 不把它和文献的 OR preference 符号混为一谈。
###############################################################################

fig2e_source[
  ,
  FDR_significance :=
    data.table::fcase(

      is.na(
        wilcox_FDR
      ),
      "NA",

      wilcox_FDR <
        0.001,
      "***",

      wilcox_FDR <
        0.01,
      "**",

      wilcox_FDR <
        0.05,
      "*",

      default =
        "ns"
    )
]


###############################################################################
# 74. Preference 分类
###############################################################################

fig2e_source[
  ,
  preference_status :=
    data.table::fcase(

      is.na(
        OR
      ),
      "NA",

      OR >
        1.5,
      "OR > 1.5",

      OR >
        1,
      "1 < OR <= 1.5",

      default =
        "OR <= 1"
    )
]


###############################################################################
# 75. Fig2E 图内 preference 符号
###############################################################################

fig2e_source[
  ,
  preference_symbol :=
    data.table::fcase(

      is.na(
        OR
      ),
      "",

      OR >
        1.5,
      "++",

      OR >
        1,
      "+",

      default =
        "+/-"
    )
]


###############################################################################
# 76. 限制绘图颜色上限
#
# 原始 OR 不修改。
###############################################################################

fig2e_source[
  ,
  OR_plot :=
    pmin(
      OR,
      OR_PLOT_CAP
    )
]


###############################################################################
# 77. Fig2E 使用与 Fig2C 相同的行序
###############################################################################

fig2e_source[
  ,
  cluster_row_order :=
    match(
      Cluster,
      PANEL_C_CLUSTER_ORDER
    )
]


fig2e_source[
  ,
  group_column_order :=
    match(
      Disease_Group3,
      GROUP_ORDER
    )
]


data.table::setorder(
  fig2e_source,
  cluster_row_order,
  group_column_order
)


###############################################################################
# 78. 正式 Fig2E source CSV
###############################################################################

data.table::fwrite(
  fig2e_source,
  file.path(
    SOURCE_DIR,
    "Fig2e_Tumor_Neurological_Ctrl_preference_OR_source.csv"
  )
)


###############################################################################
# 79. Fig2E 绘图数据
###############################################################################

fig2e_plot_data <- data.table::copy(
  fig2e_source
)


fig2e_plot_data[
  ,
  Disease_Group3 :=
    factor(
      Disease_Group3,
      levels =
        GROUP_ORDER
    )
]


fig2e_plot_data[
  ,
  Cluster :=
    factor(
      Cluster,
      levels =
        y_levels
    )
]


###############################################################################
# 80. 保存 Fig2E 真正绘图数据
###############################################################################

fig2e_plot_export <- data.table::copy(
  fig2e_plot_data
)


fig2e_plot_export[
  ,
  Disease_Group3 :=
    as.character(
      Disease_Group3
    )
]


fig2e_plot_export[
  ,
  Cluster :=
    as.character(
      Cluster
    )
]


data.table::setorder(
  fig2e_plot_export,
  cluster_row_order,
  group_column_order
)


data.table::fwrite(
  fig2e_plot_export,
  file.path(
    DATA_DIR,
    "07_Fig2e_plot_data.csv"
  )
)


###############################################################################
# 81. Fig2E 布局参数
###############################################################################

x_expand_group <- ggplot2::expansion(
  add = c(
    2.2,
    0.3
  )
)


x_dash_r_group <-
  length(
    GROUP_ORDER
  ) +
  0.4


###############################################################################
# 82. Fig2E
###############################################################################

p_fig2e <-

  ggplot2::ggplot(
    fig2e_plot_data,
    ggplot2::aes(
      x =
        Disease_Group3,

      y =
        Cluster,

      fill =
        OR_plot
    )
  ) +

  ggplot2::geom_tile(
    colour =
      "white",

    linewidth =
      0.40
  ) +

  ggplot2::geom_text(
    ggplot2::aes(
      label =
        preference_symbol
    ),

    colour =
      "black",

    family =
      FONT,

    size =
      2.8
  ) +

  ggplot2::geom_segment(

    data =
      bounds,

    ggplot2::aes(
      x =
        x_dash_l,

      xend =
        x_dash_r_group,

      y =
        y,

      yend =
        y
    ),

    inherit.aes =
      FALSE,

    linetype =
      "dashed",

    colour =
      "grey30",

    linewidth =
      0.35
  ) +

  ggplot2::annotate(
    "text",

    x =
      x_gutter,

    y =
      tier_mid$y_mid,

    label =
      tier_mid$tier,

    angle =
      90,

    size =
      8 /
      .pt,

    colour =
      "black"
  ) +

  ggplot2::scale_fill_gradientn(
    name =
      "Preference OR",

    colours =
      OR_COLORS,

    limits = c(
      0,
      OR_PLOT_CAP
    ),

    breaks = c(
      0,
      1,
      1.5,
      2,
      3
    ),

    na.value =
      "#F5F5F5"
  ) +

  ggplot2::scale_x_discrete(
    expand =
      x_expand_group
  ) +

  ggplot2::coord_cartesian(
    clip =
      "off"
  ) +

  ggplot2::labs(
    title =
      "Subcluster preference",

    subtitle =
      "Preference OR across Tumor, Neurological and Ctrl"
  ) +

  theme_landscape


###############################################################################
# 83. 保存 Fig2E
###############################################################################

save_formal_plot(
  plot_object =
    p_fig2e,

  basename =
    "Fig2e_Tumor_Neurological_Ctrl_preference_OR",

  width =
    5.3,

  height =
    8.3
)


###############################################################################
# 84. 为 Agent 创建数据字典
###############################################################################

data_dictionary <- data.table::data.table(

  file = c(
    "01_Disease17_group_mapping.csv",
    "02_sample_qc.csv",
    "03_sample_cluster_usage.csv",
    "04_cluster_disease_presence.csv",
    "05_cluster_row_order.csv",
    "06_Fig2d_plot_data.csv",
    "07_Fig2e_plot_data.csv"
  ),

  level = c(
    "Disease",
    "Sample",
    "Sample x Cluster",
    "Disease x Cluster",
    "Cluster",
    "Disease x Cluster",
    "Disease group x Cluster"
  ),

  purpose = c(

    "Disease17 到 Tumor / Neurological / Ctrl 的固定映射",

    "每个 SampleID_Global 的细胞数及 >=20 cells 纳入状态",

    "每个合格 SampleID_Global × Cluster 的样本级 usage，是 Fig2D/E 的核心基础数据",

    "Cluster × Disease 的全细胞计数、pooled proportion 和 >=20 cells presence",

    "Fig2C / Fig2D / Fig2E 共用的冻结 cluster 行顺序",

    "Fig2D 实际传给 ggplot 的最终绘图数据",

    "Fig2E 实际传给 ggplot 的最终绘图数据"
  ),

  statistical_unit = c(
    "Disease",
    "SampleID_Global",
    "SampleID_Global",
    "Descriptive cell count",
    "Cluster",
    "SampleID_Global -> Disease equal-weighted",
    "SampleID_Global"
  ),

  important_columns = c(

    "Disease, Disease_Group3",

    "SampleID_Global, DatasetID, Disease, Disease_Group3, n_cells, eligible",

    "SampleID_Global, Disease, Disease_Group3, Cluster, cluster_cells, n_cells, usage_frac, usage_pct",

    "Cluster, Disease, n_cells, disease_total_cells, cell_pooled_pct, present",

    "Cluster, tier, n_dis_present, total_cells, row_order",

    "Cluster, Disease, n_cells, present, n_samples, mean_usage_pct, cell_pooled_pct, pct_bin, cluster_row_order, disease_column_order",

    "Disease_Group3, Cluster, n_samples_in, n_samples_out, mean_usage_pct_in, mean_usage_pct_out, OR, log2_OR, wilcox_p, wilcox_FDR, FDR_significance, preference_status, cluster_row_order"
  ),

  interpretation = c(

    "Tumor = LGG/GBM/MB/LCBM；Neurological = 12种神经系统疾病；Ctrl = Ctrl",

    "eligible=TRUE 表示样本至少包含20个小胶质细胞",

    "usage_frac 为0-1；usage_pct 为0-100；每个样本30个Cluster的usage_frac之和为1",

    "present=1 表示该Cluster在该Disease中全细胞计数 >=20；cell_pooled_pct仅用于描述与Fig2C存在性对照",

    "row_order=1 表示Fig2C/D/E最上方；该顺序直接继承已经完成的Fig2C",

    "Fig2D气泡大小和颜色使用mean_usage_pct；mean_usage_pct是Disease内对SampleID_Global等权平均，不是所有细胞直接pool",

    "OR>1表示相对偏向该疾病大类；OR>1.5为更强偏好；Wilcoxon/FDR为额外样本级统计，不等同于原论文cell-level Fisher检验"
  )
)


data.table::fwrite(
  data_dictionary,
  file.path(
    DATA_DIR,
    "08_analysis_data_dictionary.csv"
  )
)


###############################################################################
# 85. 保存一个 RDS 中间数据合集
#
# 方便以后 R / Agent 工作流整体加载
###############################################################################

intermediate_bundle <- list(

  disease_group_mapping =
    group_mapping,

  sample_qc =
    sample_qc_export,

  sample_cluster_usage =
    sample_cluster_usage_export,

  cluster_disease_presence =
    cluster_disease_presence_export,

  cluster_row_order =
    cluster_order_table,

  fig2d_plot_data =
    fig2d_plot_export,

  fig2e_plot_data =
    fig2e_plot_export,

  data_dictionary =
    data_dictionary
)


saveRDS(
  intermediate_bundle,
  file.path(
    DATA_DIR,
    "09_Fig2DE_intermediate_bundle.rds"
  )
)


###############################################################################
# 86. 分析方法说明
###############################################################################

method_note <- c(

  "Fig2D / Fig2E analysis note",

  "",

  "ROW ORDER",

  "Fig2D and Fig2E use exactly the same cluster row order as the completed Fig2C.",

  paste0(
    "Authoritative Fig2C row-order file: ",
    PANEL_C_TIER_FILE
  ),

  "The current object was independently re-audited against the Fig2C tier, disease coverage and total cell counts before plotting.",

  "",

  "FIG2D",

  "Only SampleID_Global with >=20 cells were included in sample-level usage analysis.",

  "Cluster usage was first calculated independently within every eligible SampleID_Global.",

  "Disease-level mean_usage_pct was calculated by equally averaging SampleID_Global values within the Disease.",

  "Cells from different samples were not directly pooled to calculate the formal Disease-level bubble value.",

  "The Fig2C >=20 pooled-cell presence definition is retained only to determine whether a bubble is displayed.",

  "Both sample_equal_mean_pct and cell_pooled_pct are retained in the intermediate/source data for audit.",

  "Fig2D color bins are <1%, 1-5%, 5-10%, >10%.",

  paste0(
    "Fig2D colors: ",
    paste(
      names(COL_BLUES),
      COL_BLUES,
      sep = "=",
      collapse = "; "
    )
  ),

  "",

  "FIG2E",

  "Disease17 was grouped into Tumor, Neurological and Ctrl.",

  "For every cluster, p_in is the equal-weighted mean sample-level cluster proportion within the target group.",

  "p_out is the equal-weighted mean sample-level cluster proportion among samples belonging to the other two groups.",

  "Preference OR = odds(p_in) / odds(p_out).",

  "OR > 1 indicates preference for the target group and OR > 1.5 indicates stronger preference.",

  "A sample-level two-sided Wilcoxon rank-sum test was additionally performed and all 90 tests were corrected together using BH-FDR.",

  "Wilcoxon/FDR is an additional sample-level inference and is not represented as the original publication's cell-count Fisher exact test.",

  "",

  paste0(
    "Minimum cells per sample: ",
    MIN_CELLS_PER_SAMPLE
  ),

  paste0(
    "Fig2C/D presence threshold: ",
    MIN_CELLS_PRESENT
  ),

  paste0(
    "Pan-disease threshold: >= ",
    PAN_MIN,
    " diseases"
  ),

  paste0(
    "Disease-specific threshold: <= ",
    SPEC_MAX,
    " diseases"
  ),

  paste0(
    "Formal cluster field: ",
    CLUSTER_COL
  ),

  paste0(
    "Random seed: ",
    SEED
  )
)


writeLines(
  method_note,
  con = file.path(
    TABLE_DIR,
    "Fig2DE_method_note.txt"
  )
)


###############################################################################
# 87. sessionInfo
###############################################################################

capture.output(
  sessionInfo(),
  file = file.path(
    TABLE_DIR,
    "Fig2DE_sessionInfo.txt"
  )
)


###############################################################################
# 88. 控制台输出
###############################################################################

message("")
message("============================================================")
message("Fig2D / Fig2E 分析完成")
message("============================================================")

message(
  "总细胞数：",
  nrow(md)
)

message(
  "总样本数：",
  data.table::uniqueN(
    md$SampleID_Global
  )
)

message(
  ">=20 cells 合格样本：",
  length(
    eligible_samples
  )
)

message(
  "<20 cells 排除样本：",
  length(
    excluded_samples
  )
)

message(
  "Disease 数：",
  data.table::uniqueN(
    md$Disease
  )
)

message(
  "Cluster 数：",
  data.table::uniqueN(
    md$Cluster
  )
)

message(
  "三大类：",
  paste(
    GROUP_ORDER,
    collapse = " | "
  )
)

message("")
message("Fig2C/D/E cluster行序已验证一致。")

message("")
message(
  "Figures: ",
  FIG_DIR
)

message(
  "Source data: ",
  SOURCE_DIR
)

message(
  "Tables: ",
  TABLE_DIR
)

message(
  "Intermediate data: ",
  DATA_DIR
)

message("")
message(
  "Fig2D：Fig2d_subcluster_percentage.png / .pdf"
)

message(
  "Fig2E：Fig2e_Tumor_Neurological_Ctrl_preference_OR.png / .pdf"
)

message("")
message("Agent 优先读取：")

message(
  file.path(
    DATA_DIR,
    "06_Fig2d_plot_data.csv"
  )
)

message(
  file.path(
    DATA_DIR,
    "07_Fig2e_plot_data.csv"
  )
)

message(
  file.path(
    DATA_DIR,
    "03_sample_cluster_usage.csv"
  )
)

message(
  file.path(
    DATA_DIR,
    "08_analysis_data_dictionary.csv"
  )
)

message("============================================================")

