# =============================================================================
# 20_findallmarkers_res06_from_merged_obj.R
# 目的：
#   1. 基于已冻结的 merged_obj.rds，在 RNA data layer 对全部细胞和全部
#      protein-coding genes 运行 resolution 0.6 的 FindAllMarkers；
#   2. 生成全部 marker、每簇 top10、簇统计、参数/输入/完成审计、日志和
#      sessionInfo；
#   3. 新结果验收通过后，将旧 resolution 0.4 的三个正式结果安全归档。
#
# 冻结参数：
#   assay = RNA, slot = data, test.use = wilcox, only.pos = TRUE,
#   min.pct = 0.1, logfc.threshold = 0.25, densify = FALSE,
#   全 424222 个细胞、全 18687 个 protein-coding genes，不抽样。
#
# 说明：本脚本不修改 merged_obj.rds，不重算聚类，不修改 nCount_RNA/
#   nFeature_RNA，也不修改疾病 metadata。disease_analyse 必须存在且与
#   Disease_Abbrev_Analysis 逐细胞一致。
# =============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

# -----------------------------------------------------------------------------
# 1. 冻结路径与参数
# -----------------------------------------------------------------------------
INPUT_PATH <- "E:/micro_output/v3_dims16/objects/merged_obj.rds"
INPUT_AUDIT_PATH <- "E:/micro_output/v3_dims16/objects/merged_obj_audit.csv"
SAMPLE_DISEASE_AUDIT_PATH <- "E:/micro_output/v3_dims16/objects/merged_obj_sample_disease_authority_audit.csv"
OUT_DIR <- "E:/micro_output/v3_dims16/findallmarkers"
INTERPRETATION_PATH <- file.path(
  "E:/micro_output/v3_dims16",
  paste0(intToUtf8(c(35299L, 35835L)), ".md")
)

RESOLUTION <- 0.6
RESOLUTION_LABEL <- "0.6"
EXPECTED_RESOLUTION_COLUMN <- "RNA_snn_res.0.6"
EXPECTED_ACTIVE_IDENT <- "seurat_clusters"
EXPECTED_CELLS <- 424222L
EXPECTED_GENES <- 18687L
EXPECTED_CLUSTERS <- 30L
EXPECTED_SAMPLES <- 453L

ASSAY <- "RNA"
DATA_LAYER <- "data"
TEST_USE <- "wilcox"
ONLY_POS <- TRUE
MIN_PCT <- 0.1
LOGFC_THRESHOLD <- 0.25
DENSIFY <- FALSE

OLD_RESULT_NAMES <- c(
  "all_markers_res0.4.csv",
  "top10_per_cluster_res0.4.csv",
  "cluster_stats_res0.4.csv"
)
ARCHIVE_DIR <- file.path(OUT_DIR, "archive", "res0.4_20260831")

required_input_paths <- c(INPUT_PATH, INPUT_AUDIT_PATH, SAMPLE_DISEASE_AUDIT_PATH)
missing_input_paths <- required_input_paths[!file.exists(required_input_paths)]
if (length(missing_input_paths) > 0L) {
  stop("缺少必需输入文件：", paste(missing_input_paths, collapse = "; "))
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("缺少 digest 包，无法完成文件完整性审计。")
}

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
STAGING_DIR <- file.path(OUT_DIR, paste0(".staging_res0.6_", stamp))
if (dir.exists(STAGING_DIR)) {
  stop("staging 目录已存在，为避免覆盖而停止：", STAGING_DIR)
}
dir.create(STAGING_DIR, recursive = TRUE, showWarnings = FALSE)

RESULT_NAMES <- c(
  "all_markers_res0.6.csv",
  "top10_per_cluster_res0.6.csv",
  "cluster_stats_res0.6.csv",
  "findallmarkers_res0.6_parameters.csv",
  "findallmarkers_res0.6_input_audit.csv",
  "findallmarkers_res0.6_archive_audit.csv",
  "findallmarkers_res0.6_completion_audit.csv",
  "findallmarkers_res0.6.log",
  "findallmarkers_res0.6_sessionInfo.txt"
)
FINAL_RESULT_PATHS <- file.path(OUT_DIR, RESULT_NAMES)
existing_final_results <- FINAL_RESULT_PATHS[file.exists(FINAL_RESULT_PATHS)]
if (length(existing_final_results) > 0L) {
  stop("新结果目标文件已存在，为避免覆盖而停止：",
       paste(existing_final_results, collapse = "; "))
}

old_paths <- file.path(OUT_DIR, OLD_RESULT_NAMES)
missing_old <- old_paths[!file.exists(old_paths)]
if (length(missing_old) > 0L) {
  stop("待归档的旧 res0.4 文件不完整，停止且不移动任何文件：",
       paste(missing_old, collapse = "; "))
}
archive_paths <- file.path(ARCHIVE_DIR, OLD_RESULT_NAMES)
existing_archive_targets <- archive_paths[file.exists(archive_paths)]
if (length(existing_archive_targets) > 0L) {
  stop("归档目标文件已存在，停止且保留原文件：",
       paste(existing_archive_targets, collapse = "; "))
}

log_path <- file.path(STAGING_DIR, "findallmarkers_res0.6.log")
log_con <- file(log_path, open = "wt", encoding = "UTF-8")

log_msg <- function(...) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                 paste0(..., collapse = ""))
  cat(line, "\n")
  writeLines(line, log_con, useBytes = TRUE)
  flush(log_con)
}

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE)
}

sort_cluster_levels <- function(x) {
  lev <- unique(as.character(x))
  lev_num <- suppressWarnings(as.numeric(lev))
  if (all(is.finite(lev_num))) lev[order(lev_num)] else sort(lev)
}

log_msg("=== FindAllMarkers res=0.6 from merged_obj start ===")
log_msg("staging directory: ", STAGING_DIR)

# -----------------------------------------------------------------------------
# 2. 读取并核验 merged_obj 及输入审计
# -----------------------------------------------------------------------------
log_msg("读取输入审计：", INPUT_AUDIT_PATH)
input_audit <- fread(INPUT_AUDIT_PATH)
if (!all(c("item", "value") %in% names(input_audit))) {
  stop("merged_obj_audit.csv 缺少 item/value 列。")
}
audit_map <- setNames(as.character(input_audit$value), as.character(input_audit$item))
if (!identical(unname(audit_map[["output_path"]]), INPUT_PATH)) {
  stop("输入审计 output_path 与正式 merged_obj 路径不一致。")
}
if (!identical(unname(audit_map[["sample_disease_authority_status"]]), "PASS")) {
  stop("输入对象的样本疾病权威审计不是 PASS。")
}
if (!identical(as.integer(audit_map[["n_cells"]]), EXPECTED_CELLS)) {
  stop("输入审计中的 n_cells 不符合冻结值。")
}
if (!identical(as.integer(audit_map[["n_genes_protein_coding"]]), EXPECTED_GENES)) {
  stop("输入审计中的 protein-coding gene 数不符合冻结值。")
}
if (!identical(as.integer(audit_map[["n_clusters_res0.6"]]), EXPECTED_CLUSTERS)) {
  stop("输入审计中的 res0.6 簇数不符合冻结值。")
}
frozen_input_sha256 <- tolower(trimws(unname(audit_map[["sha256"]])))
if (length(frozen_input_sha256) != 1L || is.na(frozen_input_sha256) ||
    !nzchar(frozen_input_sha256)) {
  stop("merged_obj_audit.csv 缺少冻结 sha256。")
}
observed_input_sha256 <- tolower(sha256_file(INPUT_PATH))
if (!identical(observed_input_sha256, frozen_input_sha256)) {
  stop("INPUT_PATH 当前 SHA256 与 merged_obj_audit.csv 冻结值不一致。")
}

sample_disease_audit <- fread(SAMPLE_DISEASE_AUDIT_PATH)
if (!all(c("SampleID_Global", "Status") %in% names(sample_disease_audit))) {
  stop("样本疾病权威审计缺少 SampleID_Global/Status 列。")
}
if (nrow(sample_disease_audit) != EXPECTED_SAMPLES ||
    anyDuplicated(sample_disease_audit$SampleID_Global) ||
    any(sample_disease_audit$Status != "PASS")) {
  stop("453 个样本的疾病权威审计未全部 PASS。")
}

log_msg("读取正式对象：", INPUT_PATH)
obj <- readRDS(INPUT_PATH)
if (!inherits(obj, "Seurat")) stop("输入对象不是 Seurat 对象。")
if (ncol(obj) != EXPECTED_CELLS) stop("输入细胞数不符合冻结值。")
if (!ASSAY %in% Assays(obj)) stop("输入对象缺少 RNA assay。")
if (!all(c("counts", DATA_LAYER) %in% Layers(obj[[ASSAY]]))) {
  stop("RNA assay 必须同时包含 counts 和 data layer。")
}
if (nrow(obj[[ASSAY]]) != EXPECTED_GENES) {
  stop("RNA assay 基因数不符合冻结值。")
}
if (anyDuplicated(rownames(obj[[ASSAY]])) || anyDuplicated(colnames(obj))) {
  stop("输入对象存在重复基因名或重复 CellID。")
}

required_meta <- c(
  "nCount_RNA", "nFeature_RNA", "RNA_snn_res.0.6", "seurat_clusters",
  "Disease", "Disease_Abbrev_Analysis", "disease_analyse"
)
missing_meta <- setdiff(required_meta, colnames(obj[[]]))
if (length(missing_meta) > 0L) {
  stop("输入对象缺少必须 metadata 列：", paste(missing_meta, collapse = ", "))
}

disease_analyse_values <- as.character(obj[["disease_analyse", drop = TRUE]])
disease_abbrev_values <- as.character(obj[["Disease_Abbrev_Analysis", drop = TRUE]])
if (length(disease_analyse_values) != EXPECTED_CELLS ||
    anyNA(disease_analyse_values) || any(trimws(disease_analyse_values) == "") ||
    anyNA(disease_abbrev_values) || any(trimws(disease_abbrev_values) == "")) {
  stop("disease_analyse 或 Disease_Abbrev_Analysis 含 NA/空值。")
}
if (!identical(disease_analyse_values, disease_abbrev_values)) {
  stop("disease_analyse 与 Disease_Abbrev_Analysis 未逐细胞 identical。")
}
log_msg("metadata disease_analyse: PRESENT and identical to Disease_Abbrev_Analysis")

res06_labels <- as.character(obj[[EXPECTED_RESOLUTION_COLUMN, drop = TRUE]])
seurat_labels <- as.character(obj[["seurat_clusters", drop = TRUE]])
if (length(res06_labels) != EXPECTED_CELLS || anyNA(res06_labels) ||
    length(seurat_labels) != EXPECTED_CELLS || anyNA(seurat_labels)) {
  stop("res0.6 或 seurat_clusters 标签长度错误或含 NA。")
}
if (!identical(res06_labels, seurat_labels)) {
  stop("seurat_clusters 与 RNA_snn_res.0.6 不完全一致。")
}
cluster_levels <- sort_cluster_levels(res06_labels)
if (length(cluster_levels) != EXPECTED_CLUSTERS) {
  stop("实际 res0.6 簇数不符合冻结值：", length(cluster_levels))
}
original_active_labels <- as.character(Idents(obj))
if (!identical(original_active_labels, seurat_labels)) {
  stop("输入对象的原 active.ident 与 seurat_clusters/res0.6 不一致；不覆盖且停止。")
}
if (!identical(DefaultAssay(obj), ASSAY)) {
  stop("输入对象的原 DefaultAssay 不是 RNA；不覆盖且停止。")
}

log_msg("输入核验 PASS：cells=", ncol(obj), "; genes=", nrow(obj[[ASSAY]]),
        "; clusters=", length(cluster_levels), "; assay=RNA; layer=data")

input_audit_out <- data.table(
  item = c(
    "input_path", "input_size_bytes", "input_sha256", "input_audit_path",
    "sample_disease_audit_path", "sample_disease_audit_status",
    "n_cells", "n_genes", "n_clusters_res0.6", "n_samples",
    "assay", "data_layer", "counts_layer", "active_ident",
    "metadata_disease_analyse_alias", "metadata_qc_note"
  ),
  value = c(
    INPUT_PATH,
    as.character(file.info(INPUT_PATH)$size),
    observed_input_sha256,
    INPUT_AUDIT_PATH,
    SAMPLE_DISEASE_AUDIT_PATH,
    "PASS",
    as.character(ncol(obj)),
    as.character(nrow(obj[[ASSAY]])),
    as.character(length(cluster_levels)),
    as.character(EXPECTED_SAMPLES),
    ASSAY,
    DATA_LAYER,
    "counts",
    "seurat_clusters",
    "PRESENT_AND_IDENTICAL",
    "nCount_RNA/nFeature_RNA retained unchanged; no recalculation performed"
  )
)
fwrite(input_audit_out, file.path(STAGING_DIR, "findallmarkers_res0.6_input_audit.csv"))

parameters_out <- data.table(
  parameter = c(
    "input_path", "assay", "data_layer", "active_ident", "n_cells",
    "n_genes", "n_clusters", "test.use", "only.pos", "min.pct",
    "logfc.threshold", "densify", "sampling", "random_seed_for_markers"
  ),
  value = c(
    INPUT_PATH, ASSAY, DATA_LAYER, EXPECTED_ACTIVE_IDENT,
    as.character(EXPECTED_CELLS), as.character(EXPECTED_GENES),
    as.character(EXPECTED_CLUSTERS), TEST_USE, as.character(ONLY_POS),
    as.character(MIN_PCT), as.character(LOGFC_THRESHOLD), as.character(DENSIFY),
    "none; all cells and all genes", "not applicable"
  )
)
fwrite(parameters_out, file.path(STAGING_DIR, "findallmarkers_res0.6_parameters.csv"))

# -----------------------------------------------------------------------------
# 3. 计算并验收 FindAllMarkers（不抽样、RNA data layer、densify=FALSE）
# -----------------------------------------------------------------------------
features_all <- rownames(obj[[ASSAY]])
if (length(features_all) != EXPECTED_GENES) stop("FindAllMarkers features 数量异常。")
log_msg("开始 FindAllMarkers：全 cells、全 genes；test.use=wilcox；densify=FALSE")
markers <- FindAllMarkers(
  object = obj,
  assay = ASSAY,
  features = features_all,
  slot = DATA_LAYER,
  test.use = TEST_USE,
  only.pos = ONLY_POS,
  min.pct = MIN_PCT,
  logfc.threshold = LOGFC_THRESHOLD,
  densify = DENSIFY,
  verbose = FALSE
)
markers <- as.data.table(markers)
required_marker_cols <- c("cluster", "gene", "p_val", "p_val_adj")
missing_marker_cols <- setdiff(required_marker_cols, names(markers))
if (length(missing_marker_cols) > 0L) {
  stop("FindAllMarkers 输出缺少列：", paste(missing_marker_cols, collapse = ", "))
}
if (nrow(markers) == 0L || anyNA(markers$cluster) || anyNA(markers$gene)) {
  stop("FindAllMarkers 未产生有效 marker 结果。")
}
marker_clusters <- sort_cluster_levels(markers$cluster)
if (!setequal(marker_clusters, cluster_levels)) {
  stop("FindAllMarkers 结果未覆盖全部 30 个 cluster。")
}

lfc_candidates <- c("avg_log2FC", "avg_logFC", "avg_log2", "avg_log1P", "logFC")
lfc_present <- lfc_candidates[lfc_candidates %in% names(markers)]
if (length(lfc_present) != 1L) {
  stop("FindAllMarkers 未能恰好识别一个有效 logFC 列；实际列：",
       paste(lfc_present, collapse = ", "))
}
lfc_col <- lfc_present[[1L]]
setorderv(markers, c("cluster", lfc_col, "p_val_adj"), c(1L, -1L, 1L))
all_markers_path <- file.path(STAGING_DIR, "all_markers_res0.6.csv")
fwrite(markers, all_markers_path)

markers_top <- copy(markers)
markers_top[, cluster := factor(as.character(cluster), levels = cluster_levels)]
top10 <- markers_top[, head(.SD, 10L), by = cluster]
top10[, cluster := as.character(cluster)]
top10_path <- file.path(STAGING_DIR, "top10_per_cluster_res0.6.csv")
fwrite(top10, top10_path)

cluster_stats <- data.table(
  seurat_clusters = cluster_levels,
  n_cells = as.integer(table(factor(seurat_labels, levels = cluster_levels)))
)
cluster_stats_path <- file.path(STAGING_DIR, "cluster_stats_res0.6.csv")
fwrite(cluster_stats, cluster_stats_path)

if (nrow(cluster_stats) != EXPECTED_CLUSTERS || nrow(top10) == 0L ||
    length(unique(top10$cluster)) != EXPECTED_CLUSTERS ||
    nrow(top10) != EXPECTED_CLUSTERS * 10L ||
    any(top10[, .N, by = cluster]$N != 10L)) {
  stop("cluster_stats 或 top10 结果不满足每个 cluster 恰好 10 行。")
}
marker_rows <- nrow(markers)
marker_genes <- length(unique(markers$gene))
marker_clusters_n <- length(unique(markers$cluster))
top10_rows <- nrow(top10)
log_msg("FindAllMarkers 验收 PASS：", nrow(markers), " rows; ",
        length(unique(markers$gene)), " marker genes; ",
        length(unique(markers$cluster)), " clusters")

# -----------------------------------------------------------------------------
# 4. 新结果文件完整性核验（仍在 staging，尚未发布）
# -----------------------------------------------------------------------------
staged_result_names <- c(
  "all_markers_res0.6.csv", "top10_per_cluster_res0.6.csv",
  "cluster_stats_res0.6.csv", "findallmarkers_res0.6_parameters.csv",
  "findallmarkers_res0.6_input_audit.csv"
)
staged_result_paths <- file.path(STAGING_DIR, staged_result_names)
if (any(!file.exists(staged_result_paths)) || any(file.size(staged_result_paths) <= 0)) {
  stop("staging 中存在缺失或空的新结果文件。")
}
if (file.size(all_markers_path) <= 0 || file.size(top10_path) <= 0 ||
    file.size(cluster_stats_path) <= 0) {
  stop("marker 输出文件写入失败。")
}

rm(markers_top, markers, features_all, input_audit, input_audit_out,
   parameters_out, sample_disease_audit)
invisible(gc())

# -----------------------------------------------------------------------------
# 5. 事务性归档旧 res0.4（三个确切文件）
# -----------------------------------------------------------------------------
log_msg("新结果 staging 验收完成，开始归档旧 res0.4 文件。")
if (!dir.exists(ARCHIVE_DIR)) {
  dir.create(ARCHIVE_DIR, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(ARCHIVE_DIR)) stop("无法创建归档目录：", ARCHIVE_DIR)
if (any(file.exists(archive_paths))) {
  stop("归档目标在移动前已存在，停止且不移动原文件。")
}

old_hashes <- vapply(old_paths, sha256_file, character(1))
old_sizes <- as.numeric(file.info(old_paths)$size)
moved_indices <- integer(0)
published <- character(0)
interpretation_existed <- FALSE
interpretation_backup_ready <- FALSE
interpretation_backup <- file.path(STAGING_DIR, "解读_original.md")
completion_stage_path <- file.path(STAGING_DIR, "findallmarkers_res0.6_completion_audit.csv")
completion_final_path <- file.path(OUT_DIR, basename(completion_stage_path))
completion_pass_stage_path <- file.path(STAGING_DIR, "findallmarkers_res0.6_completion_PASS.tmp.csv")

clear_staged_pass <- function() {
  for (p in c(completion_stage_path, completion_pass_stage_path)) {
    if (file.exists(p)) file.remove(p)
  }
}

restore_interpretation <- function() {
  if (!interpretation_backup_ready) return(character(0))
  if (interpretation_existed) {
    ok <- tryCatch(
      isTRUE(file.copy(interpretation_backup, INTERPRETATION_PATH, overwrite = TRUE)) &&
        identical(sha256_file(interpretation_backup), sha256_file(INTERPRETATION_PATH)),
      error = function(e) FALSE
    )
    if (!ok) return(paste0("解读.md 恢复失败：", INTERPRETATION_PATH))
  } else if (file.exists(INTERPRETATION_PATH)) {
    if (!isTRUE(tryCatch(file.remove(INTERPRETATION_PATH), error = function(e) FALSE))) {
      return(paste0("解读.md 回滚失败：无法移除本次新建文件：", INTERPRETATION_PATH))
    }
  }
  character(0)
}

rollback_archive <- function() {
  errors <- character(0)
  if (length(moved_indices) > 0L) {
    for (j in rev(moved_indices)) {
      if (file.exists(old_paths[j]) && file.exists(archive_paths[j])) {
        copied_hash <- tryCatch(sha256_file(archive_paths[j]), error = function(e) NA_character_)
        if (!is.na(copied_hash) && identical(copied_hash, old_hashes[j]) &&
            isTRUE(tryCatch(file.remove(archive_paths[j]), error = function(e) FALSE))) {
          next
        }
        errors <- c(errors, paste0("归档回滚路径同时存在且无法清理：", archive_paths[j]))
      } else if (file.exists(old_paths[j]) || !file.exists(archive_paths[j])) {
        errors <- c(errors, paste0("归档回滚路径异常：", archive_paths[j]))
      } else if (!isTRUE(tryCatch(file.rename(archive_paths[j], old_paths[j]),
                                  error = function(e) FALSE))) {
        errors <- c(errors, paste0("归档回滚失败：", archive_paths[j]))
      }
    }
  }
  restored_hashes <- tryCatch(
    vapply(old_paths, sha256_file, character(1)),
    error = function(e) rep(NA_character_, length(old_paths))
  )
  if (any(!file.exists(old_paths)) || any(file.exists(archive_paths)) ||
      !identical(restored_hashes, old_hashes)) {
    errors <- c(errors, "归档回滚后原文件存在性或 SHA256 核验失败")
  }
  if (length(errors) == 0L) moved_indices <<- integer(0)
  errors
}

rollback_published <- function() {
  errors <- character(0)
  if (length(published) > 0L) {
    for (nm in rev(published)) {
      src <- file.path(OUT_DIR, nm)
      dst <- file.path(STAGING_DIR, nm)
      if (file.exists(src)) {
        if (file.exists(dst)) {
          errors <- c(errors, paste0("新结果回滚目标已存在：", dst))
          next
        }
        src_hash <- tryCatch(sha256_file(src), error = function(e) NA_character_)
        if (is.na(src_hash) || !isTRUE(file.rename(src, dst))) {
          errors <- c(errors, paste0("新结果回滚失败：", src))
          next
        }
        dst_hash <- tryCatch(sha256_file(dst), error = function(e) NA_character_)
        if (is.na(dst_hash) || !identical(src_hash, dst_hash)) {
          errors <- c(errors, paste0("新结果回滚后 SHA256 核验失败：", dst))
        }
      }
    }
  }
  if (length(published) > 0L && any(file.exists(file.path(OUT_DIR, published)))) {
    errors <- c(errors, "新结果回滚后仍有文件留在正式目录")
  }
  if (length(errors) == 0L) published <<- character(0)
  errors
}

rollback_all <- function() {
  errors <- character(0)
  interpretation_errors <- tryCatch(restore_interpretation(),
                                    error = function(e) paste0("解读回滚异常：", conditionMessage(e)))
  published_errors <- tryCatch(rollback_published(),
                               error = function(e) paste0("新结果回滚异常：", conditionMessage(e)))
  archive_errors <- tryCatch(rollback_archive(),
                             error = function(e) paste0("旧归档回滚异常：", conditionMessage(e)))
  errors <- c(errors, interpretation_errors, published_errors, archive_errors)
  unique(errors[nzchar(errors)])
}

failure_handler <- function(reason) {
  rollback_errors <- tryCatch(rollback_all(),
                              error = function(e) paste0("统一回滚异常：", conditionMessage(e)))
  clear_errors <- tryCatch({ clear_staged_pass(); character(0) },
                           error = function(e) paste0("暂存 PASS 清理异常：", conditionMessage(e)))
  all_errors <- unique(c(rollback_errors, clear_errors))
  if (length(all_errors) > 0L) {
    stop(reason, "；回滚汇总：", paste(all_errors, collapse = " | "))
  }
  stop(reason, "；统一回滚 PASS")
}

publish_one <- function(nm) {
  src <- file.path(STAGING_DIR, nm)
  dst <- file.path(OUT_DIR, nm)
  if (!file.exists(src) || file.exists(dst)) return(FALSE)
  ok <- isTRUE(file.rename(src, dst))
  # 极少数文件系统可能在返回 FALSE 时已经完成移动；先纳入回滚清单。
  if (ok || (file.exists(dst) && !file.exists(src))) {
    published <<- c(published, nm)
  }
  ok
}

for (i in seq_along(old_paths)) {
  ok <- tryCatch(isTRUE(file.rename(old_paths[i], archive_paths[i])),
                 error = function(e) FALSE)
  if (!isTRUE(ok)) {
    if (file.exists(archive_paths[i])) {
      moved_indices <- c(moved_indices, i)
    }
    failure_handler(paste0("归档移动失败：", old_paths[i]))
  }
  moved_indices <- c(moved_indices, i)
}

archive_hashes <- tryCatch(
  vapply(archive_paths, sha256_file, character(1)),
  error = function(e) rep(NA_character_, length(archive_paths))
)
archive_audit <- tryCatch(data.table(
  source_path = old_paths,
  archive_path = archive_paths,
  source_size_bytes = old_sizes,
  source_sha256 = old_hashes,
  archive_exists = file.exists(archive_paths),
  source_absent_after_move = !file.exists(old_paths),
  archive_sha256 = archive_hashes,
  status = "PASS"
), error = function(e) NULL)
if (is.null(archive_audit) || any(!archive_audit$archive_exists) ||
    anyNA(archive_audit$archive_sha256) ||
    !all(archive_audit$source_sha256 == archive_audit$archive_sha256) ||
    any(!archive_audit$source_absent_after_move)) {
  failure_handler("旧 res0.4 归档完整性核验失败")
}
archive_audit_path <- file.path(STAGING_DIR, "findallmarkers_res0.6_archive_audit.csv")
archive_audit_write_ok <- isTRUE(tryCatch({
  fwrite(archive_audit, archive_audit_path)
  TRUE
}, error = function(e) FALSE))
if (!archive_audit_write_ok) failure_handler("archive_audit 写入失败")
log_ok <- isTRUE(tryCatch({
  log_msg("旧 res0.4 三个文件归档核验 PASS：", ARCHIVE_DIR)
  TRUE
}, error = function(e) FALSE))
if (!log_ok) failure_handler("归档日志写入失败")

# -----------------------------------------------------------------------------
# 6. 发布事务：先发布结果，再发布日志；失败时新结果和旧归档均回滚
# -----------------------------------------------------------------------------
session_ok <- isTRUE(tryCatch({
  capture.output(sessionInfo(), file = file.path(STAGING_DIR, "findallmarkers_res0.6_sessionInfo.txt"))
  TRUE
}, error = function(e) FALSE))
if (!session_ok) failure_handler("sessionInfo 写入失败")
log_ok <- isTRUE(tryCatch({
  log_msg("归档和 staging 结果均已核验，准备发布主结果文件。")
  TRUE
}, error = function(e) FALSE))
if (!log_ok) failure_handler("发布前日志写入失败")
close_ok <- isTRUE(tryCatch({ close(log_con); TRUE }, error = function(e) FALSE))
if (!close_ok) failure_handler("日志连接关闭失败")

base_publish_names <- c(
  "all_markers_res0.6.csv", "top10_per_cluster_res0.6.csv",
  "cluster_stats_res0.6.csv", "findallmarkers_res0.6_parameters.csv",
  "findallmarkers_res0.6_input_audit.csv", "findallmarkers_res0.6_archive_audit.csv",
  "findallmarkers_res0.6_sessionInfo.txt"
)
for (nm in base_publish_names) {
  publish_ok <- isTRUE(tryCatch(publish_one(nm), error = function(e) FALSE))
  if (!publish_ok) failure_handler(paste0("发布新结果失败：", nm))
}

# 在修改既有解读前做可验证备份；失败时恢复原文件。
interpretation_existed <- file.exists(INTERPRETATION_PATH)
interpretation_backup <- file.path(STAGING_DIR, "解读_original.md")
if (interpretation_existed) {
  interpretation_backup_ok <- isTRUE(tryCatch({
    file.copy(INTERPRETATION_PATH, interpretation_backup, overwrite = FALSE) &&
      identical(sha256_file(INTERPRETATION_PATH), sha256_file(interpretation_backup))
  }, error = function(e) FALSE))
  if (!interpretation_backup_ok) failure_handler("解读.md 备份失败")
}
interpretation_backup_ready <- TRUE

# -----------------------------------------------------------------------------
# 7. 结果发布后追加解读；先发布 READY 审计和日志，最后才形成 completion PASS
# -----------------------------------------------------------------------------
interpretation_section <- tryCatch(c(
  "", "## FindAllMarkers resolution 0.6（merged_obj；全 protein-coding genes）", "",
  "| 项目 | 内容 |", "|---|---|",
  "| 分析目的 | 识别全局 resolution 0.6 的 30 个 cluster 的阳性差异表达 marker，支持后续细胞类型和状态解释。 |",
  paste0("| 输入数据 | `", INPUT_PATH, "`；", EXPECTED_CELLS, " 个细胞、", EXPECTED_GENES, " 个 protein-coding genes、", EXPECTED_SAMPLES, " 个样本。 |"),
  "| 方法 | Seurat `FindAllMarkers`；RNA assay 的 `data` layer；Wilcoxon；`only.pos=TRUE`；`min.pct=0.1`；`logfc.threshold=0.25`；`densify=FALSE`；不抽样。 |",
  paste0("| 结果 | 输出全部 marker（", marker_rows, " 行、", marker_genes, " 个 marker gene）、每簇 top10（", top10_rows, " 行）和 cluster 统计；结果覆盖全部 ", EXPECTED_CLUSTERS, " 个 cluster。 |"),
  "| 生物学意义 | marker 集合用于判断各 cluster 的谱系、状态和潜在混合信号，并为后续注释与疾病亚型分析提供依据。 |",
  "| 限制 | 单细胞层面的 Wilcoxon marker 主要用于探索和注释；正式疾病差异仍需在 donor/sample 层面使用未整合表达建模。 |",
  paste0("| 旧结果归档 | 旧 res0.4 三个正式文件已完整移动至 `", ARCHIVE_DIR, "`，文件未删除，仅移动归档；`cluster_stats_res0.3.csv` 未移动。 |"),
  ""
), error = function(e) NULL)
if (is.null(interpretation_section)) failure_handler("阶段解读内容形成失败")
if (!isTRUE(tryCatch({
  cat(paste0(interpretation_section, collapse = "\n"), "\n",
      file = INTERPRETATION_PATH, append = TRUE, sep = "")
  TRUE
}, error = function(e) FALSE))) {
  failure_handler("解读.md 发布失败")
}

completion_audit <- tryCatch(data.table(
  item = c(
    "status", "completed_at", "input_path", "input_sha256", "n_cells", "n_genes",
    "n_clusters_res0.6", "marker_rows", "marker_genes", "marker_clusters", "top10_rows",
    "data_layer", "test.use", "only.pos", "min.pct", "logfc.threshold",
    "densify", "sampling", "archive_status", "archive_dir",
    "disease_analyse_alias_note", "object_modified"
  ),
  value = c(
    "READY", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), INPUT_PATH, observed_input_sha256,
    as.character(EXPECTED_CELLS), as.character(EXPECTED_GENES),
    as.character(EXPECTED_CLUSTERS), as.character(marker_rows), as.character(marker_genes),
    as.character(marker_clusters_n), as.character(top10_rows), DATA_LAYER,
    TEST_USE, as.character(ONLY_POS), as.character(MIN_PCT), as.character(LOGFC_THRESHOLD),
    as.character(DENSIFY), "none", "PASS", ARCHIVE_DIR,
    "present; identical to Disease_Abbrev_Analysis; no object change",
    "FALSE"
  )
), error = function(e) NULL)
if (is.null(completion_audit)) failure_handler("completion READY 审计对象形成失败")
if (!isTRUE(tryCatch({
  fwrite(completion_audit, completion_stage_path)
  TRUE
}, error = function(e) FALSE))) {
  failure_handler("completion READY 审计形成失败")
}

# READY 审计先随结果发布；只有日志也成功发布后才把 READY 更新为 PASS。
if (!file.exists(completion_stage_path) || file.exists(completion_final_path)) {
  failure_handler("completion READY 发布前检查失败")
}
if (!isTRUE(publish_one(basename(completion_stage_path)))) {
  failure_handler("completion READY 审计发布失败")
}

# 最终日志只在所有主结果和 completion READY 审计已发布后进入正式目录。
final_log_path <- file.path(OUT_DIR, "findallmarkers_res0.6.log")
stage_log_path <- file.path(STAGING_DIR, "findallmarkers_res0.6.log")
if (!isTRUE(publish_one(basename(final_log_path)))) {
  failure_handler("最终日志发布失败")
}

# 到此为止全部正式结果、READY 审计、sessionInfo、归档审计和最终日志均已发布；
# 现在才形成并覆盖 completion PASS。任何失败都会回滚且不留下 PASS 文件。
completion_pass_update_ok <- isTRUE(tryCatch({
  completion_audit[item == "status", value := "PASS"]
  TRUE
}, error = function(e) FALSE))
if (!completion_pass_update_ok) failure_handler("completion PASS 状态更新失败")
completion_pass_write_ok <- isTRUE(tryCatch({
  fwrite(completion_audit, completion_pass_stage_path)
  TRUE
}, error = function(e) FALSE))
completion_pass_copy_ok <- isTRUE(tryCatch(
  file.copy(completion_pass_stage_path, completion_final_path, overwrite = TRUE),
  error = function(e) FALSE
))
completion_pass_hash_ok <- isTRUE(tryCatch(
  identical(sha256_file(completion_pass_stage_path), sha256_file(completion_final_path)),
  error = function(e) FALSE
))
if (!completion_pass_write_ok || !completion_pass_copy_ok || !completion_pass_hash_ok) {
  failure_handler("completion PASS 最终形成失败")
}

if (!isTRUE(tryCatch({
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
             "=== FindAllMarkers res=0.6 ALL DONE; transaction PASS ===\n"),
      file = final_log_path, append = TRUE, sep = "")
  TRUE
}, error = function(e) FALSE))) {
  failure_handler("最终日志收尾失败")
}

cleanup_ok <- isTRUE(tryCatch({
  if (dir.exists(STAGING_DIR) &&
      identical(basename(STAGING_DIR), paste0(".staging_res0.6_", stamp))) {
    unlink(STAGING_DIR, recursive = TRUE, force = TRUE)
  }
  !dir.exists(STAGING_DIR)
}, error = function(e) FALSE))
if (!cleanup_ok) failure_handler("本次精确 staging 目录清理失败")
