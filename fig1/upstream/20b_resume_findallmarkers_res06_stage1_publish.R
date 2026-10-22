# =============================================================================
# 20b_resume_findallmarkers_res06_stage1_publish.R
# 目的：从已验收的 stage1 结果恢复 res0.6 发布，不重新读取对象或运行
# FindAllMarkers。只有全部输入、CSV 结构、参数、哈希和旧结果归档核验通过后，
# 才进行旧 res0.4 归档、res0.6 发布、解读和 completion PASS。
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
})
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("缺少 digest 包，无法进行 SHA256 核验。")
}

INPUT_PATH <- "E:/micro_output/v3_dims16/objects/merged_obj.rds"
INPUT_AUDIT_PATH <- "E:/micro_output/v3_dims16/objects/merged_obj_audit.csv"
SAMPLE_DISEASE_AUDIT_PATH <- "E:/micro_output/v3_dims16/objects/merged_obj_sample_disease_authority_audit.csv"
OUT_DIR <- "E:/micro_output/v3_dims16/findallmarkers"
STAGING_DIR <- file.path(OUT_DIR, ".staging_res0.6_20260831_174739")
ARCHIVE_DIR <- file.path(OUT_DIR, "archive", "res0.4_20260831")
INTERPRETATION_PATH <- file.path(
  "E:/micro_output/v3_dims16",
  paste0(intToUtf8(c(35299L, 35835L)), ".md")
)

OLD_NAMES <- c(
  "all_markers_res0.4.csv",
  "top10_per_cluster_res0.4.csv",
  "cluster_stats_res0.4.csv"
)
OLD_PATHS <- file.path(OUT_DIR, OLD_NAMES)
ARCHIVE_PATHS <- file.path(ARCHIVE_DIR, OLD_NAMES)
BASE_NAMES <- c(
  "all_markers_res0.6.csv",
  "top10_per_cluster_res0.6.csv",
  "cluster_stats_res0.6.csv",
  "findallmarkers_res0.6_parameters.csv",
  "findallmarkers_res0.6_input_audit.csv",
  "findallmarkers_res0.6_archive_audit.csv",
  "findallmarkers_res0.6_sessionInfo.txt"
)
COMPLETION_NAME <- "findallmarkers_res0.6_completion_audit.csv"
LOG_NAME <- "findallmarkers_res0.6.log"
FINAL_NAMES <- c(BASE_NAMES, COMPLETION_NAME, LOG_NAME)

EXPECTED_CELLS <- 424222L
EXPECTED_GENES <- 18687L
EXPECTED_CLUSTERS <- 30L
EXPECTED_SAMPLES <- 453L
EXPECTED_CLUSTER_LEVELS <- as.character(0:29)
EXPECTED_MARKER_ROWS <- 64007L
EXPECTED_MARKER_GENES <- 11357L
EXPECTED_TOP10_ROWS <- 300L

stop_preflight <- function(...) stop(paste0(...))
need_file <- function(path, label = path) {
  if (!file.exists(path) || isTRUE(file.info(path)$size <= 0)) {
    stop_preflight("缺少或为空：", label, "；", path)
  }
}
read_csv <- function(path, label) {
  out <- tryCatch(fread(path, na.strings = c("", "NA")),
                  error = function(e) NULL)
  if (is.null(out)) stop_preflight("CSV 读取失败：", label, "；", path)
  out
}
audit_value <- function(x, item, label) {
  if (!("item" %in% names(x)) || !("value" %in% names(x))) {
    stop_preflight("审计缺少 item/value：", label)
  }
  y <- x[x[["item"]] == item, x[["value"]]]
  if (length(y) != 1L || is.na(y[[1L]])) {
    stop_preflight("审计项缺失或不唯一：", label, " / ", item)
  }
  as.character(y[[1L]])
}
sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE)
}
same_sha <- function(path_a, path_b) {
  isTRUE(tryCatch(identical(sha256_file(path_a), sha256_file(path_b)),
                  error = function(e) FALSE))
}
cluster_set_ok <- function(x) {
  vals <- suppressWarnings(as.integer(as.character(x)))
  length(vals) == EXPECTED_CLUSTERS &&
    identical(sort(unique(vals)), seq_len(EXPECTED_CLUSTERS) - 1L)
}
param_value <- function(x, parameter) audit_value(x, parameter, "parameters")

# -----------------------------------------------------------------------------
# 1. 只读预检：stage1、输入审计、参数和结果结构
# -----------------------------------------------------------------------------
for (p in c(INPUT_PATH, INPUT_AUDIT_PATH, SAMPLE_DISEASE_AUDIT_PATH,
            STAGING_DIR, OLD_PATHS)) {
  if (length(p) == 1L && dir.exists(p)) next
  need_file(p)
}
if (!dir.exists(STAGING_DIR)) stop_preflight("stage1 目录不存在：", STAGING_DIR)
if (!dir.exists(ARCHIVE_DIR)) {
  if (!isTRUE(tryCatch(dir.create(ARCHIVE_DIR, recursive = TRUE),
                       error = function(e) FALSE)) && !dir.exists(ARCHIVE_DIR)) {
    stop_preflight("无法创建归档目录：", ARCHIVE_DIR)
  }
}
if (any(file.exists(ARCHIVE_PATHS))) {
  stop_preflight("归档目标已存在，停止：", paste(ARCHIVE_PATHS[file.exists(ARCHIVE_PATHS)], collapse = "; "))
}
if (any(file.exists(file.path(OUT_DIR, FINAL_NAMES)))) {
  stop_preflight("res0.6 正式目标已有文件，停止：",
                 paste(file.path(OUT_DIR, FINAL_NAMES)[file.exists(file.path(OUT_DIR, FINAL_NAMES))], collapse = "; "))
}

stage_paths <- file.path(STAGING_DIR, c(BASE_NAMES, LOG_NAME))
for (p in stage_paths) need_file(p, "stage1 文件")
stage_marker <- read_csv(file.path(STAGING_DIR, "all_markers_res0.6.csv"), "stage1 markers")
stage_top10 <- read_csv(file.path(STAGING_DIR, "top10_per_cluster_res0.6.csv"), "stage1 top10")
stage_stats <- read_csv(file.path(STAGING_DIR, "cluster_stats_res0.6.csv"), "stage1 cluster stats")
stage_input <- read_csv(file.path(STAGING_DIR, "findallmarkers_res0.6_input_audit.csv"), "stage1 input audit")
stage_param <- read_csv(file.path(STAGING_DIR, "findallmarkers_res0.6_parameters.csv"), "stage1 parameters")
stage_archive <- read_csv(file.path(STAGING_DIR, "findallmarkers_res0.6_archive_audit.csv"), "stage1 archive audit")

if (nrow(stage_marker) != EXPECTED_MARKER_ROWS ||
    !all(c("cluster", "gene", "p_val", "p_val_adj") %in% names(stage_marker)) ||
    length(unique(as.character(stage_marker$gene))) != EXPECTED_MARKER_GENES ||
    !cluster_set_ok(stage_marker$cluster)) {
  stop_preflight("stage1 markers 结构/行数/基因数/cluster 覆盖不符合冻结值。")
}
lfc_candidates <- c("avg_log2FC", "avg_logFC", "avg_log2", "avg_log1P", "logFC")
lfc_columns <- intersect(lfc_candidates, names(stage_marker))
if (length(lfc_columns) != 1L) {
  stop_preflight("stage1 marker 必须恰好含一个有效 logFC 列；实际：",
                 paste(lfc_columns, collapse = ", "))
}
if (nrow(stage_top10) != EXPECTED_TOP10_ROWS ||
    !("cluster" %in% names(stage_top10)) || !cluster_set_ok(stage_top10$cluster) ||
    any(stage_top10[, .N, by = cluster]$N != 10L)) {
  stop_preflight("stage1 top10 未覆盖 30 簇或每簇不恰好 10 行。")
}
if (!all(c("seurat_clusters", "n_cells") %in% names(stage_stats)) ||
    nrow(stage_stats) != EXPECTED_CLUSTERS ||
    !cluster_set_ok(stage_stats$seurat_clusters) ||
    sum(as.numeric(stage_stats$n_cells)) != EXPECTED_CELLS) {
  stop_preflight("stage1 cluster_stats 结构或细胞总数不符合冻结值。")
}

merged_audit <- read_csv(INPUT_AUDIT_PATH, "merged_obj audit")
sample_audit <- read_csv(SAMPLE_DISEASE_AUDIT_PATH, "sample disease audit")
if (audit_value(merged_audit, "output_path", "merged_obj audit") != INPUT_PATH ||
    audit_value(merged_audit, "n_cells", "merged_obj audit") != as.character(EXPECTED_CELLS) ||
    audit_value(merged_audit, "n_genes_protein_coding", "merged_obj audit") != as.character(EXPECTED_GENES) ||
    audit_value(merged_audit, "n_clusters_res0.6", "merged_obj audit") != as.character(EXPECTED_CLUSTERS)) {
  stop_preflight("merged_obj audit 冻结值不匹配。")
}
frozen_input_sha <- audit_value(merged_audit, "sha256", "merged_obj audit")
current_input_sha <- tryCatch(sha256_file(INPUT_PATH), error = function(e) NA_character_)
if (is.na(current_input_sha) || !identical(tolower(current_input_sha), tolower(frozen_input_sha))) {
  stop_preflight("INPUT_PATH SHA256 与 merged_obj_audit.csv 不一致。")
}
if (nrow(sample_audit) != EXPECTED_SAMPLES ||
    !all(c("SampleID_Global", "Status") %in% names(sample_audit)) ||
    any(sample_audit$Status != "PASS") ||
    anyDuplicated(as.character(sample_audit$SampleID_Global))) {
  stop_preflight("权威疾病审计不是 453 个唯一样本且全部 PASS。")
}
if (audit_value(stage_input, "input_path", "stage1 input audit") != INPUT_PATH ||
    !identical(tolower(audit_value(stage_input, "input_sha256", "stage1 input audit")),
               tolower(current_input_sha)) ||
    audit_value(stage_input, "n_cells", "stage1 input audit") != as.character(EXPECTED_CELLS) ||
    audit_value(stage_input, "n_genes", "stage1 input audit") != as.character(EXPECTED_GENES) ||
    audit_value(stage_input, "n_clusters_res0.6", "stage1 input audit") != as.character(EXPECTED_CLUSTERS) ||
    audit_value(stage_input, "n_samples", "stage1 input audit") != as.character(EXPECTED_SAMPLES) ||
    audit_value(stage_input, "metadata_disease_analyse_alias", "stage1 input audit") != "PRESENT_AND_IDENTICAL") {
  stop_preflight("stage1 input audit 未通过输入 SHA/规模/疾病别名核验。")
}

param_expected <- c(
  input_path = INPUT_PATH, assay = "RNA", data_layer = "data",
  active_ident = "seurat_clusters", n_cells = as.character(EXPECTED_CELLS),
  n_genes = as.character(EXPECTED_GENES), n_clusters = as.character(EXPECTED_CLUSTERS),
  test.use = "wilcox", only.pos = "TRUE", min.pct = "0.1",
  logfc.threshold = "0.25", densify = "FALSE",
  sampling = "none; all cells and all genes"
)
for (nm in names(param_expected)) {
  if (param_value(stage_param, nm) != param_expected[[nm]]) {
    stop_preflight("stage1 参数不符合冻结值：", nm)
  }
}
if (nrow(stage_archive) != length(OLD_NAMES) ||
    !all(stage_archive$status == "PASS") ||
    !identical(as.character(stage_archive$source_path), OLD_PATHS) ||
    !identical(as.character(stage_archive$archive_path), ARCHIVE_PATHS) ||
    any(!stage_archive$archive_exists) || any(!stage_archive$source_absent_after_move) ||
    anyNA(stage_archive$source_sha256) || anyNA(stage_archive$archive_sha256) ||
    !all(stage_archive$source_sha256 == stage_archive$archive_sha256)) {
  stop_preflight("stage1 archive audit 结构或状态不符合要求。")
}
old_hashes <- vapply(OLD_PATHS, sha256_file, character(1))
if (!all(old_hashes == as.character(stage_archive$source_sha256))) {
  stop_preflight("当前旧 res0.4 源文件 SHA256 与 stage1 archive audit 不一致。")
}
if (any(as.character(stage_archive$archive_sha256) != old_hashes)) {
  stop_preflight("stage1 archive audit 中 archive SHA256 与旧源文件不一致。")
}
if (file.exists(file.path(OUT_DIR, "cluster_stats_res0.3.csv")) == FALSE) {
  stop_preflight("保护性检查失败：cluster_stats_res0.3.csv 不在原位。")
}

# -----------------------------------------------------------------------------
# 2. 发布事务。此处开始所有可能失败的操作均走统一 failure_handler。
# -----------------------------------------------------------------------------
moved_indices <- integer(0)
published <- character(0)
interpretation_existed <- FALSE
interpretation_backup_ready <- FALSE
interpretation_backup <- file.path(STAGING_DIR, "interpretation_original_backup.tmp")
completion_stage <- file.path(STAGING_DIR, COMPLETION_NAME)
completion_pass_stage <- file.path(STAGING_DIR, "findallmarkers_res0.6_completion_PASS.tmp.csv")
completion_final <- file.path(OUT_DIR, COMPLETION_NAME)

clear_completion_files <- function() {
  paths <- c(completion_stage, completion_pass_stage)
  for (p in paths) if (file.exists(p)) file.remove(p)
  invisible(character(0))
}
restore_interpretation <- function() {
  if (!interpretation_backup_ready) return(character(0))
  if (interpretation_existed) {
    ok <- tryCatch(isTRUE(file.copy(interpretation_backup, INTERPRETATION_PATH, overwrite = TRUE)) &&
                     same_sha(interpretation_backup, INTERPRETATION_PATH),
                   error = function(e) FALSE)
    if (!ok) return(paste0("解读恢复失败：", INTERPRETATION_PATH))
  } else if (file.exists(INTERPRETATION_PATH)) {
    ok <- tryCatch(isTRUE(file.remove(INTERPRETATION_PATH)), error = function(e) FALSE)
    if (!ok) return(paste0("解读回滚失败：", INTERPRETATION_PATH))
  }
  character(0)
}
rollback_published <- function() {
  errors <- character(0)
  for (nm in rev(published)) {
    src <- file.path(OUT_DIR, nm)
    dst <- file.path(STAGING_DIR, nm)
    if (!file.exists(src)) next
    if (file.exists(dst)) {
      errors <- c(errors, paste0("新结果回滚目标已存在：", dst))
      next
    }
    ok <- tryCatch(isTRUE(file.rename(src, dst)), error = function(e) FALSE)
    if (!ok || !file.exists(dst) || file.exists(src)) {
      errors <- c(errors, paste0("新结果回滚失败：", src))
    }
  }
  if (length(errors) == 0L) published <<- character(0)
  errors
}
rollback_archive <- function() {
  errors <- character(0)
  for (j in rev(moved_indices)) {
    if (file.exists(OLD_PATHS[j]) && file.exists(ARCHIVE_PATHS[j])) {
      ah <- tryCatch(sha256_file(ARCHIVE_PATHS[j]), error = function(e) NA_character_)
      if (!is.na(ah) && identical(ah, old_hashes[j])) {
        removed <- tryCatch(isTRUE(file.remove(ARCHIVE_PATHS[j])), error = function(e) FALSE)
        if (!removed) errors <- c(errors, paste0("归档回滚清理失败：", ARCHIVE_PATHS[j]))
      } else errors <- c(errors, paste0("归档回滚同时存在且 SHA 不符：", ARCHIVE_PATHS[j]))
    } else if (!file.exists(OLD_PATHS[j]) && file.exists(ARCHIVE_PATHS[j])) {
      ok <- tryCatch(isTRUE(file.rename(ARCHIVE_PATHS[j], OLD_PATHS[j])), error = function(e) FALSE)
      if (!ok) errors <- c(errors, paste0("归档回滚失败：", ARCHIVE_PATHS[j]))
    } else if (!file.exists(OLD_PATHS[j])) {
      errors <- c(errors, paste0("归档回滚缺少原文件：", OLD_PATHS[j]))
    }
  }
  final_old <- vapply(OLD_PATHS, function(p) tryCatch(sha256_file(p), error = function(e) NA_character_), character(1))
  if (any(!file.exists(OLD_PATHS)) || any(file.exists(ARCHIVE_PATHS)) || !all(final_old == old_hashes)) {
    errors <- c(errors, "归档回滚后旧 res0.4 文件存在性或 SHA256 核验失败")
  }
  if (length(errors) == 0L) moved_indices <<- integer(0)
  unique(errors)
}
rollback_all <- function() {
  out <- character(0)
  out <- c(out, tryCatch(restore_interpretation(), error = function(e) paste0("解读回滚异常：", conditionMessage(e))))
  out <- c(out, tryCatch(rollback_published(), error = function(e) paste0("新结果回滚异常：", conditionMessage(e))))
  out <- c(out, tryCatch(rollback_archive(), error = function(e) paste0("归档回滚异常：", conditionMessage(e))))
  unique(out[nzchar(out)])
}
failure_handler <- function(reason) {
  rollback_errors <- tryCatch(rollback_all(), error = function(e) paste0("统一回滚异常：", conditionMessage(e)))
  clear_errors <- tryCatch({ clear_completion_files(); character(0) },
                           error = function(e) paste0("completion 清理异常：", conditionMessage(e)))
  errors <- unique(c(rollback_errors, clear_errors))
  stop(paste0(reason, "；统一回滚状态：", if (length(errors) == 0L) "PASS" else paste(errors, collapse = " | ")))
}
publish_one <- function(nm) {
  src <- file.path(STAGING_DIR, nm)
  dst <- file.path(OUT_DIR, nm)
  ok <- isTRUE(file.rename(src, dst))
  if (ok || (file.exists(dst) && !file.exists(src))) published <<- c(published, nm)
  ok && file.exists(dst) && !file.exists(src)
}

for (i in seq_along(OLD_PATHS)) {
  moved <- tryCatch(isTRUE(file.rename(OLD_PATHS[i], ARCHIVE_PATHS[i])), error = function(e) FALSE)
  if (moved || (file.exists(ARCHIVE_PATHS[i]) && !file.exists(OLD_PATHS[i]))) {
    moved_indices <- c(moved_indices, i)
  }
  if (!moved) failure_handler(paste0("旧 res0.4 移动失败：", OLD_PATHS[i]))
}
archive_check <- tryCatch(
  all(file.exists(ARCHIVE_PATHS), !file.exists(OLD_PATHS),
      vapply(ARCHIVE_PATHS, sha256_file, character(1)) == old_hashes),
  error = function(e) FALSE
)
if (!isTRUE(archive_check)) failure_handler("旧 res0.4 归档后存在性或 SHA256 核验失败")

for (nm in BASE_NAMES) {
  ok <- tryCatch(publish_one(nm), error = function(e) FALSE)
  if (!ok) failure_handler(paste0("stage1 新结果发布失败：", nm))
}
interpretation_existed <- file.exists(INTERPRETATION_PATH)
interpretation_backup_ok <- tryCatch({
  if (interpretation_existed) {
    file.copy(INTERPRETATION_PATH, interpretation_backup, overwrite = FALSE) &&
      same_sha(INTERPRETATION_PATH, interpretation_backup)
  } else TRUE
}, error = function(e) FALSE)
if (!interpretation_backup_ok) failure_handler("解读备份失败")
interpretation_backup_ready <- TRUE

section <- tryCatch(c(
  "", "## FindAllMarkers resolution 0.6（merged_obj；stage1 恢复发布）", "",
  "| 项目 | 内容 |", "|---|---|",
  "| 分析目的 | 识别 30 个 cluster 的阳性差异表达 marker，支持细胞类型和状态解释。 |",
  paste0("| 输入数据 | `", INPUT_PATH, "`；", EXPECTED_CELLS, " 个细胞、", EXPECTED_GENES, " 个 protein-coding genes、", EXPECTED_SAMPLES, " 个样本。 |"),
  "| 方法 | RNA assay 的 data layer；Wilcoxon；only.pos=TRUE；min.pct=0.1；logfc.threshold=0.25；densify=FALSE；全细胞、全基因、不抽样。 |",
  paste0("| 结果 | stage1 已验收全部 marker（", EXPECTED_MARKER_ROWS, " 行、", EXPECTED_MARKER_GENES, " 个 marker gene）、每簇 top10（", EXPECTED_TOP10_ROWS, " 行），覆盖 30 个 cluster。 |"),
  "| 生物学意义 | marker 集合用于 cluster 注释、状态解释及后续疾病亚型分析。 |",
  "| 限制 | 单细胞 Wilcoxon marker 主要用于探索和注释；正式疾病差异应在 donor/sample 层面建模。 |",
  paste0("| 旧结果归档 | 旧 res0.4 三个正式文件已完整移动至 `", ARCHIVE_DIR, "`；未删除；cluster_stats_res0.3.csv 未移动。 |"),
  ""
), error = function(e) NULL)
if (is.null(section) || !isTRUE(tryCatch({ cat(paste0(section, collapse = "\n"), "\n", file = INTERPRETATION_PATH, append = TRUE); TRUE }, error = function(e) FALSE))) {
  failure_handler("解读发布失败")
}

completion <- tryCatch(data.table(
  item = c("status", "completed_at", "input_path", "input_sha256", "n_cells", "n_genes", "n_clusters_res0.6", "marker_rows", "marker_genes", "marker_clusters", "top10_rows", "data_layer", "test.use", "only.pos", "min.pct", "logfc.threshold", "densify", "sampling", "archive_status", "archive_dir", "disease_analyse_alias_note", "object_modified"),
  value = c("READY", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), INPUT_PATH, current_input_sha, as.character(EXPECTED_CELLS), as.character(EXPECTED_GENES), as.character(EXPECTED_CLUSTERS), as.character(EXPECTED_MARKER_ROWS), as.character(EXPECTED_MARKER_GENES), as.character(EXPECTED_CLUSTERS), as.character(EXPECTED_TOP10_ROWS), "data", "wilcox", "TRUE", "0.1", "0.25", "FALSE", "none; all cells and all genes", "PASS", ARCHIVE_DIR, "present; identical to Disease_Abbrev_Analysis", "FALSE")
), error = function(e) NULL)
if (is.null(completion) || !isTRUE(tryCatch({ fwrite(completion, completion_stage); TRUE }, error = function(e) FALSE))) {
  failure_handler("completion READY 审计写入失败")
}
if (!isTRUE(tryCatch(publish_one(COMPLETION_NAME), error = function(e) FALSE))) {
  failure_handler("completion READY 发布失败")
}
if (!isTRUE(tryCatch(publish_one(LOG_NAME), error = function(e) FALSE))) {
  failure_handler("最终日志发布失败")
}
pass_ok <- isTRUE(tryCatch({
  completion[item == "status", value := "PASS"]
  fwrite(completion, completion_pass_stage)
  file.copy(completion_pass_stage, completion_final, overwrite = TRUE) &&
    identical(sha256_file(completion_pass_stage), sha256_file(completion_final))
}, error = function(e) FALSE))
if (!pass_ok) failure_handler("completion PASS 形成失败")
log_ok <- isTRUE(tryCatch({
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] === RESUME stage1 ALL DONE; transaction PASS ===\n"), file = file.path(OUT_DIR, LOG_NAME), append = TRUE)
  TRUE
}, error = function(e) FALSE))
if (!log_ok) failure_handler("最终日志收尾失败")
cleanup_ok <- isTRUE(tryCatch({
  if (!identical(basename(STAGING_DIR), ".staging_res0.6_20260831_174739")) stop("staging 路径不精确")
  unlink(STAGING_DIR, recursive = TRUE, force = TRUE)
  !dir.exists(STAGING_DIR)
}, error = function(e) FALSE))
if (!cleanup_ok) failure_handler("精确 stage1 staging 清理失败")
invisible(NULL)
