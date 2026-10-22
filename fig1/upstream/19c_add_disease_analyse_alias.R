# =============================================================================
# 19c_add_disease_analyse_alias.R
# 目的：为已生成的 merged_obj.rds 补充用户冻结的 disease_analyse 别名列。
# 定义：disease_analyse 与权威细分列 Disease_Abbrev_Analysis 逐细胞完全一致。
# 除 meta.data 新增一列外，不改变细胞、基因、表达层、聚类、降维或图。
# =============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

OBJECT_PATH <- "E:/micro_output/v3_dims16/objects/merged_obj.rds"
AUDIT_PATH <- "E:/micro_output/v3_dims16/objects/merged_obj_audit.csv"
REPAIR_AUDIT_PATH <- "E:/micro_output/v3_dims16/objects/merged_obj_disease_analyse_alias_repair.csv"

required_paths <- c(OBJECT_PATH, AUDIT_PATH)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("缺少必需文件：", paste(missing_paths, collapse = "; "))
}
if (!requireNamespace("digest", quietly = TRUE)) stop("缺少 digest 包。")
if (file.exists(REPAIR_AUDIT_PATH)) {
  stop("修复审计已存在，为避免重复运行而停止：", REPAIR_AUDIT_PATH)
}

log_msg <- function(...) {
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
             paste0(..., collapse = ""), "\n"))
}

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE)
}

log_msg("读取并核验原对象审计。")
old_audit <- fread(AUDIT_PATH)
if (!all(c("item", "value") %in% names(old_audit))) {
  stop("merged_obj_audit.csv 缺少 item/value 列。")
}
audit_map <- setNames(as.character(old_audit$value), as.character(old_audit$item))
required_audit_items <- c(
  "output_path", "n_cells", "n_genes_protein_coding", "n_clusters_res0.6",
  "metadata_columns_n", "metadata_columns", "sha256"
)
missing_items <- setdiff(required_audit_items, names(audit_map))
if (length(missing_items) > 0L) {
  stop("对象审计缺少条目：", paste(missing_items, collapse = ", "))
}
if (!identical(unname(audit_map[["output_path"]]), OBJECT_PATH)) {
  stop("对象审计中的 output_path 与正式对象不一致。")
}
old_sha256 <- sha256_file(OBJECT_PATH)
if (!identical(old_sha256, unname(audit_map[["sha256"]]))) {
  stop("正式对象当前 SHA256 与冻结审计不一致。")
}

log_msg("读取正式对象：", OBJECT_PATH)
obj <- readRDS(OBJECT_PATH)
if (!inherits(obj, "Seurat") || !isTRUE(validObject(obj))) {
  stop("输入不是有效 Seurat 对象。")
}
if (!all(c("Disease", "Disease_Abbrev_Analysis") %in% colnames(obj[[]]))) {
  stop("对象缺少 Disease 或 Disease_Abbrev_Analysis。")
}
if ("disease_analyse" %in% colnames(obj[[]])) {
  stop("对象已包含 disease_analyse；本修复脚本不应重复运行。")
}

old_meta <- colnames(obj[[]])
audit_meta <- strsplit(
  unname(audit_map[["metadata_columns"]]), ";", fixed = TRUE
)[[1L]]
if (!identical(old_meta, audit_meta) ||
    !identical(length(old_meta), as.integer(audit_map[["metadata_columns_n"]]))) {
  stop("正式对象的实际 metadata 列与 merged_obj_audit.csv 不一致。")
}

before <- list(
  cells = colnames(obj),
  genes = rownames(obj),
  layers = Layers(obj[["RNA"]]),
  reductions = Reductions(obj),
  graphs = Graphs(obj),
  identities = as.character(Idents(obj)),
  clusters = as.character(obj$seurat_clusters)
)

obj$disease_analyse <- as.character(obj$Disease_Abbrev_Analysis)
meta_order <- colnames(obj[[]])
meta_order <- c(
  meta_order[seq_len(match("Disease", meta_order))],
  "disease_analyse",
  meta_order[seq.int(match("Disease", meta_order) + 1L, length(meta_order) - 1L)]
)
obj@meta.data <- obj@meta.data[, meta_order, drop = FALSE]

if (!identical(as.character(obj$disease_analyse),
               as.character(obj$Disease_Abbrev_Analysis)) ||
    anyNA(obj$disease_analyse)) {
  stop("disease_analyse 别名建立失败。")
}
if (!identical(colnames(obj), before$cells) ||
    !identical(rownames(obj), before$genes) ||
    !identical(Layers(obj[["RNA"]]), before$layers) ||
    !identical(Reductions(obj), before$reductions) ||
    !identical(Graphs(obj), before$graphs) ||
    !identical(as.character(Idents(obj)), before$identities) ||
    !identical(as.character(obj$seurat_clusters), before$clusters) ||
    !isTRUE(validObject(obj))) {
  stop("新增别名后，对象的细胞/基因/层/聚类或结构发生意外变化。")
}

object_dir <- dirname(OBJECT_PATH)
tmp_object <- tempfile("merged_obj_alias_", tmpdir = object_dir, fileext = ".tmp.rds")
tmp_audit <- tempfile("merged_obj_audit_alias_", tmpdir = object_dir, fileext = ".tmp.csv")
backup_object <- file.path(object_dir, "merged_obj_before_disease_analyse_alias.tmp.rds")
backup_audit <- file.path(object_dir, "merged_obj_audit_before_disease_analyse_alias.tmp.csv")
temp_paths <- c(tmp_object, tmp_audit, backup_object, backup_audit)
failed_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
failed_object <- file.path(object_dir, paste0("merged_obj_failed_alias_", failed_stamp, ".rds"))
failed_audit <- file.path(object_dir, paste0("merged_obj_audit_failed_alias_", failed_stamp, ".csv"))
on.exit({
  # 只清理尚未发布的临时新文件；失败时保留备份供人工恢复。
  for (path_i in c(tmp_object, tmp_audit)[file.exists(c(tmp_object, tmp_audit))]) {
    unlink(path_i)
  }
}, add = TRUE)
if (any(file.exists(c(backup_object, backup_audit)))) {
  stop("发现旧的修复备份文件，请先人工检查：",
       paste(c(backup_object, backup_audit)[file.exists(c(backup_object, backup_audit))],
             collapse = "; "))
}

restore_original_or_stop <- function(reason) {
  move_new_object_ok <- TRUE
  move_new_audit_ok <- TRUE
  if (file.exists(OBJECT_PATH)) {
    move_new_object_ok <- isTRUE(file.rename(OBJECT_PATH, failed_object))
  }
  if (file.exists(AUDIT_PATH)) {
    move_new_audit_ok <- isTRUE(file.rename(AUDIT_PATH, failed_audit))
  }

  restore_object_ok <- file.exists(backup_object) &&
    isTRUE(file.rename(backup_object, OBJECT_PATH)) && file.exists(OBJECT_PATH)
  restore_audit_ok <- file.exists(backup_audit) &&
    isTRUE(file.rename(backup_audit, AUDIT_PATH)) && file.exists(AUDIT_PATH)

  restored_integrity_ok <- FALSE
  if (restore_object_ok && restore_audit_ok) {
    restored_audit <- tryCatch(fread(AUDIT_PATH), error = function(e) NULL)
    restored_sha <- tryCatch(sha256_file(OBJECT_PATH), error = function(e) NA_character_)
    if (!is.null(restored_audit) && all(c("item", "value") %in% names(restored_audit))) {
      restored_map <- setNames(
        as.character(restored_audit$value), as.character(restored_audit$item)
      )
      restored_integrity_ok <- identical(restored_sha, old_sha256) &&
        identical(unname(restored_map[["sha256"]]), old_sha256)
    }
  }

  if (restore_object_ok && restore_audit_ok && restored_integrity_ok) {
    if (file.exists(failed_object)) unlink(failed_object)
    if (file.exists(failed_audit)) unlink(failed_audit)
    stop(reason, "；原对象和原审计已恢复并通过旧SHA256核验。")
  }

  stop(
    reason,
    "；自动恢复未完整成功。请勿删除以下仍存在的文件并人工恢复：",
    paste(c(backup_object, backup_audit, failed_object, failed_audit)[
      file.exists(c(backup_object, backup_audit, failed_object, failed_audit))
    ], collapse = "; "),
    "；move_new_object_ok=", move_new_object_ok,
    "；move_new_audit_ok=", move_new_audit_ok,
    "；restore_object_ok=", restore_object_ok,
    "；restore_audit_ok=", restore_audit_ok,
    "；restored_integrity_ok=", restored_integrity_ok
  )
}

log_msg("保存并重读临时对象。")
saveRDS(obj, tmp_object, compress = "gzip")
rm(obj)
invisible(gc())

check_obj <- readRDS(tmp_object)
if (!inherits(check_obj, "Seurat") || !isTRUE(validObject(check_obj)) ||
    !identical(colnames(check_obj), before$cells) ||
    !identical(rownames(check_obj), before$genes) ||
    !identical(Layers(check_obj[["RNA"]]), before$layers) ||
    !identical(Reductions(check_obj), before$reductions) ||
    !identical(Graphs(check_obj), before$graphs) ||
    !identical(as.character(check_obj$disease_analyse),
               as.character(check_obj$Disease_Abbrev_Analysis)) ||
    !identical(as.character(Idents(check_obj)), before$identities) ||
    !identical(as.character(check_obj$seurat_clusters), before$clusters)) {
  stop("临时对象重读核验失败。")
}
rm(check_obj)
invisible(gc())

new_sha256 <- sha256_file(tmp_object)
new_size <- file.info(tmp_object)$size
new_meta <- meta_order

new_audit <- copy(old_audit)
new_audit[item == "metadata_columns_n", value := as.character(length(new_meta))]
new_audit[item == "metadata_columns", value := paste(new_meta, collapse = ";")]
new_audit[item == "output_size_bytes", value := as.character(new_size)]
new_audit[item == "sha256", value := new_sha256]
new_audit <- rbind(
  new_audit,
  data.table(
    item = c("disease_analyse_definition", "disease_analyse_alias_added_at"),
    value = c(
      "exact copy of Disease_Abbrev_Analysis; subtype-level analysis column",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
  )
)
fwrite(new_audit, tmp_audit)

# 同目录事务替换：任一步失败均恢复原对象和原审计。
if (!file.rename(OBJECT_PATH, backup_object)) stop("无法备份原对象。")
if (!file.rename(AUDIT_PATH, backup_audit)) {
  restore_back_ok <- isTRUE(file.rename(backup_object, OBJECT_PATH)) &&
    file.exists(OBJECT_PATH) && identical(sha256_file(OBJECT_PATH), old_sha256)
  if (!restore_back_ok) {
    stop("无法备份原审计，且原对象自动恢复失败；备份保留于：", backup_object)
  }
  stop("无法备份原对象审计；原对象已恢复并通过旧SHA256核验。")
}
if (!file.rename(tmp_object, OBJECT_PATH)) {
  restore_original_or_stop("无法发布新增别名后的对象")
}
if (!file.rename(tmp_audit, AUDIT_PATH)) {
  restore_original_or_stop("无法发布新审计")
}

final_sha256 <- sha256_file(OBJECT_PATH)
final_audit <- fread(AUDIT_PATH)
final_map <- setNames(as.character(final_audit$value), as.character(final_audit$item))
if (!identical(final_sha256, new_sha256) ||
    !identical(unname(final_map[["sha256"]]), new_sha256) ||
    !identical(as.integer(final_map[["metadata_columns_n"]]), length(new_meta)) ||
    !"disease_analyse" %in% strsplit(
      unname(final_map[["metadata_columns"]]), ";", fixed = TRUE
    )[[1L]]) {
  restore_original_or_stop("最终对象或审计核验失败")
}

repair_audit <- data.table(
  item = c(
    "status", "object_path", "old_sha256", "new_sha256", "n_cells",
    "n_genes", "metadata_columns_before", "metadata_columns_after",
    "alias_definition", "expression_or_cluster_changed"
  ),
  value = c(
    "PASS", OBJECT_PATH, old_sha256, new_sha256, length(before$cells),
    length(before$genes), length(new_meta) - 1L, length(new_meta),
    "disease_analyse = Disease_Abbrev_Analysis", "FALSE"
  )
)
fwrite(repair_audit, REPAIR_AUDIT_PATH)

unlink(backup_object)
unlink(backup_audit)
log_msg("完成：已新增 disease_analyse，细胞/基因/聚类/表达不变；新SHA256=", new_sha256)
