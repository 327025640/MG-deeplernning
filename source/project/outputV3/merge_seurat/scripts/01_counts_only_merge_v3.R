options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(Matrix)
})

root <- "D:/data/microglia_new"
out <- file.path(root, "outputV3/merge_seurat")
manifest_dir <- file.path(out, "00_manifest_and_preflight")
merge_dir <- file.path(out, "01_counts_merge")
obj_dir <- file.path(merge_dir, "objects")
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)

frozen_path <- file.path(manifest_dir, "v3_input_manifest_FROZEN.csv")
if (!file.exists(frozen_path)) stop("Frozen V3 manifest missing: ", frozen_path)
manifest <- fread(frozen_path)
if (!nrow(manifest)) stop("Frozen V3 manifest is empty")

# Step 02b 注入表（用户 2026-08-28 裁决：Merge_Include 门控、Treatment_Status 归一化、
# Disease_Abbrev_Analysis 均以本表为准，在 merge 时注入细胞 meta.data，不修改源对象）
injection_path <- file.path(manifest_dir, "sample_metadata_injection_v3.csv")
if (!file.exists(injection_path)) stop("Metadata injection table missing: ", injection_path)
inj <- fread(injection_path)
req_inj <- c("DatasetID", "SampleID", "SampleID_Global", "Merge_Include", "Treatment_Status", "Disease_Abbrev_Analysis")
if (!all(req_inj %in% names(inj))) stop("Injection table missing columns: ", paste(setdiff(req_inj, names(inj)), collapse = ", "))
if (any(!as.logical(inj$Merge_Include))) stop("Injection table contains Merge_Include != TRUE")
if (anyDuplicated(inj$SampleID_Global)) stop("Duplicate SampleID_Global in injection table")
conflict_log <- list()

get_counts <- function(obj) {
  if (!("RNA" %in% Assays(obj))) stop("RNA assay missing")
  ass <- obj[["RNA"]]
  lyr <- Layers(ass)
  cnt <- lyr[grepl("^counts", lyr)]
  if (!length(cnt)) stop("RNA assay has no counts layer")
  if (length(cnt) > 1L) ass <- JoinLayers(ass)
  lyr <- Layers(ass)
  nm <- if ("counts" %in% lyr) "counts" else lyr[grepl("^counts", lyr)][1L]
  # LayerData 兼容 Assay5 与旧版 Assay（GetAssayData slot= 在 SeuratObject 5 已弃用）
  LayerData(ass, layer = nm)
}

standardize_one <- function(i) {
  ds <- manifest$DatasetID[i]
  path <- manifest$Canonical_RDS_Path[i]
  obj <- readRDS(path)
  counts <- get_counts(obj)
  md <- obj[[]]
  rm(obj); invisible(gc())
  if (!identical(colnames(counts), rownames(md))) {
    if (all(colnames(counts) %in% rownames(md))) md <- md[colnames(counts), , drop = FALSE]
    else stop(ds, ": counts/metadata cell names do not match")
  }
  # ---- Step 02b 注入：按 DatasetID+SampleID 以注入表覆盖/添加 metadata ----
  smap <- inj[DatasetID == ds]
  if (!nrow(smap)) stop(ds, ": no rows in injection table")
  if (!all(as.character(md$SampleID) %in% smap$SampleID)) {
    stop(ds, ": samples missing from injection table: ",
         paste(setdiff(unique(as.character(md$SampleID)), smap$SampleID), collapse = ","))
  }
  mi <- match(as.character(md$SampleID), smap$SampleID)
  for (fld in c("Merge_Include", "SampleID_Global", "Treatment_Status", "Disease_Abbrev_Analysis")) {
    if (fld %in% names(md)) {
      old_val <- as.character(md[[fld]])
      new_val <- as.character(smap[[fld]])[mi]
      diff_cells <- !is.na(old_val) & old_val != new_val
      if (any(diff_cells)) {
        key <- paste(ds, as.character(md$SampleID[diff_cells]), fld, sep = "|")
        agg <- aggregate(rep(1L, sum(diff_cells)), by = list(SampleID = as.character(md$SampleID[diff_cells])), FUN = sum)
        conflict_log[[length(conflict_log) + 1L]] <<- data.table(
          DatasetID = ds, SampleID = agg$SampleID, Field = fld,
          Object_Value = "non-matching values present", Injected_Value = unique(new_val[diff_cells])[1],
          n_cells = agg$x)
      }
    }
    md[[fld]] <- as.character(smap[[fld]])[mi]
    if (fld == "Merge_Include") md[[fld]] <- as.logical(md[[fld]])
  }
  if ("Treatment_Status_Original" %in% names(smap)) {
    md$Treatment_Status_Original <- as.character(smap$Treatment_Status_Original)[mi]
  }
  if (!all(md$Merge_Include)) stop(ds, ": Merge_Include FALSE after injection")
  old <- colnames(counts)
  new <- paste(ds, old, sep = "__")
  if (anyDuplicated(new)) stop(ds, ": duplicate global cell IDs after prefixing")
  colnames(counts) <- new
  rownames(md) <- new
  md$CellID <- new
  md$DatasetID <- ds
  md$V3_Input_DatasetID <- ds
  # Rebuild a counts-only object. This strips stale data/scale.data/reductions/graphs
  # while preserving the complete per-cell metadata used for traceability.
  out_obj <- CreateSeuratObject(counts = counts, meta.data = md, project = ds)
  DefaultAssay(out_obj) <- "RNA"
  out_obj
}

# 内存敏感的渐进合并：每套标准化后立即释放源对象，合并后释放旧 atlas。
# 本机 32 GB RAM / D 盘 15 GB 空闲，故合并产物用 compress=TRUE 保存（稀疏矩阵压缩率高）。
atlas <- standardize_one(1L)
n_cells_by_ds <- ncol(atlas)
cat(sprintf("[%s] merged 1/%d %s (n=%d)\n", format(Sys.time(), "%H:%M:%S"), nrow(manifest), manifest$DatasetID[1], n_cells_by_ds[1]))
invisible(gc())
for (j in seq_len(nrow(manifest))[-1L]) {
  o <- standardize_one(j)
  n_cells_by_ds <- c(n_cells_by_ds, ncol(o))
  m <- merge(x = atlas, y = o, merge.data = FALSE)
  rm(atlas, o); invisible(gc())
  atlas <- m
  cat(sprintf("[%s] merged %d/%d %s (n=%d)\n", format(Sys.time(), "%H:%M:%S"), j, nrow(manifest), manifest$DatasetID[j], n_cells_by_ds[j]))
}
names(n_cells_by_ds) <- manifest$DatasetID
if (inherits(atlas[["RNA"]], "Assay5")) {
  atlas <- JoinLayers(atlas, assay = "RNA")
}
if (!("RNA" %in% Assays(atlas))) stop("Merged atlas has no RNA assay")
ass <- atlas[["RNA"]]
if (inherits(ass, "Assay5")) {
  lyr <- Layers(ass)
  if (!identical(lyr, "counts")) stop("Merged RNA layers are not counts-only: ", paste(lyr, collapse = ","))
} else {
  if (!all(c("counts") %in% slotNames(ass))) stop("Merged legacy assay has no counts slot")
}
if (ncol(atlas) != sum(n_cells_by_ds)) stop("Cell count is not conserved")
if (anyDuplicated(colnames(atlas))) stop("Global CellID collision after merge")
if (!all(c("CellID", "DatasetID", "SampleID", "DonorID") %in% names(atlas[[]]))) {
  stop("Required traceability metadata missing after merge")
}
md_chk <- as.data.table(atlas[[]])
if (any(!as.logical(md_chk$Merge_Include))) stop("Merged atlas contains Merge_Include FALSE")
bad_treat <- setdiff(unique(na.omit(as.character(md_chk$Treatment_Status))), c("No", "Not_reported", "Unknown"))
if (length(bad_treat)) stop("Treatment_Status values outside allowed set after injection: ", paste(bad_treat, collapse = ","))
if (any(is.na(md_chk$Disease_Abbrev_Analysis)) || any(!nzchar(md_chk$Disease_Abbrev_Analysis))) {
  stop("Disease_Abbrev_Analysis contains NA/empty after injection")
}

saveRDS(atlas, file.path(obj_dir, "microglia_atlas_counts_merged_v3.rds"), compress = TRUE)

# 注入冲突留痕（对象原值与注入表不一致的细胞数，按样本/字段汇总）
if (length(conflict_log)) fwrite(rbindlist(conflict_log, fill = TRUE), file.path(merge_dir, "injection_conflicts_log.csv"))

md_atlas <- as.data.table(atlas[[]])
summary <- md_atlas[, .(
  Cells = .N,
  Samples = uniqueN(as.character(SampleID)),
  Donors = uniqueN(as.character(DonorID))
), by = DatasetID]
setorder(summary, DatasetID)
fwrite(summary, file.path(merge_dir, "v3_merged_summary_by_dataset.csv"))
fwrite(data.table(
  Metric = c("Datasets", "Cells", "Genes", "RNA_layers", "CellID_duplicates"),
  Value = c(uniqueN(atlas$DatasetID), ncol(atlas), nrow(atlas),
            if (inherits(ass, "Assay5")) paste(Layers(ass), collapse = ";") else "counts_slot",
            sum(duplicated(colnames(atlas))))
), file.path(merge_dir, "v3_merged_validation.csv"))
writeLines(c(
  "V3 counts-only merge completed from the frozen canonical registry.",
  paste0("Datasets: ", uniqueN(atlas$DatasetID)),
  paste0("Cells: ", ncol(atlas)),
  paste0("Genes: ", nrow(atlas)),
  "All cell IDs were prefixed DatasetID__original_cell_id.",
  "Merge-time metadata injection applied from sample_metadata_injection_v3.csv;",
  "conflicts (object value vs injection) are logged in injection_conflicts_log.csv.",
  "Memory adaptation 2026-08-28: incremental merge with rm+gc; saved with compress=TRUE",
  "(local machine 32 GB RAM, D: 15 GB free). No source RDS or raw input was modified."
), file.path(merge_dir, "README.md"))
cat(sprintf("[%s] ALL DONE - atlas: %d cells x %d genes\n", format(Sys.time(), "%H:%M:%S"), ncol(atlas), nrow(atlas)))
