options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
})

# 可选命令行参数：registry 路径（缺省仍为 V2 registry，旧行为不变）
args <- commandArgs(trailingOnly = TRUE)

root <- "D:/data/microglia_new"
out <- file.path(root, "outputV3/merge_seurat")
manifest_dir <- file.path(out, "00_manifest_and_preflight")
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

registry_path <- if (length(args)) args[1] else file.path(root, "outputV2/02_atlas_merge/tables/canonical_input_registry.csv")
if (!file.exists(registry_path)) stop("Canonical input registry not found: ", registry_path)
reg <- fread(registry_path)
required_registry <- c("DatasetID", "Canonical_RDS_Path")
if (!all(required_registry %in% names(reg))) stop("Registry missing columns: ", paste(setdiff(required_registry, names(reg)), collapse = ", "))
if (anyDuplicated(reg$DatasetID)) stop("Duplicate DatasetID in canonical registry")
reg[, Canonical_RDS_Path := normalizePath(Canonical_RDS_Path, winslash = "/", mustWork = FALSE)]
reg[, Exists := file.exists(Canonical_RDS_Path)]
if (any(!reg$Exists)) stop("Missing canonical RDS: ", paste(reg[!Exists, DatasetID], collapse = ", "))

get_counts <- function(obj) {
  if (!("RNA" %in% Assays(obj))) stop("RNA assay missing")
  ass <- obj[["RNA"]]
  lyr <- Layers(ass)
  cnt <- lyr[grepl("^counts", lyr)]
  if (!length(cnt)) stop("RNA assay has no counts layer")
  if (length(cnt) > 1L) ass <- JoinLayers(ass)
  lyr <- Layers(ass)
  nm <- if ("counts" %in% lyr) "counts" else lyr[grepl("^counts", lyr)][1L]
  # LayerData 同时兼容 Assay5 与旧版 Assay（GetAssayData slot= 在 SeuratObject 5 已弃用）
  LayerData(ass, layer = nm)
}

one <- function(i) {
  ds <- reg$DatasetID[i]
  path <- reg$Canonical_RDS_Path[i]
  ans <- list(DatasetID = ds, Canonical_RDS_Path = path, Status = "PASS", Error = NA_character_)
  obj <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(obj, "error")) { ans$Status <- "FAIL"; ans$Error <- conditionMessage(obj); return(as.data.table(ans)) }
  # 2026-08-28 规则：整体 validObject 失败不阻断（例如重处理对象残留无效 reductions），
  # 只要 RNA assay 本身 valid 且 counts 可提取即放行；失败原因记录在 Whole_Object_ValidObject_Note。
  ans$Whole_Object_ValidObject_Note <- ""
  whole <- tryCatch({ validObject(obj); "" }, error = function(e) conditionMessage(e))
  if (nzchar(whole)) ans$Whole_Object_ValidObject_Note <- whole
  ok <- tryCatch(validObject(obj[["RNA"]]), error = function(e) e)
  if (inherits(ok, "error")) { ans$Status <- "FAIL"; ans$Error <- paste("RNA assay validObject:", conditionMessage(ok)); return(as.data.table(ans)) }
  md <- obj[[]]
  # 2026-08-28 用户裁决：对象内必需列放宽为 SampleID/DonorID/DatasetID；
  # Merge_Include 门控与 Treatment_Status 归一化移至 Step 03 merge 时按
  # sample_metadata_injection_v3.csv 注入，不在预检中 FAIL。
  req <- c("SampleID", "DonorID", "DatasetID")
  miss <- setdiff(req, names(md))
  if (length(miss)) { ans$Status <- "FAIL"; ans$Error <- paste("Missing metadata:", paste(miss, collapse = ", ")); return(as.data.table(ans)) }
  if (any(as.character(md$DatasetID) != ds, na.rm = TRUE)) { ans$Status <- "FAIL"; ans$Error <- "DatasetID metadata mismatch"; return(as.data.table(ans)) }
  counts <- tryCatch(get_counts(obj), error = function(e) e)
  if (inherits(counts, "error")) { ans$Status <- "FAIL"; ans$Error <- paste("counts:", conditionMessage(counts)); return(as.data.table(ans)) }
  if (!identical(colnames(counts), rownames(md))) {
    if (all(colnames(counts) %in% rownames(md))) md <- md[colnames(counts), , drop = FALSE]
    else { ans$Status <- "FAIL"; ans$Error <- "counts columns and metadata rows do not match"; return(as.data.table(ans)) }
  }
  if (anyDuplicated(colnames(counts))) { ans$Status <- "FAIL"; ans$Error <- "duplicate cell barcodes"; return(as.data.table(ans)) }
  ans$n_cells <- ncol(counts); ans$n_genes <- nrow(counts)
  ans$n_samples <- uniqueN(as.character(md$SampleID)); ans$n_donors <- uniqueN(as.character(md$DonorID))
  ans$Treatment_values <- if ("Treatment_Status" %in% names(md)) paste(sort(unique(na.omit(as.character(md$Treatment_Status)))), collapse = ";") else NA_character_
  ans$Perturbation_values <- if ("Experimental_Perturbation_Status" %in% names(md)) paste(sort(unique(na.omit(as.character(md$Experimental_Perturbation_Status)))), collapse = ";") else NA_character_
  as.data.table(ans)
}

audit <- rbindlist(lapply(seq_len(nrow(reg)), one), fill = TRUE)
fwrite(audit, file.path(manifest_dir, "v3_input_preflight_audit.csv"))
if (any(audit$Status != "PASS")) stop("V3 input preflight failed; see v3_input_preflight_audit.csv")

frozen <- merge(reg, audit[, c("DatasetID", "n_cells", "n_genes", "n_samples", "n_donors"), with = FALSE], by = "DatasetID", all.x = TRUE, sort = FALSE)
frozen[, c("Exists", "Source_Status") := list(TRUE, "PASS")]
fwrite(frozen, file.path(manifest_dir, "v3_input_manifest_FROZEN.csv"))
fwrite(frozen[, .(DatasetID, n_cells, n_genes, n_samples, n_donors)], file.path(manifest_dir, "v3_input_summary_by_dataset.csv"))

sha_cmd <- function(p) {
if (.Platform$OS.type == "windows") {
  out <- system2("certutil", c("-hashfile", p, "SHA256"), stdout = TRUE, stderr = TRUE)
  hex <- trimws(out)
  hex <- hex[grepl("^[0-9A-Fa-f]{64}$", hex)]
  if (!length(hex)) stop("SHA256 unavailable for ", p)
  return(hex[1L])
}
  trimws(system2("sha256sum", p, stdout = TRUE))
}
sha <- data.table(DatasetID = frozen$DatasetID, Path = frozen$Canonical_RDS_Path, SHA256 = vapply(frozen$Canonical_RDS_Path, sha_cmd, character(1)))
fwrite(sha, file.path(manifest_dir, "sha256_checksums.csv"))
writeLines(c(
  "V3 input manifest built from user-specified registry (default: outputV2 canonical registry).",
  "The registry path is authoritative; no ambiguous directory search was used.",
  paste0("Registry: ", registry_path),
  paste0("Datasets: ", nrow(frozen), "; preflight: PASS."),
  "2026-08-28 rule change (user-approved): object-required metadata relaxed to SampleID/DonorID/DatasetID;",
  "Merge_Include gating and Treatment_Status normalization are injected at Step 03 merge from",
  "sample_metadata_injection_v3.csv (per-sample authoritative lookup). Original RDS files are not modified.",
  "Whole-object validObject failures are recorded but non-blocking when the RNA assay itself is valid;",
  "stale invalid reductions in reprocessed objects are dropped by the counts-only rebuild at Step 03.",
  "Raw downloads and source RDS files were not modified."
), file.path(manifest_dir, "README.md"))
