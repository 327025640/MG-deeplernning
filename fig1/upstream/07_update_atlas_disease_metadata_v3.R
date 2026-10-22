# 07_update_atlas_disease_metadata_v3.R
# 用户 2026-08-29 裁决：权威表外 20 样本按对象内疾病标签归入对应疾病（注入表已更新）。
# 直接在三个既有 atlas 对象的 meta.data 上按注入表覆盖 Disease_Abbrev_Analysis，不重跑 03-05。
# counts 与所有 reduction 不动；文件以 compress=TRUE 回存。

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(Seurat); library(data.table) })

root <- "D:/data/microglia_new/outputV3/merge_seurat"
inj <- fread(file.path(root, "00_manifest_and_preflight/sample_metadata_injection_v3.csv"))
smap <- inj[, .(SampleID_Global, Disease_Abbrev_Analysis, Disease_Source)]
if (anyDuplicated(smap$SampleID_Global)) stop("Duplicate SampleID_Global in injection table")

targets <- c(
  merged = file.path(root, "01_counts_merge/objects/microglia_atlas_counts_merged_v3.rds"),
  pca = file.path(root, "02_pca/objects/microglia_atlas_pca_v3.rds"),
  harmony = file.path(root, "03_harmony_preview/objects/microglia_atlas_harmony_preview_v3.rds")
)
for (nm in names(targets)) {
  p <- targets[[nm]]
  if (!file.exists(p)) stop("Missing atlas: ", p)
  cat(sprintf("[%s] %s loading...\n", format(Sys.time(), "%H:%M:%S"), nm))
  atlas <- readRDS(p)
  md <- as.data.table(atlas[[]])
  md[, CellID := rownames(atlas[[]])]
  if (!"SampleID_Global" %in% names(md)) stop(nm, ": SampleID_Global missing in meta.data")
  if (!all(md$SampleID_Global %in% smap$SampleID_Global)) {
    stop(nm, ": SampleID_Global values missing from injection table")
  }
  mi <- match(md$SampleID_Global, smap$SampleID_Global)
  old <- as.character(md$Disease_Abbrev_Analysis)
  new <- as.character(smap$Disease_Abbrev_Analysis)[mi]
  n_change <- sum(old != new, na.rm = TRUE) + sum(is.na(old))
  atlas$Disease_Abbrev_Analysis <- new
  atlas$Disease_Source <- as.character(smap$Disease_Source)[mi]
  cat(sprintf("[%s] %s: %d cells, %d labels changed\n",
              format(Sys.time(), "%H:%M:%S"), nm, ncol(atlas), n_change))
  saveRDS(atlas, p, compress = TRUE)
  rm(atlas, md); invisible(gc())
}
cat("ALL DONE\n")
