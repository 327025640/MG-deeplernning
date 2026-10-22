# 08_multires_scan_v3.R
# Step 06 证据包：在 harmony_v3（50 维）与 harmony_dims16（16 维）两个 embedding 上
# 扫描 Louvain resolution = 0.2/0.4/0.6/0.8/1.0/1.2，输出簇数、簇大小分布与逐分辨率标签表。
# 用法: Rscript 08_multires_scan_v3.R   （顺序处理两套，输出到各自 benchmark 目录）

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(Seurat); library(data.table) })

tasks <- list(
  list(atlas = "D:/data/microglia_new/outputV3/merge_seurat/03_harmony_preview/objects/microglia_atlas_harmony_preview_v3.rds",
       outdir = "D:/data/microglia_new/outputV3/merge_seurat/06_multires_scan/harmony_dims50",
       red = "harmony_v3", dims = 50L, tag = "dims50"),
  list(atlas = "E:/micro_output/v3_dims16/objects/microglia_atlas_harmony_dims16_v3.rds",
       outdir = "E:/micro_output/v3_dims16/06_multires_scan/harmony_dims16",
       red = "harmony_dims16", dims = 16L, tag = "dims16")
)
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)

for (tk in tasks) {
  cat(sprintf("[%s] === %s ===\n", format(Sys.time(), "%H:%M:%S"), tk$tag))
  dir.create(tk$outdir, recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(tk$atlas)) stop("Missing atlas: ", tk$atlas)
  atlas <- readRDS(tk$atlas)
  md <- as.data.table(atlas[[]])
  md[, CellID := rownames(atlas[[]])]
  cat(sprintf("[%s] FindNeighbors on %s 1:%d...\n", format(Sys.time(), "%H:%M:%S"), tk$red, tk$dims))
  atlas <- FindNeighbors(atlas, reduction = tk$red, dims = seq_len(tk$dims), verbose = FALSE)
  summary_rows <- list()
  for (r in resolutions) {
    atlas <- FindClusters(atlas, resolution = r, random.seed = 42, verbose = FALSE)
    cl <- as.character(atlas$seurat_clusters)
    sz <- table(cl)
    lab_file <- file.path(tk$outdir, paste0("cluster_labels_res", r, ".csv"))
    if (!file.exists(lab_file)) fwrite(data.table(CellID = md$CellID, cluster = cl), lab_file)
    summary_rows[[length(summary_rows) + 1L]] <- data.table(
      Tag = tk$tag, Reduction = tk$red, Dims = tk$dims, Resolution = r,
      N_clusters = length(sz),
      Min_cells = as.integer(min(sz)), Median_cells = as.integer(median(sz)),
      Max_cells = as.integer(max(sz)),
      N_clusters_lt100 = sum(sz < 100), N_clusters_lt20 = sum(sz < 20))
    cat(sprintf("[%s] res=%.1f -> %d clusters\n", format(Sys.time(), "%H:%M:%S"), r, length(sz)))
  }
  fwrite(rbindlist(summary_rows), file.path(tk$outdir, paste0("multires_summary_", tk$tag, ".csv")), append = file.exists(file.path(tk$outdir, paste0("multires_summary_", tk$tag, ".csv"))))
  rm(atlas, md); invisible(gc())
}
cat("ALL DONE\n")
