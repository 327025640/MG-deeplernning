options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
})

root <- "D:/data/microglia_new/outputV3/merge_seurat"
infile <- file.path(root, "01_counts_merge/objects/microglia_atlas_counts_merged_v3.rds")
outdir <- file.path(root, "02_pca")
dir.create(file.path(outdir, "objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "source_data"), recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(infile))

atlas <- readRDS(infile)
stopifnot(inherits(atlas, "Seurat"), "RNA" %in% Assays(atlas))

atlas <- NormalizeData(
  atlas,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)
atlas <- FindVariableFeatures(
  atlas,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)
hvg <- VariableFeatures(atlas)
stopifnot(length(hvg) == 3000L)
# 内存适配（2026-08-28，用户要求保留 counts）：counts 层全程保留；
# 仅在 PCA 后丢弃 scale.data（= data 层 + HVG 均值/方差的确定性换算，可随时无损重算）。
# 实测 counts 内存约 5-6 GB，ScaleData 峰值约 22 GB，本机 32 GB 可容纳。
invisible(gc())
cat(sprintf("[%s] ScaleData on %d HVG...\n", format(Sys.time(), "%H:%M:%S"), length(hvg)))
atlas <- ScaleData(atlas, features = hvg, verbose = FALSE)
cat(sprintf("[%s] RunPCA...\n", format(Sys.time(), "%H:%M:%S")))
atlas <- RunPCA(atlas, features = hvg, npcs = 50, verbose = FALSE)
# 仅移除 scale_data 层（Seurat5 DietSeurat 的 scale.data 参数弃用且会误删 reduction）；
# scale.data = data 层 + HVG 均值/方差的确定性换算，可随时无损重算。counts/data/pca 全部保留。
ass <- atlas[["RNA"]]
for (L in intersect(c("scale_data", "scale.data"), Layers(ass))) ass@layers[[L]] <- NULL
invisible(gc())

saveRDS(
  atlas,
  file.path(outdir, "objects/microglia_atlas_pca_v3.rds"),
  compress = TRUE
)

fwrite(
  data.table(Gene = hvg, Rank = seq_along(hvg)),
  file.path(outdir, "variable_features_3000.csv")
)
stdev <- Stdev(atlas[["pca"]])
pca_tbl <- data.table(
  PC = seq_along(stdev),
  StandardDeviation = stdev,
  Variance = stdev^2,
  VariancePercent = 100 * stdev^2 / sum(stdev^2)
)
fwrite(pca_tbl, file.path(outdir, "pca_variance_50pcs.csv"))

p <- ggplot(pca_tbl, aes(PC, StandardDeviation)) +
  geom_point(size = 1.4) +
  geom_line() +
  theme_classic(base_family = "Arial") +
  labs(title = "PCA elbow - V3", x = "PC", y = "Standard deviation")
ggsave(file.path(outdir, "figures/pca_elbow_v3.png"), p, width = 7, height = 5, dpi = 600, bg = "white")
ggsave(file.path(outdir, "figures/pca_elbow_v3.pdf"), p, width = 7, height = 5, device = cairo_pdf, bg = "white")
fwrite(pca_tbl, file.path(root, "source_data/pca_elbow_v3_source.csv"))

writeLines(
  c(
    "# V3 PCA stage",
    "",
    paste0("- Input: ", infile),
    "- LogNormalize, 3,000 HVG (vst), ScaleData on HVG only, PCA 50 PCs.",
    "- Memory adaptation 2026-08-28 (user-approved): counts and data layers retained;",
    "  scale.data dropped after PCA (deterministic recompute from data+HVG, lossless).",
    "  Saved with compress=TRUE (local machine 32 GB RAM, D: 15 GB free).",
    "- Final dims and resolution are intentionally not fixed."
  ),
  file.path(outdir, "README.md")
)
cat(sprintf("[%s] ALL DONE\n", format(Sys.time(), "%H:%M:%S")))
