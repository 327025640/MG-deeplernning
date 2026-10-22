# 03d_harmony_dims16_preview_v3.R
# dims=16 对比运行：Harmony 用 PCA 1:16（官网矩阵 API，DatasetID 批次），前后 UMAP 同参对比。
# 输出全部写入 E:/micro_output/v3_dims16/（D 盘空间不足，用户 2026-08-29 指定）。
# 预览图规格与 03b 一致（点径 0.6、图例放大、映射表配色、ASCII 标题）。

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(data.table)
  library(ggplot2)
})

root <- "D:/data/microglia_new/outputV3/merge_seurat"
infile <- file.path(root, "02_pca/objects/microglia_atlas_pca_v3.rds")
outroot <- "E:/micro_output/v3_dims16"
outdir_pre <- file.path(outroot, "pre_harmony16")
outdir_post <- file.path(outroot, "post_harmony16")
srcdir <- file.path(outroot, "source_data")
dir.create(file.path(outroot, "objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_pre, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_post, recursive = TRUE, showWarnings = FALSE)
dir.create(srcdir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(infile))

dims_use <- 16L
cat(sprintf("[%s] loading pca atlas...\n", format(Sys.time(), "%H:%M:%S")))
atlas <- readRDS(infile)
stopifnot(inherits(atlas, "Seurat"), "pca" %in% Reductions(atlas))
md <- atlas[[]]
stopifnot("DatasetID" %in% names(md))
if (!"Disease_Abbrev_Analysis" %in% names(md)) md$Disease_Abbrev_Analysis <- "Unspecified"
if (!"Brain_Region_Group" %in% names(md)) md$Brain_Region_Group <- "Not_reported"

pc <- Embeddings(atlas, "pca")
V <- pc[, seq_len(dims_use), drop = FALSE]
cat(sprintf("[%s] Harmony on dims 1:%d (official matrix API)...\n", format(Sys.time(), "%H:%M:%S"), dims_use))
h16 <- harmony::RunHarmony(V, md, "DatasetID", verbose = FALSE)
if (is.list(h16) && !is.null(h16$Embeddings)) h16 <- h16$Embeddings
h16 <- as.matrix(h16)
stopifnot(nrow(h16) == ncol(atlas))
atlas[["harmony_dims16"]] <- CreateDimReducObject(embeddings = h16, key = "HARM16_", assay = DefaultAssay(atlas))

cat(sprintf("[%s] UMAP pre (pca 1:%d)...\n", format(Sys.time(), "%H:%M:%S"), dims_use))
atlas <- RunUMAP(atlas, reduction = "pca", dims = seq_len(dims_use),
                 n.neighbors = 30, min.dist = 0.3, metric = "cosine",
                 seed.use = 42, reduction.name = "umap_pca16", reduction.key = "UMAPPC16_", verbose = FALSE)
cat(sprintf("[%s] UMAP post (harmony16 1:%d)...\n", format(Sys.time(), "%H:%M:%S"), dims_use))
atlas <- RunUMAP(atlas, reduction = "harmony_dims16", dims = seq_len(ncol(h16)),
                 n.neighbors = 30, min.dist = 0.3, metric = "cosine",
                 seed.use = 42, reduction.name = "umap_harmony16", reduction.key = "UMAPH16_", verbose = FALSE)

cat(sprintf("[%s] saving atlas...\n", format(Sys.time(), "%H:%M:%S")))
saveRDS(atlas, file.path(outroot, "objects/microglia_atlas_harmony_dims16_v3.rds"), compress = TRUE)

# ---------- 预览图（规格同 03b）----------
palette_file <- "D:/data/microglia_new/authority/03_reference_configs/disease_palette_contract.csv"
palette_tbl <- if (file.exists(palette_file)) fread(palette_file) else data.table()
dataset_palette_file <- "D:/data/microglia_new/authority/03_reference_configs/dataset_palette_contract.csv"
dataset_palette_tbl <- if (file.exists(dataset_palette_file)) fread(dataset_palette_file) else data.table()
make_palette <- function(x, variable) {
  u <- sort(unique(as.character(x[!is.na(x)])))
  if (!length(u)) return(character())
  if (variable == "Disease_Abbrev_Analysis" && nrow(palette_tbl) > 0L &&
      all(c("Disease_Abbrev", "Disease_Color_Hex") %in% names(palette_tbl))) {
    z <- setNames(as.character(palette_tbl$Disease_Color_Hex), toupper(as.character(palette_tbl$Disease_Abbrev)))
    p <- z[toupper(u)]
    miss <- is.na(p) | !nzchar(p)
    if (any(miss)) p[miss] <- hcl.colors(sum(miss), "Dark 3")
    names(p) <- u
    return(p)
  }
  if (variable == "DatasetID" && nrow(dataset_palette_tbl) > 0L &&
      all(c("DatasetID", "Color_Hex") %in% names(dataset_palette_tbl))) {
    z <- setNames(as.character(dataset_palette_tbl$Color_Hex), as.character(dataset_palette_tbl$DatasetID))
    p <- z[u]
    miss <- is.na(p) | !nzchar(p)
    if (any(miss)) p[miss] <- hcl.colors(sum(miss), "Dark 3")
    names(p) <- u
    return(p)
  }
  setNames(hcl.colors(length(u), "Dark 3"), u)
}

variables <- c("DatasetID", "SampleID", "DonorID", "Disease_Abbrev_Analysis",
               "Brain_Region_Group", "Treatment_Status")
set.seed(42)
keep <- sample(seq_len(ncol(atlas)), min(100000L, ncol(atlas)))
plot_meta <- atlas[[]][keep, , drop = FALSE]
plot_meta$CellID <- rownames(plot_meta)
pre_xy <- Embeddings(atlas, "umap_pca16")[keep, 1:2, drop = FALSE]
post_xy <- Embeddings(atlas, "umap_harmony16")[keep, 1:2, drop = FALSE]
xlim <- range(c(pre_xy[, 1], post_xy[, 1]), finite = TRUE)
ylim <- range(c(pre_xy[, 2], post_xy[, 2]), finite = TRUE)

plot_one <- function(xy, stage, variable) {
  vals <- if (variable %in% names(plot_meta)) as.character(plot_meta[[variable]]) else rep("Not_reported", nrow(plot_meta))
  nonmiss <- vals[!is.na(vals) & nzchar(vals)]
  if (variable %in% c("Treatment_Status", "Experimental_Perturbation_Status") && length(unique(nonmiss)) <= 1L) {
    return(data.table(variable = variable, stage = stage, status = "skipped_single_category"))
  }
  vals[is.na(vals) | !nzchar(vals)] <- "Not_reported"
  d <- data.table(CellID = plot_meta$CellID, UMAP_1 = xy[, 1], UMAP_2 = xy[, 2], Value = vals)
  fwrite(d, file.path(srcdir, paste0(stage, "_", variable, ".csv")))
  stage_lab <- if (stage == "pre_harmony16") "Origin (PCA 1-16)" else "Harmony (dims 1-16)"
  p <- ggplot(d, aes(UMAP_1, UMAP_2, colour = Value)) +
    geom_point(size = 0.6, alpha = 0.5, stroke = 0) +
    coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_classic(base_family = "Arial") +
    labs(title = paste0(variable, " - ", stage_lab), colour = variable, x = "UMAP 1", y = "UMAP 2")
  pal <- make_palette(vals, variable)
  if (variable %in% c("SampleID", "DonorID")) {
    p <- p + guides(colour = "none")
  } else {
    p <- p + scale_colour_manual(values = pal, na.value = "grey70") +
      guides(colour = guide_legend(
        override.aes = list(size = 3.5, alpha = 1, stroke = 0.4),
        title.theme = element_text(size = 11, face = "bold"),
        text.theme = element_text(size = 9)))
  }
  odir <- if (stage == "pre_harmony16") outdir_pre else outdir_post
  ggsave(file.path(odir, paste0(variable, ".png")), p, width = 8, height = 6, dpi = 600, bg = "white")
  ggsave(file.path(odir, paste0(variable, ".pdf")), p, width = 8, height = 6, device = cairo_pdf, bg = "white")
  data.table(variable = variable, stage = stage, status = "written", n_cells = nrow(d))
}

audit <- rbindlist(lapply(variables, function(v) plot_one(pre_xy, "pre_harmony16", v)), fill = TRUE)
audit <- rbindlist(list(audit, rbindlist(lapply(variables, function(v) plot_one(post_xy, "post_harmony16", v)), fill = TRUE)), fill = TRUE)
fwrite(audit, file.path(srcdir, "dims16_plot_audit.csv"))
writeLines(c(
  "# V3 Harmony dims=16 comparison",
  paste0("Input: ", infile),
  "Harmony input: PCA dims 1-16; batch variable: DatasetID; official matrix API; seed 42.",
  "UMAPs: pca 1-16 vs harmony_dims16 1-16, identical params (nn=30, min.dist=0.3, cosine, seed 42).",
  "Output staged on E: (D: disk full). Final dims/resolution remain unfrozen."
), file.path(outroot, "README.md"))
cat(sprintf("[%s] ALL DONE: %d figures\n", format(Sys.time(), "%H:%M:%S"), sum(audit$status == "written")))
