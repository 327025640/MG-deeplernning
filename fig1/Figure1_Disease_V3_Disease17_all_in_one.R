# =============================================================================
# Figure1_Disease_V3_Disease17_all_in_one.R
#
# 目的：基于当前 merged_obj.rds，从头生成 Figure 1 的 Panel A-F。
# DEG 契约：使用 Seurat::AggregateExpression() 按 SampleID_Global 聚合原始 RNA
#            counts；随后使用 edgeR quasi-likelihood 官方流程检验每个 Disease vs Ctrl。
#            统计单位为 sample；不使用 donor，不加入 Dataset、年龄、性别或其他协变量。
#            DEG 定义为 edgeR FDR < 0.05 且 |log2FC| >= 2。
#
# 输出：每个 panel 独立 PNG + PDF + source CSV；另存 DEG、GO、QA、sessionInfo 和解读。
# 本地 R：D:/R/R-4.4.3/bin/Rscript.exe
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(data.table)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggrepel)
  library(msigdbr)
  library(clusterProfiler)
})

# ----------------------------- 0. 路径与常量 -------------------------------
ATLAS <- "E:/micro_output/v3_dims16/objects/merged_obj.rds"
INPUT_AUDIT <- "E:/micro_output/v3_dims16/objects/merged_obj_audit.csv"
ROOT <- "E:/micro_output/Figure1_Disease_V3_Disease17"
FIG <- file.path(ROOT, "figures")
SRC <- file.path(ROOT, "source_data")
TBL <- file.path(ROOT, "tables")
AUD <- file.path(ROOT, "audit")
LOG <- file.path(ROOT, "logs")
SCRIPT_OUT <- file.path(ROOT, "script")
for (d in c(ROOT, FIG, SRC, TBL, AUD, LOG, SCRIPT_OUT)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

PALETTE_DISEASE <- "D:/data/microglia_new/authority/03_reference_configs/disease_palette_contract.csv"
PALETTE_DATASET <- "D:/data/microglia_new/authority/03_reference_configs/dataset_palette_contract.csv"
MIN_CELLS_PER_SAMPLE <- 20L
TOP_K_MAIN <- 3L
TOP_K_TABLE <- 5L
RANDOM_SEED <- 20260901L

DATASET_GSE <- c(
  "ABDELFATTAH_GBM_2022" = "GSE182109",
  "BOSTRAND_HD_GLIA_2024" = "GSE281069",
  "GSE103224" = "GSE103224",
  "GSE141383" = "GSE141383",
  "KISS_2023" = "GSE227781",
  "LIM_2022" = "GSE180928",
  "LIU_SEOAD_2025" = "GSE272082",
  "LUNG_BRAIN_MET_2022" = "GSE202371",
  "MARTIROSYAN_2024" = "GSE243639",
  "MEDULLOBLASTOMA_2021" = "GSE155446",
  "SMITH_2022" = "GSE160936",
  "TAGORE_2025_NSCLC_BM" = "GSE223499",
  "TUDDENHAM_2024" = "GSE204702",
  "WANG_C9ORF72_ALS_FTD_2024" = "GSE219280",
  "WANG_C9_DLPFC_PNAS_2025" = "GSE212630",
  "ZHU_2024" = "GSE202210"
)

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
}

save_gg <- function(plot, stem, width, height, dpi = 600) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
         dpi = dpi, bg = "white")
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
         device = cairo_pdf, bg = "white")
}

required_packages <- c(
  "Seurat", "edgeR", "data.table", "ggplot2", "ComplexHeatmap", "circlize",
  "ggrepel", "msigdbr", "clusterProfiler", "digest"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  stop("缺少 R 包：", paste(missing_packages, collapse = ", "))
}
for (p in c(ATLAS, INPUT_AUDIT, PALETTE_DISEASE, PALETTE_DATASET)) {
  if (!file.exists(p)) stop("输入文件不存在：", p)
}

# -------------------------- 1. 对象与 metadata 审计 -------------------------
log_msg("读取 Seurat 对象：", ATLAS)
obj <- readRDS(ATLAS)
if (!inherits(obj, "Seurat")) stop("ATLAS 不是 Seurat 对象")
if (!identical(Assays(obj), "RNA")) stop("对象必须只包含 RNA assay")
if (!all(c("counts", "data") %in% Layers(obj[["RNA"]]))) {
  stop("RNA assay 必须包含 counts 和 data layer")
}
if (!"umap_harmony16" %in% Reductions(obj)) stop("缺少 umap_harmony16 reduction")

md <- as.data.table(obj[[]])
md[, CellID := colnames(obj)]
required_meta <- c(
  "CellID", "DatasetID", "SampleID", "SampleID_Global",
  "Disease", "nCount_RNA", "nFeature_RNA", "seurat_clusters"
)
missing_meta <- setdiff(required_meta, names(md))
if (length(missing_meta)) stop("metadata 缺少列：", paste(missing_meta, collapse = ", "))
if (anyDuplicated(md$CellID)) stop("CellID 不唯一")
if (anyNA(md$Disease) || any(!nzchar(as.character(md$Disease)))) stop("Disease 存在缺失值")

input_audit <- fread(INPUT_AUDIT)
frozen_sha256 <- input_audit[item == "sha256", value][1]
actual_sha256 <- digest::digest(ATLAS, algo = "sha256", file = TRUE)
if (is.na(frozen_sha256) || tolower(actual_sha256) != tolower(frozen_sha256)) {
  stop("merged_obj.rds 的 SHA256 与审计文件不一致")
}

dpal <- fread(PALETTE_DISEASE)
spal <- fread(PALETTE_DATASET)
needed_dpal <- c("Disease_Abbrev", "Disease_Color_Hex", "Disease_Class", "Disease_Class_Color_Hex")
if (!all(needed_dpal %in% names(dpal))) stop("疾病颜色表列不完整")
if (!all(c("DatasetID", "Color_Hex") %in% names(spal))) stop("数据集颜色表列不完整")

# 在本脚本内部派生主 Disease 色表，不修改权威文件。
pal <- copy(dpal[Disease_Abbrev != "ctrl"])
if (!"GBM" %in% pal$Disease_Abbrev) {
  gbm_row <- copy(dpal[Disease_Abbrev == "rGBM"][1])
  gbm_row[, `:=`(Disease_Abbrev = "GBM", Disease_Color_Hex = "#B6362D", Display_Label = "GBM")]
  pal <- rbind(pal, gbm_row, fill = TRUE)
}
pal <- pal[!duplicated(Disease_Abbrev)]
dis_color <- setNames(pal$Disease_Color_Hex, pal$Disease_Abbrev)
dis_color["Ctrl"] <- "#4DAF4A"
dis_color["GBM"] <- "#B6362D"
dis_color["LGG"] <- "#8C1D18"
ds_color <- setNames(spal$Color_Hex, spal$DatasetID)
class_of_all <- setNames(pal$Disease_Class, pal$Disease_Abbrev)
class_of_all["Ctrl"] <- "Control"
class_col <- setNames(
  dpal[!duplicated(Disease_Class), Disease_Class_Color_Hex],
  dpal[!duplicated(Disease_Class), Disease_Class]
)

actual_diseases <- sort(unique(as.character(md$Disease)))
missing_colors <- setdiff(actual_diseases, names(dis_color))
if (length(missing_colors)) stop("以下 Disease 缺少冻结颜色：", paste(missing_colors, collapse = ", "))
if (anyDuplicated(unname(dis_color[actual_diseases]))) {
  stop("实际 Disease 的冻结颜色存在重复，无法保证类别可辨识")
}

md[, elig_sample := .N >= MIN_CELLS_PER_SAMPLE, by = SampleID_Global]
samp_all <- md[, .(
  n_cells = .N,
  DatasetID = first(DatasetID),
  SampleID = first(SampleID),
  Disease = first(Disease),
  elig_sample = first(elig_sample)
), by = SampleID_Global]
sample_disease_check <- md[, .(n_disease = uniqueN(Disease)), by = SampleID_Global]
if (any(sample_disease_check$n_disease > 1L)) {
  stop("同一 SampleID_Global 对应多个 Disease")
}
samp <- samp_all[elig_sample == TRUE]

elig <- samp_all[, .(
  n_cells_total = sum(n_cells),
  n_samples = .N,
  n_elig_samples = sum(elig_sample),
  n_elig_cells = sum(n_cells[elig_sample])
), by = Disease]
elig[, `:=`(
  pass_1000cells = n_cells_total >= 1000,
  pass_3elig_samples = n_elig_samples >= 3
)]
elig[, formal_eligible := pass_1000cells & pass_3elig_samples]
setorder(elig, -n_cells_total)
fwrite(elig, file.path(TBL, "disease_eligibility_audit.csv"))

formal_diseases <- elig[formal_eligible == TRUE, Disease]
if (!"Ctrl" %in% formal_diseases) stop("Ctrl 未通过既定 Figure 1 资格门")
formal_non_ctrl <- sort(setdiff(formal_diseases, "Ctrl"))

em <- Embeddings(obj, "umap_harmony16")
if (!setequal(rownames(em), md$CellID)) stop("UMAP 行名与对象 CellID 不一致")
em <- em[match(md$CellID, rownames(em)), , drop = FALSE]
cu <- data.table(
  CellID = md$CellID,
  SampleID_Global = md$SampleID_Global,
  DatasetID = md$DatasetID,
  Disease = as.character(md$Disease),
  UMAP1 = em[, 1],
  UMAP2 = em[, 2]
)

fwrite(samp[, .(SampleID_Global, DatasetID, SampleID, Disease, n_cells)],
       file.path(TBL, "sample_metadata_eligible.csv"))
fwrite(cu, file.path(SRC, "cell_level_umap_annotations_dims16.csv"))
fwrite(data.table(
  input_path = ATLAS,
  input_sha256 = actual_sha256,
  n_cells = ncol(obj),
  n_genes = nrow(obj),
  assay = "RNA",
  layers = paste(Layers(obj[["RNA"]]), collapse = ";"),
  reduction = "umap_harmony16",
  disease_column = "Disease",
  pseudobulk_method = "Seurat::AggregateExpression on RNA counts grouped by SampleID_Global",
  minimum_cells_per_sample = MIN_CELLS_PER_SAMPLE,
  deg_method = "edgeR quasi-likelihood; each Disease vs Ctrl; FDR < 0.05 and |log2FC| >= 2",
  statistical_unit = "SampleID_Global",
  donor_used = FALSE,
  covariates = "none"
), file.path(AUD, "figure1_input_and_method_audit.csv"))

# ------------------------------- 2. Panel A -------------------------------
log_msg("生成 Panel A")
ds_all <- sort(unique(samp$DatasetID))
if (!all(ds_all %in% names(DATASET_GSE))) {
  stop("DATASET_GSE 未覆盖：", paste(setdiff(ds_all, names(DATASET_GSE)), collapse = ", "))
}
dis_all <- intersect(actual_diseases, names(dis_color))

coverage <- samp[, .(n_cells = sum(n_cells), n_samples = .N), by = .(DatasetID, Disease)]
M <- dcast(coverage, DatasetID ~ Disease, value.var = c("n_cells", "n_samples"), fill = 0)
M <- M[match(ds_all, DatasetID)]
cellm <- as.matrix(M[, grep("^n_cells_", names(M)), with = FALSE])
sampm <- as.matrix(M[, grep("^n_samples_", names(M)), with = FALSE])
colnames(cellm) <- sub("^n_cells_", "", colnames(cellm))
colnames(sampm) <- sub("^n_samples_", "", colnames(sampm))
for (dd in setdiff(dis_all, colnames(cellm))) {
  cellm <- cbind(cellm, 0)
  sampm <- cbind(sampm, 0)
  colnames(cellm)[ncol(cellm)] <- dd
  colnames(sampm)[ncol(sampm)] <- dd
}
cellm <- cellm[, dis_all, drop = FALSE]
sampm <- sampm[, dis_all, drop = FALSE]
cellm[is.na(cellm)] <- 0
sampm[is.na(sampm)] <- 0

luminance <- function(hex) {
  cc <- grDevices::col2rgb(hex) / 255
  cc <- ifelse(cc <= 0.03928, cc / 12.92, ((cc + 0.055) / 1.055)^2.4)
  0.2126 * cc[1] + 0.7152 * cc[2] + 0.0722 * cc[3]
}
panelA_cell_fun <- function(j, i, x, y, w, h, fill) {
  value <- sampm[i, j]
  if (value <= 0) {
    grid.rect(x, y, w, h, gp = gpar(fill = "#F2F2F2", col = "white"))
  } else {
    color <- unname(dis_color[colnames(sampm)[j]])
    grid.rect(x, y, w, h, gp = gpar(fill = color, col = "white"))
    grid.text(as.character(value), x, y,
              gp = gpar(fontsize = 7.5, fontface = "bold", fontfamily = "Arial",
                        col = if (luminance(color) > 0.5) "black" else "white"))
  }
}

dis_cells <- colSums(cellm)
dis_samples <- colSums(sampm)
ds_cells <- rowSums(cellm)
ds_samples <- rowSums(sampm)
format_cell_count <- function(x) {
  ifelse(x >= 1e6, paste0(format(x / 1e6, trim = TRUE), "M"),
         ifelse(x >= 1e3, paste0(format(x / 1e3, trim = TRUE), "k"),
                format(x, trim = TRUE, scientific = FALSE)))
}
make_log10_cell_axis <- function(values) {
  max_value <- max(values, na.rm = TRUE)
  tick_values <- 10^(0:floor(log10(max_value)))
  list(at = log10(tick_values), labels = format_cell_count(tick_values))
}
dis_cell_axis <- make_log10_cell_axis(dis_cells)
ds_cell_axis <- make_log10_cell_axis(ds_cells)
ha_top <- HeatmapAnnotation(
  Cells = anno_barplot(log10(pmax(dis_cells, 1)), gp = gpar(fill = "#2C7FB8"), border = TRUE,
                       ylim = c(0, max(log10(pmax(dis_cells, 1))) * 1.15),
                       axis_param = dis_cell_axis),
  Samples = anno_barplot(dis_samples, gp = gpar(fill = "#238B45"), border = TRUE,
                         ylim = c(0, max(dis_samples) * 1.15)),
  annotation_height = unit(c(1.25, 1.1), "cm"),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 7, fontfamily = "Arial")
)
ha_right <- rowAnnotation(
  Cells = anno_barplot(log10(pmax(ds_cells, 1)), gp = gpar(fill = "#2C7FB8"),
                       border = TRUE,
                       ylim = c(0, max(log10(pmax(ds_cells, 1))) * 1.15),
                       axis_param = ds_cell_axis),
  Samples = anno_barplot(ds_samples, gp = gpar(fill = "#238B45"), border = TRUE),
  annotation_width = unit(c(1.9, 1.6), "cm"),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 7, fontfamily = "Arial")
)
legend_A <- Legend(
  title = "Disease",
  at = dis_all,
  legend_gp = gpar(fill = dis_color[dis_all]),
  ncol = 1,
  grid_width = unit(3.8, "mm"),
  grid_height = unit(3.8, "mm"),
  title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Arial"),
  labels_gp = gpar(fontsize = 7, fontfamily = "Arial")
)
row_labels_A <- unname(DATASET_GSE[ds_all])
htA <- Heatmap(
  cellm,
  name = "coverage",
  col = colorRamp2(c(0, 1), c("white", "white")),
  width = unit(ncol(cellm) * 7.5, "mm"),
  height = unit(nrow(cellm) * 7.5, "mm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cell_fun = panelA_cell_fun,
  rect_gp = gpar(col = "grey90", lwd = 0.4),
  row_labels = row_labels_A,
  row_names_side = "left",
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 8, fontfamily = "Arial"),
  column_names_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Arial",
                          col = "black"),
  column_names_rot = 45,
  top_annotation = ha_top,
  right_annotation = ha_right,
  show_heatmap_legend = FALSE
)
title_A <- sprintf(
  "Dataset x Disease coverage (%s cells; numbers in tiles = eligible samples)",
  format(sum(cellm), big.mark = ",")
)
png(file.path(FIG, "PanelA_dataset_disease_coverage.png"),
    width = 12, height = 8.2, units = "in", res = 600, type = "cairo")
draw(htA, column_title = title_A,
     column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial"),
     annotation_legend_list = list(legend_A), padding = unit(c(2, 12, 2, 2), "mm"))
dev.off()
cairo_pdf(file.path(FIG, "PanelA_dataset_disease_coverage.pdf"),
          width = 12, height = 8.2, family = "Arial")
draw(htA, column_title = title_A,
     column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial"),
     annotation_legend_list = list(legend_A), padding = unit(c(2, 12, 2, 2), "mm"))
dev.off()
coverage[, `:=`(
  GSE = unname(DATASET_GSE[DatasetID]),
  Statistical_Unit = "eligible sample",
  Minimum_cells_per_sample = MIN_CELLS_PER_SAMPLE
)]
fwrite(coverage, file.path(SRC, "panelA_coverage_data.csv"))

# ------------------------------- 3. Panel B -------------------------------
log_msg("生成 Panel B")
set.seed(RANDOM_SEED)
cu_plot <- cu[sample(.N)]
pB <- ggplot(cu_plot, aes(UMAP1, UMAP2, colour = Disease)) +
  geom_point(size = 0.18, alpha = 0.6, stroke = 0) +
  scale_colour_manual(values = dis_color[actual_diseases], breaks = actual_diseases,
                      name = "Disease") +
  coord_fixed(expand = FALSE) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.30, "cm")
  ) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.5, alpha = 1))) +
  labs(title = "Disease distribution in the integrated microglia atlas")
save_gg(pB, file.path(FIG, "PanelB_disease_umap"), 8.2, 6)
fwrite(cu, file.path(SRC, "panelB_umap_data.csv"))
rm(cu_plot)
invisible(gc())

# ------------------------------- 4. Panel C -------------------------------
log_msg("生成 Panel C")
eligible_sample_ids <- samp$SampleID_Global
cent <- cu[SampleID_Global %in% eligible_sample_ids,
           .(center_UMAP1 = mean(UMAP1), center_UMAP2 = mean(UMAP2)),
           by = .(SampleID_Global, Disease)]
ctrl_center <- cent[Disease == "Ctrl", .(cx = mean(center_UMAP1), cy = mean(center_UMAP2))]
if (nrow(ctrl_center) != 1L) stop("无法计算统一 Ctrl 中心")
dis_cent <- cent[Disease %in% formal_non_ctrl,
                 .(center_UMAP1 = mean(center_UMAP1), center_UMAP2 = mean(center_UMAP2),
                   n_samples = uniqueN(SampleID_Global)),
                 by = Disease]
dis_cent[, `:=`(
  ctrl_UMAP1 = ctrl_center$cx,
  ctrl_UMAP2 = ctrl_center$cy,
  displacement_UMAP1 = ctrl_center$cx - center_UMAP1,
  displacement_UMAP2 = ctrl_center$cy - center_UMAP2
)]
setorder(dis_cent, Disease)
fwrite(dis_cent, file.path(SRC, "panelC_displacement_data.csv"))

dis_cent[, Disease_plot := factor(Disease, levels = formal_non_ctrl)]
pC <- ggplot(dis_cent, aes(center_UMAP1, center_UMAP2, colour = Disease_plot)) +
  geom_segment(
    aes(xend = ctrl_UMAP1, yend = ctrl_UMAP2),
    arrow = arrow(length = unit(0.14, "cm")), linewidth = 0.55,
    alpha = 0.82, show.legend = FALSE
  ) +
  geom_point(size = 3.2) +
  geom_text_repel(
    aes(label = Disease), size = 3, family = "Arial",
    show.legend = FALSE, max.overlaps = Inf, box.padding = 0.35,
    point.padding = 0.18, segment.color = "grey65"
  ) +
  annotate("point", x = ctrl_center$cx, y = ctrl_center$cy,
           colour = "#4DAF4A", size = 3.8, shape = 17) +
  annotate("text", x = ctrl_center$cx, y = ctrl_center$cy,
           label = "Ctrl", colour = "#4DAF4A", size = 3.2,
           family = "Arial", fontface = "bold", vjust = 1.8) +
  scale_colour_manual(
    values = dis_color[formal_non_ctrl], breaks = formal_non_ctrl,
    name = "Disease", drop = FALSE
  ) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  coord_equal() +
  theme_classic(base_family = "Arial", base_size = 9) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.28, "cm")
  ) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2.5))) +
  labs(
    title = "Disease sample-centroid displacement relative to Ctrl",
    subtitle = "Descriptive UMAP geometry; arrows terminate at the unified Ctrl centroid",
    x = "UMAP 1", y = "UMAP 2", colour = "Disease"
  )
save_gg(pC, file.path(FIG, "PanelC_disease_sample_centroid_displacement"), 8.4, 5.4)

# ---------------- 5. Seurat pseudobulk + sample-level edgeR DEG ------------
log_msg("使用 Seurat::AggregateExpression 构建 SampleID_Global pseudobulk counts")
DefaultAssay(obj) <- "RNA"
gene_universe_all <- rownames(obj)

# 只有细胞数 >=20 且属于正式疾病集合的 sample 才允许进入 pseudobulk。
pb_samples <- copy(samp[Disease %in% formal_diseases])
if (!nrow(pb_samples) || any(pb_samples$n_cells < MIN_CELLS_PER_SAMPLE)) {
  stop("pseudobulk 样本门槛失败：存在少于 20 个细胞的 sample")
}
setorder(pb_samples, SampleID_Global)
pb_samples[, PseudobulkID := sprintf("PB%04d", seq_len(.N))]
pb_id_by_sample <- setNames(pb_samples$PseudobulkID, pb_samples$SampleID_Global)
eligible_cells <- md[SampleID_Global %in% pb_samples$SampleID_Global, CellID]
obj_pb <- subset(obj, cells = eligible_cells)
obj_pb$PseudobulkID <- unname(pb_id_by_sample[as.character(obj_pb$SampleID_Global)])
if (anyNA(obj_pb$PseudobulkID)) stop("PseudobulkID 映射存在缺失")
pb_cell_audit <- as.data.table(table(PseudobulkID = obj_pb$PseudobulkID))
setnames(pb_cell_audit, "N", "n_cells_from_object")
if (any(pb_cell_audit$n_cells_from_object < MIN_CELLS_PER_SAMPLE)) {
  stop("Seurat 子对象中存在少于 20 个细胞的 pseudobulk sample")
}

# 正式 pseudobulk 生成函数：Seurat::AggregateExpression()。
pb_counts_all <- AggregateExpression(
  object = obj_pb,
  assays = "RNA",
  group.by = "PseudobulkID",
  return.seurat = FALSE,
  verbose = TRUE
)$RNA
if (!all(pb_samples$PseudobulkID %in% colnames(pb_counts_all))) {
  stop("AggregateExpression 输出缺少合格 PseudobulkID")
}
pb_counts <- pb_counts_all[, pb_samples$PseudobulkID, drop = FALSE]
if (anyDuplicated(colnames(pb_counts))) stop("sample pseudobulk 列名不唯一")
pb_values <- if (inherits(pb_counts, "sparseMatrix")) pb_counts@x else as.numeric(pb_counts)
if (any(!is.finite(pb_values)) || any(pb_values < 0) ||
    any(abs(pb_values - round(pb_values)) > 1e-8)) {
  stop("AggregateExpression 输出不是非负整数 pseudobulk counts")
}
saveRDS(pb_counts, file.path(TBL, "sample_pseudobulk_counts.rds"), compress = FALSE)

pb_meta <- pb_samples[match(colnames(pb_counts), PseudobulkID),
                      .(PseudobulkID, SampleID_Global, DatasetID, SampleID, Disease, n_cells)]
if (anyNA(pb_meta$Disease) || !identical(pb_meta$PseudobulkID, colnames(pb_counts))) {
  stop("sample pseudobulk metadata 对齐失败")
}
fwrite(pb_meta,
       file.path(TBL, "sample_pseudobulk_metadata.csv"))
fwrite(merge(pb_meta, pb_cell_audit, by = "PseudobulkID", all.x = TRUE),
       file.path(AUD, "sample_pseudobulk_min20cells_audit.csv"))

# edgeR 官方 count-native quasi-likelihood 流程；设计矩阵只含 Disease。
log_msg("运行 edgeR QL：每个 Disease vs Ctrl；统计单位 = sample；无协变量")
group <- factor(pb_meta$Disease, levels = c("Ctrl", formal_non_ctrl))
if (anyNA(group)) stop("edgeR group 存在未登记 Disease")
pb_meta_df <- as.data.frame(pb_meta)
rownames(pb_meta_df) <- pb_meta_df$PseudobulkID
dge <- DGEList(counts = pb_counts, samples = pb_meta_df)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
rownames(design) <- pb_meta$PseudobulkID
keep <- filterByExpr(dge, design = design)
if (!any(keep)) stop("edgeR::filterByExpr 未保留任何基因")
fwrite(data.table(gene = rownames(dge), keep_by_filterByExpr = keep),
       file.path(AUD, "edgeR_filterByExpr_audit.csv"))
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
dge <- estimateDisp(dge, design = design, robust = TRUE)
fit <- glmQLFit(dge, design = design, robust = TRUE)
gene_universe_tested <- rownames(dge)

deg_list <- vector("list", length(formal_non_ctrl))
names(deg_list) <- formal_non_ctrl
for (disease in formal_non_ctrl) {
  contrast <- rep(0, ncol(design))
  names(contrast) <- colnames(design)
  contrast[disease] <- 1
  contrast["Ctrl"] <- -1
  qlf <- glmQLFTest(fit, contrast = contrast)
  tt <- as.data.table(topTags(qlf, n = Inf, sort.by = "none")$table,
                      keep.rownames = "gene")
  tt[, `:=`(
    Disease = disease,
    Comparison = paste0(disease, " vs Ctrl"),
    Significant = FDR < 0.05 & abs(logFC) >= 2,
    Direction = fifelse(FDR < 0.05 & logFC >= 2, "UP",
                        fifelse(FDR < 0.05 & logFC <= -2, "DOWN", "NS")),
    Method = "edgeR quasi-likelihood",
    Design = "~0 + Disease",
    Statistical_Unit = "SampleID_Global",
    Donor_Used = FALSE,
    Covariates = "none",
    DEG_Threshold = "FDR < 0.05 and |log2FC| >= 2"
  )]
  deg_list[[disease]] <- tt
}
deg_all <- rbindlist(deg_list, use.names = TRUE, fill = TRUE)
deg_sig <- deg_all[Significant == TRUE & Direction %in% c("UP", "DOWN")]
if (!nrow(deg_sig)) stop("edgeR 未检出 FDR < 0.05 且 |log2FC| >= 2 的 DEG")
fwrite(deg_all, file.path(TBL, "disease_vs_ctrl_edgeR_all_tested_genes.csv"))
fwrite(deg_sig, file.path(TBL, "disease_vs_ctrl_edgeR_FDR05_abslogFC2_DEGs.csv"))
fwrite(data.table(
  Disease = formal_non_ctrl,
  n_disease_samples = vapply(formal_non_ctrl, function(x) sum(pb_meta$Disease == x), integer(1)),
  n_ctrl_samples = sum(pb_meta$Disease == "Ctrl"),
  design = "~0 + Disease",
  comparison = paste0(formal_non_ctrl, " vs Ctrl"),
  covariates = "none",
  donor_used = FALSE,
  minimum_cells_per_sample = MIN_CELLS_PER_SAMPLE,
  deg_threshold = "FDR < 0.05 and |log2FC| >= 2"
), file.path(AUD, "edgeR_comparison_audit.csv"))

# DEG 结束后释放大对象；后续只用 DEG 表、UMAP 表和 tested universe。
rm(obj, obj_pb, md, em, pb_counts_all, dge, fit)
invisible(gc())

# ------------------------------- 6. Panel D -------------------------------
log_msg("生成 Panel D")
disease_order <- formal_non_ctrl
burden <- deg_sig[, .(
  n_up = sum(Direction == "UP"),
  n_down = sum(Direction == "DOWN"),
  n_total_DEGs = .N
), by = Disease]
burden <- merge(data.table(Disease = disease_order), burden, by = "Disease", all.x = TRUE)
burden[is.na(n_up), `:=`(n_up = 0L, n_down = 0L, n_total_DEGs = 0L)]
burden[, Disease_Class := class_of_all[Disease]]
burden[, Disease_plot := factor(Disease, levels = disease_order)]
fwrite(burden[, .(
  Disease, Disease_Class, n_up, n_down, n_total_DEGs,
  Comparison = paste0(Disease, " vs Ctrl"),
  Pseudobulk_Method = "Seurat::AggregateExpression",
  Method = "edgeR quasi-likelihood; TMM; filterByExpr; FDR < 0.05 and |log2FC| >= 2",
  Statistical_Unit = "SampleID_Global",
  Donor_Used = FALSE,
  Covariates = "none"
)], file.path(SRC, "panelD_deg_burden_data.csv"))

pD <- ggplot(burden, aes(Disease_plot)) +
  geom_col(aes(y = n_up), fill = "#D7301F", width = 0.75) +
  geom_col(aes(y = -n_down), fill = "#2166AC", width = 0.75) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "black") +
  geom_text(aes(y = n_up, label = format(n_up, big.mark = ",")),
            vjust = -0.35, size = 2.4, family = "Arial") +
  geom_text(aes(y = -n_down, label = format(n_down, big.mark = ",")),
            vjust = 1.35, size = 2.4, family = "Arial") +
  scale_y_continuous(labels = function(x) abs(x), expand = expansion(mult = c(0.08, 0.12))) +
  theme_classic(base_family = "Arial", base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,
                               colour = "black"),
    legend.position = "none"
  ) +
  labs(
    title = "Disease-associated microglial DEG burden",
    subtitle = "Seurat AggregateExpression sample pseudobulk; edgeR QL; Disease vs Ctrl; FDR < 0.05 and |log2FC| >= 2",
    x = NULL, y = "Number of DEGs"
  )
save_gg(pD, file.path(FIG, "PanelD_edgeR_sample_pseudobulk_burden"), 8.2, 4.2)

# ------------------------------- 7. Panel E -------------------------------
log_msg("生成 Panel E：上三角 UP、下三角 DOWN 的拼接 Jaccard 矩阵")
sets <- deg_sig[, .(genes = list(unique(gene))), by = .(Disease, Direction)]
get_set <- function(disease, direction) {
  ans <- sets[Disease == disease & Direction == direction, genes]
  if (!length(ans)) character() else ans[[1]]
}
jaccard_matrix <- function(direction) {
  mat <- matrix(NA_real_, length(disease_order), length(disease_order),
                dimnames = list(disease_order, disease_order))
  for (i in seq_along(disease_order)) {
    for (j in seq_along(disease_order)) {
      a <- get_set(disease_order[i], direction)
      b <- get_set(disease_order[j], direction)
      union_ab <- union(a, b)
      if (length(union_ab) == 0L) {
        mat[i, j] <- NA_real_
      } else {
        mat[i, j] <- length(intersect(a, b)) / length(union_ab)
      }
    }
  }
  mat
}
j_up <- jaccard_matrix("UP")
j_down <- jaccard_matrix("DOWN")

j_long <- rbindlist(list(
  as.data.table(as.table(j_up))[, `:=`(Direction = "UP")],
  as.data.table(as.table(j_down))[, `:=`(Direction = "DOWN")]
), use.names = TRUE)
setnames(j_long, names(j_long)[1:3], c("Disease_A", "Disease_B", "Jaccard"))
j_long[, `:=`(
  Disease_A = as.character(Disease_A),
  Disease_B = as.character(Disease_B),
  DEG_definition = "Seurat AggregateExpression sample pseudobulk; edgeR FDR < 0.05 and |log2FC| >= 2",
  Comparison = "Disease vs Ctrl",
  Statistical_Unit = "SampleID_Global",
  Donor_Used = FALSE,
  Covariates = "none",
  Display_Triangle = fifelse(Direction == "UP", "upper", "lower")
)]
fwrite(j_long, file.path(SRC, "panelE_jaccard.csv"))

mean_similarity <- (j_up + j_down) / 2
mean_similarity[is.na(mean_similarity)] <- 0
mean_similarity <- (mean_similarity + t(mean_similarity)) / 2
diag(mean_similarity) <- 1
if (length(disease_order) > 2L && any(mean_similarity[upper.tri(mean_similarity)] > 0)) {
  common_order <- disease_order[hclust(as.dist(1 - mean_similarity), method = "complete")$order]
} else {
  common_order <- disease_order
}
j_up <- j_up[common_order, common_order, drop = FALSE]
j_down <- j_down[common_order, common_order, drop = FALSE]

all_j <- c(j_up[upper.tri(j_up)], j_down[lower.tri(j_down)])
all_j <- all_j[is.finite(all_j)]
color_max <- if (length(all_j) && max(all_j) > 0) max(all_j) else 1
color_mid <- color_max / 2
col_up <- colorRamp2(c(0, color_mid, color_max), c("#FFF5F0", "#FB6A4A", "#99000D"))
col_down <- colorRamp2(c(0, color_mid, color_max), c("#F7FBFF", "#6BAED6", "#08306B"))
signed_col_fun <- colorRamp2(
  c(-color_max, -color_mid, 0, color_mid, color_max),
  c("#08306B", "#6BAED6", "#FFFFFF", "#FB6A4A", "#99000D")
)
J <- matrix(NA_real_, length(common_order), length(common_order),
            dimnames = list(common_order, common_order))
upper_cells <- upper.tri(J)
lower_cells <- lower.tri(J)
J[upper_cells] <- j_up[upper_cells]
J[lower_cells] <- -j_down[lower_cells]
diag(J) <- NA_real_

panelE_cell_fun <- function(j, i, x, y, w, h, fill) {
  value <- J[i, j]
  if (i == j) {
    grid.text("—", x, y,
              gp = gpar(fontsize = 7, fontfamily = "Arial", col = "grey45"))
  } else if (is.na(value)) {
    grid.text("NA", x, y,
              gp = gpar(fontsize = 5.8, fontfamily = "Arial", col = "grey55"))
  } else {
    text_col <- if (abs(value) >= color_mid) "white" else "black"
    grid.text(sprintf("%.2f", abs(value)), x, y,
              gp = gpar(fontsize = 6.2, fontfamily = "Arial", col = text_col))
  }
}
legend_ticks <- unique(round(c(0, color_mid, color_max), 2))
legend_up <- Legend(
  title = "UP Jaccard\n(upper triangle)", col_fun = col_up, at = legend_ticks,
  title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Arial"),
  labels_gp = gpar(fontsize = 7, fontfamily = "Arial")
)
legend_down <- Legend(
  title = "DOWN Jaccard\n(lower triangle)", col_fun = col_down, at = legend_ticks,
  title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Arial"),
  labels_gp = gpar(fontsize = 7, fontfamily = "Arial")
)
htE <- Heatmap(
  J,
  name = "signed_jaccard",
  col = signed_col_fun,
  na_col = "#F2F2F2",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  width = unit(length(common_order) * 8.2, "mm"),
  height = unit(length(common_order) * 8.2, "mm"),
  row_names_gp = gpar(fontsize = 8, fontfamily = "Arial", col = "black"),
  column_names_gp = gpar(fontsize = 8, fontfamily = "Arial", col = "black"),
  column_names_rot = 45,
  rect_gp = gpar(col = "grey88", lwd = 0.45),
  cell_fun = panelE_cell_fun,
  show_heatmap_legend = FALSE
)
panelE_title <- "DEG overlap Jaccard (upper: UP; lower: DOWN; edgeR Disease vs Ctrl)"
panelE_stem <- file.path(FIG, "PanelE_split_triangle_Jaccard")
png(paste0(panelE_stem, ".png"), width = 8.8, height = 8.2,
    units = "in", res = 600, type = "cairo")
draw(htE, column_title = panelE_title,
     column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial"),
     heatmap_legend_list = list(legend_up, legend_down),
     heatmap_legend_side = "right", padding = unit(c(2, 16, 2, 2), "mm"))
dev.off()
cairo_pdf(paste0(panelE_stem, ".pdf"), width = 8.8, height = 8.2, family = "Arial")
draw(htE, column_title = panelE_title,
     column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial"),
     heatmap_legend_list = list(legend_up, legend_down),
     heatmap_legend_side = "right", padding = unit(c(2, 16, 2, 2), "mm"))
dev.off()

# ------------------------------- 8. Panel F -------------------------------
log_msg("运行 GO Biological Process ORA 并生成 Panel F")
gobp_gs <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
TERM2GENE <- unique(data.frame(term = gobp_gs$gs_name, gene = gobp_gs$gene_symbol))

ora_list <- list()
ora_audit <- list()
for (disease in disease_order) {
  for (direction in c("UP", "DOWN")) {
    genes <- unique(deg_sig[Disease == disease & Direction == direction, gene])
    run_ora <- length(genes) >= 5L
    ora_audit[[length(ora_audit) + 1L]] <- data.table(
      Disease = disease,
      Direction = direction,
      n_input_genes = length(genes),
      ora_run = run_ora,
      DEG_definition = "Seurat AggregateExpression sample pseudobulk; edgeR FDR < 0.05 and |log2FC| >= 2",
      Statistical_Unit = "SampleID_Global",
      Donor_Used = FALSE,
      ORA_background = "all protein-coding genes in merged_obj.rds"
    )
    if (!run_ora) next
    er <- tryCatch(
      enricher(
        gene = genes,
        TERM2GENE = TERM2GENE,
        universe = gene_universe_tested,
        pAdjustMethod = "BH",
        minGSSize = 5,
        maxGSSize = 500
      ),
      error = function(e) {
        log_msg("GO ORA 失败：", disease, " / ", direction, "：", conditionMessage(e))
        NULL
      }
    )
    if (is.null(er) || nrow(as.data.frame(er)) == 0L) next
    dd <- as.data.table(as.data.frame(er))
    dd[, `:=`(Disease = disease, Direction = direction)]
    ora_list[[length(ora_list) + 1L]] <- dd
  }
}
ora_all <- if (length(ora_list)) rbindlist(ora_list, fill = TRUE) else data.table()
ora_audit <- rbindlist(ora_audit, fill = TRUE)
fwrite(ora_all, file.path(SRC, "panelF_go_bp_ora_all.csv"))
fwrite(ora_audit, file.path(SRC, "panelF_go_bp_selection_audit.csv"))

readable_term <- function(x) {
  x <- sub("^GOBP_", "", x)
  x <- gsub("_", " ", x)
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}

make_go_bubble <- function(direction, top_k, stem) {
  title <- paste0("Top ", top_k, " GO biological processes per disease - ", direction)
  if (!nrow(ora_all)) {
    selected <- data.table()
  } else {
    selected <- copy(ora_all[Direction == direction & p.adjust < 0.05])
  }
  if (nrow(selected)) {
    setorder(selected, Disease, p.adjust)
    selected <- selected[, head(.SD, top_k), by = Disease]
    selected[, Term_Readable := readable_term(ID)]
    selected[, Disease_Class := class_of_all[Disease]]
  }
  fwrite(selected, paste0(stem, "_source.csv"))
  if (!nrow(selected)) {
    p <- ggplot() + theme_void(base_family = "Arial") +
      labs(title = paste0(title, " - no significant terms"))
    save_gg(p, stem, 10, 4)
    return(invisible(NULL))
  }
  selected[, Disease_plot := factor(Disease, levels = disease_order)]
  disease_class_order <- unique(na.omit(class_of_all[disease_order]))
  selected[, Disease_Class := factor(Disease_Class, levels = disease_class_order)]
  term_order <- selected[, .(n_diseases = uniqueN(Disease), min_fdr = min(p.adjust)),
                         by = Term_Readable]
  setorder(term_order, -n_diseases, min_fdr)
  selected[, Term_plot := factor(Term_Readable, levels = rev(term_order$Term_Readable))]
  p <- ggplot(selected, aes(Disease_plot, Term_plot,
                            size = -log10(p.adjust), colour = Disease_plot)) +
    geom_point(alpha = 0.9) +
    scale_colour_manual(values = dis_color[disease_order], guide = "none") +
    scale_size_continuous(range = c(2.2, 7), name = expression(-log[10](FDR))) +
    facet_grid(. ~ Disease_Class, scales = "free_x", space = "free_x", drop = TRUE) +
    theme_classic(base_family = "Arial", base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1,
                                 colour = "black"),
      strip.background = element_rect(fill = "grey95", colour = "grey70", linewidth = 0.3),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.position = "right",
      plot.margin = margin(6, 8, 6, 6)
    ) +
    labs(
      title = title,
      subtitle = "Seurat AggregateExpression sample pseudobulk; edgeR Disease vs Ctrl; FDR < 0.05 and |log2FC| >= 2",
      x = NULL, y = "GO Biological Process"
    )
  height <- max(6, uniqueN(selected$Term_Readable) * 0.26 + 2.5)
  save_gg(p, stem, 11.5, height)
}

make_go_bubble("UP", TOP_K_MAIN, file.path(FIG, "PanelF_GO_BP_UP_top3"))
make_go_bubble("DOWN", TOP_K_MAIN, file.path(FIG, "PanelF_GO_BP_DOWN_top3"))
make_go_bubble("UP", TOP_K_TABLE, file.path(FIG, "PanelF_GO_BP_UP_top5"))
make_go_bubble("DOWN", TOP_K_TABLE, file.path(FIG, "PanelF_GO_BP_DOWN_top5"))

# ------------------------------- 9. QA 与解读 ------------------------------
expected_stems <- c(
  "PanelA_dataset_disease_coverage",
  "PanelB_disease_umap",
  "PanelC_disease_sample_centroid_displacement",
  "PanelD_edgeR_sample_pseudobulk_burden",
  "PanelE_split_triangle_Jaccard",
  "PanelF_GO_BP_UP_top3",
  "PanelF_GO_BP_DOWN_top3",
  "PanelF_GO_BP_UP_top5",
  "PanelF_GO_BP_DOWN_top5"
)
qa <- rbindlist(lapply(expected_stems, function(stem) {
  data.table(
    panel = stem,
    png_exists = file.exists(file.path(FIG, paste0(stem, ".png"))),
    pdf_exists = file.exists(file.path(FIG, paste0(stem, ".pdf")))
  )
}))
qa[, passed := png_exists & pdf_exists]
fwrite(qa, file.path(AUD, "figure1_output_qa.csv"))
if (!all(qa$passed)) stop("部分 Figure 1 图件未成功导出，详见 figure1_output_qa.csv")

deg_summary <- burden[, paste0(Disease, ": UP=", n_up, ", DOWN=", n_down)]
interpretation <- c(
  "# Figure 1 解读",
  "",
  "## 分析目的",
  "展示整合小胶质细胞图谱的疾病覆盖、Disease UMAP 分布、样本质心位移、疾病相关 marker 数量、marker 集合共享程度及 GO Biological Process 富集。",
  "",
  "## 输入数据",
  paste0("- Seurat 对象：`", ATLAS, "`"),
  paste0("- 细胞数：", format(nrow(cu), big.mark = ",")),
  paste0("- 蛋白编码基因数：", format(length(gene_universe_all), big.mark = ",")),
  paste0("- Dataset 数：", uniqueN(cu$DatasetID)),
  paste0("- Disease 标签数：", uniqueN(cu$Disease)),
  "",
  "## 方法",
  "| 模块 | 方法 |",
  "|---|---|",
  "| Panel A | Dataset x Disease 合格样本覆盖；每个样本至少 20 个细胞 |",
  "| Panel B-C | `umap_harmony16`；Panel C 为样本质心的描述性几何位移 |",
  "| DEG | 仅保留每个 sample >=20 个细胞；`Seurat::AggregateExpression()` 按 `SampleID_Global` 聚合 RNA counts；edgeR TMM + filterByExpr + quasi-likelihood；每病 vs Ctrl；FDR < 0.05 且 |log2FC| >= 2；不使用 donor 或协变量 |",
  "| Panel E | edgeR 显著 DEG 集合的 Jaccard；上三角为 UP、下三角为 DOWN，拼接为同一热图 |",
  "| Panel F | GO:BP ORA；BH 校正；背景为 edgeR 实际检验的蛋白编码基因；按 Disease_Class 分面 |",
  "",
  "## 关键结果",
  paste0("- sample pseudobulk 数：", ncol(pb_counts), "。"),
  paste0("- edgeR 实际检验基因数：", format(length(gene_universe_tested), big.mark = ","), "。"),
  paste0("- edgeR FDR < 0.05 且 |log2FC| >= 2 DEG 总数：", format(nrow(deg_sig), big.mark = ","), "。"),
  paste0("- 正式非 Ctrl 疾病：", paste(disease_order, collapse = ", "), "。"),
  paste0("- ", deg_summary),
  "",
  "## 生物学意义",
  "Panel D-F 描述各 Disease samples 相对于 Ctrl samples 的差异表达。Jaccard 高值表示显著 DEG 集合共享较多，GO 富集用于概括这些集合的功能主题。",
  "",
  "## 限制",
  "edgeR 以 SampleID_Global pseudobulk 为统计观测，每个疾病与 Ctrl 比较；未控制 Dataset、年龄、性别或其他混杂因素，也不使用 donor。因此结果是未校正混杂因素的 sample-level 疾病相关差异，不应表述为独立于数据来源或临床变量的效应。"
)
writeLines(interpretation, file.path(ROOT, "解读.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(LOG, "sessionInfo_all_in_one.txt"), useBytes = TRUE)

log_msg("Figure 1 全流程完成。输出目录：", ROOT)
