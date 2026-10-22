windows_locale <- tryCatch(
  Sys.setlocale("LC_CTYPE", "Chinese_China.utf8"),
  warning = function(w) NA_character_,
  error = function(e) NA_character_
)
if (is.na(windows_locale)) stop("无法设置 Windows 中文 LC_CTYPE，无法安全处理中文输出路径。")

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

atlas_path <- "E:/micro_output/v3_dims16/objects/merged_obj.rds"
palette_path <- "D:/data/microglia_new/authority/03_reference_configs/disease_palette_contract.csv"
output_root <- "E:/micro_output/分析审阅中心/待审阅/Fig1A_全部样本_排除两个混合GSM"
excluded_gsm <- c("GSM6186080", "GSM6186087")
sample_note <- "纳入当前 Seurat 对象全部样本，排除两个混合疾病文库 GSM6186080、GSM6186087；样本按 SampleID_Global 计数。"

dataset_gse <- c(
  ABDELFATTAH_GBM_2022 = "GSE182109",
  BOSTRAND_HD_GLIA_2024 = "GSE281069",
  GSE103224 = "GSE103224",
  GSE141383 = "GSE141383",
  KISS_2023 = "GSE227781",
  LIM_2022 = "GSE180928",
  LIU_SEOAD_2025 = "GSE272082",
  LUNG_BRAIN_MET_2022 = "GSE202371",
  MARTIROSYAN_2024 = "GSE243639",
  MEDULLOBLASTOMA_2021 = "GSE155446",
  SMITH_2022 = "GSE160936",
  TAGORE_2025_NSCLC_BM = "GSE223499",
  TUDDENHAM_2024 = "GSE204702",
  WANG_C9ORF72_ALS_FTD_2024 = "GSE219280",
  WANG_C9_DLPFC_PNAS_2025 = "GSE212630",
  ZHU_2024 = "GSE202210"
)

for (p in c(atlas_path, palette_path)) {
  if (!file.exists(p)) stop("输入文件不存在：", p)
}
if (!dir.exists(dirname(output_root)) && !dir.create(dirname(output_root), recursive = TRUE)) {
  stop("无法创建输出父目录：", dirname(output_root))
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "figures"), showWarnings = FALSE)
dir.create(file.path(output_root, "tables"), showWarnings = FALSE)
dir.create(file.path(output_root, "source_data"), showWarnings = FALSE)

log_path <- file.path(output_root, "Fig1A_run.log")
log_msg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_path, append = TRUE)
}

log_msg("读取当前 Seurat 对象：", atlas_path)
obj <- readRDS(atlas_path)
if (!inherits(obj, "Seurat")) stop("输入对象不是 Seurat 对象")
assay_names <- Assays(obj)
default_assay <- DefaultAssay(obj)
rna_layers <- if ("RNA" %in% assay_names) Layers(obj[["RNA"]]) else character()
md <- as.data.table(obj[[]], keep.rownames = "CellID")
if (nrow(md) != ncol(obj)) stop("meta.data 行数与 Seurat 细胞数不一致")
if (!identical(md$CellID, colnames(obj))) {
  md <- md[match(colnames(obj), CellID)]
}
if (anyNA(md$CellID) || anyDuplicated(md$CellID)) stop("CellID 缺失或重复")

required_cols <- c("DatasetID", "SampleID_Global", "SampleID", "Disease", "GSM")
missing_cols <- setdiff(required_cols, names(md))
if (length(missing_cols)) stop("meta.data 缺少字段：", paste(missing_cols, collapse = ", "))
for (nm in c("DatasetID", "SampleID_Global", "Disease")) {
  md[, (nm) := as.character(get(nm))]
  if (anyNA(md[[nm]]) || any(!nzchar(trimws(md[[nm]])))) {
    stop("关键字段 ", nm, " 存在缺失/空值；已停止汇总，未静默丢弃。")
  }
}
for (nm in c("SampleID", "GSM")) md[, (nm) := as.character(get(nm))]

sample_map <- md[, .(
  n_DatasetID = uniqueN(DatasetID),
  n_Disease = uniqueN(Disease),
  n_SampleID = uniqueN(SampleID),
  n_GSM = uniqueN(GSM),
  DatasetID_values = paste(sort(unique(DatasetID)), collapse = ";"),
  Disease_values = paste(sort(unique(Disease)), collapse = ";"),
  SampleID_values = paste(sort(unique(na.omit(SampleID))), collapse = ";"),
  GSM_values = paste(sort(unique(na.omit(GSM))), collapse = ";")
), by = SampleID_Global]
sample_map[, mapping_status := fifelse(n_DatasetID == 1L & n_Disease == 1L,
                                       "unique DatasetID and Disease", "CONFLICT")]
fwrite(sample_map, file.path(output_root, "tables", "sample_label_mapping_audit.csv"))
conflicts <- sample_map[mapping_status == "CONFLICT"]
if (nrow(conflicts)) {
  fwrite(conflicts, file.path(output_root, "tables", "sample_label_conflicts.csv"))
  stop("发现 SampleID_Global 对应多个 DatasetID 或 Disease；冲突已输出，未任选标签。")
}

meta_before <- copy(md)
before_cells <- nrow(meta_before)
before_samples <- uniqueN(meta_before$SampleID_Global)
gsm_audit <- meta_before[GSM %chin% excluded_gsm, .(
  DatasetID = paste(sort(unique(DatasetID)), collapse = ";"),
  SampleID_Global = paste(sort(unique(SampleID_Global)), collapse = ";"),
  SampleID = paste(sort(unique(na.omit(SampleID))), collapse = ";"),
  Disease = paste(sort(unique(Disease)), collapse = ";"),
  n_samples = uniqueN(SampleID_Global),
  n_cells = .N
), by = GSM]
missing_exclusions <- setdiff(excluded_gsm, gsm_audit$GSM)
if (length(missing_exclusions)) {
  missing_rows <- data.table(
    GSM = missing_exclusions, DatasetID = NA_character_,
    SampleID_Global = NA_character_, SampleID = NA_character_,
    Disease = NA_character_, n_samples = 0L, n_cells = 0L
  )
  gsm_audit <- rbind(gsm_audit, missing_rows, fill = TRUE)
}
gsm_audit[, exclusion_order := match(GSM, excluded_gsm)]
setorder(gsm_audit, exclusion_order)
gsm_audit[, exclusion_order := NULL]
fwrite(gsm_audit, file.path(output_root, "tables", "excluded_gsm_audit.csv"))
if (any(gsm_audit$n_cells == 0L)) {
  stop("指定 GSM 未在对象 GSM 字段中精确匹配到细胞；排除审计表已写出。")
}
md <- meta_before[!GSM %chin% excluded_gsm]
after_cells <- nrow(md)
after_samples <- uniqueN(md$SampleID_Global)
removed_cells <- before_cells - after_cells
removed_samples <- before_samples - after_samples
if (any(md$GSM %chin% excluded_gsm)) stop("排除 GSM 后仍有目标细胞残留")

sample_summary <- md[, .(
  DatasetID = first(DatasetID),
  SampleID = first(SampleID),
  GSM = paste(sort(unique(na.omit(GSM))), collapse = ";"),
  Disease = first(Disease),
  n_cells = .N
), by = SampleID_Global]
setorder(sample_summary, DatasetID, Disease, SampleID_Global)
fwrite(sample_summary, file.path(output_root, "tables", "sample_labels_and_cell_counts.csv"))

dataset_disease <- sample_summary[, .(
  n_samples = uniqueN(SampleID_Global),
  n_cells = sum(n_cells)
), by = .(DatasetID, Disease)]
if (!all(dataset_disease$DatasetID %in% names(dataset_gse))) {
  stop("DatasetID 缺少原 Panel A 的 GSE 对照：",
       paste(setdiff(unique(dataset_disease$DatasetID), names(dataset_gse)), collapse = ", "))
}
dataset_disease[, GSE := unname(dataset_gse[DatasetID])]
setcolorder(dataset_disease, c("DatasetID", "GSE", "Disease", "n_samples", "n_cells"))
setorder(dataset_disease, DatasetID, Disease)
fwrite(dataset_disease, file.path(output_root, "tables", "dataset_disease_summary.csv"))

disease_summary <- sample_summary[, .(
  n_samples = uniqueN(SampleID_Global),
  n_cells = sum(n_cells)
), by = Disease]
setorder(disease_summary, Disease)
fwrite(disease_summary, file.path(output_root, "tables", "disease_summary.csv"))

expected_gse202371 <- data.table(
  Disease = c("GBM", "LGG", "LCBM"),
  expected_n_samples = c(2L, 2L, 10L)
)
actual_gse202371 <- dataset_disease[
  DatasetID == "LUNG_BRAIN_MET_2022", .(Disease, actual_n_samples = n_samples)
]
gse202371_audit <- merge(expected_gse202371, actual_gse202371,
                          by = "Disease", all = TRUE)
gse202371_audit[is.na(actual_n_samples), actual_n_samples := 0L]
gse202371_audit[, check_status := fifelse(
  !is.na(expected_n_samples) & expected_n_samples == actual_n_samples,
  "PASS", "MISMATCH")]
fwrite(gse202371_audit, file.path(output_root, "tables", "GSE202371_sample_count_check.csv"))

actual_diseases <- sort(unique(sample_summary$Disease))
actual_datasets <- sort(unique(sample_summary$DatasetID))
if (!all(actual_datasets %in% names(dataset_gse))) stop("实际 DatasetID 未全部映射到 GSE")
dpal <- fread(palette_path)
needed_palette_cols <- c("Disease_Abbrev", "Disease_Color_Hex")
if (!all(needed_palette_cols %in% names(dpal))) stop("冻结疾病配色表缺少必需字段")
disease_colors <- setNames(dpal$Disease_Color_Hex, dpal$Disease_Abbrev)
disease_colors["Ctrl"] <- "#4DAF4A"
disease_colors["GBM"] <- "#B6362D"
disease_colors["LGG"] <- "#8C1D18"
missing_colors <- setdiff(actual_diseases, names(disease_colors))
if (length(missing_colors)) stop("实际疾病标签缺少冻结颜色：", paste(missing_colors, collapse = ", "))

coverage <- dataset_disease[, .(DatasetID, GSE, Disease, n_samples, n_cells)]
sample_matrix <- dcast(coverage, DatasetID ~ Disease, value.var = "n_samples", fill = 0)
cell_matrix <- dcast(coverage, DatasetID ~ Disease, value.var = "n_cells", fill = 0)
dataset_order <- actual_datasets
disease_order <- actual_diseases
sample_matrix <- sample_matrix[match(dataset_order, DatasetID)]
cell_matrix <- cell_matrix[match(dataset_order, DatasetID)]
sample_mat <- as.matrix(sample_matrix[, ..disease_order])
cell_mat <- as.matrix(cell_matrix[, ..disease_order])
storage.mode(sample_mat) <- "integer"
storage.mode(cell_mat) <- "integer"
rownames(sample_mat) <- dataset_order
rownames(cell_mat) <- dataset_order
colnames(sample_mat) <- disease_order
colnames(cell_mat) <- disease_order

source_long <- melt(coverage, id.vars = c("DatasetID", "GSE", "Disease"),
                    variable.name = "measure", value.name = "value")
source_long[, statistical_unit := "SampleID_Global"]
source_long[, note := sample_note]
fwrite(source_long, file.path(output_root, "source_data", "panelA_source.csv"))

disease_table_md <- c(
  "| Disease | SampleID_Global 数 | 细胞数 |",
  "|---|---:|---:|",
  vapply(seq_len(nrow(disease_summary)), function(i) {
    sprintf("| %s | %s | %s |", disease_summary$Disease[i],
            format(disease_summary$n_samples[i], big.mark = ","),
            format(disease_summary$n_cells[i], big.mark = ","))
  }, character(1))
)

# Marginal totals are calculated from the same post-exclusion object summaries.
disease_cells <- colSums(cell_mat)
disease_samples <- colSums(sample_mat)
dataset_cells <- rowSums(cell_mat)
dataset_samples <- rowSums(sample_mat)
historical_crosscheck <- data.table(
  metric = c("samples_before", "cells_before", "excluded_samples", "excluded_cells",
             "samples_after", "cells_after"),
  expected_for_comparison = c(453L, 424222L, 2L, 8982L, 451L, 415240L),
  observed = c(before_samples, before_cells, removed_samples, removed_cells,
              after_samples, after_cells)
)
historical_crosscheck[, comparison_status := fifelse(
  expected_for_comparison == observed, "MATCH", "DIFFERS; inspect metadata audit")]
fwrite(historical_crosscheck,
       file.path(output_root, "tables", "historical_count_crosscheck.csv"))

format_cell_count <- function(x) {
  ifelse(x >= 1e6,
         paste0(formatC(x / 1e6, format = "f", digits = 1), "M"),
         ifelse(x >= 1e3,
                paste0(formatC(x / 1e3, format = "f", digits = 0), "k"),
                formatC(x, format = "f", digits = 0, big.mark = ",")))
}
make_log10_cell_axis <- function(values, labels_rot = 0) {
  max_value <- max(values, na.rm = TRUE)
  max_exponent <- floor(log10(max_value))
  exponents <- unique(c(0L, seq(2L, max_exponent, by = 2L), max_exponent))
  tick_values <- 10^exponents
  list(at = log10(tick_values), labels = format_cell_count(tick_values),
       labels_rot = labels_rot, gp = gpar(fontsize = 5.5, fontfamily = "Arial"))
}
cell_luminance <- function(hex) {
  rgb <- grDevices::col2rgb(hex) / 255
  linear <- ifelse(rgb <= 0.03928, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
  0.2126 * linear[1] + 0.7152 * linear[2] + 0.0722 * linear[3]
}
cell_fun <- function(j, i, x, y, width, height, fill) {
  n <- sample_mat[i, j]
  disease <- colnames(sample_mat)[j]
  if (n == 0L) {
    grid.rect(x, y, width, height, gp = gpar(fill = "#F2F2F2", col = "white"))
  } else {
    clr <- unname(disease_colors[disease])
    grid.rect(x, y, width, height, gp = gpar(fill = clr, col = "white"))
    grid.text(as.character(n), x, y,
              gp = gpar(fontsize = 7.5, fontface = "bold", fontfamily = "Arial",
                        col = if (cell_luminance(clr) > 0.5) "black" else "white"))
  }
}

top_anno <- HeatmapAnnotation(
  Cells = anno_barplot(log10(pmax(disease_cells, 1)),
                       gp = gpar(fill = "#2C7FB8"), border = TRUE,
                       ylim = c(0, max(log10(pmax(disease_cells, 1))) * 1.15),
                       axis_param = make_log10_cell_axis(disease_cells, labels_rot = 0)),
  Samples = anno_barplot(disease_samples, gp = gpar(fill = "#238B45"), border = TRUE,
                         ylim = c(0, max(disease_samples) * 1.15)),
  annotation_height = unit(c(1.25, 1.1), "cm"),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 7, fontfamily = "Arial")
)
right_anno <- rowAnnotation(
  Cells = anno_barplot(log10(pmax(dataset_cells, 1)),
                       gp = gpar(fill = "#2C7FB8"), border = TRUE,
                       ylim = c(0, max(log10(pmax(dataset_cells, 1))) * 1.15),
                       axis_param = make_log10_cell_axis(dataset_cells, labels_rot = 90)),
  Samples = anno_barplot(dataset_samples, gp = gpar(fill = "#238B45"), border = TRUE,
                         ylim = c(0, max(dataset_samples) * 1.15)),
  annotation_width = unit(c(1.9, 1.6), "cm"),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 7, fontfamily = "Arial")
)
legend <- Legend(
  title = "Disease",
  at = disease_order,
  legend_gp = gpar(fill = disease_colors[disease_order]),
  ncol = 1,
  grid_width = unit(3.8, "mm"),
  grid_height = unit(3.8, "mm"),
  title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Arial"),
  labels_gp = gpar(fontsize = 7, fontfamily = "Arial")
)
row_labels <- unname(dataset_gse[dataset_order])
heatmap <- Heatmap(
  cell_mat,
  name = "cells",
  col = colorRamp2(c(0, 1), c("white", "white")),
  width = unit(ncol(cell_mat) * 7.5, "mm"),
  height = unit(nrow(cell_mat) * 7.5, "mm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cell_fun = cell_fun,
  rect_gp = gpar(col = "grey90", lwd = 0.4),
  row_labels = row_labels,
  row_names_side = "left",
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 8, fontfamily = "Arial", col = "black"),
  column_names_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Arial", col = "black"),
  column_names_rot = 45,
  top_annotation = top_anno,
  right_annotation = right_anno,
  show_heatmap_legend = FALSE
)
figure_title <- paste0(
  "Dataset × Disease coverage (", format(sum(cell_mat), big.mark = ","),
  " cells; tile numbers = samples)"
)
draw_panel <- function() {
  draw(heatmap, column_title = figure_title,
       column_title_gp = gpar(fontsize = 8.5, fontface = "bold", fontfamily = "Arial"),
       annotation_legend_list = list(legend),
       padding = unit(c(2, 12, 12, 2), "mm"))
  upViewport(0)
  grid.text(sample_note, x = unit(0.5, "npc"), y = unit(3, "mm"),
            just = "bottom",
            gp = gpar(fontsize = 7, fontfamily = "Microsoft YaHei", col = "black"))
}
png_path <- file.path(output_root, "figures", "Fig1A_dataset_disease_all_samples_excluding_mixed_GSM.png")
pdf_path <- file.path(output_root, "figures", "Fig1A_dataset_disease_all_samples_excluding_mixed_GSM.pdf")
png(png_path, width = 12, height = 8.2, units = "in", res = 600, type = "cairo", bg = "white")
draw_panel()
dev.off()
cairo_pdf(pdf_path, width = 12, height = 8.2, family = "Arial", bg = "white")
draw_panel()
dev.off()

metrics <- data.table(
  metric = c("cells_before_exclusion", "samples_before_exclusion",
             "cells_excluded", "samples_excluded",
             "cells_after_exclusion", "samples_after_exclusion",
             "n_disease_labels", "n_datasets"),
  value = c(before_cells, before_samples, removed_cells, removed_samples,
            after_cells, after_samples, length(actual_diseases), length(actual_datasets))
)
fwrite(metrics, file.path(output_root, "tables", "filtering_and_scope_summary.csv"))

object_summary <- data.table(
  input_path = atlas_path,
  default_assay = default_assay,
  assays = paste(assay_names, collapse = ";"),
  RNA_layers = paste(rna_layers, collapse = ";"),
  metadata_source = "Seurat meta.data",
  sample_unit = "SampleID_Global",
  exclusion_GSM = paste(excluded_gsm, collapse = ";"),
  cells_before = before_cells,
  samples_before = before_samples,
  cells_excluded = removed_cells,
  samples_excluded = removed_samples,
  cells_after = after_cells,
  samples_after = after_samples,
  diseases_after = paste(actual_diseases, collapse = ";"),
  n_datasets_after = length(actual_datasets)
)
fwrite(object_summary, file.path(output_root, "tables", "input_and_method_audit.csv"))

readme <- c(
  "# Fig1A：全部对象样本，排除两个混合GSM",
  "",
  "## 分析目的",
  "",
  "展示当前 Seurat 对象中各 DatasetID × Disease 组合的样本覆盖与细胞覆盖；样本按 SampleID_Global 计数，不代表独立供体数。",
  "",
  "## 输入与过滤",
  "",
  paste0("- 输入对象：", atlas_path),
  paste0("- Seurat assays：", paste(assay_names, collapse = ", "),
         "；DefaultAssay：", default_assay,
         "；RNA layers：", paste(rna_layers, collapse = ", ")),
  paste0("- 过滤前：", before_samples, " 个 SampleID_Global，", format(before_cells, big.mark = ","), " 个细胞。"),
  paste0("- 排除 GSM：", paste(excluded_gsm, collapse = ", "),
         "；移除 ", removed_samples, " 个样本、", format(removed_cells, big.mark = ","), " 个细胞。"),
  paste0("- 过滤后：", after_samples, " 个样本，", format(after_cells, big.mark = ","), " 个细胞。"),
  paste0("- Disease 标签（来自过滤后对象）：", paste(actual_diseases, collapse = ", "), "。"),
  "",
  "## 方法",
  "",
  "直接读取 Seurat meta.data。核对每个 SampleID_Global 唯一对应 DatasetID 和 Disease；精确匹配 GSM 字段排除指定混合文库。未设置每样本最少细胞数、供体数或疾病组门槛。按 DatasetID × Disease 汇总 n_samples 与 n_cells，图中方格数字为样本数，顶部及右侧边际柱显示疾病和数据集的细胞数、样本数。",
  "",
  "图件沿用原 Panel A 疾病填色和边际统计布局；冻结疾病颜色来自权威色表，疾病轴文字黑色，Arial，白底。标题说明纳入全部当前对象样本、排除两个混合 GSM，并按 SampleID_Global 计数。",
  "",
  "## 关键结果",
  "",
  paste0("- 过滤后共 ", length(actual_diseases), " 个 Disease 标签、", after_samples,
         " 个样本和 ", format(after_cells, big.mark = ","), " 个细胞。"),
  "",
  disease_table_md,
  "",
  paste0("- GSE202371 核对：",
         paste(dataset_disease[DatasetID == "LUNG_BRAIN_MET_2022",
                               paste0(Disease, "=", n_samples)], collapse = ", "), "。"),
  paste0("- 历史数值交叉核对：", paste(historical_crosscheck$metric,
       historical_crosscheck$comparison_status, sep = "=", collapse = "; "), "。"),
  "",
  "DatasetID × Disease 和 Disease 汇总数字分别见 tables/dataset_disease_summary.csv 与 tables/disease_summary.csv。",
  "",
  "## 生物学解释与限制",
  "",
  "Panel A 是对象内细胞和样本覆盖描述，用于识别疾病在各数据集中的组成；不同数据集的样本数和细胞数不等，图中不代表疾病间统计检验，也不能把样本数解释为独立供体数。两个排除 GSM 在原对象中的疾病标签为混合文库的暂存标签，按用户指定全部排除；没有进行 HTO 去卷积，也没有修改原对象。",
  "",
  "## 输出",
  "",
  "- figures/Fig1A_dataset_disease_all_samples_excluding_mixed_GSM.png（600 dpi）",
  "- figures/Fig1A_dataset_disease_all_samples_excluding_mixed_GSM.pdf（Arial 可编辑文字）",
  "- source_data/panelA_source.csv",
  "- tables/sample_labels_and_cell_counts.csv",
  "- tables/dataset_disease_summary.csv",
  "- tables/disease_summary.csv",
  "- tables/excluded_gsm_audit.csv",
  "- tables/filtering_and_scope_summary.csv",
  "- tables/input_and_method_audit.csv"
)
writeLines(enc2utf8(readme), file.path(output_root, "解读.md"), useBytes = TRUE)
log_msg("完成：", after_samples, " samples / ", after_cells, " cells; diseases=", paste(actual_diseases, collapse = ", "))

