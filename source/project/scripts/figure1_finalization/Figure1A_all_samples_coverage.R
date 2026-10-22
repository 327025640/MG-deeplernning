# Windows UTF-8路径兼容，避免继承无效的C.UTF-8设置。
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
# Figure 1A：复用冻结逐细胞源表，统计全部样本；仅导出A图。
suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})
ROOT <- "E:/micro_output/分析审阅中心/待审阅/Figure1_图件收尾"
FIG <- file.path(ROOT, "figures")
SRC <- file.path(ROOT, "source_data")
INPUT <- "E:/micro_output/Figure1_Disease_V3_Disease17/source_data/panelB_umap_data.csv"
PALETTE_DISEASE <- "D:/data/microglia_new/authority/03_reference_configs/disease_palette_contract.csv"
stopifnot(file.exists(INPUT), file.exists(PALETTE_DISEASE))
# 使用base连接兼容Windows中文输出路径。
write_csv <- function(x, path) write.csv(as.data.frame(x), path, row.names = FALSE, fileEncoding = "UTF-8")
md <- fread(INPUT)
stopifnot(all(c("CellID", "SampleID_Global", "DatasetID", "Disease") %in% names(md)))
stopifnot(nrow(md) == 424222L, !anyDuplicated(md$CellID))
for (v in c("CellID", "SampleID_Global", "DatasetID", "Disease")) {
  stopifnot(!anyNA(md[[v]]), all(nzchar(md[[v]])))
}
keys <- md[, .(n_dataset = uniqueN(DatasetID), n_disease = uniqueN(Disease)), by = SampleID_Global]
stopifnot(all(keys$n_dataset == 1L), all(keys$n_disease == 1L))
samp <- md[, .(n_cells = .N), by = .(SampleID_Global, DatasetID, Disease)]
stopifnot(nrow(samp) == 453L, uniqueN(samp$DatasetID) == 16L, uniqueN(samp$Disease) == 17L)
actual_diseases <- sort(unique(md$Disease))
dpal <- fread(PALETTE_DISEASE)
stopifnot(all(c("Disease_Abbrev", "Disease_Color_Hex") %in% names(dpal)))
dis_color <- setNames(dpal$Disease_Color_Hex, dpal$Disease_Abbrev)
# 与现有A图完全一致的主疾病映射，不修改共享色表。
dis_color["Ctrl"] <- "#4DAF4A"
dis_color["GBM"] <- "#B6362D"
dis_color["LGG"] <- "#8C1D18"
stopifnot(all(actual_diseases %in% names(dis_color)), !anyNA(dis_color[actual_diseases]))
for (p in c(ROOT, FIG, SRC, file.path(ROOT, "logs"))) dir.create(p, recursive = TRUE, showWarnings = FALSE)
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
ds_all <- sort(unique(samp$DatasetID))
# 显式锁定原图中两个WANG数据集的行顺序，避免locale改变排序。
w <- which(ds_all %in% c("WANG_C9_DLPFC_PNAS_2025", "WANG_C9ORF72_ALS_FTD_2024"))
ds_all[w] <- c("WANG_C9_DLPFC_PNAS_2025", "WANG_C9ORF72_ALS_FTD_2024")
if (!all(ds_all %in% names(DATASET_GSE))) {
  stop("DATASET_GSE 未覆盖：", paste(setdiff(ds_all, names(DATASET_GSE)), collapse = ", "))
}
# 2026-09-28 用户裁决：热图列序 = 图例顺序，按疾病大类排列，类内按亚类聚。
# 疾病到"大类/亚类"的归属从权威色表契约读取，本脚本不重复维护；
# 大类的先后与亚类的先后由下面两个常量显式冻结。
# 仅 GBM 为契约外别名（契约写作 GBM_IDHwt），与上方颜色特例同源处理。
CLASS_ORDER <- c(
  "Control",
  "Neurodegenerative_disease",
  "Neuroinflammatory_immune_disease",
  "Primary_brain_tumor",
  "Brain_metastasis"
)
SUBCLASS_ORDER <- c(
  "Control",
  "Alzheimer_spectrum",
  "Synucleinopathy",
  "Motor_neuron_disease",
  "FTD_PSP_spectrum",
  "Huntington_disease",
  "Multiple_sclerosis",
  "Glioblastoma",
  "Diffuse_glioma",
  "Medulloblastoma",
  "Lung_cancer_brain_metastasis"
)
DISPLAY_ORDER <- c(
  "Ctrl",
  "AD", "EOAD",                        # Alzheimer_spectrum
  "PD", "DLBD_PD", "DLB",              # Synucleinopathy
  "ALS", "ALS_FTD", "FTD_MND",         # Motor_neuron_disease
  "FTD", "PSP",                        # FTD_PSP_spectrum
  "HD",                                # Huntington_disease
  "MS",                                # Multiple_sclerosis
  "GBM", "LGG", "MB",                  # Primary_brain_tumor
  "LCBM"                               # Brain_metastasis
)
dis_all <- DISPLAY_ORDER[DISPLAY_ORDER %in% actual_diseases]
stopifnot(
  length(dis_all) == length(actual_diseases),
  setequal(dis_all, actual_diseases),
  anyDuplicated(dis_all) == 0L
)
# 结构校验：展现顺序的大类分块与亚类分块必须与冻结顺序一致，防手工换序。
# GBM 别名只做新增；契约若已自带 GBM，必须是同一大类/亚类，避免别名静默覆盖契约。
stopifnot("GBM_IDHwt" %in% dpal$Disease_Abbrev)
class_of <- setNames(dpal$Disease_Class, dpal$Disease_Abbrev)
sub_of <- setNames(dpal$Disease_Subclass, dpal$Disease_Abbrev)
if ("GBM" %in% names(class_of)) {
  stopifnot(identical(class_of[["GBM"]], class_of[["GBM_IDHwt"]]),
            identical(sub_of[["GBM"]], sub_of[["GBM_IDHwt"]]))
} else {
  class_of["GBM"] <- class_of[["GBM_IDHwt"]]
  sub_of["GBM"] <- sub_of[["GBM_IDHwt"]]
}
dis_class <- unname(class_of[dis_all])
dis_subclass <- unname(sub_of[dis_all])
stopifnot(
  !anyNA(dis_class), !anyNA(dis_subclass),
  all(dis_class %in% CLASS_ORDER), all(dis_subclass %in% SUBCLASS_ORDER),
  all(diff(match(dis_class, CLASS_ORDER)) >= 0),
  all(diff(match(dis_subclass, SUBCLASS_ORDER)) >= 0)
)

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
# 只改刻度文字；保留原有log10坐标变换及刻度位置。
format_cell_count <- function(x) format(x, scientific = FALSE, big.mark = ",", trim = TRUE)
make_log10_cell_axis <- function(values) {
  max_value <- max(values, na.rm = TRUE)
  tick_values <- 10^(0:floor(log10(max_value)))
  list(at = log10(tick_values), labels = format_cell_count(tick_values), gp = gpar(fontsize = 7, fontfamily = "Arial"))
}
dis_cell_axis <- make_log10_cell_axis(dis_cells)
ds_cell_axis <- make_log10_cell_axis(ds_cells)
ha_top <- HeatmapAnnotation(
  `Cells (log10)` = anno_barplot(log10(pmax(dis_cells, 1)), gp = gpar(fill = "#2C7FB8"), border = TRUE,
                       ylim = c(0, max(log10(pmax(dis_cells, 1))) * 1.15),
                       axis_param = dis_cell_axis),
  Samples = anno_barplot(dis_samples, gp = gpar(fill = "#238B45"), border = TRUE,
                         ylim = c(0, max(dis_samples) * 1.15)),
  annotation_height = unit(c(2.0, 1.1), "cm"),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 7, fontfamily = "Arial")
)
ha_right <- rowAnnotation(
  `Cells (log10)` = anno_barplot(log10(pmax(ds_cells, 1)), gp = gpar(fill = "#2C7FB8"),
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
  "Dataset x Disease coverage (%s cells; numbers in tiles = all samples)",
  format(sum(cellm), big.mark = ",")
)
png(file.path(FIG, "Figure1A_all_samples_coverage.png"),
    width = 12, height = 8.2, units = "in", res = 600, type = "cairo")
draw(htA, column_title = title_A,
     column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial"),
     annotation_legend_list = list(legend_A), padding = unit(c(2, 12, 2, 2), "mm"))
dev.off()
cairo_pdf(file.path(FIG, "Figure1A_all_samples_coverage.pdf"),
          width = 12, height = 8.2, family = "Arial")
draw(htA, column_title = title_A,
     column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial"),
     annotation_legend_list = list(legend_A), padding = unit(c(2, 12, 2, 2), "mm"))
dev.off()
coverage[, `:=`(
  GSE = unname(DATASET_GSE[DatasetID]),
  Statistical_Unit = "all atlas samples; no minimum cell-count filter",
  Minimum_cells_per_sample = 0L
)]
stopifnot(sum(coverage$n_cells) == 424222L, sum(coverage$n_samples) == 453L)
write_csv(coverage, file.path(SRC, "Figure1A_all_samples_coverage_source.csv"))
write_csv(samp, file.path(SRC, "Figure1A_all_samples_sample_counts.csv"))
writeLines(capture.output(sessionInfo()), file.path(ROOT, "logs", "Figure1A_sessionInfo.txt"))
cat("COMPLETE: Figure1A; cells=424222; samples=453; datasets=16; diseases=17\n")
