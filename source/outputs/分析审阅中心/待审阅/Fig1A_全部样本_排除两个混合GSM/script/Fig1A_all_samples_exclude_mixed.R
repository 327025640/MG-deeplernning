# 独立重绘 Fig1A：当前对象全部样本，精确排除两个混合疾病 GSM。
# 仅统计元数据；不修改输入对象，不运行 Panel B–F。
invisible(Sys.setlocale("LC_CTYPE", "Chinese_China.utf8"))
options(stringsAsFactors = FALSE)
required_packages <- c("SeuratObject", "data.table", "ComplexHeatmap", "circlize", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                             quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages)) stop("缺少依赖：", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(SeuratObject)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

ATLAS <- "E:/micro_output/v3_dims16/objects/merged_obj.rds"
INPUT_AUDIT <- "E:/micro_output/v3_dims16/objects/merged_obj_audit.csv"
PALETTE <- "D:/data/microglia_new/authority/03_reference_configs/disease_palette_contract.csv"
ROOT <- "E:/micro_output/分析审阅中心/待审阅/Fig1A_全部样本_排除两个混合GSM"
EXCLUDE_GSM <- c("GSM6186080", "GSM6186087")
FIG <- file.path(ROOT, "figures")
SRC <- file.path(ROOT, "source_data")
AUD <- file.path(ROOT, "audit")
LOG <- file.path(ROOT, "logs")
for (p in c(ATLAS, INPUT_AUDIT, PALETTE)) {
  if (!file.exists(p)) stop("输入文件不存在：", p)
}
for (p in c(ROOT, FIG, SRC, AUD, LOG)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), paste0(..., collapse = ""))
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = file.path(LOG, "run_log.txt"), append = TRUE)
}
DATASET_GSE <- c(
  ABDELFATTAH_GBM_2022 = "GSE182109", BOSTRAND_HD_GLIA_2024 = "GSE281069",
  GSE103224 = "GSE103224", GSE141383 = "GSE141383", KISS_2023 = "GSE227781",
  LIM_2022 = "GSE180928", LIU_SEOAD_2025 = "GSE272082", LUNG_BRAIN_MET_2022 = "GSE202371",
  MARTIROSYAN_2024 = "GSE243639", MEDULLOBLASTOMA_2021 = "GSE155446",
  SMITH_2022 = "GSE160936", TAGORE_2025_NSCLC_BM = "GSE223499",
  TUDDENHAM_2024 = "GSE204702", WANG_C9ORF72_ALS_FTD_2024 = "GSE219280",
  WANG_C9_DLPFC_PNAS_2025 = "GSE212630", ZHU_2024 = "GSE202210"
)

# 输入指纹同时用于缓存验证，防止复用其他对象的元数据。
log_msg("核对输入指纹")
input_sha <- digest::digest(ATLAS, algo = "sha256", file = TRUE)
input_audit <- fread(INPUT_AUDIT)
stopifnot(all(c("item", "value") %in% names(input_audit)))
recorded_sha <- input_audit[item == "sha256", value]
if (length(recorded_sha) != 1L || !identical(tolower(input_sha), tolower(recorded_sha))) {
  stop("当前对象指纹与其审计文件不一致")
}
cache_file <- file.path(AUD, "metadata_cache.rds")
use_cache <- "--metadata-cache" %in% commandArgs(trailingOnly = TRUE)
if (use_cache) {
  if (!file.exists(cache_file)) stop("未找到本分析的元数据缓存")
  cache <- readRDS(cache_file)
  if (!identical(cache$input_sha, input_sha)) stop("缓存对象指纹不一致")
  md <- copy(cache$metadata)
  object_info <- copy(cache$object_info)
  rm(cache)
  log_msg("使用指纹一致的当前对象元数据缓存")
} else {
  log_msg("直接读取 Seurat 对象")
  obj <- readRDS(ATLAS)
  if (!inherits(obj, "Seurat")) stop("输入不是 Seurat 对象")
  assay_names <- Assays(obj)
  if (!"RNA" %in% assay_names) stop("输入对象没有 RNA assay")
  layer_names <- Layers(obj[["RNA"]])
  log_msg("assays=", paste(assay_names, collapse = ";"),
          "；RNA layers=", paste(layer_names, collapse = ";"))
  meta_columns <- names(obj[[]])
  required_meta <- c("DatasetID", "SampleID", "SampleID_Global", "Disease")
  if (!all(required_meta %in% meta_columns)) {
    stop("缺少元数据字段：", paste(setdiff(required_meta, meta_columns), collapse = ", "))
  }
  if (!identical(rownames(obj[[]]), colnames(obj))) stop("细胞名与 metadata 行顺序不一致")
  md <- as.data.table(obj[[]][, intersect(c(required_meta, "GSM"), meta_columns), drop = FALSE])
  md[, CellID := colnames(obj)]
  object_info <- data.table(
    input_path = ATLAS, input_sha256 = input_sha, n_cells = ncol(obj), n_genes = nrow(obj),
    assays = paste(assay_names, collapse = ";"), rna_layers = paste(layer_names, collapse = ";"),
    metadata_columns = paste(meta_columns, collapse = ";"),
    metadata_extracted_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  rm(obj)
  invisible(gc())
  saveRDS(list(input_sha = input_sha, metadata = md, object_info = object_info), cache_file)
}
fwrite(object_info, file.path(AUD, "object_input_audit.csv"))
required_meta <- c("CellID", "DatasetID", "SampleID", "SampleID_Global", "Disease")
for (nm in required_meta) {
  set(md, j = nm, value = as.character(md[[nm]]))
  bad <- is.na(md[[nm]]) | !nzchar(trimws(md[[nm]]))
  if (any(bad)) {
    fwrite(md[bad], file.path(AUD, paste0("invalid_", nm, ".csv")))
    stop(nm, " 存在缺失或空值，详见 audit")
  }
}
if (anyDuplicated(md$CellID)) stop("CellID 不唯一")
sample_check <- md[, .(n_dataset = uniqueN(DatasetID), n_disease = uniqueN(Disease),
                       n_sample_id = uniqueN(SampleID)), by = SampleID_Global]
fwrite(sample_check, file.path(AUD, "sample_label_uniqueness.csv"))
if (sample_check[n_dataset != 1L | n_disease != 1L | n_sample_id != 1L, .N] > 0L) {
  stop("同一样本存在数据集、疾病或原始样本编号冲突，详见 sample_label_uniqueness.csv")
}

# 在汇总之前排除；两个字段存在时交叉核对，而非仅匹配名字中的子串。
md[, Exclusion_GSM := fifelse(SampleID %in% EXCLUDE_GSM, SampleID, NA_character_)]
if ("GSM" %in% names(md)) {
  md[, GSM := as.character(GSM)]
  conflict <- md[SampleID %in% EXCLUDE_GSM & !is.na(GSM) & nzchar(GSM) & GSM != SampleID]
  if (nrow(conflict)) {
    fwrite(conflict, file.path(AUD, "exclusion_identifier_conflict.csv"))
    stop("目标样本的 GSM 与 SampleID 不一致")
  }
  md[GSM %in% EXCLUDE_GSM, Exclusion_GSM := GSM]
}
md[, Excluded := !is.na(Exclusion_GSM)]
if (!setequal(md[Excluded == TRUE, unique(Exclusion_GSM)], EXCLUDE_GSM)) {
  stop("没有完整匹配到两个指定排除 GSM")
}
if (md[, any(uniqueN(Excluded) > 1L), by = SampleID_Global]$V1 |> any()) {
  stop("排除规则只匹配到某个样本的部分细胞")
}
sample_all <- md[, .(
  DatasetID = first(DatasetID), SampleID = first(SampleID), Disease = first(Disease),
  GSM = if ("GSM" %in% names(md)) paste(sort(unique(na.omit(GSM))), collapse = ";") else NA_character_,
  n_cells = .N, Excluded = first(Excluded), Exclusion_GSM = first(Exclusion_GSM)
), by = SampleID_Global]
sample_all[, GSE := unname(DATASET_GSE[DatasetID])]
if (anyNA(sample_all$GSE)) stop("存在未覆盖的数据集到 GSE 映射")
sample_all[, Exclusion_reason := fifelse(Excluded, "混合疾病文库，按本次指令排除", "")]
setorder(sample_all, DatasetID, Disease, SampleID_Global)
included <- sample_all[Excluded == FALSE]
excluded <- sample_all[Excluded == TRUE]
fwrite(sample_all, file.path(SRC, "sample_metadata_all_with_exclusion.csv"), bom = TRUE)
fwrite(included, file.path(SRC, "sample_metadata_included.csv"), bom = TRUE)
fwrite(excluded, file.path(AUD, "excluded_samples.csv"), bom = TRUE)
fwrite(md[Excluded == TRUE, .(CellID, DatasetID, SampleID, SampleID_Global, Disease, Exclusion_GSM)],
       file.path(AUD, "excluded_cells.csv"))

summary_counts <- data.table(
  Stage = c("before_exclusion", "excluded", "included"),
  n_samples = c(nrow(sample_all), nrow(excluded), nrow(included)),
  n_cells = c(sum(sample_all$n_cells), sum(excluded$n_cells), sum(included$n_cells))
)
fwrite(summary_counts, file.path(AUD, "before_after_counts.csv"))
crosscheck <- data.table(
  Metric = c("before_samples", "before_cells", "GSM6186080_cells", "GSM6186087_cells",
             "included_samples", "included_cells"),
  Previously_reported = c(453L, 424222L, 3072L, 5910L, 451L, 415240L),
  Actual = c(nrow(sample_all), sum(sample_all$n_cells),
             excluded[Exclusion_GSM == "GSM6186080", sum(n_cells)],
             excluded[Exclusion_GSM == "GSM6186087", sum(n_cells)],
             nrow(included), sum(included$n_cells))
)
crosscheck[, Matches := Previously_reported == Actual]
fwrite(crosscheck, file.path(AUD, "previous_count_crosscheck.csv"))
if (any(!crosscheck$Matches)) warning("与此前核查数字不同，详见 previous_count_crosscheck.csv")

coverage <- included[, .(n_samples = .N, n_cells = sum(n_cells)), by = .(DatasetID, Disease)]
ds_all <- sort(unique(included$DatasetID))
dis_all <- sort(unique(included$Disease))
coverage <- merge(CJ(DatasetID = ds_all, Disease = dis_all), coverage,
                  by = c("DatasetID", "Disease"), all.x = TRUE, sort = TRUE)
coverage[is.na(n_samples), `:=`(n_samples = 0L, n_cells = 0L)]
coverage[, `:=`(GSE = unname(DATASET_GSE[DatasetID]), Statistical_Unit = "SampleID_Global",
               Sample_scope = "All samples present in object after GSM exclusion",
               Minimum_cells_present = 1L, Excluded_GSM = paste(EXCLUDE_GSM, collapse = ";"))]
fwrite(coverage, file.path(SRC, "panelA_coverage_data.csv"))
disease_summary <- included[, .(n_samples = .N, n_cells = sum(n_cells),
                                n_datasets = uniqueN(DatasetID)), by = Disease][order(Disease)]
dataset_summary <- included[, .(n_samples = .N, n_cells = sum(n_cells),
                                n_diseases = uniqueN(Disease)), by = .(DatasetID, GSE)][order(DatasetID)]
fwrite(disease_summary, file.path(SRC, "disease_summary.csv"))
fwrite(dataset_summary, file.path(SRC, "dataset_summary.csv"))
gse202371 <- coverage[GSE == "GSE202371" & n_samples > 0, .(GSE, Disease, n_samples, n_cells)]
fwrite(gse202371, file.path(AUD, "GSE202371_current_labels.csv"))
gse_expected <- c(GBM = 2L, LGG = 2L, LCBM = 10L)
gse_ok <- setequal(gse202371$Disease, names(gse_expected)) &&
  all(gse202371$n_samples == unname(gse_expected[gse202371$Disease]))
if (!gse_ok) warning("GSE202371 与此前核查不一致，详见 GSE202371_current_labels.csv")
stopifnot(nrow(included) == uniqueN(md[Excluded == FALSE, SampleID_Global]),
          sum(coverage$n_samples) == nrow(included),
          sum(coverage$n_cells) == md[Excluded == FALSE, .N],
          sum(disease_summary$n_cells) == sum(included$n_cells),
          sum(dataset_summary$n_samples) == nrow(included))
log_msg("排除后：", nrow(included), " 个样本，", sum(included$n_cells), " 个细胞；",
        length(ds_all), " 个数据集，", length(dis_all), " 个疾病/对照标签")

# 配色与旧 Panel A 相同；类别列表始终来自当前纳入元数据。
dpal <- fread(PALETTE)
stopifnot(all(c("Disease_Abbrev", "Disease_Color_Hex") %in% names(dpal)))
dis_color <- setNames(dpal$Disease_Color_Hex, dpal$Disease_Abbrev)
dis_color[c("Ctrl", "GBM", "LGG")] <- c("#4DAF4A", "#B6362D", "#8C1D18")
if (length(setdiff(dis_all, names(dis_color)))) stop("当前疾病标签缺少冻结颜色")
if (anyDuplicated(unname(dis_color[dis_all]))) stop("疾病颜色存在重复")
fwrite(data.table(Disease = dis_all, Color_Hex = unname(dis_color[dis_all])),
       file.path(AUD, "disease_palette_used.csv"))
cellm <- sampm <- matrix(0, nrow = length(ds_all), ncol = length(dis_all),
                        dimnames = list(ds_all, dis_all))
matrix_index <- cbind(match(coverage$DatasetID, ds_all), match(coverage$Disease, dis_all))
cellm[matrix_index] <- coverage$n_cells
sampm[matrix_index] <- coverage$n_samples
stopifnot(sum(cellm) == sum(included$n_cells), sum(sampm) == nrow(included))
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
format_cell_count <- function(x) {
  ifelse(x >= 1e6, paste0(format(x / 1e6, trim = TRUE, scientific = FALSE), "M"),
         ifelse(x >= 1e3, paste0(format(x / 1e3, trim = TRUE, scientific = FALSE), "k"),
                format(x, trim = TRUE, scientific = FALSE)))
}
make_log10_cell_axis <- function(values) {
  ticks <- 10^(0:floor(log10(max(values))))
  list(at = log10(ticks), labels = format_cell_count(ticks),
       gp = gpar(fontfamily = "Arial", fontsize = 6))
}
dis_cells <- colSums(cellm)
dis_samples <- colSums(sampm)
ds_cells <- rowSums(cellm)
ds_samples <- rowSums(sampm)
# 在构建注释前开启 Cairo 设备，Arial 的文本尺寸计算与实际输出使用相同字体。
log_msg("输出 PNG 与 PDF")
png(file.path(FIG, "PanelA_dataset_disease_coverage.png"),
    width = 12, height = 8.2, units = "in", res = 600, type = "cairo", bg = "white")
ha_top <- HeatmapAnnotation(
  Cells = anno_barplot(log10(pmax(dis_cells, 1)), gp = gpar(fill = "#2C7FB8"), border = TRUE,
                      ylim = c(0, max(log10(pmax(dis_cells, 1))) * 1.15),
                      axis_param = make_log10_cell_axis(dis_cells)),
  Samples = anno_barplot(dis_samples, gp = gpar(fill = "#238B45"), border = TRUE,
                        ylim = c(0, max(dis_samples) * 1.15),
                        axis_param = list(gp = gpar(fontfamily = "Arial", fontsize = 6))),
  annotation_height = unit(c(1.25, 1.1), "cm"), show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 7, fontfamily = "Arial")
)
ha_right <- rowAnnotation(
  Cells = anno_barplot(log10(pmax(ds_cells, 1)), gp = gpar(fill = "#2C7FB8"), border = TRUE,
                      ylim = c(0, max(log10(pmax(ds_cells, 1))) * 1.15),
                      axis_param = make_log10_cell_axis(ds_cells)),
  Samples = anno_barplot(ds_samples, gp = gpar(fill = "#238B45"), border = TRUE,
                        axis_param = list(gp = gpar(fontfamily = "Arial", fontsize = 6))),
  annotation_width = unit(c(1.9, 1.6), "cm"), show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 7, fontfamily = "Arial")
)
legend_A <- Legend(title = "Disease", at = dis_all, legend_gp = gpar(fill = dis_color[dis_all]),
                   ncol = 1, grid_width = unit(3.8, "mm"), grid_height = unit(3.8, "mm"),
                   title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Arial"),
                   labels_gp = gpar(fontsize = 7, fontfamily = "Arial"))
htA <- Heatmap(
  cellm, name = "coverage", col = colorRamp2(c(0, 1), c("white", "white")),
  width = unit(ncol(cellm) * 7.5, "mm"), height = unit(nrow(cellm) * 7.5, "mm"),
  cluster_rows = FALSE, cluster_columns = FALSE, cell_fun = panelA_cell_fun,
  rect_gp = gpar(col = "grey90", lwd = 0.4), row_labels = unname(DATASET_GSE[ds_all]),
  row_names_side = "left", column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 8, fontfamily = "Arial"),
  column_names_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Arial", col = "black"),
  column_names_rot = 45, top_annotation = ha_top, right_annotation = ha_right,
  show_heatmap_legend = FALSE
)
title_A <- sprintf("Dataset x Disease coverage (%s cells; %s samples)",
                   format(sum(cellm), big.mark = ","), format(sum(sampm), big.mark = ","))
caption <- paste("All samples in the current Seurat object; unit: SampleID_Global; tile numbers: samples.",
                 "Excluded mixed-disease libraries: GSM6186080 and GSM6186087.", sep = "\n")
draw_panel <- function() {
  draw(htA, column_title = title_A,
       column_title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial"),
       annotation_legend_list = list(legend_A), padding = unit(c(14, 12, 2, 2), "mm"))
  grid.text(caption, x = unit(0.5, "npc"), y = unit(5, "mm"),
            gp = gpar(fontfamily = "Arial", fontsize = 7, col = "black"))
}
draw_panel()
invisible(dev.off())
cairo_pdf(file.path(FIG, "PanelA_dataset_disease_coverage.pdf"),
          width = 12, height = 8.2, family = "Arial", bg = "white")
draw_panel()
invisible(dev.off())

fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
disease_lines <- sprintf("| %s | %s | %s | %s |", disease_summary$Disease,
                         fmt(disease_summary$n_samples), fmt(disease_summary$n_cells),
                         disease_summary$n_datasets)
report <- c(
  "# Fig1A：全部样本，排除两个混合疾病 GSM", "",
  "## 分析目的", "",
  "按当前 Seurat 对象重新展示数据集与疾病覆盖，纳入所有对象内样本，并排除尚未拆分疾病来源的两个混合文库。", "",
  "## 输入与方法", "",
  paste0("- 对象：`", ATLAS, "`；SHA256：`", input_sha, "`。"),
  "- 直接读取 meta.data 的 DatasetID、SampleID、SampleID_Global、Disease 和 GSM（如存在）。",
  "- 样本统计单位为 SampleID_Global，保留对象内至少有一个细胞的所有样本，不设置额外细胞数、供体数或疾病组门槛。",
  "- 在样本和细胞汇总前，精确排除 GSM6186080 与 GSM6186087；未进行 HTO 去卷积。",
  "- 格内数字为样本数，色块为冻结疾病颜色；蓝色边际柱为 log10 细胞数（刻度显示原始数量），绿色边际柱为样本数。",
  "- 独立输出 600 DPI PNG、保留可编辑文字的 PDF 和源数据 CSV；所有疾病轴文字为黑色。", "",
  "## 结果", "",
  "| 阶段 | 样本数 | 细胞数 |", "|---|---:|---:|",
  sprintf("| 原对象 | %s | %s |", fmt(nrow(sample_all)), fmt(sum(sample_all$n_cells))),
  sprintf("| 排除 | %s | %s |", fmt(nrow(excluded)), fmt(sum(excluded$n_cells))),
  sprintf("| 新 Fig1A | %s | %s |", fmt(nrow(included)), fmt(sum(included$n_cells))), "",
  sprintf("新图覆盖 %s 个数据集、%s 个疾病/对照标签；纳入细胞数不足 20 的样本 %s 个（%s 个细胞）。",
          length(ds_all), length(dis_all), included[n_cells < 20, .N],
          included[n_cells < 20, sum(n_cells)]), "",
  "| Disease | 样本数 | 细胞数 | 数据集数 |", "|---|---:|---:|---:|", disease_lines, "",
  "## 核查", "",
  "- 样本的 DatasetID、Disease 与 SampleID 唯一性，以及缺失值、重复细胞名均经实际输入检查。",
  "- 格内、行列边际与总体样本数/细胞数加和一致；两个指定 GSM 的全部细胞均被排除。",
  paste0("- 此前计数交叉核查：", if (all(crosscheck$Matches)) "全部一致。" else "存在差异，详见 audit/previous_count_crosscheck.csv。"),
  paste0("- GSE202371 当前标签核查：", if (gse_ok) "GBM=2、LGG=2、LCBM=10，与当前对象一致。" else "与此前预期不同，详见 audit/GSE202371_current_labels.csv。"),
  "- 本次只生成独立 Panel A，原始 Seurat 对象与旧版图件未被覆盖。", "",
  "## 生物学意义与限制", "",
  "该图描述当前对象实际包含的疾病覆盖范围。排除混合疾病文库后，避免将其所有细胞统一归入 ALS_FTD 或 AD。",
  "样本数不是独立供体数；同一供体的多个组织或文库可能对应多个样本。低细胞数样本只用于覆盖描述，不意味着满足后续推断分析的纳入标准。",
  "本次未全量重建其他数据集的供体/脑区身份，也未改变其他下游分析的冻结输入。", "",
  "## 输出", "",
  "- figures/PanelA_dataset_disease_coverage.png 与 .pdf",
  "- source_data/sample_metadata_all_with_exclusion.csv、sample_metadata_included.csv",
  "- source_data/panelA_coverage_data.csv、disease_summary.csv、dataset_summary.csv",
  "- audit/：输入指纹、样本唯一性、排除清单及细胞列表、计数核查、配色和元数据缓存。"
)
writeLines(enc2utf8(report), file.path(ROOT, "解读.md"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(LOG, "sessionInfo.txt"))
log_msg("运行完成")
print(disease_summary)
