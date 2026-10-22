# Figure 2 formal default: consensus 12 MP usage, res=0.6, Disease17.
suppressPackageStartupMessages({
  library(data.table); library(Seurat); library(ggplot2)
  library(ComplexHeatmap); library(circlize); library(grid)
})

OBJ <- "E:/micro_output/v3_dims16/objects/merged_obj.rds"
USAGE <- "E:/micro_output/v3_dims16/05_nmf/03_projection/primary_MP_v2_usage_dims16.rds"
ROOT <- "E:/micro_output/NMF"
FIG <- file.path(ROOT, "figures"); SRC <- file.path(ROOT, "source_data"); TBL <- file.path(ROOT, "tables")
MIN_CELLS <- 20L
for (p in c(OBJ, USAGE)) if (!file.exists(p)) stop("Missing input: ", p)
for (p in c(FIG, SRC, TBL)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))

# Archive only files replaced by this update; copy rather than delete.
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
archive_root <- file.path(ROOT, "archive", paste0("before_12MP_Disease17_usage_", stamp))
archive_files <- c(
  file.path(FIG, c(
    "PanelD_pseudobulk_MP_usage_heatmap.png", "PanelD_pseudobulk_MP_usage_heatmap.pdf",
    "PanelD_cluster_MP_usage.png", "PanelD_cluster_MP_usage.pdf",
    "PanelD_cluster_MP_usage_rowZ.png", "PanelD_cluster_MP_usage_rowZ.pdf",
    "PanelE_MP_disease_mean_usage.png", "PanelE_MP_disease_mean_usage.pdf",
    "PanelE_MP_disease_mean_usage_raw.png", "PanelE_MP_disease_mean_usage_raw.pdf"
  )),
  file.path(SRC, c(
    "PanelD_pseudobulk_MP_usage_raw_source.csv", "PanelD_pseudobulk_MP_usage_rowZ_source.csv",
    "PanelD_cluster_MP_usage_source.csv", "PanelD_cluster_MP_usage_rowZ_source.csv",
    "PanelE_MP_disease_mean_usage_source.csv", "PanelE_MP_disease_mean_usage_rowZ_source.csv"
  )),
  file.path(TBL, c(
    "PanelD_pseudobulk_metadata.csv", "PanelD_cluster_MP_usage_wide.csv",
    "PanelE_disease_MP_usage_wide.csv", "Figure2_12MP_usage_update_audit.csv"
  ))
)
archive_files <- archive_files[file.exists(archive_files)]
if (length(archive_files)) {
  dir.create(archive_root, recursive = TRUE, showWarnings = FALSE)
  if (!all(file.copy(archive_files, archive_root, overwrite = FALSE))) stop("Archive failed.")
  log_msg("Archived ", length(archive_files), " files -> ", archive_root)
}

disease_palette <- data.table(
  Disease = c("Ctrl", "LGG", "GBM", "MB", "LCBM", "AD", "EOAD", "PD", "DLB",
              "DLBD_PD", "ALS", "ALS_FTD", "FTD", "FTD_MND", "PSP", "HD", "MS"),
  Color_Hex = c("#4DAF4A", "#8C1D18", "#B6362D", "#D85A50", "#F39C72", "#073B77",
                "#3F88C5", "#007F8B", "#1B9AA3", "#4DBBC5", "#2F6FB0", "#92BCE0",
                "#5D4E9A", "#5F8FC9", "#A6A2D4", "#6B7C93", "#008C9E")
)
disease_order <- disease_palette$Disease
disease_colors <- setNames(disease_palette$Color_Hex, disease_palette$Disease)
fwrite(disease_palette, file.path(TBL, "Figure2_12MP_Disease17_palette.csv"))
ds_pal <- fread("D:/data/microglia_new/authority/03_reference_configs/dataset_palette_contract.csv")
dataset_colors <- setNames(ds_pal$Color_Hex, ds_pal$DatasetID)

log_msg("Loading frozen 12 MP usage ...")
u <- readRDS(USAGE); H <- u$H_norm; rm(u)
mp_order <- sprintf("MP%02d", 1:12)
stopifnot(is.matrix(H), identical(rownames(H), mp_order), !anyNA(H), all(is.finite(H)), all(H >= 0))
if (any(abs(colSums(H) - 1) > 1e-8)) stop("Per-cell usage does not sum to one.")

log_msg("Loading current merged_obj metadata ...")
obj <- readRDS(OBJ); md <- as.data.table(obj[[]]); md[, CellID := colnames(obj)]
need <- c("CellID", "SampleID_Global", "DatasetID", "Disease", "disease_analyse",
          "RNA_snn_res.0.6", "seurat_clusters")
if (length(setdiff(need, names(md)))) stop("Missing metadata: ", paste(setdiff(need, names(md)), collapse = ", "))
stopifnot(
  identical(as.character(md$`RNA_snn_res.0.6`), as.character(md$seurat_clusters)),
  identical(as.character(Idents(obj)), as.character(md$`RNA_snn_res.0.6`)),
  setequal(colnames(H), md$CellID), setequal(unique(md$Disease), disease_order)
)
rm(obj); invisible(gc())
md <- md[match(colnames(H), CellID)]; stopifnot(identical(md$CellID, colnames(H)))

dt <- as.data.table(t(H)); setnames(dt, mp_order)
dt[, `:=`(CellID = md$CellID, SampleID_Global = md$SampleID_Global, DatasetID = md$DatasetID,
          Disease = md$Disease, disease_analyse = md$disease_analyse,
          cluster_res0.6 = as.character(md$`RNA_snn_res.0.6`))]
setcolorder(dt, c("CellID", "SampleID_Global", "DatasetID", "Disease", "disease_analyse", "cluster_res0.6", mp_order))
rm(H, md); invisible(gc())

# Sample-level usage; disease labels must be unique within SampleID_Global.
meta_unique <- unique(dt[, .(SampleID_Global, DatasetID, Disease, disease_analyse)])
if (meta_unique[, anyDuplicated(SampleID_Global)] > 0L) stop("A sample maps to multiple metadata combinations.")
pb_all <- dt[, c(list(n_cells = .N), lapply(.SD, mean)),
             by = .(SampleID_Global, DatasetID, Disease, disease_analyse), .SDcols = mp_order]
pb_all[, included_ge20 := n_cells >= MIN_CELLS]
pb <- pb_all[included_ge20 == TRUE]
pb[, disease_rank := match(Disease, disease_order)]
setorder(pb, disease_rank, DatasetID, SampleID_Global); pb[, disease_rank := NULL]
pb_meta_out <- pb_all[, .(
  PseudobulkID = SampleID_Global, SampleID_Global, DatasetID, Disease,
  disease_analyse, Disease_Abbrev_Analysis = disease_analyse,
  n_cells, included_ge20
)]
fwrite(pb_meta_out, file.path(TBL, "PanelD_pseudobulk_metadata.csv"))
pb_source <- pb[, c("SampleID_Global", mp_order), with = FALSE]
setnames(pb_source, "SampleID_Global", "PseudobulkID")
fwrite(pb_source, file.path(SRC, "PanelD_pseudobulk_MP_usage_raw_source.csv"))

pb_raw <- as.matrix(pb[, ..mp_order]); rownames(pb_raw) <- pb$SampleID_Global
pb_z <- scale(pb_raw); pb_z[!is.finite(pb_z)] <- 0
fwrite(as.data.table(pb_z, keep.rownames = "PseudobulkID"),
       file.path(SRC, "PanelD_pseudobulk_MP_usage_rowZ_source.csv"))
missing_ds <- setdiff(unique(pb$DatasetID), names(dataset_colors))
if (length(missing_ds)) dataset_colors <- c(dataset_colors, setNames(hcl.colors(length(missing_ds), "Dark 3"), missing_ds))

ha_pb <- HeatmapAnnotation(
  Disease = pb$Disease, Dataset = pb$DatasetID,
  n_cells = anno_barplot(pb$n_cells, gp = gpar(fill = "grey45", col = NA), border = FALSE),
  col = list(Disease = disease_colors[unique(pb$Disease)], Dataset = dataset_colors[unique(pb$DatasetID)]),
  show_legend = c(TRUE, TRUE, FALSE),
  annotation_name_gp = gpar(fontfamily = "Arial", fontsize = 7, fontface = "bold")
)
ht_pb <- Heatmap(t(pb_z), name = "Row Z",
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
  row_names_gp = gpar(fontfamily = "Arial", fontsize = 8, fontface = "bold", col = "black"),
  top_annotation = ha_pb, use_raster = TRUE, raster_quality = 3, rect_gp = gpar(col = NA))
png(file.path(FIG, "PanelD_pseudobulk_MP_usage_heatmap.png"), 3600, 2500, res = 300, type = "cairo")
draw(ht_pb, heatmap_legend_side = "right", annotation_legend_side = "right"); dev.off()
cairo_pdf(file.path(FIG, "PanelD_pseudobulk_MP_usage_heatmap.pdf"), 12, 8.3, family = "Arial")
draw(ht_pb, heatmap_legend_side = "right", annotation_legend_side = "right"); dev.off()

# Cluster-level descriptive usage at res=0.6.
cl <- dt[, c(list(n_cells = .N), lapply(.SD, mean)), by = cluster_res0.6, .SDcols = mp_order]
cl[, ord := as.integer(cluster_res0.6)]; setorder(cl, ord); cluster_order <- cl$cluster_res0.6
cl_raw <- as.matrix(cl[, ..mp_order]); rownames(cl_raw) <- cluster_order
cl_z <- scale(cl_raw); cl_z[!is.finite(cl_z)] <- 0
cl_long <- melt(cl[, c("cluster_res0.6", "n_cells", mp_order), with = FALSE],
                id.vars = c("cluster_res0.6", "n_cells"), variable.name = "MP", value.name = "mean_usage")
cl_long[, `:=`(
  cluster_res0.6 = factor(cluster_res0.6, levels = cluster_order), MP = factor(MP, levels = mp_order),
  row_zscore = as.numeric(cl_z[cbind(match(as.character(cluster_res0.6), rownames(cl_z)),
                                     match(as.character(MP), colnames(cl_z)))])
)]
fwrite(cl[, ord := NULL], file.path(TBL, "PanelD_cluster_MP_usage_wide.csv"))
fwrite(cl_long, file.path(SRC, "PanelD_cluster_MP_usage_source.csv"))
fwrite(as.data.table(cl_z, keep.rownames = "cluster_res0.6"), file.path(SRC, "PanelD_cluster_MP_usage_rowZ_source.csv"))

theme_usage <- theme_minimal(base_family = "Arial", base_size = 9) +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black"), legend.text = element_text(colour = "black"))
p_cl_raw <- ggplot(cl_long, aes(cluster_res0.6, MP, fill = mean_usage)) + geom_tile() +
  scale_fill_gradient(low = "white", high = "#B2182B") +
  labs(x = "Seurat cluster (resolution 0.6)", y = NULL, fill = "Mean usage") + theme_usage
p_cl_z <- ggplot(cl_long, aes(cluster_res0.6, MP, fill = row_zscore)) + geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(x = "Seurat cluster (resolution 0.6)", y = NULL, fill = "Z-score") + theme_usage
save_gg <- function(p, stem) {
  ggsave(file.path(FIG, paste0(stem, ".png")), p, width = 8.4, height = 5.4, dpi = 600, bg = "white")
  ggsave(file.path(FIG, paste0(stem, ".pdf")), p, width = 8.4, height = 5.4, device = cairo_pdf, bg = "white")
}
save_gg(p_cl_raw, "PanelD_cluster_MP_usage"); save_gg(p_cl_z, "PanelD_cluster_MP_usage_rowZ")

# Disease-level equal-weight sample means.
pb_long <- melt(pb[, c("SampleID_Global", "DatasetID", "Disease", "disease_analyse", "n_cells", mp_order), with = FALSE],
                id.vars = c("SampleID_Global", "DatasetID", "Disease", "disease_analyse", "n_cells"),
                variable.name = "MP", value.name = "mean_usage")
pb_long[, `:=`(Disease = factor(Disease, disease_order), MP = factor(MP, mp_order))]
dis <- pb_long[, .(mean_usage = mean(mean_usage), sd_usage = sd(mean_usage),
                   n_sample = uniqueN(SampleID_Global), n_cells = sum(n_cells)), by = .(Disease, MP)]
setorder(dis, Disease, MP)
dis_raw_dt <- dcast(dis, MP ~ Disease, value.var = "mean_usage")
dis_mat <- as.matrix(dis_raw_dt[, -1]); rownames(dis_mat) <- dis_raw_dt$MP
dis_z <- t(scale(t(dis_mat))); dis_z[!is.finite(dis_z)] <- 0
dis[, row_zscore := as.numeric(dis_z[cbind(match(as.character(MP), rownames(dis_z)),
                                           match(as.character(Disease), colnames(dis_z)))])]
fwrite(dis, file.path(SRC, "PanelE_MP_disease_mean_usage_source.csv"))
fwrite(as.data.table(dis_z, keep.rownames = "MP"), file.path(SRC, "PanelE_MP_disease_mean_usage_rowZ_source.csv"))
fwrite(dcast(dis, Disease + n_sample + n_cells ~ MP, value.var = "mean_usage"),
       file.path(TBL, "PanelE_disease_MP_usage_wide.csv"))
n_by_dis <- dis[, .(n_sample = unique(n_sample)), by = Disease][match(disease_order, as.character(Disease))]
stopifnot(!anyNA(n_by_dis$n_sample))

draw_dis <- function(mat, legend_name, col_fun, stem) {
  ha <- HeatmapAnnotation(
    Disease = colnames(mat), n_sample = anno_text(n_by_dis$n_sample, gp = gpar(fontfamily = "Arial", fontsize = 7)),
    col = list(Disease = disease_colors[colnames(mat)]), show_legend = FALSE,
    annotation_name_gp = gpar(fontfamily = "Arial", fontsize = 8, fontface = "bold"))
  ht <- Heatmap(mat, name = legend_name, col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_gp = gpar(fontfamily = "Arial", fontsize = 8, fontface = "bold", col = "black"),
    column_names_gp = gpar(fontfamily = "Arial", fontsize = 8, col = "black"), column_names_rot = 45,
    top_annotation = ha, rect_gp = gpar(col = "white", lwd = 0.4))
  png(file.path(FIG, paste0(stem, ".png")), 3000, 2100, res = 300, type = "cairo")
  draw(ht, heatmap_legend_side = "right"); dev.off()
  cairo_pdf(file.path(FIG, paste0(stem, ".pdf")), 10, 7, family = "Arial")
  draw(ht, heatmap_legend_side = "right"); dev.off()
}
draw_dis(dis_mat, "Mean usage", colorRamp2(c(0, max(dis_mat)), c("white", "#2166AC")),
         "PanelE_MP_disease_mean_usage_raw")
draw_dis(dis_z, "Row Z", colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
         "PanelE_MP_disease_mean_usage")

audit <- data.table(
  item = c("default_MP_model", "object", "usage", "n_cells", "n_MP", "n_clusters_res0.6",
           "n_diseases", "n_samples_total", "n_samples_ge20", "n_samples_excluded_lt20",
           "disease_aggregation", "cluster_identity", "cell_id_alignment"),
  value = c("consensus_12MP", OBJ, USAGE, nrow(dt), 12L, length(cluster_order), 17L,
            nrow(pb_all), nrow(pb), sum(!pb_all$included_ge20),
            "sample mean -> equal-weight disease mean", "RNA_snn_res.0.6", "PASS")
)
fwrite(audit, file.path(TBL, "Figure2_12MP_usage_update_audit.csv"))
interpretation <- c(
  "# Figure 2 NMF 12 MP usage 更新解读", "",
  "## 分析目的",
  "以后续正式分析默认采用的 consensus 12 MP 为基础，按照当前 res=0.6 聚类和合并后的 Disease17 更新 usage 图；不重新拟合 NMF。", "",
  "## 输入数据",
  paste0("输入对象：", OBJ),
  paste0("全细胞 usage：", USAGE),
  paste0("共 ", nrow(dt), " 个细胞、12 个 MP、", length(cluster_order), " 个 res=0.6 簇、17 类疾病。"), "",
  "## 方法",
  paste0("以 SampleID_Global 为统计单位，仅纳入至少 ", MIN_CELLS, " 个细胞的样本。先在样本内对细胞 usage 求均值，再在 Disease 内对样本等权求均值；不使用 donor，也不按疾病细胞数加权。"),
  "Panel D 同时展示样本层级和 cluster 层级 usage；Panel E 展示 Disease17 的原始均值与每个 MP 跨疾病的行 Z 分数。", "",
  "## 关键数字",
  paste0("合格样本：", nrow(pb), "/", nrow(pb_all), "；排除低于门槛样本：", sum(!pb_all$included_ge20), "；CellID 对齐：PASS。"), "",
  "## 结果",
  "12 MP 在不同 res=0.6 簇和疾病间呈现不同的相对 usage 模式；原始 usage 与行 Z 分数均已输出，可分别用于观察绝对组成和相对富集。", "",
  "## 生物学意义",
  "该结果为后续跨疾病 MP 活性比较、疾病相关状态解释及样本级统计建模提供统一的12 MP输入。", "",
  "## 限制",
  "usage 热图属于描述性汇总，不等同于正式疾病效应检验；样本数很少的疾病必须通过独立的 sample-level 模型谨慎推断。"
)
writeLines(interpretation, file.path(ROOT, "Figure2_12MP_Disease17_interpretation.md"), useBytes = TRUE)
writeLines(interpretation, file.path(ROOT, "\u89e3\u8bfb.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(ROOT, "sessionInfo_12MP_usage_Disease17.txt"))
log_msg("DONE: 12 MP; ", length(cluster_order), " clusters; 17 diseases; ", nrow(pb), "/", nrow(pb_all), " samples.")
quit(save = "no", status = 0L, runLast = FALSE)
