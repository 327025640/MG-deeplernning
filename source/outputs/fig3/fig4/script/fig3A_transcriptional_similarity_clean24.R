# fig3A_transcriptional_similarity_clean24.R
# 正式 Figure 3A：仅使用 24 个干净 MG 亚簇（排除用户固定污染簇 C10/C15/C20/C24/C27/C28）
# 的转录相似性层次聚类树。
#
# 与 30 簇 QC 版（fig3A_transcriptional_similarity.R，保留为 QC/supplementary）的差异：
#   1) 只纳入 24 个干净亚簇；叶标签格式 "C00 | MG_TLR6"（簇 ID + 正式名）；
#   2) 术语统一：sample-averaged log-normalized expression profile
#      （对 log-normalized 表达先按样本求均值、再在簇内对合格样本等权平均；
#      未聚合原始 counts，不得称为 pseudobulk）；
#   3) 距离 Euclidean、linkage Ward.D2（论文未报告 linkage，审计留痕）；
#   4) 不运行 cutree(k=10)，也不因论文得到 10 分支而强制切 10；
#   5) summary 使用 n_eligible_cells（通过阈值的细胞）与 n_cells_total（簇内真实总细胞），
#      新增 n_eligible_samples 与 stability_warning；
#   6) 修复底部副标题与横轴标题重叠；
#   7) 不根据树状图反向修改亚簇命名。
# 不覆盖 Fig3A.png / Fig3A.pdf（30 簇 QC 版原样保留）。

options(stringsAsFactors = FALSE)
set.seed(42)

INPUT_RDS <- "E:/micro_output/v3_dims16/objects/merged_obj.rds"
MAPPING_CSV <- "E:/micro_output/fig3/fig4/data/cluster_name_mapping_used.csv"
OUT_ROOT <- "E:/micro_output/fig3/fig4"

CLUSTER_COL <- "RNA_snn_res.0.6"
CONTAMINANT_CLUSTERS <- sprintf("C%02d", c(10, 15, 20, 24, 27, 28))
MIN_CELLS_PER_SAMPLE_CLUSTER <- 20L
STABILITY_MIN_SAMPLES <- 10L   # 合格样本数 < 10 时给出 stability_warning（适配阈值，写入审计）
DIST_METHOD <- "euclidean"
LINKAGE <- "ward.D2"           # 论文未报告 linkage；适配选择，写入审计

FIG_DIR <- file.path(OUT_ROOT, "figures")
TABLE_DIR <- file.path(OUT_ROOT, "tables")
SOURCE_DIR <- file.path(OUT_ROOT, "source_data")
DATA_DIR <- file.path(OUT_ROOT, "data")

required_pkgs <- c("Seurat", "SeuratObject", "Matrix", "data.table")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs)) {
  stop("缺少R包（脚本不会自动安装）：", paste(missing_pkgs, collapse = ", "))
}
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(data.table)
})

for (p in c(INPUT_RDS, MAPPING_CSV)) {
  if (!file.exists(p)) stop("输入不存在：", p)
}
for (d in c(FIG_DIR, TABLE_DIR, SOURCE_DIR, DATA_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

###############################################################################
# 1. 命名映射表：仅保留 24 个干净亚簇
###############################################################################

mapping <- fread(MAPPING_CSV)
stopifnot(nrow(mapping) == 30L, !anyDuplicated(mapping$Cluster))
stopifnot(setequal(mapping$Cluster, sprintf("C%02d", 0:29)))
stopifnot(all(mapping[Cluster %in% CONTAMINANT_CLUSTERS, ContaminantFlag] == TRUE))

clean_mapping <- mapping[ContaminantFlag != TRUE][
  order(as.integer(sub("C", "", Cluster)))
]
stopifnot(nrow(clean_mapping) == 24L, !anyDuplicated(clean_mapping$Cluster))
clean_clusters <- clean_mapping$Cluster
name_lab <- setNames(
  paste0(clean_mapping$Cluster, " | ", clean_mapping$FinalName),
  clean_mapping$Cluster
)   # 叶标签："C00 | MG_TLR6"
cat(sprintf("[1] 干净亚簇 %d 个（已排除 %s）\n", length(clean_clusters),
            paste(CONTAMINANT_CLUSTERS, collapse = ",")))

###############################################################################
# 2. 正式对象校验（只读；不重算任何降维/聚类；不写回）
###############################################################################

obj <- readRDS(INPUT_RDS)
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "RNA_snn_res.0.6"
stopifnot(
  ncol(obj) == 424222,
  nrow(obj) == 18687,
  "RNA_snn_res.0.6" %in% colnames(obj[[]]),
  "seurat_clusters" %in% colnames(obj[[]]),
  "Disease" %in% colnames(obj[[]]),
  "SampleID_Global" %in% colnames(obj[[]]),
  identical(as.character(obj$seurat_clusters), as.character(obj$RNA_snn_res.0.6))
)

md <- obj[[]]
cl_vec <- sprintf("C%02d", as.integer(as.character(md[[CLUSTER_COL]])))
samp_vec <- as.character(md$SampleID_Global)
stopifnot(setequal(unique(cl_vec), sprintf("C%02d", 0:29)))

data_layer <- LayerData(obj, assay = "RNA", layer = "data")
gene_names <- rownames(data_layer)
rm(obj); invisible(gc())

###############################################################################
# 3. sample-averaged log-normalized expression profile（每 sample-cluster >=20 细胞）
###############################################################################

prof_cluster <- matrix(NA_real_, nrow = length(gene_names),
                       ncol = length(clean_clusters),
                       dimnames = list(gene_names, clean_clusters))
elig_list <- list(); meta_list <- list()

for (cc in clean_clusters) {
  idx <- which(cl_vec == cc)
  n_cells_total <- length(idx)
  sub <- data_layer[, idx, drop = FALSE]
  svec <- samp_vec[idx]
  tb <- table(svec)
  ok_samples <- names(tb)[tb >= MIN_CELLS_PER_SAMPLE_CLUSTER]
  pm <- matrix(0, nrow = nrow(sub), ncol = length(ok_samples),
               dimnames = list(gene_names, ok_samples))
  for (s in ok_samples) {
    pm[, s] <- Matrix::rowMeans(sub[, svec == s, drop = FALSE])
  }
  # 簇谱 = 合格样本的等权平均 = sample-averaged log-normalized expression profile
  prof_cluster[, cc] <- rowMeans(pm)
  elig_list[[cc]] <- data.table(
    Cluster = cc, SampleID_Global = ok_samples,
    n_cells = as.integer(tb[ok_samples]), Eligible = TRUE
  )
  excl <- setdiff(names(tb), ok_samples)
  if (length(excl)) {
    elig_list[[cc]] <- rbind(
      elig_list[[cc]],
      data.table(Cluster = cc, SampleID_Global = excl,
                 n_cells = as.integer(tb[excl]), Eligible = FALSE)
    )
  }
  meta_list[[cc]] <- data.table(
    Cluster = cc,
    n_cells_total = n_cells_total,
    n_eligible_cells = as.integer(sum(tb[ok_samples])),
    n_eligible_samples = length(ok_samples),
    n_sample_cluster_pairs_total = length(tb)
  )
  cat(sprintf("  %s: n_cells_total=%d | n_eligible_cells=%d | n_eligible_samples=%d\n",
              cc, n_cells_total, sum(tb[ok_samples]), length(ok_samples)))
  rm(sub); invisible(gc())
}
elig_all <- rbindlist(elig_list)
cluster_meta <- rbindlist(meta_list)
cluster_meta[, stability_warning := fifelse(
  n_eligible_samples < STABILITY_MIN_SAMPLES,
  paste0("only ", n_eligible_samples, " eligible samples (<", STABILITY_MIN_SAMPLES,
         "); leaf position may be unstable"),
  ""
)]
fwrite(elig_all, file.path(SOURCE_DIR, "Fig3A_clean24_sample_cluster_eligibility.csv"))

###############################################################################
# 4. 基因过滤 + Euclidean + Ward.D2
###############################################################################

keep_genes <- rowMeans(prof_cluster) > 0
prof_used <- prof_cluster[keep_genes, , drop = FALSE]
cat(sprintf("[4] 基因过滤：%d/%d 基因在至少一个簇的 sample-averaged profile 中表达>0\n",
            sum(keep_genes), nrow(prof_cluster)))

d_mat <- dist(t(prof_used), method = DIST_METHOD)
hc <- hclust(d_mat, method = LINKAGE)
leaf_order_clusters <- hc$labels[hc$order]     # 簇 ID 的叶顺序

###############################################################################
# 5. 树状图（叶标签 "C00 | MG_TLR6"；正式结果不运行 cutree(k=10)）
###############################################################################

dend <- as.dendrogram(hc)
dend_col <- dendrapply(dend, function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    cl <- a$label
    attr(n, "label") <- unname(name_lab[cl])
    attr(n, "nodePar") <- c(a$nodePar, list(lab.col = "black", lab.cex = 0.72, pch = NA))
  }
  n
})

fig_png <- file.path(FIG_DIR, "Fig3A_clean24.png")
fig_pdf <- file.path(FIG_DIR, "Fig3A_clean24.pdf")
draw_tree <- function() {
  par(mar = c(4.6, 4.5, 3.5, 13.2), family = "Arial", xpd = TRUE)
  plot(dend_col, horiz = TRUE, axes = FALSE,
       main = "Transcriptional similarity of the 24 clean microglia subclusters",
       sub = "", xlab = "")
  title(sub = paste0("sample-averaged log-normalized expression profile (",
                     sum(keep_genes), " genes); Euclidean distance; Ward.D2 linkage ",
                     "(linkage not reported in the original paper)"),
        line = 0.4, cex.sub = 0.72, col.sub = "grey30")
  mtext("Cluster dissimilarity (Euclidean)", side = 1, line = 2.6, cex = 0.8)
}
png(fig_png, width = 11, height = 8, units = "in", res = 600,
    bg = "white", type = "cairo", family = "Arial")
draw_tree()
dev.off()
cairo_pdf(fig_pdf, width = 11, height = 8, family = "Arial")
draw_tree()
dev.off()
cat("[5] 图形输出完成：", fig_png, "|", fig_pdf, "\n")

###############################################################################
# 6. 叶顺序（供 Fig3B/3C 使用）+ source data + summary + method audit
###############################################################################

leaf_order_dt <- data.table(
  leaf_order = seq_along(leaf_order_clusters),
  Cluster = leaf_order_clusters,
  FinalNameWithStatus = clean_mapping$FinalNameWithStatus[
    match(leaf_order_clusters, clean_mapping$Cluster)
  ]
)
fwrite(leaf_order_dt, file.path(DATA_DIR, "Fig3A_clean24_leaf_order.csv"))

dist_out <- as.data.table(as.matrix(d_mat), keep.rownames = "Cluster")
fwrite(dist_out, file.path(SOURCE_DIR, "Fig3A_clean24_distance_matrix.csv"))
prof_out <- as.data.table(prof_used, keep.rownames = "Gene")
fwrite(prof_out, file.path(SOURCE_DIR, "Fig3A_clean24_sample_averaged_profiles.csv"))

summary_tab <- clean_mapping[, .(Cluster, FinalName, FinalNameWithStatus,
                                 RepresentativeGene, ContaminantFlag = FALSE,
                                 Confidence, ReviewRequired)]
summary_tab <- merge(
  summary_tab, cluster_meta[, .(Cluster, n_eligible_cells, n_cells_total,
                                n_eligible_samples, n_sample_cluster_pairs_total,
                                stability_warning)],
  by = "Cluster", all.x = TRUE
)
summary_tab[, n_eligible_samples_warning_threshold := STABILITY_MIN_SAMPLES]
setorder(summary_tab, Cluster)
fwrite(summary_tab, file.path(TABLE_DIR, "Fig3A_clean24_summary.csv"))

method_audit <- data.table(
  Item = c(
    "Reference_panel", "Formal_scope", "QC30_version_status",
    "Paper_input", "Paper_method", "Paper_linkage",
    "Adaptation_clusters", "Adaptation_expression",
    "Adaptation_terminology", "Adaptation_statistical_unit",
    "Min_cells_per_sample_cluster", "Adaptation_gene_filter",
    "Adaptation_distance", "Adaptation_linkage", "Adaptation_linkage_note",
    "Cutree_k10", "Leaf_label_format", "Leaf_label_source",
    "n_eligible_cells_definition", "n_cells_total_definition",
    "Stability_warning_rule", "Downstream_leaf_order",
    "Figure_devices", "Statistical_interpretation"
  ),
  Value = c(
    "Nie et al. Cell Discovery 2026 Figure 3A (transcriptional similarity dendrogram)",
    "formal main figure: 24 clean MG subclusters only",
    "30-cluster version (Fig3A.png/Fig3A.pdf) retained as QC/supplementary only",
    "subtype x gene transcriptional profiles (per-subtype mean expression)",
    "Euclidean distance matrix + unsupervised hierarchical clustering dendrogram",
    "not reported in the paper",
    paste0(length(clean_clusters), " clean subclusters (contaminants ",
           paste(CONTAMINANT_CLUSTERS, collapse = ","), " excluded)"),
    "Seurat RNA data layer (log-normalized)",
    "sample-averaged log-normalized expression profile (NOT pseudobulk; raw counts not aggregated)",
    paste0("per sample-cluster mean of log-normalized expression (>= ",
           MIN_CELLS_PER_SAMPLE_CLUSTER, " cells), then equal-weight mean across ",
           "eligible samples within each cluster"),
    MIN_CELLS_PER_SAMPLE_CLUSTER,
    "genes with mean profile > 0 in at least one cluster",
    DIST_METHOD, LINKAGE,
    "linkage not reported in the original paper; Ward.D2 recorded as an adaptation choice",
    "NOT run for the formal result; no forced 10-branch cut",
    "C00 | MG_TLR6 (cluster id + final name)",
    "fig4/data/cluster_name_mapping_used.csv via FinalName",
    "cells inside eligible sample-cluster pairs (>= 20 cells) that enter the profile",
    "all cells assigned to the cluster regardless of eligibility",
    paste0("n_eligible_samples < ", STABILITY_MIN_SAMPLES, " triggers stability_warning ",
           "(project adaptation threshold)"),
    "fig4/data/Fig3A_clean24_leaf_order.csv (consumed by Fig3B/3C)",
    "PNG 600 dpi (cairo, white bg); PDF cairo_pdf (editable text); Arial",
    "Descriptive similarity tree; not proof of lineage hierarchy or developmental ordering"
  )
)
fwrite(method_audit, file.path(TABLE_DIR, "Fig3A_clean24_method_audit.csv"))
writeLines(capture.output(sessionInfo()),
           file.path(TABLE_DIR, "Fig3A_clean24_sessionInfo.txt"), useBytes = TRUE)

message("Fig3A_clean24 完成：", length(clean_clusters), " clean clusters x ",
        sum(keep_genes), " genes; linkage=", LINKAGE, "; 无 cutree(k=10)。")
