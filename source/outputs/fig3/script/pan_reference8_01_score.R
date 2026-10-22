suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
  library(UCell)
  library(BiocParallel)
})

options(stringsAsFactors = FALSE)
set.seed(42)
write_utf8 <- function(text, path) {
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw(text), con)
}

FORMAL_OBJECT <- "E:/micro_output/v3_dims16/objects/merged_obj.rds"
PAN_ROOT <- "E:/micro_output/fig3/MP_FUN_PTG02/pan_disease_allcells_v3"
OUT_ROOT <- file.path(PAN_ROOT, "reference8")
DIRS <- file.path(OUT_ROOT, c("data", "tables", "source_data", "figures", "logs"))
META_FILE <- file.path(PAN_ROOT, "data/pan14_metadata_allcells.csv.gz")
CLUSTER_ORDER_FILE <- "E:/micro_output/fig3/data/05_cluster_row_order.csv"
CLUSTER_MAP_FILE <- "E:/micro_output/fig3/fig4/tables/cluster_name_mapping.csv"
DISEASE_GROUP_FILE <- "E:/micro_output/fig3/data/cluster_disease_group_map.csv"
DISEASE_PALETTE_FILE <- "E:/micro_output/fig3/data/disease_palette.csv"
TRAJECTORY_MAP_FILE <- file.path(PAN_ROOT, "data/P4_frozen_trajectory_allcells_score_mapping.csv.gz")

PAN14 <- c("C00", "C09", "C01", "C03", "C04", "C07", "C02", "C05", "C06", "C13", "C08", "C19", "C17", "C12")
EXCLUDED <- c("C10", "C15", "C20", "C24", "C27", "C28")
EXPECTED_PAN14 <- 365563L
BLOCK_SIZE <- 20000L

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))
required_paths <- c(FORMAL_OBJECT, META_FILE, CLUSTER_ORDER_FILE, CLUSTER_MAP_FILE)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) stop("缺少必要输入：", paste(missing_paths, collapse = "; "))

signature_list <- list(
  F8_Surveillance = c("OLFML3", "P2RY12", "TGFBR1", "TMEM119"),
  F8_Neuroprotection = c("BDNF", "GDNF", "IGF1", "LIF", "TGFB1"),
  F8_Phagocytosis = c("APOE", "AXL", "C3AR1", "CD63", "CD68", "CD9", "CLEC7A", "LGALS3", "MYO1E", "SPP1", "TREM2", "TYROBP"),
  F8_Inflammation = c("AOAH", "APOD", "C5AR1", "CCL8", "FGR", "MSR1", "RNF169"),
  F8_Cytokine_production = c("CCL2", "CXCL10", "IL10", "IL1B", "NFKBIA", "TNF"),
  F8_Antigen_presentation = c("B2M", "HLA-A", "HLA-DMA", "HLA-DQA1", "HLA-DQB1", "HLA-DRB5", "NLRC5"),
  F8_IFN_signature = c("CD69", "CXCL10", "IFI16", "IFIT2", "IFIT3", "IFITM3", "IRF7", "ISG15", "MNDA", "MX1", "USP18"),
  F8_Proliferation = c("BRIP1", "RAD51B", "TOP2A")
)
signature_labels <- c(
  F8_Surveillance = "监视 / Surveillance",
  F8_Neuroprotection = "神经保护 / Neuroprotection",
  F8_Phagocytosis = "吞噬 / Phagocytosis",
  F8_Inflammation = "炎症 / Inflammation",
  F8_Cytokine_production = "细胞因子产生 / Cytokine production",
  F8_Antigen_presentation = "抗原呈递 / Antigen presentation",
  F8_IFN_signature = "干扰素反应 / IFN signature",
  F8_Proliferation = "增殖 / Proliferation"
)

registry <- rbindlist(lapply(names(signature_list), function(sig) {
  data.table(signature = sig, signature_label = unname(signature_labels[[sig]]), gene = signature_list[[sig]], gene_position = seq_along(signature_list[[sig]]))
}), use.names = TRUE)
registry[, gene_unique_id := match(gene, unique(gene))]
fwrite(registry, file.path(OUT_ROOT, "tables/FUN8_T1_only_registry.csv"))

obj <- readRDS(FORMAL_OBJECT)
if (!inherits(obj, "Seurat")) stop("正式对象不是Seurat对象")
if (!"RNA" %in% Assays(obj)) stop("正式对象缺少RNA assay")
if (DefaultAssay(obj) != "RNA") DefaultAssay(obj) <- "RNA"
meta <- as.data.table(obj[[]], keep.rownames = "CellID")
required_meta <- c("CellID", "SampleID_Global", "DonorID", "DatasetID", "Disease", "RNA_snn_res.0.6")
missing_meta <- setdiff(required_meta, names(meta))
if (length(missing_meta)) stop("对象metadata缺少列：", paste(missing_meta, collapse = ", "))
if (anyDuplicated(meta$CellID)) stop("对象metadata存在重复CellID")
if (!identical(meta$CellID, colnames(obj))) stop("metadata CellID与对象列名顺序/集合不一致")
meta[, Cluster_raw := as.character(RNA_snn_res.0.6)]
meta[, Cluster := ifelse(grepl("^C", Cluster_raw), Cluster_raw, paste0("C", sprintf("%02d", as.integer(Cluster_raw))))]
pan_meta <- meta[Cluster %in% PAN14]
pan_meta[, cluster_order := match(Cluster, PAN14)]
setorder(pan_meta, cluster_order, CellID)
pan_ids <- pan_meta$CellID
if (anyDuplicated(pan_ids)) stop("pan14 CellID重复")

expr <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data")
expr_genes <- rownames(expr)
expr_cells <- colnames(expr)
if (is.null(expr_genes) || is.null(expr_cells)) stop("RNA data layer缺少行名或列名")
if (!identical(expr_cells, colnames(obj))) stop("RNA data layer列名与对象列名不一致")
expr_pan <- expr[, pan_ids, drop = FALSE]

gene_universe <- expr_genes
gene_audit <- registry[, .(
  n_requested = .N,
  n_in_universe = sum(gene %in% gene_universe),
  in_RNA_data_universe = all(gene %in% gene_universe),
  missing_genes = paste(gene[!gene %in% gene_universe], collapse = ";")
), by = .(signature, signature_label)]
gene_audit[, coverage_fraction := n_in_universe / n_requested]
gene_audit[, missing_gene_reason := fifelse(in_RNA_data_universe, "none", "requested gene absent from RNA data layer")]
fwrite(gene_audit, file.path(OUT_ROOT, "tables/FUN8_gene_mapping_audit.csv"))
if (any(!gene_audit$in_RNA_data_universe)) stop("至少一个8功能签名存在RNA data layer缺失基因；已输出审计，未静默缩短签名")
sig_obs <- lapply(signature_list, function(x) x)

meta_out <- pan_meta[, .(CellID, object_cell_name = CellID, SampleID_Global, DonorID, DatasetID, Disease, Cluster)]
meta_out[, sample_filter_applied := FALSE]
meta_out[, trajectory_sampling_used := FALSE]
fwrite(meta_out, file.path(OUT_ROOT, "data/FUN8_pan14_metadata_allcells.csv.gz"))

blocks <- split(seq_along(pan_ids), ceiling(seq_along(pan_ids) / BLOCK_SIZE))
score_blocks <- vector("list", length(blocks))
for (b in seq_along(blocks)) {
  idx <- blocks[[b]]
  message(sprintf("Scoring block %d/%d: cells %d-%d", b, length(blocks), min(idx), max(idx)))
  block_expr <- expr_pan[, idx, drop = FALSE]
  sc <- UCell::ScoreSignatures_UCell(
    matrix = block_expr,
    features = sig_obs,
    maxRank = 1500,
    ties.method = "average",
    chunk.size = 1000,
    ncores = 1,
    BPPARAM = BiocParallel::SerialParam(),
    force.gc = TRUE,
    name = "_UCell"
  )
  sc <- as.data.table(as.data.frame(sc), keep.rownames = "CellID")
  score_cols <- grep("_UCell$", names(sc), value = TRUE)
  if (length(score_cols) != length(signature_list)) stop("UCell输出列数异常：", length(score_cols))
  setnames(sc, score_cols, sub("_UCell$", "", score_cols))
  sc <- sc[, c("CellID", names(signature_list)), with = FALSE]
  if (!setequal(sc$CellID, pan_ids[idx])) stop("UCell block CellID集合不匹配")
  sc <- sc[match(pan_ids[idx], CellID)]
  score_blocks[[b]] <- sc
  rm(block_expr, sc)
  gc(verbose = FALSE)
}

scores <- rbindlist(score_blocks, use.names = TRUE)
if (nrow(scores) != length(pan_ids) || anyDuplicated(scores$CellID) || !identical(scores$CellID, pan_ids)) stop("合并后评分CellID覆盖不完整")
score_matrix <- as.matrix(scores[, names(signature_list), with = FALSE])
if (any(!is.finite(score_matrix))) stop("评分结果包含非有限值")
saveRDS(scores, file.path(OUT_ROOT, "data/FUN8_UCell_scores_allcells.rds"), compress = "xz")
fwrite(scores, file.path(OUT_ROOT, "data/FUN8_UCell_scores_allcells.csv.gz"))

coverage <- data.table(
  metric = c("formal_object_cells", "formal_pan14_cells", "FUN8_scored_cells", "FUN8_finite_cells", "FUN8_unique_CellID", "requested_signatures", "requested_gene_rows", "unique_requested_genes"),
  value = c(ncol(obj), length(pan_ids), nrow(scores), sum(complete.cases(score_matrix)), uniqueN(scores$CellID), length(signature_list), nrow(registry), uniqueN(registry$gene))
)
coverage <- rbind(coverage, data.table(
  metric = c("pan14_expected_reference", "pan14_coverage_fraction", "score_coverage_fraction", "sample_filter_applied", "trajectory_sampling_used", "duplicate_CellID"),
  value = c(EXPECTED_PAN14, length(pan_ids) / EXPECTED_PAN14, nrow(scores) / length(pan_ids), FALSE, FALSE, anyDuplicated(scores$CellID) > 0)
))
fwrite(coverage, file.path(OUT_ROOT, "tables/FUN8_cell_coverage_audit.csv"))

param_table <- data.table(
  parameter = c("algorithm", "maxRank", "ties.method", "chunk.size", "ncores", "BPPARAM", "force.gc", "external_block_size", "seed", "assay", "layer", "cell_scope", "sample_filter", "trajectory_sampling"),
  value = c("UCell::ScoreSignatures_UCell", "1500", "average", "1000", "1", "BiocParallel::SerialParam()", "TRUE", as.character(BLOCK_SIZE), "42", "RNA", "data", "complete pan14", "FALSE", "FALSE")
)
fwrite(param_table, file.path(OUT_ROOT, "tables/FUN8_UCell_parameters.csv"))

signature_summary <- registry[, .(signature_label = signature_label[1], n_requested = .N, n_unique_genes = uniqueN(gene), genes = paste(gene, collapse = ",")), by = signature]
fwrite(signature_summary, file.path(OUT_ROOT, "tables/FUN8_signature_summary.csv"))

hash_inputs <- c(FORMAL_OBJECT, META_FILE, CLUSTER_ORDER_FILE, CLUSTER_MAP_FILE)
hash_manifest <- data.table(path = hash_inputs, exists = file.exists(hash_inputs), md5 = vapply(hash_inputs, function(x) if (file.exists(x)) unname(tools::md5sum(x)) else NA_character_, character(1)))
fwrite(hash_manifest, file.path(OUT_ROOT, "tables/FUN8_input_hash_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(OUT_ROOT, "logs/FUN8_sessionInfo.txt"), useBytes = TRUE)

interpretation <- paste0(
  "# FUN8 T1-only全细胞评分\n\n",
  "## 分析目的\n",
  "在正式对象的完整pan14细胞范围内，按冻结的原始8个T1功能签名重新计算UCell评分，作为当前规划阶段A的基础结果。\n\n",
  "## 输入与范围\n",
  "- 正式对象：`", FORMAL_OBJECT, "`\n",
  "- pan14细胞：", format(length(pan_ids), big.mark = ","), "个；参考核对值：", format(EXPECTED_PAN14, big.mark = ","), "个\n",
  "- 亚簇顺序：", paste(PAN14, collapse = "、"), "\n",
  "- RNA assay既有data layer；未重算PCA、Harmony、UMAP、聚类或轨迹。\n\n",
  "## 方法\n",
  "使用`UCell::ScoreSignatures_UCell`，maxRank=1500、ties.method=average、chunk.size=1000、ncores=1、SerialParam、force.gc=TRUE；外部分块大小", BLOCK_SIZE, "，seed=42。未施加sample或sample×cluster门槛，未使用轨迹抽样细胞限制评分。\n\n",
  "## 结果\n",
  "8个签名、", nrow(registry), "条签名基因记录、", uniqueN(registry$gene), "个唯一请求基因均完成RNA data layer映射；评分覆盖", nrow(scores), "/", length(pan_ids), "个pan14细胞（", format(nrow(scores) / length(pan_ids), digits = 6), "），有限值覆盖", format(sum(complete.cases(score_matrix)) / nrow(scores), digits = 6), "。\n\n",
  "## 限制与边界\n",
  "本结果只完成T1-only 8功能评分和审计，不替代既有12FUN、12MP、PTG或轨迹结果，也不启动规划中的B–F阶段。"
)
interpretation_file <- file.path(OUT_ROOT, paste0(intToUtf8(c(35299L, 35835L)), ".md"))
write_utf8(interpretation, interpretation_file)

writeLines(capture.output(list(
  completed = Sys.time(), formal_object = FORMAL_OBJECT, pan14_cells = length(pan_ids), scored_cells = nrow(scores), signatures = names(signature_list), sample_filter = FALSE, trajectory_sampling = FALSE
)), file.path(OUT_ROOT, "logs/FUN8_score_run.log"), useBytes = TRUE)

rm(obj, expr, expr_pan, score_blocks, score_matrix)
gc(verbose = FALSE)
cat("Stage A FUN8 T1-only scoring complete.\n")
