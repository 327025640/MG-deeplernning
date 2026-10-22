suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
  library(UCell)
  library(BiocParallel)
})

options(stringsAsFactors = FALSE)
set.seed(42)

FIG3 <- "D:/data/microglia_new/analysis_scripts/runtime_io_aliases/fig3"
FORMAL_OBJECT <- "E:/micro_output/v3_dims16/objects/merged_obj.rds"
CELL_MANIFEST <- file.path(FIG3, "data/cell_manifest_clean24.csv.gz")
PAN14_FUN8 <- "D:/data/microglia_new/analysis_scripts/runtime_io_aliases/fun8_reviewed/source_data/FUN8_heatmap_source.csv"
MP12_CACHE <- "E:/micro_output/fig3/shared_function_framework/tables/MP_cluster_usage_mean.csv"
PARAM_LOG <- file.path(FIG3, "tables/S1_nonpan10_FUN8_parameters.csv")
SCORE_FILE <- file.path(FIG3, "data/nonpan10_FUN8_scores_with_metadata.csv.gz")
RUN_LOG <- file.path(FIG3, "logs/S1_nonpan10_FUN8_run.log")
BLOCK_SIZE <- 20000L
NONPAN10 <- c("C11", "C14", "C16", "C18", "C21", "C22", "C23", "C25", "C26", "C29")
PAN14 <- c("C00", "C09", "C01", "C03", "C04", "C07", "C02", "C05", "C06", "C13", "C08", "C19", "C17", "C12")
CLUSTER_ORDER <- c(PAN14, "C11", "C14", "C16", "C18", "C21", "C22", "C23", "C25", "C26", "C29")

required_files <- c(FORMAL_OBJECT, CELL_MANIFEST, PAN14_FUN8, MP12_CACHE)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing frozen input(s): ", paste(missing_files, collapse = "; "))
if (file.exists(SCORE_FILE)) stop("Refusing to overwrite existing nonpan10 FUN8 output: ", SCORE_FILE)
dir.create(dirname(SCORE_FILE), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(PARAM_LOG), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(RUN_LOG), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIG3, "source_data"), recursive = TRUE, showWarnings = FALSE)

started <- Sys.time()
phase_times <- data.table(phase = character(), seconds = numeric(), detail = character())
log_phase <- function(phase, phase_start, detail) {
  phase_times <<- rbind(phase_times, data.table(phase = phase, seconds = as.numeric(difftime(Sys.time(), phase_start, units = "secs")), detail = detail))
  fwrite(phase_times, RUN_LOG)
}

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
registry <- rbindlist(lapply(names(signature_list), function(sig) data.table(signature = sig, gene = signature_list[[sig]])))

manifest <- fread(CELL_MANIFEST, na.strings = c("", "NA"))
required_manifest <- c("CellID", "Cluster", "Disease", "DatasetID", "SampleID_Global", "sample_mapping_status")
if (length(setdiff(required_manifest, names(manifest)))) stop("clean24 manifest lacks required columns")
if (anyDuplicated(manifest$CellID)) stop("clean24 manifest has duplicate CellID")
if (!setequal(unique(manifest$Cluster), CLUSTER_ORDER)) stop("clean24 manifest does not contain the exact frozen 24-cluster scope")
nonpan <- manifest[Cluster %in% NONPAN10]
if (nrow(nonpan) == 0L || anyDuplicated(nonpan$CellID) || !setequal(unique(nonpan$Cluster), NONPAN10)) stop("nonpan10 scope is empty, duplicated, or incomplete")
if (!all(nonpan$sample_mapping_status %in% c("matched_frozen_sample_map", "unmapped_in_frozen_sample_map"))) stop("Unexpected sample mapping status")

read_start <- Sys.time()
obj <- readRDS(FORMAL_OBJECT)
if (!inherits(obj, "Seurat")) stop("Formal object is not Seurat")
if (!"RNA" %in% Assays(obj)) stop("Formal object lacks RNA assay")
meta <- as.data.table(obj[[]], keep.rownames = "CellID")
if (anyDuplicated(meta$CellID) || !identical(meta$CellID, colnames(obj))) stop("Object metadata CellID order does not exactly match object columns")
if (!"seurat_clusters" %in% names(meta)) stop("Formal object metadata lacks frozen seurat_clusters")
if (!setequal(nonpan$CellID, intersect(nonpan$CellID, meta$CellID))) stop("Some clean24 nonpan cells are absent from the formal object")
object_cluster_raw <- as.character(meta$seurat_clusters[match(nonpan$CellID, meta$CellID)])
if (anyNA(object_cluster_raw) || any(!grepl("^[0-9]{1,2}$", object_cluster_raw))) stop("Unexpected raw seurat_clusters labels in the formal object")
object_cluster <- paste0("C", sprintf("%02d", as.integer(object_cluster_raw)))
if (anyNA(object_cluster) || !identical(object_cluster, nonpan$Cluster)) stop("Clean24 cluster labels do not match formal object seurat_clusters by CellID")
ordered_idx <- match(colnames(obj), nonpan$CellID)
ordered_idx <- ordered_idx[!is.na(ordered_idx)]
nonpan <- nonpan[ordered_idx]
if (!identical(nonpan$CellID, colnames(obj)[colnames(obj) %in% nonpan$CellID])) stop("nonpan manifest/object cell order reconciliation failed")
DefaultAssay(obj) <- "RNA"
expr <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data")
if (is.null(rownames(expr)) || is.null(colnames(expr))) stop("RNA data layer is missing gene or cell names")
if (!identical(colnames(expr), colnames(obj))) stop("RNA data layer cell order differs from formal object")
if (!is(expr, "sparseMatrix")) stop("RNA data layer is not sparse; stop before creating a dense full matrix")
if (!all(registry$gene %in% rownames(expr))) {
  missing_genes <- unique(registry[!gene %in% rownames(expr), gene])
  fwrite(data.table(missing_gene = missing_genes), file.path(FIG3, "tables/S1_FUN8_missing_signature_genes.csv"))
  stop("Frozen FUN8 signature genes are missing from RNA data layer")
}
log_phase("load_object_and_validate_RNA_data", read_start,
          sprintf("object_cells=%d; object_genes=%d; RNA_assay=RNA; layer=data; sparse=%s; nonpan10_cells=%d",
                  ncol(obj), nrow(expr), is(expr, "sparseMatrix"), nrow(nonpan)))

score_start <- Sys.time()
expr_nonpan <- expr[, nonpan$CellID, drop = FALSE]
if (!identical(colnames(expr_nonpan), nonpan$CellID)) stop("RNA nonpan cell ordering changed after subsetting")
blocks <- split(seq_len(ncol(expr_nonpan)), ceiling(seq_len(ncol(expr_nonpan)) / BLOCK_SIZE))
score_blocks <- vector("list", length(blocks))
for (b in seq_along(blocks)) {
  idx <- blocks[[b]]
  message(sprintf("Scoring nonpan FUN8 block %d/%d: %d cells", b, length(blocks), length(idx)))
  block_expr <- expr_nonpan[, idx, drop = FALSE]
  scored <- UCell::ScoreSignatures_UCell(
    matrix = block_expr,
    features = signature_list,
    maxRank = 1500,
    ties.method = "average",
    chunk.size = 1000,
    ncores = 1,
    BPPARAM = BiocParallel::SerialParam(),
    force.gc = TRUE,
    name = "_UCell"
  )
  scored <- as.data.table(as.data.frame(scored), keep.rownames = "CellID")
  score_cols <- grep("_UCell$", names(scored), value = TRUE)
  if (length(score_cols) != length(signature_list)) stop("UCell output signature count mismatch")
  setnames(scored, score_cols, sub("_UCell$", "", score_cols))
  scored <- scored[, c("CellID", names(signature_list)), with = FALSE]
  if (!setequal(scored$CellID, nonpan$CellID[idx])) stop("UCell output CellID set differs from the current block")
  scored <- scored[match(nonpan$CellID[idx], CellID)]
  if (any(!is.finite(as.matrix(scored[, names(signature_list), with = FALSE])))) stop("UCell block has non-finite values")
  score_blocks[[b]] <- scored
  rm(block_expr, scored)
  gc(verbose = FALSE)
}
scores <- rbindlist(score_blocks, use.names = TRUE)
if (nrow(scores) != nrow(nonpan) || anyDuplicated(scores$CellID) || !identical(scores$CellID, nonpan$CellID)) stop("Nonpan FUN8 CellID coverage/order is incomplete")
if (any(as.matrix(scores[, names(signature_list), with = FALSE]) < 0) || any(as.matrix(scores[, names(signature_list), with = FALSE]) > 1)) stop("UCell scores are outside [0,1]")
log_phase("UCell_scoring", score_start, sprintf("cells=%d; signatures=%d; external_block=%d; internal_chunk=1000", nrow(scores), length(signature_list), BLOCK_SIZE))

write_start <- Sys.time()
score_cols <- setdiff(names(scores), "CellID")
score_output <- merge(nonpan, scores, by = "CellID", all = FALSE, sort = FALSE)
score_output <- score_output[match(nonpan$CellID, CellID)]
if (nrow(score_output) != nrow(nonpan) || anyDuplicated(score_output$CellID)) stop("Joined FUN8 metadata changed CellID coverage")
score_output[, `:=`(assay = "RNA", layer = "data", score_method = "UCell::ScoreSignatures_UCell",
                    frozen_pan14_cache_overlap = FALSE)]
fwrite(score_output, SCORE_FILE)
score_summary <- score_output[, c(lapply(.SD, mean), list(n_cells = .N)), by = Cluster, .SDcols = score_cols]
pan <- fread(PAN14_FUN8)
if (!all(c("Cluster", "Function", "raw_mean") %in% names(pan))) stop("Frozen pan14 FUN8 cache lacks Cluster/Function/raw_mean")
pan <- pan[, .(Cluster, Function, raw_mean)]
if (anyDuplicated(pan[, .(Cluster, Function)]) || nrow(pan) != length(PAN14) * length(signature_list) || !setequal(unique(pan$Cluster), PAN14)) stop("Frozen pan14 FUN8 cache is not exactly 14x8")
pan[, n_cells := NA_integer_]
nonpan_long <- melt(score_summary, id.vars = c("Cluster", "n_cells"), variable.name = "Function", value.name = "raw_mean")
nonpan_long[, origin := "new_nonpan10_full_RNA_UCell_cell_mean"]
pan[, origin := "frozen_pan14_cluster_raw_mean_cache"]
fun8_long <- rbindlist(list(pan, nonpan_long), use.names = TRUE, fill = TRUE)
fun8_long[, cluster_order := match(Cluster, CLUSTER_ORDER)]
fun8_long[, function_order := match(Function, names(signature_list))]
if (nrow(fun8_long) != 24L * 8L || anyDuplicated(fun8_long[, .(Cluster, Function)]) || anyNA(fun8_long$raw_mean)) stop("Combined FUN8 summary is incomplete or contains missing raw means")
fun8_long[, `:=`(display_z = NA_real_, z_status = NA_character_)]
for (f in names(signature_list)) {
  idx <- which(fun8_long$Function == f)
  vals <- fun8_long$raw_mean[idx]
  scale <- sd(vals)
  if (!is.finite(scale) || scale == 0) {
    fun8_long$z_status[idx] <- "no_between_cluster_variation"
  } else {
    fun8_long$display_z[idx] <- (vals - mean(vals)) / scale
    fun8_long$z_status[idx] <- "standardized_across_24_clusters_for_display"
  }
}
setorder(fun8_long, cluster_order, function_order)
fwrite(fun8_long, file.path(FIG3, "tables/cluster24_FUN8_summary.csv"))
fwrite(fun8_long, file.path(FIG3, "source_data/cluster24_FUN8_source.csv"))
fun8_wide <- dcast(fun8_long, Cluster + cluster_order ~ Function, value.var = "raw_mean")
fwrite(fun8_wide, file.path(FIG3, "tables/cluster24_FUN8_summary_wide.csv"))

mp <- fread(MP12_CACHE)
if (!"cluster_c" %in% names(mp) || !"n_samples" %in% names(mp)) stop("Frozen MP12 cache lacks cluster_c/n_samples")
mp_cols <- paste0("MP", sprintf("%02d", 1:12))
if (!all(mp_cols %in% names(mp)) || nrow(mp) != 24L || anyDuplicated(mp$cluster_c) || !setequal(mp$cluster_c, CLUSTER_ORDER)) stop("Frozen MP12 cache is not exactly 24x12")
mp_summary <- copy(mp)
setnames(mp_summary, "cluster_c", "Cluster")
mp_summary[, cluster_order := match(Cluster, CLUSTER_ORDER)]
setorder(mp_summary, cluster_order)
fwrite(mp_summary, file.path(FIG3, "tables/cluster24_MP12_summary.csv"))
mp_long <- melt(mp_summary, id.vars = c("Cluster", "cluster_order", "n_samples"), measure.vars = mp_cols,
                variable.name = "Function", value.name = "raw_mean")
mp_long[, `:=`(function_order = match(Function, mp_cols), display_z = NA_real_, z_status = NA_character_)]
for (f in mp_cols) {
  idx <- which(mp_long$Function == f)
  vals <- mp_long$raw_mean[idx]
  scale <- sd(vals)
  if (!is.finite(scale) || scale == 0) mp_long$z_status[idx] <- "no_between_cluster_variation" else {
    mp_long$display_z[idx] <- (vals - mean(vals)) / scale
    mp_long$z_status[idx] <- "standardized_across_24_clusters_for_display"
  }
}
fwrite(mp_long, file.path(FIG3, "source_data/cluster24_MP12_source.csv"))
param <- data.table(
  parameter = c("algorithm", "maxRank", "ties.method", "chunk.size", "ncores", "BPPARAM", "force.gc", "external_block_size", "seed", "assay", "layer", "input_gene_ranking_scope", "cell_scope", "pan14_cache_handling"),
  value = c("UCell::ScoreSignatures_UCell", "1500", "average", "1000", "1", "BiocParallel::SerialParam()", "TRUE", as.character(BLOCK_SIZE), "42", "RNA", "data", "all RNA data-layer genes retained before UCell rank calculation", "all clean24 nonpan10 cells including Ctrl", "reuse frozen pan14 raw_mean; no rescoring")
)
fwrite(param, PARAM_LOG)
log_phase("write_scores_and_24_cluster_summaries", write_start,
          sprintf("nonpan_rows=%d; FUN8_grid=%d; MP12_grid=%d", nrow(score_output), nrow(fun8_long), nrow(mp_long)))
fwrite(data.table(metric = c("formal_object_cells", "clean24_cells", "nonpan10_scored_cells", "pan14_clusters_reused", "clean24_clusters_FUN8", "clean24_clusters_MP12", "FUN8_signatures", "FUN8_unique_signature_genes"),
                  value = c(ncol(obj), nrow(manifest), nrow(score_output), length(unique(pan$Cluster)), length(unique(fun8_long$Cluster)), nrow(mp_summary), length(signature_list), uniqueN(registry$gene))),
      file.path(FIG3, "tables/S1_function_score_coverage.csv"))
writeLines(capture.output(sessionInfo()), file.path(FIG3, "logs/S1_FUN8_sessionInfo.txt"), useBytes = TRUE)
cat(sprintf("S1 nonpan10 FUN8 completed: %d cells, %.1f elapsed seconds\n", nrow(score_output), as.numeric(difftime(Sys.time(), started, units = "secs"))))
