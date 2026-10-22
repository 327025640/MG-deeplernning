#!/usr/bin/env Rscript

project_lib <- "D:/data/microglia_new/Rlibs/Fig4C_R4.4.3"
.libPaths(c(project_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
  library(decoupleR)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 04_cluster_mean_ulm.R <formal_merged_obj.rds|--resume-checkpoints> <staging_dir>")
}
resume_from_checkpoints <- identical(args[[1L]], "--resume-checkpoints")
object_path <- if (resume_from_checkpoints) {
  "E:/micro_output/v3_dims16/objects/merged_obj.rds (not accessed; checkpoint resume)"
} else {
  normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
}
stage_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
manifest_path <- file.path(stage_dir, "cell_manifest_clean24.csv.gz")
order_path <- file.path(stage_dir, "cluster_display_order.csv")
scope_path <- file.path(stage_dir, "cluster_scope.csv")
prior_path <- file.path(stage_dir, "collectri_human_viper_regulon_20260624.csv")
anchor_geometry_path <- file.path(stage_dir, "anchor_row_center_geometry.csv")
cache_audit_path <- file.path(stage_dir, "expression_cache_search_audit.csv")
required_paths <- c(manifest_path, order_path, scope_path, prior_path,
                    anchor_geometry_path, cache_audit_path)
if (!all(file.exists(required_paths))) {
  stop("One or more frozen Fig4C inputs are missing from the staging directory.")
}

output_dir <- file.path(stage_dir, "ulm_outputs")
checkpoint_dir <- file.path(output_dir, "cluster_mean_checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

expected_order <- c(
  "C11", "C21", "C14", "C18", "C16", "C22",
  "C23", "C25", "C26", "C29",
  "C00", "C09", "C01", "C03", "C04", "C07", "C02", "C05",
  "C06", "C13", "C08", "C19", "C17", "C12"
)

manifest <- fread(manifest_path)
order_table <- fread(order_path)
scope_table <- fread(scope_path)
prior <- fread(prior_path)
anchor_geometry <- fread(anchor_geometry_path)
cache_audit <- fread(cache_audit_path)
required_manifest <- c("CellID", "Cluster")
if (!all(required_manifest %in% names(manifest))) {
  stop("Frozen clean24 manifest lacks CellID or Cluster.")
}
if (!all(c("row_top_to_bottom", "Cluster") %in% names(order_table)) ||
    !identical(as.character(order_table[order(row_top_to_bottom), Cluster]), expected_order)) {
  stop("cluster_display_order.csv differs from the user-confirmed 24-cluster order.")
}
if (!"scope_group" %in% names(order_table)) {
  stop("cluster_display_order.csv lacks its frozen cluster-group annotation.")
}
if (!all(c("row_top_to_bottom", "Cluster", "center_mm_from_top", "tile_height_mm") %in% names(anchor_geometry)) ||
    !identical(as.character(anchor_geometry[order(row_top_to_bottom), Cluster]), expected_order)) {
  stop("Anchor row-center geometry is incomplete or differs from the confirmed FUN8 order.")
}
if (!"compatible_cache_found" %in% names(cache_audit) || any(cache_audit$compatible_cache_found)) {
  stop("A compatible expression-mean cache was found or cache audit schema is invalid; resolve cache identity before aggregation.")
}
if (!all(c("Cluster", "analysis_scope") %in% names(scope_table))) {
  stop("Frozen cluster_scope.csv lacks the required group mapping.")
}
if (!all(c("source", "target", "mor") %in% names(prior))) {
  stop("The audited CollecTRI regulon lacks source, target, or mor columns.")
}
if (nrow(manifest) != 400026L || anyDuplicated(manifest$CellID) != 0L ||
    !setequal(unique(manifest$Cluster), expected_order)) {
  stop("Frozen clean24 manifest identity or size does not match the approved scope.")
}
if (anyNA(prior$source) || anyNA(prior$target) || any(!prior$mor %in% c(-1, 1)) ||
    anyDuplicated(prior[, .(source, target)]) > 0L) {
  stop("The frozen CollecTRI network has missing/invalid signs or duplicate TF-target pairs.")
}
prior[, likelihood := 1]

available_mib <- suppressWarnings(as.numeric(Sys.getenv("FIG4C_AVAILABLE_MIB", NA_character_)))
if (length(available_mib) != 1L || !is.finite(available_mib) || available_mib <= 0) {
  stop("PowerShell must pass current available memory as FIG4C_AVAILABLE_MIB; analysis was not started.")
}
available_gib <- available_mib / 1024
if (resume_from_checkpoints && available_gib < 4) {
  stop("Checkpoint-resume mode requires at least 4 GiB available system memory.")
}
manifest[, CellID := as.character(CellID)]
manifest[, Cluster := as.character(Cluster)]
cluster_cell_counts <- manifest[, .(n_cells = .N), by = Cluster]
cluster_cell_counts[, Cluster := factor(Cluster, levels = expected_order)]
setorder(cluster_cell_counts, Cluster)

if (resume_from_checkpoints) {
  checkpoint_paths <- file.path(checkpoint_dir, paste0("cluster_mean_", expected_order, ".rds"))
  if (!all(file.exists(checkpoint_paths))) {
    stop("Checkpoint resume requires one checkpoint for each of the 24 frozen clusters.")
  }
  checkpoint_files <- list.files(checkpoint_dir, pattern = "^cluster_mean_.*[.]rds$", full.names = TRUE)
  if (length(checkpoint_files) != length(expected_order) || !setequal(checkpoint_files, checkpoint_paths)) {
    stop("Checkpoint directory must contain exactly the 24 expected cluster-mean files.")
  }
  checkpoint_records <- lapply(checkpoint_paths, readRDS)
  if (any(!vapply(checkpoint_records, function(x) {
    is.list(x) && all(c("cluster", "n_cells", "mean_log_normalized_expression") %in% names(x))
  }, logical(1)))) {
    stop("A checkpoint is missing required identity or expression-mean fields.")
  }
  if (!identical(vapply(checkpoint_records, `[[`, character(1), "cluster"), expected_order)) {
    stop("Checkpoint cluster IDs or ordering differ from the frozen 24-cluster order.")
  }
  checkpoint_names <- lapply(checkpoint_records, function(x) names(x$mean_log_normalized_expression))
  if (!all(vapply(checkpoint_records, function(x) {
    is.numeric(x$mean_log_normalized_expression) &&
      length(x$mean_log_normalized_expression) > 0L &&
      all(is.finite(x$mean_log_normalized_expression)) &&
      !is.null(names(x$mean_log_normalized_expression))
  }, logical(1)))) {
    stop("A checkpoint contains missing names or non-finite expression means.")
  }
  if (!all(vapply(checkpoint_names[-1L], identical, logical(1), checkpoint_names[[1L]])) ||
      anyDuplicated(checkpoint_names[[1L]]) || anyNA(checkpoint_names[[1L]]) ||
      any(!nzchar(checkpoint_names[[1L]]))) {
    stop("Checkpoint gene identifiers or order are inconsistent across clusters.")
  }
  checkpoint_counts <- vapply(checkpoint_records, `[[`, integer(1), "n_cells")
  expected_counts <- cluster_cell_counts$n_cells[match(expected_order, as.character(cluster_cell_counts$Cluster))]
  if (!identical(as.integer(checkpoint_counts), as.integer(expected_counts)) ||
      sum(checkpoint_counts) != nrow(manifest)) {
    stop("Checkpoint cell counts differ from the frozen clean24 manifest.")
  }
  features <- checkpoint_names[[1L]]
  cluster_means <- do.call(cbind, lapply(checkpoint_records, `[[`, "mean_log_normalized_expression"))
  dimnames(cluster_means) <- list(features, expected_order)
  if (nrow(cluster_means) != 18687L || ncol(cluster_means) != 24L || any(!is.finite(cluster_means))) {
    stop("Resumed checkpoint matrix dimensions or numerical integrity do not match the audited 18,687 x 24 matrix.")
  }
  cluster_map <- NULL
  fwrite(data.table(
    resume_mode = "cluster_mean_checkpoints", formal_object_read = FALSE,
    expression_aggregation_repeated = FALSE, checkpoint_count = length(checkpoint_records),
    cell_count = sum(checkpoint_counts), feature_count = length(features),
    finite_mean_values = sum(is.finite(cluster_means)),
    available_system_memory_GiB = available_gib, status = "validated_and_resumed"
  ), file.path(output_dir, "checkpoint_resume_audit.csv"))
  fwrite(data.table(
    Cluster = expected_order, n_cells = as.integer(checkpoint_counts),
    evidence = "previously generated cluster mean checkpoint; no object or metadata reread"
  ), file.path(output_dir, "checkpoint_cluster_identity.csv"))
  rm(checkpoint_records, checkpoint_names)
  gc(verbose = FALSE)
  message("Resumed validated expression means from 24 checkpoints; the formal Seurat object was not read.")
} else {
  serialized_gib <- file.info(object_path)$size / 1024^3
  minimum_headroom_gib <- max(12, 2 * serialized_gib + 6)
  resource_preflight <- data.table(
    formal_RDS_serialized_GiB = serialized_gib,
    available_system_memory_GiB = available_gib,
    minimum_start_headroom_GiB = minimum_headroom_gib,
    matched_expression_cache = FALSE,
    sparse_cluster_block_aggregation = TRUE,
    status = if (available_gib >= minimum_headroom_gib) "passed_preload_guard" else "blocked_before_object_load"
  )
  fwrite(resource_preflight, file.path(output_dir, "resource_preflight.csv"))
  if (available_gib < minimum_headroom_gib) {
    stop("Available memory is below the conservative object-load headroom; inspect cache/chunking rather than loading the formal object.")
  }
  message("Loading formal Seurat object; serialized GiB: ", round(serialized_gib, 2),
          "; available GiB: ", round(available_gib, 2),
          "; required headroom GiB: ", round(minimum_headroom_gib, 2))
  gc(verbose = FALSE)
  obj <- readRDS(object_path)
  if (!inherits(obj, "Seurat")) stop("The formal RDS did not deserialize to a Seurat object.")
  if (!"RNA" %in% Assays(obj)) stop("The formal Seurat object lacks the RNA assay.")
  meta <- obj@meta.data
  if (!all(c("RNA_snn_res.0.6", "seurat_clusters") %in% names(meta))) {
    stop("The formal object lacks the frozen res0.6 cluster metadata fields.")
  }
  if (!identical(as.character(meta$RNA_snn_res.0.6), as.character(meta$seurat_clusters))) {
    stop("RNA_snn_res.0.6 and seurat_clusters disagree in the formal object.")
  }
  rna_assay <- obj[["RNA"]]
  rna_layers <- Layers(rna_assay)
  if (!"data" %in% rna_layers) {
    stop("The formal RNA assay does not contain a single explicit data layer; no layer joining was attempted.")
  }
  rna <- GetAssayData(obj, assay = "RNA", layer = "data")
  if (!inherits(rna, "sparseMatrix")) {
    stop("RNA data is not sparse; refusing to create or process a dense cell-level matrix.")
  }
  if (anyDuplicated(rownames(rna)) || anyNA(rownames(rna)) || any(!nzchar(rownames(rna)))) {
    stop("RNA feature identifiers must be unique, present, and non-empty.")
  }
  object_cell_index <- match(manifest$CellID, colnames(rna))
  meta_cell_index <- match(manifest$CellID, rownames(meta))
  if (anyNA(object_cell_index) || anyNA(meta_cell_index)) {
    stop("At least one frozen clean24 CellID is absent from the current formal object/RNA data layer.")
  }
  raw_object_cluster <- as.character(meta$RNA_snn_res.0.6[meta_cell_index])
  cluster_map <- unique(data.table(
    raw_RNA_snn_res0.6 = raw_object_cluster,
    frozen_cluster = manifest$Cluster
  ))
  if (nrow(cluster_map) != length(expected_order) ||
      anyDuplicated(cluster_map$raw_RNA_snn_res0.6) > 0L ||
      anyDuplicated(cluster_map$frozen_cluster) > 0L ||
      !setequal(cluster_map$frozen_cluster, expected_order)) {
    stop("Frozen CellID-to-object res0.6 mapping is not one-to-one across the 24 approved clusters.")
  }
  manifest[, object_column := object_cell_index]
  rm(meta, rna_assay, obj)
  gc(verbose = FALSE)

  features <- rownames(rna)
  cluster_means <- matrix(
    NA_real_, nrow = length(features), ncol = length(expected_order),
    dimnames = list(features, expected_order)
  )
  for (cluster_id in expected_order) {
    rows <- manifest[Cluster == cluster_id]
    if (nrow(rows) < 1L) stop("An approved cluster has zero frozen cells: ", cluster_id)
    cell_columns <- rows$object_column
    block <- rna[, cell_columns, drop = FALSE]
    cluster_means[, cluster_id] <- as.numeric(Matrix::rowSums(block) / nrow(rows))
    saveRDS(
      list(cluster = cluster_id, n_cells = nrow(rows), mean_log_normalized_expression = cluster_means[, cluster_id]),
      file.path(checkpoint_dir, paste0("cluster_mean_", cluster_id, ".rds")),
      compress = FALSE
    )
    rm(block, rows, cell_columns)
    gc(verbose = FALSE)
    message("Completed sparse cell-equal mean for ", cluster_id)
  }
}

finite_gene <- apply(cluster_means, 1L, function(x) all(is.finite(x)))
centered_means <- matrix(NA_real_, nrow = nrow(cluster_means), ncol = ncol(cluster_means),
                         dimnames = dimnames(cluster_means))
if (any(finite_gene)) {
  gene_background_mean <- rowMeans(cluster_means[finite_gene, , drop = FALSE])
  centered_means[finite_gene, ] <- cluster_means[finite_gene, , drop = FALSE] - gene_background_mean
} else {
  stop("No finite RNA gene means remain across all 24 clusters.")
}
finite_indices <- which(finite_gene)
raw_finite <- cluster_means[finite_gene, , drop = FALSE]
centered_finite <- centered_means[finite_gene, , drop = FALSE]
max_abs <- apply(abs(raw_finite), 1L, max)
max_centered_abs <- apply(abs(centered_finite), 1L, max)
zero_tol <- .Machine$double.eps * 16 * pmax(1, max_abs)
zero_across_clusters <- max_abs <= zero_tol
constant_across_clusters <- !zero_across_clusters & max_centered_abs <= zero_tol
variable_across_clusters <- !zero_across_clusters & !constant_across_clusters
gene_status <- data.table(
  feature_id = features,
  finite_in_all_clusters = finite_gene,
  all_zero_across_clusters = FALSE,
  no_between_cluster_variation_nonzero = FALSE,
  variable_finite = FALSE,
  included_in_24_cluster_matrix = finite_gene,
  included_as_informative_target = FALSE,
  exclusion_or_status = "nonfinite_in_at_least_one_cluster"
)
gene_status[finite_gene, all_zero_across_clusters := zero_across_clusters]
gene_status[finite_gene, no_between_cluster_variation_nonzero := constant_across_clusters]
gene_status[finite_gene, variable_finite := variable_across_clusters]
gene_status[finite_gene, exclusion_or_status := "retained_in_common_finite_background"]
gene_status[finite_gene & all_zero_across_clusters, exclusion_or_status := "all_zero; retained_as_zero_row_in_common_background"]
gene_status[finite_gene & no_between_cluster_variation_nonzero, exclusion_or_status := "no_between_cluster_variation; retained_as_centered_zero_row"]
gene_status[finite_gene & variable_finite, included_as_informative_target := TRUE]

valid_features <- features[finite_gene]
informative_features <- gene_status[variable_finite == TRUE, feature_id]
effective_network <- prior[target %in% informative_features]
prior_tf_targets <- prior[, .(prior_unique_target_n = uniqueN(target)), by = source]
finite_target_counts <- prior[target %in% valid_features,
  .(finite_RNA_target_n = uniqueN(target)), by = source]
informative_target_counts <- prior[target %in% informative_features,
  .(effective_target_n = uniqueN(target)), by = source]
coverage <- merge(prior_tf_targets, finite_target_counts, by = "source", all.x = TRUE)
coverage <- merge(coverage, informative_target_counts, by = "source", all.x = TRUE)
coverage[is.na(finite_RNA_target_n), finite_RNA_target_n := 0L]
coverage[is.na(effective_target_n), effective_target_n := 0L]
coverage[, finite_RNA_coverage := finite_RNA_target_n / prior_unique_target_n]
coverage[, effective_target_coverage := effective_target_n / prior_unique_target_n]
coverage[, eligible_min10 := effective_target_n >= 10L]
coverage[, exclusion_reason := fifelse(eligible_min10, "",
  fifelse(effective_target_n == 0L, "no informative prior targets in RNA feature background",
          "fewer than 10 informative prior targets in RNA feature background"))]
eligible_tfs <- coverage[eligible_min10 == TRUE, source]
if (length(eligible_tfs) == 0L) stop("No TF meets the user-specified 10-target threshold.")
effective_network <- effective_network[source %in% eligible_tfs]
if (!setequal(unique(effective_network$source), eligible_tfs)) {
  stop("Filtered effective network TF list differs from the audited eligible TF list.")
}

ulm_matrix <- centered_means[valid_features, , drop = FALSE]
storage.mode(ulm_matrix) <- "double"
if (any(!is.finite(ulm_matrix))) stop("The common ULM gene-by-cluster matrix contains non-finite values.")
score_long <- decoupleR::run_ulm(
  mat = ulm_matrix,
  network = as.data.frame(effective_network),
  .source = "source",
  .target = "target",
  .mor = "mor",
  .likelihood = "likelihood",
  minsize = 10L,
  center = FALSE,
  na.rm = FALSE,
  sparse = FALSE
)
score_long <- as.data.table(score_long)
if (!all(c("source", "condition", "score", "p_value", "statistic") %in% names(score_long))) {
  stop("run_ulm returned an unexpected result schema.")
}
if (anyDuplicated(score_long[, .(source, condition)]) > 0L) {
  stop("run_ulm returned duplicate TF-cluster activity rows.")
}
score_long[, condition := as.character(condition)]
score_long[, source := as.character(source)]
if (!setequal(unique(score_long$condition), expected_order)) {
  stop("run_ulm cluster conditions differ from the frozen 24-cluster input matrix.")
}
score_stats <- score_long[, .(
  estimated_cluster_n = uniqueN(condition),
  all_24_finite = uniqueN(condition) == 24L && all(is.finite(score)),
  raw_score_mean = mean(score[is.finite(score)]),
  raw_score_sd = sd(score[is.finite(score)]),
  raw_score_IQR = IQR(score[is.finite(score)])
), by = source]
coverage <- merge(coverage, score_stats, by = "source", all.x = TRUE)
coverage[is.na(estimated_cluster_n), estimated_cluster_n := 0L]
coverage[is.na(all_24_finite), all_24_finite := FALSE]
coverage[, has_between_cluster_score_variation := is.finite(raw_score_sd) & raw_score_sd > 0]
coverage[, display_eligible := eligible_min10 & all_24_finite & has_between_cluster_score_variation]
coverage[, display_exclusion_reason := fifelse(!eligible_min10, exclusion_reason,
  fifelse(!all_24_finite, "ULM score missing/nonfinite in one or more clusters",
          fifelse(!has_between_cluster_score_variation, "ULM score has no between-cluster variation", "")))]

all_tf_cluster_grid <- CJ(source = coverage$source, condition = expected_order, unique = TRUE)
all_tf_scores <- merge(
  all_tf_cluster_grid,
  score_long[, .(source, condition, statistic, score, p_value)],
  by = c("source", "condition"), all.x = TRUE, sort = FALSE
)
all_tf_scores <- merge(
  all_tf_scores,
  coverage[, .(source, eligible_min10, exclusion_reason, display_exclusion_reason)],
  by = "source", all.x = TRUE, sort = FALSE
)
all_tf_scores[, score_status := fifelse(
  !eligible_min10, "not_run_below_10_informative_targets",
  fifelse(is.finite(score), "scored", "missing_or_nonfinite_from_run_ulm")
)]

ranking <- coverage[display_eligible == TRUE][order(-raw_score_IQR, -raw_score_sd, source)]
ranking[, display_rank := .I]
ranking[, selected_for_figure := display_rank <= 20L]
coverage <- merge(coverage, ranking[, .(source, display_rank, selected_for_figure)], by = "source", all.x = TRUE)
coverage[is.na(selected_for_figure), selected_for_figure := FALSE]
coverage[is.na(display_rank), display_rank := NA_integer_]
selected_tfs <- ranking[selected_for_figure == TRUE, source]
if (length(selected_tfs) == 0L) stop("No complete, variable TF score is available for the figure.")

selected_long <- score_long[source %in% selected_tfs]
selected_long[, Cluster := factor(condition, levels = expected_order)]
selected_long[, TF := factor(source, levels = selected_tfs)]
selected_long[, score_mean := mean(score), by = source]
selected_long[, score_sd := sd(score), by = source]
selected_long[, z_score := (score - score_mean) / score_sd]
selected_long[, TF_display_rank := match(source, selected_tfs)]
selected_long[, cluster_display_order := match(condition, expected_order)]
if (any(!is.finite(selected_long$z_score))) stop("Selected TF Z-score calculation returned non-finite values.")

write_matrix_csv <- function(x, path) {
  out <- as.data.table(x, keep.rownames = "feature_id")
  fwrite(out, path, na = "NA")
}
write_matrix_csv(cluster_means, file.path(output_dir, "cluster_mean_log_normalized_expression_wide.csv"))
write_matrix_csv(centered_means, file.path(output_dir, "cluster_mean_expression_centered_wide.csv"))
fwrite(cluster_cell_counts, file.path(output_dir, "cluster_cell_counts.csv"))
if (!is.null(cluster_map)) {
  fwrite(cluster_map, file.path(output_dir, "object_to_frozen_cluster_mapping.csv"))
}
fwrite(gene_status, file.path(output_dir, "RNA_gene_background_status.csv"))
fwrite(coverage[order(-effective_target_n, source)], file.path(output_dir, "TF_target_coverage_and_ULM_status.csv"))
fwrite(all_tf_scores, file.path(output_dir, "TF_activity_all_raw_long.csv"))
fwrite(coverage[display_eligible == TRUE][order(display_rank)], file.path(output_dir, "TF_display_ranking_and_selection.csv"))
fwrite(selected_long[order(cluster_display_order, TF_display_rank)],
       file.path(output_dir, "Fig4C_cluster_by_TF_source.csv"))
fwrite(order_table[order(row_top_to_bottom)], file.path(output_dir, "cluster_display_order_used.csv"))
fwrite(anchor_geometry[order(row_top_to_bottom)], file.path(output_dir, "anchor_row_center_geometry_used.csv"))

scope_map <- unique(scope_table[, .(Cluster, analysis_scope)])
scope_map <- scope_map[Cluster %in% expected_order]
scope_map[, analysis_scope := as.character(analysis_scope)]
scope_colors <- c(Shared6 = "#44718C", Restricted4 = "#B87555", pan14 = "#CFCFCF")
if (nrow(scope_map) != 24L || anyNA(scope_map$analysis_scope) ||
    any(!scope_map$analysis_scope %in% names(scope_colors))) {
  stop("Frozen cluster group annotation is incomplete or contains an unknown group.")
}
order_scope_map <- order_table[order(row_top_to_bottom), .(Cluster, anchor_scope = as.character(scope_group))]
scope_check <- merge(order_scope_map, scope_map, by = "Cluster", all = TRUE, sort = FALSE)
if (nrow(scope_check) != 24L || anyNA(scope_check$anchor_scope) || anyNA(scope_check$analysis_scope) ||
    any(scope_check$anchor_scope != scope_check$analysis_scope)) {
  stop("FUN8 row-order scope labels disagree with the frozen cluster_scope mapping.")
}

anchor_geometry <- anchor_geometry[order(row_top_to_bottom)]
plot_height_mm <- unique(anchor_geometry$anchor_page_height_mm)
if (length(plot_height_mm) != 1L || !is.finite(plot_height_mm) || plot_height_mm != 490) {
  stop("The anchor geometry height is not the measured FUN8 page height of 490 mm.")
}
if (nrow(anchor_geometry) != 24L || any(!is.finite(anchor_geometry$center_mm_from_top)) ||
    any(!is.finite(anchor_geometry$tile_height_mm)) || any(anchor_geometry$tile_height_mm <= 0) ||
    anyDuplicated(anchor_geometry$center_mm_from_top) > 0L ||
    any(diff(anchor_geometry$center_mm_from_top) <= 0) ||
    any(anchor_geometry$center_mm_from_top - anchor_geometry$tile_height_mm / 2 < 0) ||
    any(anchor_geometry$center_mm_from_top + anchor_geometry$tile_height_mm / 2 > plot_height_mm)) {
  stop("Measured anchor row-center geometry is invalid or out of page bounds.")
}
row_centres_from_top <- setNames(anchor_geometry$center_mm_from_top, anchor_geometry$Cluster)
row_centres_from_bottom <- plot_height_mm - row_centres_from_top
row_tile_heights_mm <- setNames(anchor_geometry$tile_height_mm, anchor_geometry$Cluster)
plot_width_mm <- max(210, 52 + 10.8 * length(selected_tfs))
heatmap_left_mm <- 45
heatmap_right_mm <- 8
heatmap_width_mm <- plot_width_mm - heatmap_left_mm - heatmap_right_mm
tf_cell_width_mm <- heatmap_width_mm / length(selected_tfs)
tf_x_centres <- heatmap_left_mm + (seq_along(selected_tfs) - 0.5) * tf_cell_width_mm
z_limit <- max(abs(selected_long$z_score))
if (!is.finite(z_limit) || z_limit <= 0) stop("Invalid Z-score display range.")

layout_audit <- data.table(
  row_top_to_bottom = seq_along(expected_order),
  Cluster = expected_order,
  anchor_center_mm_from_top = as.numeric(row_centres_from_top[expected_order]),
  Fig4C_center_mm_from_top = as.numeric(row_centres_from_top[expected_order]),
  center_delta_mm = 0,
  tile_height_mm = as.numeric(row_tile_heights_mm[expected_order]),
  anchor_source = "FUN8 anchor PDF raster geometry; 1500px rendering, 490mm page height"
)
fwrite(layout_audit, file.path(output_dir, "Fig4C_anchor_row_center_alignment.csv"))

map_z_colour <- function(z) {
  result <- rep("#BDBDBD", length(z))
  valid <- is.finite(z)
  if (any(valid)) {
    position <- pmin(1, pmax(0, (z[valid] / z_limit + 1) / 2))
    rgb_values <- grDevices::colorRamp(
      c("#2166AC", "#F7F7F7", "#B2182B"), space = "Lab"
    )(position)
    result[valid] <- grDevices::rgb(
      rgb_values[, 1], rgb_values[, 2], rgb_values[, 3], maxColorValue = 255
    )
  }
  result
}

draw_panel <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    xscale = c(0, plot_width_mm), yscale = c(0, plot_height_mm),
    width = grid::unit(1, "npc"), height = grid::unit(1, "npc"), clip = "off"
  ))
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  grid::grid.text("C. Inferred TF activity across 24 clusters",
                  x = grid::unit(9, "mm"), y = grid::unit(486, "mm"),
                  just = c("left", "centre"), gp = grid::gpar(fontfamily = "Arial", fontsize = 10, fontface = "bold", col = "black"))
  grid::grid.text("Centered cluster-mean RNA expression; ULM scores are relative to the 24-cluster background",
                  x = grid::unit(9, "mm"), y = grid::unit(478, "mm"),
                  just = c("left", "centre"), gp = grid::gpar(fontfamily = "Arial", fontsize = 7.2, col = "#333333"))

  strip_x <- heatmap_left_mm - 5
  for (i in seq_len(nrow(anchor_geometry))) {
    cl <- anchor_geometry$Cluster[i]
    y <- row_centres_from_bottom[[cl]]
    group <- scope_map$analysis_scope[match(cl, scope_map$Cluster)]
    grid::grid.rect(
      x = grid::unit(strip_x, "mm"), y = grid::unit(y, "mm"),
      width = grid::unit(3.7, "mm"), height = grid::unit(row_tile_heights_mm[[cl]], "mm"),
      gp = grid::gpar(fill = unname(scope_colors[group]), col = NA)
    )
    grid::grid.text(
      cl, x = grid::unit(heatmap_left_mm - 8, "mm"), y = grid::unit(y, "mm"),
      just = c("right", "centre"), gp = grid::gpar(fontfamily = "Arial", fontsize = 7.2, col = "black")
    )
  }

  selected_long[, draw_colour := map_z_colour(z_score)]
  for (i in seq_len(nrow(selected_long))) {
    row <- selected_long[i]
    x <- tf_x_centres[row$TF_display_rank]
    y <- row_centres_from_bottom[[as.character(row$Cluster)]]
    grid::grid.rect(
      x = grid::unit(x, "mm"), y = grid::unit(y, "mm"),
      width = grid::unit(tf_cell_width_mm * 0.94, "mm"),
      height = grid::unit(row_tile_heights_mm[[as.character(row$Cluster)]], "mm"),
      gp = grid::gpar(fill = row$draw_colour, col = "white", lwd = 0.5)
    )
  }
  for (boundary in c(6L, 10L)) {
    top_y <- mean(row_centres_from_top[expected_order[c(boundary, boundary + 1L)]])
    line_y <- plot_height_mm - top_y
    grid::grid.lines(
      x = grid::unit(c(heatmap_left_mm, plot_width_mm - heatmap_right_mm), "mm"),
      y = grid::unit(c(line_y, line_y), "mm"),
      gp = grid::gpar(col = "#555555", lwd = 0.6)
    )
  }

  label_y <- plot_height_mm - max(row_centres_from_top) - max(row_tile_heights_mm) / 2 - 7
  for (i in seq_along(selected_tfs)) {
    tf <- selected_tfs[i]
    grid::grid.text(
      tf, x = grid::unit(tf_x_centres[i], "mm"), y = grid::unit(label_y, "mm"),
      rot = 45, just = c("right", "centre"),
      gp = grid::gpar(fontfamily = "Arial", fontsize = 7, col = "black")
    )
  }

  cb_width <- min(64, plot_width_mm * 0.28)
  cb_left <- plot_width_mm / 2 - cb_width / 2
  cb_bottom <- 15.5
  cb_height <- 2.5
  n_steps <- 101L
  ramp_colours <- map_z_colour(seq(-z_limit, z_limit, length.out = n_steps))
  for (i in seq_len(n_steps)) {
    grid::grid.rect(
      x = grid::unit(cb_left + (i - 0.5) * cb_width / n_steps, "mm"),
      y = grid::unit(cb_bottom + cb_height / 2, "mm"),
      width = grid::unit(cb_width / n_steps + 0.02, "mm"),
      height = grid::unit(cb_height, "mm"),
      gp = grid::gpar(fill = ramp_colours[i], col = NA)
    )
  }
  grid::grid.text("Relative TF activity (Z-score)",
                  x = grid::unit(plot_width_mm / 2, "mm"), y = grid::unit(21.5, "mm"),
                  gp = grid::gpar(fontfamily = "Arial", fontsize = 7.2, col = "black"))
  cb_breaks <- c(-z_limit, 0, z_limit)
  cb_x <- cb_left + (cb_breaks + z_limit) / (2 * z_limit) * cb_width
  cb_labels <- c(formatC(-z_limit, digits = 2, format = "f"), "0", formatC(z_limit, digits = 2, format = "f"))
  for (i in seq_along(cb_breaks)) {
    grid::grid.text(cb_labels[i], x = grid::unit(cb_x[i], "mm"), y = grid::unit(12.3, "mm"),
                    gp = grid::gpar(fontfamily = "Arial", fontsize = 6.6, col = "black"))
  }

  group_order <- c("Shared6", "Restricted4", "pan14")
  group_spacing <- 36
  legend_start <- plot_width_mm / 2 - group_spacing
  for (i in seq_along(group_order)) {
    group <- group_order[i]
    x <- legend_start + (i - 1) * group_spacing
    grid::grid.rect(x = grid::unit(x, "mm"), y = grid::unit(6.2, "mm"),
                    width = grid::unit(4, "mm"), height = grid::unit(3, "mm"),
                    gp = grid::gpar(fill = unname(scope_colors[group]), col = NA))
    grid::grid.text(group, x = grid::unit(x + 3, "mm"), y = grid::unit(6.2, "mm"),
                    just = c("left", "centre"),
                    gp = grid::gpar(fontfamily = "Arial", fontsize = 6.8, col = "black"))
  }
  grid::grid.text("Within-TF Z-scores show relative cluster patterns, not absolute activation.",
                  x = grid::unit(9, "mm"), y = grid::unit(2.2, "mm"), just = c("left", "centre"),
                  gp = grid::gpar(fontfamily = "Arial", fontsize = 6.2, col = "#333333"))
  grid::popViewport()
}

png_path <- file.path(output_dir, "Fig4C_cluster24_TF_activity.png")
pdf_path <- file.path(output_dir, "Fig4C_cluster24_TF_activity.pdf")
grDevices::png(
  filename = png_path, width = round(plot_width_mm / 25.4 * 600),
  height = round(plot_height_mm / 25.4 * 600), units = "px", res = 600,
  bg = "white", family = "Arial", type = "cairo"
)
draw_panel()
grDevices::dev.off()
grDevices::cairo_pdf(
  filename = pdf_path, width = plot_width_mm / 25.4,
  height = plot_height_mm / 25.4, family = "Arial", bg = "white"
)
draw_panel()
grDevices::dev.off()

package_versions <- data.table(
  package = c("R", "SeuratObject", "Matrix", "decoupleR", "data.table", "grid"),
  version = c(
    as.character(getRversion()),
    as.character(packageVersion("SeuratObject")),
    as.character(packageVersion("Matrix")),
    as.character(packageVersion("decoupleR")),
    as.character(packageVersion("data.table")),
    as.character(packageVersion("grid"))
  )
)
fwrite(package_versions, file.path(output_dir, "package_versions.csv"))
input_parameters <- data.table(
  item = c(
    "analysis_route", "formal_object", "RNA_assay_layer", "frozen_scope", "cell_count",
    "available_RNA_feature_count", "cluster_mean_definition", "centering",
    "score_method", "decoupleR_run_ulm_center", "minimum_effective_targets",
    "ULM_sparse", "ULM_na_rm", "TF_display_rule", "Figure4A_anchor_row_order",
    "group_boundaries_after_rows", "figure_dimensions_mm", "figure_dpi", "interpretation_boundary"
  ),
  value = c(
    "cluster mean expression + human CollecTRI + decoupleR ULM",
    object_path, "RNA/data", "frozen clean24; Ctrl included", as.character(nrow(manifest)),
    as.character(length(features)), "Matrix rowSums / all cells in each frozen cluster; mean log-normalized expression",
    "gene-wise subtract equal-weight mean across the 24 cluster means; no SD scaling",
    "decoupleR::run_ulm", "FALSE; matrix was pre-centered", "10 informative signed target genes",
    "FALSE", "FALSE", "raw ULM IQR descending; ties SD descending, TF ascending; up to 20",
    paste(expected_order, collapse = ","), "6,10", paste(plot_width_mm, plot_height_mm, sep = " x "),
    "600", "descriptive cluster-level inferred relative activity; not causal or donor-level validation"
  )
)
fwrite(input_parameters, file.path(output_dir, "input_parameter_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
writeLines(c(
  "Renderer: grid vector primitives; editable PDF text via cairo_pdf.",
  paste0("Canvas: ", plot_width_mm, " x ", plot_height_mm, " mm."),
  "Cluster row centers use measured Fig4A FUN8 PDF row-center coordinates from anchor_row_center_geometry.csv.",
  paste0("Heatmap rows: ", length(expected_order), "; group boundaries after rows 6 and 10."),
  paste0("Displayed TFs: ", length(selected_tfs), "; each TF column follows raw ULM IQR rank.")
), file.path(output_dir, "plot_render_audit.txt"))

message("Finished Fig4C cluster-mean ULM pipeline. Features: ", length(features),
        "; finite background: ", sum(finite_gene), "; variable genes: ", sum(variable_across_clusters),
        "; eligible TFs: ", length(eligible_tfs), "; displayed TFs: ", length(selected_tfs),
        "; output: ", output_dir)
