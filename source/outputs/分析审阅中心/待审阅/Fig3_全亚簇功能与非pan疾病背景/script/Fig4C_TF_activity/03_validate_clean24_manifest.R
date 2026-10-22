manifest_file <- "D:/data/microglia_new/_staging_Fig4C_TF_activity_20261007/cell_manifest_clean24.csv.gz"
cluster_file <- "D:/data/microglia_new/_staging_Fig4C_TF_activity_20261007/cluster_labels_res0.6.csv"
audit_file <- "D:/data/microglia_new/_staging_Fig4C_TF_activity_20261007/clean24_object_cluster_match_audit.csv"
cluster_count_file <- "D:/data/microglia_new/_staging_Fig4C_TF_activity_20261007/clean24_cluster_match_counts.csv"

stopifnot(file.exists(manifest_file), file.exists(cluster_file))
if (file.exists(audit_file) || file.exists(cluster_count_file)) {
  stop("Refusing to overwrite input-match audit outputs")
}

manifest <- readr::read_csv(manifest_file, show_col_types = FALSE, progress = FALSE)
object_labels <- readr::read_csv(cluster_file, show_col_types = FALSE, progress = FALSE)
required_manifest <- c("CellID", "Cluster", "analysis_scope", "Disease")
required_labels <- c("CellID", "cluster")
if (!all(required_manifest %in% names(manifest)) ||
    !all(required_labels %in% names(object_labels))) {
  stop("Unexpected schema in frozen manifest or object-derived cluster labels")
}

label_position <- match(manifest$CellID, object_labels$CellID)
observed_cluster <- sprintf("C%02d", as.integer(object_labels$cluster[label_position]))
cluster_counts <- as.data.frame(table(manifest$Cluster), stringsAsFactors = FALSE)
names(cluster_counts) <- c("Cluster", "frozen_manifest_cells")
matched_counts <- as.data.frame(table(observed_cluster), stringsAsFactors = FALSE)
names(matched_counts) <- c("Cluster", "object_label_matched_cells")
cluster_counts <- merge(cluster_counts, matched_counts, by = "Cluster", all = TRUE, sort = FALSE)
cluster_counts$frozen_manifest_cells[is.na(cluster_counts$frozen_manifest_cells)] <- 0L
cluster_counts$object_label_matched_cells[is.na(cluster_counts$object_label_matched_cells)] <- 0L
cluster_counts$exact_count_match <-
  cluster_counts$frozen_manifest_cells == cluster_counts$object_label_matched_cells

checks <- data.frame(
  check = c(
    "frozen_manifest_rows",
    "frozen_manifest_unique_CellID",
    "object_label_rows",
    "object_label_unique_CellID",
    "all_frozen_CellID_present_in_object_label_export",
    "all_24_frozen_cluster_labels_match_object_export",
    "all_frozen_cluster_counts_match_object_export",
    "Ctrl_included_in_frozen_manifest",
    "frozen_cluster_count"
  ),
  value = c(
    nrow(manifest),
    length(unique(manifest$CellID)),
    nrow(object_labels),
    length(unique(object_labels$CellID)),
    sum(!is.na(label_position)),
    sum(observed_cluster == manifest$Cluster, na.rm = TRUE),
    sum(cluster_counts$exact_count_match),
    any(manifest$Disease == "Ctrl", na.rm = TRUE),
    length(unique(manifest$Cluster))
  ),
  expected = c(
    400026, 400026, 424222, 424222, 400026, 400026, 24, TRUE, 24
  ),
  stringsAsFactors = FALSE
)
checks$pass <- as.character(checks$value) == as.character(checks$expected)
utils::write.csv(checks, audit_file, row.names = FALSE)
utils::write.csv(cluster_counts, cluster_count_file, row.names = FALSE)
cat(paste(checks$check, checks$value, sep = "=", collapse = "\n"), "\n")
