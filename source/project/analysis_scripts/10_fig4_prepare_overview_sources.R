suppressPackageStartupMessages({
  library(data.table)
})
options(stringsAsFactors = FALSE)

ROOT <- "D:/data/microglia_new/analysis_scripts/runtime_io_aliases/fig3"
SCOPE <- file.path(ROOT, "tables/cluster_scope.csv")
FUN_SOURCE <- file.path(ROOT, "source_data/cluster24_FUN8_source.csv")
MP_SOURCE <- file.path(ROOT, "source_data/cluster24_MP12_source.csv")
ELIGIBILITY <- "E:/micro_output/fig3/fig4/source_data/Fig3A_clean24_sample_cluster_eligibility.csv"
SAMPLE_MAP <- "E:/micro_output/fig3/shared_function_framework/tables/disease_sample_coverage.csv"
OUT_FUN <- file.path(ROOT, "source_data/fig4_cluster24_FUN8_source.csv")
OUT_MP <- file.path(ROOT, "source_data/fig4_cluster24_MP12_source.csv")

required_files <- c(SCOPE, FUN_SOURCE, MP_SOURCE, ELIGIBILITY, SAMPLE_MAP)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing frozen source inputs: ", paste(missing_files, collapse = "; "))

scope <- fread(SCOPE)
if (!all(c("Cluster", "n_cells_total", "n_eligible_samples", "analysis_scope", "scope_order", "row_order") %in% names(scope))) {
  stop("cluster_scope.csv does not contain the expected clean24 identity fields")
}
if (nrow(scope) != 24L || anyDuplicated(scope$Cluster) || sum(scope$n_cells_total) != 400026L) {
  stop("cluster_scope.csv does not match the frozen clean24 24-cluster scope")
}
scope[, Cluster := as.character(Cluster)]

fun <- fread(FUN_SOURCE)
if (nrow(fun) != 192L || anyDuplicated(fun[, .(Cluster, Function)]) ||
    !all(c("raw_mean", "display_z", "cluster_order", "function_order") %in% names(fun))) {
  stop("FUN8 source is not a complete unique 24 x 8 matrix")
}
fun <- merge(fun, scope[, .(Cluster, n_cells_total, n_eligible_samples, analysis_scope, scope_order)],
             by = "Cluster", all.x = TRUE, sort = FALSE)
if (nrow(fun) != 192L || anyNA(fun$n_cells_total)) stop("FUN8 source did not join to all clean24 clusters")
fun[, `:=`(n_cells = n_cells_total,
            z_unclipped = display_z,
            z_display = pmax(-3, pmin(3, display_z)),
            z_display_clipped = abs(display_z) > 3)]
setorder(fun, scope_order, cluster_order, function_order)
fwrite(fun, OUT_FUN, na = "")

mp <- fread(MP_SOURCE)
mp_ids <- paste0("MP", sprintf("%02d", 1:12))
if (nrow(mp) != 288L || anyDuplicated(mp[, .(Cluster, Function)]) ||
    !all(c("raw_mean", "display_z", "cluster_order", "function_order", "n_samples") %in% names(mp)) ||
    !setequal(unique(mp$Function), mp_ids)) {
  stop("MP12 source is not a complete unique 24 x 12 matrix with frozen MP identities")
}
elig <- fread(ELIGIBILITY)
if (!all(c("Cluster", "SampleID_Global", "n_cells", "Eligible") %in% names(elig)) || anyDuplicated(elig[, .(Cluster, SampleID_Global)])) {
  stop("Clean24 sample x cluster eligibility source has an unexpected schema or duplicate keys")
}
map <- fread(SAMPLE_MAP)
if (!all(c("SampleID_Global", "DatasetID") %in% names(map)) || anyDuplicated(map$SampleID_Global)) {
  stop("Frozen sample map lacks unique SampleID_Global to DatasetID mapping")
}
eligible_map <- merge(elig[Eligible == TRUE, .(Cluster, SampleID_Global)],
                      map[, .(SampleID_Global, DatasetID)], by = "SampleID_Global", all.x = TRUE, sort = FALSE)
if (anyNA(eligible_map$DatasetID)) stop("At least one eligible MP sample lacks DatasetID mapping")
eligible_coverage <- eligible_map[, .(n_eligible_samples_reconstructed = uniqueN(SampleID_Global),
                                      n_eligible_datasets = uniqueN(DatasetID)), by = Cluster]
mp <- merge(mp, scope[, .(Cluster, n_cells_total, analysis_scope, scope_order)],
            by = "Cluster", all.x = TRUE, sort = FALSE)
mp <- merge(mp, eligible_coverage, by = "Cluster", all.x = TRUE, sort = FALSE)
if (nrow(mp) != 288L || anyNA(mp$n_cells_total) || anyNA(mp$n_eligible_datasets) ||
    any(mp$n_samples != mp$n_eligible_samples_reconstructed)) {
  stop("MP12 source coverage disagrees with the frozen summary or clean24 cluster scope")
}
mp[, `:=`(z_unclipped = display_z,
           z_display = pmax(-3, pmin(3, display_z)),
           z_display_clipped = abs(display_z) > 3,
           n_cells_total = n_cells_total)]
setnames(mp, "n_samples", "n_eligible_samples")
setorder(mp, scope_order, cluster_order, function_order)
fwrite(mp, OUT_MP, na = "")

cat(sprintf("FUN8 source: %d rows; MP12 source: %d rows; clean24 cells: %d\n",
            nrow(fun), nrow(mp), sum(scope$n_cells_total)))
