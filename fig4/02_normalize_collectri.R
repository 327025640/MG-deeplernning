.libPaths(c("D:/data/microglia_new/Rlibs/Fig4C_R4.4.3", .libPaths()))

input_file <- "D:/data/microglia_new/_staging_Fig4C_TF_activity_20261007/collectri_human_raw_with_evidence_baseparser_20260624.tsv"
output_file <- "D:/data/microglia_new/_staging_Fig4C_TF_activity_20261007/collectri_human_viper_regulon_20260624.csv"
audit_file <- "D:/data/microglia_new/_staging_Fig4C_TF_activity_20261007/collectri_normalization_audit_20260624.csv"
duplicate_file <- "D:/data/microglia_new/_staging_Fig4C_TF_activity_20261007/collectri_duplicate_direction_audit_20260624.csv"

stopifnot(file.exists(input_file))
if (file.exists(output_file) || file.exists(audit_file) || file.exists(duplicate_file)) {
  stop("Refusing to overwrite existing normalized prior outputs")
}

raw <- read.delim(
  input_file,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "\"",
  comment.char = ""
)
stopifnot(
  all(c("source", "target", "source_genesymbol", "target_genesymbol",
        "is_stimulation", "is_inhibition", "evidences") %in% names(raw)),
  all(raw$ncbi_tax_id_source == 9606L),
  all(raw$ncbi_tax_id_target == 9606L)
)

decode_evidences <- get("deserialize_evidences", asNamespace("OmnipathR"))
filter_to_dataset <- get("only_from", asNamespace("OmnipathR"))
raw <- decode_evidences(raw, simplifyVector = FALSE)
strict <- filter_to_dataset(raw, datasets = "collectri", .keep = FALSE)

complex_source <- !is.na(strict$source) & grepl("COMPLEX", strict$source, fixed = TRUE)
simple_edges <- strict[!complex_source, , drop = FALSE]
complex_edges <- strict[complex_source, , drop = FALSE]
complex_edges$source_genesymbol <- dplyr::case_when(
  grepl("JUN", complex_edges$source_genesymbol, fixed = TRUE) |
    grepl("FOS", complex_edges$source_genesymbol, fixed = TRUE) ~ "AP1",
  grepl("REL", complex_edges$source_genesymbol, fixed = TRUE) |
    grepl("NFKB", complex_edges$source_genesymbol, fixed = TRUE) ~ "NFKB",
  TRUE ~ complex_edges$source_genesymbol
)

signed_edges <- rbind(simple_edges, complex_edges)
signed_edges$mor <- ifelse(
  signed_edges$is_stimulation == 1,
  1,
  ifelse(signed_edges$is_stimulation == 0, -1, NA_real_)
)

missing_symbol <- is.na(signed_edges$source_genesymbol) |
  signed_edges$source_genesymbol == "" |
  is.na(signed_edges$target_genesymbol) |
  signed_edges$target_genesymbol == ""
usable_edges <- signed_edges[!missing_symbol & !is.na(signed_edges$mor), , drop = FALSE]

key <- paste(usable_edges$source_genesymbol, usable_edges$target_genesymbol, sep = "\r")
duplicate_key <- duplicated(key) | duplicated(key, fromLast = TRUE)
key_direction <- split(usable_edges$mor, key)
conflict_keys <- names(Filter(function(x) length(unique(x)) > 1L, key_direction))
conflict_edges <- key %in% conflict_keys
keep_first <- !duplicated(key)
regulon <- data.frame(
  source = usable_edges$source_genesymbol[keep_first],
  target = usable_edges$target_genesymbol[keep_first],
  mor = usable_edges$mor[keep_first],
  stringsAsFactors = FALSE
)

audit <- data.frame(
  stage = c(
    "archive_rows_before_strict_evidence_filter",
    "after_strict_collectri_evidence_filter",
    "complex_source_rows_collapsed_to_AP1_or_NFKB_or_preserved",
    "rows_missing_gene_symbol_or_direction",
    "duplicate_TF_target_rows_before_first_record_dedup",
    "TF_target_keys_with_conflicting_directions",
    "unique_signed_gene_symbol_edges_after_first_record_dedup"
  ),
  n = c(
    nrow(raw), nrow(strict), sum(complex_source), sum(missing_symbol),
    sum(duplicate_key), length(conflict_keys), nrow(regulon)
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(regulon, output_file, row.names = FALSE, na = "")
utils::write.csv(audit, audit_file, row.names = FALSE)
duplicate_audit <- usable_edges[duplicate_key, c(
  "source", "target", "source_genesymbol", "target_genesymbol",
  "is_stimulation", "is_inhibition", "mor", "sources", "references"
), drop = FALSE]
duplicate_audit$conflicting_direction_within_TF_target <- conflict_edges[duplicate_key]
utils::write.csv(duplicate_audit, duplicate_file, row.names = FALSE, na = "")
cat(paste(audit$stage, audit$n, sep = "=", collapse = "\n"), "\n")
cat("Regulon file:", output_file, "\n")
