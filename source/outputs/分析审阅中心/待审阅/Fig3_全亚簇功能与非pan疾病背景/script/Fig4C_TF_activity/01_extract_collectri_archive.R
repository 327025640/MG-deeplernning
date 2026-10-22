input_file <- "D:/data/microglia_new/_staging_Fig4C_TF_activity_20261007/omnipath_interactions_archive_20260624.tsv.xz"
output_file <- "D:/data/microglia_new/_staging_Fig4C_TF_activity_20261007/collectri_human_raw_with_evidence_baseparser_20260624.tsv"
audit_file <- "D:/data/microglia_new/_staging_Fig4C_TF_activity_20261007/collectri_archive_selection_audit_baseparser_20260624.csv"

stopifnot(file.exists(input_file))
if (file.exists(output_file) || file.exists(audit_file)) {
  stop("Refusing to overwrite existing extraction outputs")
}

selected_columns <- c(
  "source", "target", "source_genesymbol", "target_genesymbol",
  "is_directed", "is_stimulation", "is_inhibition", "consensus_direction",
  "consensus_stimulation", "consensus_inhibition", "sources", "references",
  "curation_effort", "evidences", "collectri", "ncbi_tax_id_source",
  "entity_type_source", "ncbi_tax_id_target", "entity_type_target"
)

column_classes <- c(
  source = "character", target = "character",
  source_genesymbol = "character", target_genesymbol = "character",
  is_directed = "numeric", is_stimulation = "numeric", is_inhibition = "numeric",
  consensus_direction = "numeric", consensus_stimulation = "numeric",
  consensus_inhibition = "numeric", sources = "character",
  references = "character", curation_effort = "numeric",
  evidences = "character", collectri = "logical",
  ncbi_tax_id_source = "integer", entity_type_source = "character",
  ncbi_tax_id_target = "integer", entity_type_target = "character"
)

connection <- xzfile(input_file, open = "rt")
header <- strsplit(readLines(connection, n = 1L), "\t", fixed = TRUE)[[1L]]
if (!all(c(selected_columns, names(column_classes)) %in% header)) {
  close(connection)
  stop("Official archive snapshot has an unexpected schema")
}

col_types <- rep("NULL", length(header))
names(col_types) <- header
col_types[names(column_classes)] <- column_classes

counts <- c(total_rows = 0, collectri_rows = 0, human_collectri_rows = 0)
has_output <- FALSE
repeat {
  chunk <- read.delim(
    connection,
    header = FALSE,
    sep = "\t",
    quote = "\"",
    comment.char = "",
    col.names = header,
    colClasses = col_types,
    nrows = 100000L,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (nrow(chunk) == 0L) break
  counts[["total_rows"]] <- counts[["total_rows"]] + nrow(chunk)
  is_collectri <- !is.na(chunk$collectri) & chunk$collectri
  counts[["collectri_rows"]] <- counts[["collectri_rows"]] + sum(is_collectri)
  is_human <- !is.na(chunk$ncbi_tax_id_source) &
    !is.na(chunk$ncbi_tax_id_target) &
    chunk$ncbi_tax_id_source == 9606L &
    chunk$ncbi_tax_id_target == 9606L
  selected <- chunk[is_collectri & is_human, selected_columns, drop = FALSE]
  counts[["human_collectri_rows"]] <-
    counts[["human_collectri_rows"]] + nrow(selected)
  if (nrow(selected) > 0L) {
    write.table(
      selected,
      file = output_file,
      sep = "\t",
      quote = TRUE,
      qmethod = "double",
      row.names = FALSE,
      col.names = !has_output,
      append = has_output,
      na = ""
    )
    has_output <- TRUE
  }
}
close(connection)

audit <- data.frame(
  source_snapshot = "OmniPath legacy webservice archive, 2026-06-24",
  organism_tax_id = 9606L,
  total_rows = unname(counts[["total_rows"]]),
  collectri_rows = unname(counts[["collectri_rows"]]),
  human_collectri_rows = unname(counts[["human_collectri_rows"]])
)
utils::write.csv(audit, audit_file, row.names = FALSE)

cat(paste(names(counts), counts, sep = "=", collapse = "\n"), "\n")
cat("Selected edge file:", output_file, "\n")
