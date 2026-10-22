args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[[1]]
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})
orgdb_version <- as.character(packageVersion("org.Hs.eg.db"))
orgdb_metadata <- AnnotationDbi::metadata(org.Hs.eg.db)
official_keys <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")
alias_keys <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "ALIAS")
official <- AnnotationDbi::select(org.Hs.eg.db, keys = official_keys,
                                   keytype = "SYMBOL", columns = c("SYMBOL", "ENTREZID"))
official <- unique(official[!is.na(official$SYMBOL) & !is.na(official$ENTREZID), c("SYMBOL", "ENTREZID")])
aliases <- AnnotationDbi::select(org.Hs.eg.db, keys = alias_keys,
                                  keytype = "ALIAS", columns = c("ALIAS", "ENTREZID"))
aliases <- unique(aliases[!is.na(aliases$ALIAS) & !is.na(aliases$ENTREZID), c("ALIAS", "ENTREZID")])
utils::write.table(official, file.path(out_dir, "orgHs_official_symbol_entrez.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE, na = "")
utils::write.table(aliases, file.path(out_dir, "orgHs_alias_entrez.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE, na = "")
manifest <- list(resource = "Bioconductor org.Hs.eg.db", version = orgdb_version,
                 official_mapping_rows = nrow(official), alias_mapping_rows = nrow(aliases),
                 metadata = as.data.frame(orgdb_metadata, stringsAsFactors = FALSE))
jsonlite::write_json(manifest, file.path(out_dir, "orgHs_mapping_resource_manifest.json"),
                     auto_unbox = TRUE, pretty = TRUE, dataframe = "rows", na = "null")
cat(sprintf("org.Hs.eg.db %s: %d official mapping rows; %d alias rows\n",
            orgdb_version, nrow(official), nrow(aliases)))
