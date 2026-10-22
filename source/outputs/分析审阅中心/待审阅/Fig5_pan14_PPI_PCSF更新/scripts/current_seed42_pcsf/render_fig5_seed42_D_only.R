#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 999)
suppressWarnings(invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8")))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
  library(jsonlite)
  library(scales)
})

ROOT <- "E:/micro_output/分析审阅中心/待审阅/Fig5_pan14_PPI_PCSF更新"
SRC <- file.path(ROOT, "source_data", "current_seed42_pcsf")
SRC_BASE <- file.path(ROOT, "source_data")
FIG <- file.path(ROOT, "figures", "current_seed42_pcsf")
LOG <- file.path(ROOT, "logs", "current_seed42_pcsf", "render_fig5_seed42_D_only.log")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG), recursive = TRUE, showWarnings = FALSE)

sha256_file <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) stop("Cannot hash missing file: ", path)
  out <- suppressWarnings(system2("certutil.exe", c("-hashfile", shQuote(normalizePath(path, winslash = "\\", mustWork = TRUE)), "SHA256"), stdout = TRUE, stderr = TRUE))
  h <- tolower(out[grepl("^[A-Fa-f0-9]{64}$", out, perl = TRUE)])
  if (length(h) != 1L) stop("SHA-256 could not be parsed for ", path)
  h
}
log_line <- function(x) cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", x, "\n"), file = LOG, append = TRUE)
write_parameters <- function(path, x) fwrite(data.table(parameter = names(x), value = as.character(unlist(x))), path, bom = TRUE)

required <- c(
  file.path(SRC, "seed42_model_evidence_76_registry.csv"),
  file.path(SRC_BASE, "all_gene_mapping_audit.csv"),
  file.path(ROOT, "data", "orgHs_mapping_resource_manifest.json"),
  file.path(ROOT, "logs", "current_seed42_pcsf", "D_display_name_diagnostic_summary.csv"),
  file.path(ROOT, "logs", "current_seed42_pcsf", "D_display_name_diagnostic_displayed_terms.csv"),
  file.path(ROOT, "source_data", "current_seed42_pcsf", "QuickGO_134_missing_terms_metadata_20261022.json")
)
if (!all(file.exists(required))) stop("Missing D-only input: ", paste(required[!file.exists(required)], collapse = "; "))
out_png <- file.path(FIG, "Fig5D_seed42_PCSF_GO_BP_barplot.png")
out_pdf <- file.path(FIG, "Fig5D_seed42_PCSF_GO_BP_barplot.pdf")
if (file.exists(out_png) || file.exists(out_pdf)) stop("Refusing to overwrite existing D output.")
mapping_meta <- fromJSON(file.path(ROOT, "data", "orgHs_mapping_resource_manifest.json"), simplifyVector = TRUE)
if (as.character(packageVersion("org.Hs.eg.db")) != as.character(mapping_meta$version) ||
    as.character(packageVersion("org.Hs.eg.db")) != "3.19.1" ||
    as.character(packageVersion("GO.db")) != "3.20.0") {
  stop("Frozen GO mapping resources are not available at the accepted versions; no substitute will be used.")
}
registry <- fread(required[1], na.strings = c("", "NA"), encoding = "UTF-8")
mapping <- fread(required[2], na.strings = c("", "NA"), encoding = "UTF-8")
if (!all(c("gene_id", "core_10_of_10") %in% names(registry)) || anyDuplicated(registry$gene_id) || nrow(registry) != 76L) {
  stop("The frozen 76-gene seed42 model-evidence registry is invalid.")
}
if (!all(c("seed", "gene_id", "mapping_status", "chosen_Entrez", "in_qualified_PPI", "mapping_source") %in% names(mapping))) {
  stop("Frozen mapping audit fields do not match the accepted D input.")
}
map42 <- mapping[seed == 42L]
mapped_background <- unique(map42[mapping_status %chin% c("official_unique", "alias_unique") & in_qualified_PPI == TRUE,
  .(gene_id = as.character(gene_id), ENTREZID = as.character(chosen_Entrez))])
if (anyDuplicated(mapped_background$gene_id)) stop("GO background has duplicate model-gene mapping keys.")
entrez_keys <- unique(mapped_background$ENTREZID)
annotation_raw <- as.data.table(suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
  keys = entrez_keys, keytype = "ENTREZID", columns = c("GOALL", "ONTOLOGYALL"))))
bp_entrez <- unique(annotation_raw[ONTOLOGYALL == "BP" & !is.na(GOALL),
  .(ENTREZID = as.character(ENTREZID), GOID = as.character(GOALL))])
bp_map <- unique(merge(mapped_background, bp_entrez, by = "ENTREZID", all = FALSE, sort = FALSE)[, .(gene_id, GOID)])
candidate_ids <- unique(as.character(registry$gene_id))
candidate_annotated <- intersect(candidate_ids, unique(bp_map$gene_id))
background_annotated <- unique(bp_map$gene_id)
N_background <- length(background_annotated)
n_candidate <- length(candidate_annotated)
term_table <- bp_map[, .(K_term = uniqueN(gene_id)), by = GOID]
go_terms <- as.data.table(suppressMessages(AnnotationDbi::select(GO.db::GO.db,
  keys = term_table$GOID, keytype = "GOID", columns = "TERM")))
go_terms <- unique(go_terms[!is.na(TERM), .(GOID = as.character(GOID), term_name_local = as.character(TERM))])
term_table <- merge(term_table, go_terms, by = "GOID", all.x = TRUE, sort = FALSE)
candidate_map <- bp_map[gene_id %in% candidate_annotated,
  .(candidate_genes = paste(sort(unique(gene_id)), collapse = ";"), k_hit = uniqueN(gene_id)), by = GOID]
tests <- merge(term_table, candidate_map, by = "GOID", all.x = TRUE, sort = FALSE)
tests[is.na(k_hit), `:=`(k_hit = 0L, candidate_genes = "")]
tests[, `:=`(N_background = as.integer(N_background), n_candidate_annotated = as.integer(n_candidate),
  K_term = as.integer(K_term), k_hit = as.integer(k_hit))]
tests[, P := stats::phyper(k_hit - 1L, K_term, N_background - K_term, n_candidate_annotated, lower.tail = FALSE)]
tests[, q_BH_all_BP := stats::p.adjust(P, method = "BH")]
tests[, significant_BH_005 := q_BH_all_BP < 0.05]
setorder(tests, q_BH_all_BP, P, GOID)
tests[, `:=`(displayed_on_main_figure = FALSE, display_basis = "not selected for main figure")]
sig <- tests[significant_BH_005 == TRUE][order(q_BH_all_BP, P, GOID)]
if (nrow(sig)) {
  display_raw <- head(sig, 20L)
  mode <- "BH-significant GO BP terms"
  basis <- "Top 20 q<0.05 terms ordered by BH q, raw P, GOID; complete family retained"
} else {
  display_raw <- head(tests[k_hit > 0L][order(-k_hit, GOID)], 20L)
  mode <- "Descriptive GO BP terms; no term met BH q<0.05"
  basis <- "Top 20 by candidate hit count then GOID; complete family retained"
}
if (nrow(display_raw)) tests[GOID %in% display_raw$GOID, `:=`(displayed_on_main_figure = TRUE, display_basis = basis)]
tests[, neg_log10_BH_q := -log10(pmax(q_BH_all_BP, .Machine$double.xmin))]

# Save all statistical products before trying to resolve term names.
stats_path <- file.path(SRC, "Fig5D_seed42_PCSF_GO_BP_full_test_table_statistics_before_name_completion.csv")
fwrite(tests, stats_path, na = "NA", bom = TRUE)
candidate_mapping_path <- file.path(SRC, "Fig5D_seed42_PCSF_GO_BP_candidate_gene_term_mapping.csv")
fwrite(bp_map[gene_id %in% candidate_annotated], candidate_mapping_path, na = "NA", bom = TRUE)
candidate_hits_path <- file.path(SRC, "Fig5D_seed42_PCSF_GO_BP_candidate_term_hits.csv")
fwrite(candidate_map, candidate_hits_path, na = "NA", bom = TRUE)
top20_raw_path <- file.path(SRC, "Fig5D_seed42_PCSF_GO_BP_top20_raw_selection.csv")
display_raw[, display_rank := seq_len(.N)]
fwrite(display_raw, top20_raw_path, na = "NA", bom = TRUE)
background_audit <- copy(map42[, .(gene_id = as.character(gene_id), mapping_status, mapping_source,
  chosen_Entrez = as.character(chosen_Entrez), in_qualified_PPI)])
background_audit[, `:=`(in_GO_BP_background = gene_id %in% background_annotated,
  in_seed42_model_evidence_set = gene_id %in% candidate_ids,
  excluded_reason = fifelse(!in_qualified_PPI, "not_uniquely_mapped_to_qualified_PPI",
    fifelse(!gene_id %in% background_annotated, "no_GO_BP_annotation", "included_in_GO_BP_background")))]
background_path <- file.path(SRC, "Fig5D_seed42_PCSF_GO_BP_background_audit.csv")
fwrite(background_audit, background_path, na = "NA", bom = TRUE)
log_line(paste0("Saved complete statistical family before name resolution: candidates=", nrow(registry),
  "; annotated candidates=", n_candidate, "; background=", N_background, "; tested=", nrow(tests),
  "; significant=", sum(tests$significant_BH_005)))

# Official QuickGO metadata supplements only the missing display names; GO-to-gene relations and tests remain frozen.
quickgo_path <- required[6]
quickgo <- fromJSON(quickgo_path, simplifyVector = TRUE)
quickgo <- as.data.table(quickgo)
if (!all(c("id", "name", "aspect", "isObsolete", "comment") %in% names(quickgo)) || anyDuplicated(quickgo$id)) {
  stop("QuickGO metadata response is malformed.")
}
missing_ids <- unique(as.character(tests[is.na(term_name_local) | !nzchar(term_name_local), GOID]))
if (length(missing_ids) != 134L || !setequal(missing_ids, quickgo$id)) stop("QuickGO metadata does not exactly cover the 134 missing local GO.db names.")
if (any(quickgo$aspect != "biological_process")) stop("QuickGO aspect differs from the frozen GO BP test family.")
if (any(is.na(quickgo$name) | !nzchar(quickgo$name))) stop("Official QuickGO name is missing for at least one tested GO ID.")
replacement <- c("GO:0044271" = "GO:0009058", "GO:0019438" = "GO:0009058", "GO:0018130" = "GO:0009058")
name_map <- data.table(GOID = as.character(tests$GOID), term_name_local = tests$term_name_local)
name_map <- unique(name_map)
name_map[, quick_idx := match(GOID, quickgo$id)]
missing_local <- is.na(name_map$term_name_local) | !nzchar(name_map$term_name_local)
name_map[, `:=`(official_name = quickgo$name[quick_idx], ontology_aspect = quickgo$aspect[quick_idx],
  is_obsolete = as.logical(quickgo$isObsolete[quick_idx]), obsoletion_comment = quickgo$comment[quick_idx],
  source_url = ifelse(missing_local, paste0("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/", GOID), "https://bioconductor.org/packages/GO.db/"),
  source_path = ifelse(missing_local, quickgo_path, system.file(package = "GO.db")),
  source_type = ifelse(missing_local, "QuickGO official REST API", "Frozen local GO.db package"),
  resource_version = ifelse(missing_local, "Current API response; ontology release identifier not returned", "GO.db 3.20.0"),
  access_date = "2026-10-22")]
name_map[, display_name := fifelse(!is.na(term_name_local) & nzchar(term_name_local), term_name_local, official_name)]
name_map[, name_source := fifelse(!is.na(term_name_local) & nzchar(term_name_local), "Frozen local GO.db 3.20.0",
  "QuickGO official REST API; 2026-10-22")]
name_map[, term_status := fifelse(is.na(is_obsolete), "present in frozen GO.db 3.20.0; current status not separately reported",
  fifelse(is_obsolete, "obsolete", "current"))]
name_map[, official_replacement_id := unname(replacement[GOID])]
if (any(is.na(name_map$display_name) | !nzchar(name_map$display_name)) || anyDuplicated(name_map$GOID) || nrow(name_map) != nrow(tests)) {
  stop("GO ID display-name mapping is incomplete or non-unique.")
}
name_map_path <- file.path(SRC, "Fig5D_seed42_PCSF_GO_ID_display_name_mapping.csv")
fwrite(name_map, name_map_path, na = "NA", bom = TRUE)
name_idx <- match(tests$GOID, name_map$GOID)
tests[, `:=`(term_name = name_map$display_name[name_idx],
  name_source = name_map$name_source[name_idx], term_status = name_map$term_status[name_idx])]
if (anyNA(tests$term_name) || any(!nzchar(tests$term_name))) stop("An official display label remains unresolved.")

# Cross-check the previous diagnostic summary and exact diagnostic display order/statistics.
diag_summary <- fread(required[4], na.strings = c("", "NA"))
diag_display <- fread(required[5], na.strings = c("", "NA"))
get_diag <- function(key) as.integer(diag_summary[metric == key, value][1])
if (nrow(tests) != 11444L || get_diag("tested_terms") != 11444L || nrow(tests) != get_diag("tested_terms") || N_background != 4036L || n_candidate != 76L ||
    sum(tests$significant_BH_005) != get_diag("significant_terms") ||
    sum(is.na(name_map$term_name_local)) != get_diag("missing_official_term_names") ||
    nrow(display_raw) != get_diag("display_terms")) stop("Recomputed D counts disagree with the preserved diagnostic summary.")
display <- tests[match(display_raw$GOID, GOID)]
if (!identical(as.character(display$GOID), as.character(diag_display$GOID))) stop("Original Top20 term ordering differs from the preserved diagnostic.")
for (nm in c("P", "q_BH_all_BP")) {
  if (max(abs(display[[nm]] - diag_display[[nm]]), na.rm = TRUE) > 1e-12) stop("Recomputed Top20 ", nm, " differs from preserved diagnostic.")
}
if (sum(is.na(display_raw$term_name_local)) != get_diag("missing_display_names")) stop("Count of missing original Top20 local names differs from diagnostic.")
full_test_path <- file.path(SRC, "Fig5D_seed42_PCSF_GO_BP_full_test_table.csv")
fwrite(tests, full_test_path, na = "NA", bom = TRUE)
display[, `:=`(display_rank = seq_len(.N), term_label = term_name,
  term_label_wrapped = vapply(term_name, function(x) paste(strwrap(x, width = 38L), collapse = "\n"), character(1)))]
display_source_path <- file.path(SRC, "Fig5D_seed42_PCSF_GO_BP_barplot_source.csv")
fwrite(display, display_source_path, na = "NA", bom = TRUE)

if (anyNA(display$term_name) || any(!nzchar(display$term_name))) stop("Displayed GO term name unresolved; statistics are saved, drawing is blocked.")
display[, term_axis := factor(paste0("term_", display_rank), levels = rev(paste0("term_", display_rank)))]
axis_labels <- setNames(display$term_label_wrapped, as.character(display$term_axis))
pD <- ggplot(display, aes(x = neg_log10_BH_q, y = term_axis)) +
  geom_col(width = 0.7, fill = "#2166AC", colour = "#184F83", linewidth = 0.25) +
  scale_x_continuous(name = expression(-log[10]("BH-adjusted P")), breaks = scales::pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0.025))) +
  scale_y_discrete(name = NULL, labels = axis_labels, drop = FALSE) +
  labs(title = "GO Biological Process enrichment",
    caption = "Official term labels are shown as recorded; obsolete GO IDs retain their original statistics and are not replaced.") +
  theme_minimal(base_family = "Arial", base_size = 9) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    axis.text.y = element_text(colour = "black", size = 8, lineheight = 0.95),
    axis.text.x = element_text(colour = "black", size = 8),
    axis.title.x = element_text(colour = "black", size = 9, margin = margin(t = 7)),
    plot.title = element_text(colour = "black", size = 12, face = "plain"),
    plot.caption = element_text(colour = "#444444", size = 7.2, hjust = 0),
    plot.margin = margin(9, 12, 9, 9, unit = "pt"))
heightD <- max(90, 40 + 7.5 * nrow(display))
ggsave(out_png, pD, width = 265, height = heightD, units = "mm", dpi = 600, bg = "white", limitsize = FALSE)
ggsave(out_pdf, pD, width = 265, height = heightD, units = "mm", device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)

write_parameters(file.path(SRC, "Fig5D_seed42_PCSF_GO_BP_parameters.csv"), list(
  candidate_count = length(candidate_ids), candidate_annotated = n_candidate,
  background = "all pan14 model features uniquely mapped to qualified PPI Entrez nodes and annotated with GO BP before PCSF selection",
  background_genes = N_background, tested_terms = nrow(tests), significant_q_lt_005 = sum(tests$significant_BH_005),
  test = "one-sided hypergeometric over-representation; BH over all testable BP terms",
  mapping = "frozen org.Hs.eg.db 3.19.1 and GO.db 3.20.0; QuickGO metadata supplements only missing display names",
  display_mode = mode, display_basis = basis, display_missing_original_local_names = sum(is.na(display_raw$term_name_local)),
  obsolete_displayed_terms = paste(display[term_status == "obsolete", GOID], collapse = ";"),
  png_dpi = 600, font = "Arial", pdf_device = "cairo_pdf"))
input_paths <- c(required, file.path(ROOT, "scripts", "current_seed42_pcsf", "render_fig5_seed42_BCD.R"),
  file.path(ROOT, "scripts", "current_seed42_pcsf", "render_fig5_seed42_D_only.R"))
input_manifest <- data.table(path = input_paths, sha256 = vapply(input_paths, sha256_file, character(1)),
  size_bytes = file.info(input_paths)$size, role = c(rep("input_or_diagnostic", length(required)), "parent_renderer", "D_only_renderer"))
fwrite(input_manifest, file.path(SRC, "Fig5D_seed42_PCSF_GO_BP_input_sha256.csv"), na = "NA", bom = TRUE)
output_paths <- c(out_png, out_pdf, stats_path, full_test_path, background_path, candidate_mapping_path,
  candidate_hits_path, top20_raw_path, name_map_path, display_source_path,
  file.path(SRC, "Fig5D_seed42_PCSF_GO_BP_parameters.csv"), file.path(SRC, "Fig5D_seed42_PCSF_GO_BP_input_sha256.csv"))
if (!all(file.exists(output_paths))) stop("D-only output inventory is incomplete.")
fwrite(data.table(path = output_paths, sha256 = vapply(output_paths, sha256_file, character(1)),
  size_bytes = file.info(output_paths)$size), file.path(SRC, "Fig5D_seed42_PCSF_GO_BP_output_sha256.csv"), na = "NA", bom = TRUE)
log_line(paste0("D-only complete: candidates=", length(candidate_ids), "; background=", N_background,
  "; tests=", nrow(tests), "; BH q<0.05=", sum(tests$significant_BH_005),
  "; displayed=", nrow(display), "; displayed obsolete=", sum(display$term_status == "obsolete")))
