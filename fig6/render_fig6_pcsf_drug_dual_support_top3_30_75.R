#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
  library(grid)
})
options(stringsAsFactors = FALSE, scipen = 999)
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"))

PROJECT_ROOT <- "E:/micro_output/分析审阅中心/待审阅/Fig6_疾病vsCtrl_SHAP_PCSF模块"
RUN_ROOT <- file.path(PROJECT_ROOT, "robustness/20261021_disease_compact_30_75_parameter_selection_01")
OUT <- file.path(RUN_ROOT, "drug_accessibility")
NETWORK_CSV <- file.path(OUT, "tables/all_disease_drug_network_proximity.csv")
NETWORK_ROOT <- file.path(RUN_ROOT, "selected_networks")
REGISTRY_CSV <- file.path(RUN_ROOT, "summaries/selected_network_delivery_registry.csv")
FINAL_SUMMARY_CSV <- file.path(RUN_ROOT, "summaries/final_disease_summary.csv")
SELECTED_LOCK_JSON <- file.path(RUN_ROOT, "config/selected_parameters.lock.json")
RUN_MANIFEST_JSON <- file.path(RUN_ROOT, "summaries/run_manifest.json")
RUN_CONFIG_JSON <- file.path(RUN_ROOT, "config/run_config.json")
INPUT_MANIFEST_CSV <- file.path(RUN_ROOT, "inputs/input_sha256_manifest.csv")
REVERSAL_CSV <- "E:/micro_output/分析审阅中心/待审阅/38核心疾病背景药物再定位/tables/v3_main_results_all_disease_drug.csv"
PALETTE_CSV <- "D:/data/microglia_new/authority/03_reference_configs/disease_palette_contract.csv"
LAYOUT_TEMPLATE <- "D:/data/microglia_new/analysis_scripts/render_Fig6_Top10_PPI_compact.R"
DUAL_RENDERER <- "E:/micro_output/分析审阅中心/待审阅/Fig6_疾病vsCtrl_SHAP_PCSF模块/drug_accessibility/scripts/render_fig6_pcsf_drug_dual_support_top3_12d.R"
ALL_DISEASES <- c("AD", "ALS", "ALS_FTD", "DLB", "DLBD_PD", "EOAD", "FTD", "FTD_MND", "GBM", "HD", "LCBM", "LGG", "MB", "MS", "PD", "PSP")
DISEASE_ORDER <- c("LGG", "GBM", "MB", "LCBM", "AD", "EOAD", "PD", "DLB",
                   "DLBD_PD", "ALS", "ALS_FTD", "FTD", "FTD_MND", "PSP", "HD", "MS")
MAX_DISPLAY <- 3L

FIG_PNG <- file.path(OUT, "figures/Fig6_PCSF_drug_dual_support_top3_30_75.png")
FIG_PDF <- file.path(OUT, "figures/Fig6_PCSF_drug_dual_support_top3_30_75.pdf")
CANDIDATE_OUT <- file.path(OUT, "tables/dual_support_candidates_top3_30_75_all.csv")
DISPLAY_OUT <- file.path(OUT, "tables/dual_support_display_top3_30_75.csv")
SUMMARY_OUT <- file.path(OUT, "tables/displayed_drug_summary_by_disease_dual_support_30_75.csv")
EXCLUDED_OUT <- file.path(OUT, "tables/excluded_diseases_dual_support_30_75.csv")
SOURCE_OUT <- file.path(OUT, "source_data/Fig6_PCSF_drug_dual_support_top3_30_75_source.csv")
KEY_AUDIT_OUT <- file.path(OUT, "tables/reversal_disease_drug_key_audit_30_75.csv")
HASH_OUT <- file.path(OUT, "config/Fig6_PCSF_drug_dual_support_top3_30_75_input_sha256.csv")
CONFIG_OUT <- file.path(OUT, "config/Fig6_PCSF_drug_dual_support_top3_30_75_config.json")
MANIFEST_OUT <- file.path(OUT, "config/Fig6_PCSF_drug_dual_support_top3_30_75_manifest.json")
LOG_OUT <- file.path(OUT, "logs/Fig6_PCSF_drug_dual_support_top3_30_75.log")
SCRIPT_COPY <- file.path(OUT, "scripts/render_fig6_pcsf_drug_dual_support_top3_30_75.R")

for (p in c(file.path(OUT, "tables"), file.path(OUT, "source_data"),
            file.path(OUT, "figures"), file.path(OUT, "config"),
            file.path(OUT, "scripts"), file.path(OUT, "logs")))
  dir.create(p, recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) as.character(digest::digest(file = path, algo = "sha256", serialize = FALSE))
read_csv <- function(path) data.table::fread(path, na.strings = c("", "NA"), encoding = "UTF-8")
write_csv <- function(x, path) data.table::fwrite(x, path, na = "NA", bom = TRUE)
logmsg <- function(...) cat(paste0("[", format(Sys.time(), "%F %T %z"), "] ",
                                  paste(..., collapse = ""), "\n"), file = LOG_OUT, append = TRUE)
stop_if_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Required input missing: ", paste(missing, collapse = "; "))
}
archive_once <- function(path) {
  if (!file.exists(path)) return(invisible(NULL))
  archive_dir <- file.path(OUT, "archive", "pre_dual_support_top3_30_75")
  dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(archive_dir, basename(path))
  if (file.exists(dest)) dest <- paste0(dest, ".", substr(sha256(path), 1L, 10L))
  if (!file.exists(dest) && !file.copy(path, dest, overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE))
    stop("Could not preserve existing output: ", path)
  if (!identical(sha256(path), sha256(dest))) stop("Archive SHA256 mismatch: ", dest)
  invisible(dest)
}

registry <- read_csv(REGISTRY_CSV)
final_status <- read_csv(FINAL_SUMMARY_CSV)
selected_lock <- jsonlite::fromJSON(SELECTED_LOCK_JSON, simplifyVector = TRUE)
run_manifest <- jsonlite::fromJSON(RUN_MANIFEST_JSON, simplifyVector = TRUE)
selected_parameters_path <- as.character(selected_lock$selected_parameters_path)
if (!file.exists(selected_parameters_path)) stop("锁定参数表缺失")
selected_parameters <- read_csv(selected_parameters_path)
if (!identical(as.character(run_manifest$status), "complete") ||
    !identical(as.character(run_manifest$selected_parameters_sha256), as.character(selected_lock$selected_parameters_sha256)) ||
    !identical(sha256(selected_parameters_path), as.character(selected_lock$selected_parameters_sha256)) ||
    !identical(sha256(RUN_CONFIG_JSON), as.character(selected_lock$run_config_sha256)) ||
    !identical(sha256(RUN_CONFIG_JSON), as.character(run_manifest$run_config_sha256)) ||
    !identical(sha256(INPUT_MANIFEST_CSV), as.character(selected_lock$input_manifest_sha256)) ||
    !identical(sha256(INPUT_MANIFEST_CSV), as.character(run_manifest$input_manifest_sha256)))
  stop("30-75 run manifest或参数锁SHA未通过核验")
if (anyDuplicated(registry$Disease) || anyDuplicated(final_status$Disease) || anyDuplicated(selected_parameters$Disease) ||
    !setequal(final_status$Disease, ALL_DISEASES) || !setequal(registry$Disease, selected_parameters$Disease) ||
    !setequal(registry$Disease, final_status[status == "complete_selected_and_evaluated", Disease]))
  stop("30-75 Disease登记表、状态表和参数表集合不一致")
TARGET_DISEASES <- DISEASE_ORDER[DISEASE_ORDER %chin% registry$Disease]
if (!length(TARGET_DISEASES)) stop("30-75 run没有可展示的已选网络")
LEFT_LAYOUT_DISEASES <- c("LGG", "GBM", "MB", "LCBM", "AD", "EOAD", "PD", "DLB")
RIGHT_LAYOUT_DISEASES <- c("DLBD_PD", "ALS", "ALS_FTD", "FTD", "FTD_MND", "PSP", "HD", "MS")
LEFT_DISEASES <- LEFT_LAYOUT_DISEASES[LEFT_LAYOUT_DISEASES %chin% TARGET_DISEASES]
RIGHT_DISEASES <- RIGHT_LAYOUT_DISEASES[RIGHT_LAYOUT_DISEASES %chin% TARGET_DISEASES]
node_paths <- setNames(file.path(NETWORK_ROOT, TARGET_DISEASES, "selected_nodes.csv"), paste0("selected_nodes_", TARGET_DISEASES))
manifest_paths <- setNames(file.path(NETWORK_ROOT, TARGET_DISEASES, "selected_parameter_manifest.json"), paste0("selected_parameter_manifest_", TARGET_DISEASES))
network <- read_csv(NETWORK_CSV)
reversal <- read_csv(REVERSAL_CSV)
palette <- read_csv(PALETTE_CSV)
input_paths <- c(network_proximity = NETWORK_CSV, historical_reversal = REVERSAL_CSV,
                 disease_palette = PALETTE_CSV, dual_renderer_reference = DUAL_RENDERER,
                 compact_layout_reference = LAYOUT_TEMPLATE, selected_registry = REGISTRY_CSV,
                 final_disease_status = FINAL_SUMMARY_CSV, selected_lock = SELECTED_LOCK_JSON,
                 run_manifest = RUN_MANIFEST_JSON, run_config = RUN_CONFIG_JSON,
                 input_manifest = INPUT_MANIFEST_CSV, selected_parameters = selected_parameters_path,
                 node_paths, manifest_paths)
stop_if_missing(input_paths)
key_cols <- c("Disease", "Drug")
required_network <- c("Disease", "Drug", "network_status", "network_graph_node_count",
  "candidate_unique_entrez_count", "candidate_in_max_component_n", "n_drug_targets_in_max_component",
  "observed_mean_min_distance", "Z_dual_random", "p_empirical", "q_global_BH")
required_reversal <- c("Disease", "Drug", "p_drug", "q_reversal", "formal_deg_status",
  "local_signature_status", "reversal_status", "best_instance_id", "analysis_scope", "direction_status")
required_palette <- c("Disease_Abbrev", "Disease_Color_Hex", "Status")
for (check in list(network = list(x = network, cols = required_network),
                   reversal = list(x = reversal, cols = required_reversal),
                   palette = list(x = palette, cols = required_palette))) {
  missing <- setdiff(check$cols, names(check$x))
  if (length(missing)) stop("Input schema mismatch; missing: ", paste(missing, collapse = ", "))
}
if (anyNA(network$Disease) || anyNA(network$Drug) || anyNA(reversal$Disease) || anyNA(reversal$Drug))
  stop("Disease/Drug join keys contain missing values")
if (anyDuplicated(network, by = key_cols) || anyDuplicated(reversal, by = key_cols))
  stop("Disease + Drug keys must be unique in both frozen input tables")
if (!setequal(unique(network$Disease), ALL_DISEASES) || !setequal(unique(reversal$Disease), ALL_DISEASES) ||
    !setequal(unique(network$Drug), unique(reversal$Drug)) || nrow(network) != 12224L || nrow(reversal) != 12224L)
  stop("Frozen Disease x Drug universe differs from the expected complete key universe")
network_keys <- unique(network[, ..key_cols])
reversal_keys <- unique(reversal[, ..key_cols])
key_audit <- merge(network_keys[, network_key_present := TRUE],
                   reversal_keys[, reversal_key_present := TRUE],
                   by = key_cols, all = TRUE, sort = TRUE)
key_audit[is.na(network_key_present), network_key_present := FALSE]
key_audit[is.na(reversal_key_present), reversal_key_present := FALSE]
network_audit <- rbindlist(lapply(ALL_DISEASES, function(d) {
  summary_row <- final_status[Disease == d]
  reg <- registry[Disease == d]
  stats <- unique(network[Disease == d, .(network_graph_node_count, candidate_unique_entrez_count)])
  disease_network_status <- unique(as.character(network[Disease == d, network_status]))
  if (nrow(summary_row) != 1L || nrow(stats) != 1L || !length(disease_network_status))
    stop("Disease summary or proximity result is missing/duplicated: ", d)
  if (summary_row$status == "no_eligible_candidate") {
    if (nrow(reg) != 0L || stats$network_graph_node_count != 0L || stats$candidate_unique_entrez_count != 0L ||
        !identical(disease_network_status, "empty_candidate_set"))
      stop("No-candidate disease does not have an empty-network result row: ", d)
    return(data.table(Disease = d, selected_node_count = 0L, manifest_status = "no_eligible_candidate",
      network_present = FALSE, exclusion_reason = "no_eligible_candidate_in_30_75_selection"))
  }
  if (summary_row$status != "complete_selected_and_evaluated" || nrow(reg) != 1L) stop("Selected disease status/registry mismatch: ", d)
  locked <- selected_parameters[Disease == d]
  if (nrow(locked) != 1L || as.character(locked$parameter_id) != as.character(reg$parameter_id) ||
      as.numeric(locked$b) != as.numeric(reg$b) || as.integer(locked$w) != as.integer(reg$w) || as.integer(locked$g) != as.integer(reg$g))
    stop("Per-disease parameter id/b/w/g differs across locked table and registry: ", d)
  nodes <- read_csv(node_paths[[paste0("selected_nodes_", d)]])
  manifest <- jsonlite::fromJSON(manifest_paths[[paste0("selected_parameter_manifest_", d)]], simplifyVector = TRUE)
  if (!all(c("Disease", "parameter_id", "Entrez") %in% names(nodes)) ||
      any(nodes$Disease != d) || anyNA(nodes$Entrez) || anyDuplicated(nodes$Entrez))
    stop("selected_nodes identity/Entrez uniqueness failed: ", d)
  norm <- function(x) tolower(gsub("\\\\", "/", as.character(x)))
  if (!identical(norm(reg$selected_nodes_path), norm(node_paths[[paste0("selected_nodes_", d)]])) ||
      !identical(norm(reg$manifest_path), norm(manifest_paths[[paste0("selected_parameter_manifest_", d)]])) ||
      !identical(sha256(node_paths[[paste0("selected_nodes_", d)]]), as.character(reg$selected_nodes_sha256)) ||
      !identical(sha256(manifest_paths[[paste0("selected_parameter_manifest_", d)]]), as.character(reg$manifest_sha256)) ||
      !identical(as.character(manifest$Disease), d) || !identical(as.character(manifest$status), "awaiting_user_review") ||
      !identical(as.character(manifest$selected_parameters_sha256), as.character(selected_lock$selected_parameters_sha256)) ||
      as.character(manifest$parameter_id) != as.character(locked$parameter_id) ||
      as.numeric(manifest$b) != as.numeric(locked$b) || as.integer(manifest$w) != as.integer(locked$w) || as.integer(manifest$g) != as.integer(locked$g) ||
      nrow(nodes) != as.integer(reg$node_count) || stats$network_graph_node_count != nrow(nodes) ||
      stats$candidate_unique_entrez_count != nrow(nodes)) stop("Node table, manifest, registry, and proximity input mismatch: ", d)
  data.table(Disease = d, selected_node_count = nrow(nodes), manifest_status = as.character(manifest$status),
    network_present = TRUE, exclusion_reason = NA_character_)
}), use.names = TRUE)
if (!setequal(network_audit[network_present == TRUE, Disease], TARGET_DISEASES) ||
    !setequal(network_audit[network_present == FALSE, Disease], ALL_DISEASES[!ALL_DISEASES %chin% TARGET_DISEASES]))
  stop("Dynamic active/excluded disease set does not match the final 30-75 registry")
excluded <- network_audit[network_present == FALSE]
# Exact Disease + Drug join; historical p/q and assessment labels are copied, never recomputed.
reverse_fields <- c("p_drug", "q_reversal", "formal_deg_status", "local_signature_status",
                    "reversal_status", "best_instance_id", "analysis_scope", "direction_status")
net <- copy(network)
net[, source_network_row := .I]
rev <- reversal[, c(key_cols, reverse_fields), with = FALSE]
master <- merge(net, rev, by = key_cols, all.x = TRUE, sort = FALSE)
setorder(master, source_network_row)
if (nrow(master) != nrow(network) || anyDuplicated(master, by = key_cols))
  stop("Exact reversal join changed network key cardinality")
master[, `:=`(Z = suppressWarnings(as.numeric(Z_dual_random)),
              d_obs = suppressWarnings(as.numeric(observed_mean_min_distance)),
              q_network_global = suppressWarnings(as.numeric(q_global_BH)),
              p_drug = suppressWarnings(as.numeric(p_drug)),
              q_reversal = suppressWarnings(as.numeric(q_reversal)))]
master[, reversal_record_found := !is.na(formal_deg_status) | !is.na(local_signature_status) |
  !is.na(reversal_status) | is.finite(q_reversal) | is.finite(p_drug)]
master[, reversal_assessment_status := fcase(
  !reversal_record_found, "not_assessed_missing_reversal_key",
  !is.na(formal_deg_status) & formal_deg_status == "no_formal_DEG", "not_assessed_no_formal_DEG",
  !is.na(local_signature_status) & local_signature_status == "no_local_signature", "not_assessed_no_local_signature",
  !is.na(reversal_status) & reversal_status == "evaluated" & is.finite(p_drug) & is.finite(q_reversal), "evaluated",
  default = "not_assessed_missing_or_invalid_reversal_result")]
master[, reversal_cell_status := fcase(
  reversal_assessment_status == "evaluated" & q_reversal < 0.05, "reversal_support",
  reversal_assessment_status == "evaluated", "assessed_no_support",
  default = "not_assessed")]
master[, network_cell_status := fcase(
  is.finite(Z) & Z < -1.96, "exploratory_close",
  is.finite(Z) & Z < -1.645, "marginal_close",
  is.finite(Z), "other_finite_Z",
  suppressWarnings(as.numeric(n_drug_targets_in_max_component)) == 0, "no_mapped_targets",
  default = "undefined_z")]

# Dual-support candidate selection: formal PCSF set present, valid network Z, and assessed historical reversal q<0.05.
valid_network_status <- "complete"
master[, network_result_valid := network_status %chin% valid_network_status]
candidates <- master[Disease %chin% TARGET_DISEASES & network_result_valid &
  network_graph_node_count > 0 & is.finite(Z) & Z < -1.96 &
  reversal_assessment_status == "evaluated" & is.finite(q_reversal) & q_reversal < 0.05]
setorderv(candidates, c("Disease", "q_reversal", "p_drug", "Z", "d_obs", "Drug"),
          c(1L, 1L, 1L, 1L, 1L, 1L), na.last = TRUE)
candidates[, qualifying_rank := seq_len(.N), by = Disease]
candidates[, displayed := qualifying_rank <= MAX_DISPLAY]
display <- candidates[displayed == TRUE]
display[, selection_rank := qualifying_rank]
if (any(display[, .N, by = Disease]$N > MAX_DISPLAY) || any(display$Z >= -1.96) ||
    any(display$q_reversal >= 0.05) || any(display$reversal_assessment_status != "evaluated"))
  stop("Dual-support display violates the frozen eligibility or Top3 rule")

summary_rows <- lapply(TARGET_DISEASES, function(d) {
  x <- master[Disease == d]
  c <- candidates[Disease == d]
  sh <- display[Disease == d]
  n_rev <- sum(x$reversal_assessment_status == "evaluated")
  data.table(Disease = d, Disease_order = match(d, TARGET_DISEASES),
    n_network_keys = nrow(x), n_valid_network_results = sum(x$network_result_valid),
    n_finite_Z = sum(is.finite(x$Z)), n_reversal_assessed = n_rev,
    n_dual_support_qualified = nrow(c), n_displayed = nrow(sh),
    group_status = if (nrow(x) > 0L && sum(x$network_result_valid & is.finite(x$Z)) == 0L)
      "Network proximity not evaluable; dual-support selection unavailable" else if (nrow(x) > 0L && n_rev == 0L)
      "Reversal not assessed; dual-support selection unavailable" else if (nrow(x) > 0L && nrow(c) == 0L)
      "No qualifying dual-support candidate" else NA_character_,
    display_label = sprintf("%s    %d/%d dual-support candidates shown", d, nrow(sh), nrow(c)))
})
display_summary <- rbindlist(summary_rows)

# Frozen palette mapping; GBM uses the authority's GBM_IDHwt color alias.
palette_keys <- unique(c(TARGET_DISEASES, "GBM_IDHwt"))
palette_lookup <- palette[Status %chin% c("frozen_2026-08-29", "frozen_2026-08-04") &
                            Disease_Abbrev %chin% palette_keys,
                          .(Disease_Abbrev, Disease_Color_Hex)]
if (anyDuplicated(palette_lookup$Disease_Abbrev)) stop("Frozen disease palette abbreviations are duplicated")
colors <- setNames(palette_lookup$Disease_Color_Hex, palette_lookup$Disease_Abbrev)
palette_alias <- c(GBM = "GBM_IDHwt")
for (d in TARGET_DISEASES) {
  key <- if (d %in% names(palette_alias)) unname(palette_alias[[d]]) else d
  if (length(colors[key]) != 1L || is.na(colors[key])) stop("Frozen palette is missing disease: ", d)
  colors[d] <- colors[key]
}
colors <- colors[TARGET_DISEASES]

blend <- function(a, b, w) {
  aa <- grDevices::col2rgb(a); bb <- grDevices::col2rgb(b)
  grDevices::rgb((aa[1, ] * (1 - w) + bb[1, ] * w) / 255,
                 (aa[2, ] * (1 - w) + bb[2, ] * w) / 255,
                 (aa[3, ] * (1 - w) + bb[3, ] * w) / 255)
}
master <- master[Disease %chin% TARGET_DISEASES]
master[, target_color_hex := unname(colors[Disease])]
master[, `:=`(network_cell_fill_hex = "#FFFFFF", network_cell_border_hex = target_color_hex,
              network_cell_glyph = "none")]
master[network_cell_status == "exploratory_close", network_cell_fill_hex := blend(target_color_hex, "#202020", 0.30)]
master[network_cell_status == "marginal_close", network_cell_fill_hex := blend(target_color_hex, "#FFFFFF", 0.46)]
master[network_cell_status == "no_mapped_targets", `:=`(network_cell_fill_hex = "#F2F2F2",
  network_cell_border_hex = "#858585", network_cell_glyph = "diagonal_hatch")]
master[network_cell_status == "undefined_z", `:=`(network_cell_fill_hex = "#B8B8B8",
  network_cell_border_hex = "#858585")]
master[, `:=`(reversal_glyph = fifelse(reversal_cell_status == "reversal_support", "filled_up_triangle",
  fifelse(reversal_cell_status == "assessed_no_support", "open_circle", "horizontal_dash")),
  reversal_glyph_hex = fifelse(reversal_cell_status == "reversal_support", target_color_hex, "#68717A"))]

# Expand each displayed source-drug row over every disease with a selected 30-75 network.
display_keys <- display[, .(source_disease = Disease, Drug, selection_rank, qualifying_rank,
  source_Z = Z, source_d_obs = d_obs, source_network_status = network_status,
  source_reversal_status = reversal_assessment_status, source_q_reversal = q_reversal,
  source_p_drug = p_drug)]
target_cells <- master[, .(target_disease = Disease, Drug, network_status, Z, d_obs, p_empirical,
  q_network_global, p_drug, q_reversal, formal_deg_status, local_signature_status, reversal_status,
  best_instance_id, analysis_scope, direction_status, reversal_assessment_status, network_cell_status,
  reversal_cell_status, target_color_hex, network_cell_fill_hex, network_cell_border_hex,
  network_cell_glyph, reversal_glyph, reversal_glyph_hex)]
source <- merge(display_keys, target_cells, by = "Drug", all.x = TRUE, allow.cartesian = TRUE, sort = FALSE)
source[, target_disease_order := match(target_disease, TARGET_DISEASES)]
source[, source_disease_order := match(source_disease, TARGET_DISEASES)]
setorder(source, source_disease_order, selection_rank, target_disease_order)
if (nrow(source) != nrow(display) * length(TARGET_DISEASES) ||
    anyDuplicated(source, by = c("source_disease", "Drug", "target_disease")) ||
    anyNA(source$target_disease_order) ||
    any(source$target_disease_order != match(source$target_disease, TARGET_DISEASES)))
  stop("Figure source does not equal displayed drugs x current selected-disease target universe")
source[, source_disease_order := NULL]

# Verify every source target cell still matches the frozen exact-key joined master.
source_check <- merge(source, master[, .(target_disease = Disease, Drug, check_Z = Z,
  check_d_obs = d_obs, check_p_empirical = p_empirical, check_q_network_global = q_network_global,
  check_p_drug = p_drug, check_q_reversal = q_reversal, check_network_status = network_status,
  check_reversal_status = reversal_status)], by = c("target_disease", "Drug"), all.x = TRUE, sort = FALSE)
if (nrow(source_check) != nrow(source)) stop("Source table cannot be keyed back to the frozen master")
for (field in c("Z", "d_obs", "p_empirical", "q_network_global", "p_drug", "q_reversal",
                "network_status", "reversal_status")) {
  lhs <- source_check[[field]]; rhs <- source_check[[paste0("check_", field)]]
  bad <- xor(is.na(lhs), is.na(rhs)); both <- !is.na(lhs) & !is.na(rhs)
  bad[both] <- as.character(lhs[both]) != as.character(rhs[both])
  if (any(bad)) stop("Drawing source differs from frozen result in field: ", field)
}

input_manifest <- data.table(input = names(input_paths), path = normalizePath(input_paths, winslash = "/", mustWork = TRUE),
  sha256 = vapply(input_paths, sha256, character(1)), size_bytes = as.numeric(file.info(input_paths)$size))

# Compact two-column vector layout adapted from the accepted 12-disease renderer.
W <- 300; GROUPS_PER_COLUMN <- max(length(LEFT_DISEASES), length(RIGHT_DISEASES)); H <- max(225, 123 + 22 * GROUPS_PER_COLUMN)
cell <- 2.15; gap <- 0.42; step <- cell + gap; strip_w <- length(TARGET_DISEASES) * cell + (length(TARGET_DISEASES) - 1L) * gap
panel_x <- c(6, 151)
network_offset <- 52; reversal_offset <- network_offset + strip_w + 3.6
row_top <- H - 40; row_floor <- 48
row_items <- lapply(list(LEFT_DISEASES, RIGHT_DISEASES), function(side_diseases) {
  items <- list()
  for (d in side_diseases) {
    ss <- display_summary[Disease == d]
    items[[length(items) + 1L]] <- list(type = "header", disease = d, label = ss$display_label)
    selected <- display[Disease == d]
    if (nrow(selected)) {
      setorder(selected, selection_rank)
      for (drug in selected$Drug) items[[length(items) + 1L]] <- list(type = "drug", disease = d, Drug = drug, label = drug)
    } else if (!is.na(ss$group_status) && nzchar(ss$group_status)) {
      status_label <- switch(as.character(ss$group_status),
        "Network proximity not evaluable; dual-support selection unavailable" = "Network not evaluable;\ndual-support unavailable",
        "Reversal not assessed; dual-support selection unavailable" = "Reversal not assessed;\ndual-support unavailable",
        "No qualifying dual-support candidate" = "No qualifying dual-support\ncandidate",
        as.character(ss$group_status))
      items[[length(items) + 1L]] <- list(type = "status", disease = d, label = status_label)
    }
    items[[length(items) + 1L]] <- list(type = "gap", disease = d)
  }
  items
})
row_weight <- function(items) sum(vapply(items, function(x) if (x$type == "gap") 0.30 else if (x$type == "status") 1.35 else 1, numeric(1)))
row_step <- (row_top - row_floor) / max(vapply(row_items, row_weight, numeric(1)))
if (row_step < 3.7) stop("Calculated row pitch is too small for the required readable labels")

draw_figure <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(xscale = c(0, W), yscale = c(0, H), clip = "off"))
  grid::grid.rect(x = unit(W / 2, "native"), y = unit(H / 2, "native"),
    width = unit(W, "mm"), height = unit(H, "mm"), gp = gpar(fill = "white", col = NA))
  txt <- function(label, x, y, size = 6, just = "left", face = "plain", col = "black")
    grid.text(label, x = unit(x, "native"), y = unit(y, "native"), just = just,
      gp = gpar(fontfamily = "Arial", fontsize = size, fontface = face, col = col))
  rect <- function(x, y, w, h, fill, border = "#777777", lwd = 0.4)
    grid.rect(x = unit(x, "native"), y = unit(y, "native"), width = unit(w, "mm"), height = unit(h, "mm"),
      gp = gpar(fill = fill, col = border, lwd = lwd))
  line <- function(x1, y1, x2, y2, col = "#68717A", lwd = 0.55)
    grid.lines(x = unit(c(x1, x2), "native"), y = unit(c(y1, y2), "native"), gp = gpar(col = col, lwd = lwd))
  txt("Drug candidates supported by network proximity and expression reversal", 6, H - 9, 10.5, "left", "bold")
  for (side in 1:2) {
    px <- panel_x[side]
    txt("Source disease / drug", px, H - 19, 6.8, "left", "bold")
    txt("Network proximity", px + network_offset + strip_w / 2, H - 19, 6.8, "centre", "bold")
    txt("Reversal support", px + reversal_offset + strip_w / 2, H - 19, 6.8, "centre", "bold")
    for (j in seq_along(TARGET_DISEASES)) {
      nx <- px + network_offset + (j - 1L) * step + cell / 2
      rx <- px + reversal_offset + (j - 1L) * step + cell / 2
      txt(sprintf("%02d", j), nx, H - 25, 5.8, "centre", col = "#303030")
      txt(sprintf("%02d", j), rx, H - 25, 5.8, "centre", col = "#303030")
    }
    y <- row_top
    for (item in row_items[[side]]) {
      if (item$type == "gap") { y <- y - 0.30 * row_step; next }
      yy <- y - row_step / 2
      if (item$type == "header") {
        rect(px + 0.65, yy, 1.35, row_step * 0.70, colors[[item$disease]], colors[[item$disease]], 0.2)
        txt(item$label, px + 2.8, yy, 6.1, "left", "bold")
      } else if (item$type == "status") {
        txt(item$label, px + 2.8, yy, 5.7, "left", col = "#555555")
      } else {
        txt(item$label, px + 2.8, yy, 8.0, "left")
        cells <- source[source_disease == item$disease & Drug == item$Drug]
        setorder(cells, target_disease_order)
        if (nrow(cells) != length(TARGET_DISEASES) ||
            !identical(as.integer(cells$target_disease_order), seq_along(TARGET_DISEASES)))
          stop("Figure cells do not preserve the 01-to-current-target-count disease order")
        for (j in seq_along(TARGET_DISEASES)) {
          cstat <- cells[j]
          nx <- px + network_offset + (j - 1L) * step + cell / 2
          rx <- px + reversal_offset + (j - 1L) * step + cell / 2
          rect(nx, yy, cell, cell, cstat$network_cell_fill_hex, cstat$network_cell_border_hex, 0.34)
          if (cstat$network_cell_glyph == "diagonal_hatch")
            line(nx - cell * 0.34, yy - cell * 0.34, nx + cell * 0.34, yy + cell * 0.34, "#777777", 0.45)
          if (cstat$reversal_glyph == "filled_up_triangle") {
            grid.polygon(x = unit(c(rx - 0.72, rx + 0.72, rx), "native"),
              y = unit(c(yy - 0.52, yy - 0.52, yy + 0.76), "native"),
              gp = gpar(fill = cstat$reversal_glyph_hex, col = cstat$reversal_glyph_hex, lwd = 0.35))
          } else if (cstat$reversal_glyph == "open_circle") {
            grid.circle(x = unit(rx, "native"), y = unit(yy, "native"), r = unit(0.67, "mm"),
              gp = gpar(fill = "white", col = "#68717A", lwd = 0.52))
          } else line(rx - 0.7, yy, rx + 0.7, yy, "#68717A", 0.65)
        }
      }
      y <- y - row_step
    }
  }

  # Three compact legend rows; all network and reversal evidence is displayed independently.
  legend_disease_color <- unname(colors[[TARGET_DISEASES[1L]]])
  txt("Network proximity", 6, 43, 6.5, "left", "bold")
  net_leg <- list(
    list(x=48, fill=blend(legend_disease_color, "#202020", 0.30), border=legend_disease_color, label="Z < −1.96"),
    list(x=102, fill=blend(legend_disease_color, "#FFFFFF", 0.46), border=legend_disease_color, label="−1.96 ≤ Z < −1.645"),
    list(x=185, fill="#FFFFFF", border=legend_disease_color, label="other finite Z"),
    list(x=244, fill="#F2F2F2", border="#858585", label="no mapped drug targets")
  )
  for (it in net_leg) {
    rect(it$x, 43, 2.6, 2.6, it$fill, it$border, 0.42)
    if (grepl("no mapped", it$label)) line(it$x - .75, 42.25, it$x + .75, 43.75, "#777777", .45)
    txt(it$label, it$x + 2.3, 43, 5.8)
  }
  rect(244, 36.5, 2.6, 2.6, "#B8B8B8", "#858585", 0.42)
  txt("undefined Z / other missing", 246.3, 36.5, 5.8)
  txt("Reversal support", 6, 36.5, 6.5, "left", "bold")
  grid.polygon(x = unit(c(57 - .8, 57 + .8, 57), "native"), y = unit(c(36.0, 36.0, 37.4), "native"),
    gp = gpar(fill = legend_disease_color, col = legend_disease_color, lwd = .35))
  txt("q_reversal < 0.05", 59, 36.5, 5.8)
  grid.circle(x = unit(116, "native"), y = unit(36.5, "native"), r = unit(.67, "mm"),
    gp = gpar(fill = "white", col = "#68717A", lwd = .52))
  txt("assessed, no support", 118.3, 36.5, 5.8)
  line(212, 36.5, 214, 36.5, "#68717A", .65); txt("not assessed", 216.3, 36.5, 5.8)

  key_x0 <- 6; KEYS_PER_ROW <- ceiling(length(TARGET_DISEASES) / 2); key_step <- (W - 12) / KEYS_PER_ROW
  for (j in seq_along(TARGET_DISEASES)) {
    row <- if (j <= KEYS_PER_ROW) 1L else 2L
    col <- if (j <= KEYS_PER_ROW) j else j - KEYS_PER_ROW
    yy <- if (row == 1L) 27.5 else 22.0
    xx <- key_x0 + (col - 1L) * key_step
    rect(xx + 1.1, yy, 2.2, 2.2, colors[[TARGET_DISEASES[j]]], colors[[TARGET_DISEASES[j]]], .2)
    txt(sprintf("%02d %s", j, TARGET_DISEASES[j]), xx + 3.4, yy, 5.8)
  }
  txt("Candidates require network proximity (Z < −1.96) and historical expression-reversal support (q < 0.05).", W / 2, 13.7, 6.1, "centre")
  excluded_label <- if (nrow(excluded)) paste(excluded$Disease, collapse = ", ") else "none"
  txt(sprintf("Excluded without an eligible selected network: %s; computational prioritization is not efficacy evidence.", excluded_label),
      W / 2, 8.7, 5.5, "centre")
  txt("Reversal reuses historical frozen expression evidence; it was not recomputed on the selected PCSF nodes.",
      W / 2, 4.2, 5.2, "centre")
  grid::popViewport()
}

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot identify renderer script path")
SCRIPT_PATH <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
for (p in c(FIG_PNG, FIG_PDF, CANDIDATE_OUT, DISPLAY_OUT, SUMMARY_OUT, EXCLUDED_OUT, SOURCE_OUT, KEY_AUDIT_OUT,
            HASH_OUT, CONFIG_OUT, MANIFEST_OUT, SCRIPT_COPY)) archive_once(p)

write_csv(candidates, CANDIDATE_OUT)
write_csv(display, DISPLAY_OUT)
write_csv(display_summary, SUMMARY_OUT)
write_csv(excluded, EXCLUDED_OUT)
write_csv(source, SOURCE_OUT)
write_csv(key_audit, KEY_AUDIT_OUT)
write_csv(input_manifest, HASH_OUT)
if (!file.copy(SCRIPT_PATH, SCRIPT_COPY, overwrite = TRUE) || !identical(sha256(SCRIPT_PATH), sha256(SCRIPT_COPY)))
  stop("Could not save and verify the executable renderer copy")

config <- list(
  figure = "Fig6_PCSF_drug_dual_support_top3_30_75",
  selection_rule = "formal PCSF candidate set nonempty AND valid network result AND finite Z < -1.96 AND historical reversal assessed with finite q_reversal < 0.05; rank by q_reversal, p_drug, Z, d_obs, Drug ascending; top 3 per source disease",
  reversal_use = "exact Disease + Drug join from historical frozen source; used for dual-support eligibility and displayed independently by target disease; no new signatures or BH",
  candidate_network_statuses_accepted = valid_network_status,
  all_diseases = ALL_DISEASES, included_source_and_target_order = TARGET_DISEASES,
  left_source_diseases = LEFT_DISEASES, right_source_diseases = RIGHT_DISEASES,
  excluded_diseases = excluded[, .(Disease, selected_node_count, manifest_status, exclusion_reason)],
  max_display_per_source_disease = MAX_DISPLAY,
  counts = list(network_rows = nrow(network), reversal_rows = nrow(reversal),
    included_network_diseases = length(TARGET_DISEASES), excluded_formal_empty_diseases = nrow(excluded),
    eligible_candidates = nrow(candidates), displayed_rows = nrow(display),
    source_grid_rows = nrow(source), expected_grid_rows = nrow(display) * length(TARGET_DISEASES),
    disease_drug_key_audit = list(network_keys = nrow(network_keys), reversal_keys = nrow(reversal_keys),
      network_without_reversal_key = sum(key_audit$network_key_present & !key_audit$reversal_key_present),
      reversal_without_network_key = sum(!key_audit$network_key_present & key_audit$reversal_key_present)),
    actual_diseases_with_display = uniqueN(display$Disease)),
  target_reversal_unassessed_diseases = display_summary[n_reversal_assessed == 0L, Disease],
  provenance = "Network proximity comes from the frozen current Disease-vs-Ctrl PCSF full-node analysis. Reversal is copied from the historical frozen disease-expression result and was not recomputed on current PCSF nodes; no signatures or test families were changed.",
  layout = list(width_mm = W, height_mm = H, dpi = 600, columns = 2, source_groups_per_column = c(length(LEFT_DISEASES), length(RIGHT_DISEASES)),
    target_disease_order = TARGET_DISEASES, editable_pdf_family = "Arial", row_step_mm = row_step),
  inputs = input_manifest,
  outputs = c(FIG_PNG, FIG_PDF, CANDIDATE_OUT, DISPLAY_OUT, SUMMARY_OUT, EXCLUDED_OUT, SOURCE_OUT, KEY_AUDIT_OUT,
    HASH_OUT, SCRIPT_COPY))
jsonlite::write_json(config, CONFIG_OUT, pretty = TRUE, auto_unbox = TRUE, na = "null")

grDevices::png(FIG_PNG, width = W / 25.4, height = H / 25.4, units = "in", res = 600, type = "cairo", bg = "white")
draw_figure(); grDevices::dev.off()
grDevices::cairo_pdf(FIG_PDF, width = W / 25.4, height = H / 25.4, family = "Arial", onefile = TRUE, bg = "white")
draw_figure(); grDevices::dev.off()

if (!file.exists(FIG_PNG) || !file.exists(FIG_PDF) || file.info(FIG_PNG)$size <= 0 || file.info(FIG_PDF)$size <= 0)
  stop("Figure export missing or empty")
manifest <- list(status = "complete", finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  exit_code = 0L, figure_config_path = CONFIG_OUT, script_path = SCRIPT_COPY,
  script_sha256 = sha256(SCRIPT_COPY), input_hash_table = HASH_OUT,
  inputs_unchanged_check = input_manifest,
  counts = list(included_diseases = length(TARGET_DISEASES), displayed_drug_rows = nrow(display),
    eligible_candidates = nrow(candidates), source_cells = nrow(source), source_cells_expected = nrow(display) * length(TARGET_DISEASES)),
  outputs = data.table(path = c(FIG_PNG, FIG_PDF, CANDIDATE_OUT, DISPLAY_OUT, SUMMARY_OUT, EXCLUDED_OUT, SOURCE_OUT, KEY_AUDIT_OUT,
    HASH_OUT, CONFIG_OUT, SCRIPT_COPY), sha256 = vapply(c(FIG_PNG, FIG_PDF, CANDIDATE_OUT, DISPLAY_OUT, SUMMARY_OUT,
    EXCLUDED_OUT, SOURCE_OUT, KEY_AUDIT_OUT, HASH_OUT, CONFIG_OUT, SCRIPT_COPY), sha256, character(1)),
    size_bytes = as.numeric(file.info(c(FIG_PNG, FIG_PDF, CANDIDATE_OUT, DISPLAY_OUT, SUMMARY_OUT, EXCLUDED_OUT,
      SOURCE_OUT, KEY_AUDIT_OUT, HASH_OUT, CONFIG_OUT, SCRIPT_COPY))$size)))
jsonlite::write_json(manifest, MANIFEST_OUT, pretty = TRUE, auto_unbox = TRUE, na = "null")
logmsg("COMPLETE included_diseases=", length(TARGET_DISEASES), " eligible_dual=", nrow(candidates),
       " displayed=", nrow(display), " source_cells=", nrow(source), " row_step_mm=", round(row_step, 3))
cat("DUAL_SUPPORT_TOP3_30_75_COMPLETE\nDISEASES=", length(TARGET_DISEASES), "\nELIGIBLE=", nrow(candidates),
    "\nDISPLAYED=", nrow(display), "\nSOURCE_CELLS=", nrow(source), "\nPNG=", FIG_PNG,
    "\nPDF=", FIG_PDF, "\n", sep = "")

