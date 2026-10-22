#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(Rcpp)
  library(digest)
  library(jsonlite)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

if (!nzchar(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.65001"))) {
  stop("无法设置UTF-8字符区域", call. = FALSE)
}
options(stringsAsFactors = FALSE, scipen = 999)

PROJECT_ROOT <- "E:/micro_output/分析审阅中心/待审阅/Fig6_疾病vsCtrl_SHAP_PCSF模块"
RUN_ROOT <- file.path(PROJECT_ROOT, "robustness/20261021_disease_compact_30_75_parameter_selection_01")
OUT <- file.path(RUN_ROOT, "drug_accessibility")
NETWORK_ROOT <- file.path(RUN_ROOT, "selected_networks")
REGISTRY_FILE <- file.path(RUN_ROOT, "summaries/selected_network_delivery_registry.csv")
FINAL_SUMMARY_FILE <- file.path(RUN_ROOT, "summaries/final_disease_summary.csv")
SELECTED_PARAMETERS_FILE <- file.path(RUN_ROOT, "summaries/selected_parameters.csv")
SELECTED_LOCK_FILE <- file.path(RUN_ROOT, "config/selected_parameters.lock.json")
RUN_MANIFEST_FILE <- file.path(RUN_ROOT, "summaries/run_manifest.json")
PPI_FILE <- "D:/data/microglia_new/input/drug_repurposing_38core/mentor_20260923/Network-1.txt"
DRUG_FILE <- "E:/micro_output/fig3/input/DrugTarget-final (1).txt"
DISEASES <- c("AD", "ALS", "ALS_FTD", "DLB", "DLBD_PD", "EOAD", "FTD", "FTD_MND",
              "GBM", "HD", "LCBM", "LGG", "MB", "MS", "PD", "PSP")
B <- 1000L
MAIN_SEED <- 42L
MIN_BIN_SIZE <- 100L
PPI_MIN_SUPPORT <- 2L
IMPL_VERSION <- "Fig6_PCSF_30_75_selected_nodes_dual_random_ge2_v1"

for (d in c("data", "tables", "checkpoints", "scripts", "logs"))
  dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(OUT, "logs", "drug_accessibility_run.log")
writeLines(as.character(Sys.getpid()), file.path(OUT, "logs", "drug_accessibility_pid.txt"), useBytes = TRUE)
log_line <- function(...) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "] ", paste(..., collapse = ""))
  cat(line, "\n", file = LOG_FILE, append = TRUE, sep = "")
  flush.console()
}
write_csv <- function(x, path) data.table::fwrite(x, path, na = "NA", bom = TRUE)
sha256 <- function(path) as.character(digest::digest(file = path, algo = "sha256", serialize = FALSE))
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
SCRIPT_PATH <- if (length(script_arg) == 1L) normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE) else ""
if (!nzchar(SCRIPT_PATH)) stop("无法确定执行脚本路径，停止以保留可追溯代码")
SCRIPT_COPY <- file.path(OUT, "scripts", basename(SCRIPT_PATH))
if (!file.copy(SCRIPT_PATH, SCRIPT_COPY, overwrite = TRUE) || !identical(sha256(SCRIPT_PATH), sha256(SCRIPT_COPY)))
  stop("无法保存当前已审查R脚本的副本")
atomic_save_rds <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(x, tmp, compress = "gzip")
  if (!file.rename(tmp, path)) stop("无法原子完成检查点：", path)
}
assert_cols <- function(x, need, label) {
  miss <- setdiff(need, names(x))
  if (length(miss)) stop(label, " 缺少字段：", paste(miss, collapse = ", "))
}
stable_seed <- function(...) substr(digest::digest(paste(c(MAIN_SEED, ...), collapse = "|"),
                                                   algo = "sha256", serialize = FALSE), 1L, 16L)

required <- c(PPI_FILE, DRUG_FILE, REGISTRY_FILE, FINAL_SUMMARY_FILE, SELECTED_PARAMETERS_FILE, SELECTED_LOCK_FILE, RUN_MANIFEST_FILE)
if (!all(file.exists(required))) stop("缺少必要输入：", paste(required[!file.exists(required)], collapse = "; "))
selection_registry <- fread(REGISTRY_FILE, na.strings = c("", "NA"), encoding = "UTF-8")
final_disease_status <- fread(FINAL_SUMMARY_FILE, na.strings = c("", "NA"), encoding = "UTF-8")
selected_parameters <- fread(SELECTED_PARAMETERS_FILE, na.strings = c("", "NA"), encoding = "UTF-8")
selection_lock <- fromJSON(SELECTED_LOCK_FILE, simplifyVector = TRUE)
run_manifest <- fromJSON(RUN_MANIFEST_FILE, simplifyVector = TRUE)
run_config_file <- file.path(RUN_ROOT, "config/run_config.json")
input_manifest_file <- file.path(RUN_ROOT, "inputs/input_sha256_manifest.csv")
if (!file.exists(run_config_file) || !file.exists(input_manifest_file) || !identical(as.character(run_manifest$status), "complete") ||
    !identical(as.character(run_manifest$selected_parameters_sha256), as.character(selection_lock$selected_parameters_sha256)) ||
    !identical(sha256(run_config_file), as.character(selection_lock$run_config_sha256)) ||
    !identical(sha256(input_manifest_file), as.character(selection_lock$input_manifest_sha256)))
  stop("30-75运行总manifest、配置或输入SHA未通过锁定核验")
assert_cols(selection_registry, c("Disease", "status", "b", "w", "g", "parameter_id", "node_count", "edge_count", "selected_nodes_path", "selected_nodes_sha256", "selected_edges_path", "selected_edges_sha256", "manifest_path", "manifest_sha256"), "30-75 selected network registry")
assert_cols(final_disease_status, c("Disease", "status"), "30-75 final disease summary")
assert_cols(selected_parameters, c("Disease", "parameter_id", "b", "w", "g"), "30-75 locked selected parameters")
if (anyDuplicated(selection_registry, by = "Disease") || anyDuplicated(final_disease_status, by = "Disease") ||
    anyDuplicated(selected_parameters, by = "Disease") || !setequal(final_disease_status$Disease, DISEASES) ||
    !setequal(selection_registry$Disease, selected_parameters$Disease) ||
    !setequal(selection_registry$Disease, final_disease_status[status == "complete_selected_and_evaluated", Disease]) ||
    !identical(as.character(selection_lock$selected_parameters_sha256), sha256(SELECTED_PARAMETERS_FILE)))
  stop("30-75 parameter lock, disease summary, and selected network registry do not agree")

ppi_raw <- fread(PPI_FILE, na.strings = c("", "NA"), encoding = "UTF-8")
drug_raw <- fread(DRUG_FILE, na.strings = c("", "NA"), encoding = "UTF-8")
assert_cols(ppi_raw, c("Protein_A_Entrez_ID", "Protein_B_Entrez_ID", "data_source.s.", "number"), "导师PPI")
assert_cols(drug_raw, c("Drug", "Target"), "药物靶点表")
drug_raw[, `:=`(Drug = trimws(as.character(Drug)), Target = trimws(as.character(Target)))]
drug_raw <- unique(drug_raw[!is.na(Drug) & nzchar(Drug) & !is.na(Target) & nzchar(Target)],
                   by = c("Drug", "Target"))
DRUGS <- sort(unique(drug_raw$Drug), method = "radix")

number_numeric <- suppressWarnings(as.numeric(ppi_raw$number))
if (any(!is.finite(number_numeric)) || any(number_numeric != trunc(number_numeric)) || any(number_numeric < 0) || any(number_numeric > .Machine$integer.max))
  stop("PPI number字段含非有限、非整数或超出整数范围的值")
ppi_raw[, `:=`(a = trimws(as.character(Protein_A_Entrez_ID)),
               b = trimws(as.character(Protein_B_Entrez_ID)),
               sources = as.character(`data_source.s.`),
               number_num = as.integer(number_numeric))]
if (anyNA(ppi_raw$number_num) || anyNA(ppi_raw$a) || anyNA(ppi_raw$b) || any(!nzchar(ppi_raw$a)) || any(!nzchar(ppi_raw$b)))
  stop("PPI端点或number字段包含缺失/非法值")
ppi_raw[, source_n_recounted := vapply(strsplit(sources, ",", fixed = TRUE), function(x) {
  x <- trimws(x); length(unique(x[nzchar(x)]))
}, integer(1))]
ppi_source_audit <- ppi_raw[, .(n_records = .N, n_source_number_mismatch = sum(source_n_recounted != number_num),
                               number_min = min(number_num), number_max = max(number_num))]
if (ppi_source_audit$n_source_number_mismatch != 0L) {
  write_csv(ppi_raw[source_n_recounted != number_num,
                    .(Protein_A_Entrez_ID = a, Protein_B_Entrez_ID = b, data_source.s. = sources,
                      number = number_num, source_n_recounted)],
            file.path(OUT, "tables", "ppi_source_number_mismatches.csv"))
  stop("PPI来源去重计数与number存在不一致；保留证据并停止，不静默改筛边口径")
}

ppi_ge2 <- ppi_raw[number_num >= PPI_MIN_SUPPORT]
ppi_ge2[, `:=`(from = pmin(a, b), to = pmax(a, b))]
n_self_loop <- ppi_ge2[from == to, .N]
ppi_nonloop <- ppi_ge2[from != to]
ppi_edges <- ppi_nonloop[, {
  src <- trimws(unlist(strsplit(sources, ",", fixed = TRUE)))
  src <- sort(unique(src[nzchar(src)]), method = "radix")
  .(number = length(src), data_sources = paste(src, collapse = ","), source_record_n = .N)
}, by = .(from, to)]
setorder(ppi_edges, from, to)
ppi_graph <- graph_from_data_frame(ppi_edges[, .(from, to)], directed = FALSE)
ppi_graph <- simplify(ppi_graph, remove.multiple = TRUE, remove.loops = TRUE)
cc <- components(ppi_graph)
largest_n <- max(cc$csize)
largest_ids <- which(cc$csize == largest_n)
largest_id <- largest_ids[which.min(vapply(largest_ids, function(i) {
  min(suppressWarnings(as.numeric(V(ppi_graph)$name[cc$membership == i])))
}, numeric(1)))]
all_ppi_nodes <- data.table(Entrez = as.character(V(ppi_graph)$name),
                            component_id = as.integer(cc$membership),
                            component_size = as.integer(cc$csize[cc$membership]),
                            in_max_component = cc$membership == largest_id,
                            degree = as.integer(degree(ppi_graph)))
all_ppi_nodes[, Entrez_numeric := suppressWarnings(as.numeric(Entrez))]
setorderv(all_ppi_nodes, c("component_id", "Entrez_numeric"))
all_ppi_nodes[, Entrez_numeric := NULL]
main_nodes <- all_ppi_nodes[in_max_component == TRUE]
main_graph <- induced_subgraph(ppi_graph, vids = V(ppi_graph)[name %chin% main_nodes$Entrez])
if (!is_connected(main_graph)) stop("按最大连通分量定义诱导出的距离背景不连通")
main_nodes <- data.table(Entrez = as.character(V(main_graph)$name), degree = as.integer(degree(main_graph)))
if (!identical(main_nodes$Entrez, as.character(V(main_graph)$name))) stop("距离背景节点索引顺序不一致")
node_index <- setNames(seq_len(nrow(main_nodes)), main_nodes$Entrez)
main_ids <- main_nodes$Entrez
all_graph_ids <- all_ppi_nodes$Entrez

candidate_nodes_list <- vector("list", length(DISEASES)); names(candidate_nodes_list) <- DISEASES
forest_audit <- vector("list", length(DISEASES)); names(forest_audit) <- DISEASES
input_paths <- c(PPI = PPI_FILE, DrugTarget = DRUG_FILE, selected_network_registry = REGISTRY_FILE,
                 final_disease_summary = FINAL_SUMMARY_FILE, selected_parameters = SELECTED_PARAMETERS_FILE,
                 selected_parameters_lock = SELECTED_LOCK_FILE, run_manifest = RUN_MANIFEST_FILE, run_config = run_config_file, input_manifest = input_manifest_file)
for (d in DISEASES) {
  status_row <- final_disease_status[Disease == d]
  if (nrow(status_row) != 1L) stop("30-75疾病状态记录缺失或重复：", d)
  if (status_row$status == "no_eligible_candidate") {
    if (d %chin% selection_registry$Disease) stop("无候选疾病错误出现在网络登记中：", d)
    candidate_nodes_list[[d]] <- data.table(Disease = character(), Entrez = character())
    forest_audit[[d]] <- data.table(Disease = d, n_forest_nodes = 0L, n_selected_nodes = 0L,
                                    n_forest_edges = 0L, forest_endpoints_subset_selected_nodes = TRUE,
                                    network_status = "empty_candidate_set")
    next
  }
  if (status_row$status != "complete_selected_and_evaluated") stop("当前疾病状态不是完整参数选择及评价：", d)
  reg <- selection_registry[Disease == d]
  if (nrow(reg) != 1L || reg$status != "awaiting_user_review") stop("30-75网络登记状态不符：", d)
  ddir <- file.path(NETWORK_ROOT, d)
  node_path <- file.path(ddir, "selected_nodes.csv")
  forest_path <- file.path(ddir, "selected_forest_edges.csv")
  manifest_path <- file.path(ddir, "selected_parameter_manifest.json")
  expected_paths <- c(node_path, forest_path, manifest_path)
  if (!all(file.exists(expected_paths))) stop("锁定疾病网络文件缺失：", d)
  if (!identical(normalizePath(reg$selected_nodes_path, winslash = "/", mustWork = TRUE), normalizePath(node_path, winslash = "/", mustWork = TRUE)) ||
      !identical(normalizePath(reg$selected_edges_path, winslash = "/", mustWork = TRUE), normalizePath(forest_path, winslash = "/", mustWork = TRUE)) ||
      !identical(normalizePath(reg$manifest_path, winslash = "/", mustWork = TRUE), normalizePath(manifest_path, winslash = "/", mustWork = TRUE)))
    stop("网络登记路径不指向当前30-75 run：", d)
  if (!identical(sha256(node_path), as.character(reg$selected_nodes_sha256)) ||
      !identical(sha256(forest_path), as.character(reg$selected_edges_sha256)) ||
      !identical(sha256(manifest_path), as.character(reg$manifest_sha256))) stop("当前网络文件SHA与登记不一致：", d)
  nodes <- fread(node_path, na.strings = c("", "NA"), encoding = "UTF-8")
  forest <- fread(forest_path, na.strings = c("", "NA"), encoding = "UTF-8")
  manifest <- fromJSON(manifest_path, simplifyVector = TRUE)
  assert_cols(nodes, c("Disease", "parameter_id", "Entrez", "model_SHAP_evidence", "mean_abs_SHAP", "normalized_prize", "solver_prize"), paste0(d, " selected_nodes"))
  assert_cols(forest, c("Disease", "Entrez_A", "Entrez_B"), paste0(d, " selected_forest_edges"))
  locked_disease_row <- selected_parameters[Disease == d]
  if (nrow(locked_disease_row) != 1L || as.character(locked_disease_row$parameter_id) != as.character(reg$parameter_id) ||
      as.numeric(locked_disease_row$b) != as.numeric(reg$b) || as.integer(locked_disease_row$w) != as.integer(reg$w) || as.integer(locked_disease_row$g) != as.integer(reg$g))
    stop("网络登记与锁定参数表不一致：", d)
  if (!identical(as.character(manifest$Disease), d) || !identical(as.character(manifest$status), "awaiting_user_review") ||
      !identical(as.character(manifest$selected_parameters_sha256), as.character(selection_lock$selected_parameters_sha256)) ||
      any(nodes$Disease != d) || any(forest$Disease != d) || any(nodes$parameter_id != reg$parameter_id) ||
      nrow(nodes) != as.integer(reg$node_count) || nrow(forest) != as.integer(reg$edge_count) ||
      as.numeric(manifest$b) != as.numeric(reg$b) || as.integer(manifest$w) != as.integer(reg$w) || as.integer(manifest$g) != as.integer(reg$g))
    stop("疾病身份、参数锁或网络计数与manifest不一致：", d)
  nodes[, Entrez := as.character(Entrez)]
  if (anyNA(nodes$Entrez) || any(!nzchar(nodes$Entrez)) || anyDuplicated(nodes$Entrez)) stop("PCSF候选Entrez非唯一或无效：", d)
  nodes[, `:=`(network_graph_node_count = .N, in_ge2_ppi = Entrez %chin% all_graph_ids,
               in_max_component = Entrez %chin% main_ids)]
  nodes[, exclusion_reason := fifelse(in_max_component, "entered_max_component",
                             fifelse(in_ge2_ppi, "outside_max_component", "not_in_number_ge2_ppi"))]
  candidate_nodes_list[[d]] <- nodes
  forest_endpoints <- unique(c(as.character(forest$Entrez_A), as.character(forest$Entrez_B)))
  forest_audit[[d]] <- data.table(Disease = d, n_forest_nodes = length(forest_endpoints),
                                  n_selected_nodes = nrow(nodes), n_forest_edges = nrow(forest),
                                  forest_endpoints_subset_selected_nodes = all(forest_endpoints %chin% nodes$Entrez),
                                  network_status = "complete")
  if (!forest_audit[[d]]$forest_endpoints_subset_selected_nodes) stop("PCSF骨架含selected_nodes之外的端点：", d)
  input_paths <- c(input_paths, setNames(expected_paths, paste(d, c("nodes", "forest", "manifest"), sep = "_")))
}
candidate_nodes <- rbindlist(candidate_nodes_list, use.names = TRUE, fill = TRUE)
candidate_count_check <- candidate_nodes[, .(n_rows = .N, n_unique_ids = uniqueN(Entrez),
                                              n_annotated_values = uniqueN(network_graph_node_count),
                                              annotated_count = network_graph_node_count[[1L]]), by = Disease]
if (any(candidate_count_check$n_rows != candidate_count_check$n_unique_ids) ||
    any(candidate_count_check$n_annotated_values != 1L) || any(candidate_count_check$n_rows != candidate_count_check$annotated_count) ||
    anyDuplicated(candidate_nodes, by = c("Disease", "Entrez"))) stop("疾病PCSF候选全集不是疾病内唯一节点表")
write_csv(rbindlist(forest_audit, use.names = TRUE, fill = TRUE), file.path(OUT, "tables", "network_graph_forest_node_audit.csv"))
write_csv(candidate_nodes, file.path(OUT, "tables", "all_disease_candidate_nodes_and_distance_coverage.csv"))

coverage <- candidate_nodes[, .(
  network_graph_node_count = .N,
  candidate_unique_entrez_count = uniqueN(Entrez),
  entered_distance_max_component_n = uniqueN(Entrez[in_max_component == TRUE]),
  not_in_number_ge2_ppi_n = uniqueN(Entrez[exclusion_reason == "not_in_number_ge2_ppi"]),
  outside_max_component_n = uniqueN(Entrez[exclusion_reason == "outside_max_component"]),
  distance_input_count = uniqueN(Entrez[in_max_component == TRUE]),
  status = if (.N == 0L) "empty_candidate_set" else if (all(in_max_component)) "all_candidates_in_distance_background" else "candidate_set_intersected_with_max_component"
), by = Disease]
coverage <- merge(data.table(Disease = DISEASES), coverage, by = "Disease", all.x = TRUE, sort = FALSE)
for (nm in c("network_graph_node_count", "candidate_unique_entrez_count", "entered_distance_max_component_n",
             "not_in_number_ge2_ppi_n", "outside_max_component_n", "distance_input_count"))
  set(coverage, which(is.na(coverage[[nm]])), nm, 0L)
coverage[is.na(status), status := "empty_candidate_set"]
coverage[, disease_order := match(Disease, DISEASES)]
setorderv(coverage, "disease_order")
coverage[, disease_order := NULL]
write_csv(coverage, file.path(OUT, "tables", "disease_candidate_distance_coverage.csv"))

ppi_edges_out <- copy(ppi_edges)
write.table(ppi_ge2[, .(Protein_A_Entrez_ID = a, Protein_B_Entrez_ID = b,
                        data_source.s. = sources, number = number_num, source_n_recounted, from, to)],
            file.path(OUT, "data", "ppi_number_ge2_source_records.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
fwrite(ppi_edges_out, file.path(OUT, "data", "ppi_number_ge2_unique_undirected_edges.tsv"), sep = "\t", na = "NA")
write_csv(all_ppi_nodes, file.path(OUT, "data", "ppi_number_ge2_component_node_manifest.csv"))
write_csv(main_nodes, file.path(OUT, "data", "ppi_number_ge2_max_component_nodes.csv"))
main_edgelist <- as_edgelist(main_graph, names = TRUE)
main_edge_keys <- paste(pmin(main_edgelist[, 1L], main_edgelist[, 2L]),
                        pmax(main_edgelist[, 1L], main_edgelist[, 2L]), sep = "|")
edge_keys_all <- paste(ppi_edges$from, ppi_edges$to, sep = "|")
main_edge_idx <- match(main_edge_keys, edge_keys_all)
if (anyNA(main_edge_idx)) stop("最大连通分量边未能映射回number>=2批准PPI边")
fwrite(ppi_edges[main_edge_idx], file.path(OUT, "data", "ppi_number_ge2_max_component_edges.tsv"), sep = "\t", na = "NA")
edge_component <- cc$membership[as_edgelist(ppi_graph, names = FALSE)[, 1L]]
component_edges <- tabulate(edge_component, nbins = length(cc$csize))
write_csv(data.table(component_id = seq_along(cc$csize), n_nodes = as.integer(cc$csize),
                     n_edges = as.integer(component_edges),
                     in_max_component = seq_along(cc$csize) == largest_id),
          file.path(OUT, "tables", "ppi_component_summary.csv"))
write_csv(data.table(
  step = c("raw_source_records", "source_number_consistency_mismatches", "number_ge2_source_records",
           "number_ge2_self_loop_records_removed", "nonloop_records_before_dedup", "unique_undirected_edges",
           "background_nodes", "background_components", "maximum_component_nodes", "maximum_component_edges"),
  n = c(nrow(ppi_raw), ppi_source_audit$n_source_number_mismatch, nrow(ppi_ge2), n_self_loop,
        nrow(ppi_nonloop), nrow(ppi_edges), vcount(ppi_graph), length(cc$csize),
        vcount(main_graph), ecount(main_graph)),
  max_component_id = c(rep(NA_integer_, 9L), as.integer(largest_id))),
  file.path(OUT, "tables", "ppi_filter_and_component_audit.csv"))

map_symbols <- function(symbols) {
  keys <- sort(unique(trimws(as.character(symbols[!is.na(symbols) & nzchar(trimws(as.character(symbols)))]))), method = "radix")
  if (!length(keys)) return(data.table(symbol = character(), map_status = character(), Entrez = character(), n_entrez = integer()))
  mapped <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = keys, keytype = "SYMBOL", columns = "ENTREZID"))
  mapped <- unique(as.data.table(mapped)[!is.na(ENTREZID), .(symbol = as.character(SYMBOL), Entrez = as.character(ENTREZID))])
  counts <- mapped[, .(n_entrez = uniqueN(Entrez), unique_entrez = if (uniqueN(Entrez) == 1L) unique(Entrez) else NA_character_), by = symbol]
  out <- data.table(symbol = keys)[counts, on = "symbol"]
  out[is.na(n_entrez), n_entrez := 0L]
  out[, map_status := fifelse(n_entrez == 1L, "unique_official_symbol_to_entrez",
                     fifelse(n_entrez > 1L, "ambiguous_one_to_many", "unmapped_symbol"))]
  out[, Entrez := unique_entrez]
  out[, unique_entrez := NULL]
  out
}
target_map <- map_symbols(drug_raw$Target)
drug_raw <- target_map[drug_raw, on = c(symbol = "Target")]
setnames(drug_raw, "symbol", "target_symbol")
drug_raw[, `:=`(Target = target_symbol,
                in_ge2_ppi = !is.na(Entrez) & Entrez %chin% all_graph_ids,
                in_max_component = !is.na(Entrez) & Entrez %chin% main_ids)]
drug_raw[, mapping_network_status := map_status]
drug_raw[map_status == "unique_official_symbol_to_entrez" & !in_ge2_ppi, mapping_network_status := "unique_mapped_not_in_ge2_ppi"]
drug_raw[map_status == "unique_official_symbol_to_entrez" & in_ge2_ppi & !in_max_component, mapping_network_status := "unique_mapped_outside_max_component"]
drug_raw[map_status == "unique_official_symbol_to_entrez" & in_max_component, mapping_network_status := "unique_mapped_in_max_component"]
write_csv(drug_raw, file.path(OUT, "tables", "drug_target_mapping_audit.csv"))
write_csv(drug_raw[, .(n_target_relations = .N, n_unique_target_symbols = uniqueN(Target),
                       n_unique_entrez = uniqueN(Entrez[!is.na(Entrez)]),
                       n_unique_entrez_in_max_component = uniqueN(Entrez[in_max_component]),
                       n_unmapped_symbols = uniqueN(Target[map_status == "unmapped_symbol"]),
                       n_ambiguous_symbols = uniqueN(Target[map_status == "ambiguous_one_to_many"])), by = Drug],
          file.path(OUT, "tables", "drug_target_coverage_by_drug.csv"))

# Degree bins preserve exact-degree blocks and accumulate until >=100 nodes; a final tail joins the prior bin.
ord <- order(main_nodes$degree, suppressWarnings(as.numeric(main_nodes$Entrez)), method = "radix")
degree_blocks <- split(ord, factor(main_nodes$degree[ord], levels = sort(unique(main_nodes$degree))))
bin_members <- list(); current <- integer()
for (block in degree_blocks) {
  current <- c(current, block)
  if (length(current) >= MIN_BIN_SIZE) { bin_members[[length(bin_members) + 1L]] <- current; current <- integer() }
}
if (length(current)) {
  if (length(bin_members)) bin_members[[length(bin_members)]] <- c(bin_members[[length(bin_members)]], current)
  else bin_members[[1L]] <- current
}
bin_id <- integer(nrow(main_nodes))
bins <- rbindlist(lapply(seq_along(bin_members), function(i) {
  bin_id[bin_members[[i]]] <<- i
  data.table(bin_id = i - 1L, degree_min = min(main_nodes$degree[bin_members[[i]]]),
             degree_max = max(main_nodes$degree[bin_members[[i]]]), n_nodes = length(bin_members[[i]]),
             node_ids = paste(main_nodes$Entrez[bin_members[[i]]], collapse = ";"))
}))
write_csv(bins, file.path(OUT, "tables", "ppi_degree_bins.csv"))
bin_pools <- lapply(bin_members, as.integer)
match_counts <- function(idx) tabulate(bin_id[idx], nbins = length(bin_members))

CPP_FILE <- file.path(OUT, "scripts", "fig6_dual_random_bfs_sampler.cpp")
cpp_source <- '#include <Rcpp.h>\n#include <cstdint>\n#include <random>\n#include <string>\n#include <unordered_set>\n#include <vector>\nusing namespace Rcpp;\n// [[Rcpp::export]]\nIntegerMatrix fig6_sample_stratified(List pools, IntegerVector counts, int B, std::string seed_hex) { int nb=pools.size(); if(counts.size()!=nb||B<1) stop("invalid dimensions"); int total=sum(counts); IntegerMatrix out(B,total); std::mt19937_64 rng((uint64_t)std::stoull(seed_hex,nullptr,16)); for(int b=0;b<B;++b){int col=0; for(int j=0;j<nb;++j){IntegerVector pool=pools[j]; int n=pool.size(),k=counts[j]; if(k<0||k>n) stop("degree-bin pool too small"); std::unordered_set<int> used; used.reserve((size_t)k*2+1); for(int t=n-k;t<n;++t){std::uniform_int_distribution<int> draw(0,t); int x=draw(rng); int y=used.find(x)==used.end()?x:t; used.insert(y); out(b,col++)=pool[y];}}} return out;}\n// [[Rcpp::export]]\nIntegerMatrix fig6_batched_multisource_bfs(IntegerVector from, IntegerVector to, int n, IntegerMatrix sources) { if(from.size()!=to.size()) stop("edge vectors differ"); std::vector<std::vector<int>> adj(n); for(int e=0;e<from.size();++e){int a=from[e]-1,b=to[e]-1; if(a<0||b<0||a>=n||b>=n) stop("edge index outside graph"); adj[a].push_back(b); adj[b].push_back(a);} int B=sources.nrow(),k=sources.ncol(); IntegerMatrix out(B,n); std::vector<int> q(n); for(int r=0;r<B;++r){for(int v=0;v<n;++v) out(r,v)=-1; int head=0,tail=0; for(int j=0;j<k;++j){int s=sources(r,j)-1; if(s<0||s>=n) stop("source index outside graph"); if(out(r,s)<0){out(r,s)=0;q[tail++]=s;}} while(head<tail){int u=q[head++],du=out(r,u)+1; for(int v:adj[u]) if(out(r,v)<0){out(r,v)=du;q[tail++]=v;}} if(tail!=n) stop("random source BFS did not reach every maximum-component node");} return out;}\n// [[Rcpp::export]]\nNumericVector fig6_target_null_mean(IntegerMatrix min_dist, List pools, IntegerVector counts, int B, std::string seed_hex) {int nb=pools.size(); if(counts.size()!=nb||min_dist.nrow()!=B||B<1) stop("invalid dimensions"); int total=sum(counts); if(total<1) stop("empty target set"); NumericVector out(B); std::mt19937_64 rng((uint64_t)std::stoull(seed_hex,nullptr,16)); for(int b=0;b<B;++b){double sumd=0.0;for(int j=0;j<nb;++j){IntegerVector pool=pools[j];int n=pool.size(),k=counts[j];if(k<0||k>n) stop("target degree-bin pool too small");std::unordered_set<int> used;used.reserve((size_t)k*2+1);for(int t=n-k;t<n;++t){std::uniform_int_distribution<int> draw(0,t);int x=draw(rng);int y=used.find(x)==used.end()?x:t;used.insert(y);sumd+=min_dist(b,pool[y]-1);}}out[b]=sumd/total;}return out;}\n'
writeLines(cpp_source, CPP_FILE, useBytes = TRUE)
Rcpp::sourceCpp(CPP_FILE, rebuild = TRUE, showOutput = TRUE, verbose = FALSE)

input_paths <- c(input_paths, CPP = CPP_FILE, script_source = SCRIPT_PATH, script_preserved_copy = SCRIPT_COPY)
input_paths <- input_paths[file.exists(input_paths)]
input_manifest <- rbindlist(lapply(input_paths, function(p) data.table(path = normalizePath(p, winslash = "/", mustWork = TRUE),
                                                                       size_bytes = file.info(p)$size, sha256 = sha256(p))))
write_csv(input_manifest, file.path(OUT, "tables", "input_sha256_manifest.csv"))
environment_identity <- paste(c(as.character(getRversion()), vapply(c("data.table", "igraph", "Rcpp", "digest", "jsonlite", "AnnotationDbi", "org.Hs.eg.db"), function(p) paste0(p, "=", as.character(packageVersion(p))), character(1))), collapse = "|")
fingerprint <- digest::digest(paste(IMPL_VERSION, B, MAIN_SEED, MIN_BIN_SIZE, PPI_MIN_SUPPORT, environment_identity,
                                    paste(input_manifest$path, input_manifest$sha256, sep = "="), collapse = "\n"),
                              algo = "sha256", serialize = FALSE)

available_mem <- suppressWarnings(as.numeric(Sys.getenv("FIG6_AVAILABLE_RAM_BYTES", unset = NA_character_)))
if (!is.finite(available_mem)) {
  mem_cmd <- "Add-Type -AssemblyName Microsoft.VisualBasic; [Microsoft.VisualBasic.Devices.ComputerInfo]::new().AvailablePhysicalMemory"
  mem_out <- tryCatch(system2("powershell.exe", c("-NoProfile", "-Command", shQuote(mem_cmd)), stdout = TRUE, stderr = FALSE),
                      error = function(e) character())
  available_mem <- suppressWarnings(as.numeric(tail(mem_out, 1L)))
}
estimated_peak_bytes <- 4 * B * nrow(main_nodes) + 8 * min(250L, nrow(main_nodes)) * nrow(main_nodes) + 768 * 1024^2
if (!is.finite(available_mem)) stop("无法测量当前可用物理内存，停止距离零模型")
if (estimated_peak_bytes > 0.70 * available_mem) stop("预计峰值内存超过当前可用物理内存的70%")

cfg <- list(implementation = IMPL_VERSION, diseases = DISEASES, candidate_gene_set = "all unique Entrez in frozen selected_nodes.csv; no reward/status/topK filtering",
            ppi_source = PPI_FILE, ppi_filter = "number >= 2", edge_semantics = "undirected, unweighted; unique endpoints; self-loops removed",
            distance_background = "largest connected component; candidate input retained in full audit; source set is its intersection with this component",
            n_raw_ppi_records = nrow(ppi_raw), n_ppi_records_ge2 = nrow(ppi_ge2), n_unique_edges_ge2 = nrow(ppi_edges),
            n_ppi_nodes_ge2 = vcount(ppi_graph), n_components_ge2 = length(cc$csize), n_max_component_nodes = vcount(main_graph),
            n_max_component_edges = ecount(main_graph), max_component_id = as.integer(largest_id),
            drug_target_input = DRUG_FILE, n_drugs = length(DRUGS), n_unique_drug_target_relations = nrow(drug_raw),
            mapper = paste0("org.Hs.eg.db ", packageVersion("org.Hs.eg.db"), "; AnnotationDbi SYMBOL -> ENTREZID; ambiguous mappings excluded"),
            statistic = "mean over unique mapped drug targets in maximum component of min unweighted shortest-path distance to all PCSF-selected disease nodes intersected with maximum component",
            null_model = "B=1000 dual random; independent disease-source and drug-target sets sampled without replacement, exact matched counts per frozen degree bin; cross-set overlap allowed; entire maximum component eligible",
            degree_bins = paste0("exact-degree ties kept intact; sorted degree then Entrez; accumulate to >=", MIN_BIN_SIZE,
                                 "; final remainder joins previous bin; if whole graph smaller than threshold, one bin"),
            seed = MAIN_SEED, stable_seed = "first 16 hex chars of SHA256('42|Fig6allnodes|set_type|Disease|Drug') parsed as mt19937_64 seed",
            empirical_p = "lower-tail (1 + sum(P_null <= P_obs))/(B+1)", p_adjustment = "single global Benjamini-Hochberg over every finite empirical P across all 16 diseases and all drugs",
            R = as.character(getRversion()), packages = vapply(c("data.table", "igraph", "Rcpp", "digest", "jsonlite", "AnnotationDbi", "org.Hs.eg.db"),
                                                               function(p) as.character(packageVersion(p)), character(1)),
            B = B, main_seed = MAIN_SEED, fingerprint = fingerprint, estimated_peak_bytes = estimated_peak_bytes,
            measured_available_physical_bytes = available_mem, created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
write_json(cfg, file.path(OUT, "data", "run_config.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
log_line("Start; fingerprint=", fingerprint, "; max_component=", vcount(main_graph), " nodes/", ecount(main_graph),
         " edges; disease candidate instances=", nrow(candidate_nodes), "; drugs=", length(DRUGS),
         "; estimated peak bytes=", estimated_peak_bytes)

drug_targets <- split(drug_raw[in_max_component == TRUE & !is.na(Entrez), Entrez],
                      drug_raw[in_max_component == TRUE & !is.na(Entrez), Drug])
for (dr in DRUGS) if (is.null(drug_targets[[dr]])) drug_targets[[dr]] <- character()
results_by_disease <- vector("list", length(DISEASES)); names(results_by_disease) <- DISEASES
for (d in DISEASES) {
  log_line("Disease start: ", d)
  cp_path <- file.path(OUT, "checkpoints", paste0("checkpoint_", d, "_", substr(fingerprint, 1L, 16L), ".rds"))
  if (file.exists(cp_path)) {
    cp <- readRDS(cp_path)
    if (!identical(cp$fingerprint, fingerprint) || !identical(cp$disease, d) || nrow(cp$rows) != length(DRUGS) ||
        anyDuplicated(cp$rows, by = c("Disease", "Drug")) || !identical(as.character(cp$rows$Drug), DRUGS) ||
        !identical(dim(cp$null), c(length(DRUGS), B))) stop("Disease checkpoint schema/fingerprint mismatch: ", d)
    results_by_disease[[d]] <- cp$rows
    log_line("Disease resumed from validated checkpoint: ", d)
    next
  }
  full_ids <- unique(candidate_nodes[Disease == d, Entrez])
  source_ids <- full_ids[full_ids %chin% main_ids]
  source_idx <- unname(node_index[source_ids])
  rows <- vector("list", length(DRUGS))
  null_mat <- matrix(NA_real_, nrow = length(DRUGS), ncol = B,
                     dimnames = list(DRUGS, paste0("b", seq_len(B))))
  source_seed <- stable_seed("Fig6allnodes", "source", d)
  obs_min <- core_draws <- null_min <- NULL
  core_counts <- integer(length(bin_pools))
  if (length(full_ids) > 0L && length(source_idx) > 0L) {
    obs_d <- igraph::distances(main_graph, v = source_idx, to = seq_len(nrow(main_nodes)), mode = "all", weights = NA)
    obs_min <- apply(obs_d, 2L, min)
    core_counts <- match_counts(source_idx)
    if (any(core_counts > lengths(bin_pools))) stop("Disease source degree-bin count unavailable: ", d)
    core_draws <- fig6_sample_stratified(bin_pools, as.integer(core_counts), B, source_seed)
    ed <- as_edgelist(main_graph, names = FALSE)
    null_min <- fig6_batched_multisource_bfs(as.integer(ed[, 1L]), as.integer(ed[, 2L]), vcount(main_graph), core_draws)
    rm(obs_d, ed)
  }
  for (ji in seq_along(DRUGS)) {
    dr <- DRUGS[[ji]]
    target_ids <- unique(drug_targets[[dr]])
    target_idx <- unname(node_index[target_ids])
    d_obs <- null_mean <- null_sd <- z <- p_emp <- NA_real_
    status <- if (!length(full_ids)) "empty_candidate_set" else
      if (!length(source_idx)) "no_candidate_in_max_component" else
      if (!nrow(drug_raw[Drug == dr])) "no_target_relations" else
      if (!length(target_idx)) "no_mapped_target_in_max_component" else "complete"
    target_seed <- NA_character_
    complete_null <- FALSE
    target_counts <- integer(length(bin_pools))
    if (status == "complete") {
      d_obs <- mean(obs_min[target_idx])
      target_counts <- match_counts(target_idx)
      if (any(target_counts > lengths(bin_pools))) {
        status <- "degree_bin_sample_unavailable"
      } else {
        target_seed <- stable_seed("Fig6allnodes", "target", d, dr)
        null <- fig6_target_null_mean(null_min, bin_pools, as.integer(target_counts), B, target_seed)
        null_mat[ji, ] <- null
        complete_null <- length(null) == B && all(is.finite(null)) && is.finite(d_obs)
        if (!complete_null) {
          status <- "null_nonfinite"
        } else {
          null_mean <- mean(null); null_sd <- sd(null)
          p_emp <- (1 + sum(null <= d_obs)) / (B + 1)
          z <- if (is.finite(null_sd) && null_sd > 0) (d_obs - null_mean) / null_sd else NA_real_
          status <- if (is.finite(z)) "complete" else "null_sd_zero_z_undefined"
        }
      }
    }
    rows[[ji]] <- data.table(Disease = d, Drug = dr,
      network_graph_node_count = length(full_ids), candidate_unique_entrez_count = length(full_ids),
      candidate_in_max_component_n = length(source_idx), candidate_excluded_n = length(full_ids) - length(source_idx),
      n_drug_targets_in_max_component = length(target_idx), network_status = status,
      observed_mean_min_distance = d_obs, null_mean = null_mean, null_sd = null_sd, Z_dual_random = z,
      p_empirical = p_emp, q_global_BH = NA_real_, null_B = if (complete_null) B else 0L,
      null_complete_finite = complete_null, source_seed = if (length(source_idx)) source_seed else NA_character_,
      target_seed = target_seed, n_source_degree_bins = sum(core_counts > 0L),
      n_target_degree_bins = sum(target_counts > 0L),
      z_category = if (is.finite(z) && z < -1.96) "exploratory_dual_random_Z_below_minus_1.96" else
        if (is.finite(z) && z < -1.645) "exploratory_marginal_dual_random_Z" else
        if (is.finite(z)) "other_finite_Z" else "undefined_Z")
  }
  disease_rows <- rbindlist(rows, use.names = TRUE)
  if (nrow(disease_rows) != length(DRUGS) || anyDuplicated(disease_rows, by = c("Disease", "Drug")) ||
      !identical(as.character(disease_rows$Drug), DRUGS)) stop("Disease×Drug keys incomplete/duplicated: ", d)
  cp <- list(fingerprint = fingerprint, disease = d, rows = disease_rows, null = null_mat,
             core_draws = core_draws, source_seed = source_seed,
             candidate_ids_full = full_ids, candidate_ids_distance = source_ids,
             created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
  atomic_save_rds(cp, cp_path)
  results_by_disease[[d]] <- disease_rows
  rm(cp, null_mat, null_min, core_draws, obs_min); gc(verbose = FALSE)
  log_line("Disease complete: ", d, "; graph_nodes=", length(full_ids), "; distance_nodes=", length(source_idx),
           "; finite_P=", sum(is.finite(disease_rows$p_empirical)), "; checkpoint=", basename(cp_path))
}

result <- rbindlist(results_by_disease, use.names = TRUE)
if (nrow(result) != length(DISEASES) * length(DRUGS) || anyDuplicated(result, by = c("Disease", "Drug")))
  stop("最终Disease×Drug结果键不完整或重复")
valid_p <- is.finite(result$p_empirical)
result[valid_p, q_global_BH := p.adjust(p_empirical, method = "BH")]
write_csv(result, file.path(OUT, "tables", "all_disease_drug_network_proximity.csv"))
write_csv(result[network_status == "complete" | network_status == "null_sd_zero_z_undefined"],
          file.path(OUT, "tables", "all_evaluable_network_proximity.csv"))
write_csv(result[is.finite(Z_dual_random) & Z_dual_random < -1.96],
          file.path(OUT, "tables", "exploratory_dual_random_Z_below_minus_1_96.csv"))
write_csv(result[is.finite(q_global_BH) & q_global_BH < 0.05],
          file.path(OUT, "tables", "global_BH_q_below_0_05.csv"))
write_csv(result[is.finite(p_empirical), .(Disease, Drug, p_empirical)],
          file.path(OUT, "tables", "global_BH_valid_P_manifest.csv"))
write_csv(result[is.finite(Z_dual_random), .(Disease, Drug, Z_dual_random)],
          file.path(OUT, "tables", "global_Z_assessable_manifest.csv"))
write_csv(result[, .(n_drug_keys = .N, n_valid_P = sum(is.finite(p_empirical)),
                     n_valid_Z = sum(is.finite(Z_dual_random)), n_q_below_0_05 = sum(is.finite(q_global_BH) & q_global_BH < 0.05),
                     n_empty_candidate_set = sum(network_status == "empty_candidate_set"),
                     n_drug_pairs_with_distance_exclusions = sum(candidate_excluded_n > 0L),
                     candidate_nodes_excluded_from_distance = max(candidate_excluded_n)), by = Disease],
          file.path(OUT, "tables", "disease_drug_result_summary.csv"))
write_csv(result[, .(Disease, Drug, p_empirical, q_global_BH)][is.finite(p_empirical)][order(p_empirical)],
          file.path(OUT, "tables", "global_P_BH_order.csv"))

log_line("Global BH completed; keys=", nrow(result), "; valid P=", sum(valid_p),
         "; valid Z=", sum(is.finite(result$Z_dual_random)), "; q<0.05=", sum(result$q_global_BH < 0.05, na.rm = TRUE),
         "; empty candidate disease keys=", sum(result$network_status == "empty_candidate_set"))

# Rebuild manifest after all outputs are written so hashes point to exact inputs/code.
out_files <- list.files(OUT, recursive = TRUE, full.names = TRUE)
out_files <- out_files[file.info(out_files)$isdir %in% FALSE & !grepl("\\.tmp$", out_files)]
manifest <- list(status = "complete", fingerprint = fingerprint, disease_count = length(DISEASES), drug_count = length(DRUGS),
                 disease_drug_key_count = nrow(result), valid_empirical_p_count = sum(valid_p),
                 valid_z_count = sum(is.finite(result$Z_dual_random)), global_bh_q_lt_0_05_count = sum(result$q_global_BH < 0.05, na.rm = TRUE),
                 network_candidate_coverage = coverage, input_manifest = input_manifest,
                 outputs = lapply(out_files, function(p) list(path = normalizePath(p, winslash = "/", mustWork = TRUE),
                                                            size_bytes = file.info(p)$size, sha256 = sha256(p))),
                 finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
write_json(manifest, file.path(OUT, "run_manifest.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
log_line("Run manifest written; complete.")
