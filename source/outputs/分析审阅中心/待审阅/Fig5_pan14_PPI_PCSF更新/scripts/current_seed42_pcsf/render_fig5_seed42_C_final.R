#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 999)
suppressWarnings(invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8")))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

ROOT <- "E:/micro_output/分析审阅中心/待审阅/Fig5_pan14_PPI_PCSF更新"
MODEL_ROOT <- "E:/micro_output/分析审阅中心/待审阅/PTG02_CD_分支疾病与调控候选/ptg14_training_v2"
FIG5_OLD <- "E:/micro_output/分析审阅中心/待审阅/Fig5_主分析与流程图收尾"
SRC_BASE <- file.path(ROOT, "source_data")
SRC <- file.path(ROOT, "source_data", "current_seed42_pcsf")
FIG <- file.path(ROOT, "figures", "current_seed42_pcsf")
LOG <- file.path(ROOT, "logs", "current_seed42_pcsf", "render_fig5_seed42_BCD.log")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG), recursive = TRUE, showWarnings = FALSE)
PTG_COLORS <- c(PTG01 = "#2166AC", PTG02 = "#B2182B", PTG03 = "#1B9E77", PTG04 = "#7570B3")
CORE_RED <- "#D73027"

need <- c(file.path(SRC, "seed42_model_evidence_76_registry.csv"),
          file.path(SRC, "seed42_B_weights_76x10.csv"),
          file.path(SRC_BASE, "all_seed_full_gene_weights.csv"),
          file.path(SRC, "seed42_C_nodes_82.csv"),
          file.path(SRC, "seed42_C_backbone_edges_81.csv"),
          file.path(SRC_BASE, "all_gene_mapping_audit.csv"),
          file.path(ROOT, "data", "ppi_edges_ge2.tsv"))
if (!all(file.exists(need))) stop("缺少PCSF后候选或网络输入：", paste(need[!file.exists(need)], collapse = "; "))
if (!requireNamespace("ggrepel", quietly = TRUE) || !requireNamespace("ggtext", quietly = TRUE) ||
    !requireNamespace("digest", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
    !requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace("GO.db", quietly = TRUE) ||
    !requireNamespace("jsonlite", quietly = TRUE) || !requireNamespace("scales", quietly = TRUE)) {
  stop("缺少已安装绘图或GO注释依赖；不自动安装。")
}
sha256_file <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) stop("Cannot hash missing file: ", path)
  normalized <- normalizePath(path, winslash = "\\", mustWork = TRUE)
  out <- suppressWarnings(system2("certutil.exe",
    c("-hashfile", shQuote(normalized), "SHA256"), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  hash_lines <- grepl("^[A-Fa-f0-9]{64}$", out, perl = TRUE, useBytes = TRUE)
  hashes <- tolower(out[hash_lines])
  if (!is.null(status) && status != 0L) stop("certutil SHA-256 failed for ", path, ": ", paste(out, collapse = " | "))
  if (length(hashes) != 1L) stop("Could not parse streaming SHA-256 for ", path, ": ", paste(out, collapse = " | "))
  hashes[[1L]]
}
log_line <- function(x) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", x)
  cat(line, "\n", file = LOG, append = TRUE)
}
write_parameters <- function(path, data) data.table::fwrite(data.table(parameter = names(data), value = as.character(unlist(data))), path, bom = TRUE)
for (path in c(file.path(FIG, "Fig5C_seed42_PCSF_candidate_network_final.png"), file.path(FIG, "Fig5C_seed42_PCSF_candidate_network_final.pdf"))) {
  if (file.exists(path)) stop("拒绝覆盖既有输出：", path)
}
registry <- fread(need[1], na.strings = c("", "NA"), encoding = "UTF-8")
weights <- fread(need[2], na.strings = c("", "NA"), encoding = "UTF-8")
full_weights <- fread(need[3], na.strings = c("", "NA"), encoding = "UTF-8")
nodes <- fread(need[4], na.strings = c("", "NA"), encoding = "UTF-8")
edges <- fread(need[5], na.strings = c("", "NA"), encoding = "UTF-8")
mapping <- fread(need[6], na.strings = c("", "NA"), encoding = "UTF-8")
if (nrow(registry) && (!all(c("gene_id", "Entrez", "selected_runs", "core_10_of_10", "median_weight", "PTG") %in% names(registry)) || anyDuplicated(registry$gene_id))) stop("seed42模型证据基因表字段或唯一性错误。")
if (nrow(weights) && (!all(c("gene_id", "seed", "weight", "mean_abs_SHAP", "selected_this_seed") %in% names(weights)) ||
    anyDuplicated(weights[, .(gene_id, seed)]))) stop("B图权重长表字段或 gene×seed 唯一键错误。")
if (nrow(weights) && (any(!is.finite(weights$weight)) || any(weights$weight < 0) || any(weights$seed < 42 | weights$seed > 51))) stop("B图含无效权重或seed。")
if (nrow(weights) != nrow(registry) * 10L) stop("每个模型候选必须包含十个真实seed权重。")
if (nrow(full_weights) != 48430L || anyDuplicated(full_weights[, .(gene_id, seed)])) stop("十seed全基因权重表必须有48,430个唯一gene×seed观测。")
for (seed in 42:51) {
  seed_value <- as.integer(seed)
  fw <- full_weights[full_weights$seed == seed_value]
  expected_weight <- fw$mean_abs_SHAP / sum(fw$mean_abs_SHAP)
  if (nrow(fw) != 4843L || any(!is.finite(expected_weight)) || max(abs(expected_weight - fw$weight)) > 1e-12 || abs(sum(fw$weight) - 1) > 1e-9) {
    stop("seed ", seed, "全4843特征mean_abs_SHAP归一化权重公式核对失败。")
  }
}
full_candidate_weights <- full_weights[gene_id %in% registry$gene_id]
if (nrow(full_candidate_weights) != nrow(weights) || (nrow(weights) > 0L && max(abs(full_candidate_weights[weights, on = .(gene_id, seed), x.weight - i.weight])) > 1e-12)) {
  stop("B图绘图表权重与全基因冻结权重不一致。")
}
if (nrow(nodes) && (!all(c("Entrez", "gene_id", "official_symbols", "PTG", "selected_runs", "core_10_of_10", "median_weight", "node_role", "model_evaluable") %in% names(nodes)) || anyDuplicated(nodes$Entrez))) stop("C图节点表字段或Entrez唯一性错误。")
if (nrow(edges) && (!all(c("Protein_A_Entrez_ID", "Protein_B_Entrez_ID", "selected_in_seed_count", "display_category") %in% names(edges)) ||
    any(edges$Protein_A_Entrez_ID == edges$Protein_B_Entrez_ID) || anyDuplicated(edges[, .(Protein_A_Entrez_ID, Protein_B_Entrez_ID)]))) stop("C图边表存在自环或重复无向边。")

# C: established ring organization for the frozen seed42 reference PCSF forest only.
if (nrow(nodes) != 82L || nrow(edges) != 81L) stop("C-only renderer requires the frozen seed42 82-node/81-backbone input")
if (nrow(nodes)) {
  nodes[, `:=`(core_10_of_10 = as.logical(core_10_of_10), model_evaluable = as.logical(model_evaluable))]
  if (anyNA(nodes$core_10_of_10) || anyNA(nodes$model_evaluable) || anyNA(nodes$node_role) || any(!nzchar(nodes$node_role)) ||
      !all(nodes$model_evaluable == (nodes$node_role == "model_evidence")) ||
      any(!nodes$node_role %in% c("model_evidence", "zero_reward_PPI_connector")) || anyNA(nodes$official_symbols)) stop("Node evidence/core fields are incomplete or inconsistent.")
  ptg_order <- c("PTG01", "PTG02", "PTG03", "PTG04")
  nodes[, ptg_sort := match(PTG, ptg_order)]
  nodes[is.na(ptg_sort), ptg_sort := length(ptg_order) + 1L]
  nodes <- nodes[order(ptg_sort, -selected_runs, Entrez)]
  nodes[, order_idx := seq_len(.N)]
  nnode <- nrow(nodes)
  nodes[, angle_rad := pi / 2 - 2 * pi * (order_idx - 1) / nnode]
  nodes[, `:=`(x = cos(angle_rad), y = sin(angle_rad))]
  if (any(!is.finite(nodes[model_evaluable == TRUE, median_weight])) || any(nodes[model_evaluable == TRUE, median_weight] < 0)) stop("Model-evidence median weights must be finite and non-negative.")
  valid_weight <- nodes$model_evaluable & is.finite(nodes$median_weight)
  node_colors <- ifelse(nodes$model_evaluable, unname(PTG_COLORS[match(nodes$PTG, names(PTG_COLORS))]), "#999999")
  node_colors[is.na(node_colors)] <- "#999999"
  nodes[, `:=`(node_color = node_colors, label_color = ifelse(core_10_of_10, CORE_RED, "#000000"), label_nudge_x = 0.22 * x, label_nudge_y = 0.22 * y)]
  if (any(nodes$model_evaluable & is.na(nodes$PTG))) stop("Model-evidence node has no frozen PTG class for the legend.")
  nodes[, node_fill_group := fifelse(model_evaluable, as.character(PTG), "Connector / no model SHAP")]
  # Fixed local label offsets for the crowded lower arc; node coordinates and ring layout remain unchanged.
  offset_map <- data.table(gene_id = c("ATP5F1E", "ATP5MC2", "EEF1G", "FOS", "CCL3", "CCL4", "VIM", "ATP5MG"),
    label_nudge_x = c(-0.10, 0.16, 0.13, 0.10, -0.10, -0.08, -0.03, 0),
    label_nudge_y = c(0.12, -0.08, 0.06, -0.10, 0.13, -0.13, 0.14, -0.15))
  if (anyDuplicated(offset_map$gene_id) || !all(offset_map$gene_id %in% nodes[core_10_of_10 == TRUE, gene_id])) stop("Label offset IDs must be unique seed42 core nodes")
  offset_map[, label_nudge_x := label_nudge_x + 0.22 * nodes$x[match(gene_id, nodes$gene_id)]]
  offset_map[, label_nudge_y := label_nudge_y + 0.22 * nodes$y[match(gene_id, nodes$gene_id)]]
  nodes[offset_map, on = "gene_id", `:=`(label_nudge_x = i.label_nudge_x, label_nudge_y = i.label_nudge_y)]
  if (anyNA(offset_map$label_nudge_x) || anyNA(offset_map$label_nudge_y)) stop("Label offset vectors must be finite and complete")
  offset_rows <- nodes[match(offset_map$gene_id, gene_id)]
  if (nrow(offset_rows) != nrow(offset_map) || !identical(as.character(offset_rows$gene_id), as.character(offset_map$gene_id)) || any(offset_rows$label_nudge_x != offset_map$label_nudge_x) || any(offset_rows$label_nudge_y != offset_map$label_nudge_y)) stop("Label offsets were not aligned to their intended nodes")
  offset_map[, `:=`(Entrez = offset_rows$Entrez, label = ifelse(!is.na(offset_rows$gene_id) & nzchar(as.character(offset_rows$gene_id)), as.character(offset_rows$gene_id), ifelse(!is.na(offset_rows$official_symbols) & nzchar(as.character(offset_rows$official_symbols)), as.character(offset_rows$official_symbols), as.character(offset_rows$Entrez))), core_10_of_10 = offset_rows$core_10_of_10, label_color = offset_rows$label_color)]
  if (anyDuplicated(offset_map$Entrez) || any(offset_map$core_10_of_10 != TRUE) || any(offset_map$label_color != CORE_RED)) stop("Every offset must map uniquely to a red 10/10 core label")
  if (any(valid_weight)) {
    weight_range <- range(nodes$median_weight[valid_weight], na.rm = TRUE)
    if (diff(weight_range) == 0) nodes[, plot_size := ifelse(valid_weight, 5.5, 2)] else
      nodes[, plot_size := ifelse(valid_weight, 2 + 7 * (median_weight - weight_range[1]) / diff(weight_range), 2)]
  } else nodes[, plot_size := 2]
  if (nnode == 1L) nodes[, `:=`(x = 0, y = 0)]
  if (sum(nodes$core_10_of_10) != 15L || sum(nodes$model_evaluable) != 76L || sum(!nodes$model_evaluable) != 6L || anyNA(nodes$label_color)) stop("Frozen core/evidence/connector counts differ from seed42 contract")
  if (nrow(edges) != 81L || any(edges$selected_in_seed_count < 1L)) stop("C must contain exactly the 81 frozen seed42 backbone edges; support frequency is annotation only")
  edges[, `:=`(display_category = "seed42_reference_PCSF_backbone", display_color = "#555555", display_linetype = "solid")]
  edge_key <- function(a,b) paste(pmin(as.numeric(a),as.numeric(b)),pmax(as.numeric(a),as.numeric(b)),sep="-")
  if (nrow(edges)!=81L || uniqueN(edge_key(edges$Protein_A_Entrez_ID,edges$Protein_B_Entrez_ID))!=81L ||
      !all(c(edges$Protein_A_Entrez_ID,edges$Protein_B_Entrez_ID) %in% nodes$Entrez)) stop("C diagram must use exactly 81 unique frozen seed42 backbone edges with valid node endpoints")
  edge_i <- match(edges$Protein_A_Entrez_ID, nodes$Entrez); edge_j <- match(edges$Protein_B_Entrez_ID, nodes$Entrez)
  draw_edges <- edges[!is.na(edge_i) & !is.na(edge_j)]
  if (nrow(draw_edges)!=81L) stop("A seed42 backbone edge was lost during coordinate joining")
  ei <- edge_i[!is.na(edge_i) & !is.na(edge_j)]; ej <- edge_j[!is.na(edge_i) & !is.na(edge_j)]
  if (any(valid_weight)) {
    size_range_weight <- range(nodes$median_weight[valid_weight], na.rm = TRUE)
    weight_breaks <- if (diff(size_range_weight) > 0) pretty(size_range_weight, n = 3) else size_range_weight[1]
    weight_breaks <- weight_breaks[weight_breaks >= size_range_weight[1] & weight_breaks <= size_range_weight[2]]
    if (!length(weight_breaks)) weight_breaks <- size_range_weight[1]
    weight_positions <- if (diff(size_range_weight) > 0) 2 + 7 * (weight_breaks - size_range_weight[1]) / diff(size_range_weight) else rep(5.5, length(weight_breaks))
  } else {
    weight_breaks <- numeric(); weight_positions <- numeric()
  }
  pC <- ggplot() +
    geom_segment(data = data.frame(x = nodes$x[ei], y = nodes$y[ei], xend = nodes$x[ej], yend = nodes$y[ej]),
      aes(x = x, y = y, xend = xend, yend = yend), color = "#555555", alpha = 0.78, linewidth = 0.28) +
    geom_point(data = nodes, aes(x = x, y = y, size = plot_size, fill = node_fill_group), shape = 21, color = "#FFFFFF", stroke = 0.5) +
    scale_fill_manual(name = "PTG / node role",
      values = c(PTG_COLORS, "Connector / no model SHAP" = "#999999"),
      breaks = c("PTG01", "PTG02", "PTG03", "PTG04", "Connector / no model SHAP"),
      labels = c("PTG01", "PTG02", "PTG03", "PTG04", "PPI connector (no model SHAP)"),
      drop = FALSE, guide = guide_legend(override.aes = list(shape = 21, size = 3, color = "#FFFFFF"))) + scale_size_identity(name = "Median normalized mean |SHAP|", breaks = weight_positions,
      labels = format(weight_breaks, trim = TRUE, scientific = FALSE), guide = guide_legend(override.aes = list(fill = "#777777", color = "#777777"))) +
    ggrepel::geom_text_repel(data = nodes, aes(x = x, y = y, label = ifelse(!is.na(gene_id) & nzchar(as.character(gene_id)), as.character(gene_id), ifelse(!is.na(official_symbols) & nzchar(as.character(official_symbols)), as.character(official_symbols), as.character(Entrez))), color = label_color),
      family = "Arial", size = 3.1, seed = 42, max.overlaps = Inf, max.iter = 100000, max.time = 30, point.padding = grid::unit(2, "mm"), nudge_x = nodes$label_nudge_x, nudge_y = nodes$label_nudge_y,
      box.padding = 0.24, force = 5, force_pull = 0.08, min.segment.length = 0, segment.color = "#777777", segment.size = 0.25, show.legend = FALSE) +
    scale_color_identity() + coord_equal(clip = "off") +
    labs(title = "Seed42 reference PCSF network",
      subtitle = paste0(nnode, " nodes; ", nrow(draw_edges), " real PCSF backbone edges; 76 model-evidence nodes and 6 connectors"),
      caption = "Fill: PTG class; grey nodes: PPI connectors without model SHAP evidence. Size: median normalized mean |SHAP|; red labels: 10/10 nodes across frozen seeds.") +
    theme_void(base_family = "Arial") + theme(plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA), plot.title = element_text(size = 11, hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = 8, hjust = 0.5, color = "#444444"),
      plot.caption = element_text(size = 7, hjust = 0.5, color = "#333333"),
      legend.position = "bottom", legend.box = "vertical",
      legend.title = element_text(size = 8, color = "black"), legend.text = element_text(size = 7.5, color = "black"),
      plot.margin = margin(10, 25, 10, 25, unit = "mm"))
  ggsave(file.path(FIG, "Fig5C_seed42_PCSF_candidate_network_final.pdf"), pC, device = grDevices::cairo_pdf, width = 250, height = 210, units = "mm", family = "Arial", bg = "white")
  ggsave(file.path(FIG, "Fig5C_seed42_PCSF_candidate_network_final.png"), pC, width = 250, height = 210, units = "mm", dpi = 600, bg = "white")
} else {
  pC <- ggplot() + annotate("text", x = 0, y = 0, label = "No seed42 PCSF node available", family = "Arial", size = 4) +
    labs(title = "Seed42 PCSF candidate network") + theme_void(base_family = "Arial")
  ggsave(file.path(FIG, "Fig5C_seed42_PCSF_candidate_network_final.pdf"), pC, device = grDevices::cairo_pdf, width = 250, height = 210, units = "mm", family = "Arial", bg = "white")
  ggsave(file.path(FIG, "Fig5C_seed42_PCSF_candidate_network_final.png"), pC, width = 250, height = 210, units = "mm", dpi = 600, bg = "white")
}
fwrite(nodes, file.path(SRC, "Fig5C_seed42_PCSF_candidate_network_final_nodes_source.csv"), na = "NA", bom = TRUE)
fwrite(offset_map, file.path(SRC, "Fig5C_seed42_PCSF_candidate_network_final_label_offsets.csv"), na = "NA", bom = TRUE)
label_positions <- nodes[, .(Entrez, gene_id, official_symbols, label = ifelse(!is.na(gene_id) & nzchar(as.character(gene_id)), as.character(gene_id), ifelse(!is.na(official_symbols) & nzchar(as.character(official_symbols)), as.character(official_symbols), as.character(Entrez))),
  label_color, core_10_of_10, label_nudge_x, label_nudge_y, node_x = x, node_y = y,
  label_target_x = x + label_nudge_x, label_target_y = y + label_nudge_y)]
label_positions[, label_rule := "All 82 network nodes labeled; v3 gene labels preserved; outward radial nudge plus frozen local core offsets; ggrepel final placement is deterministic at seed 42."]
if (nrow(label_positions) != 82L || anyNA(label_positions$label) || any(!nzchar(label_positions$label)) || anyDuplicated(label_positions$Entrez) || anyDuplicated(label_positions$label) ||
    sum(label_positions$core_10_of_10, na.rm = TRUE) != 15L) stop("Label position table must cover 82 unique labels and preserve the 15-node core.")
fwrite(label_positions, file.path(SRC, "Fig5C_seed42_PCSF_candidate_network_final_label_positions.csv"), na = "NA", bom = TRUE)
fwrite(edges, file.path(SRC, "Fig5C_seed42_PCSF_candidate_network_final_edges_source.csv"), na = "NA", bom = TRUE)
fwrite(label_positions[, .(Entrez, gene_id, label, label_color, core_10_of_10, label_rule)],
  file.path(SRC, "Fig5C_seed42_PCSF_candidate_network_final_label_selection.csv"), na = "NA", bom = TRUE)
write_parameters(file.path(SRC, "Fig5C_seed42_PCSF_candidate_network_final_parameters.csv"), list(
  node_rule = "All 82 nodes from the unperturbed seed42 b10/w3/g3 PCSF network; evidence genes and six zero-reward connectors retained separately",
  edges = "Exactly the 81 real selected_forest_edges from the seed42 PCSF forest; no induced edges or virtual-root edges",
  selected_edge = "all 81 displayed edges are the actual seed42 frozen PCSF backbone; cross-seed frequency is retained in source only and does not imply additional displayed PPI edges",
  node_size = "median per-seed sum-normalized mean_abs_SHAP; connector fixed minimum",
  layout = "Established ring layout retained exactly; only label outward nudge/local core label offsets changed; ggrepel max.overlaps=Inf, seed=42",
  label_rule = "All 82 nodes labeled; 15 10/10 network nodes red; connectors remain neutral gray and are not treated as model evidence",
  legend = "PTG fill colors; gray PPI connectors with no model SHAP; median normalized SHAP size; red names mean selected in all 10 frozen seed networks",
  seeds = "42-51", font = "Arial", png_dpi = 600))

log_line("C-only label-layout revision complete; fixed seed42 nodes/edges and ring coordinates retained; lower-arc label offsets are recorded in source; ggrepel max.iter=100000, max.time=30, force=5; output canvas 250x210 mm")
output_files <- c(file.path(FIG, c("Fig5C_seed42_PCSF_candidate_network_final.png", "Fig5C_seed42_PCSF_candidate_network_final.pdf")),
  file.path(SRC, c("Fig5C_seed42_PCSF_candidate_network_final_nodes_source.csv", "Fig5C_seed42_PCSF_candidate_network_final_edges_source.csv",
    "Fig5C_seed42_PCSF_candidate_network_final_parameters.csv", "Fig5C_seed42_PCSF_candidate_network_final_label_offsets.csv",
    "Fig5C_seed42_PCSF_candidate_network_final_label_positions.csv", "Fig5C_seed42_PCSF_candidate_network_final_label_selection.csv")))
if (!all(file.exists(output_files))) stop("C-only renderer output inventory incomplete")
hash_inventory <- data.table(path=output_files, sha256=vapply(output_files, sha256_file, character(1)), size_bytes=file.info(output_files)$size)
fwrite(hash_inventory, file.path(SRC, "Fig5C_seed42_PCSF_candidate_network_final_output_sha256.csv"), na="NA", bom=TRUE)



