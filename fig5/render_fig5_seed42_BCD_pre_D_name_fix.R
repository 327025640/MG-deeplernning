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
for (path in c(file.path(FIG, "Fig5B_seed42_PCSF_weight_distribution.png"), file.path(FIG, "Fig5B_seed42_PCSF_weight_distribution.pdf"),
               file.path(FIG, "Fig5C_seed42_PCSF_candidate_network.png"), file.path(FIG, "Fig5C_seed42_PCSF_candidate_network.pdf"),
               file.path(FIG, "Fig5D_seed42_PCSF_GO_BP_barplot.png"), file.path(FIG, "Fig5D_seed42_PCSF_GO_BP_barplot.pdf"))) {
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
if (nrow(nodes) && (!all(c("Entrez", "gene_id", "PTG", "selected_runs", "core_10_of_10", "median_weight", "node_role") %in% names(nodes)) || anyDuplicated(nodes$Entrez))) stop("C图节点表字段或Entrez唯一性错误。")
if (nrow(edges) && (!all(c("Protein_A_Entrez_ID", "Protein_B_Entrez_ID", "selected_in_seed_count", "display_category") %in% names(edges)) ||
    any(edges$Protein_A_Entrez_ID == edges$Protein_B_Entrez_ID) || anyDuplicated(edges[, .(Protein_A_Entrez_ID, Protein_B_Entrez_ID)]))) stop("C图边表存在自环或重复无向边。")

# B: preserve the established ten-seed boxplot, using the original sum-normalized weights.
registry[, core_10_of_10 := as.logical(core_10_of_10)]
registry[, selected_runs := as.integer(selected_runs)]
registry <- registry[order(-selected_runs, gene_id)]
if (nrow(registry)) {
  weights[, gene_id := factor(gene_id, levels = registry$gene_id)]
  registry[, label := ifelse(core_10_of_10, paste0("<span style='color:", CORE_RED, "'>", gene_id, "</span>"), gene_id)]
  labels <- setNames(registry$label, registry$gene_id)
  max_y <- max(weights$weight)
  pB <- ggplot(weights, aes(x = gene_id, y = weight, color = PTG, group = gene_id)) +
    geom_boxplot(width = 0.62, fill = "#F2F2F2", linewidth = 0.48, outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.10, height = 0), shape = 21, fill = "white", size = 1.65, stroke = 0.48, alpha = 0.95) +
    scale_color_manual(values = PTG_COLORS, drop = FALSE, name = "PTG module") +
    scale_y_continuous(name = "Normalized mean absolute SHAP weight", expand = expansion(mult = c(0, 0.08))) +
    scale_x_discrete(name = NULL, labels = labels, expand = expansion(add = c(0.5, 0.5))) +
    coord_cartesian(ylim = c(0, max_y * 1.12)) +
    labs(title = "Seed42 PCSF model-evidence genes: SHAP weights across ten frozen models",
      subtitle = "All ten frozen seed weights are shown, including seeds without PCSF selection.",
      caption = "Red labels: seed42 model-evidence nodes selected in 10/10 frozen networks.") +
    theme_classic(base_family = "Arial") +
    theme(plot.background = element_rect(fill = "white", color = NA), panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major.y = element_line(color = "#E6E6E6", linewidth = 0.3), panel.grid.minor = element_blank(),
      axis.line = element_line(color = "#333333", linewidth = 0.35),
      axis.text.x = ggtext::element_markdown(family = "Arial", size = 7.2, color = "#000000", angle = 58, hjust = 1, vjust = 1),
      axis.ticks.x = element_blank(), axis.text.y = element_text(family = "Arial", size = 8, color = "#111111"),
      axis.title.y = element_text(family = "Arial", size = 9, color = "#111111", margin = margin(r = 5)),
      plot.title = element_text(family = "Arial", size = 11, color = "#111111", hjust = 0.5),
      plot.subtitle = element_text(family = "Arial", size = 8, color = "#444444", hjust = 0.5),
      plot.caption = element_text(family = "Arial", size = 7.5, color = "#333333", hjust = 0.98),
      legend.position = "bottom", legend.title = element_text(family = "Arial", size = 8), legend.text = element_text(family = "Arial", size = 7.5),
      plot.margin = margin(4, 5, 8, 5, unit = "mm"))
  widthB <- max(250, 90 + nrow(registry) * 5.5)
  ggsave(file.path(FIG, "Fig5B_seed42_PCSF_weight_distribution.pdf"), pB, device = grDevices::cairo_pdf,
         width = widthB, height = 150, units = "mm", family = "Arial", bg = "white", limitsize = FALSE)
  ggsave(file.path(FIG, "Fig5B_seed42_PCSF_weight_distribution.png"), pB,
         width = widthB, height = 150, units = "mm", dpi = 600, bg = "white", limitsize = FALSE)
} else {
  pB <- ggplot() + annotate("text", x = 0, y = 0, label = "No model-evaluable gene met the ≥6/10 PCSF selection rule", family = "Arial", size = 4) +
    labs(title = "Seed42 PCSF model-evidence genes: SHAP weights across ten frozen models") + theme_void(base_family = "Arial")
  ggsave(file.path(FIG, "Fig5B_seed42_PCSF_weight_distribution.pdf"), pB, device = grDevices::cairo_pdf, width = 250, height = 100, units = "mm", family = "Arial", bg = "white")
  ggsave(file.path(FIG, "Fig5B_seed42_PCSF_weight_distribution.png"), pB, width = 250, height = 100, units = "mm", dpi = 600, bg = "white")
}
fwrite(weights, file.path(SRC, "Fig5B_seed42_PCSF_weight_distribution_source.csv"), na = "NA", bom = TRUE)
write_parameters(file.path(SRC, "Fig5B_seed42_PCSF_weight_distribution_parameters.csv"), list(
  candidate_count = nrow(registry), core_10_of_10_count = sum(registry$core_10_of_10, na.rm = TRUE),
  weight_definition = "per-seed mean_abs_SHAP / sum over all 4,843 model features; no selection-to-zero substitution",
  seeds = "42-51", font = "Arial", png_dpi = 600, white_background = TRUE))

# C: established ring organization for the frozen seed42 reference PCSF forest only.
if (nrow(nodes)) {
  nodes[, core_10_of_10 := as.logical(core_10_of_10)]
  ptg_order <- c("PTG01", "PTG02", "PTG03", "PTG04")
  nodes[, ptg_sort := match(PTG, ptg_order)]
  nodes[is.na(ptg_sort), ptg_sort := length(ptg_order) + 1L]
  nodes <- nodes[order(ptg_sort, -selected_runs, Entrez)]
  nodes[, order_idx := seq_len(.N)]
  nnode <- nrow(nodes)
  nodes[, angle_rad := pi / 2 - 2 * pi * (order_idx - 1) / nnode]
  nodes[, `:=`(x = cos(angle_rad), y = sin(angle_rad))]
  valid_weight <- nodes$model_evaluable & is.finite(nodes$median_weight)
  node_colors <- ifelse(nodes$model_evaluable, unname(PTG_COLORS[match(nodes$PTG, names(PTG_COLORS))]), "#999999")
  node_colors[is.na(node_colors)] <- "#999999"
  nodes[, `:=`(node_color = node_colors, label_color = ifelse(core_10_of_10, CORE_RED, "#000000"))]
  if (any(valid_weight)) {
    weight_range <- range(nodes$median_weight[valid_weight], na.rm = TRUE)
    if (diff(weight_range) == 0) nodes[, plot_size := ifelse(valid_weight, 5.5, 2)] else
      nodes[, plot_size := ifelse(valid_weight, 2 + 7 * (median_weight - weight_range[1]) / diff(weight_range), 2)]
  } else nodes[, plot_size := 2]
  if (nnode == 1L) nodes[, `:=`(x = 0, y = 0)]
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
    geom_point(data = nodes, aes(x = x, y = y, size = plot_size, fill = node_color), shape = 21, color = "#FFFFFF", stroke = 0.5) +
    scale_fill_identity() + scale_size_identity(name = "Median normalized mean |SHAP|", breaks = weight_positions,
      labels = format(weight_breaks, trim = TRUE, scientific = FALSE), guide = guide_legend(override.aes = list(fill = "#777777", color = "#777777"))) +
    ggrepel::geom_text_repel(data = nodes, aes(x = x, y = y, label = ifelse(nzchar(gene_id), gene_id, ifelse(nzchar(official_symbols), official_symbols, Entrez)), color = label_color),
      family = "Arial", size = 3.1, seed = 42, max.overlaps = Inf, point.padding = grid::unit(3, "mm"),
      box.padding = 0.16, force = 2, min.segment.length = 0, segment.color = "#777777", segment.size = 0.25, show.legend = FALSE) +
    scale_color_identity() + coord_equal(clip = "off") +
    labs(title = "Seed42 reference PCSF network",
      subtitle = paste0(nnode, " nodes; ", nrow(draw_edges), " real PCSF backbone edges; 76 model-evidence nodes and 6 connectors")) +
    theme_void(base_family = "Arial") + theme(plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA), plot.title = element_text(size = 11, hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = 8, hjust = 0.5, color = "#444444"), legend.position = "bottom",
      legend.title = element_text(size = 8, color = "black"), legend.text = element_text(size = 7.5, color = "black"),
      plot.margin = margin(10, 25, 10, 25, unit = "mm"))
  ggsave(file.path(FIG, "Fig5C_seed42_PCSF_candidate_network.pdf"), pC, device = grDevices::cairo_pdf, width = 220, height = 190, units = "mm", family = "Arial", bg = "white")
  ggsave(file.path(FIG, "Fig5C_seed42_PCSF_candidate_network.png"), pC, width = 220, height = 190, units = "mm", dpi = 600, bg = "white")
} else {
  pC <- ggplot() + annotate("text", x = 0, y = 0, label = "No seed42 PCSF node available", family = "Arial", size = 4) +
    labs(title = "Seed42 PCSF candidate network") + theme_void(base_family = "Arial")
  ggsave(file.path(FIG, "Fig5C_seed42_PCSF_candidate_network.pdf"), pC, device = grDevices::cairo_pdf, width = 220, height = 150, units = "mm", family = "Arial", bg = "white")
  ggsave(file.path(FIG, "Fig5C_seed42_PCSF_candidate_network.png"), pC, width = 220, height = 150, units = "mm", dpi = 600, bg = "white")
}
fwrite(nodes, file.path(SRC, "Fig5C_seed42_PCSF_candidate_network_nodes_source.csv"), na = "NA", bom = TRUE)
fwrite(edges, file.path(SRC, "Fig5C_seed42_PCSF_candidate_network_edges_source.csv"), na = "NA", bom = TRUE)
write_parameters(file.path(SRC, "Fig5C_seed42_PCSF_candidate_network_parameters.csv"), list(
  node_rule = "All 82 nodes from the unperturbed seed42 b10/w3/g3 PCSF network; evidence genes and six zero-reward connectors retained separately",
  edges = "Exactly the 81 real selected_forest_edges from the seed42 PCSF forest; no induced edges or virtual-root edges",
  selected_edge = "all 81 displayed edges are the actual seed42 frozen PCSF backbone; cross-seed frequency is retained in source only and does not imply additional displayed PPI edges",
  node_size = "median per-seed sum-normalized mean_abs_SHAP; connector fixed minimum",
  layout = "Established ring layout; PTG order with connectors neutral; no directed arrows", seeds = "42-51", font = "Arial", png_dpi = 600))

# D: recompute GO BP over-representation using the frozen mapping and pre-PCSF PPI-mapped feature background.
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
mapping_meta <- jsonlite::fromJSON(file.path(ROOT, "data", "orgHs_mapping_resource_manifest.json"), simplifyVector = TRUE)
if (as.character(utils::packageVersion("org.Hs.eg.db")) != as.character(mapping_meta$version)) {
  stop("Installed org.Hs.eg.db version differs from the frozen mapping resource; no substitute annotation will be used.")
}
if (as.character(utils::packageVersion("org.Hs.eg.db")) != "3.19.1" || as.character(utils::packageVersion("GO.db")) != "3.20.0") {
  stop("GO resource versions differ from the accepted Fig5 contract (org.Hs.eg.db 3.19.1; GO.db 3.20.0); no substitute resource will be used.")
}
map42 <- mapping[seed == 42L]
mapped_background <- unique(map42[mapping_status %chin% c("official_unique", "alias_unique") & in_qualified_PPI == TRUE,
  .(gene_id = as.character(gene_id), ENTREZID = as.character(chosen_Entrez))])
if (nrow(mapped_background) && anyDuplicated(mapped_background$gene_id)) stop("GO背景模型Symbol一对一映射键重复。")
entrez_keys <- unique(mapped_background$ENTREZID)
annotation_raw <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = entrez_keys,
  keytype = "ENTREZID", columns = c("GOALL", "ONTOLOGYALL")))
annotation_raw <- as.data.table(annotation_raw)
bp_entrez <- unique(annotation_raw[ONTOLOGYALL == "BP" & !is.na(GOALL), .(ENTREZID = as.character(ENTREZID), GOID = as.character(GOALL))])
bp_map <- merge(mapped_background, bp_entrez, by = "ENTREZID", all = FALSE, sort = FALSE)
bp_map <- unique(bp_map[, .(gene_id, GOID)])
candidate_ids <- unique(as.character(registry$gene_id))
candidate_annotated <- intersect(candidate_ids, unique(bp_map$gene_id))
background_annotated <- unique(bp_map$gene_id)
N_background <- length(background_annotated)
n_candidate <- length(candidate_annotated)
term_table <- bp_map[, .(K_term = uniqueN(gene_id)), by = GOID]
go_terms <- if (nrow(term_table)) as.data.table(suppressMessages(AnnotationDbi::select(GO.db::GO.db,
  keys = term_table$GOID, keytype = "GOID", columns = "TERM"))) else data.table(GOID = character(), TERM = character())
go_terms <- unique(go_terms[!is.na(TERM), .(GOID = as.character(GOID), term_name = as.character(TERM))])
term_table <- merge(term_table, go_terms, by = "GOID", all.x = TRUE, sort = FALSE)
if (anyNA(term_table$term_name) || any(!nzchar(term_table$term_name))) {
  missing_go <- term_table[is.na(term_name) | !nzchar(term_name), GOID]
  stop("GO.db is missing official labels for tested terms; refusing placeholder labels: ", paste(missing_go, collapse = ","))
}
candidate_map <- bp_map[gene_id %in% candidate_annotated,
  .(candidate_genes = paste(sort(unique(gene_id)), collapse = ";"), k_hit = uniqueN(gene_id)), by = GOID]
tests <- merge(term_table, candidate_map, by = "GOID", all.x = TRUE, sort = FALSE)
tests[is.na(k_hit), `:=`(k_hit = 0L, candidate_genes = "")]
tests[, `:=`(N_background = as.integer(N_background), n_candidate_annotated = as.integer(n_candidate), K_term = as.integer(K_term), k_hit = as.integer(k_hit))]
tests[, P := stats::phyper(k_hit - 1L, K_term, N_background - K_term, n_candidate_annotated, lower.tail = FALSE)]
tests[, q_BH_all_BP := stats::p.adjust(P, method = "BH")]
tests[, significant_BH_005 := q_BH_all_BP < 0.05]
setorder(tests, q_BH_all_BP, P, GOID)
tests[, `:=`(displayed_on_main_figure = FALSE, display_basis = "not selected for main figure")]
sig <- tests[significant_BH_005 == TRUE][order(q_BH_all_BP, P, GOID)]
if (nrow(sig)) {
  display <- head(sig, 20L); mode <- "BH-significant GO BP terms"
  basis <- "Top 20 q<0.05 terms ordered by BH q, raw P, GOID; complete family retained"
} else {
  display <- head(tests[k_hit > 0L][order(-k_hit, GOID)], 20L)
  mode <- "Descriptive GO BP terms; no term met BH q<0.05"
  basis <- "Top 20 by candidate hit count then GOID; complete family retained"
}
if (nrow(display)) tests[GOID %in% display$GOID, `:=`(displayed_on_main_figure = TRUE, display_basis = basis)]
tests[, neg_log10_BH_q := -log10(pmax(q_BH_all_BP, .Machine$double.xmin))]
fwrite(tests, file.path(SRC, "Fig5D_seed42_PCSF_GO_BP_full_test_table.csv"), na = "NA", bom = TRUE)
background_audit <- copy(map42[, .(gene_id = as.character(gene_id), mapping_status, mapping_source,
  chosen_Entrez = as.character(chosen_Entrez), in_qualified_PPI)])
background_audit[, `:=`(in_GO_BP_background = gene_id %in% background_annotated,
  in_seed42_model_evidence_set = gene_id %in% candidate_ids,
  excluded_reason = fifelse(!in_qualified_PPI, "not_uniquely_mapped_to_qualified_PPI",
    fifelse(!gene_id %in% background_annotated, "no_GO_BP_annotation", "included_in_GO_BP_background")))]
fwrite(background_audit,
  file.path(SRC, "Fig5D_seed42_PCSF_GO_BP_background_audit.csv"), na = "NA", bom = TRUE)
display[, `:=`(display_rank = seq_len(.N), term_label = term_name,
  term_label_wrapped = vapply(term_name, function(x) paste(strwrap(x, width = 38L), collapse = "\n"), character(1)))]
if (nrow(display)) {
  display[, term_axis := factor(paste0("term_", display_rank), levels = rev(paste0("term_", display_rank)))]
  axis_labels <- setNames(display$term_label_wrapped, as.character(display$term_axis))
  pD <- ggplot(display, aes(x = neg_log10_BH_q, y = term_axis)) +
    geom_col(width = 0.7, fill = "#2166AC", colour = "#184F83", linewidth = 0.25) +
    scale_x_continuous(name = expression(-log[10]("BH-adjusted P")), breaks = scales::pretty_breaks(n = 6), expand = expansion(mult = c(0, 0.025))) +
    scale_y_discrete(name = NULL, labels = axis_labels, drop = FALSE) +
    labs(title = "GO Biological Process enrichment") + theme_minimal(base_family = "Arial", base_size = 9) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(colour = "black", size = 8, lineheight = 0.95),
      axis.text.x = element_text(colour = "black", size = 8), axis.title.x = element_text(colour = "black", size = 9, margin = margin(t = 7)),
      plot.title = element_text(colour = "black", size = 12, face = "plain"), plot.margin = margin(9, 12, 9, 9, unit = "pt"))
} else {
  pD <- ggplot() + annotate("text", x = 0, y = 0, label = "No GO BP terms among seed42 model-evidence genes", family = "Arial", size = 4) +
    labs(title = "GO Biological Process enrichment") + theme_void(base_family = "Arial")
}
heightD <- max(90, 40 + 7.5 * nrow(display))
ggsave(file.path(FIG, "Fig5D_seed42_PCSF_GO_BP_barplot.png"), pD, width = 265, height = heightD, units = "mm", dpi = 600, bg = "white", limitsize = FALSE)
ggsave(file.path(FIG, "Fig5D_seed42_PCSF_GO_BP_barplot.pdf"), pD, width = 265, height = heightD, units = "mm", device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
fwrite(display, file.path(SRC, "Fig5D_seed42_PCSF_GO_BP_barplot_source.csv"), na = "NA", bom = TRUE)
write_parameters(file.path(SRC, "Fig5D_seed42_PCSF_GO_BP_parameters.csv"), list(
  candidate_count = length(candidate_ids), candidate_annotated = n_candidate,
  background = "all 4,843 model features uniquely mapped to eligible PPI Entrez nodes and annotated with GO BP before PCSF selection",
  background_genes = N_background, tested_terms = nrow(tests), significant_q_lt_005 = sum(tests$significant_BH_005),
  test = "one-sided hypergeometric over-representation; BH over all testable BP terms",
  mapping = "same frozen local Org.Hs.eg.db 3.19.1-derived official/alias mappings used for PCSF input; GO annotations retrieved from installed org.Hs.eg.db GOALL/ONTOLOGYALL",
  display_mode = mode, display_basis = basis, font = "Arial", png_dpi = 600))
log_line(paste0("B/C/D complete; seed42 model-evidence genes=", nrow(registry), "; network display nodes=", nrow(nodes),
  "; GO background=", N_background, "; tested terms=", nrow(tests), "; significant=", sum(tests$significant_BH_005)))

output_files <- c(file.path(FIG, c("Fig5B_seed42_PCSF_weight_distribution.png","Fig5B_seed42_PCSF_weight_distribution.pdf",
  "Fig5C_seed42_PCSF_candidate_network.png","Fig5C_seed42_PCSF_candidate_network.pdf",
  "Fig5D_seed42_PCSF_GO_BP_barplot.png","Fig5D_seed42_PCSF_GO_BP_barplot.pdf")),
  file.path(SRC, c("Fig5B_seed42_PCSF_weight_distribution_source.csv","Fig5B_seed42_PCSF_weight_distribution_parameters.csv",
  "Fig5C_seed42_PCSF_candidate_network_nodes_source.csv","Fig5C_seed42_PCSF_candidate_network_edges_source.csv",
  "Fig5C_seed42_PCSF_candidate_network_parameters.csv","Fig5D_seed42_PCSF_GO_BP_full_test_table.csv",
  "Fig5D_seed42_PCSF_GO_BP_background_audit.csv","Fig5D_seed42_PCSF_GO_BP_barplot_source.csv",
  "Fig5D_seed42_PCSF_GO_BP_parameters.csv")))
if (!all(file.exists(output_files))) stop("Renderer output inventory incomplete: ",paste(output_files[!file.exists(output_files)],collapse="; "))
hash_inventory <- data.table(path=output_files,sha256=vapply(output_files,sha256_file,character(1)),size_bytes=file.info(output_files)$size)
fwrite(hash_inventory,file.path(SRC,"rendered_output_sha256_manifest.csv"),na="NA",bom=TRUE)


