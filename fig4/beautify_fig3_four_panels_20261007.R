# Local layout refinement for four collected Fig3 panels.
# Reuses frozen sources only; does not rerun trajectory, clustering, scoring,
# module assignment, or any statistical analysis.

options(stringsAsFactors = FALSE, scipen = 999)
invisible(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(monocle)
  library(igraph)
  library(dplyr)
  library(scales)
})

OUT <- Sys.getenv("FIG3_FOUR_OUT",
                  unset = "D:/data/microglia_new/_staging_fig3_four_panel_20261007")
FIG <- file.path(OUT, "figures")
SRC <- file.path(OUT, "source_data")
LOG <- file.path(OUT, "logs")
for (d in c(OUT, FIG, SRC, LOG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

FONT <- "Arial"
GROUP_COLS <- c(Ctrl = "#4DAF4A", Neurological = "#4477AA", Tumor = "#B6362D")
PTG_COLS <- c(PTG01 = "#4477AA", PTG02 = "#EE6677", PTG03 = "#228833", PTG04 = "#AA3377")

save_panel <- function(plot_obj, stem, width_mm, height_mm) {
  png_path <- file.path(FIG, paste0(stem, ".png"))
  pdf_path <- file.path(FIG, paste0(stem, ".pdf"))
  grDevices::png(png_path, width = width_mm / 25.4, height = height_mm / 25.4,
                 units = "in", res = 600, type = "cairo", family = FONT,
                 bg = "white")
  print(plot_obj)
  grDevices::dev.off()
  grDevices::cairo_pdf(pdf_path, width = width_mm / 25.4,
                        height = height_mm / 25.4, family = FONT,
                        bg = "white", onefile = TRUE)
  print(plot_obj)
  grDevices::dev.off()
  invisible(c(png = png_path, pdf = pdf_path))
}

###############################################################################
# A. C2 segment composition: frozen cell-level denominators and proportions.
###############################################################################
c2_in <- file.path(OUT, "input_source", "C2_segment_group3_composition_source.csv")
stopifnot(file.exists(c2_in))
c2 <- fread(c2_in)
stopifnot(identical(names(c2), c("Segment", "Group3", "n_cells", "segment_total",
                                 "proportion", "percentage")))
segment_order <- c("shared", "arm1", "arm2")
group_order <- c("Ctrl", "Neurological", "Tumor")
stopifnot(setequal(unique(c2$Segment), segment_order),
          setequal(unique(c2$Group3), group_order),
          nrow(c2) == 9L,
          all(abs(c2$proportion - c2$n_cells / c2$segment_total) < 1e-12))
c2[, Segment := factor(Segment, levels = segment_order)]
# Keep the source PDF's bottom-to-top stack order: Tumor, Neurological, Ctrl.
c2[, Group3 := factor(Group3, levels = c("Tumor", "Neurological", "Ctrl"))]
c2[, x := match(as.character(Segment), segment_order)]
c2[, ymin := vapply(seq_len(.N), function(i) {
  prior <- c2[Segment == c2$Segment[i] & as.integer(Group3) < as.integer(c2$Group3[i]),
             proportion]
  sum(prior)
}, numeric(1))]
c2[, ymax := ymin + proportion]
c2[, ymid := (ymin + ymax) / 2]
c2[, pct_label := sprintf("%.1f%%", percentage)]
c2[, label_inside := !(Group3 == "Ctrl" & Segment != "shared")]
c2[, label_colour := fifelse(Group3 == "Ctrl", "#111111", "#FFFFFF")]
# Wider bars reduce the visual gap while preserving the three segment centers.
c2[, `:=`(x0 = x - 0.35, x1 = x + 0.35)]
c2[, segment_axis := sprintf("%s\n%s cells",
  c(shared = "Shared", arm1 = "Arm 1", arm2 = "Arm 2")[as.character(Segment)],
  format(segment_total, big.mark = ","))]
segment_axis_map <- c2[!duplicated(Segment),
                       setNames(segment_axis, as.character(Segment))]

inside <- c2[label_inside == TRUE]
callouts <- c2[Group3 == "Ctrl" & Segment != "shared"]
callouts[, `:=`(
  x_start = ifelse(Segment == "arm1", x0, x1),
  x_end = ifelse(Segment == "arm1", x0 - 0.09, x1 + 0.09),
  x_label = ifelse(Segment == "arm1", x0 - 0.12, x1 + 0.12),
  hjust = ifelse(Segment == "arm1", 1, 0),
  y_label = 0.052
)]
p_c2 <- ggplot() +
  geom_rect(data = c2,
    aes(xmin = x0, xmax = x1, ymin = ymin, ymax = ymax, fill = Group3),
    colour = "white", linewidth = 0.45) +
  geom_text(data = inside,
    aes(x = x, y = ymid, label = pct_label, colour = label_colour),
    family = FONT, size = 2.45, fontface = "plain", show.legend = FALSE) +
  geom_segment(data = callouts,
    aes(x = x_start, xend = x_end, y = ymid, yend = y_label),
    linewidth = 0.35, colour = "#222222", lineend = "butt") +
  geom_text(data = callouts,
    aes(x = x_label, y = y_label, label = pct_label, hjust = hjust),
    family = FONT, size = 2.45, colour = "#111111", show.legend = FALSE) +
  scale_fill_manual(values = GROUP_COLS, breaks = group_order, name = NULL) +
  scale_colour_identity() +
  scale_x_continuous(breaks = 1:3,
    labels = unname(segment_axis_map[segment_order]),
    limits = c(0.65, 3.72), expand = expansion(mult = 0)) +
  scale_y_continuous(labels = label_percent(accuracy = 1),
    breaks = seq(0, 1, 0.25), limits = c(0, 1.015),
    expand = expansion(mult = c(0, 0))) +
  labs(title = "Disease-group composition", x = NULL,
       y = "Share of segment cells") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE,
    override.aes = list(colour = NA))) +
  theme_classic(base_family = FONT, base_size = 8) +
  theme(plot.title = element_text(size = 9.5, face = "bold", colour = "black",
                                  margin = margin(b = 2)),
        axis.text = element_text(size = 7.5, colour = "black"),
        axis.title.y = element_text(size = 8, colour = "black"),
        axis.text.x = element_text(lineheight = 0.96, margin = margin(t = 4)),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", legend.text = element_text(size = 7.5),
        legend.key.size = grid::unit(3.2, "mm"),
        plot.margin = margin(3, 4, 2, 3, unit = "mm"))
fwrite(c2[, .(Segment = as.character(Segment), Group3 = as.character(Group3),
              n_cells, segment_total, proportion, percentage, ymin, ymax,
              label_inside, pct_label)],
       file.path(SRC, "C2_segment_group3_composition_display_source.csv"))
save_panel(p_c2, "C2_segment_group3_composition", 100, 82)

###############################################################################
# B. Monocle2 trajectory: render frozen CDS/tree; update only GB04 Disease.
###############################################################################
traj_csv <- "E:/micro_output/fig3/source_data/Fig2i_monocle2_trajectory_source.csv"
cds_rds <- "E:/micro_output/fig3/data/Fig2i_monocle2_cds_default.rds"
changed_csv <- paste0(
  "E:/micro_output/v3_dims16/objects/",
  "merged_obj_GB04_Disease_patch_20261006_214443_changed_cells.csv"
)
palette_csv <- "D:/data/microglia_new/authority/03_reference_configs/disease_palette_contract.csv"
stopifnot(all(file.exists(c(traj_csv, cds_rds, changed_csv, palette_csv))))
traj <- fread(traj_csv)
changed <- fread(changed_csv)
stopifnot(all(c("CellID", "DDRTree_1", "DDRTree_2", "Cluster", "Disease",
                "SampleID_Global", "DatasetID", "State", "Pseudotime") %in% names(traj)),
          "CellID" %in% names(changed), !anyDuplicated(traj$CellID),
          nrow(traj) == 16800L)
gb04_ids <- intersect(changed$CellID, traj$CellID)
stopifnot(length(gb04_ids) == 61L,
          all(traj[CellID %in% gb04_ids, Disease] == "LGG"),
          uniqueN(traj[CellID %in% gb04_ids, SampleID_Global]) == 1L,
          grepl("GB04", unique(traj[CellID %in% gb04_ids, SampleID_Global])),
          all(traj[CellID %in% gb04_ids, DatasetID] == "LUNG_BRAIN_MET_2022"))
traj[, Disease_before_update := Disease]
traj[CellID %in% gb04_ids, Disease := "GBM"]
stopifnot(sum(traj$Disease == "GBM") == 3612L,
          sum(traj$Disease == "LGG") == 433L,
          all(traj[CellID %in% gb04_ids, Disease] == "GBM"))

cds <- readRDS(cds_rds)
stopifnot(ncol(cds) == 16800L,
          setequal(colnames(cds), traj$CellID),
          !is.null(pData(cds)$State), !is.null(pData(cds)$Pseudotime))
coords_before <- reducedDimS(cds)
tree_before <- minSpanningTree(cds)
tree_signature <- paste(igraph::V(tree_before)$name,
                         paste(igraph::as_edgelist(tree_before), collapse = ","),
                         sep = "|")
pd <- Biobase::pData(cds)
pd_source_idx <- match(rownames(pd), traj$CellID)
coord_source_idx <- match(colnames(cds), traj$CellID)
stopifnot(!anyNA(pd_source_idx), !anyNA(coord_source_idx),
          isTRUE(all.equal(as.numeric(coords_before[1, ]),
                           traj$DDRTree_1[coord_source_idx], tolerance = 1e-12)),
          isTRUE(all.equal(as.numeric(coords_before[2, ]),
                           traj$DDRTree_2[coord_source_idx], tolerance = 1e-12)),
          identical(as.character(pd$State), as.character(traj$State[pd_source_idx])),
          isTRUE(all.equal(as.numeric(pd$Pseudotime),
                           as.numeric(traj$Pseudotime[pd_source_idx]), tolerance = 1e-12)))
pd_disease_before <- as.character(pd$Disease)
stopifnot(all(pd_disease_before[match(gb04_ids, rownames(pd))] == "LGG"))
pd_disease_after <- pd_disease_before
pd_disease_after[match(gb04_ids, rownames(pd))] <- "GBM"
if (is.factor(pd$Disease)) {
  pData(cds)$Disease <- factor(pd_disease_after,
    levels = unique(c(levels(pd$Disease), "GBM")))
} else {
  pData(cds)$Disease <- pd_disease_after
}
stopifnot(identical(as.character(pData(cds)$Disease), pd_disease_after),
          identical(coords_before, reducedDimS(cds)))
tree_after <- minSpanningTree(cds)
tree_signature_after <- paste(igraph::V(tree_after)$name,
  paste(igraph::as_edgelist(tree_after), collapse = ","), sep = "|")
stopifnot(identical(tree_signature, tree_signature_after))

pal <- fread(palette_csv)
stopifnot(all(c("Disease_Abbrev", "Disease_Color_Hex", "Display_Label") %in% names(pal)))
pal <- pal[!duplicated(Disease_Abbrev, fromLast = TRUE)]
disease_cols <- setNames(pal$Disease_Color_Hex, pal$Disease_Abbrev)
disease_labels <- setNames(pal$Display_Label, pal$Disease_Abbrev)
# Preserve the disease palette compatibility mapping already declared in the
# frozen Fig2i generator; the current formal label GBM is absent from the older
# palette CSV, while the generator explicitly mapped GBM to the Tumor color.
disease_cols["GBM"] <- "#B6362D"
disease_labels["GBM"] <- "Glioblastoma (GBM)"
present_disease <- sort(unique(traj$Disease))
missing_disease <- setdiff(present_disease, names(disease_cols))
if (length(missing_disease)) {
  stop("Disease abbreviations missing from frozen palette: ",
       paste(missing_disease, collapse = ", "),
       " | palette keys: ", paste(names(disease_cols), collapse = ", "))
}
present_levels <- pal[Disease_Abbrev %chin% present_disease, Disease_Abbrev]
traj[, Disease := factor(Disease, levels = present_levels)]
graph_coords <- reducedDimK(cds)
edge_src <- as.data.table(igraph::as_data_frame(tree_after, what = "edges"))
vertex_xy <- data.table(
  node = colnames(graph_coords),
  graph_x = as.numeric(graph_coords[1, ]),
  graph_y = as.numeric(graph_coords[2, ])
)
edge_plot <- merge(edge_src,
  vertex_xy[, .(from = node, x = graph_x, y = graph_y)],
  by = "from", all.x = TRUE, sort = FALSE)
edge_plot <- merge(edge_plot,
  vertex_xy[, .(to = node, xend = graph_x, yend = graph_y)],
  by = "to", all.x = TRUE, sort = FALSE)
branch_nodes <- as.character(cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points)
branch_plot <- merge(data.table(node = branch_nodes,
                                branch_point_idx = seq_along(branch_nodes)),
                     vertex_xy, by = "node", all.x = TRUE, sort = FALSE)
stopifnot(nrow(edge_plot) == igraph::ecount(tree_after),
          all(is.finite(edge_plot$x)), all(is.finite(edge_plot$y)),
          all(is.finite(edge_plot$xend)), all(is.finite(edge_plot$yend)),
          nrow(branch_plot) == length(branch_nodes),
          all(is.finite(branch_plot$graph_x)), all(is.finite(branch_plot$graph_y)))
set.seed(42)
traj_plot <- traj[sample(.N)]
p_traj <- ggplot() +
  geom_segment(data = edge_plot,
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "#4A4A4A", linewidth = 0.38, lineend = "round") +
  geom_point(data = traj_plot,
    aes(x = DDRTree_1, y = DDRTree_2, colour = Disease),
    size = 0.25, alpha = 0.62, stroke = 0) +
  geom_point(data = branch_plot,
    aes(x = graph_x, y = graph_y),
    inherit.aes = FALSE, size = 2.4, colour = "#222222") +
  geom_text(data = branch_plot,
    aes(x = graph_x, y = graph_y, label = branch_point_idx),
    inherit.aes = FALSE, size = 1.8, colour = "white", family = FONT) +
  scale_color_manual(values = disease_cols, breaks = present_levels,
                     labels = disease_labels[present_levels], drop = FALSE,
                     name = "Disease") +
  labs(title = "Disease distribution along the trajectory",
       subtitle = NULL, x = "Component 1", y = "Component 2") +
  guides(colour = guide_legend(ncol = 3, byrow = TRUE,
    override.aes = list(size = 1.45, alpha = 1))) +
  theme_classic(base_family = FONT, base_size = 8) +
  theme(plot.title = element_text(size = 9.5, face = "bold", colour = "black",
                                  margin = margin(b = 3)),
        axis.title = element_text(size = 8, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"),
        legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(size = 6.5, colour = "black"),
        legend.title = element_text(size = 7, colour = "black"),
        legend.key.size = grid::unit(2.8, "mm"),
        legend.spacing.x = grid::unit(1.2, "mm"),
        plot.margin = margin(3, 4, 2, 3, unit = "mm"))
traj_width_mm <- 145
coord_range_x <- diff(range(traj$DDRTree_1))
coord_range_y <- diff(range(traj$DDRTree_2))
# Preserve the original trajectory panel's physical aspect while moving its
# legend below the plot; this does not rotate or alter any saved coordinates.
source_panel_aspect <- 0.73
source_coord_ratio <- source_panel_aspect / (coord_range_y / coord_range_x)
p_traj <- p_traj + coord_fixed(ratio = source_coord_ratio)
traj_height_mm <- traj_width_mm * source_panel_aspect + 34
fwrite(traj[, .(CellID, DDRTree_1, DDRTree_2, Cluster, Disease,
                SampleID_Global, DatasetID, State, Pseudotime)],
       file.path(SRC, "Fig2i_monocle2_trajectory_GB04_updated_source.csv"))
fwrite(data.table(CellID = traj_plot$CellID,
                  plot_order = seq_len(nrow(traj_plot)), seed = 42L),
       file.path(SRC, "Fig2i_disease_trajectory_plot_order_seed42.csv"))
fwrite(data.table(
  change = "Disease-only synchronization from authorized current Seurat metadata",
  matched_trajectory_cells = length(gb04_ids),
  before = "LGG", after = "GBM",
  sample = paste(unique(traj[CellID %in% gb04_ids, SampleID_Global]), collapse = ";"),
  dataset = "LUNG_BRAIN_MET_2022",
  GBM_palette_mapping = "#B6362D, same compatibility mapping as frozen Fig2i generator",
  drawing_order = "Fixed pseudorandom order, seed 42; all 16800 cells retained",
  GBM_before = 3551L, GBM_after = 3612L,
  LGG_before = 494L, LGG_after = 433L,
  coordinates_tree_state_pseudotime_changed = FALSE
), file.path(SRC, "Fig2i_GB04_disease_display_update_audit.csv"))
save_panel(p_traj, "Fig2i_monocle2_trajectory_by_disease",
           traj_width_mm, traj_height_mm)

###############################################################################
# D. Frozen S4h matrix and row order: layout-only heatmap rendition.
###############################################################################
hm_file <- "E:/micro_output/fig3/trajectory_S4/source_data/default/S4h_no_labels_display_source.csv"
ord_file <- "E:/micro_output/fig3/trajectory_S4/source_data/default/S4h_no_labels_row_order.csv"
grid_file <- "E:/micro_output/fig3/trajectory_S4/data/default/s4h_time_grid.csv"
module_file <- "E:/micro_output/fig3/trajectory_S4/tables/default/s4h_dynamic_genes_modules.csv"
tree_file <- "E:/micro_output/fig3/trajectory_S4/data/default/s4h_row_dendrogram_hclust.rds"
stopifnot(all(file.exists(c(hm_file, ord_file, grid_file, module_file, tree_file))))
hm <- fread(hm_file, header = TRUE)
ord <- fread(ord_file)
grid <- fread(grid_file)
modules <- fread(module_file)
hc <- readRDS(tree_file)
stopifnot(nrow(hm) == 4843L, nrow(ord) == 4843L,
          ncol(hm) == 101L, nrow(grid) == 100L,
          all(c("gene_id", "PTG", "row_rank") %in% names(ord)),
          identical(hm[[1]], ord$gene_id),
          identical(ord$row_rank, seq_len(nrow(ord))),
          length(hc$labels) == 4843L,
          setequal(hc$labels, ord$gene_id),
          setequal(unique(ord$PTG), names(PTG_COLS)))
tree_leaf_order_matches_display <- identical(
  unname(hc$labels[hc$order]), as.character(ord$gene_id)
)
grid_cols <- names(hm)[-1]
hm_mat <- as.matrix(hm[, ..grid_cols])
storage.mode(hm_mat) <- "double"
stopifnot(all(is.finite(hm_mat)), min(hm_mat) >= -3 - 1e-12,
          max(hm_mat) <= 3 + 1e-12)
expected_counts <- c(PTG01 = 1614L, PTG02 = 1377L, PTG03 = 934L, PTG04 = 918L)
actual_counts <- table(factor(ord$PTG, levels = names(expected_counts)))
stopifnot(identical(as.integer(actual_counts), unname(expected_counts)),
          identical(as.integer(modules[, .N, by = PTG][
            match(names(expected_counts), PTG), N]), unname(expected_counts)))
ptg_runs <- rle(as.character(ord$PTG))
stopifnot(length(ptg_runs$values) == 4L,
          identical(ptg_runs$values, names(expected_counts)))
long <- data.table(
  gene_id = rep(ord$gene_id, each = ncol(hm_mat)),
  PTG = rep(ord$PTG, each = ncol(hm_mat)),
  row_rank = rep(ord$row_rank, each = ncol(hm_mat)),
  grid_point = rep(seq_len(ncol(hm_mat)), times = nrow(hm_mat)),
  pseudotime = rep(grid$pseudotime, times = nrow(hm_mat)),
  row_z = as.vector(t(hm_mat))
)
long[, x := grid_point]
module_spans <- ord[, .(first = min(row_rank), last = max(row_rank), n = .N), by = PTG]
module_spans[, `:=`(y_mid = (first + last) / 2,
                    y_min = first - 0.5, y_max = last + 0.5)]
module_spans[, strip_xmin := -0.75]
module_spans[, strip_xmax := -0.45]
module_spans[, label_x := -1.0]
module_spans[, module_label := sprintf("%s (n = %s)", PTG,
                                       format(n, big.mark = ","))]
grid_breaks <- c(1L, 25L, 50L, 75L, 100L)
grid_labels <- format(signif(grid$pseudotime[grid_breaks], 3), trim = TRUE)
separator_y <- module_spans$last[-nrow(module_spans)] + 0.5
p_hm <- ggplot(long, aes(x = x, y = row_rank, fill = row_z)) +
  geom_raster(interpolate = FALSE) +
  geom_rect(data = module_spans,
    aes(xmin = strip_xmin, xmax = strip_xmax, ymin = y_min, ymax = y_max),
    inherit.aes = FALSE, fill = unname(PTG_COLS[module_spans$PTG]),
    colour = NA) +
  geom_text(data = module_spans,
    aes(x = label_x, y = y_mid, label = module_label),
    inherit.aes = FALSE, hjust = 1, family = FONT, size = 2.0,
    colour = "#111111") +
  geom_hline(yintercept = separator_y, colour = "white", linewidth = 0.38) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, limits = c(-3, 3), breaks = c(-3, 0, 3),
    oob = scales::squish, name = "Row Z-score",
    guide = guide_colorbar(barheight = grid::unit(34, "mm"),
                           barwidth = grid::unit(2.8, "mm"),
                           title.position = "top")) +
  scale_x_continuous(breaks = grid_breaks, labels = grid_labels,
    limits = c(-28, 100.5), expand = expansion(mult = c(0, 0))) +
  scale_y_reverse(limits = c(nrow(ord) + 0.5, 0.5), breaks = NULL,
    expand = expansion(mult = c(0, 0))) +
  labs(title = "Dynamic gene expression across pseudotime",
       x = "Pseudotime", y = NULL) +
  theme_minimal(base_family = FONT, base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 7, colour = "black"),
        axis.title.x = element_text(size = 8, colour = "black", margin = margin(t = 3)),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.title = element_text(size = 9.5, face = "bold", colour = "black",
                                  margin = margin(b = 3)),
        legend.position = "right", legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 7),
        legend.margin = margin(0, 0, 0, 1, unit = "mm"),
        plot.margin = margin(3, 2, 2, 2, unit = "mm"))
display_source <- copy(ord)
display_source <- cbind(display_source, as.data.table(hm_mat))
setnames(display_source, names(display_source)[-(1:3)],
         sprintf("grid_%03d", seq_len(ncol(hm_mat))))
fwrite(display_source, file.path(SRC, "S4h_dynamic_gene_heatmap_display_source.csv"))
fwrite(grid, file.path(SRC, "S4h_pseudotime_grid_source.csv"))
save_panel(p_hm, "S4h_dynamic_gene_heatmap_no_labels", 112, 132)

###############################################################################
# Local provenance and quantitative invariance report.
###############################################################################
invariance <- data.table(
  figure = c("C2_segment_group3_composition", "Fig2i_monocle2_trajectory_by_disease",
             "Fig3A_clean24", "S4h_dynamic_gene_heatmap_no_labels"),
  source_verified = c(TRUE, TRUE, FALSE, TRUE),
  output_rendered = c(TRUE, TRUE, FALSE, TRUE),
  invariant_check = c(
    "Original segment x disease counts, denominators and proportions retained",
    "All 16800 CellIDs, DDRTree coordinates, tree edges, State and Pseudotime retained; 61 GB04 Disease labels updated LGG to GBM",
    "Not rendered: frozen 24-cluster hclust/dendrogram with merge heights is absent; no reclustering or PDF reconstruction permitted",
    paste0("4843 genes, 100 grid values, explicit display row order, PTG assignments and [-3,3] display range retained; serialized hclust gene set matches but leaf order ",
           if (tree_leaf_order_matches_display) "matches" else "DIFFERS from display order; no dendrogram redrawn")
  )
)
fwrite(invariance, file.path(LOG, "beautification_invariance_audit.csv"))
cat("Rendered panels A, B, and D. Fig3A remains pending missing frozen hclust source.\n")
