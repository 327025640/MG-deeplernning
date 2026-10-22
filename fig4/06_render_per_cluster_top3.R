#!/usr/bin/env Rscript

project_lib <- "D:/data/microglia_new/Rlibs/Fig4C_R4.4.3"
.libPaths(c(project_lib, .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 06_render_per_cluster_top3.R <staging_input_dir> <output_dir>")
input_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

score_path <- file.path(input_dir, "TF_activity_all_raw_long.csv")
coverage_path <- file.path(input_dir, "TF_target_coverage_and_ULM_status.csv")
order_path <- file.path(input_dir, "cluster_display_order.csv")
geometry_path <- file.path(input_dir, "anchor_row_center_geometry_used.csv")
scope_path <- file.path(input_dir, "cluster_scope.csv")
required <- c(score_path, coverage_path, order_path, geometry_path, scope_path)
if (!all(file.exists(required))) stop("A full ULM score, coverage, order, geometry, or scope table is missing.")

all_scores <- fread(score_path)
coverage <- fread(coverage_path)
order_table <- fread(order_path)
geometry <- fread(geometry_path)
scope_table <- fread(scope_path)
expected_order <- as.character(order_table[order(row_top_to_bottom), Cluster])
if (length(expected_order) != 24L || anyDuplicated(expected_order)) stop("Frozen cluster display order is invalid.")
if (!all(c("source", "condition", "score", "eligible_min10", "score_status") %in% names(all_scores)) ||
    !all(c("source", "effective_target_n", "eligible_min10", "all_24_finite") %in% names(coverage))) {
  stop("Full score or coverage table schema is incomplete.")
}
all_scores$source <- as.character(all_scores$source)
all_scores$condition <- as.character(all_scores$condition)
coverage$source <- as.character(coverage$source)
if (anyDuplicated(all_scores[, .(source, condition)]) > 0L ||
    !setequal(unique(all_scores$condition), expected_order)) {
  stop("Full ULM score table contains duplicate TF-cluster rows or an unexpected cluster set.")
}
eligible_sources <- coverage[eligible_min10 == TRUE & all_24_finite == TRUE, source]
if (length(eligible_sources) != 570L) stop("Eligible TF set differs from the audited 570 complete regulators.")
eligible_scores <- all_scores[source %in% eligible_sources & is.finite(score) & score_status == "scored"]
eligible_scores[, cluster_display_order := match(condition, expected_order)]
if (uniqueN(eligible_scores$source) != length(eligible_sources) ||
    uniqueN(eligible_scores$condition) != 24L ||
    uniqueN(eligible_scores[, .(source, condition)]) != 570L * 24L) {
  stop("Eligible TF score matrix is not a complete 570 x 24 finite matrix.")
}
cluster_score_counts <- eligible_scores[, .(n_tf = uniqueN(source)), by = condition]
if (nrow(cluster_score_counts) != 24L || any(cluster_score_counts$n_tf != 570L)) {
  stop("At least one cluster lacks one or more eligible TF scores.")
}

setorder(eligible_scores, cluster_display_order, -score, source)
per_cluster_top3 <- eligible_scores[, head(.SD, 3L), by = .(condition, cluster_display_order)]
per_cluster_top3[, within_cluster_rank := seq_len(.N), by = condition]
setnames(per_cluster_top3, c("condition", "source", "score"),
         c("Cluster", "TF", "raw_score"))
per_cluster_top3 <- merge(
  per_cluster_top3,
  coverage[, .(TF = source, effective_target_n, prior_unique_target_n,
               effective_target_coverage)],
  by = "TF", all.x = TRUE, sort = FALSE
)
setorder(per_cluster_top3, cluster_display_order, within_cluster_rank)
selection_counts <- per_cluster_top3[, .N, by = Cluster]
if (nrow(per_cluster_top3) != 72L || nrow(selection_counts) != 24L ||
    any(selection_counts$N != 3L) ||
    anyDuplicated(per_cluster_top3[, .(Cluster, TF)]) > 0L ||
    anyNA(per_cluster_top3$effective_target_n)) {
  stop("Per-cluster Top3 selection table is not exactly 72 valid unique rows.")
}

tf_order <- per_cluster_top3[!duplicated(TF),
  .(TF, first_cluster = Cluster, first_cluster_rank = within_cluster_rank)]
union_tfs <- tf_order$TF
if (length(union_tfs) < 3L || anyDuplicated(union_tfs)) stop("Per-cluster Top3 TF union is invalid.")

full_union_scores <- eligible_scores[source %in% union_tfs]
full_union_scores$Cluster <- full_union_scores$condition
full_union_scores$TF <- full_union_scores$source
full_union_scores[, Cluster := factor(Cluster, levels = expected_order)]
full_union_scores[, TF := factor(TF, levels = union_tfs)]
if (nrow(full_union_scores) != 24L * length(union_tfs) ||
    uniqueN(full_union_scores[, .(condition, source)]) != 24L * length(union_tfs) ||
    any(!is.finite(full_union_scores$score))) {
  stop("Union TFs do not have finite scores across all 24 clusters.")
}
full_union_scores[, raw_score_mean := mean(score), by = source]
full_union_scores[, raw_score_sd := sd(score), by = source]
if (any(!is.finite(full_union_scores$raw_score_sd)) || any(full_union_scores$raw_score_sd <= 0)) {
  stop("A union TF has no cross-cluster score variation for Z standardization.")
}
full_union_scores[, z_score := (score - raw_score_mean) / raw_score_sd]
full_union_scores$Cluster <- as.character(full_union_scores$Cluster)
full_union_scores$TF <- as.character(full_union_scores$TF)
full_union_scores[, cluster_display_order := match(Cluster, expected_order)]
full_union_scores[, TF_display_order := match(TF, union_tfs)]
selection_key <- per_cluster_top3[, .(Cluster, TF, is_cluster_top3 = TRUE)]
full_union_scores <- merge(
  full_union_scores, selection_key, by = c("Cluster", "TF"), all.x = TRUE, sort = FALSE
)
full_union_scores[is.na(is_cluster_top3), is_cluster_top3 := FALSE]
if (sum(full_union_scores$is_cluster_top3) != 72L) stop("Top3 selection markers do not identify exactly 72 matrix cells.")

geometry <- geometry[order(row_top_to_bottom)]
if (nrow(geometry) != 24L || !identical(as.character(geometry$Cluster), expected_order) ||
    length(unique(geometry$anchor_page_height_mm)) != 1L || geometry$anchor_page_height_mm[1L] != 490 ||
    any(!is.finite(geometry$center_mm_from_top)) || any(!is.finite(geometry$tile_height_mm))) {
  stop("FUN8 anchor geometry is incomplete or differs from the frozen cluster order.")
}
scope_map <- unique(scope_table[, .(Cluster, analysis_scope)])
scope_map$Cluster <- as.character(scope_map$Cluster)
scope_map$analysis_scope <- as.character(scope_map$analysis_scope)
scope_colors <- c(Shared6 = "#44718C", Restricted4 = "#B87555", pan14 = "#CFCFCF")
order_scope <- order_table[order(row_top_to_bottom), .(Cluster, scope_group = as.character(scope_group))]
if (nrow(scope_map) != 24L || anyNA(scope_map$analysis_scope) ||
    any(!scope_map$analysis_scope %in% names(scope_colors))) stop("Frozen cluster group mapping is invalid.")
scope_check <- merge(order_scope, scope_map, by = "Cluster", all = TRUE, sort = FALSE)
if (nrow(scope_check) != 24L || anyNA(scope_check$analysis_scope) ||
    any(scope_check$scope_group != scope_check$analysis_scope)) stop("Order group labels disagree with frozen cluster scope.")

plot_height_mm <- 490
n_tfs <- length(union_tfs)
heatmap_left_mm <- 45
heatmap_right_mm <- 8
row_tile_heights_mm <- setNames(geometry$tile_height_mm, geometry$Cluster)
# Preserve the FUN8 row geometry and make TF tiles square to within the small
# row-specific height variation already present in the frozen anchor.
tf_cell_width_mm <- mean(row_tile_heights_mm)
heatmap_width_mm <- n_tfs * tf_cell_width_mm
plot_width_mm <- heatmap_left_mm + heatmap_width_mm + heatmap_right_mm
tf_cell_width_mm <- heatmap_width_mm / n_tfs
tf_x_centres <- heatmap_left_mm + (seq_len(n_tfs) - 0.5) * tf_cell_width_mm
row_centres_from_top <- setNames(geometry$center_mm_from_top, geometry$Cluster)
row_centres_from_bottom <- plot_height_mm - row_centres_from_top
z_limit <- max(abs(full_union_scores$z_score))
if (!is.finite(z_limit) || z_limit <= 0) stop("The per-cluster Top3 union Z-score scale is invalid.")

map_z_colour <- function(z) {
  result <- rep("#BDBDBD", length(z))
  valid <- is.finite(z)
  if (any(valid)) {
    position <- pmin(1, pmax(0, (z[valid] / z_limit + 1) / 2))
    rgb_values <- grDevices::colorRamp(c("#2166AC", "#F7F7F7", "#B2182B"), space = "Lab")(position)
    result[valid] <- grDevices::rgb(rgb_values[, 1], rgb_values[, 2], rgb_values[, 3], maxColorValue = 255)
  }
  result
}
full_union_scores[, draw_colour := map_z_colour(z_score)]

layout_audit <- data.table(
  row_top_to_bottom = seq_along(expected_order), Cluster = expected_order,
  anchor_center_mm_from_top = as.numeric(row_centres_from_top[expected_order]),
  rendered_center_mm_from_top = as.numeric(row_centres_from_top[expected_order]),
  center_delta_mm = 0,
  tile_height_mm = as.numeric(row_tile_heights_mm[expected_order]),
  group = order_scope$scope_group
)

draw_panel <- function() {
  grid.newpage()
  pushViewport(viewport(xscale = c(0, plot_width_mm), yscale = c(0, plot_height_mm),
                        width = unit(1, "npc"), height = unit(1, "npc"), clip = "off"))
  grid.rect(gp = gpar(fill = "white", col = NA))
  grid.text("C. Per-cluster top-three inferred regulator activity",
            x = unit(9, "mm"), y = unit(486, "mm"), just = c("left", "centre"),
            gp = gpar(fontfamily = "Arial", fontsize = 10, fontface = "bold", col = "black"))
  grid.text("Columns are the union of raw-score Top3; full scores shown across all clusters",
            x = unit(9, "mm"), y = unit(478, "mm"), just = c("left", "centre"),
            gp = gpar(fontfamily = "Arial", fontsize = 7.2, col = "#333333"))

  for (i in seq_len(nrow(geometry))) {
    cl <- geometry$Cluster[i]
    y <- row_centres_from_bottom[[cl]]
    group <- scope_map$analysis_scope[match(cl, scope_map$Cluster)]
    grid.rect(x = unit(heatmap_left_mm - 5, "mm"), y = unit(y, "mm"),
              width = unit(3.7, "mm"), height = unit(row_tile_heights_mm[[cl]], "mm"),
              gp = gpar(fill = unname(scope_colors[group]), col = NA))
    grid.text(cl, x = unit(heatmap_left_mm - 8, "mm"), y = unit(y, "mm"),
              just = c("right", "centre"), gp = gpar(fontfamily = "Arial", fontsize = 7.2, col = "black"))
  }
  for (i in seq_len(nrow(full_union_scores))) {
    row <- full_union_scores[i]
    x <- tf_x_centres[row$TF_display_order]
    y <- row_centres_from_bottom[[as.character(row$Cluster)]]
    grid.rect(x = unit(x, "mm"), y = unit(y, "mm"),
              width = unit(tf_cell_width_mm, "mm"),
              height = unit(row_tile_heights_mm[[as.character(row$Cluster)]], "mm"),
              gp = gpar(fill = row$draw_colour, col = "white", lwd = 0.5))
  }
  for (boundary in c(6L, 10L)) {
    top_y <- mean(row_centres_from_top[expected_order[c(boundary, boundary + 1L)]])
    y <- plot_height_mm - top_y
    grid.lines(x = unit(c(heatmap_left_mm, plot_width_mm - heatmap_right_mm), "mm"),
               y = unit(c(y, y), "mm"), gp = gpar(col = "#555555", lwd = 0.6))
  }
  label_y <- plot_height_mm - max(row_centres_from_top) - max(row_tile_heights_mm) / 2 - 7
  label_cex <- max(5.2, min(7, 400 / n_tfs))
  for (i in seq_along(union_tfs)) {
    grid.text(union_tfs[i], x = unit(tf_x_centres[i], "mm"), y = unit(label_y, "mm"),
              rot = 45, just = c("right", "centre"),
              gp = gpar(fontfamily = "Arial", fontsize = label_cex, col = "black"))
  }

  cb_width <- min(64, plot_width_mm * 0.25)
  cb_left <- plot_width_mm / 2 - cb_width / 2
  cb_bottom <- 15.5
  cb_height <- 2.5
  n_steps <- 101L
  ramp <- map_z_colour(seq(-z_limit, z_limit, length.out = n_steps))
  for (i in seq_len(n_steps)) {
    grid.rect(x = unit(cb_left + (i - 0.5) * cb_width / n_steps, "mm"),
              y = unit(cb_bottom + cb_height / 2, "mm"),
              width = unit(cb_width / n_steps + 0.02, "mm"), height = unit(cb_height, "mm"),
              gp = gpar(fill = ramp[i], col = NA))
  }
  grid.text("Relative inferred activity (Z-score)", x = unit(plot_width_mm / 2, "mm"), y = unit(21.5, "mm"),
            gp = gpar(fontfamily = "Arial", fontsize = 7.2, col = "black"))
  breaks <- c(-z_limit, 0, z_limit)
  break_x <- cb_left + (breaks + z_limit) / (2 * z_limit) * cb_width
  break_labels <- c(formatC(-z_limit, digits = 2, format = "f"), "0", formatC(z_limit, digits = 2, format = "f"))
  for (i in seq_along(breaks)) {
    grid.text(break_labels[i], x = unit(break_x[i], "mm"), y = unit(12.3, "mm"),
              gp = gpar(fontfamily = "Arial", fontsize = 6.6, col = "black"))
  }
  group_order <- c("Shared6", "Restricted4", "pan14")
  legend_spacing <- 36
  legend_start <- plot_width_mm / 2 - legend_spacing
  for (i in seq_along(group_order)) {
    group <- group_order[i]
    x <- legend_start + (i - 1) * legend_spacing
    grid.rect(x = unit(x, "mm"), y = unit(6.2, "mm"), width = unit(4, "mm"), height = unit(3, "mm"),
              gp = gpar(fill = unname(scope_colors[group]), col = NA))
    grid.text(group, x = unit(x + 3, "mm"), y = unit(6.2, "mm"), just = c("left", "centre"),
              gp = gpar(fontfamily = "Arial", fontsize = 6.8, col = "black"))
  }
  grid.text("Within-regulator Z-scores are relative patterns, not absolute activity.",
            x = unit(9, "mm"), y = unit(2.2, "mm"), just = c("left", "centre"),
            gp = gpar(fontfamily = "Arial", fontsize = 6.2, col = "#333333"))
  popViewport()
}

png_path <- file.path(output_dir, "Fig4C_cluster24_TF_activity_per_cluster_top3.png")
pdf_path <- file.path(output_dir, "Fig4C_cluster24_TF_activity_per_cluster_top3.pdf")
png(filename = png_path, width = round(plot_width_mm / 25.4 * 600),
    height = round(plot_height_mm / 25.4 * 600), units = "px", res = 600,
    bg = "white", family = "Arial", type = "cairo")
draw_panel()
dev.off()
cairo_pdf(filename = pdf_path, width = plot_width_mm / 25.4, height = plot_height_mm / 25.4,
          family = "Arial", bg = "white")
draw_panel()
dev.off()

selection_export <- per_cluster_top3[, .(Cluster, TF, within_cluster_rank, raw_score, effective_target_n,
                                         prior_unique_target_n, effective_target_coverage)]
fwrite(selection_export, file.path(output_dir, "TF_per_cluster_top3_selection.csv"))
source_export <- full_union_scores[, .(
  Cluster, TF, cluster_display_order, TF_display_order, score, raw_score_mean, raw_score_sd,
  z_score, is_cluster_top3
)]
source_export <- merge(
  source_export,
  coverage[, .(TF = source, effective_target_n)],
  by = "TF", all.x = TRUE, sort = FALSE
)
fwrite(source_export[order(cluster_display_order, TF_display_order)],
       file.path(output_dir, "Fig4C_cluster_by_TF_source_per_cluster_top3.csv"))
fwrite(tf_order, file.path(output_dir, "TF_per_cluster_top3_union_column_order.csv"))
fwrite(layout_audit, file.path(output_dir, "Fig4C_per_cluster_top3_anchor_row_center_alignment.csv"))
fwrite(data.table(
  item = c("selection_rule", "eligible_TF_n", "cluster_n", "selected_rows", "union_TF_n",
           "union_column_order", "color_scale_symmetric_z_limit", "figure_width_mm", "figure_height_mm",
           "tile_width_mm", "tile_height_range_mm", "A_anchor_alignment", "interpretation"),
  value = c("within each cluster, descending raw signed ULM score; ties TF ascending",
            as.character(length(eligible_sources)), "24", as.character(nrow(per_cluster_top3)),
            as.character(n_tfs), paste(union_tfs, collapse = ","),
            as.character(z_limit), as.character(plot_width_mm), as.character(plot_height_mm),
            as.character(tf_cell_width_mm), paste(range(geometry$tile_height_mm), collapse = " to "),
            "same 24 row centers/order/tile heights; groups after rows 6 and 10",
            "descriptive relative inferred activity; not significance, donor validation, or causality")
), file.path(output_dir, "per_cluster_top3_render_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "per_cluster_top3_sessionInfo.txt"))
message("Rendered per-cluster Top3 union; ULM scores and expression aggregation were not rerun.")
