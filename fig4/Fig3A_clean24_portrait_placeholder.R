#!/usr/bin/env Rscript

# Re-layout the frozen Fig3A clean24 dendrogram to match the portrait Fig4
# placeholder. The linkage is reconstructed only from the saved distance
# matrix and recorded Ward.D2 method; no expression profiles or distances are
# recalculated. See the input and output paths below.

suppressPackageStartupMessages({
  library(data.table)
  library(grid)
})

distance_path <- "E:/micro_output/fig3/fig4/source_data/Fig3A_clean24_distance_matrix.csv"
leaf_order_path <- "E:/micro_output/fig3/fig4/data/Fig3A_clean24_leaf_order.csv"
output_dir <- "D:/data/microglia_new/_staging_Fig3A_portrait_placeholder_v3"

if (!file.exists(distance_path) || !file.exists(leaf_order_path)) {
  stop("Frozen distance matrix or leaf-order source is missing.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

distance_dt <- fread(distance_path)
if (ncol(distance_dt) != 25L || names(distance_dt)[1L] != "Cluster") {
  stop("Distance matrix must contain one Cluster column and 24 cluster columns.")
}
distance_labels <- as.character(distance_dt$Cluster)
distance_matrix <- as.matrix(distance_dt[, -1L, with = FALSE])
storage.mode(distance_matrix) <- "double"
rownames(distance_matrix) <- distance_labels
if (!identical(distance_labels, names(distance_dt)[-1L]) ||
    !identical(rownames(distance_matrix), colnames(distance_matrix)) ||
    anyDuplicated(distance_labels) || any(!is.finite(distance_matrix)) ||
    any(abs(diag(distance_matrix)) > 1e-12) ||
    max(abs(distance_matrix - t(distance_matrix))) > 1e-10) {
  stop("Frozen distance matrix labels, values, diagonal, or symmetry are invalid.")
}

leaf_table <- fread(leaf_order_path)
if (!all(c("leaf_order", "Cluster", "FinalNameWithStatus") %in% names(leaf_table)) ||
    nrow(leaf_table) != 24L || anyDuplicated(leaf_table$Cluster) ||
    any(!is.finite(leaf_table$leaf_order)) ||
    any(leaf_table$leaf_order != as.integer(leaf_table$leaf_order)) ||
    anyDuplicated(as.integer(leaf_table$leaf_order)) > 0L ||
    !identical(sort(as.integer(leaf_table$leaf_order)), seq_len(24L))) {
  stop("Frozen leaf-order source must contain 24 unique cluster labels and names.")
}
setorder(leaf_table, leaf_order)
frozen_order_bottom_to_top <- as.character(leaf_table$Cluster)
if (!setequal(frozen_order_bottom_to_top, distance_labels) ||
    anyNA(leaf_table$FinalNameWithStatus) || any(!nzchar(leaf_table$FinalNameWithStatus))) {
  stop("Frozen leaf order and distance matrix cluster identities do not match.")
}

# Reconstruct the linkage object from the frozen distance matrix with the
# recorded method, then require exact agreement with the frozen leaf order.
hc <- stats::hclust(stats::as.dist(distance_matrix), method = "ward.D2")
if (!identical(as.character(hc$labels[hc$order]), frozen_order_bottom_to_top) ||
    nrow(hc$merge) != 23L || any(!is.finite(hc$height)) ||
    any(diff(hc$height) < -1e-12)) {
  stop("Ward.D2 reconstruction does not reproduce the frozen leaf order or valid merge heights.")
}
saveRDS(hc, file.path(output_dir, "Fig3A_clean24_hclust_reconstructed_from_frozen_distance_matrix.rds"))

# Pixel geometry measured from fig-4.jpg at 96 dpi: black placeholder x=22..162,
# y=598..860 (inclusive). This fixes the panel page dimensions and aspect ratio.
page_width_mm <- 141 * 25.4 / 96
page_height_mm <- 263 * 25.4 / 96
page_ratio <- page_width_mm / page_height_mm
if (abs(page_ratio - 141 / 263) > 1e-12) stop("Placeholder page aspect ratio changed unexpectedly.")

leaf_table[, label := paste(Cluster, "|", FinalNameWithStatus)]
leaf_table[, display_order_top_to_bottom := 25L - as.integer(leaf_order)]
setorder(leaf_table, display_order_top_to_bottom)

tree_x_left_mm <- 2.0
tree_x_right_mm <- 15.1
labels_x_mm <- 15.6
tree_y_top_mm <- page_height_mm - 7.0
tree_y_bottom_mm <- 10.2
axis_y_mm <- 6.1
tick_label_y_mm <- 4.5
axis_title_y_mm <- 1.8
title_y_mm <- page_height_mm - 2.1
label_font_pt <- 4.2
title_font_pt <- 5.7
title_text <- "Transcriptional similarity"
line_width_pt <- 0.38

leaf_table[, y_mm := seq(tree_y_top_mm, tree_y_bottom_mm, length.out = .N)]
leaf_y <- setNames(leaf_table$y_mm, leaf_table$Cluster)
leaf_name <- setNames(leaf_table$label, leaf_table$Cluster)
max_height <- max(hc$height)
if (!is.finite(max_height) || max_height <= 0) stop("Ward.D2 merge height range is invalid.")
x_from_height <- function(h) {
  tree_x_left_mm + (max_height - h) / max_height * (tree_x_right_mm - tree_x_left_mm)
}

segments <- list()
merge_table <- data.table(
  merge_id = seq_len(nrow(hc$merge)),
  child_left = hc$merge[, 1L], child_right = hc$merge[, 2L],
  merge_height = hc$height
)
node_positions <- vector("list", nrow(hc$merge))
resolve_node <- function(ref) {
  if (ref < 0L) {
    cluster <- hc$labels[-ref]
    return(list(y = unname(leaf_y[[cluster]]), height = 0,
                clusters = cluster))
  }
  i <- as.integer(ref)
  if (!is.null(node_positions[[i]])) return(node_positions[[i]])
  left <- resolve_node(hc$merge[i, 1L])
  right <- resolve_node(hc$merge[i, 2L])
  parent_height <- hc$height[i]
  parent_x <- x_from_height(parent_height)
  left_x <- x_from_height(left$height)
  right_x <- x_from_height(right$height)
  segments[[length(segments) + 1L]] <<- data.table(
    merge_id = i, segment_type = "child_left_horizontal",
    x0_mm = parent_x, y0_mm = left$y, x1_mm = left_x, y1_mm = left$y,
    height_from = parent_height, height_to = left$height
  )
  segments[[length(segments) + 1L]] <<- data.table(
    merge_id = i, segment_type = "child_right_horizontal",
    x0_mm = parent_x, y0_mm = right$y, x1_mm = right_x, y1_mm = right$y,
    height_from = parent_height, height_to = right$height
  )
  segments[[length(segments) + 1L]] <<- data.table(
    merge_id = i, segment_type = "merge_vertical",
    x0_mm = parent_x, y0_mm = min(left$y, right$y),
    x1_mm = parent_x, y1_mm = max(left$y, right$y),
    height_from = parent_height, height_to = parent_height
  )
  node <- list(y = mean(c(left$y, right$y)), height = parent_height,
               clusters = c(left$clusters, right$clusters))
  node_positions[[i]] <<- node
  node
}
invisible(resolve_node(nrow(hc$merge)))
branch_source <- rbindlist(segments)
if (nrow(branch_source) != 69L || any(!is.finite(as.matrix(branch_source[, .(x0_mm, y0_mm, x1_mm, y1_mm)])))) {
  stop("Reconstructed dendrogram branch geometry is incomplete or non-finite.")
}

png_path <- file.path(output_dir, "Fig3A_clean24_portrait_placeholder.png")
pdf_path <- file.path(output_dir, "Fig3A_clean24_portrait_placeholder.pdf")
draw_panel <- function() {
  grid.newpage()
  pushViewport(viewport(xscale = c(0, page_width_mm), yscale = c(0, page_height_mm),
                        width = unit(1, "npc"), height = unit(1, "npc"), clip = "on"))
  grid.rect(gp = gpar(fill = "white", col = NA))
  title_gp <- gpar(fontfamily = "Arial", fontsize = title_font_pt,
                   fontface = "bold", col = "black")
  title_width_mm <- convertWidth(grobWidth(textGrob(title_text, gp = title_gp)),
                                 unitTo = "mm", valueOnly = TRUE)
  label_width_mm <- convertWidth(
    grobWidth(textGrob(leaf_table$label[which.max(nchar(leaf_table$label))],
                       gp = gpar(fontfamily = "Arial", fontsize = label_font_pt))),
    unitTo = "mm", valueOnly = TRUE
  )
  if (1.1 + title_width_mm > page_width_mm - 0.5 ||
      labels_x_mm + label_width_mm > page_width_mm - 0.5) {
    stop("Title or cluster labels exceed the measured portrait placeholder width.")
  }
  grid.text(title_text,
            x = unit(1.1, "mm"), y = unit(title_y_mm, "mm"),
            just = c("left", "centre"),
            gp = title_gp)
  for (i in seq_len(nrow(branch_source))) {
    s <- branch_source[i]
    grid.lines(x = unit(c(s$x0_mm, s$x1_mm), "mm"),
               y = unit(c(s$y0_mm, s$y1_mm), "mm"),
               gp = gpar(col = "#222222", lwd = line_width_pt, lineend = "butt"))
  }
  for (i in seq_len(nrow(leaf_table))) {
    grid.text(leaf_table$label[i], x = unit(labels_x_mm, "mm"),
              y = unit(leaf_table$y_mm[i], "mm"), just = c("left", "centre"),
              gp = gpar(fontfamily = "Arial", fontsize = label_font_pt, col = "black"))
  }
  tick_values <- c(max_height, max_height / 2, 0)
  tick_x <- vapply(tick_values, x_from_height, numeric(1))
  grid.lines(x = unit(c(tree_x_left_mm, tree_x_right_mm), "mm"),
             y = unit(c(axis_y_mm, axis_y_mm), "mm"),
             gp = gpar(col = "#333333", lwd = 0.35))
  for (i in seq_along(tick_values)) {
    grid.lines(x = unit(c(tick_x[i], tick_x[i]), "mm"),
               y = unit(c(axis_y_mm, axis_y_mm + 0.55), "mm"),
               gp = gpar(col = "#333333", lwd = 0.3))
    grid.text(formatC(tick_values[i], digits = 1, format = "f"),
              x = unit(tick_x[i], "mm"), y = unit(tick_label_y_mm, "mm"),
              gp = gpar(fontfamily = "Arial", fontsize = 3.4, col = "black"))
  }
  grid.text("Ward.D2 merge height (decreases toward leaves)",
            x = unit(page_width_mm / 2, "mm"), y = unit(axis_title_y_mm, "mm"),
            gp = gpar(fontfamily = "Arial", fontsize = 3.5, col = "black"))
  popViewport()
}

grDevices::png(filename = png_path, width = page_width_mm / 25.4,
               height = page_height_mm / 25.4, units = "in", res = 600,
               type = "cairo", family = "Arial", bg = "white")
draw_panel()
grDevices::dev.off()
grDevices::cairo_pdf(filename = pdf_path, width = page_width_mm / 25.4,
                     height = page_height_mm / 25.4, family = "Arial", bg = "white")
draw_panel()
grDevices::dev.off()

fwrite(branch_source, file.path(output_dir, "Fig3A_clean24_portrait_branch_source.csv"))
fwrite(leaf_table[, .(display_order_top_to_bottom, leaf_order,
                      Cluster, FinalNameWithStatus, label, y_mm)],
       file.path(output_dir, "Fig3A_clean24_portrait_leaf_source.csv"))
fwrite(merge_table, file.path(output_dir, "Fig3A_clean24_reconstructed_merge_source.csv"))
fwrite(data.table(
  item = c("distance_source", "linkage", "leaf_order_match", "cluster_n",
           "page_width_mm", "page_height_mm", "page_aspect_ratio",
           "placeholder_pixel_bounds_at_96dpi", "branch_segment_n",
           "source_resolution", "source_interpretation"),
  value = c(distance_path, "Ward.D2", "exactly matches saved bottom-to-top leaf_order",
            "24", as.character(page_width_mm), as.character(page_height_mm),
            as.character(page_ratio), "x=22..162; y=598..860 inclusive",
            as.character(nrow(branch_source)), "600 dpi PNG; editable vector PDF",
            "same frozen distances, topology and merge heights; portrait re-layout only")
), file.path(output_dir, "Fig3A_clean24_portrait_render_manifest.csv"))
message("Portrait Fig3A rendered from frozen distances; input profiles and distances were not recalculated.")
