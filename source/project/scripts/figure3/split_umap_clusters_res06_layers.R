# Split the frozen res0.6 UMAP into an Illustrator-aligned point layer and text layer.
# This reads the existing point/label sources only; it never loads or recomputes Seurat, UMAP, or clusters.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(grid)
  library(ragg)
  library(pdftools)
})

INPUT_FIG3_ROOT <- "E:/micro_output/fig3"
OUTPUT_ROOT <- "D:/data/microglia_new/tmp/figure3_umap_layer_probe/output_stage_no_ticks"
SOURCE_ROOT <- "E:/micro_output/v3_dims16/resolution_determine/source_data"
POINT_SOURCE <- file.path(SOURCE_ROOT, "umap_source_res0.6.csv")
LABEL_SOURCE <- file.path(SOURCE_ROOT, "umap_label_positions_res0.6.csv")
COLOR_SOURCE <- file.path(INPUT_FIG3_ROOT, "tables/cluster_color_unified_res0.6.csv")
OUTPUT_DIR <- OUTPUT_ROOT
SOURCE_OUT_DIR <- file.path(OUTPUT_ROOT, "source_data")
TABLE_OUT_DIR <- file.path(OUTPUT_ROOT, "tables")
POINTS_OUT <- file.path(OUTPUT_DIR, "umap_clusters_res0.6_points_only.png")
TEXT_OUT <- file.path(OUTPUT_DIR, "umap_clusters_res0.6_text_only.pdf")
POINT_SOURCE_OUT <- file.path(SOURCE_OUT_DIR, "umap_clusters_res0.6_points_only_source.csv")
TEXT_POSITION_OUT <- file.path(SOURCE_OUT_DIR, "umap_clusters_res0.6_text_positions_source.csv")
AXIS_ARROW_SOURCE_OUT <- file.path(SOURCE_OUT_DIR, "umap_clusters_res0.6_text_axis_arrows_source.csv")
AUDIT_OUT <- file.path(TABLE_OUT_DIR, "umap_clusters_res0.6_layer_export_audit.csv")
ORIGINAL_PNG <- file.path(INPUT_FIG3_ROOT, "figures/umap_clusters_res0.6.png")
TEMP_DIR <- "D:/data/microglia_new/tmp/figure3_umap_layer_probe"
TEMP_REFERENCE <- file.path(TEMP_DIR, "umap_clusters_res0.6_original_device_reference.png")
WIDTH_IN <- 7.4
HEIGHT_IN <- 6.0
DPI <- 600L
CLUSTER_LEVELS <- as.character(0:29)

required_files <- c(POINT_SOURCE, LABEL_SOURCE, COLOR_SOURCE)
if (!all(file.exists(required_files))) {
  stop("Missing frozen input(s): ", paste(required_files[!file.exists(required_files)], collapse = "; "))
}
if (any(file.exists(c(POINTS_OUT, TEXT_OUT, POINT_SOURCE_OUT, TEXT_POSITION_OUT,
                      AXIS_ARROW_SOURCE_OUT, AUDIT_OUT)))) {
  stop("A requested output already exists; refusing to overwrite.")
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TEMP_DIR, recursive = TRUE, showWarnings = FALSE)

points <- fread(POINT_SOURCE)
label_anchors <- fread(LABEL_SOURCE)
color_table <- fread(COLOR_SOURCE)
point_required <- c("CellID", "UMAP_1", "UMAP_2", "Cluster")
label_required <- c("Cluster", "UMAP_1", "UMAP_2", "CellNumber")
color_required <- c("Cluster", "Cluster_Color_Hex")
if (!all(point_required %in% names(points)) || !all(label_required %in% names(label_anchors)) ||
    !all(color_required %in% names(color_table))) {
  stop("Frozen source schema mismatch.")
}
if (nrow(points) != 424222L || anyDuplicated(points$CellID) ||
    any(!is.finite(points$UMAP_1)) || any(!is.finite(points$UMAP_2))) {
  stop("Frozen point source failed the expected 424,222-cell integrity check.")
}
if (!setequal(unique(as.character(points$Cluster)), CLUSTER_LEVELS) ||
    !setequal(unique(as.character(label_anchors$Cluster)), CLUSTER_LEVELS) ||
    nrow(label_anchors) != 30L || nrow(color_table) != 30L ||
    anyDuplicated(color_table$Cluster)) {
  stop("Frozen res0.6 cluster universe or palette does not match clusters 0-29.")
}
if (!all(grepl("^#[0-9A-Fa-f]{6}$", color_table$Cluster_Color_Hex))) {
  stop("Invalid cluster color mapping.")
}

# Preserve source row/point draw order and use the frozen cluster-to-color table.
points[, Cluster := factor(as.character(Cluster), levels = CLUSTER_LEVELS)]
label_anchors[, Cluster := factor(as.character(Cluster), levels = CLUSTER_LEVELS)]
cluster_colors <- setNames(color_table$Cluster_Color_Hex, as.character(color_table$Cluster))

# Reconstruct only the existing plot from frozen coordinates and the original plotting contract.
p_umap <- ggplot(points, aes(UMAP_1, UMAP_2, colour = Cluster)) +
  geom_point(size = 0.22, alpha = 0.55, stroke = 0) +
  scale_colour_manual(values = cluster_colors, drop = FALSE) +
  coord_fixed(expand = FALSE) +
  labs(
    title = "Harmony UMAP - resolution 0.6",
    subtitle = "30 clusters; labels shown at cluster medians",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme_classic(base_family = "Arial", base_size = 8) +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.text = element_text(colour = "black"),
    plot.title = element_text(size = 9, face = "bold"),
    plot.subtitle = element_text(size = 7.5, colour = "#4D4D4D"),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  geom_label_repel(
    data = label_anchors,
    aes(UMAP_1, UMAP_2, label = Cluster),
    inherit.aes = FALSE,
    family = "Arial",
    size = 2.7,
    fontface = "bold",
    colour = "black",
    fill = scales::alpha("white", 0.88),
    label.size = 0.18,
    box.padding = 0.25,
    point.padding = 0.1,
    min.segment.length = 0,
    seed = 42,
    max.overlaps = Inf
  )

# Draw once at the exact original 7.4 x 6.0 inch, 600 dpi ragg geometry,
# force ggrepel, then capture the laid-out grob before removing non-target layers.
ragg::agg_png(TEMP_REFERENCE, width = WIDTH_IN, height = HEIGHT_IN,
              units = "in", res = DPI, background = "white")
grid.newpage()
grid.draw(ggplotGrob(p_umap))
grid.force()
laid_out_grob <- grid.grab(wrap = TRUE, warn = FALSE)
grDevices::dev.off()
if (!file.exists(TEMP_REFERENCE) || file.info(TEMP_REFERENCE)$size <= 0) {
  stop("Original-device layout capture did not complete.")
}
if (!file.exists(ORIGINAL_PNG)) stop("Original reference PNG is missing.")
reference_png_byte_identical <- identical(
  unname(tools::md5sum(TEMP_REFERENCE)),
  unname(tools::md5sum(ORIGINAL_PNG))
)

filter_grob <- function(grob, keep) {
  if (keep(grob)) return(grob)
  if (!inherits(grob, "gTree")) return(NULL)
  children <- as.list(grob$children)
  if (!length(children)) return(NULL)
  children <- lapply(children, filter_grob, keep = keep)
  children <- Filter(Negate(is.null), children)
  if (!length(children)) return(NULL)
  grob$children <- do.call(grid::gList, children)
  grob
}

leaf_classes <- function(grob) {
  if (!inherits(grob, "gTree")) return(class(grob)[[1L]])
  children <- as.list(grob$children)
  if (!length(children)) return(character())
  unlist(lapply(children, leaf_classes), use.names = FALSE)
}

points_grob <- filter_grob(laid_out_grob, function(x) inherits(x, "points"))
text_grob <- filter_grob(laid_out_grob, function(x) {
  is_repel_label <- inherits(x, "text") && grepl("^textrepelgrob[0-9]+$", x$name)
  is_plot_heading <- inherits(x, "gTree") && grepl("^plot\\.(title|subtitle)\\.\\.titleGrob", x$name)
  is_repel_label || is_plot_heading
})
point_leaf_classes <- unique(leaf_classes(points_grob))
text_leaf_classes <- unique(leaf_classes(text_grob))
if (!length(point_leaf_classes) || !all(point_leaf_classes == "points")) {
  stop("Points layer contains a non-point leaf grob.")
}
if (!length(text_leaf_classes) || !all(text_leaf_classes == "text")) {
  stop("Text layer contains a non-text leaf grob.")
}

# The two deliverables are drawn from the same captured layout/viewports.
ragg::agg_png(POINTS_OUT, width = WIDTH_IN, height = HEIGHT_IN,
              units = "in", res = DPI, background = "white")
grid.newpage()
grid.draw(points_grob)
grDevices::dev.off()

grDevices::cairo_pdf(TEXT_OUT, width = WIDTH_IN, height = HEIGHT_IN,
                     bg = "transparent", family = "Arial")
grid.newpage()
grid.draw(text_grob)
arrow_lwd <- 0.7
arrow_head_in <- 0.07
grid.lines(
  x = unit(c(0.14, 0.24), "npc"), y = unit(c(0.14, 0.14), "npc"),
  arrow = grid::arrow(type = "closed", ends = "last", length = unit(arrow_head_in, "in")),
  gp = gpar(col = "black", lwd = arrow_lwd)
)
grid.lines(
  x = unit(c(0.14, 0.14), "npc"), y = unit(c(0.14, 0.24), "npc"),
  arrow = grid::arrow(type = "closed", ends = "last", length = unit(arrow_head_in, "in")),
  gp = gpar(col = "black", lwd = arrow_lwd)
)
grid.text("X", x = unit(0.25, "npc"), y = unit(0.14, "npc"), just = "left",
          gp = gpar(fontfamily = "Arial", fontsize = 8, fontface = "bold"))
grid.text("Y", x = unit(0.14, "npc"), y = unit(0.25, "npc"), just = "centre",
          gp = gpar(fontfamily = "Arial", fontsize = 8, fontface = "bold"))
grid.text("UMAP1", x = unit(0.19, "npc"), y = unit(0.105, "npc"), just = "centre",
          gp = gpar(fontfamily = "Arial", fontsize = 8))
grid.text("UMAP2", x = unit(0.105, "npc"), y = unit(0.19, "npc"), rot = 90,
          just = "centre", gp = gpar(fontfamily = "Arial", fontsize = 8))
grDevices::dev.off()

# Retain auditable point input and actual PDF text positions (including repelled labels).
point_source_out <- points[, .(
  CellID,
  UMAP_1,
  UMAP_2,
  Cluster = as.character(Cluster),
  Cluster_Color_Hex = unname(cluster_colors[as.character(Cluster)])
)]
fwrite(point_source_out, POINT_SOURCE_OUT)
text_boxes <- pdftools::pdf_data(TEXT_OUT)
if (length(text_boxes) != 1L || !all(c("x", "y", "width", "height", "text") %in% names(text_boxes[[1L]]))) {
  stop("Could not extract text positions from the one-page PDF.")
}
text_positions <- as.data.table(text_boxes[[1L]])
text_positions[, `:=`(
  Text_Order = .I,
  Page_Width_pt = WIDTH_IN * 72,
  Page_Height_pt = HEIGHT_IN * 72,
  Page_Width_in = WIDTH_IN,
  Page_Height_in = HEIGHT_IN,
  Position_Origin = "PDF page top-left; x/y/width/height in PDF points"
)]
setcolorder(text_positions, c("Text_Order", "text", "x", "y", "width", "height",
                              "space", "Page_Width_pt", "Page_Height_pt",
                              "Page_Width_in", "Page_Height_in", "Position_Origin"))
fwrite(text_positions, TEXT_POSITION_OUT)

axis_arrow_source <- data.table(
  axis = c("UMAP1", "UMAP2"), direction = c("X", "Y"),
  x_start_npc = c(0.14, 0.14), y_start_npc = c(0.14, 0.14),
  x_end_npc = c(0.24, 0.14), y_end_npc = c(0.14, 0.24),
  arrow_head = "closed", line_width_pt = arrow_lwd,
  arrow_head_length_in = arrow_head_in,
  canvas_width_in = WIDTH_IN, canvas_height_in = HEIGHT_IN
)
fwrite(axis_arrow_source, AXIS_ARROW_SOURCE_OUT)

pdf_text <- paste(pdftools::pdf_text(TEXT_OUT), collapse = "\n")
required_text <- c("Harmony UMAP - resolution 0.6",
                   "30 clusters; labels shown at cluster medians", "UMAP1", "UMAP2")
if (!all(vapply(required_text, grepl, logical(1), x = pdf_text, fixed = TRUE))) {
  stop("The PDF is missing one or more required title/axis text strings.")
}
pdf_atoms <- trimws(text_positions$text)
if (!all(CLUSTER_LEVELS %in% pdf_atoms) || !all(c("X", "Y", "UMAP1", "UMAP2") %in% pdf_atoms)) {
  stop("One or more cluster labels 0-29 are absent from the text PDF.")
}

audit <- data.table(
  item = c("frozen_cells", "frozen_clusters", "canvas_width_in", "canvas_height_in",
           "png_width_px", "png_height_px", "dpi", "point_leaf_classes",
           "text_leaf_classes", "pdf_text_spans", "cluster_labels_present",
           "source_point_rows", "text_positions_source_rows", "axis_arrows",
           "axis_arrow_source_rows"),
  value = c(nrow(points), length(CLUSTER_LEVELS), WIDTH_IN, HEIGHT_IN,
            WIDTH_IN * DPI, HEIGHT_IN * DPI, DPI,
            paste(point_leaf_classes, collapse = ";"),
            paste(text_leaf_classes, collapse = ";"), nrow(text_positions),
            sum(CLUSTER_LEVELS %in% pdf_atoms), nrow(point_source_out), nrow(text_positions),
            "X/Y arrows with UMAP1/UMAP2 labels", nrow(axis_arrow_source)),
  status = rep("verified", 15L)
)
audit <- rbind(audit, data.table(
  item = c("full_original_device_png_byte_identical", "original_png_path"),
  value = c(reference_png_byte_identical, ORIGINAL_PNG),
  status = c(if (reference_png_byte_identical) "verified" else "not_identical", "reference")
))
fwrite(audit, AUDIT_OUT)

unlink(TEMP_REFERENCE)
cat("COMPLETE: point PNG and text-only editable PDF; cells=", nrow(points),
    "; clusters=", length(CLUSTER_LEVELS),
    "; canvas=", WIDTH_IN, "x", HEIGHT_IN, "in; dpi=", DPI,
    "; point leaves=", paste(point_leaf_classes, collapse = ","),
    "; text leaves=", paste(text_leaf_classes, collapse = ","),
    "; PDF text spans=", nrow(text_positions), "\n", sep = "")
