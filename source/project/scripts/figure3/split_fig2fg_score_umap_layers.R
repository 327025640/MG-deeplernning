# Split the frozen Fig2F/G score UMAPs into a point-only PNG and a non-point PDF layer.
# Reads only existing Fig2 source CSVs; no Seurat object or score recalculation.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(ragg)
  library(pdftools)
})

INPUT_ROOT <- "E:/micro_output/fig3"
OUTPUT_ROOT <- "D:/data/microglia_new/tmp/fig2fg_score_umap_layers_stage_v2"
FIGURE_OUT <- file.path(OUTPUT_ROOT, "figures")
SOURCE_OUT <- file.path(OUTPUT_ROOT, "source_data")
TABLE_OUT <- file.path(OUTPUT_ROOT, "tables")
WIDTH_IN <- 5.5
HEIGHT_IN <- 5.2
DPI <- 600L

jobs <- list(
  list(
    stem = "Fig2f_microglia_M1like_AddModuleScore_UMAP",
    source = file.path(INPUT_ROOT, "source_data/Fig2f_microglia_M1like_UMAP_source.csv"),
    raw_score = "M1_like_score",
    plot_score = "M1_like_score_plot",
    title = "M1-like score",
    low = "#F2F2F2",
    high = "#2166AC"
  ),
  list(
    stem = "Fig2g_microglia_M2like_AddModuleScore_UMAP",
    source = file.path(INPUT_ROOT, "source_data/Fig2g_microglia_M2like_UMAP_source.csv"),
    raw_score = "M2_like_score",
    plot_score = "M2_like_score_plot",
    title = "M2-like score",
    low = "#F2F2F2",
    high = "#B2182B"
  )
)

dir.create(FIGURE_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_OUT, recursive = TRUE, showWarnings = FALSE)

filter_grob <- function(grob, keep_leaf) {
  if (inherits(grob, "gTree")) {
    children <- as.list(grob$children)
    if (!length(children)) return(NULL)
    children <- lapply(children, filter_grob, keep_leaf = keep_leaf)
    children <- Filter(Negate(is.null), children)
    if (!length(children)) return(NULL)
    grob$children <- do.call(grid::gList, children)
    return(grob)
  }
  if (isTRUE(keep_leaf(grob))) grob else NULL
}

leaf_classes <- function(grob) {
  if (!inherits(grob, "gTree")) return(class(grob)[[1L]])
  children <- as.list(grob$children)
  if (!length(children)) return(character())
  unlist(lapply(children, leaf_classes), use.names = FALSE)
}

theme_score_umap <- theme_classic(base_family = "Arial", base_size = 9) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.35),
    axis.title = element_text(size = 8, colour = "black"),
    plot.title = element_text(face = "bold", size = 10, hjust = 0, colour = "black"),
    plot.subtitle = element_text(size = 7, colour = "grey30"),
    legend.title = element_text(size = 8, colour = "black"),
    legend.text = element_text(size = 7, colour = "black"),
    legend.key.height = unit(0.65, "cm"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    aspect.ratio = 1
  )

audit_rows <- list()
for (job in jobs) {
  if (!file.exists(job$source)) stop("Missing frozen source CSV: ", job$source)
  points <- fread(job$source)
  required <- c("CellID", "UMAP_1", "UMAP_2", job$raw_score, job$plot_score)
  if (!all(required %in% names(points))) stop("Source schema mismatch: ", job$source)
  if (!nrow(points) || anyDuplicated(points$CellID) ||
      any(!is.finite(points$UMAP_1)) || any(!is.finite(points$UMAP_2)) ||
      any(!is.finite(points[[job$plot_score]]))) {
    stop("Invalid cell coordinates or score values: ", job$source)
  }

  score_limits <- range(points[[job$plot_score]])
  score_breaks <- c(score_limits[1L], mean(score_limits), score_limits[2L])
  score_labels <- round(score_breaks, 2)
  p <- ggplot(points, aes(x = UMAP_1, y = UMAP_2, colour = .data[[job$plot_score]])) +
    geom_point(size = 0.10, alpha = 0.85, shape = 16, stroke = 0) +
    scale_colour_gradient(
      name = gsub(" ", "\n", job$title, fixed = TRUE),
      low = job$low, high = job$high, limits = score_limits,
      breaks = score_breaks, labels = score_labels
    ) +
    coord_equal() +
    labs(title = job$title, x = "UMAP 1", y = "UMAP 2") +
    theme_score_umap

  full_png <- file.path(OUTPUT_ROOT, paste0(job$stem, "_reference.png"))
  points_png <- file.path(FIGURE_OUT, paste0(job$stem, "_points_only.png"))
  nonpoint_pdf <- file.path(FIGURE_OUT, paste0(job$stem, "_text_legend_only.pdf"))
  point_source <- file.path(SOURCE_OUT, paste0(job$stem, "_points_source.csv"))
  text_source <- file.path(SOURCE_OUT, paste0(job$stem, "_text_positions_source.csv"))
  audit_file <- file.path(TABLE_OUT, paste0(job$stem, "_layer_export_audit.csv"))
  outputs <- c(full_png, points_png, nonpoint_pdf, point_source, text_source, audit_file)
  if (any(file.exists(outputs))) stop("Refusing to overwrite staged output for ", job$stem)

  ragg::agg_png(full_png, width = WIDTH_IN, height = HEIGHT_IN,
                units = "in", res = DPI, background = "white")
  grid.newpage()
  grid.draw(ggplotGrob(p))
  grid.force()
  captured <- grid.grab(wrap = TRUE, warn = FALSE)
  grDevices::dev.off()

  point_grob <- filter_grob(captured, function(x) inherits(x, "points"))
  nonpoint_grob <- filter_grob(captured, function(x) {
    if (inherits(x, "points")) return(FALSE)
    if (inherits(x, "rect") && grepl("background", x$name, ignore.case = TRUE)) return(FALSE)
    TRUE
  })
  point_classes <- unique(leaf_classes(point_grob))
  nonpoint_classes <- unique(leaf_classes(nonpoint_grob))
  if (!length(point_classes) || !all(point_classes == "points")) stop("Point layer contains non-point grobs.")
  if (!length(nonpoint_classes) || any(nonpoint_classes == "points")) stop("Non-point layer still contains points.")

  ragg::agg_png(points_png, width = WIDTH_IN, height = HEIGHT_IN,
                units = "in", res = DPI, background = "white")
  grid.newpage()
  grid.draw(point_grob)
  grDevices::dev.off()

  cairo_pdf(nonpoint_pdf, width = WIDTH_IN, height = HEIGHT_IN,
            bg = "transparent", family = "Arial")
  grid.newpage()
  grid.draw(nonpoint_grob)
  grDevices::dev.off()

  keep_cols <- intersect(
    c("CellID", "UMAP_1", "UMAP_2", job$raw_score, job$plot_score,
      "Cluster", "Disease", "DatasetID", "SampleID_Global"),
    names(points)
  )
  fwrite(points[, ..keep_cols], point_source)

  pdf_boxes <- pdf_data(nonpoint_pdf)
  if (length(pdf_boxes) != 1L) stop("Expected a single-page PDF: ", nonpoint_pdf)
  text_dt <- as.data.table(pdf_boxes[[1L]])
  text_dt[, `:=`(Text_Order = .I, Page_Width_in = WIDTH_IN, Page_Height_in = HEIGHT_IN)]
  setcolorder(text_dt, c("Text_Order", "text", "x", "y", "width", "height", "space",
                         "Page_Width_in", setdiff(names(text_dt),
                           c("Text_Order", "text", "x", "y", "width", "height", "space",
                             "Page_Width_in"))))
  fwrite(text_dt, text_source)

  pdf_content <- paste(pdf_text(nonpoint_pdf), collapse = "\n")
  required_strings <- c(job$title, "UMAP 1", "UMAP 2", "score")
  if (!all(vapply(required_strings, grepl, logical(1), x = pdf_content, fixed = TRUE))) {
    stop("Non-point PDF is missing title, axis text, or score legend: ", nonpoint_pdf)
  }
  if (length(nonpoint_classes) == 0L || !nrow(text_dt)) stop("Empty non-point layer.")
  font_data <- pdf_fonts(nonpoint_pdf)
  if (!nrow(font_data) || !all(font_data$embedded)) stop("PDF fonts are not embedded: ", nonpoint_pdf)

  audit <- data.table(
    item = c("cells", "canvas_width_in", "canvas_height_in", "png_width_px", "png_height_px",
             "dpi", "score_column", "score_plot_min", "score_plot_max", "point_grob_classes",
             "nonpoint_grob_classes", "pdf_text_spans", "pdf_pages", "embedded_fonts"),
    value = c(nrow(points), WIDTH_IN, HEIGHT_IN, WIDTH_IN * DPI, HEIGHT_IN * DPI, DPI,
              job$plot_score, score_limits[1L], score_limits[2L],
              paste(point_classes, collapse = ";"), paste(nonpoint_classes, collapse = ";"),
              nrow(text_dt), pdf_info(nonpoint_pdf)$pages, paste(font_data$name, collapse = ";")),
    status = rep("verified", 14L)
  )
  fwrite(audit, audit_file)
  audit_rows[[job$stem]] <- audit

  cat("COMPLETE ", job$stem, ": cells=", nrow(points),
      "; score_limits=", paste(signif(score_limits, 5), collapse = ","),
      "; PDF text spans=", nrow(text_dt),
      "; nonpoint grobs=", paste(nonpoint_classes, collapse = ","), "\n", sep = "")
  rm(points, p, captured, point_grob, nonpoint_grob, text_dt)
  invisible(gc())
}
