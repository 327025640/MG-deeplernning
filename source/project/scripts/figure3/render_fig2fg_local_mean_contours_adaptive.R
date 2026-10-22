# Re-render Fig2F/G from saved local-mean grids only.
# Does not recalculate scores, smoothing, coordinates, or support masks.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(ragg)
  library(pdftools)
  library(digest)
})

INPUT_ROOT <- "D:/data/microglia_new/tmp/fig2fg_local_mean_contours_stage_v6/source_data"
OUTPUT_ROOT <- "D:/data/microglia_new/tmp/fig2fg_contrast_render_v9"
FIGURE_OUT <- file.path(OUTPUT_ROOT, "figures")
SOURCE_OUT <- file.path(OUTPUT_ROOT, "source_data")
TABLE_OUT <- file.path(OUTPUT_ROOT, "tables")
WIDTH_IN <- 5.5
HEIGHT_IN <- 5.2
DPI <- 600L
GRID_N <- 400L
GRID_HASHES <- c(
  Fig2f_microglia_M1like_AddModuleScore_UMAP =
    "BA8E95E0B95D81A2817967A9783CE31DFE13905E454E6D9F319CAB0AB255AF09",
  Fig2g_microglia_M2like_AddModuleScore_UMAP =
    "9B020A3A5AB920C2B8D192543C9852987D9DC38866D4FE037C101BF01B723420"
)
LINE_ALPHA <- 0.20
LINE_WIDTH_MM <- 0.08

jobs <- list(
  list(
    stem = "Fig2f_microglia_M1like_AddModuleScore_UMAP",
    module = "M1-like",
    source = file.path(INPUT_ROOT,
      "Fig2f_microglia_M1like_AddModuleScore_UMAP_local_mean_contour_grid_smoothing_source.csv"),
    limits = c(-0.40, 1.70),
    ticks = c(-0.4, 0.2, 0.8, 1.4),
    palette = c("#F1F7FC", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
    line_color = "#174A78",
    title = "M1-like local mean score"
  ),
  list(
    stem = "Fig2g_microglia_M2like_AddModuleScore_UMAP",
    module = "M2-like",
    source = file.path(INPUT_ROOT,
      "Fig2g_microglia_M2like_AddModuleScore_UMAP_local_mean_contour_grid_smoothing_source.csv"),
    limits = c(-0.30, 0.40),
    ticks = c(-0.3, -0.1, 0.1, 0.3),
    palette = c("#FFF5F0", "#FDC3AD", "#FB8060", "#D73027", "#8E1029"),
    line_color = "#6F1323",
    title = "M2-like local mean score"
  )
)

dir.create(FIGURE_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_OUT, recursive = TRUE, showWarnings = FALSE)

filter_leaf_grob <- function(grob, keep_leaf) {
  if (inherits(grob, "gTree")) {
    children <- as.list(grob$children)
    if (!length(children)) return(NULL)
    children <- lapply(children, filter_leaf_grob, keep_leaf = keep_leaf)
    children <- Filter(Negate(is.null), children)
    if (!length(children)) return(NULL)
    grob$children <- do.call(grid::gList, children)
    return(grob)
  }
  if (isTRUE(keep_leaf(grob))) grob else NULL
}

filter_named_subtrees <- function(grob, pattern, drop_background = FALSE) {
  clean_tree <- function(x) {
    if (drop_background && inherits(x, "rect") &&
        grepl("background", x$name, ignore.case = TRUE)) return(NULL)
    if (!inherits(x, "gTree")) return(x)
    children <- as.list(x$children)
    if (!length(children)) return(NULL)
    children <- lapply(children, clean_tree)
    children <- Filter(Negate(is.null), children)
    if (!length(children)) return(NULL)
    x$children <- do.call(grid::gList, children)
    x
  }
  if (inherits(grob, "gTree") && grepl(pattern, grob$name)) return(clean_tree(grob))
  if (!inherits(grob, "gTree")) return(NULL)
  children <- as.list(grob$children)
  if (!length(children)) return(NULL)
  children <- lapply(children, filter_named_subtrees, pattern = pattern,
                     drop_background = drop_background)
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

format_tick <- function(x) formatC(x, format = "f", digits = 1L)

audit_rows <- list()
config_rows <- list()

for (job in jobs) {
  if (!file.exists(job$source)) stop("Missing frozen smoothed-grid CSV: ", job$source)
  actual_file_hash <- toupper(digest::digest(file = job$source, algo = "sha256"))
  expected_file_hash <- unname(GRID_HASHES[job$stem])
  if (!identical(actual_file_hash, expected_file_hash)) {
    stop("Smoothed-grid source hash changed; refusing to render: ", job$source)
  }

  grid_dt <- fread(job$source)
  required <- c("UMAP_1", "UMAP_2", "score_value", "Neff",
                "nearest_cell_distance", "support_mask")
  if (!all(required %in% names(grid_dt)) || nrow(grid_dt) != GRID_N^2L ||
      anyNA(grid_dt$support_mask)) stop("Smoothed-grid source schema/dimensions failed: ", job$source)
  if (!is.logical(grid_dt$support_mask)) {
    grid_dt[, support_mask := tolower(as.character(support_mask)) == "true"]
  }
  if (grid_dt[, uniqueN(UMAP_1)] != GRID_N ||
      grid_dt[, uniqueN(UMAP_2)] != GRID_N ||
      grid_dt[support_mask == TRUE, any(!is.finite(score_value))] ||
      grid_dt[support_mask == FALSE, any(is.finite(score_value))]) {
    stop("Grid coordinates, supported values, or blank mask values failed integrity checks.")
  }

  mask_hash_before <- digest::digest(grid_dt$support_mask, algo = "sha256")
  score_hash_before <- digest::digest(grid_dt$score_value, algo = "sha256")
  support_n <- sum(grid_dt$support_mask)
  values <- grid_dt[support_mask == TRUE, score_value]
  actual_range <- range(values)
  if (job$limits[1L] > actual_range[1L] || job$limits[2L] < actual_range[2L]) {
    stop("Panel-specific scale would clip supported values for ", job$module)
  }
  if (length(job$ticks) < 4L || length(job$ticks) > 5L ||
      any(job$ticks < job$limits[1L] | job$ticks > job$limits[2L])) {
    stop("Invalid panel-specific colorbar ticks for ", job$module)
  }

  # Four subtle contours span the full observed range; raster interpolation is
  # display-only antialiasing of the saved grid. Unsupported cells remain NA.
  line_levels <- seq(actual_range[1L], actual_range[2L], length.out = 7L)[2:6]
  ramp <- grDevices::colorRamp(job$palette, space = "Lab")
  ramp_rgb <- ramp(seq(0, 1, length.out = 256L))
  palette_256 <- grDevices::rgb(ramp_rgb[, 1L], ramp_rgb[, 2L], ramp_rgb[, 3L],
                                maxColorValue = 255)
  p <- ggplot() +
    geom_raster(data = grid_dt,
      aes(x = UMAP_1, y = UMAP_2, fill = score_value),
      interpolate = TRUE, na.rm = TRUE) +
    geom_contour(data = grid_dt,
      aes(x = UMAP_1, y = UMAP_2, z = score_value),
      breaks = line_levels, colour = scales::alpha(job$line_color, LINE_ALPHA),
      linewidth = LINE_WIDTH_MM, na.rm = TRUE, show.legend = FALSE) +
    scale_fill_gradientn(
      name = paste0(job$module, " local mean\nAddModuleScore\nPanel-specific scale"),
      colours = palette_256, values = seq(0, 1, length.out = 256L), space = "Lab",
      limits = job$limits, breaks = job$ticks, labels = format_tick,
      oob = scales::oob_censor, na.value = "white",
      guide = guide_colorbar(frame.colour = NA, ticks.colour = "black")
    ) +
    coord_equal(
      xlim = range(grid_dt$UMAP_1) + c(0.30, -0.30),
      ylim = range(grid_dt$UMAP_2) + c(0.30, -0.30),
      expand = TRUE
    ) +
    labs(title = job$title, x = "UMAP1", y = "UMAP2") +
    theme_classic(base_family = "Arial", base_size = 9) +
    theme(
      axis.text = element_blank(), axis.ticks = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.title = element_text(size = 8, colour = "black"),
      plot.title = element_text(face = "bold", size = 10, hjust = 0, colour = "black"),
      legend.title = element_text(size = 8, colour = "black"),
      legend.text = element_text(size = 7, colour = "black"),
      legend.key.height = unit(0.65, "cm"),
      legend.key.width = unit(0.42, "cm"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      aspect.ratio = 1
    )

  prefix <- paste0(job$stem, "_local_mean_contour_contrast")
  contour_png <- file.path(FIGURE_OUT, paste0(prefix, "_only.png"))
  text_pdf <- file.path(FIGURE_OUT, paste0(prefix, "_text_only.pdf"))
  legend_pdf <- file.path(FIGURE_OUT, paste0(prefix, "_legend_only.pdf"))
  legend_csv <- file.path(SOURCE_OUT, paste0(prefix, "_legend_source.csv"))
  audit_csv <- file.path(TABLE_OUT, paste0(prefix, "_audit.csv"))
  config_csv <- file.path(TABLE_OUT, paste0(prefix, "_config.csv"))
  if (any(file.exists(c(contour_png, text_pdf, legend_pdf, legend_csv, audit_csv, config_csv)))) {
    stop("Refusing to overwrite adaptive staging outputs: ", prefix)
  }

  ragg::agg_png(contour_png, width = WIDTH_IN, height = HEIGHT_IN,
                units = "in", res = DPI, background = "white")
  grid.newpage()
  grid.draw(ggplotGrob(p))
  grid.force()
  full_grob <- grid.grab(wrap = TRUE, warn = FALSE)
  grDevices::dev.off()

  text_pattern <- "(plot\\.title\\.\\.titleGrob|axis\\.title\\.[xy]\\.(bottom|left)\\.\\.titleGrob)"
  text_grob <- filter_named_subtrees(full_grob, text_pattern, drop_background = TRUE)
  text_classes <- unique(leaf_classes(text_grob))
  if (!length(text_classes) || !all(text_classes == "text")) {
    stop("Text-only PDF would contain non-text grobs for ", job$module)
  }
  legend_grob <- filter_named_subtrees(full_grob, "^guide-box-right", drop_background = TRUE)
  legend_classes <- unique(leaf_classes(legend_grob))
  if (!length(legend_classes) || any(legend_classes == "points")) {
    stop("Legend-only PDF missing or contains point data for ", job$module)
  }
  panel_grob <- filter_named_subtrees(full_grob, "^panel-1")
  contour_grob <- filter_leaf_grob(panel_grob, function(x) {
    grepl("rastergrob", x$name, fixed = TRUE) || inherits(x, "polyline")
  })
  contour_classes <- unique(leaf_classes(contour_grob))
  if (!all(c("rastergrob", "polyline") %in% contour_classes)) {
    stop("Contour-only layer does not contain smoothed raster and contour lines: ",
         paste(contour_classes, collapse = ";"))
  }

  ragg::agg_png(contour_png, width = WIDTH_IN, height = HEIGHT_IN,
                units = "in", res = DPI, background = "white")
  grid.newpage(); grid.draw(contour_grob); grDevices::dev.off()
  grDevices::cairo_pdf(text_pdf, width = WIDTH_IN, height = HEIGHT_IN,
                       bg = "transparent", family = "Arial")
  grid.newpage(); grid.draw(text_grob); grDevices::dev.off()
  grDevices::cairo_pdf(legend_pdf, width = WIDTH_IN, height = HEIGHT_IN,
                       bg = "transparent", family = "Arial")
  grid.newpage(); grid.draw(legend_grob); grDevices::dev.off()

  text_content <- paste(pdftools::pdf_text(text_pdf), collapse = "\n")
  legend_content <- paste(pdftools::pdf_text(legend_pdf), collapse = "\n")
  if (!all(vapply(c(job$title, "UMAP1", "UMAP2"), grepl, logical(1),
                  x = text_content, fixed = TRUE))) stop("Text labels missing for ", job$module)
  if (!grepl("Panel-specific scale", legend_content, fixed = TRUE) ||
      !grepl("AddModuleScore", legend_content, fixed = TRUE)) {
    stop("Independent-scale wording missing from legend for ", job$module)
  }
  legend_data <- as.data.table(pdftools::pdf_data(legend_pdf)[[1L]])
  legend_data[, text_order := .I]
  fwrite(legend_data, legend_csv)
  text_fonts <- pdftools::pdf_fonts(text_pdf)
  legend_fonts <- pdftools::pdf_fonts(legend_pdf)
  if (!nrow(text_fonts) || !all(text_fonts$embedded) ||
      !nrow(legend_fonts) || !all(legend_fonts$embedded)) {
    stop("Arial embedding check failed for ", job$module)
  }

  mask_hash_after <- digest::digest(grid_dt$support_mask, algo = "sha256")
  score_hash_after <- digest::digest(grid_dt$score_value, algo = "sha256")
  if (!identical(mask_hash_before, mask_hash_after) ||
      !identical(score_hash_before, score_hash_after)) {
    stop("Input grid values or support mask changed during render for ", job$module)
  }
  audit <- data.table(
    item = c("input_grid_csv", "input_grid_sha256", "support_mask_sha256",
      "score_vector_sha256", "grid_nx", "grid_ny", "supported_grid_cells",
      "supported_grid_fraction", "actual_supported_min", "actual_supported_max",
      "panel_scale_min", "panel_scale_max", "panel_ticks", "palette_low_mid_high",
      "palette_interpolation", "palette_sample_count",
      "fill_interpolation", "contour_line_levels", "contour_line_color",
      "contour_line_alpha", "contour_line_width_mm", "contour_layer_grobs",
      "text_layer_grobs", "legend_layer_grobs", "text_pdf_pages", "legend_pdf_pages",
      "text_pdf_embedded_fonts", "legend_pdf_embedded_fonts", "grid_values_unchanged",
      "mask_unchanged"),
    value = c(job$source, actual_file_hash, mask_hash_before, score_hash_before,
      GRID_N, GRID_N, support_n, support_n / nrow(grid_dt), actual_range[1L],
      actual_range[2L], job$limits[1L], job$limits[2L],
      paste(format_tick(job$ticks), collapse = ";"), paste(job$palette, collapse = ";"),
      "Lab", 256L,
      "bilinear raster interpolation for rendering only; source grid unchanged; unsupported cells NA",
      paste(format(signif(line_levels, 6), trim = TRUE), collapse = ";"),
      job$line_color, LINE_ALPHA, LINE_WIDTH_MM, paste(contour_classes, collapse = ";"),
      paste(text_classes, collapse = ";"), paste(legend_classes, collapse = ";"),
      pdftools::pdf_info(text_pdf)$pages, pdftools::pdf_info(legend_pdf)$pages,
      paste(text_fonts$name, collapse = ";"), paste(legend_fonts$name, collapse = ";"),
      identical(score_hash_before, score_hash_after), identical(mask_hash_before, mask_hash_after)),
    status = "recorded"
  )
  fwrite(audit, audit_csv)
  audit_rows[[job$stem]] <- audit
  config <- data.table(
    parameter = c("render_source", "score_range", "panel_scale_limits", "ticks",
      "palette_low_mid_high", "palette_interpolation", "palette_samples",
      "interpolation", "contour_levels", "contour_line_color",
      "contour_line_alpha", "contour_line_width_mm", "grid_nx", "grid_ny",
      "support_mask_rule", "grid_values_and_mask_modified"),
    value = c(job$source, paste(signif(actual_range, 8), collapse = ";"),
      paste(job$limits, collapse = ";"), paste(job$ticks, collapse = ";"),
      paste(job$palette, collapse = ";"), "Lab", 256L,
      "bilinear raster interpolation for rendering only; source grid unchanged",
      paste(signif(line_levels, 8), collapse = ";"), job$line_color, LINE_ALPHA,
      LINE_WIDTH_MM, GRID_N, GRID_N, "unchanged saved support_mask; unsupported cells remain NA",
      "no; render-only")
  )
  fwrite(config, config_csv)
  cat("RENDERED ", job$module, " range=", paste(signif(actual_range, 6), collapse = ","),
      " ticks=", paste(format_tick(job$ticks), collapse = ","),
      " mask=", support_n, "/", nrow(grid_dt), "\n", sep = "")
}
cat("ADAPTIVE_RENDER_COMPLETE\n")
