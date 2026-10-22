# 依照 Fig3 预留黑框尺寸导出 res0.6 UMAP 独立图例；只读取冻结的颜色映射。
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(ragg)
  library(svglite)
  library(pdftools)
})

COLOR_SOURCE <- "E:/micro_output/fig3/tables/cluster_color_unified_res0.6.csv"
OUTPUT_ROOT <- "D:/data/microglia_new/tmp/umap_cluster_legend_res06_stage_20260930_v3"
FIGURE_DIR <- file.path(OUTPUT_ROOT, "figures")
SOURCE_DIR <- file.path(OUTPUT_ROOT, "source_data")
TABLE_DIR <- file.path(OUTPUT_ROOT, "tables")

PNG_OUT <- file.path(FIGURE_DIR, "umap_clusters_res0.6_legend_blackbox_fit.png")
PDF_OUT <- file.path(FIGURE_DIR, "umap_clusters_res0.6_legend_blackbox_fit.pdf")
SVG_OUT <- file.path(FIGURE_DIR, "umap_clusters_res0.6_legend_blackbox_fit.svg")
SOURCE_OUT <- file.path(SOURCE_DIR, "umap_clusters_res0.6_legend_blackbox_fit_source.csv")
AUDIT_OUT <- file.path(TABLE_DIR, "umap_clusters_res0.6_legend_blackbox_fit_audit.csv")

# Fig3.jpg 的黑框约 97 × 227 px；按其 72 ppi 输出尺寸固定图例画布。
REFERENCE_DPI <- 72
TARGET_WIDTH_PX <- 97
TARGET_HEIGHT_PX <- 227
WIDTH_IN <- TARGET_WIDTH_PX / REFERENCE_DPI
HEIGHT_IN <- TARGET_HEIGHT_PX / REFERENCE_DPI
PNG_DPI <- 600
CLUSTER_LEVELS <- 0:29

if (!file.exists(COLOR_SOURCE)) stop("Missing frozen cluster color table: ", COLOR_SOURCE)
if (any(file.exists(c(PNG_OUT, PDF_OUT, SVG_OUT, SOURCE_OUT, AUDIT_OUT)))) {
  stop("At least one staged output already exists; refusing to overwrite.")
}

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

color_table <- fread(COLOR_SOURCE)
required_columns <- c("Cluster", "Cluster_Color_Hex")
if (!all(required_columns %in% names(color_table))) {
  stop("Frozen cluster color table schema mismatch.")
}
if (nrow(color_table) != length(CLUSTER_LEVELS) ||
    anyDuplicated(as.character(color_table$Cluster)) ||
    !setequal(as.integer(as.character(color_table$Cluster)), CLUSTER_LEVELS) ||
    any(!grepl("^#[0-9A-Fa-f]{6}$", color_table$Cluster_Color_Hex))) {
  stop("Frozen cluster color table must contain exactly clusters 0-29 with valid hex colors.")
}

legend_data <- color_table[, .(
  Cluster = as.integer(as.character(Cluster)),
  Cluster_Color_Hex = toupper(as.character(Cluster_Color_Hex))
)]
setorder(legend_data, Cluster)
legend_data[, `:=`(
  Label = as.character(Cluster),
  Legend_Row = 15L - (Cluster %/% 2L),
  Legend_Column = (Cluster %% 2L) + 1L
)]
legend_data[, `:=`(
  Swatch_X = fifelse(Legend_Column == 1L, 0.16, 1.12),
  Label_X = fifelse(Legend_Column == 1L, 0.30, 1.26),
  Y = as.numeric(Legend_Row)
)]

legend_plot <- ggplot(legend_data, aes(x = Swatch_X, y = Y)) +
  geom_point(
    aes(colour = Cluster_Color_Hex),
    size = 2.5,
    shape = 16,
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = Label_X, label = Label),
    family = "Arial",
    size = 6.2 * 25.4 / 72,
    colour = "black",
    hjust = 0,
    vjust = 0.5,
    show.legend = FALSE
  ) +
  scale_colour_identity() +
  annotate(
    "text", x = 1, y = 15.92, label = "Cluster",
    family = "Arial", fontface = "bold", size = 7 * 25.4 / 72,
    colour = "black", hjust = 0.5, vjust = 0.5
  ) +
  coord_cartesian(xlim = c(0, 2), ylim = c(0.35, 16.35), expand = FALSE, clip = "off") +
  theme_void(base_family = "Arial", base_size = 6.2) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(3.5, 4, 3.5, 4, unit = "pt")
  )

# 导出均使用同一 R 图对象和固定画布尺寸，不裁边、不重新排列簇色。
ragg::agg_png(
  filename = PNG_OUT, width = WIDTH_IN, height = HEIGHT_IN,
  units = "in", res = PNG_DPI, background = "white"
)
print(legend_plot)
grDevices::dev.off()

grDevices::cairo_pdf(
  filename = PDF_OUT, width = WIDTH_IN, height = HEIGHT_IN,
  family = "Arial", bg = "white", onefile = TRUE
)
print(legend_plot)
grDevices::dev.off()

svglite::svglite(
  file = SVG_OUT, width = WIDTH_IN, height = HEIGHT_IN,
  bg = "white", pointsize = 6.2
)
print(legend_plot)
grDevices::dev.off()

fwrite(legend_data[, .(
  Cluster, Label, Cluster_Color_Hex, Legend_Column,
  Legend_Row, Swatch_X, Label_X, Y
)], SOURCE_OUT)

pdf_info <- pdftools::pdf_info(PDF_OUT)
pdf_page_size <- pdftools::pdf_pagesize(PDF_OUT)
pdf_labels <- paste(pdftools::pdf_text(PDF_OUT), collapse = "\n")
required_labels <- c("Cluster", as.character(CLUSTER_LEVELS))
if (pdf_info$pages != 1L ||
    abs(pdf_page_size$width[1] - WIDTH_IN * 72) > 0.2 ||
    abs(pdf_page_size$height[1] - HEIGHT_IN * 72) > 0.2 ||
    !all(required_labels %in% unlist(strsplit(pdf_labels, "\\s+")))) {
  stop("PDF QA failed: page geometry or one or more cluster labels are missing.")
}

audit <- data.table(
  item = c(
    "source_color_table", "cluster_count", "cluster_ids", "layout",
    "target_box_px_at_72ppi", "target_pdf_page_inches", "png_pixel_dimensions",
    "color_mapping_unchanged", "pdf_page_count", "pdf_labels_present"
  ),
  value = c(
    COLOR_SOURCE, as.character(nrow(legend_data)),
    paste(legend_data$Cluster, collapse = ";"), "2 columns × 15 rows; row-major IDs",
    paste0(TARGET_WIDTH_PX, " × ", TARGET_HEIGHT_PX),
    sprintf("%.6f × %.6f", WIDTH_IN, HEIGHT_IN),
    sprintf("%d × %d", floor(WIDTH_IN * PNG_DPI), floor(HEIGHT_IN * PNG_DPI)),
    "TRUE", as.character(pdf_info$pages), "TRUE"
  )
)
fwrite(audit, AUDIT_OUT)

cat(sprintf(
  "LEGEND_RENDER_COMPLETE clusters=%d page=%.4f x %.4f in png=%d x %d px\n",
  nrow(legend_data), WIDTH_IN, HEIGHT_IN,
  floor(WIDTH_IN * PNG_DPI), floor(HEIGHT_IN * PNG_DPI)
))
