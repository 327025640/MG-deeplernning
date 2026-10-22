suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})
options(stringsAsFactors = FALSE)

ROOT <- "D:/data/microglia_new/analysis_scripts/runtime_io_aliases/fig3"
FIG_DIR <- file.path(ROOT, "figures")
SCOPE_FILE <- file.path(ROOT, "tables/cluster_scope.csv")
FUN_FILE <- file.path(ROOT, "source_data/fig4_cluster24_FUN8_source.csv")
OUTPUT_DIR <- Sys.getenv("FIG4_B_OUTPUT_DIR", unset = FIG_DIR)
PNG_OUT <- file.path(OUTPUT_DIR, "Fig4_B_cluster24_FUN8.png")
PDF_OUT <- file.path(OUTPUT_DIR, "Fig4_B_cluster24_FUN8.pdf")

required_files <- c(SCOPE_FILE, FUN_FILE)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing existing Fig4B inputs/outputs: ", paste(missing_files, collapse = "; "))
if (!dir.exists(OUTPUT_DIR) && !dir.create(OUTPUT_DIR, recursive = TRUE)) stop("Cannot create output directory: ", OUTPUT_DIR)

scope <- fread(SCOPE_FILE)
fun <- fread(FUN_FILE, na.strings = c("", "NA"))
if (nrow(scope) != 24L || anyDuplicated(scope$Cluster) || sum(scope$n_cells_total) != 400026L) {
  stop("Cluster scope is not the verified clean24 set")
}
scope_order <- c("Shared6", "Restricted4", "pan14")
order_by_frozen_row <- function(cluster_ids) {
  scope[Cluster %in% cluster_ids][order(row_order), as.character(Cluster)]
}
scope_map <- list(
  Shared6 = order_by_frozen_row(c("C11", "C14", "C16", "C18", "C21", "C22")),
  Restricted4 = order_by_frozen_row(c("C23", "C25", "C26", "C29")),
  pan14 = scope[analysis_scope == "pan14"][order(row_order), Cluster]
)
ordered_clusters <- unlist(scope_map[scope_order], use.names = FALSE)
scope_colors <- c(Shared6 = "#44718C", Restricted4 = "#B87555", pan14 = "#CFCFCF")
scope_lookup <- rbindlist(lapply(scope_order, function(group) {
  data.table(Cluster = scope_map[[group]], scope_group = group)
}))
scope_lookup[, Cluster := factor(Cluster, levels = rev(ordered_clusters))]
scope_lookup[, scope_group := factor(scope_group, levels = scope_order)]
scope_legend <- scope_lookup[, .SD[1], by = scope_group]

fun_labels <- c(
  F8_Surveillance = "Surveillance",
  F8_Neuroprotection = "Neuroprotective\nsignaling",
  F8_Phagocytosis = "Phagocytosis",
  F8_Inflammation = "Inflammation",
  F8_Cytokine_production = "Cytokine\nproduction",
  F8_Antigen_presentation = "Antigen\npresentation",
  F8_IFN_signature = "Interferon\nresponse",
  F8_Proliferation = "Proliferation"
)
if (nrow(fun) != 192L || anyDuplicated(fun[, .(Cluster, Function)]) ||
    !setequal(unique(fun$Cluster), ordered_clusters) || !setequal(unique(fun$Function), names(fun_labels))) {
  stop("Panel B source does not match the complete frozen 24 x 8 FUN8 matrix")
}
fun[, Cluster := factor(Cluster, levels = rev(ordered_clusters))]
fun[, feature_index := match(Function, names(fun_labels))]
if (sum(fun$z_display_clipped) != 1L) stop("FUN8 clipping marker count differs from its source")

base_theme <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", colour = "#111111"),
    axis.text = element_text(colour = "#111111", size = 7),
    axis.title = element_text(colour = "#111111", size = 8),
    plot.title = element_text(face = "bold", size = 10, hjust = 0),
    plot.subtitle = element_text(size = 7, colour = "#333333", hjust = 0),
    plot.caption = element_text(size = 6.1, colour = "#333333", hjust = 0, lineheight = 1.05),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6.5),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.height = grid::unit(3.4, "mm"),
    legend.key.width = grid::unit(13, "mm"),
    axis.ticks = element_line(linewidth = 0.25, colour = "#555555"),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(5, 5, 4, 4, unit = "mm")
  )
wrap_caption <- function(x) paste(strwrap(x, width = 118L), collapse = "\n")

plot_b <- ggplot(fun, aes(x = feature_index, y = Cluster, fill = z_display)) +
  geom_tile(width = 0.94, height = 0.94, colour = "white", linewidth = 0.18) +
  geom_point(data = fun[z_display_clipped == TRUE], aes(x = feature_index + 0.33, y = Cluster),
             inherit.aes = FALSE, shape = 24, size = 1.05, fill = "#202020",
             colour = "white", stroke = 0.16) +
  scale_fill_gradient2(name = "Z across 24 clusters", low = "#3B6FB6", mid = "#F7F7F7",
                       high = "#C44E52", midpoint = 0, limits = c(-3, 3),
                       oob = scales::squish, breaks = c(-3, -1.5, 0, 1.5, 3),
                       na.value = "#D9D9D9") +
  scale_colour_manual(name = "Cluster group", values = scope_colors,
                      breaks = scope_order, drop = FALSE) +
  guides(colour = guide_legend(order = 2,
                               override.aes = list(shape = 15, size = 4, alpha = 1)),
         fill = guide_colourbar(order = 1)) +
  scale_x_continuous(breaks = seq_along(fun_labels), labels = unname(fun_labels),
                     limits = c(0, length(fun_labels) + 0.5), expand = c(0, 0)) +
  scale_y_discrete(limits = rev(ordered_clusters), expand = expansion(add = 0.02)) +
  coord_fixed(ratio = 1, clip = "off") +
  labs(title = "B. FUN8 profiles across all 24 microglial clusters",
       subtitle = "Cell-equal cluster means; each feature is standardized across the 24 clusters",
       x = NULL, y = NULL,
       caption = wrap_caption("Z scores are relative signature-score summaries, not experimental evidence of functional activation. FUN8 cluster means are descriptive and cell-equal. Triangles mark unclipped Z values outside the -3 to +3 display range.")) +
  base_theme +
  theme(axis.text.x = element_text(size = 6.5, lineheight = 0.95),
        axis.text.y = element_text(size = 6.8))

for (group in scope_order) {
  gd <- scope_lookup[scope_group == group]
  plot_b <- plot_b + geom_tile(data = gd, aes(x = 0.25, y = Cluster), inherit.aes = FALSE,
                               width = 0.38, height = 0.9, fill = scope_colors[group])
}
plot_b <- plot_b +
  geom_point(data = scope_legend,
             aes(x = 0.25, y = Cluster, colour = scope_group),
             inherit.aes = FALSE, shape = 15, size = 0.1, alpha = 0,
             show.legend = TRUE)
for (boundary in c(6L, 10L)) {
  plot_b <- plot_b + geom_hline(yintercept = length(ordered_clusters) - boundary + 0.5,
                                colour = "#555555", linewidth = 0.3)
}

png_device <- function(filename, width, height, bg, dpi = 600, ...) {
  grDevices::png(filename = filename, width = round(width * dpi), height = round(height * dpi),
                 bg = bg, res = dpi, family = "Arial", type = "cairo", ...)
}
pdf_device <- function(filename, width, height, bg, ...) {
  grDevices::cairo_pdf(filename = filename, width = width, height = height,
                       family = "Arial", bg = bg, ...)
}
ggsave(PNG_OUT, plot_b, device = png_device, width = 180, height = 490,
       units = "mm", dpi = 600, bg = "white", limitsize = FALSE)
ggsave(PDF_OUT, plot_b, device = pdf_device, width = 180 / 25.4, height = 490 / 25.4,
       units = "in", bg = "white", limitsize = FALSE)
cat("Updated Fig4B with square cells and an external cluster-group legend; output size 180 x 490 mm.\n")
