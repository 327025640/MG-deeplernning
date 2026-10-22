# Fig2F/G local-mean score contours from frozen UMAP and AddModuleScore CSVs.
# No Seurat loading, score recalculation, or UMAP/clustering recomputation.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(ragg)
  library(pdftools)
  library(RANN)
})

INPUT_ROOT <- "E:/micro_output/fig3"
OUTPUT_ROOT <- "D:/data/microglia_new/tmp/fig2fg_local_mean_contours_stage_v6"
FIGURE_OUT <- file.path(OUTPUT_ROOT, "figures")
SOURCE_OUT <- file.path(OUTPUT_ROOT, "source_data")
TABLE_OUT <- file.path(OUTPUT_ROOT, "tables")

N_GRID <- 400L
H <- 0.10
KERNEL_CUTOFF_H <- 3
MASK_NEFF_MIN <- 30
MASK_NEAREST_MAX <- 2 * H
WIDTH_IN <- 5.5
HEIGHT_IN <- 5.2
DPI <- 600L
SEED_NN <- 20260929L
CONTOUR_TARGET_INTERVALS <- 12L
CONTOUR_COLORS <- c("#4575B4", "#74ADD1", "#FFFFBF", "#FDAE61", "#D73027")
CONTOUR_LINE_COLOR <- "#666666"
CONTOUR_LINE_ALPHA <- 0.35
CONTOUR_LINE_WIDTH_MM <- 0.18

jobs <- list(
  list(
    stem = "Fig2f_microglia_M1like_AddModuleScore_UMAP",
    source = file.path(INPUT_ROOT, "source_data/Fig2f_microglia_M1like_UMAP_source.csv"),
    score = "M1_like_score",
    plot_score = "M1_like_score_plot",
    title = "M1-like local mean score"
  ),
  list(
    stem = "Fig2g_microglia_M2like_AddModuleScore_UMAP",
    source = file.path(INPUT_ROOT, "source_data/Fig2g_microglia_M2like_UMAP_source.csv"),
    score = "M2_like_score",
    plot_score = "M2_like_score_plot",
    title = "M2-like local mean score"
  )
)

dir.create(FIGURE_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_OUT, recursive = TRUE, showWarnings = FALSE)

read_cells <- function(job) {
  if (!file.exists(job$source)) stop("Missing frozen input: ", job$source)
  dt <- fread(job$source, select = c("CellID", "UMAP_1", "UMAP_2", job$score,
                                     job$plot_score, "Cluster", "Disease",
                                     "DatasetID", "SampleID_Global"))
  required <- c("CellID", "UMAP_1", "UMAP_2", job$score, job$plot_score)
  if (!all(required %in% names(dt))) stop("Frozen source schema mismatch: ", job$source)
  if (nrow(dt) != 424222L || anyNA(dt$CellID) || anyDuplicated(dt$CellID) ||
      any(!is.finite(dt$UMAP_1)) || any(!is.finite(dt$UMAP_2)) ||
      any(!is.finite(dt[[job$score]]))) {
    stop("Frozen source row count, ID uniqueness, coordinates, or raw scores failed integrity checks: ",
         job$source)
  }
  dt
}

cells_m1 <- read_cells(jobs[[1L]])
cells_m2 <- read_cells(jobs[[2L]])
if (!setequal(cells_m1$CellID, cells_m2$CellID)) stop("M1/M2 cell universes differ.")
idx_m2 <- match(cells_m2$CellID, cells_m1$CellID)
if (anyNA(idx_m2) || !identical(cells_m1$UMAP_1[idx_m2], cells_m2$UMAP_1) ||
    !identical(cells_m1$UMAP_2[idx_m2], cells_m2$UMAP_2)) {
  stop("M1/M2 frozen UMAP coordinates do not match after CellID alignment.")
}

coord_x <- range(cells_m1$UMAP_1)
coord_y <- range(cells_m1$UMAP_2)
grid_x <- seq(coord_x[1L] - KERNEL_CUTOFF_H * H,
              coord_x[2L] + KERNEL_CUTOFF_H * H, length.out = N_GRID)
grid_y <- seq(coord_y[1L] - KERNEL_CUTOFF_H * H,
              coord_y[2L] + KERNEL_CUTOFF_H * H, length.out = N_GRID)
step_x <- diff(range(grid_x)) / (N_GRID - 1L)
step_y <- diff(range(grid_y)) / (N_GRID - 1L)
n_total <- N_GRID * N_GRID

# Assign cells to the nearest common grid center and retain per-bin cell counts
# and raw-score sums. This is the only spatial binning approximation.
make_binned_sums <- function(dt, score_col) {
  ix <- pmax(1L, pmin(N_GRID, as.integer(round((dt$UMAP_1 - grid_x[1L]) / step_x)) + 1L))
  iy <- pmax(1L, pmin(N_GRID, as.integer(round((dt$UMAP_2 - grid_y[1L]) / step_y)) + 1L))
  bin <- iy + (ix - 1L) * N_GRID
  counts <- tabulate(bin, nbins = n_total)
  grouped <- data.table(bin = bin, score = dt[[score_col]])[
    , .(score_sum = sum(score)), by = bin
  ]
  sums <- numeric(n_total)
  sums[grouped$bin] <- grouped$score_sum
  list(counts = matrix(counts, nrow = N_GRID, ncol = N_GRID),
       score_sum = matrix(sums, nrow = N_GRID, ncol = N_GRID),
       ix = ix, iy = iy)
}

binned_m1 <- make_binned_sums(cells_m1, jobs[[1L]]$score)
binned_m2 <- make_binned_sums(cells_m2, jobs[[2L]]$score)
if (!identical(binned_m1$counts, binned_m2$counts)) stop("Binned common coordinate counts differ.")
count_matrix <- binned_m1$counts

radius <- KERNEL_CUTOFF_H * H
offsets <- CJ(dr = seq.int(-ceiling(radius / step_y), ceiling(radius / step_y)),
              dc = seq.int(-ceiling(radius / step_x), ceiling(radius / step_x)))
offsets[, `:=`(dy = dr * step_y, dx = dc * step_x)]
offsets[, dist2 := dx^2 + dy^2]
offsets <- offsets[dist2 <= radius^2 + .Machine$double.eps]
offsets[, `:=`(w = exp(-dist2 / (2 * H^2)), w2 = exp(-dist2 / H^2))]

convolve_grid <- function(mat, weight_col) {
  out <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  nr <- nrow(mat)
  nc <- ncol(mat)
  for (k in seq_len(nrow(offsets))) {
    dr <- offsets$dr[k]
    dc <- offsets$dc[k]
    rs <- max(1L, 1L + dr):min(nr, nr + dr)
    cs <- max(1L, 1L + dc):min(nc, nc + dc)
    rt <- rs - dr
    ct <- cs - dc
    out[rt, ct] <- out[rt, ct] + mat[rs, cs] * offsets[[weight_col]][k]
  }
  out
}

denominator <- convolve_grid(count_matrix, "w")
sumw2 <- convolve_grid(count_matrix, "w2")
effective_n <- denominator^2 / sumw2
effective_n[!is.finite(effective_n)] <- NA_real_

# Exact nearest-observed-cell support check on the shared grid (not a density
# estimate); keeps external regions and unsupported gaps blank.
grid_query <- cbind(rep(grid_x, each = N_GRID), rep(grid_y, times = N_GRID))
cell_xy <- as.matrix(cells_m1[, .(UMAP_1, UMAP_2)])
set.seed(SEED_NN)
nearest_distance <- RANN::nn2(data = cell_xy, query = grid_query, k = 1L)$nn.dists[, 1L]
nearest_matrix <- matrix(nearest_distance, nrow = N_GRID, ncol = N_GRID)
support_mask <- effective_n >= MASK_NEFF_MIN & nearest_matrix <= MASK_NEAREST_MAX & denominator > 0

local_m1 <- convolve_grid(binned_m1$score_sum, "w") / denominator
local_m2 <- convolve_grid(binned_m2$score_sum, "w") / denominator
local_m1[!support_mask] <- NA_real_
local_m2[!support_mask] <- NA_real_
finite_values <- c(local_m1[is.finite(local_m1)], local_m2[is.finite(local_m2)])
if (!length(finite_values)) stop("No supported grid cells remain under the approved mask; no contours were fabricated.")
actual_score_range <- range(finite_values)
if (!all(is.finite(actual_score_range)) || diff(actual_score_range) <= 0) {
  stop("Supported local means are constant/degenerate; no artificial contrast or levels will be created.")
}

nice_breaks <- function(values, target_n = CONTOUR_TARGET_INTERVALS) {
  rg <- range(values)
  target_step <- diff(rg) / target_n
  mag <- 10^floor(log10(target_step))
  options <- c(1, 2, 2.5, 5, 10) * mag
  step <- options[which(options >= target_step)[1L]]
  lo <- floor(rg[1L] / step) * step
  hi <- ceiling(rg[2L] / step) * step
  if (lo == hi) {
    lo <- lo - step
    hi <- hi + step
  }
  seq(lo, hi, by = step)
}
contour_breaks <- nice_breaks(finite_values)
legend_limits <- range(contour_breaks)
zero_position <- if (legend_limits[1L] <= 0 && legend_limits[2L] >= 0) {
  (0 - legend_limits[1L]) / diff(legend_limits)
} else NA_real_

# Direct point-wise checks at grid cells spanning both module-score ranges.
# These are numeric approximation checks, not an inferential test suite.
select_check_indices <- function(z) {
  valid_idx <- which(support_mask)
  vals <- as.vector(z)[valid_idx]
  probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  qv <- as.numeric(quantile(vals, probs = probs, names = FALSE, na.rm = TRUE))
  unique(vapply(qv, function(q) valid_idx[which.min(abs(vals - q))], integer(1)))
}
check_indices <- sort(unique(c(select_check_indices(local_m1), select_check_indices(local_m2))))
check_ix <- ((check_indices - 1L) %/% N_GRID) + 1L
check_iy <- ((check_indices - 1L) %% N_GRID) + 1L
check_x <- grid_x[check_ix]
check_y <- grid_y[check_iy]
direct_check_module <- function(dt, score_col, module, approx_surface) {
  rbindlist(lapply(seq_along(check_indices), function(i) {
    dx <- dt$UMAP_1 - check_x[i]
    dy <- dt$UMAP_2 - check_y[i]
    d2 <- dx^2 + dy^2
    in_kernel <- d2 <= radius^2
    w <- exp(-d2[in_kernel] / (2 * H^2))
    scores <- dt[[score_col]][in_kernel]
    direct_mean <- sum(w * scores) / sum(w)
    direct_neff <- sum(w)^2 / sum(w^2)
    approx_mean <- as.vector(approx_surface)[check_indices[i]]
    approx_neff <- as.vector(effective_n)[check_indices[i]]
    data.table(
      module = module, grid_index = check_indices[i], grid_ix = check_ix[i], grid_iy = check_iy[i],
      UMAP_1 = check_x[i], UMAP_2 = check_y[i],
      binned_mean = approx_mean, direct_mean = direct_mean,
      mean_abs_error = abs(approx_mean - direct_mean),
      binned_Neff = approx_neff, direct_Neff = direct_neff,
      Neff_abs_error = abs(approx_neff - direct_neff),
      direct_cells_within_3h = sum(in_kernel)
    )
  }))
}
direct_checks <- rbindlist(list(
  direct_check_module(cells_m1, jobs[[1L]]$score, "M1", local_m1),
  direct_check_module(cells_m2, jobs[[2L]]$score, "M2", local_m2)
))

grid_dt <- as.data.table(expand.grid(UMAP_1 = grid_x, UMAP_2 = grid_y))
grid_dt[, `:=`(
  M1_local_mean = as.vector(t(local_m1)),
  M2_local_mean = as.vector(t(local_m2)),
  Neff = as.vector(t(effective_n)),
  nearest_cell_distance = as.vector(t(nearest_matrix)),
  support_mask = as.vector(t(support_mask))
)]

palette_position <- seq(0, 1, length.out = length(CONTOUR_COLORS))
theme_contour <- theme_classic(base_family = "Arial", base_size = 9) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.35),
    axis.title = element_text(size = 8, colour = "black"),
    plot.title = element_text(face = "bold", size = 10, hjust = 0, colour = "black"),
    legend.title = element_text(size = 8, colour = "black"),
    legend.text = element_text(size = 7, colour = "black"),
    legend.key.height = unit(0.65, "cm"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    aspect.ratio = 1
  )

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

format_break <- function(x) format(signif(x, 5), trim = TRUE, scientific = FALSE)
audit_list <- list()
config_rows <- list()

for (job in jobs) {
  surface_col <- if (job$score == "M1_like_score") "M1_local_mean" else "M2_local_mean"
  plot_data <- copy(grid_dt)
  plot_data[, score_value := get(surface_col)]
  p <- ggplot() +
    geom_blank(data = cells_m1, aes(x = UMAP_1, y = UMAP_2), inherit.aes = FALSE) +
    geom_contour_filled(
      data = plot_data,
      aes(x = UMAP_1, y = UMAP_2, z = score_value, fill = after_stat(level_mid)),
      breaks = contour_breaks,
      na.rm = TRUE,
      colour = scales::alpha(CONTOUR_LINE_COLOR, CONTOUR_LINE_ALPHA),
      linewidth = CONTOUR_LINE_WIDTH_MM,
      show.legend = TRUE
    ) +
    scale_fill_gradientn(
      name = "Local mean\nAddModuleScore",
      colours = CONTOUR_COLORS,
      values = palette_position,
      limits = legend_limits,
      breaks = contour_breaks,
      labels = format_break,
      guide = guide_colorbar(frame.colour = NA, ticks.colour = "black")
    ) +
    coord_equal(xlim = coord_x, ylim = coord_y, expand = TRUE) +
    labs(title = job$title, x = "UMAP1", y = "UMAP2") +
    theme_contour

  prefix <- paste0(job$stem, "_local_mean_contour")
  full_ref <- file.path(OUTPUT_ROOT, paste0(prefix, "_reference.png"))
  contour_png <- file.path(FIGURE_OUT, paste0(prefix, "_only.png"))
  text_pdf <- file.path(FIGURE_OUT, paste0(prefix, "_text_only.pdf"))
  legend_pdf <- file.path(FIGURE_OUT, paste0(prefix, "_legend_only.pdf"))
  cell_csv <- file.path(SOURCE_OUT, paste0(prefix, "_cell_scores_source.csv"))
  grid_csv <- file.path(SOURCE_OUT, paste0(prefix, "_grid_smoothing_source.csv"))
  text_csv <- file.path(SOURCE_OUT, paste0(prefix, "_text_positions_source.csv"))
  legend_csv <- file.path(SOURCE_OUT, paste0(prefix, "_legend_source.csv"))
  direct_check_csv <- file.path(SOURCE_OUT, paste0(prefix, "_approximation_check_source.csv"))
  audit_csv <- file.path(TABLE_OUT, paste0(prefix, "_audit.csv"))
  config_csv <- file.path(TABLE_OUT, paste0(prefix, "_config.csv"))
  expected <- c(full_ref, contour_png, text_pdf, legend_pdf, cell_csv, grid_csv,
                text_csv, legend_csv, direct_check_csv, audit_csv, config_csv)
  if (any(file.exists(expected))) stop("Refusing to overwrite staged output: ", prefix)

  ragg::agg_png(full_ref, width = WIDTH_IN, height = HEIGHT_IN,
                units = "in", res = DPI, background = "white")
  grid.newpage()
  grid.draw(ggplotGrob(p))
  grid.force()
  full_grob <- grid.grab(wrap = TRUE, warn = FALSE)
  grDevices::dev.off()

  grob_inventory <- grid.ls(full_grob, grobs = TRUE, viewports = FALSE,
                            fullNames = TRUE, recursive = TRUE, print = FALSE)
  fwrite(as.data.table(unclass(as.list(grob_inventory))),
         file.path(OUTPUT_ROOT, paste0(prefix, "_grob_inventory.csv")))
  cat("LEAF_CLASSES ", prefix, ": ",
      paste(unique(leaf_classes(full_grob)), collapse = ";"), "\n", sep = "")
  contour_grob <- filter_leaf_grob(full_grob, function(x) {
    inherits(x, "pathgrob") && grepl("^geom_polygon\\.pathgrob", x$name)
  })
  contour_classes <- unique(leaf_classes(contour_grob))
  if (!length(contour_classes) || !all(contour_classes == "pathgrob")) {
    stop("Contour-only layer did not resolve to polygon grobs for ", prefix,
         "; observed: ", paste(contour_classes, collapse = ","))
  }
  text_pattern <- "(plot\\.title\\.\\.titleGrob|axis\\.title\\.[xy]\\.(bottom|left)\\.\\.titleGrob)"
  text_grob <- filter_named_subtrees(full_grob, text_pattern, drop_background = TRUE)
  text_classes <- unique(leaf_classes(text_grob))
  if (!length(text_classes) || !all(text_classes == "text")) {
    stop("Text-only layer contains non-text grobs for ", prefix, "; observed: ",
         paste(text_classes, collapse = ","))
  }
  legend_grob <- filter_named_subtrees(full_grob, "^guide-box-right", drop_background = TRUE)
  legend_classes <- unique(leaf_classes(legend_grob))
  if (!length(legend_classes) || any(legend_classes == "points")) {
    stop("Legend-only layer missing or contains point data for ", prefix)
  }

  ragg::agg_png(contour_png, width = WIDTH_IN, height = HEIGHT_IN,
                units = "in", res = DPI, background = "white")
  grid.newpage()
  grid.draw(contour_grob)
  grDevices::dev.off()

  cairo_pdf(text_pdf, width = WIDTH_IN, height = HEIGHT_IN,
            bg = "transparent", family = "Arial")
  grid.newpage()
  grid.draw(text_grob)
  grDevices::dev.off()

  cairo_pdf(legend_pdf, width = WIDTH_IN, height = HEIGHT_IN,
            bg = "transparent", family = "Arial")
  grid.newpage()
  grid.draw(legend_grob)
  grDevices::dev.off()

  cell_out <- copy(if (job$score == "M1_like_score") cells_m1 else cells_m2)
  cell_out[, grid_ix := if (job$score == "M1_like_score") binned_m1$ix else binned_m2$ix]
  cell_out[, grid_iy := if (job$score == "M1_like_score") binned_m1$iy else binned_m2$iy]
  setnames(cell_out, job$score, "raw_AddModuleScore")
  setnames(cell_out, job$plot_score, "original_1_99pct_plot_clipped_score")
  fwrite(cell_out, cell_csv)
  fwrite(grid_dt[, .(UMAP_1, UMAP_2, score_value = get(surface_col),
                     Neff, nearest_cell_distance, support_mask)], grid_csv)

  text_data <- as.data.table(pdf_data(text_pdf)[[1L]])
  text_data[, Text_Order := .I]
  setcolorder(text_data, c("Text_Order", names(text_data)[names(text_data) != "Text_Order"]))
  fwrite(text_data, text_csv)
  legend_data <- as.data.table(pdf_data(legend_pdf)[[1L]])
  legend_data[, Text_Order := .I]
  setcolorder(legend_data, c("Text_Order", names(legend_data)[names(legend_data) != "Text_Order"]))
  fwrite(legend_data, legend_csv)
  fwrite(direct_checks[module == if (job$score == "M1_like_score") "M1" else "M2"],
         direct_check_csv)

  text_content <- paste(pdf_text(text_pdf), collapse = "\n")
  legend_content <- paste(pdf_text(legend_pdf), collapse = "\n")
  if (!all(vapply(c(job$title, "UMAP1", "UMAP2"), grepl, logical(1),
                   x = text_content, fixed = TRUE))) stop("Required text missing from ", text_pdf)
  if (!grepl("Local mean", legend_content, fixed = TRUE) ||
      !grepl("AddModuleScore", legend_content, fixed = TRUE) || !nrow(legend_data)) {
    stop("Legend title/text missing from ", legend_pdf)
  }
  text_fonts <- pdf_fonts(text_pdf)
  legend_fonts <- pdf_fonts(legend_pdf)
  if (!nrow(text_fonts) || !all(text_fonts$embedded) ||
      !nrow(legend_fonts) || !all(legend_fonts$embedded)) {
    stop("Text or legend PDF font embedding failed for ", prefix)
  }

  mask_n <- sum(support_mask)
  module_values <- if (job$score == "M1_like_score") local_m1[support_mask] else local_m2[support_mask]
  input_md5 <- unname(tools::md5sum(job$source))
  audit <- data.table(
    item = c("input_csv", "input_md5", "cells", "finite_raw_scores", "umap1_min", "umap1_max",
             "umap2_min", "umap2_max", "grid_nx", "grid_ny", "grid_step_x", "grid_step_y",
             "kernel_h", "kernel_cutoff_h", "kernel_cutoff_distance", "mask_neff_min",
             "mask_nearest_max", "supported_grid_cells", "supported_grid_fraction",
             "module_supported_min", "module_supported_max", "combined_supported_min",
             "combined_supported_max", "common_legend_min", "common_legend_max",
             "contour_breaks", "zero_position_in_legend", "contour_line_color",
             "contour_line_alpha", "contour_line_width_mm", "contour_layer_grobs",
             "text_layer_grobs", "legend_layer_grobs", "text_pdf_pages", "legend_pdf_pages",
             "text_pdf_embedded_fonts", "legend_pdf_embedded_fonts", "direct_check_n",
             "direct_mean_max_abs_error", "direct_neff_max_abs_error"),
    value = c(job$source, input_md5, nrow(cell_out), sum(is.finite(cell_out$raw_AddModuleScore)),
              coord_x[1L], coord_x[2L], coord_y[1L], coord_y[2L], N_GRID, N_GRID,
              step_x, step_y, H, KERNEL_CUTOFF_H, KERNEL_CUTOFF_H * H,
              MASK_NEFF_MIN, MASK_NEAREST_MAX, mask_n, mask_n / n_total,
              min(module_values), max(module_values), actual_score_range[1L], actual_score_range[2L],
              legend_limits[1L], legend_limits[2L], paste(format_break(contour_breaks), collapse = ";"),
              zero_position, CONTOUR_LINE_COLOR, CONTOUR_LINE_ALPHA, CONTOUR_LINE_WIDTH_MM,
              paste(contour_classes, collapse = ";"),
              paste(text_classes, collapse = ";"), paste(legend_classes, collapse = ";"),
              pdf_info(text_pdf)$pages, pdf_info(legend_pdf)$pages,
              paste(text_fonts$name, collapse = ";"), paste(legend_fonts$name, collapse = ";"),
              nrow(direct_checks[module == if (job$score == "M1_like_score") "M1" else "M2"]),
              max(direct_checks[module == if (job$score == "M1_like_score") "M1" else "M2", mean_abs_error]),
              max(direct_checks[module == if (job$score == "M1_like_score") "M1" else "M2", Neff_abs_error])),
    status = rep("recorded", 40L)
  )
  fwrite(audit, audit_csv)
  audit_list[[job$stem]] <- audit
  config_rows[[job$stem]] <- data.table(
    parameter = c("method", "score_column", "score_source", "bandwidth_h", "kernel", "kernel_cutoff_h",
                  "grid_nx", "grid_ny", "grid_x_min", "grid_x_max", "grid_y_min", "grid_y_max",
                  "grid_step_x", "grid_step_y", "mask_neff_min", "mask_nearest_distance_max",
                  "common_actual_score_min", "common_actual_score_max", "common_legend_min",
                  "common_legend_max", "contour_breaks", "palette_low_to_high", "contour_line_color",
                  "contour_line_alpha", "contour_line_width_mm", "binning_approximation",
                  "replicate_inference"),
    value = c("binned Gaussian local weighted mean", job$score, job$source, H,
              "exp(-d^2/(2h^2)); truncate beyond 3h", KERNEL_CUTOFF_H,
              N_GRID, N_GRID, grid_x[1L], grid_x[N_GRID], grid_y[1L], grid_y[N_GRID],
              step_x, step_y, MASK_NEFF_MIN, MASK_NEAREST_MAX,
              actual_score_range[1L], actual_score_range[2L], legend_limits[1L],
              legend_limits[2L], paste(format_break(contour_breaks), collapse = ";"),
              paste(CONTOUR_COLORS, collapse = ";"), CONTOUR_LINE_COLOR,
              CONTOUR_LINE_ALPHA, CONTOUR_LINE_WIDTH_MM,
              "cells assigned to nearest grid center; numerator/denominator and Neff computed from binned cell counts/scores",
              "descriptive spatial summary only; cells are not independent inferential replicates")
  )
  fwrite(config_rows[[job$stem]], config_csv)

  cat("COMPLETE ", prefix, ": cells=", nrow(cell_out),
      "; supported_grid=", mask_n, "/", n_total,
      "; module_range=", paste(signif(range(module_values), 6), collapse = ","),
      "; common_range=", paste(signif(actual_score_range, 6), collapse = ","),
      "; levels=", length(contour_breaks) - 1L,
      "; direct_check_max_abs_error=", signif(max(direct_checks[
        module == if (job$score == "M1_like_score") "M1" else "M2", mean_abs_error
      ]), 5),
      "; PDF text/legend spans=", nrow(text_data), "/", nrow(legend_data), "\n", sep = "")
  rm(p, plot_data, full_grob, contour_grob, text_grob, legend_grob, cell_out)
  invisible(gc())
}

cat("ALL_MODULES_COMPLETE\n")
