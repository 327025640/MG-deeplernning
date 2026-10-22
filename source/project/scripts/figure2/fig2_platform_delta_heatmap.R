# Platform-specific MP x Disease descriptive delta heatmaps (scRNA and snRNA).
# Descriptive only: no regression, random effects, Welch test, BH or any inference;
# no donor-level analysis. Same-platform Ctrl only (DatasetID and brain region NOT matched).
# Reuses the frozen 12 MP sample usage; does not refit NMF, projection, clustering or scoring.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales) })

NMF_ROOT <- "E:/micro_output/NMF"
OUT_ROOT <- file.path(NMF_ROOT, "figure2_platform_delta_heatmap")
FIG_DIR  <- OUT_ROOT
LOG_DIR  <- file.path(OUT_ROOT, "logs")

IN_USAGE  <- file.path(NMF_ROOT, "source_data/PanelD_pseudobulk_MP_usage_raw_source.csv")
IN_META   <- file.path(NMF_ROOT, "tables/PanelD_pseudobulk_metadata.csv")
# Platform/sample mapping is taken from the verified B-panel sample source. The run04 copy
# named in the instruction was removed during the 2026-09-29 cleanup; run07 is its verified
# successor (per-sample usage identical, only S03/S04 stratum codes differ).
IN_PLAT   <- file.path(NMF_ROOT, "figure2_final_20260929_run07/tables/B_MP04_sample_usage_common_source.csv")
IN_COMP   <- file.path(NMF_ROOT, "stage2_comparability/sample_comparability.csv")
IN_MP_PAL <- file.path(NMF_ROOT, "source_data/NMFVIS_MP_palette.csv")
IN_PANEL_E<- file.path(NMF_ROOT, "source_data/PanelE_MP_disease_mean_usage_rowZ_source.csv")

INPUTS <- c(IN_USAGE, IN_META, IN_PLAT, IN_COMP, IN_MP_PAL, IN_PANEL_E)
if (!all(file.exists(INPUTS))) stop("Required frozen input missing: ", paste(INPUTS[!file.exists(INPUTS)], collapse = "; "))

# Deterministic string collation (this environment cannot set LC_COLLATE=C.UTF-8).
suppressWarnings(invisible(Sys.setlocale("LC_COLLATE", "C")))

if (dir.exists(OUT_ROOT) && length(list.files(OUT_ROOT, recursive = TRUE, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory is not empty; refusing to overwrite: ", OUT_ROOT)
}
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
sink(file.path(LOG_DIR, "fig2_platform_delta_heatmap.log"), split = TRUE)
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

# ---- frozen contracts --------------------------------------------------------
MP_LEVELS <- sprintf("MP%02d", 1:12)
# Frozen difference palette, extracted from existing project scripts
# (A panel: scripts/figure2/fig2_final_20260929_usage_panels.R;
#  NMFVIS_07_panelE_mp_disease_usage_heatmap.R uses the identical ramp).
COL_LOW <- "#2166AC"; COL_MID <- "#FFFFFF"; COL_HIGH <- "#B2182B"; COL_NA <- "#D9D9D9"
# The literal U+0394 glyph is silently DROPPED by this R/cairo setup (verified for Arial, sans,
# DejaVu Sans, Symbol, Times New Roman and Calibri in cairo_pdf), so the rendered title uses
# plotmath Delta, which keeps the glyph visually identical in both PNG and PDF. Consequence:
# pdftotext extracts the Delta as the Symbol-encoding letter "D".
COLBAR_TITLE    <- expression(Delta ~ "relative usage (percentage points)")
COLBAR_TITLE_TXT<- "\u0394 relative usage (percentage points)"

# Mandated Chinese figure legend, recorded verbatim. Chinese glyphs are silently DROPPED by
# this R/cairo setup (verified for Arial, sans, Microsoft YaHei and SimSun in both the PNG and
# cairo_pdf devices), so the figure carries an English caption and the exact Chinese wording is
# delivered in the source/record files for the manuscript legend.
LEGEND_ZH <- paste0("\u6837\u672c\u7b49\u6743\u7684\u540c\u5e73\u53f0\u75be\u75c5\u2212Ctrl\u63cf\u8ff0\u6027\u5dee\u503c\uff1b",
                    "\u672a\u6821\u6b63\u6570\u636e\u96c6\u6216\u8111\u533a\u5dee\u5f02\uff0c",
                    "\u4e0d\u80fd\u89e3\u91ca\u4e3a\u72ec\u7acb\u75be\u75c5\u6548\u5e94\u3002")
CAPTION_EN_1 <- "Sample-equal-weight same-platform disease minus Ctrl descriptive difference; dataset and brain-region differences are not adjusted, so this is not an independent disease effect."
wrap <- function(x, width) paste(strwrap(x, width = width), collapse = "\n")

# ---- inputs ------------------------------------------------------------------
usage <- fread(IN_USAGE); meta <- fread(IN_META); plat <- fread(IN_PLAT)
comp <- fread(IN_COMP);   mp_pal <- fread(IN_MP_PAL); panel_e <- fread(IN_PANEL_E)

if (!all(c("PseudobulkID", MP_LEVELS) %in% names(usage))) stop("Frozen usage table schema mismatch")
if (!all(c("MP", "Color_Hex") %in% names(mp_pal)) || !identical(mp_pal$MP, MP_LEVELS)) stop("Frozen MP palette/order mismatch")
if (!all(c("SampleID_Global", "platform", "platform_source", "Disease", "DatasetID") %in% names(plat))) stop("Platform source schema mismatch")
if (!all(c("PseudobulkID", "SampleID_Global", "Disease", "n_cells", "included_ge20") %in% names(meta))) stop("Frozen metadata schema mismatch")

# Row order follows the current formal PanelE: its source table is written in the frozen
# MP palette order, and the rendered PanelE figure shows the same MP01..MP12 order.
e_mps <- panel_e[[1]]
if (!identical(e_mps, MP_LEVELS)) stop("PanelE MP row order is not the frozen MP01..MP12 order")
# Column order: PanelE disease order, Ctrl dropped (it is the baseline, its delta is 0 by construction).
e_dis <- setdiff(names(panel_e), "MP")
if (!identical(e_dis[1], "Ctrl")) stop("PanelE first disease column is not Ctrl; refusing to guess the frozen order")
DIS_ORDER <- e_dis[-1]
cat("PanelE MP row order:", paste(e_mps, collapse = ", "), "\n")
cat("PanelE disease column order (Ctrl dropped):", paste(DIS_ORDER, collapse = ", "), "\n")

# ---- eligible samples and platform mapping -----------------------------------
elig <- meta[included_ge20 == TRUE]
if (nrow(elig) != 432L || anyDuplicated(elig$SampleID_Global)) stop("Eligible sample universe is not the frozen 432")
if (any(elig$n_cells < 20L)) stop("A sample has fewer than 20 cells")
if (anyDuplicated(plat$SampleID_Global)) stop("Duplicate SampleID_Global in the platform source")

d <- merge(elig[, .(PseudobulkID, SampleID_Global, DatasetID, Disease, n_cells)],
           usage, by = "PseudobulkID", all.x = FALSE)
d <- merge(d, plat[, .(SampleID_Global, platform, platform_source)], by = "SampleID_Global", all.x = TRUE)
if (nrow(d) != 432L || anyDuplicated(d$SampleID_Global) || anyNA(d$platform)) stop("Platform join did not preserve the frozen 432-sample universe")
if (!all(d$platform %in% c("scRNA", "snRNA"))) stop("Unresolved platform value present")

# cross-check the platform mapping against the frozen comparability table
cs <- unique(comp[, .(SampleID_Global, sc_or_sn)])
cs[, sc_or_sn := fifelse(sc_or_sn == "scRNA-seq", "scRNA", sc_or_sn)]
chk <- merge(d[, .(SampleID_Global, platform, DatasetID)], cs, by = "SampleID_Global", all.x = TRUE)
tagore <- chk$DatasetID == "TAGORE_2025_NSCLC_BM"
agree <- chk[!tagore & !is.na(sc_or_sn), platform == sc_or_sn]
cat("Platform cross-check vs sample_comparability (TAGORE excluded): ",
    sum(agree), "/", length(agree), " agree\n")
cat("TAGORE samples resolved from local series matrix as:",
    paste(unique(d$platform[d$DatasetID == "TAGORE_2025_NSCLC_BM"]), collapse = ", "),
    "| n =", sum(tagore), "\n")
cat("platform_source record:", paste(unique(d[platform_source != "sample_comparability.sc_or_sn", platform_source]), collapse = "; "), "\n")

mat <- as.matrix(d[, ..MP_LEVELS])
if (any(!is.finite(mat)) || min(mat) < 0 || max(mat) > 1) stop("Usage contains NA/non-finite or falls outside [0,1]")
if (max(abs(rowSums(mat) - 1)) > 1e-6) stop("Frozen sample usage rows do not sum to 1")
cat("Usage rows sum to 1 (max abs deviation", max(abs(rowSums(mat) - 1)), ")\n")

# ---- per-platform disease minus same-platform Ctrl ---------------------------
PLATFORMS <- c("scRNA", "snRNA")
res <- list(); k <- 0L
for (pf in PLATFORMS) {
  dp <- d[platform == pf]
  ctrl <- dp[Disease == "Ctrl"]
  cat("\n[", pf, "] samples:", nrow(dp), "| Ctrl samples:", nrow(ctrl),
      "| Ctrl datasets:", length(unique(ctrl$DatasetID)), "\n", sep = "")
  if (nrow(ctrl) < 1L) stop("Platform ", pf, " has no Ctrl sample; delta is undefined")
  ctrl_mean <- colMeans(as.matrix(ctrl[, ..MP_LEVELS]))
  for (mp in MP_LEVELS) for (dz in DIS_ORDER) {
    k <- k + 1L
    sub <- dp[Disease == dz]
    reason <- NA_character_
    if (nrow(sub) == 0L) {
      reason <- "no_sample_of_this_disease_on_this_platform"
      md <- NA_real_; delta <- NA_real_
    } else {
      md <- mean(sub[[mp]]); delta <- 100 * (md - ctrl_mean[[mp]])
    }
    res[[k]] <- data.table(platform = pf, MP = mp, Disease = dz,
                           n_disease = nrow(sub), n_ctrl = nrow(ctrl),
                           mean_disease = md, mean_ctrl = ctrl_mean[[mp]],
                           delta_percentage_points = delta, missing_reason = reason)
  }
}
res <- rbindlist(res)
if (nrow(res) != length(PLATFORMS) * length(MP_LEVELS) * length(DIS_ORDER)) stop("Unexpected number of cells")

# Per-platform independent colour scale and column set.
# User decision 2026-09-29: supersedes the shared symmetric scale and the grey NA columns of the
# first run. Each platform now uses its own symmetric limit covering all of its finite deltas, and
# diseases with no sample on that platform are dropped from the plot instead of drawn as grey NA.
PLAT_PLAN <- list()
for (pf in PLATFORMS) {
  z <- res[platform == pf]
  keep_dis <- DIS_ORDER[DIS_ORDER %in% z[n_disease > 0L, Disease]]
  lim_pf <- max(abs(z[is.finite(delta_percentage_points), delta_percentage_points]))
  if (any(abs(z$delta_percentage_points) > lim_pf + 1e-9, na.rm = TRUE)) stop("Delta outside the per-platform limit for ", pf)
  PLAT_PLAN[[pf]] <- list(dis = keep_dis, lim = lim_pf)
  cat(sprintf("\n[%s] plotted disease columns: %d -> %s\n", pf, length(keep_dis), paste(keep_dis, collapse = ", ")))
  cat(sprintf("[%s] independent symmetric colour limits: -%.6f to +%.6f\n", pf, lim_pf, lim_pf))
}
res[, plotted_column := n_disease > 0L]
cat("\nFinite delta cells:", sum(is.finite(res$delta_percentage_points)), "of", nrow(res),
    "; dropped disease columns:", sum(!res$plotted_column) / length(MP_LEVELS), "\n")

# ---- plot ---------------------------------------------------------------------
# Square tiles (coord_fixed) mean the 12-row panel is square, so the canvas has to be tall
# enough for the panel to also fill the requested width. At 160 x 100 mm the square panel
# occupied only ~76 mm and left ~55 mm of dead space on the right, so the height is set to
# the value that makes the panel fill the width. Both figures share the same canvas and the
# same panel height, so the tile side length is identical in the two panels.
FIG_W <- 160
FIG_H <- 128
make_plot <- function(dt, pf, n_ctrl, n_ctrl_ds, dis_keep, lim) {
  z <- copy(dt[platform == pf & Disease %in% dis_keep])
  z[, MP := factor(MP, levels = rev(MP_LEVELS))]
  z[, Disease := factor(Disease, levels = dis_keep)]
  cap <- paste(
    wrap(CAPTION_EN_1, 118),
    wrap(sprintf("%s baseline: %d Ctrl sample%s across %d dataset%s; Ctrl is used only as the subtraction baseline and is not drawn as a column. Diseases with no sample on this platform are omitted, not shown as zero.",
                 pf, n_ctrl, ifelse(n_ctrl == 1L, "", "s"), n_ctrl_ds, ifelse(n_ctrl_ds == 1L, "", "s")), 118),
    sep = "\n")
  ggplot(z, aes(x = Disease, y = MP, fill = delta_percentage_points)) +
    geom_tile(colour = "white", linewidth = .5) +
    scale_fill_gradient2(low = COL_LOW, mid = COL_MID, high = COL_HIGH, midpoint = 0,
                         limits = c(-lim, lim), na.value = COL_NA, name = COLBAR_TITLE,
                         breaks = scales::breaks_pretty(n = 5)) +
    scale_x_discrete(limits = dis_keep, drop = FALSE, expand = c(0, 0)) +
    scale_y_discrete(limits = rev(MP_LEVELS), drop = FALSE, expand = c(0, 0)) +
    # square tiles: one discrete category = one data unit on both axes, so with ratio = 1
    # every cell is square. Both figures are rendered at the same canvas size and both are
    # height-limited, so the tile side length is identical in the two panels.
    coord_fixed(ratio = 1, expand = FALSE) +
    labs(title = pf, x = NULL, y = NULL, caption = cap) +
    theme_minimal(base_family = "Arial", base_size = 7) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, vjust = 1, size = 6.4),
          axis.text.y = element_text(colour = "black", size = 6.4),
          plot.title = element_text(colour = "black", face = "bold", size = 8.5, hjust = 0),
          legend.title = element_text(colour = "black", size = 6.2),
          legend.text = element_text(colour = "black", size = 6),
          legend.key.height = unit(12, "mm"), legend.key.width = unit(2.6, "mm"),
          plot.caption = element_text(colour = "black", size = 4.6, hjust = 0, lineheight = 1.15),
          plot.margin = margin(4, 4, 3, 3))
}

for (pf in PLATFORMS) {
  nc <- unique(res[platform == pf, n_ctrl])
  nds <- unique(d[platform == pf & Disease == "Ctrl", .N, by = DatasetID])[, .N]
  p <- make_plot(res, pf, nc, nds, PLAT_PLAN[[pf]]$dis, PLAT_PLAN[[pf]]$lim)
  # measure on an explicit device: with no device open, ggplotGrob falls back to the default
  # pdf device, which cannot resolve the Arial/CID fonts used here
  grDevices::png(tempfile(fileext = ".png"), width = FIG_W, height = FIG_H, units = "mm", res = 600, bg = "white")
  gt <- ggplotGrob(p); print(gt)
  # The panel track is a null (flexible) unit, so it cannot be converted directly. The height
  # track resolves correctly as canvas minus the absolute tracks, and coord_fixed(ratio = 1)
  # forces one x-unit and one y-unit to the same physical length, so the tile side follows from
  # the height and the panel width is ncol x tile. (Verifying this on the exported PNG by
  # measuring the white tile borders gives 8.50 mm for scRNA and 8.47 mm for snRNA.)
  ph <- FIG_H - sum(as.numeric(grid::convertHeight(gt$heights, "mm", valueOnly = TRUE)))
  grDevices::dev.off()
  ncol_p <- length(PLAT_PLAN[[pf]]$dis); nrow_p <- length(MP_LEVELS)
  tile <- ph / nrow_p; pw <- tile * ncol_p
  cat(sprintf("[%s] panel %.1f x %.1f mm -> tile %.3f mm square (tile w/h = 1 by coord_fixed; ncol=%d, nrow=%d)\n",
              pf, pw, ph, tile, ncol_p, nrow_p))
  base <- file.path(FIG_DIR, paste0("MP_Disease_delta_", pf))
  withCallingHandlers({
    ggsave(paste0(base, ".png"), p, width = FIG_W, height = FIG_H, units = "mm", dpi = 600, bg = "white")
    ggsave(paste0(base, ".pdf"), p, width = FIG_W, height = FIG_H, units = "mm", device = cairo_pdf, bg = "white")
  }, warning = function(w) { cat("RENDER WARNING:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
}

# ---- source and records -------------------------------------------------------
lim_by_pf <- setNames(vapply(PLATFORMS, function(p) PLAT_PLAN[[p]]$lim, numeric(1)), PLATFORMS)
out <- res[, .(platform, MP, Disease, n_disease, n_ctrl, mean_disease, mean_ctrl,
               delta_percentage_points, missing_reason, plotted_column)]
out[, `:=`(delta_formula = "100 * (mean_disease - mean_ctrl) within the same platform",
           ctrl_definition = "all Ctrl samples of the same platform (DatasetID and brain region NOT matched)",
           row_order = match(MP, MP_LEVELS), column_order_in_panel_e = match(Disease, DIS_ORDER),
           colour_low = COL_LOW, colour_mid = COL_MID, colour_high = COL_HIGH,
           colour_na_defined_but_unused = COL_NA,
           platform_colour_min = -lim_by_pf[platform], platform_colour_max = lim_by_pf[platform],
           caption_figure_en = CAPTION_EN_1,
           legend_zh = LEGEND_ZH)]
setorder(out, platform, row_order, column_order_in_panel_e)
for (pf in PLATFORMS) {
  fwrite(out[platform == pf], file.path(FIG_DIR, paste0("MP_Disease_delta_", pf, "_source.csv")))
}

record <- data.table(
  item = c("difference_low", "difference_mid", "difference_high", "na_grey_status", "colourbar_title",
           "colour_scale_rule", "scRNA_colour_min", "scRNA_colour_max",
           "snRNA_colour_min", "snRNA_colour_max",
           "scRNA_plotted_columns", "snRNA_plotted_columns",
           "dropped_columns_rule", "mp_row_order", "mp_row_order_source",
           "disease_column_order", "disease_column_order_source", "ctrl_definition",
           "delta_formula", "platform_source_file", "colour_palette_source",
           "legend_zh_verbatim", "caption_in_figure_en"),
  value = c(COL_LOW, COL_MID, COL_HIGH,
            paste0(COL_NA, " (still defined in the source; no grey NA cell remains because NA columns are dropped)"),
            COLBAR_TITLE_TXT,
            "independent per-platform symmetric scale, each covering that platform's finite deltas without truncation",
            -lim_by_pf[["scRNA"]], lim_by_pf[["scRNA"]],
            -lim_by_pf[["snRNA"]], lim_by_pf[["snRNA"]],
            paste(PLAT_PLAN[["scRNA"]]$dis, collapse = " > "),
            paste(PLAT_PLAN[["snRNA"]]$dis, collapse = " > "),
            "diseases with zero samples on that platform are omitted from the plot (user decision 2026-09-29)",
            paste(MP_LEVELS, collapse = " > "),
            "E:/micro_output/NMF/source_data/PanelE_MP_disease_mean_usage_rowZ_source.csv (row order) and the rendered PanelE figure",
            paste(DIS_ORDER, collapse = " > "),
            "PanelE disease column order with Ctrl dropped; each panel keeps this relative order among its plotted diseases",
            "all Ctrl samples of the same platform; DatasetID and brain region not matched",
            "100 * (platform disease mean usage - platform Ctrl mean usage)",
            "figure2_final_20260929_run07/tables/B_MP04_sample_usage_common_source.csv (platform, verified successor of the run04 copy)",
            "scripts/figure2/fig2_final_20260929_usage_panels.R (A panel) and NMFVIS_07_panelE_mp_disease_usage_heatmap.R",
            LEGEND_ZH, CAPTION_EN_1))
fwrite(record, file.path(FIG_DIR, "delta_heatmap_color_and_order_record.csv"))

cat("\n=== per-platform summary ===\n")
print(res[, .(finite = sum(is.finite(delta_percentage_points)),
              na = sum(!is.finite(delta_percentage_points)),
              min = min(delta_percentage_points, na.rm = TRUE),
              max = max(delta_percentage_points, na.rm = TRUE)), by = platform])
cat("\nNA columns per platform:\n")
print(res[!is.finite(delta_percentage_points), .(Disease = paste(unique(Disease), collapse = ", ")), by = platform])
cat("\nscRNA Ctrl disclosure: n_ctrl =", unique(res[platform == "scRNA", n_ctrl]),
    "; snRNA n_ctrl =", unique(res[platform == "snRNA", n_ctrl]), "\n")
cat("Outputs:", OUT_ROOT, "\n")
graphics.off()
cat("Graphics devices after explicit close:", length(dev.list()), "\n")
cat("End:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
sink()
