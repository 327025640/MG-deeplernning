# PanelD variant: pseudobulk x MP row-Z heatmap, transposed, split by sequencing platform.
# MP on rows (12), pseudobulk samples on columns (432), columns split into a left scRNA block
# and a right snRNA block with one shared cell width.
# Removes the Dataset annotation bar and the n_cells bar of the original PanelD; adds a Platform
# bar. Reuses the FROZEN row-Z matrix and metadata: no usage, Z score, clustering or projection is
# recomputed.
suppressPackageStartupMessages({
  library(data.table); library(ComplexHeatmap); library(circlize); library(grid)
})

NMF_ROOT <- "E:/micro_output/NMF"
OUT_ROOT <- file.path(NMF_ROOT, "figure2_panelD_platform_split")
FIG_DIR  <- OUT_ROOT
LOG_DIR  <- file.path(OUT_ROOT, "logs")

IN_ROWZ <- file.path(NMF_ROOT, "source_data/PanelD_pseudobulk_MP_usage_rowZ_source.csv")
IN_META <- file.path(NMF_ROOT, "tables/PanelD_pseudobulk_metadata.csv")
IN_PLAT <- file.path(NMF_ROOT, "figure2_final_20260929_run07/tables/B_MP04_sample_usage_common_source.csv")
IN_COMP <- file.path(NMF_ROOT, "stage2_comparability/sample_comparability.csv")
IN_DPAL <- file.path(NMF_ROOT, "tables/Figure2_12MP_Disease17_palette.csv")
IN_MPAL <- file.path(NMF_ROOT, "source_data/NMFVIS_MP_palette.csv")
INPUTS <- c(IN_ROWZ, IN_META, IN_PLAT, IN_COMP, IN_DPAL, IN_MPAL)
if (!all(file.exists(INPUTS))) stop("Required frozen input missing: ", paste(INPUTS[!file.exists(INPUTS)], collapse = "; "))

suppressWarnings(invisible(Sys.setlocale("LC_COLLATE", "C")))
if (dir.exists(OUT_ROOT) && length(list.files(OUT_ROOT, recursive = TRUE, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory is not empty; refusing to overwrite: ", OUT_ROOT)
}
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
sink(file.path(LOG_DIR, "fig2_panelD_platform_split.log"), split = TRUE)
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

MP_LEVELS <- sprintf("MP%02d", 1:12)
COL_LOW <- "#2166AC"; COL_MID <- "#FFFFFF"; COL_HIGH <- "#B2182B"   # frozen row-Z ramp (unchanged)
PLAT_LEVELS <- c("scRNA", "snRNA")
PLAT_COL <- c(scRNA = "#969696", snRNA = "#4D4D4D")                 # frozen platform palette (as in S1)

# ---- inputs ------------------------------------------------------------------
rz   <- fread(IN_ROWZ); meta <- fread(IN_META)
plat <- fread(IN_PLAT); comp <- fread(IN_COMP); dpal <- fread(IN_DPAL); mpal <- fread(IN_MPAL)

if (!all(c("PseudobulkID", MP_LEVELS) %in% names(rz))) stop("Frozen row-Z table schema mismatch")
if (!all(c("PseudobulkID", "SampleID_Global", "DatasetID", "Disease", "Disease_Abbrev_Analysis", "n_cells") %in% names(meta))) stop("Frozen metadata schema mismatch")
if (!all(c("MP", "Color_Hex") %in% names(mpal)) || !identical(mpal$MP, MP_LEVELS)) stop("Frozen MP palette/order mismatch")
if (nrow(rz) != 432L || anyDuplicated(rz$PseudobulkID)) stop("Row-Z table is not the frozen 432-pseudobulk universe")
if (max(abs(colMeans(rz[, ..MP_LEVELS]))) > 1e-8 || max(abs(apply(rz[, ..MP_LEVELS], 2, sd) - 1)) > 1e-3) {
  stop("Frozen row-Z columns are not mean 0 / sd 1; refusing to treat them as the frozen Z")
}

# platform mapping, cross-checked against the frozen comparability table
cs <- unique(comp[, .(SampleID_Global, sc_or_sn)])
cs[, sc_or_sn := fifelse(sc_or_sn == "scRNA-seq", "scRNA", sc_or_sn)]
if (!"included_ge20" %in% names(meta)) stop("Frozen metadata lacks included_ge20")
meta <- meta[included_ge20 == TRUE]          # 453 -> the frozen 432 eligible pseudobulks
if (nrow(meta) != 432L) stop("Eligible pseudobulk set is not 432")
# Disease, not Disease_Abbrev_Analysis: the canonical column carries the 17 frozen classes
# (pGBM + rGBM collapse into GBM, Anaplastic_astrocytoma collapses into LGG). The alias column
# would draw 19 categories and split those tumours into separate legend entries.
d <- merge(meta[, .(PseudobulkID, SampleID_Global, DatasetID, Disease,
                    Disease_Abbrev_Analysis, n_cells)],
           plat[, .(SampleID_Global, platform, platform_source)], by = "SampleID_Global", all.x = TRUE)
d <- merge(d, cs, by = "SampleID_Global", all.x = TRUE)
if (nrow(d) != 432L || anyNA(d$platform)) stop("Platform join did not preserve the 432 pseudobulks")
if (!all(d$platform %in% PLAT_LEVELS)) stop("Unresolved platform value")
tagore <- d$DatasetID == "TAGORE_2025_NSCLC_BM"
agree <- d[!tagore & !is.na(sc_or_sn), platform == sc_or_sn]
cat("Platform cross-check vs sample_comparability (TAGORE excluded):", sum(agree), "/", length(agree), "agree\n")
cat("TAGORE resolved from local series matrix as snRNA, n =", sum(tagore), ", source =",
    paste(unique(d[tagore, platform_source]), collapse = "; "), "\n")

# ---- column order: platform block, then dataset, disease, sample (same convention as S1) -----
d[, platform := factor(platform, levels = PLAT_LEVELS)]
setorder(d, platform, DatasetID, Disease, SampleID_Global)
if (anyDuplicated(d$PseudobulkID)) stop("Duplicate pseudobulk after ordering")
Zs <- as.matrix(rz[match(d$PseudobulkID, PseudobulkID), ..MP_LEVELS])
rownames(Zs) <- d$PseudobulkID
stopifnot(all(is.finite(Zs)))
Zt <- t(Zs)                       # MP x pseudobulk
rownames(Zt) <- MP_LEVELS
cat("\nMatrix after transpose:", nrow(Zt), "MP x", ncol(Zt), "pseudobulks\n")
print(table(d$platform))

# ---- colours -----------------------------------------------------------------
if (!all(c("Disease", "Color_Hex") %in% names(dpal))) stop("Frozen Disease palette schema mismatch")
d_col <- setNames(dpal$Color_Hex, dpal$Disease)
dis_used <- sort(unique(as.character(d$Disease)))
miss <- setdiff(dis_used, names(d_col))
if (length(miss)) stop("Disease label without a frozen colour: ", paste(miss, collapse = ", "))
cat("\nDisease classes drawn:", length(dis_used), "->", paste(dis_used, collapse = ", "), "\n")
# report the alias collapse explicitly so the legend change is auditable
ali <- unique(d[, .(Disease, Disease_Abbrev_Analysis)])
cat("Alias -> canonical collapse:\n"); print(ali[order(Disease, Disease_Abbrev_Analysis)])
# ---- heatmap -----------------------------------------------------------------
n_by_pf <- as.integer(table(d$platform)[PLAT_LEVELS])
split_lab <- paste0(PLAT_LEVELS, "  (n = ", n_by_pf, ")")
split_f <- factor(split_lab[match(as.character(d$platform), PLAT_LEVELS)], levels = split_lab)
ha <- HeatmapAnnotation(
  Platform = as.character(d$platform),
  Disease  = as.character(d$Disease),
  col = list(Platform = PLAT_COL, Disease = d_col[dis_used]),
  show_legend = c(TRUE, TRUE),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 8))

ht <- Heatmap(
  Zt, name = "row Z",
  col = colorRamp2(c(-2, 0, 2), c(COL_LOW, COL_MID, COL_HIGH)),
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8, fontface = "bold",
                      col = setNames(mpal$Color_Hex, mpal$MP)[MP_LEVELS]),
  show_column_names = FALSE,
  column_split = split_f,
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  top_annotation = ha,
  use_raster = FALSE,
  border = TRUE, rect_gp = gpar(col = NA))

W_IN <- 14; H_IN <- 6
# Font handling: do NOT set gpar(fontfamily=) here. Requesting a family that a device cannot
# resolve makes grid fall back to the PostScript font database for text metrics and emit
# "font family 'Arial' not found in PostScript font database" once per text element. Both output
# devices already give Arial by default:
#   - ragg::agg_png -> systemfonts maps the default sans to C:/Windows/Fonts/arial.ttf
#   - cairo_pdf(family = "Arial") -> fontconfig resolves Arial and embeds ArialMT / Arial-BoldMT
#   (the grDevices Windows png device cannot be used: windowsFonts() only exposes serif/sans/mono)
if (!requireNamespace("ragg", quietly = TRUE)) stop("ragg is required for the PNG output")
ragg::agg_png(file.path(FIG_DIR, "PanelD_pseudobulk_MP_usage_heatmap_platform_split.png"),
              width = W_IN, height = H_IN, units = "in", res = 600, background = "white")
draw(ht)
dev.off()
cairo_pdf(file.path(FIG_DIR, "PanelD_pseudobulk_MP_usage_heatmap_platform_split.pdf"),
          width = W_IN, height = H_IN, family = "Arial")
draw(ht)
dev.off()

# ---- source and record -------------------------------------------------------
out <- data.table(PseudobulkID = d$PseudobulkID, SampleID_Global = d$SampleID_Global,
                  DatasetID = d$DatasetID, platform = as.character(d$platform),
                  Disease = as.character(d$Disease),
                  Disease_alias_original = as.character(d$Disease_Abbrev_Analysis),
                  n_cells = d$n_cells, plot_column_order = seq_len(nrow(d)))
out <- cbind(out, as.data.table(t(Zt)[, MP_LEVELS]))
out[, colour_scale := "frozen PanelD row Z, taken unchanged from PanelD_pseudobulk_MP_usage_rowZ_source.csv (mean 0 / sd 1 across all 432 pseudobulks)"]
fwrite(out, file.path(FIG_DIR, "PanelD_pseudobulk_MP_usage_heatmap_platform_split_source.csv"))

record <- data.table(
  item = c("orientation", "row_order", "column_order_rule", "platform_split",
           "removed_annotations", "kept_annotations", "added_annotation",
           "disease_mapping_source", "disease_alias_collapse",
           "rowZ_definition", "rowZ_recomputed",
           "rowZ_low", "rowZ_mid", "rowZ_high",
           "platform_scRNA_colour", "platform_snRNA_colour",
           "platform_mapping_source", "n_scRNA", "n_snRNA", "figure_size_in"),
  value = c("transposed: 12 MP rows x 432 pseudobulk columns",
            paste(MP_LEVELS, collapse = " > "),
            "platform block > DatasetID > Disease > SampleID_Global (same convention as S1)",
            "left = scRNA, right = snRNA, one shared cell width, separated by a gap with a block title",
            "Dataset bar; n_cells barplot (and their legends)",
            "Disease bar with legend",
            "Platform bar (SN/SC) with legend",
            paste0("canonical Disease column of the frozen PanelD metadata (", length(dis_used),
                   " classes), coloured with tables/Figure2_12MP_Disease17_palette.csv; supersedes the earlier build that used the 19-label Disease_Abbrev_Analysis alias column"),
            "pGBM + rGBM -> GBM ; Anaplastic_astrocytoma + LGG -> LGG",
            "row Z per MP across all 432 pseudobulks (frozen file)",
            "no",
            COL_LOW, COL_MID, COL_HIGH,
            PLAT_COL[["scRNA"]], PLAT_COL[["snRNA"]],
            "figure2_final_20260929_run07/tables/B_MP04_sample_usage_common_source.csv, cross-checked against stage2_comparability/sample_comparability.csv",
            sum(d$platform == "scRNA"), sum(d$platform == "snRNA"),
            paste0(W_IN, " x ", H_IN)))
fwrite(record, file.path(FIG_DIR, "panelD_platform_split_record.csv"))

cat("\nPlatform block sizes:\n"); print(table(split_f))
cat("Figure size:", W_IN, "x", H_IN, "in\nOutputs:", OUT_ROOT, "\n")
graphics.off()
cat("Graphics devices after explicit close:", length(dev.list()), "\n")
cat("End:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
sink()
