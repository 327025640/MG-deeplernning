suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})
options(stringsAsFactors = FALSE, scipen = 999)
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"))

OUT <- "E:/micro_output/分析审阅中心/待审阅/Fig6_疾病Top10_导师PPI药物可及性"
DISEASES <- c("LGG", "GBM", "MB", "LCBM", "AD", "EOAD", "PD", "DLB",
              "DLBD_PD", "ALS", "ALS_FTD", "FTD", "FTD_MND", "PSP", "HD", "MS")
SOURCE_CSV <- file.path(OUT, "source_data/Fig6_Top10_PPI_drug_dual_evidence_source.csv")
MASTER_CSV <- file.path(OUT, "tables/Fig6_all_disease_drug_master.csv")
DISPLAY_CSV <- file.path(OUT, "tables/dual_support_display_top3.csv")
SUMMARY_CSV <- file.path(OUT, "tables/displayed_drug_summary_by_disease.csv")
PALETTE_CSV <- "D:/data/microglia_new/authority/03_reference_configs/disease_palette_contract.csv"
REVISION <- "compact_reference_style_20261010"

for (p in c(SOURCE_CSV, MASTER_CSV, DISPLAY_CSV, SUMMARY_CSV, PALETTE_CSV))
  if (!file.exists(p)) stop("Missing frozen Figure input: ", p)
dir.create(file.path(OUT, "archive/Fig6_initial_layout_20261010/figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "archive/Fig6_initial_layout_20261010/source_data"), recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) as.character(digest::digest(file = path, algo = "sha256", serialize = FALSE))
read_csv <- function(path) fread(path, na.strings = c("", "NA"))
write_csv <- function(x, path) fwrite(x, path, na = "NA", bom = TRUE)
archive_once <- function(src, dst) {
  verify_destination <- FALSE
  if (file.exists(dst)) {
    if (identical(sha256(src), sha256(dst))) {
      verify_destination <- TRUE
    } else {
      backup <- paste0(dst, ".pre_redraw_", format(Sys.time(), "%Y%m%d_%H%M%S"))
      if (file.exists(backup) || !isTRUE(file.copy(src, backup, overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE)))
        stop("Could not preserve current Fig6 artifact before redraw: ", src)
      if (!identical(sha256(src), sha256(backup))) stop("Pre-redraw archive SHA256 verification failed: ", backup)
    }
  } else if (!isTRUE(file.copy(src, dst, overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE))) {
    stop("Could not preserve previous Fig6 artifact: ", src)
  } else {
    verify_destination <- TRUE
  }
  if (verify_destination && !identical(sha256(src), sha256(dst))) stop("Archive SHA256 verification failed: ", dst)
}

old_source <- read_csv(SOURCE_CSV)
master <- read_csv(MASTER_CSV)
display <- read_csv(DISPLAY_CSV)
display_summary <- read_csv(SUMMARY_CSV)
palette <- read_csv(PALETTE_CSV)
required_source <- c("source_disease", "Drug", "target_disease", "target_disease_order", "Z", "d_obs",
  "p_empirical", "q_network_global", "network_status", "network_cell_status", "reversal_cell_status",
  "q_reversal", "p_drug", "formal_deg_status", "local_signature_status", "reversal_status", "best_instance_id")
if (length(setdiff(required_source, names(old_source)))) stop("Frozen Figure source columns are incomplete")
if (nrow(old_source) != 528L || anyDuplicated(old_source, by = c("source_disease", "Drug", "target_disease")))
  stop("Frozen Figure source is not 33 x 16 unique cells")
if (!setequal(unique(old_source$target_disease), DISEASES) || any(table(old_source$target_disease_order) != 33L))
  stop("Frozen Figure source does not contain all 16 disease columns per displayed drug")
if (anyNA(old_source$target_disease_order) ||
    any(old_source$target_disease_order != match(old_source$target_disease, DISEASES)))
  stop("Frozen Figure disease labels do not match the fixed 01-16 column order")
if (nrow(master) != 12224L || nrow(display) != 33L || nrow(display_summary) != 16L)
  stop("Frozen result/display table sizes do not match the independently validated run")
effective_p_count <- sum(is.finite(master$p_empirical))
network_q_significant_count <- sum(is.finite(master$q_network_global) & master$q_network_global < 0.05)
if (effective_p_count != 11952L || network_q_significant_count != 0L)
  stop("Frozen master counts differ from the independently validated network P/BH totals")
if (any(display$dual_support != TRUE) || any(display[, .N, by = Disease]$N > 3L))
  stop("A display row does not satisfy the frozen dual-support/top-3 contract")

# Recheck the frozen figure cells against the final keyed master table; no network work is run here.
check <- merge(old_source, master, by.x = c("target_disease", "Drug"), by.y = c("Disease", "Drug"),
               all.x = TRUE, sort = FALSE, suffixes = c("_source", "_master"))
if (nrow(check) != 528L || anyNA(check$network_status_master)) stop("A frozen Figure cell cannot be joined to the result master")
for (field in c("Z", "d_obs", "p_empirical", "q_network_global", "network_status",
                "q_reversal", "p_drug", "formal_deg_status", "local_signature_status",
                "reversal_status", "best_instance_id")) {
  lhs <- check[[paste0(field, "_source")]]
  rhs <- check[[paste0(field, "_master")]]
  bad <- xor(is.na(lhs), is.na(rhs))
  both <- !is.na(lhs) & !is.na(rhs)
  bad[both] <- as.character(lhs[both]) != as.character(rhs[both])
  if (any(bad)) stop("Frozen Figure data changed during style-only revision: ", field)
}
rm(check)

palette_alias <- c(GBM = "GBM_IDHwt")
palette_lookup <- palette[Status == "frozen_2026-08-29", .(Disease_Abbrev, Disease_Color_Hex)]
palette_lookup <- rbind(palette_lookup, palette[Status == "frozen_2026-08-04", .(Disease_Abbrev, Disease_Color_Hex)])
colors <- setNames(palette_lookup$Disease_Color_Hex, palette_lookup$Disease_Abbrev)
for (d in DISEASES) {
  key <- if (d %in% names(palette_alias)) unname(palette_alias[[d]]) else d
  if (is.na(colors[key])) stop("Frozen disease color missing: ", d)
  colors[d] <- colors[key]
}
old_source[, target_color_hex := unname(colors[target_disease])]
if (anyNA(old_source$target_color_hex)) stop("Figure cells lack frozen disease colors")

blend <- function(a, b, w) {
  aa <- grDevices::col2rgb(a); bb <- grDevices::col2rgb(b)
  grDevices::rgb((aa[1, ] * (1 - w) + bb[1, ] * w) / 255,
                 (aa[2, ] * (1 - w) + bb[2, ] * w) / 255,
                 (aa[3, ] * (1 - w) + bb[3, ] * w) / 255)
}
old_source[, `:=`(network_cell_fill_hex = "#FFFFFF", network_cell_border_hex = target_color_hex)]
old_source[network_cell_status == "exploratory_close", network_cell_fill_hex := blend(target_color_hex, "#202020", 0.30)]
old_source[network_cell_status == "marginal_close", network_cell_fill_hex := blend(target_color_hex, "#FFFFFF", 0.46)]
old_source[network_cell_status == "undefined_z", `:=`(network_cell_fill_hex = "#B8B8B8", network_cell_border_hex = "#858585")]
old_source[network_cell_status == "no_mapped_candidate", `:=`(network_cell_fill_hex = "#D8D8D8", network_cell_border_hex = "#858585")]
old_source[network_cell_status == "random_failure", `:=`(network_cell_fill_hex = "#C8C8C8", network_cell_border_hex = "#858585")]
old_source[network_cell_status == "no_mapped_targets", `:=`(network_cell_fill_hex = "#EEEEEE", network_cell_border_hex = "#858585")]
old_source[, `:=`(network_cell_glyph = ifelse(network_cell_status == "no_mapped_targets", "diagonal_hatch", "none"),
                  reversal_glyph = fifelse(reversal_cell_status == "reversal_support", "filled_up_triangle",
                    fifelse(reversal_cell_status == "assessed_no_support", "open_circle", "horizontal_dash")),
                  reversal_glyph_hex = fifelse(reversal_cell_status == "reversal_support", target_color_hex, "#68717A"),
                  figure_style_revision = REVISION)]

# Preserve the original one-page figure/source before replacing the canonical deliverables.
pdf_out <- file.path(OUT, "figures/Fig6_Top10_PPI_drug_dual_evidence.pdf")
png_out <- file.path(OUT, "figures/Fig6_Top10_PPI_drug_dual_evidence.png")
archive_pdf <- file.path(OUT, "archive/Fig6_initial_layout_20261010/figures/Fig6_Top10_PPI_drug_dual_evidence_initial.pdf")
archive_png <- file.path(OUT, "archive/Fig6_initial_layout_20261010/figures/Fig6_Top10_PPI_drug_dual_evidence_initial.png")
archive_source <- file.path(OUT, "archive/Fig6_initial_layout_20261010/source_data/Fig6_Top10_PPI_drug_dual_evidence_initial_source.csv")
for (pair in list(c(pdf_out, archive_pdf), c(png_out, archive_png), c(SOURCE_CSV, archive_source)))
  archive_once(pair[[1L]], pair[[2L]])

# Tight, vector-native geometry based on the accepted Fig5 compact two-column reference.
W <- 300; H <- 225
cell <- 2.15; gap <- 0.45; step <- cell + gap; strip_w <- 16 * cell + 15 * gap
panel_x <- c(7, 153)
label_w <- 49; network_offset <- 50; reversal_offset <- network_offset + strip_w + 4
row_top <- 188; row_floor <- 47
row_items <- lapply(1:2, function(side) {
  srcs <- if (side == 1L) DISEASES[1:8] else DISEASES[9:16]
  items <- list()
  for (src in srcs) {
    candidates <- display[Disease == src]
    ss <- display_summary[Disease == src]
    if (nrow(ss) != 1L) stop("Missing display summary for ", src)
    items[[length(items) + 1L]] <- list(type = "header", disease = src,
      label = sprintf("%02d %s  %d/%d dual support", match(src, DISEASES), src, nrow(candidates), ss$n_dual_support_all))
    if (nrow(candidates)) {
      for (dr in candidates$Drug) items[[length(items) + 1L]] <- list(type = "drug", disease = src, Drug = dr, label = dr)
    } else {
      msg <- if (ss$group_status == "Reversal not assessed") "Reversal not assessed" else "No qualifying dual-support candidate"
      items[[length(items) + 1L]] <- list(type = "empty", disease = src, label = msg)
    }
    items[[length(items) + 1L]] <- list(type = "gap", disease = src)
  }
  items
})
row_weight <- function(items) sum(vapply(items, function(x) if (x$type == "gap") 0.32 else 1, numeric(1)))
max_weight <- max(vapply(row_items, row_weight, numeric(1)))
row_step <- (row_top - row_floor) / max_weight

draw_figure <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(xscale = c(0, W), yscale = c(0, H), clip = "off"))
  grid::grid.rect(x = grid::unit(W / 2, "native"), y = grid::unit(H / 2, "native"),
                  width = grid::unit(W, "mm"), height = grid::unit(H, "mm"),
                  gp = grid::gpar(fill = "white", col = NA))
  txt <- function(label, x, y, size = 6, just = "left", face = "plain", col = "black")
    grid::grid.text(label, x = grid::unit(x, "native"), y = grid::unit(y, "native"), just = just,
                    gp = grid::gpar(fontfamily = "Arial", fontsize = size, fontface = face, col = col))
  rect <- function(x, y, w, h, fill, border = "#777777", lwd = 0.4)
    grid::grid.rect(x = grid::unit(x, "native"), y = grid::unit(y, "native"),
                    width = grid::unit(w, "mm"), height = grid::unit(h, "mm"),
                    gp = grid::gpar(fill = fill, col = border, lwd = lwd))
  line <- function(x1, y1, x2, y2, col = "#68717A", lwd = 0.55)
    grid::grid.lines(x = grid::unit(c(x1, x2), "native"), y = grid::unit(c(y1, y2), "native"),
                     gp = grid::gpar(col = col, lwd = lwd))
  txt("Fig6  Top10 candidate network proximity and reversal support", W / 2, 217.5, 9, "centre", "bold")
  for (side in 1:2) {
    px <- panel_x[side]
    diseases_side <- if (side == 1L) "Source diseases 01–08" else "Source diseases 09–16"
    txt(diseases_side, px + 65, 211.2, 7.4, "centre", "bold")
    txt("Source disease / drug", px, 204.5, 6.3, "left", "bold")
    txt("Network proximity", px + network_offset + strip_w / 2, 204.5, 6.3, "centre", "bold")
    txt("Reversal support", px + reversal_offset + strip_w / 2, 204.5, 6.3, "centre", "bold")
    for (j in seq_along(DISEASES)) {
      nx <- px + network_offset + (j - 1) * step + cell / 2
      rx <- px + reversal_offset + (j - 1) * step + cell / 2
      txt(sprintf("%02d", j), nx, 198.3, 5.4, "centre", col = "#303030")
      txt(sprintf("%02d", j), rx, 198.3, 5.4, "centre", col = "#303030")
    }
    y <- row_top
    for (item in row_items[[side]]) {
      if (item$type == "gap") {
        y <- y - 0.32 * row_step
        next
      }
      yy <- y - row_step / 2
      if (item$type == "header") {
        rect(px + 0.65, yy, 1.25, row_step * 0.72, colors[item$disease], colors[item$disease], 0.2)
        txt(item$label, px + 2.5, yy, 6.2, "left", "bold")
      } else if (item$type == "empty") {
        txt(item$label, px + 2.5, yy, 5.8, "left")
      } else {
        txt(item$label, px + 2.5, yy, 5.9, "left")
        cells <- old_source[source_disease == item$disease & Drug == item$Drug]
        setorder(cells, target_disease_order)
        if (nrow(cells) != 16L || !identical(as.integer(cells$target_disease_order), 1:16))
          stop("Figure source cells do not preserve the frozen 01–16 disease order")
        for (j in seq_len(16L)) {
          cstat <- cells[j]
          nx <- px + network_offset + (j - 1) * step + cell / 2
          rx <- px + reversal_offset + (j - 1) * step + cell / 2
          rect(nx, yy, cell, cell, cstat$network_cell_fill_hex, cstat$network_cell_border_hex, 0.34)
          if (cstat$network_cell_glyph == "diagonal_hatch")
            line(nx - cell * 0.34, yy - cell * 0.34, nx + cell * 0.34, yy + cell * 0.34, "#777777", 0.45)
          if (cstat$reversal_glyph == "filled_up_triangle") {
            grid::grid.polygon(x = grid::unit(c(rx - 0.72, rx + 0.72, rx), "native"),
              y = grid::unit(c(yy - 0.52, yy - 0.52, yy + 0.76), "native"),
              gp = grid::gpar(fill = cstat$reversal_glyph_hex, col = cstat$reversal_glyph_hex, lwd = 0.35))
          } else if (cstat$reversal_glyph == "open_circle") {
            grid::grid.circle(x = grid::unit(rx, "native"), y = grid::unit(yy, "native"), r = grid::unit(0.67, "mm"),
              gp = grid::gpar(fill = "white", col = "#68717A", lwd = 0.52))
          } else {
            line(rx - 0.7, yy, rx + 0.7, yy, "#68717A", 0.65)
          }
        }
      }
      y <- y - row_step
    }
  }

  # Symbolic legends use the same fills and glyphs as the actual cells.
  txt("Network", 7, 40.5, 6.4, "left", "bold")
  legend_x <- c(31, 67, 103, 133, 164, 195, 232)
  ad <- colors[["AD"]]
  rect(legend_x[1], 40.5, 2.5, 2.5, blend(ad, "#202020", 0.30), ad, 0.4); txt("Z < −1.96 (exploratory)", legend_x[1] + 2.4, 40.5, 5.9)
  rect(legend_x[2], 40.5, 2.5, 2.5, blend(ad, "#FFFFFF", 0.46), ad, 0.4); txt("−1.96 ≤ Z < −1.645", legend_x[2] + 2.4, 40.5, 5.9)
  rect(legend_x[3], 40.5, 2.5, 2.5, "white", ad, 0.45); txt("other finite Z", legend_x[3] + 2.4, 40.5, 5.9)
  rect(legend_x[4], 40.5, 2.5, 2.5, "#B8B8B8", "#858585", 0.35); txt("Undefined Z", legend_x[4] + 2.4, 40.5, 5.8)
  rect(legend_x[5], 40.5, 2.5, 2.5, "#D8D8D8", "#858585", 0.35); txt("No mapped candidate", legend_x[5] + 2.4, 40.5, 5.8)
  rect(legend_x[6], 40.5, 2.5, 2.5, "#C8C8C8", "#858585", 0.35); txt("Random failure", legend_x[6] + 2.4, 40.5, 5.8)
  rect(legend_x[7], 40.5, 2.5, 2.5, "#EEEEEE", "#858585", 0.35)
  line(legend_x[7] - 0.75, 40.5 - 0.75, legend_x[7] + 0.75, 40.5 + 0.75, "#777777", 0.45)
  txt("No mapped targets", legend_x[7] + 2.4, 40.5, 5.8)

  txt("Reversal", 7, 33.8, 6.4, "left", "bold")
  grid::grid.polygon(x = grid::unit(c(31 - 0.8, 31 + 0.8, 31), "native"),
    y = grid::unit(c(33.3, 33.3, 34.7), "native"), gp = grid::gpar(fill = colors[["AD"]], col = colors[["AD"]], lwd = 0.35))
  txt("q_reversal < 0.05; triangle uses target-disease color", 34, 33.8, 5.9)
  grid::grid.circle(x = grid::unit(122, "native"), y = grid::unit(33.8, "native"), r = grid::unit(0.67, "mm"),
                    gp = grid::gpar(fill = "white", col = "#68717A", lwd = 0.52))
  txt("assessed, no support", 125, 33.8, 5.9)
  line(217, 33.8, 219, 33.8, "#68717A", 0.65); txt("not assessed", 222, 33.8, 5.9)

  # Disease key is two rows of eight, with the same numbering used above both strips.
  key_x0 <- 7; key_step <- (W - 14) / 8
  for (j in seq_along(DISEASES)) {
    row <- if (j <= 8L) 1 else 2
    col <- if (j <= 8L) j else j - 8L
    yy <- if (row == 1L) 25.7 else 18.8
    xx <- key_x0 + (col - 1L) * key_step
    rect(xx + 1.0, yy, 2.3, 2.3, colors[[DISEASES[j]]], colors[[DISEASES[j]]], 0.2)
    txt(sprintf("%02d %s", j, DISEASES[j]), xx + 3.0, yy, 5.8)
  }
  txt(sprintf("Network empirical-P global BH: q < 0.05 for %s / %s effective tests. Network Z thresholds are exploratory; reversal is historical signature support. Neither is drug efficacy evidence.",
      format(network_q_significant_count, big.mark = ","), format(effective_p_count, big.mark = ",")),
      W / 2, 6.5, 5.9, "centre")
  grid::popViewport()
}

# Retain the exact numeric/state columns while adding only figure color and glyph metadata.
original_columns <- names(old_source)
write_csv(old_source, SOURCE_CSV)
device_png <- file.path(OUT, "figures/Fig6_Top10_PPI_drug_dual_evidence.png")
device_pdf <- file.path(OUT, "figures/Fig6_Top10_PPI_drug_dual_evidence.pdf")
grDevices::png(device_png, width = W / 25.4, height = H / 25.4, units = "in", res = 600, type = "cairo", bg = "white")
draw_figure(); grDevices::dev.off()
grDevices::cairo_pdf(device_pdf, width = W / 25.4, height = H / 25.4, family = "Arial", onefile = TRUE, bg = "white")
draw_figure(); grDevices::dev.off()

final_source <- read_csv(SOURCE_CSV)
if (nrow(final_source) != 528L || !all(original_columns %in% names(final_source)) ||
    nrow(display) != 33L || any(final_source$figure_style_revision != REVISION))
  stop("Final figure source integrity check failed")
for (field in original_columns) {
  lhs <- old_source[[field]]; rhs <- final_source[[field]]
  bad <- xor(is.na(lhs), is.na(rhs)); both <- !is.na(lhs) & !is.na(rhs)
  bad[both] <- as.character(lhs[both]) != as.character(rhs[both])
  if (any(bad)) stop("Style revision changed original source values: ", field)
}
if (!all(file.exists(c(device_png, device_pdf, SOURCE_CSV))) || file.info(device_pdf)$size < 10000 || file.info(device_png)$size < 10000)
  stop("Compact Fig6 PDF/PNG/source output is missing or undersized")

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) == 1L) normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE) else stop("Cannot identify the standalone renderer path")
manifest_paths <- c(SOURCE_CSV, MASTER_CSV, DISPLAY_CSV, SUMMARY_CSV, PALETTE_CSV, script_path,
                    archive_pdf, archive_png, archive_source, device_pdf, device_png)
manifest <- rbindlist(lapply(manifest_paths, function(p) data.table(path = p, size_bytes = file.info(p)$size, sha256 = sha256(p))))
write_csv(manifest, file.path(OUT, "tables/Fig6_figure_revision_manifest.csv"))
cat("STYLE_REVISION_COMPLETE\nPDF=", device_pdf, "\nPNG=", device_png, "\nSOURCE=", SOURCE_CSV,
    "\nSHA256_RENDERER=", sha256(script_path), "\nMINIMUM_COMBINED_WIDTH_MM=300\n", sep = "")
