#!/usr/bin/env Rscript
# Figure 2 局部绘图修订：超几何参考签名图改为 MP 行、功能签名列。
# 只读取现有 source CSV；不重新计算 NMF、基因筛选或超几何检验。

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

ROOT <- "E:/micro_output/NMF"
FIG <- file.path(ROOT, "figures")
SRC <- file.path(ROOT, "source_data")
HYP_FILE <- file.path(SRC, "NMFVIS_reference_signature_hypergeom_source.csv")
PANELE_FILE <- file.path(SRC, "PanelE_MP_disease_mean_usage_source.csv")
PALETTE_FILE <- file.path(SRC, "NMFVIS_MP_palette.csv")
PNG_FILE <- file.path(FIG, "Panel_reference_signature_hypergeom.png")
PDF_FILE <- file.path(FIG, "Panel_reference_signature_hypergeom.pdf")

for (path in c(HYP_FILE, PANELE_FILE, PALETTE_FILE)) {
  if (!file.exists(path)) stop("Required input does not exist: ", path)
}

hyp <- fread(HYP_FILE)
panele <- fread(PANELE_FILE)
palette <- fread(PALETTE_FILE)

required_hyp <- c("MP", "Signature", "SignatureGroup", "overlap", "log2_enrichment", "Significant_FDR05")
if (length(setdiff(required_hyp, names(hyp)))) {
  stop("Hypergeometric source is missing columns: ", paste(setdiff(required_hyp, names(hyp)), collapse = ", "))
}
if (!"MP" %in% names(panele) || !"MP" %in% names(palette) || !"Color_Hex" %in% names(palette)) {
  stop("Panel E source or frozen MP palette is missing required columns.")
}

# MP order is read from the current Panel E source, then checked against its frozen palette.
mp_order <- unique(as.character(panele$MP))
if (length(mp_order) != 12L || anyDuplicated(mp_order) ||
    !identical(mp_order, as.character(palette$MP))) {
  stop("Panel E MP order does not match the frozen 12-MP palette.")
}
if (!setequal(unique(as.character(hyp$MP)), mp_order)) {
  stop("Hypergeometric source MP values differ from current Panel E.")
}
signature_order <- unique(as.character(hyp$Signature))
if (length(signature_order) != 12L || anyDuplicated(signature_order)) {
  stop("Expected 12 unique functional signatures in the existing source.")
}
if (nrow(hyp) != length(mp_order) * length(signature_order) ||
    anyDuplicated(hyp[, .(MP, Signature)])) {
  stop("Hypergeometric source is not a complete unique MP-by-signature grid.")
}
if (any(!is.finite(hyp$overlap)) || any(!is.finite(hyp$log2_enrichment)) ||
    anyNA(hyp$Significant_FDR05)) {
  stop("Hypergeometric source contains missing or non-finite plotted values.")
}

hyp[, MP := factor(as.character(MP), levels = rev(mp_order))]
hyp[, Signature := factor(as.character(Signature), levels = signature_order)]
group_order <- c("Aging", "DAM", "Homeostatic")
if (!setequal(unique(as.character(hyp$SignatureGroup)), group_order)) {
  stop("Unexpected signature groups in hypergeometric source.")
}
hyp[, SignatureGroup := factor(as.character(SignatureGroup), levels = group_order)]

p <- ggplot(hyp, aes(x = Signature, y = MP, size = overlap, fill = log2_enrichment)) +
  geom_point(data = hyp[Significant_FDR05 == FALSE], shape = 21,
             stroke = 0.3, colour = "grey85") +
  geom_point(data = hyp[Significant_FDR05 == TRUE], shape = 21,
             stroke = 0.8, colour = "black") +
  facet_grid(. ~ SignatureGroup, scales = "free_x", space = "free_x", drop = TRUE) +
  scale_size(range = c(1, 8), guide = guide_legend(title = "overlap")) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "log2 obs/exp") +
  scale_y_discrete(drop = FALSE) +
  labs(title = "MP x reference signature hypergeometric enrichment (top 100 vs 5,000 HVG)",
       subtitle = "HuMiCA - Martins-Ferreira et al., Nat Commun 2025, DOI 10.1038/s41467-025-56124-1",
       x = NULL, y = NULL) +
  theme_minimal(base_family = "Arial", base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 7),
        axis.text.y = element_text(colour = "black", size = 8),
        strip.text.x = element_text(face = "bold", size = 8),
        plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 8.5),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))

# Preserve the previous rendered files before replacing them; source data remain unchanged.
existing <- c(PNG_FILE, PDF_FILE)[file.exists(c(PNG_FILE, PDF_FILE))]
if (length(existing)) {
  archive_dir <- file.path(ROOT, "archive", paste0("reference_hypergeom_pre_transpose_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  if (!dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create archive directory: ", archive_dir)
  }
  copied <- file.copy(existing, archive_dir, overwrite = FALSE)
  if (!all(copied)) stop("Could not preserve all previous figure files before export.")
}

ggsave(PNG_FILE, p, width = 7.0, height = 5.1, units = "in", dpi = 600, bg = "white")
ggsave(PDF_FILE, p, width = 7.0, height = 5.1, units = "in", device = cairo_pdf, bg = "white")

message("Rendered transposed hypergeometric signature panel using frozen source values.")
message("MP visual order (top to bottom): ", paste(mp_order, collapse = ", "))
message("PNG: ", PNG_FILE)
message("PDF: ", PDF_FILE)
