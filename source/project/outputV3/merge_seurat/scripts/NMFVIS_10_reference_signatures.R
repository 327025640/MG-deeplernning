# =============================================================================
# NMFVIS_10_reference_signatures.R — 微胶质状态/DAM/衰老参考签名富集
# 权威参考: 08_plot_reference_signature_enrichment_transposed.R (fgsea 权重版)
#           + 90c_reference_signature_enrichment.R (top100 超几何版)
# 两法并列: A) consensus W 全权重 fgseaMultilevel; B) 每 MP top100 基因超几何右尾
# 全局 BH; 12 签名 x 3 组（Aging/Homeostatic/DAM, HuMicA v1）
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(fgsea); library(ggplot2)
})

ROOT <- "E:/micro_output/NMF"
FIG  <- file.path(ROOT, "figures"); SRC <- file.path(ROOT, "source_data")
TBL  <- file.path(ROOT, "tables")
for (d in c(FIG, SRC, TBL)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
NV <- "E:/micro_output/v3_dims16/05_nmf"
GS_CSV <- "D:/data/microglia_new/input/gene_sets/humica_dam_aging_gene_sets_v1.csv"

pal <- fread(file.path(ROOT, "source_data", "NMFVIS_MP_palette.csv"))
mp_ord <- pal$MP
mp_col <- setNames(pal$Color_Hex, pal$MP)

log_msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...), "\n")
log_msg("=== NMFVIS_10 reference signatures start ===")

cns <- readRDS(file.path(NV, "02_consensus", "consensus_programs.rds"))
W <- cns$consensus_W
stopifnot(identical(colnames(W), mp_ord), nrow(W) == 5000L)
universe <- rownames(W)

gs <- fread(GS_CSV)
sig_sets <- lapply(split(gs$Gene, gs$GeneSet), unique)
groups <- data.table(Signature = names(sig_sets))
groups[, SignatureGroup := fifelse(grepl("Homeos", Signature), "Homeostatic",
                            fifelse(grepl("DAM|Inflam|Lipo|DIM|MAC|Ribo", Signature), "DAM", "Aging"))]
sig_sets <- lapply(sig_sets, function(g) intersect(g, universe))
log_msg("signatures: ", paste(names(sig_sets), collapse = ", "))

# ---- A) 权重 GSEA（全 5000 排序）----
res_list <- list()
for (m in mp_ord) {
  ranks <- sort(setNames(W[, m], universe), decreasing = TRUE)
  er <- fgseaMultilevel(pathways = sig_sets, stats = ranks,
                        minSize = 5, maxSize = 500, eps = 1e-50)
  er <- as.data.table(er)
  er[, MP := m]
  res_list[[m]] <- er
}
gsea <- rbindlist(res_list, fill = TRUE)
gsea[, FDR_mp := padj]
gsea[, FDR := p.adjust(pval, "BH")]
gsea[, neglog10FDR := -log10(pmax(FDR, 1e-300))]
gsea <- merge(gsea, groups, by.x = "pathway", by.y = "Signature", all.x = TRUE)
fwrite(gsea[, .(MP, Signature = pathway, SignatureGroup, NES, P = pval, FDR_mp, FDR,
                neglog10FDR, size)],
       file.path(SRC, "NMFVIS_MP_reference_signature_GSEA_source.csv"))

gsea[, MP := factor(MP, levels = mp_ord)]
gsea[, Signature := factor(pathway, levels = rev(names(sig_sets)))]
p <- ggplot(gsea, aes(MP, Signature, size = neglog10FDR, fill = NES)) +
  geom_point(shape = 21, stroke = 0.4) +
  facet_grid(SignatureGroup ~ ., scales = "free_y", space = "free_y") +
  scale_size(range = c(0.8, 4.5), guide = guide_legend(title = "-log10(FDR)")) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "NES") +
  theme_minimal(base_family = "Arial", base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = mp_col)) +
  labs(title = "MP x reference signature GSEA (HuMicA v1; full consensus loadings)",
       subtitle = "HuMicA - Martins-Ferreira et al., Nat Commun 2025, DOI 10.1038/s41467-025-56124-1",
       x = NULL, y = NULL)
ggsave(file.path(FIG, "Panel_reference_signature_GSEA.png"), p,
       width = 6.2, height = 5.4, dpi = 600, bg = "white")
ggsave(file.path(FIG, "Panel_reference_signature_GSEA.pdf"), p,
       width = 6.2, height = 5.4, device = cairo_pdf, bg = "white")

# ---- B) top100 超几何 ----
hyp_list <- list()
for (m in mp_ord) {
  top100 <- names(sort(W[, m], decreasing = TRUE))[1:100]
  top100 <- intersect(top100, universe)
  for (sn in names(sig_sets)) {
    m_sz <- length(sig_sets[[sn]])
    k <- length(intersect(top100, sig_sets[[sn]]))
    P <- if (k == 0) 1 else phyper(k - 1, m_sz, length(universe) - m_sz,
                                   length(top100), lower.tail = FALSE)
    hyp_list[[length(hyp_list) + 1L]] <- data.table(
      MP = m, Signature = sn, SignatureGroup = groups[Signature == sn, SignatureGroup],
      n_top100 = length(top100), n_sig_in_universe = m_sz, overlap = k,
      expected = length(top100) * m_sz / length(universe),
      P = P)
  }
}
hyp <- rbindlist(hyp_list)
hyp[, FDR := p.adjust(P, "BH")]
hyp[, log2_enrichment := round(log2((overlap + 1e-9) / (expected + 1e-9)), 4)]
hyp[, Significant_FDR05 := FDR < 0.05]
fwrite(hyp, file.path(SRC, "NMFVIS_reference_signature_hypergeom_source.csv"))

hyp[, MP := factor(MP, levels = mp_ord)]
hyp[, Signature := factor(Signature, levels = rev(names(sig_sets)))]
p2 <- ggplot(hyp, aes(MP, Signature, size = overlap, fill = log2_enrichment)) +
  geom_point(shape = 21, stroke = ifelse(hyp$Significant_FDR05, 0.8, 0.3),
             colour = ifelse(hyp$Significant_FDR05, "black", "grey85")) +
  facet_grid(SignatureGroup ~ ., scales = "free_y", space = "free_y") +
  scale_size(range = c(1, 8), guide = guide_legend(title = "overlap")) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "log2 obs/exp") +
  theme_minimal(base_family = "Arial", base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = mp_col)) +
  labs(title = "MP x reference signature hypergeometric enrichment (top100 vs 5000 HVG)",
       subtitle = "HuMicA - Martins-Ferreira et al., Nat Commun 2025, DOI 10.1038/s41467-025-56124-1",
       x = NULL, y = NULL)
ggsave(file.path(FIG, "Panel_reference_signature_hypergeom.png"), p2,
       width = 6.6, height = 5.6, dpi = 600, bg = "white")
ggsave(file.path(FIG, "Panel_reference_signature_hypergeom.pdf"), p2,
       width = 6.6, height = 5.6, device = cairo_pdf, bg = "white")
log_msg("== NMFVIS_10 done ==")
