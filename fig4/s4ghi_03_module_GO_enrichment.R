# =============================================================================
# S4g-i Step 3: PTG module GO enrichment (S4i)
# Method: original trajectory approach = compareCluster(fun = "enrichGO")
# (DIFFERENT from the old subcluster ORA enricher+msigdbr C5 pipeline - that
# one stays untouched). Human org.Hs.eg.db, GO:BP.
#   input: PTG01-04 dynamic genes from step 2 (this project's real grouping)
#   background: genes passing the >=10% detection gate that entered the dynamic
#     test and are GO:BP-annotatable (NOT the significant union, NOT the whole DB)
#   params (project-transparent, original undisclosed in this section):
#     BH, pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500
#   display: top 5 FDR terms per module; bubble size = GeneRatio,
#     colour = -log10(BH FDR) (REVERSED vs original - declared adaptation)
#   full term results and mapping audits saved; no-significant modules are
#   shown honestly.
# Both orientations: PTG membership per orientation (identical here, verified
# by the intersection mapping in step 2) - each gets its own tables/figures.
# =============================================================================
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(data.table)
  library(ggplot2)
})
options(stringsAsFactors = FALSE)
set.seed(42)

OUT  <- "E:/micro_output/fig3/trajectory_S4"
DATA <- file.path(OUT, "data"); TAB <- file.path(OUT, "tables")
FIG  <- file.path(OUT, "figures"); SRC <- file.path(OUT, "source_data")
t0 <- Sys.time()

PARAMS <- list(pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05,
               minGSSize = 10, maxGSSize = 500)
PTG_COLS <- c(PTG01 = "#4477AA", PTG02 = "#EE6677", PTG03 = "#228833", PTG04 = "#AA3377")

run_enrichment <- function(orient) {
  cat("\n=====", orient, "=====\n")
  dgt  <- fread(file.path(DATA, orient, "dynamic_gene_test.csv"))
  memb <- fread(file.path(DATA, orient, "s4h_ptg_membership.csv"))
  cand <- fread(file.path(TAB, orient, "s4h_dynamic_test_candidates.csv"))

  # SYMBOL -> ENTREZ mapping with audit (multi-mapping kept, unmapped recorded)
  sym_all <- unique(c(cand$gene_short_name, memb$gene_id))
  bm <- clusterProfiler::bitr(sym_all, fromType = "SYMBOL", toType = "ENTREZID",
                              OrgDb = org.Hs.eg.db)
  bm_dt <- as.data.table(bm)
  multi <- bm_dt[, .N, by = SYMBOL][N > 1]
  unmapped <- setdiff(sym_all, bm_dt$SYMBOL)
  fwrite(bm_dt, file.path(TAB, orient, "s4i_symbol_entrez_map.csv"))
  fwrite(data.table(symbol = unmapped),
         file.path(TAB, orient, "s4i_unmapped_symbols.csv"),
         row.names = FALSE)
  cat(sprintf("mapping: %d/%d symbols mapped | multi-ENTREZ symbols: %d | unmapped: %d\n",
              length(unique(bm_dt$SYMBOL)), length(sym_all), nrow(multi), length(unmapped)))

  sym2ez <- bm_dt
  to_ez <- function(genes) {
    unique(sym2ez[SYMBOL %in% genes]$ENTREZID)
  }
  # background: candidates (detection gate, entered the dynamic test) with GO:BP annotation
  bg_ez <- to_ez(cand$gene_short_name)
  cat("background ENTREZ genes:", length(bg_ez), "\n")

  gene_sets <- lapply(sort(unique(memb$PTG)), function(g) {
    to_ez(memb[PTG == g]$gene_id)
  })
  names(gene_sets) <- sort(unique(memb$PTG))
  cat("module ENTREZ sizes:", paste(sprintf("%s=%d", names(gene_sets),
                                            sapply(gene_sets, length)), collapse = ", "), "\n")

  cc <- compareCluster(geneClusters = gene_sets, fun = "enrichGO",
                       OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                       universe = bg_ez,
                       pAdjustMethod = PARAMS$pAdjustMethod,
                       pvalueCutoff = PARAMS$pvalueCutoff,
                       qvalueCutoff = PARAMS$qvalueCutoff,
                       minGSSize = PARAMS$minGSSize,
                       maxGSSize = PARAMS$maxGSSize,
                       readable = TRUE)
  cc_dt <- as.data.table(cc)
  fwrite(cc_dt, file.path(TAB, orient, "s4i_compareCluster_enrichGO_full.csv"))
  cat("enrichGO significant terms:", nrow(cc_dt), "\n")

  # top 5 FDR per module (fixed rule); no hand-picking
  top5 <- cc_dt[, head(.SD[order(-(-log10(p.adjust)))], 5), by = Cluster]
  fwrite(top5, file.path(TAB, orient, "s4i_top5_per_module.csv"))
  if (nrow(cc_dt) == 0) {
    cat("No module reached the significance thresholds - figure will show the\n",
        "'no significant terms' state honestly (no invented labels).\n")
    return(invisible(list(cc = cc_dt, top5 = top5)))
  }

  # GeneRatio numeric + components (denominator = enrichGO annotated size)
  top5[, GeneRatio_num := sapply(GeneRatio, function(x) {
    a <- as.numeric(strsplit(x, "/")[[1]]); a[1] / a[2]
  })]
  top5[, log10FDR := -log10(p.adjust)]
  top5[, log10FDR_disp := pmin(log10FDR, 20)]   # display cap only

  # ------------------------------ bubble plot ------------------------------
  # y = GO term (grouped by module, FDR within module), x = module
  top5[, Cluster := factor(Cluster, levels = names(PTG_COLS))]
  top5 <- top5[order(Cluster, p.adjust)]
  top5[, term_label := factor(Description, levels = rev(unique(Description)))]
  gp <- ggplot(top5, aes(x = Cluster, y = term_label)) +
    geom_point(aes(size = GeneRatio_num, color = log10FDR_disp)) +
    scale_color_gradient(low = "grey90", high = "#B2182B", name = "-log10(BH FDR)") +
    scale_size_continuous(range = c(2.5, 7.5), name = "GeneRatio") +
    facet_grid(~ Cluster, scales = "free_x", space = "free_x") +
    labs(title = "GO:BP enrichment of pseudotime gene modules (PTG01-04)",
         subtitle = sprintf("compareCluster + enrichGO (org.Hs.eg.db); background = detection-gated genes (n=%d); BH q<0.05\npoint size = GeneRatio, colour = -log10 FDR (display cap 20)",
                            length(bg_ez)),
         x = "Gene module", y = "GO biological process term") +
    theme_bw(base_size = 10) +
    theme(plot.subtitle = element_text(size = 8.5),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = "grey92"),
          strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank())
  ggsave(file.path(FIG, orient, "S4i_module_GO_bubble.png"), gp,
         width = 11, height = 7, dpi = 600, bg = "white")
  ggsave(file.path(FIG, orient, "S4i_module_GO_bubble.pdf"), gp,
         width = 11, height = 7, bg = "white", useDingbats = FALSE)
  fwrite(top5[, .(Cluster, Description, GeneRatio, GeneRatio_num,
                  BgRatio, pvalue, p.adjust, log10FDR, log10FDR_disp,
                  geneID, Count)],
         file.path(SRC, orient, "S4i_bubble_source.csv"))
  cat("S4i bubble written:", orient, "\n")
  invisible(list(cc = cc_dt, top5 = top5, bg = bg_ez))
}

res_def <- run_enrichment("default")
res_c00 <- run_enrichment("C00_associated_state")

capture.output(sessionInfo(), file = file.path(OUT, "logs", "sessionInfo_s4ghi_03.txt"))
cat(sprintf("\nStep 3 done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
