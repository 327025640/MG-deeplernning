# =============================================================================
# S4g-i Step 2: pseudotime-dynamic genes (S4h) + module heatmaps + branch panel
# Method (contract v1, user-confirmed):
#   - candidates = genes detected in >= 10% of the 16,800 trajectory cells
#     (4,868 genes); dynamic test = monocle::differentialGeneTest,
#     full ~sm.ns(Pseudotime, df=3), reduced ~1, relative_expr=TRUE, cores=1
#     (df=3 is the local default made explicit; original states ~sm.ns only)
#   - project exploratory screen: status OK & qval < 0.05 (original final
#     threshold undisclosed - recorded in S4ghi_method_alignment.csv)
#   - smooth/cluster flow replicates monocle 2.34.0 plot_pseudotime_heatmap
#     internals (verified against deparse): 100-pt grid -> genSmoothCurves ->
#     log10(x+1) -> per-gene Z (center+scale) -> display clip +-3 ->
#     row dist (1-r)/2, ward.D2, cluster_cols = FALSE, cutree k = 4
#   - PTG01-04 named by increasing median peak pseudotime (project rule)
#   - smooth matrix / Z matrix / grid / tree / membership cached explicitly
#     (return_heatmap object is NOT used as data)
#   - both orientations computed from their own CDS; C00 orientation runs a
#     50-gene spot check after input-identity verification, then reuses the
#     default test (inputs verified identical; MST bookkeeping differs only)
#   - branch panel: paths recovered from the projected principal graph
#     (root cell -> branch point Y_48 -> leaves); shared trunk marked;
#     no extrapolation; coverage tables saved
# =============================================================================
suppressPackageStartupMessages({
  library(monocle)
  library(ComplexHeatmap)
  library(circlize)
  library(data.table)
  library(igraph)
})
options(stringsAsFactors = FALSE)
set.seed(42)

OUT  <- "E:/micro_output/fig3/trajectory_S4"
DATA <- file.path(OUT, "data"); TAB <- file.path(OUT, "tables")
FIG  <- file.path(OUT, "figures"); SRC <- file.path(OUT, "source_data")
LOG  <- file.path(OUT, "logs")
t0 <- Sys.time()

P_CDS <- c(default = "E:/micro_output/fig3/data/Fig2i_monocle2_cds_default.rds",
           C00_associated_state = "E:/micro_output/fig3/data/Fig2i_monocle2_cds_C00_root.rds")
P_MAP <- "E:/micro_output/fig3/fig4/tables/cluster_name_mapping.csv"
PTG_COLS <- c(PTG01 = "#4477AA", PTG02 = "#EE6677", PTG03 = "#228833", PTG04 = "#AA3377")
BLOCK <- 400L

map <- fread(P_MAP, encoding = "UTF-8")

# ------------------------------------------------------------------ helpers --
run_dgt_chunked <- function(cds, genes, orient) {
  blk_dir <- file.path(DATA, orient, "dgt_blocks")
  dir.create(blk_dir, showWarnings = FALSE, recursive = TRUE)
  blocks <- split(genes, ceiling(seq_along(genes) / BLOCK))
  res <- list()
  for (b in seq_along(blocks)) {
    fp <- file.path(blk_dir, sprintf("block_%03d.rds", b))
    if (file.exists(fp)) {
      res[[b]] <- readRDS(fp)
      cat(sprintf("  [%s block %d/%d] cached (%d genes)\n", orient, b,
                  length(blocks), length(blocks[[b]])))
      next
    }
    g <- blocks[[b]]
    r <- monocle::differentialGeneTest(cds[g, ],
                                       fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)",
                                       reducedModelFormulaStr = "~1",
                                       relative_expr = TRUE, cores = 1)
    r$gene_id <- rownames(r)
    res[[b]] <- as.data.table(r)
    saveRDS(res[[b]], fp)
    cat(sprintf("  [%s block %d/%d] tested %d genes (%.1f min elapsed)\n",
                orient, b, length(blocks), length(g),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  out <- rbindlist(res)
  setorder(out, gene_id)
  out
}

smooth_and_modules <- function(cds, dyn_genes, orient) {
  fp_cache <- file.path(DATA, orient, "s4h_smooth_and_modules.rds")
  if (file.exists(fp_cache)) {
    cat("smooth/modules cached for", orient, "- reusing\n")
    return(readRDS(fp_cache))
  }
  newdata <- data.frame(Pseudotime = seq(min(pData(cds)$Pseudotime),
                                         max(pData(cds)$Pseudotime),
                                         length.out = 100))
  m <- monocle:::genSmoothCurves(cds[dyn_genes, ], cores = 1,
                                 trend_formula = "~sm.ns(Pseudotime, df=3)",
                                 relative_expr = TRUE, new_data = newdata)
  m <- log10(m + 1)
  z <- Matrix::t(scale(Matrix::t(m), center = TRUE))   # per-gene Z on the grid
  z[!is.finite(z)] <- 0
  z_disp <- pmin(pmax(z, -3), 3)
  cm <- cor(t(z_disp))                 # gene x gene on the 100-pt grid
  cm[!is.finite(cm)] <- 0              # flat genes (sd=0) would give NA
  row_dist <- as.dist((1 - cm)/2)
  hc <- hclust(row_dist, method = "ward.D2")
  mod <- cutree(hc, k = 4)
  # PTG numbering: modules ordered by median peak position on the grid
  pk <- sapply(names(table(mod)), function(g) {
    ii <- which(mod == g)
    median(apply(z[ii, , drop = FALSE], 1, which.max))
  })
  ord_mod <- order(pk)
  map_mod <- setNames(sprintf("PTG%02d", seq_along(pk)), names(sort(pk)))  # sort by peak position; names(sort(pk)) = module ids in that order
  ptg <- unname(map_mod[as.character(mod)])
  names(ptg) <- dyn_genes
  list(grid = newdata$Pseudotime, smooth = m, z = z, z_disp = z_disp,
       hc = hc, ptg = ptg, peak_median = pk)
}

save_modules <- function(sm, orient) {
  fwrite(as.data.table(sm$smooth, keep.rownames = "gene_id"),
         file.path(DATA, orient, "s4h_smooth_log10_grid100.csv.gz"))
  fwrite(as.data.table(sm$z, keep.rownames = "gene_id"),
         file.path(DATA, orient, "s4h_z_matrix_grid100.csv.gz"))
  fwrite(data.table(gene_id = names(sm$ptg), PTG = sm$ptg),
         file.path(DATA, orient, "s4h_ptg_membership.csv"))
  fwrite(data.table(grid_point = seq_along(sm$grid), pseudotime = sm$grid),
         file.path(DATA, orient, "s4h_time_grid.csv"))
  saveRDS(sm$hc, file.path(DATA, orient, "s4h_row_dendrogram_hclust.rds"))
  saveRDS(sm, file.path(DATA, orient, "s4h_smooth_and_modules.rds"))
}

# ==================================================================== step 2.1
cds_d <- readRDS(P_CDS["default"])
fd <- fData(cds_d)
n_cells <- ncol(cds_d)
cand <- rownames(fd)[fd$num_cells_expressed >= 0.10 * n_cells]
cat(sprintf("Candidate genes (detected >=10%% of %d cells): %d\n", n_cells, length(cand)))
fwrite(data.table(gene_id = cand,
                  gene_short_name = fd[cand, ]$gene_short_name,
                  num_cells_expressed = fd[cand, ]$num_cells_expressed),
       file.path(TAB, "default", "s4h_dynamic_test_candidates.csv"))

dgt_fp <- file.path(DATA, "default", "dynamic_gene_test.csv")
if (file.exists(dgt_fp)) {
  dgt <- fread(dgt_fp)
  cat("dynamic gene test cached:", nrow(dgt), "genes\n")
} else {
  cat("Running differentialGeneTest (chunked, serial)...\n")
  dgt <- run_dgt_chunked(cds_d, cand, "default")
  setorder(dgt, qval, pval)
  fwrite(dgt, dgt_fp)
}
cat(sprintf("tested: %d | status OK: %d | qval<0.05: %d\n",
            nrow(dgt), sum(dgt$status == "OK"), sum(dgt$status == "OK" & dgt$qval < 0.05)))

dyn <- dgt[status == "OK" & qval < 0.05]$gene_id
cat("dynamic genes:", length(dyn), "\n")

sm_d <- smooth_and_modules(cds_d, dyn, "default")
save_modules(sm_d, "default")

mod_tab <- data.table(gene_id = dyn,
                      gene_short_name = fd[dyn, ]$gene_short_name,
                      num_cells_expressed = fd[dyn, ]$num_cells_expressed,
                      PTG = sm_d$ptg[dyn])
mod_tab[, peak_grid_point := apply(sm_d$z[dyn, , drop = FALSE], 1, which.max)]
mod_tab[, peak_pseudotime := sm_d$grid[peak_grid_point]]
mod_tab[, peak_z := apply(sm_d$z[dyn, , drop = FALSE], 1, max)]
setorder(mod_tab, PTG, -peak_z)
fwrite(mod_tab, file.path(TAB, "default", "s4h_dynamic_genes_modules.csv"))

# ==================================================================== step 2.2
# C00 orientation: verify input identity, spot-check 50 genes, then reuse
cds_c <- readRDS(P_CDS["C00_associated_state"])
id_exprs <- isTRUE(all.equal(as.numeric(exprs(cds_d)[cand[1:100], ]),
                             as.numeric(exprs(cds_c)[cand[1:100], ])))
id_sf <- isTRUE(all.equal(pData(cds_d)$Size_Factor, pData(cds_c)$Size_Factor))
id_pt <- isTRUE(all.equal(pData(cds_d)$Pseudotime, pData(cds_c)$Pseudotime))
cat(sprintf("input identity default vs C00: exprs=%s, Size_Factor=%s, Pseudotime=%s\n",
            id_exprs, id_sf, id_pt))
stopifnot(id_exprs, id_sf, id_pt)
set.seed(42)
spot <- sample(cand, 50)
dgt_spot <- monocle::differentialGeneTest(cds_c[spot, ],
                                          fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)",
                                          reducedModelFormulaStr = "~1",
                                          relative_expr = TRUE, cores = 1)
dgt_spot$gene_id <- rownames(dgt_spot)
cmp <- merge(as.data.table(dgt_spot)[, .(gene_id, pval_c = pval, qval_c = qval)],
             dgt[gene_id %in% spot, .(gene_id, pval_d = pval, qval_d = qval)],
             by = "gene_id")
# identity is judged on pvals (the test statistic); qvals differ by design
# because BH scope differs (50-gene spot set vs the 4,868-gene full set).
# RHOA-level denormal pvals (~1e-308) are treated as equal below 1e-300.
pval_diff <- abs(cmp$pval_c - cmp$pval_d)
pval_diff[cmp$pval_c < 1e-300 & cmp$pval_d < 1e-300] <- 0
qval_consistent <- max(abs(cmp$qval_c -
                             p.adjust(cmp$pval_c, method = "BH"))) < 1e-8
spot_ok <- nrow(cmp) == 50 && max(pval_diff) < 1e-12 * max(cmp$pval_d) &&
  qval_consistent
cat(sprintf("spot check: max |dpval| = %.2e (excl. denormal), qval BH-consistent = %s
",
            max(pval_diff), qval_consistent))
fwrite(cmp, file.path(TAB, "C00_associated_state", "s4h_spot_check_50genes.csv"))
cat("C00 spot check (50 genes) identical:", spot_ok, "\n")
stopifnot(spot_ok)
dgt_c <- copy(dgt)
fwrite(dgt_c, file.path(DATA, "C00_associated_state", "dynamic_gene_test.csv"))
sm_c <- smooth_and_modules(cds_c, dyn, "C00_associated_state")
save_modules(sm_c, "C00_associated_state")
fwrite(mod_tab, file.path(TAB, "C00_associated_state", "s4h_dynamic_genes_modules.csv"))

# PTG intersection mapping between orientations (no id-to-id mapping assumed)
ptg_xmap <- data.table(gene_id = dyn,
                       PTG_default = sm_d$ptg[dyn],
                       PTG_c00 = sm_c$ptg[dyn])
xmap_tab <- ptg_xmap[, .N, by = .(PTG_default, PTG_c00)]
setorder(xmap_tab, PTG_default, PTG_c00)
fwrite(xmap_tab, file.path(TAB, "s4h_PTG_intersection_mapping.csv"))
cat("\nPTG intersection mapping (default x C00):\n"); print(xmap_tab)

# ==================================================================== step 2.3
# S4h global heatmaps: rows = ward.D2 dendrogram order, cols = pseudotime grid
draw_s4h <- function(sm, orient, mod_tab) {
  g_order <- sm$hc$labels[sm$hc$order]
  zplot <- sm$z_disp[g_order, , drop = FALSE]
  ptg_vec <- factor(sm$ptg[g_order], levels = names(PTG_COLS))
  reps <- mod_tab[, head(.SD[order(-peak_z)], 4), by = PTG]
  lab <- setNames(rep("", nrow(zplot)), g_order)
  lab[match(reps$gene_id, g_order)] <-
    unname(setNames(sprintf("%s (%s)", reps$gene_id, reps$gene_short_name),
                    reps$gene_id)[reps$gene_id])

  col_fun <- colorRamp2(c(-3, 0, 3), c("#2166AC", "#F7F7F7", "#B2182B"))
  ha <- rowAnnotation(
    PTG = ptg_vec,
    col = list(PTG = PTG_COLS),
    annotation_legend_param = list(PTG = list(title = "Gene module",
                                              title_gp = gpar(fontsize = 9, fontfamily = "Arial"),
                                              labels_gp = gpar(fontsize = 8, fontfamily = "Arial"))),
    annotation_name_gp = gpar(fontsize = 9, fontfamily = "Arial"))
  ht <- Heatmap(zplot, name = "Z per gene",
                col = col_fun,
                cluster_rows = sm$hc, cluster_columns = FALSE,
                show_row_names = TRUE, row_labels = lab[g_order],
                row_names_side = "right",
                row_names_gp = gpar(fontsize = 6.5, fontfamily = "Arial"),
                show_column_names = FALSE,
                column_title = sprintf("Pseudotime-dynamic genes (q<0.05, n=%d), %s orientation - columns = 100-pt pseudotime grid",
                                       length(g_order), orient),
                column_title_gp = gpar(fontsize = 10, fontfamily = "Arial"),
                left_annotation = ha,
                heatmap_legend_param = list(title = "Z (per gene,\nclipped +/-3)",
                                            title_gp = gpar(fontsize = 9, fontfamily = "Arial"),
                                            labels_gp = gpar(fontsize = 8, fontfamily = "Arial")))
  png(file.path(FIG, orient, "S4h_dynamic_gene_heatmap.png"),
      width = 10.5, height = 14, units = "in", res = 600, family = "Arial", bg = "white")
  draw(ht, merge_legends = TRUE, padding = unit(c(10, 12, 6, 4), "mm"))
  dev.off()
  cairo_pdf(file.path(FIG, orient, "S4h_dynamic_gene_heatmap.pdf"),
            width = 10.5, height = 14, family = "Arial", bg = "white")
  draw(ht, merge_legends = TRUE, padding = unit(c(10, 12, 6, 4), "mm"))
  dev.off()
  fwrite(as.data.table(zplot, keep.rownames = "gene_id"),
         file.path(SRC, orient, "S4h_heatmap_display_source.csv"))
  fwrite(data.table(gene_id = g_order, row_rank = seq_along(g_order),
                    PTG = sm$ptg[g_order]),
         file.path(SRC, orient, "S4h_heatmap_row_order_source.csv"))
  cat("S4h heatmap written:", orient, "\n")
}
draw_s4h(sm_d, "default", mod_tab)
draw_s4h(sm_c, "C00_associated_state", mod_tab)

# ==================================================================== step 2.4
# Branch-resolved panel. Verified geometry (default CDS K-tree): root point
# Y_14, branch point Y_48 (degree 3) at pseudotime >= 19.458. The root-side
# component (6,552 State-1 cells) is the SHARED TRUNK; beyond Y_48 the tree
# splits into two terminal arms (State 2: 8,477 cells; State 3: 1,771 cells).
# Each path = shared trunk + one terminal arm (trunk re-displayed and marked
# by the segment annotation); no extrapolation; coverage tables saved.
build_branches <- function(cds, orient, mod_tab, sm) {
  ao <- cds@auxOrderingData$DDRTree
  closest_idx <- as.integer(ao$pr_graph_cell_proj_closest_vertex[, 1])
  names(closest_idx) <- colnames(cds)
  root_cell <- ao$root_cell
  cds_ref <- readRDS(P_CDS["default"])
  k_tree <- cds_ref@minSpanningTree
  stopifnot(length(igraph::V(k_tree)) == 136)
  node_names <- igraph::V(k_tree)$name
  cell_node <- setNames(node_names[closest_idx], names(closest_idx))
  root_v <- unname(cell_node[root_cell])
  bp <- intersect(ao$branch_points, node_names)
  stopifnot(length(bp) == 1)
  deg <- igraph::degree(k_tree)
  bp_idx <- match(bp, node_names)
  bp_pt <- min(pData(cds)$Pseudotime[closest_idx == bp_idx], na.rm = TRUE)
  cat(sprintf("[%s] root point: %s | branch point: %s (degree %d) | bp pseudotime >= %.3f\n",
              orient, root_v, bp, deg[bp], bp_pt))

  # components after removing the branch point: the one holding the root is
  # the shared trunk; the rest are terminal arms
  g_no_bp <- igraph::delete_vertices(k_tree, bp)
  comps <- igraph::components(g_no_bp)
  comp_of <- comps$membership
  trunk_comp <- comp_of[root_v]
  arm_comps <- setdiff(unique(comp_of), trunk_comp)
  arm_order <- order(sapply(arm_comps, function(cc) {
    nodes_cc <- names(comp_of)[comp_of == cc]
    cells_cc <- names(cell_node)[cell_node %in% nodes_cc]
    median(pData(cds)$Pseudotime[cells_cc])
  }))
  arm_comps <- arm_comps[arm_order]
  node_path <- setNames(rep("trunk_shared", length(node_names)), node_names)
  state_of_arm <- list()
  for (b in seq_along(arm_comps)) {
    nodes_b <- names(comp_of)[comp_of == arm_comps[b]]
    node_path[names(node_path) %in% nodes_b] <- sprintf("arm_%d", b)
    cells_b <- names(cell_node)[cell_node %in% nodes_b]
    state_of_arm[[b]] <- names(sort(table(pData(cds)$State[cells_b]),
                                    decreasing = TRUE))[1]
  }
  # cell-level assignment: trunk component cells AND branch-point cells stay
  # in trunk_shared; terminal-arm cells follow their component
  cell_path <- setNames(unname(node_path[cell_node]), names(cell_node))

  groups <- c("trunk_shared", sprintf("arm_%d", seq_along(arm_comps)))
  cov <- data.table(path = groups)
  cov[, n_cells := sapply(path, function(pp) sum(cell_path == pp))]
  cov[, n_samples := sapply(path, function(pp)
    uniqueN(pData(cds)$SampleID_Global[cell_path == pp]))]
  cov[, states := sapply(path, function(pp)
    paste(sort(unique(pData(cds)$State[cell_path == pp])), collapse = ","))]
  cat("group coverage:\n"); print(cov)
  fwrite(cov, file.path(TAB, orient, "s4h_branch_path_cells.csv"))

  # sample x group coverage (>=20-cell rule flags formal sample support)
  cd <- data.table(CellID = names(cell_path), path = cell_path,
                   SampleID_Global = as.character(pData(cds)$SampleID_Global),
                   Cluster = as.character(pData(cds)$Cluster))
  cov_sp <- cd[, .(n_cells = .N), by = .(path, SampleID_Global)]
  cov_sp[, ge20 := n_cells >= 20]
  setorder(cov_sp, path, -n_cells)
  fwrite(cov_sp, file.path(TAB, orient, "s4h_branch_sample_path_coverage.csv"))
  cov_sc <- cd[, .(n_cells = .N), by = .(SampleID_Global, Cluster)]
  cov_sc[, ge20 := n_cells >= 20]
  fwrite(cov_sc, file.path(TAB, orient, "s4h_branch_sample_cluster_coverage.csv"))

  # one panel per terminal path: shared trunk + that arm (trunk re-displayed,
  # marked by the segment annotation at bp_pt); rows keep the global PTG order
  plots <- list()
  pt_vec <- setNames(pData(cds)$Pseudotime, colnames(cds))
  for (b in seq_along(arm_comps)) {
    pp <- sprintf("arm_%d", b)
    cells_pp <- names(cell_path)[cell_path %in% c("trunk_shared", pp)]
    pt_pp <- unname(pt_vec[cells_pp])
    grid_pp <- seq(min(pt_pp), max(pt_pp), length.out = 100)
    m_pp <- monocle:::genSmoothCurves(cds[dyn, cells_pp], cores = 1,
                                      trend_formula = "~sm.ns(Pseudotime, df=3)",
                                      relative_expr = TRUE,
                                      new_data = data.frame(Pseudotime = grid_pp))
    m_pp <- log10(m_pp + 1)
    z_pp <- Matrix::t(scale(Matrix::t(m_pp), center = TRUE))
    z_pp[!is.finite(z_pp)] <- 0
    z_pp <- pmin(pmax(z_pp, -3), 3)
    g_order <- sm$hc$labels[sm$hc$order]
    z_pp <- z_pp[g_order, , drop = FALSE]
    seg_lab <- ifelse(grid_pp < bp_pt, "shared trunk", "terminal arm")
    ha_col <- HeatmapAnnotation(
      Segment = factor(seg_lab, levels = c("shared trunk", "terminal arm")),
      col = list(Segment = c("shared trunk" = "grey85",
                             "terminal arm" = "#F2C500")),
      annotation_name_gp = gpar(fontsize = 8, fontfamily = "Arial"),
      annotation_legend_param = list(Segment = list(
        title = "Path segment", title_gp = gpar(fontsize = 8, fontfamily = "Arial"),
        labels_gp = gpar(fontsize = 7, fontfamily = "Arial"))))
    ptg_vec <- factor(sm$ptg[g_order], levels = names(PTG_COLS))
    ha <- rowAnnotation(PTG = ptg_vec, col = list(PTG = PTG_COLS),
                        show_legend = FALSE,
                        annotation_name_gp = gpar(fontsize = 8, fontfamily = "Arial"))
    reps <- mod_tab[, head(.SD[order(-peak_z)], 1), by = PTG]  # 1/panel; full lists in CSV
    lab <- setNames(rep("", nrow(z_pp)), g_order)
    lab[match(reps$gene_id, g_order)] <-
      unname(setNames(sprintf("%s (%s)", reps$gene_id, reps$gene_short_name),
                      reps$gene_id)[reps$gene_id])
    ht <- Heatmap(z_pp, name = "Z per gene",
                  col = colorRamp2(c(-3, 0, 3), c("#2166AC", "#F7F7F7", "#B2182B")),
                  cluster_rows = FALSE, cluster_columns = FALSE,
                  show_row_names = TRUE, row_labels = lab[g_order],
                  row_names_side = "right",
                  row_names_gp = gpar(fontsize = 6, fontfamily = "Arial"),
                  show_column_names = FALSE,
                  top_annotation = ha_col,
                  left_annotation = ha,
                  column_title = sprintf("%s - path %d = shared trunk + %s (State %s; n=%d cells, %d samples); columns = pseudotime grid",
                                         orient, b, pp, state_of_arm[[b]],
                                         length(cells_pp),
                                         uniqueN(cd[path %in% c("trunk_shared", pp)]$SampleID_Global)),
                  column_title_gp = gpar(fontsize = 9.5, fontfamily = "Arial"),
                  heatmap_legend_param = list(title = "Z (per gene,\nclipped +/-3)",
                                              title_gp = gpar(fontsize = 8, fontfamily = "Arial"),
                                              labels_gp = gpar(fontsize = 7, fontfamily = "Arial")))
    plots[[pp]] <- ht
    fwrite(as.data.table(z_pp, keep.rownames = "gene_id"),
           file.path(SRC, orient, sprintf("S4h_path_%s_display_source.csv", pp)))
    fwrite(data.table(gene_id = g_order, PTG = sm$ptg[g_order]),
           file.path(SRC, orient, sprintf("S4h_path_%s_row_order_source.csv", pp)))
  }
  png(file.path(FIG, orient, "S4h_branch_resolved_heatmaps.png"),
      width = 12, height = 14, units = "in", res = 600, family = "Arial", bg = "white")
  draw(Reduce(`%v%`, plots), merge_legends = TRUE,
       padding = unit(c(10, 12, 6, 4), "mm"))
  dev.off()
  cairo_pdf(file.path(FIG, orient, "S4h_branch_resolved_heatmaps.pdf"),
            width = 12, height = 14, family = "Arial", bg = "white")
  draw(Reduce(`%v%`, plots), merge_legends = TRUE,
       padding = unit(c(10, 12, 6, 4), "mm"))
  dev.off()
  fwrite(data.table(path = names(plots),
                    composition = sprintf("shared trunk + %s", names(plots)),
                    state = unlist(state_of_arm)),
         file.path(TAB, orient, "s4h_branch_paths_used.csv"))
  cat("branch heatmaps written:", orient, "\n")
  invisible(list(cov = cov, bp = bp, bp_pt = bp_pt,
                 state_of_arm = state_of_arm))
}

br_d <- build_branches(cds_d, "default", mod_tab, sm_d)
br_c <- build_branches(cds_c, "C00_associated_state", mod_tab, sm_c)

capture.output(sessionInfo(), file = file.path(LOG, "sessionInfo_s4ghi_02.txt"))
cat(sprintf("\nStep 2 done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
