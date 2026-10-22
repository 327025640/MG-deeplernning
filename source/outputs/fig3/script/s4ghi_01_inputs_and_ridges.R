# =============================================================================
# S4g-i Step 1: inputs + S4g ridge plots (both orientations)
# Contract v1 (user-confirmed 2026-09-05):
#   x = Pseudotime; horizontal viridis gradient fills (no rainbow); rows =
#   frozen cluster order subset (14); labels "Cxx MG_*"; density grids cached;
#   descriptive only (balanced sampling of 1200 cells/cluster - ridge height is
#   NOT atlas abundance). Both orientations computed from their own CDS.
# Output: trajectory_S4/{data,tables,figures,source_data}/<orientation>
# =============================================================================
suppressPackageStartupMessages({
  library(monocle)
  library(data.table)
})
options(stringsAsFactors = FALSE)
set.seed(42)

OUT  <- "E:/micro_output/fig3/trajectory_S4"
DATA <- file.path(OUT, "data"); TAB <- file.path(OUT, "tables")
FIG  <- file.path(OUT, "figures"); SRC <- file.path(OUT, "source_data")
for (o in c("default", "C00_associated_state")) {
  for (d in c(file.path(DATA, o), file.path(TAB, o),
              file.path(FIG, o), file.path(SRC, o)))
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
P_MAP <- "E:/micro_output/fig3/fig4/tables/cluster_name_mapping.csv"
P_ORD <- "E:/micro_output/fig3/data/05_cluster_row_order.csv"

map <- fread(P_MAP, encoding = "UTF-8")
ord <- fread(P_ORD)
TRAJ_CL <- c("C00", "C01", "C02", "C03", "C04", "C05", "C06", "C07",
             "C08", "C09", "C12", "C13", "C17", "C19")
frozen14 <- ord[Cluster %in% TRAJ_CL][order(row_order)]$Cluster
stopifnot(length(frozen14) == 14)
lab14 <- setNames(sprintf("%s %s", map$Cluster, map$FinalName), map$Cluster)[frozen14]

t0 <- Sys.time()
build_orientation <- function(cds_path, orient) {
  cat("\n=====", orient, "=====\n")
  cds <- readRDS(cds_path)
  pd <- as.data.table(pData(cds))
  stopifnot(nrow(pd) == 16800, setequal(unique(pd$Cluster), TRAJ_CL))
  pd[, FinalName := map$FinalName[match(Cluster, map$Cluster)]]
  setorder(pd, CellID)
  fwrite(pd[, .(CellID, Cluster, FinalName, Disease, SampleID_Global, DatasetID,
                State, Pseudotime)],
         file.path(DATA, orient, "cells_with_pseudotime.csv"))

  st <- pd[, .(n_cells = .N,
               n_samples = uniqueN(SampleID_Global),
               pseudotime_median = median(Pseudotime),
               pseudotime_q25 = quantile(Pseudotime, .25),
               pseudotime_q75 = quantile(Pseudotime, .75),
               pseudotime_min = min(Pseudotime),
               pseudotime_max = max(Pseudotime)),
           by = .(Cluster, FinalName)]
  st[, row_order := ord$row_order[match(Cluster, ord$Cluster)]]
  setorder(st, row_order)
  fwrite(st, file.path(TAB, orient, "s4g_cluster_stats.csv"))

  # density grid with common bandwidth + common grid, per-ridge max-normalised
  grid_max <- max(pd$Pseudotime)
  bws <- sapply(frozen14, function(cl) bw.nrd0(pd[Cluster == cl]$Pseudotime))
  bw_common <- mean(bws)
  den <- rbindlist(lapply(frozen14, function(cl) {
    d <- density(pd[Cluster == cl]$Pseudotime, bw = bw_common, n = 512,
                 from = 0, to = grid_max)
    data.table(Cluster = cl, x = d$x, dens_scaled = d$y / max(d$y))
  }))
  fwrite(den, file.path(DATA, orient, "s4g_density_grid.csv"))
  cat("cells/stats/density written | common bw =", round(bw_common, 4), "\n")

  # ---------------- ridge plot: base graphics, Arial, viridis x-gradient -----
  # gradient fills are the x-axis pseudotime encoding (contract: viridis,
  # continuous, no rainbow); ridge heights equal (balanced sampling)
  vir <- function(n) grDevices::hcl.colors(n, palette = "viridis")
  ridge_h <- 0.85
  ncl <- length(frozen14)
  row_of <- setNames(seq_len(ncl), frozen14)   # 1 = bottom
  xs_seg <- seq(0, grid_max, length.out = 121) # 120 fill slices
  seg_cols <- vir(120)

  png(file.path(FIG, orient, "S4g_ridge_pseudotime.png"),
      width = 11, height = 7.5, units = "in", res = 600,
      family = "Arial", bg = "white")
  par(mar = c(4.2, 10.8, 3.4, 1.2), family = "Arial", las = 1)
  plot(NULL, xlim = c(0, grid_max * 1.01), ylim = c(0.35, ncl + 1.25),
       xlab = "Pseudotime", ylab = NA, yaxt = "n", xaxt = "n",
       main = "Pan-disease microglia subclusters along pseudotime (Monocle2)", cex.main = 0.95)
  axis(1, at = pretty(c(0, grid_max)), cex.axis = 0.9)
  axis(2, at = seq_len(ncl), labels = lab14[frozen14], tick = FALSE, cex.axis = 0.85)
  # draw top ridge last so upper ridges sit in front (ggridges look)
  for (cl in rev(frozen14)) {
    d <- den[Cluster == cl][order(x)]
    yb <- row_of[cl]
    ys <- yb + approx(d$x, d$dens_scaled, xout = xs_seg)$y * ridge_h
    for (k in seq_len(length(xs_seg) - 1)) {
      polygon(x = c(xs_seg[k], xs_seg[k + 1], xs_seg[k + 1], xs_seg[k]),
              y = c(yb, yb, ys[k + 1], ys[k]),
              col = adjustcolor(seg_cols[k], alpha.f = 0.85), border = NA)
    }
    lines(xs_seg, ys, col = "grey25", lwd = 0.6)
    abline(h = yb, col = "grey80", lwd = 0.5)
    # median tick
    med <- st[Cluster == cl]$pseudotime_median
    ytop <- yb + approx(d$x, d$dens_scaled, xout = med)$y * ridge_h
    segments(med, yb, med, ytop, col = "black", lwd = 1.4)
  }
  mtext("Density per cluster, equal height (balanced sampling; height is not atlas abundance)",
        side = 3, line = 0.3, cex = 0.72, col = "grey35")
  legend("topright", bty = "n", cex = 0.78, lwd = c(1.4, NA), pch = c(NA, 15),
         col = c("black", seg_cols[60]),
         legend = c("cluster median", "pseudotime gradient (early \u2192 late)"))
  dev.off()

  # source CSV: the plotted density grid + medians
  fwrite(den, file.path(SRC, orient, "S4g_ridge_density_source.csv"))
  fwrite(st, file.path(SRC, orient, "S4g_ridge_stats_source.csv"))

  cairo_pdf(file.path(FIG, orient, "S4g_ridge_pseudotime.pdf"),
            width = 11, height = 7.5, family = "Arial", bg = "white")
  par(mar = c(4.2, 10.8, 3.4, 1.2), family = "Arial", las = 1)
  plot(NULL, xlim = c(0, grid_max * 1.01), ylim = c(0.35, ncl + 1.25),
       xlab = "Pseudotime", ylab = NA, yaxt = "n", xaxt = "n",
       main = "Pan-disease microglia subclusters along pseudotime (Monocle2)", cex.main = 0.95)
  axis(1, at = pretty(c(0, grid_max)), cex.axis = 0.9)
  axis(2, at = seq_len(ncl), labels = lab14[frozen14], tick = FALSE, cex.axis = 0.85)
  for (cl in rev(frozen14)) {
    d <- den[Cluster == cl][order(x)]
    yb <- row_of[cl]
    ys <- yb + approx(d$x, d$dens_scaled, xout = xs_seg)$y * ridge_h
    for (k in seq_len(length(xs_seg) - 1)) {
      polygon(x = c(xs_seg[k], xs_seg[k + 1], xs_seg[k + 1], xs_seg[k]),
              y = c(yb, yb, ys[k + 1], ys[k]),
              col = adjustcolor(seg_cols[k], alpha.f = 0.85), border = NA)
    }
    lines(xs_seg, ys, col = "grey25", lwd = 0.6)
    abline(h = yb, col = "grey80", lwd = 0.5)
    med <- st[Cluster == cl]$pseudotime_median
    ytop <- yb + approx(d$x, d$dens_scaled, xout = med)$y * ridge_h
    segments(med, yb, med, ytop, col = "black", lwd = 1.4)
  }
  mtext("Density per cluster, equal height (balanced sampling; height is not atlas abundance)",
        side = 3, line = 0.3, cex = 0.72, col = "grey35")
  legend("topright", bty = "n", cex = 0.78, lwd = c(1.4, NA), pch = c(NA, 15),
         col = c("black", seg_cols[60]),
         legend = c("cluster median", "pseudotime gradient (early \u2192 late)"))
  dev.off()
  cat("S4g figure written\n")
  invisible(list(pd = pd, st = st, den = den, bw_common = bw_common))
}

r1 <- build_orientation("E:/micro_output/fig3/data/Fig2i_monocle2_cds_default.rds",
                        "default")
r2 <- build_orientation("E:/micro_output/fig3/data/Fig2i_monocle2_cds_C00_root.rds",
                        "C00_associated_state")

# recorded equivalence check between the two stored orientations (audit found
# identical Pseudotime; keep the fact in the delivered tables)
eq <- r1$pd[r2$pd, on = "CellID"][, .(
  n = .N,
  identical_pseudotime = isTRUE(all.equal(Pseudotime, i.Pseudotime, tolerance = 1e-10)),
  max_abs_diff = max(abs(Pseudotime - i.Pseudotime)))]
fwrite(eq, file.path(TAB, "orientation_equivalence_check.csv"))
cat("\norientation equivalence:\n"); print(eq)

capture.output(sessionInfo(), file = file.path(OUT, "logs", "sessionInfo_s4ghi_01.txt"))
cat(sprintf("\nStep 1 done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
