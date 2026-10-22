# =============================================================================
# NMFVIS_02_gene_stability_and_summary.R — 由既有 bootstrap draw 结果事后聚合：
#   1) gene_stability（每 MP：5000 基因跨 bootstrap consensus loadings 的 CV）
#   2) per-MP summary（n_programs / n_draws）
#   3) 相似度矩阵的排序长表索引（供 PanelA；全矩阵直接读 CSV）
# 不重新运行 NMF；对齐/consensus 逻辑与 13_nmf_consensus.R 冻结版逐行一致。
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(abind) })

NV     <- "E:/micro_output/v3_dims16/05_nmf"
ROOT   <- "E:/micro_output/NMF"
B1_DIR <- file.path(NV, "01_bootstrap")
dir.create(file.path(ROOT, "tables"), recursive = TRUE, showWarnings = FALSE)
K_RANGE <- 4:15
N_INIT  <- 20L

log_msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...), "\n")
log_msg("=== NMFVIS_02 gene stability + summary start ===")

a <- fread(file.path(NV, "02_consensus", "consensus_assignment.csv"))
a[, program_in_draw := as.integer(sub(".*_p", "", program_id))]
stopifnot(nrow(a) == 11400L)

# ---- 程序池重建：每 draw x K 的 20-inits 对齐 median consensus（与 13 逐行一致）----
# 只保留被 assignment 标记为各 MP 的 consensus 列，按 MP 分桶收集
mp_vectors <- setNames(vector("list", 12), sprintf("MP%02d", 1:12))
genes_ref <- NULL

for (i in 1:100) {
  d <- readRDS(file.path(B1_DIR, sprintf("draw_%03d.rds", i)))
  if (is.null(genes_ref)) genes_ref <- d$genes
  stopifnot(identical(d$genes, genes_ref))
  a_k <- a[draw == i]
  for (k in as.character(K_RANGE)) {
    W_k <- d$Ws[[k]]
    if (is.null(W_k)) next
    K <- as.integer(k)
    ref <- W_k[, 1:K, drop = FALSE]
    lr <- sqrt(colSums(ref^2)) + 1e-12
    ref_n <- sweep(ref, 2, lr, "/")
    cols <- vector("list", N_INIT)
    cols[[1]] <- ref
    for (j in 2:N_INIT) {
      Wj <- W_k[, ((j - 1) * K + 1):(j * K), drop = FALSE]
      lj <- sqrt(colSums(Wj^2)) + 1e-12
      C <- crossprod(ref_n, sweep(Wj, 2, lj, "/"))
      used <- rep(FALSE, K)
      perm <- integer(K)
      for (p in 1:K) {
        cand <- which(!used)
        best <- cand[which.max(C[p, cand])]
        perm[p] <- best
        used[best] <- TRUE
      }
      cols[[j]] <- Wj[, perm, drop = FALSE]
    }
    W_cons <- apply(abind::abind(cols, along = 3), c(1, 2), stats::median)
    sub <- a_k[K == k]
    for (r in seq_len(nrow(sub))) {
      mp <- sub$cluster_final[r]
      p <- sub$program_in_draw[r]
      v <- W_cons[, p]
      v <- v / (sqrt(sum(v^2)) + 1e-12)   # L2 标准化, 消除 NMF 程序尺度差异（与 consensus 流程一致）
      mp_vectors[[mp]] <- c(mp_vectors[[mp]], list(v))
    }
  }
  if (i %% 10 == 0) log_msg("draw ", i, "/100 pooled")
}

# ---- per-MP gene stability（mean/sd/CV per gene across member consensus vectors）----
stab_list <- list()
for (mp in names(mp_vectors)) {
  M <- do.call(cbind, mp_vectors[[mp]])       # 5000 x n_members
  mn <- rowMeans(M)
  sdv <- apply(M, 1, stats::sd)
  cv <- ifelse(mn > 1e-12, sdv / mn, NA_real_)
  stab_list[[mp]] <- data.table(MP = mp, gene = genes_ref,
                                mean_loading = mn, sd_loading = sdv, CV = cv)
  log_msg("  ", mp, ": ", ncol(M), " member programs")
}
stab <- rbindlist(stab_list)
fwrite(stab, file.path(ROOT, "tables", "MP_gene_stability.csv"))

# ---- per-MP summary ----
summ <- a[, .(n_programs = .N, n_draws = uniqueN(draw)), by = cluster_final]
setnames(summ, "cluster_final", "MP")
setorder(summ, MP)
gs_n <- stab[, .(n_genes = .N, median_CV = median(CV, na.rm = TRUE)), by = MP]
summ <- merge(summ, gs_n, by = "MP", all.x = TRUE)
summ[, draw_support := n_draws / 100]
fwrite(summ, file.path(ROOT, "tables", "MP_summary.csv"))
log_msg("== NMFVIS_02 done: gene stability + MP summary written ==")
