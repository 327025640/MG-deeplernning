# =============================================================================
# 03_consensus.R — 正式 NMF: 02_consensus 程序聚类
# 输入: 01_bootstrap/draw_{001..100}.rds
# 程序池构建: 每 (draw, K) 的 20 inits 对齐 (以 init1 为参照, 其余列按 cosine
#   最佳匹配重排) -> 逐基因 median consensus -> K 个程序
#   => 程序池 = 100 x 12 x K = 11,400 程序 (全 bootstrap 汇总, 非 per-init 全量)
# 相似度: top-50 Jaccard (主) + 完整 loading cosine (验证)
# 聚类: hclust(average) + 3% 规则 (minSize=max(2,floor(n*0.03)), deepSplit=2,
#   CLUSTER_CAP=30); 簇 >=3 程序有效
# 每簇签名: L2 normalize -> 簇内逐基因 median -> L1 normalize
# 输出: 02_consensus/{consensus_selection_diagnostics.csv,
#       consensus_programs.rds, consensus_programs_top50.csv,
#       consensus_assignment.csv, jaccard_similarity.csv, cosine_similarity.csv,
#       selected_K_diagnostics.{png,pdf}}
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(dynamicTreeCut)
})

PROJ     <- "D:/data/microglia_new"
NMF_ROOT <- "E:/micro_output/v3_dims16/05_nmf"
B1_DIR   <- file.path(NMF_ROOT, "01_bootstrap")
C2_DIR   <- file.path(NMF_ROOT, "02_consensus")
dir.create(C2_DIR, recursive = TRUE, showWarnings = FALSE)

K_RANGE <- 4:15
N_INIT  <- 20L

log_msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...), "\n")
log_msg("=== 03 consensus start ===")

draw_files <- file.path(B1_DIR, sprintf("draw_%03d.rds", 1:100))
stopifnot(all(file.exists(draw_files)), all(file.exists(
  file.path(B1_DIR, sprintf("draw_%03d_DONE.txt", 1:100)))))

# ============================================================================
# Step 1: 程序池 (每 draw x K 的 20-inits median consensus)
# ============================================================================
all_W <- list()        # 程序池: 5000 x n_programs
prog_meta <- list()
n_prog <- 0L
genes_ref <- NULL

for (i in 1:100) {
  d <- readRDS(draw_files[i])
  if (is.null(genes_ref)) genes_ref <- d$genes
  stopifnot(identical(d$genes, genes_ref))
  for (k in as.character(K_RANGE)) {
    W_k <- d$Ws[[k]]                     # 5000 x (20*K)
    if (is.null(W_k)) next
    K <- as.integer(k)
    # 对齐: 以 init1 的 K 列为参照, 其余 init 列按 cosine 最佳匹配重排
    ref <- W_k[, 1:K, drop = FALSE]
    lr <- sqrt(colSums(ref^2)) + 1e-12
    ref_n <- sweep(ref, 2, lr, "/")
    cols <- vector("list", N_INIT)
    cols[[1]] <- ref
    for (j in 2:N_INIT) {
      Wj <- W_k[, ((j - 1) * K + 1):(j * K), drop = FALSE]
      lj <- sqrt(colSums(Wj^2)) + 1e-12
      C <- crossprod(ref_n, sweep(Wj, 2, lj, "/"))     # K x K cosine
      # 贪心匹配: 每 ref 列选未用且余弦最大的 Wj 列
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
    # 逐基因 median consensus
    W_cons <- apply(abind::abind(cols, along = 3), c(1, 2), stats::median)
    all_W[[length(all_W) + 1]] <- W_cons
    prog_meta[[length(prog_meta) + 1]] <- data.frame(
      draw = i, K = K, n_init = N_INIT,
      program_in_draw = seq_len(K),
      program_id = sprintf("d%03d_k%02d_p%02d", i, K, seq_len(K)),
      stringsAsFactors = FALSE)
  }
  if (i %% 10 == 0) log_msg("draw", i, "/100 loaded")
}
W_all <- do.call(cbind, all_W)
rownames(W_all) <- genes_ref              # 必须设置行名, top50 提取依赖
prog_meta <- rbindlist(prog_meta)
n_programs <- ncol(W_all)
log_msg("程序池: ", n_programs, " 程序 (100 draws x 12 K x K)")

# ============================================================================
# Step 2: 相似度 (Jaccard top50 主 + cosine 验证)
# ============================================================================
top50_list <- lapply(seq_len(n_programs), function(p) {
  w <- W_all[, p]
  names(sort(w, decreasing = TRUE))[1:50]
})
log_msg("Jaccard 相似度 (", n_programs, "^2 对)...")
jaccard_mat <- matrix(0, n_programs, n_programs)
for (i in seq_len(n_programs)) {
  gi <- top50_list[[i]]
  for (j in seq_len(i)) {
    inter <- length(intersect(gi, top50_list[[j]]))
    jaccard_mat[i, j] <- inter / (100 - inter)
    jaccard_mat[j, i] <- jaccard_mat[i, j]
  }
  if (i %% 1000 == 0) log_msg("  Jaccard ", i, "/", n_programs)
}
jd <- as.data.frame(jaccard_mat)
jd$program_id <- prog_meta$program_id
jd <- jd[, c("program_id", setdiff(names(jd), "program_id"))]
fwrite(jd, file.path(C2_DIR, "jaccard_similarity.csv"))

log_msg("Cosine 相似度 (验证)...")
Wn <- sweep(W_all, 2, sqrt(colSums(W_all^2)) + 1e-12, "/")
cosine_mat <- crossprod(Wn)
cosine_mat[cosine_mat > 1] <- 1
cd <- as.data.frame(cosine_mat)
cd$program_id <- prog_meta$program_id
cd <- cd[, c("program_id", setdiff(names(cd), "program_id"))]
fwrite(cd, file.path(C2_DIR, "cosine_similarity.csv"))

# ============================================================================
# Step 3: 层次聚类 (3% 规则)
# ============================================================================
jaccard_dist <- as.dist(1 - jaccard_mat)
hc <- hclust(jaccard_dist, method = "average")
min_size <- max(2L, floor(n_programs * 0.03))
labels <- cutreeDynamic(hc, distM = as.matrix(jaccard_dist),
                        minClusterSize = min_size, method = "hybrid",
                        deepSplit = 2, verbose = 0)
if (all(labels == 0) || length(unique(labels)) <= 1 || length(unique(labels)) > 30) {
  h <- quantile(hc$height, 0.75)
  labels <- cutree(hc, h = h)
  log_msg("动态剪枝退化 -> 固定高度 h=", round(h, 4))
}
names(labels) <- prog_meta$program_id
sizes <- table(labels)
valid <- as.integer(names(sizes[sizes >= 3]))
orphan <- as.integer(names(sizes[sizes < 3]))
n_cons <- length(valid)
log_msg("簇: ", length(unique(labels)), " (min size=", min_size,
        "); 有效 (>=3): ", n_cons, "; orphan: ", length(orphan))

# ============================================================================
# Step 4: 每簇 consensus 签名 (L2 -> median -> L1)
# ============================================================================
old2new <- setNames(paste0("MP", sprintf("%02d", seq_along(valid))), as.character(valid))
assign <- data.table(
  program_id = prog_meta$program_id,
  cluster_raw = labels[prog_meta$program_id],
  cluster_final = ifelse(labels[prog_meta$program_id] %in% valid,
                         old2new[as.character(labels[prog_meta$program_id])], "orphan"),
  stringsAsFactors = FALSE)
assign <- merge(assign, prog_meta[, .(program_id, draw, K)],
                by = "program_id", all.x = TRUE)
fwrite(assign, file.path(C2_DIR, "consensus_assignment.csv"))

cons_W <- list(); cons_top50 <- list()
for (cid in valid) {
  cname <- old2new[as.character(cid)]
  members <- assign[cluster_raw == cid, program_id]
  idx <- which(prog_meta$program_id %in% members)
  Wm <- W_all[, idx, drop = FALSE]
  # L2 normalize 每程序 -> 逐基因 median -> L1 normalize
  Wn_m <- sweep(Wm, 2, sqrt(colSums(Wm^2)) + 1e-12, "/")
  w_med <- apply(Wn_m, 1, stats::median)
  w_l1 <- w_med / sum(w_med)
  cons_W[[cname]] <- w_l1
  ord <- order(w_l1, decreasing = TRUE)
  cons_top50[[cname]] <- data.frame(
    MP = cname, rank = seq_len(50),
    gene = genes_ref[ord[1:50]], weight = w_l1[ord[1:50]],
    stringsAsFactors = FALSE)
  log_msg("  ", cname, ": ", length(members), " programs")
}
consensus_W <- do.call(cbind, cons_W)
rownames(consensus_W) <- genes_ref
MOD <- colnames(consensus_W)
fwrite(as.data.frame(consensus_W), file.path(C2_DIR, "consensus_programs.csv"),
       row.names = TRUE)
fwrite(rbindlist(cons_top50), file.path(C2_DIR, "consensus_programs_top50.csv"))
saveRDS(list(consensus_W = consensus_W, assignment = assign,
             jaccard = jaccard_mat, cosine = cosine_mat, hc = hc,
             labels = labels, valid_clusters = valid,
             params = list(n_programs = n_programs, n_consensus = n_cons,
                           min_size = min_size, K_range = K_RANGE,
                           n_init = N_INIT)),
        file.path(C2_DIR, "consensus_programs.rds"))

# ============================================================================
# Step 5: K 选择诊断 (从 k_scan_all_draws.csv)
# ============================================================================
k_scan <- fread(file.path(B1_DIR, "k_scan_all_draws.csv"))
diag <- k_scan[, .(mse = mean(mse), stability = mean(stability),
                   mse_sd = sd(mse), stability_sd = sd(stability)),
               by = K][order(K)]
fwrite(diag, file.path(C2_DIR, "consensus_selection_diagnostics.csv"))
log_msg("K 诊断:"); print(diag)

if (capabilities("cairo")) {
  pdf(file.path(C2_DIR, "selected_K_diagnostics.pdf"), width = 8, height = 6)
  par(mar = c(5, 4, 4, 4))
  plot(diag$K, diag$stability, type = "b", col = "#3182BD", lwd = 2,
       xlab = "K", ylab = "stability (init consensus)", main = "K selection diagnostics")
  par(new = TRUE)
  plot(diag$K, diag$mse, type = "b", col = "#D24B40", lwd = 2, axes = FALSE,
       xlab = NA, ylab = NA)
  axis(4); mtext("MSE", side = 4, line = 3, col = "#D24B40")
  legend("topright", legend = c("stability", "MSE"), col = c("#3182BD", "#D24B40"),
         lwd = 2)
  dev.off()
  png(file.path(C2_DIR, "selected_K_diagnostics.png"), width = 1600, height = 1200,
      res = 200, type = "cairo-png")
  par(mar = c(5, 4, 4, 4))
  plot(diag$K, diag$stability, type = "b", col = "#3182BD", lwd = 2,
       xlab = "K", ylab = "stability (init consensus)", main = "K selection diagnostics")
  par(new = TRUE)
  plot(diag$K, diag$mse, type = "b", col = "#D24B40", lwd = 2, axes = FALSE,
       xlab = NA, ylab = NA)
  axis(4); mtext("MSE", side = 4, line = 3, col = "#D24B40")
  legend("topright", legend = c("stability", "MSE"), col = c("#3182BD", "#D24B40"),
         lwd = 2)
  dev.off()
}

log_msg("== 03 done: ", n_cons, " consensus programs (MP) ==")
