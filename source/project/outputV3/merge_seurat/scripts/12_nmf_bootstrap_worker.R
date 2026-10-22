# =============================================================================
# 02b_bootstrap_run.R — 正式 NMF: 每 draw K 扫描 (RcppML::nmf)
# 用法: Rscript 02b_bootstrap_run.R --draw-start 1 --draw-end 10
# 参数 (用户冻结):
#   K_RANGE = 4:15, N_INIT = 20 (每 K 固定初始化),
#   seed_init = 20260801 + i*100 + K*1000 + j
#   RcppML::nmf(x, k, tol=1e-4, maxit=100)
# 每 draw 每 K 保存: MSE (||x - wdh||^2 / nnz), 初始化稳定性 (20 W 两两 cosine
#   最佳匹配均值), 每 init 每程序 top-50 gene, 程序权重 (W 拼接 5000 x 20K)
# 输出: 01_bootstrap/draw_{001..100}.rds + DONE.txt
#       k_scan_all_draws.csv + raw_program_manifest.csv (append, 每 draw 原子写)
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  idx <- match(key, args)
  if (is.na(idx)) return(default)
  args[idx + 1]
}
DRAW_START <- as.integer(get_arg("--draw-start"))
DRAW_END   <- as.integer(get_arg("--draw-end"))
stopifnot(!is.na(DRAW_START), !is.na(DRAW_END))

PROJ     <- "D:/data/microglia_new"
NMF_ROOT <- "E:/micro_output/v3_dims16/05_nmf"
SAMP_DIR <- file.path(NMF_ROOT, "01_bootstrap", "sampling")
B1_DIR   <- file.path(NMF_ROOT, "01_bootstrap")
dir.create(B1_DIR, recursive = TRUE, showWarnings = FALSE)

K_RANGE  <- 4:15
N_INIT   <- 20L
BASE_SEED <- 20260801L
NMF_TOL  <- 1e-4
NMF_MAXIT <- 100L

log_msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...), "\n")

# ---- MSE (稀疏高效): ||x - wdh||^2 = ||x||^2 - 2*sum_nz(x*(wdh)) + ||wdh||^2 ----
calc_mse <- function(x, fit) {
  WH <- fit$w %*% (fit$d * fit$h)      # cells x genes (稠密)
  x_nz <- x@x
  row_i <- x@i + 1L                    # CSC 行索引 (1-based)
  col_i <- rep(seq_len(ncol(x)) - 1L, diff(x@p))
  wh_nz <- as.vector(WH)[row_i + nrow(WH) * col_i]
  s2 <- sum(WH^2)
  x2 <- sum(x_nz^2)
  (x2 - 2 * sum(x_nz * wh_nz) + s2) / length(x_nz)
}

# ---- 初始化稳定性: 20 个 W 两两 cosine 最佳匹配均值 ----
init_stability <- function(W_list) {
  n <- length(W_list)
  sims <- numeric(0)
  for (a in seq_len(n - 1)) {
    Wa <- W_list[[a]]
    for (b in (a + 1):n) {
      Wb <- W_list[[b]]
      la <- sqrt(colSums(Wa^2)) + 1e-12
      lb <- sqrt(colSums(Wb^2)) + 1e-12
      C <- crossprod(sweep(Wa, 2, la, "/"), sweep(Wb, 2, lb, "/"))  # K x K
      sims <- c(sims, mean(apply(C, 1, max)))   # 每程序最佳匹配, 平均
    }
  }
  mean(sims)
}

for (i in DRAW_START:DRAW_END) {
  sdir <- file.path(SAMP_DIR, sprintf("draw_%03d", i))
  out_rds <- file.path(B1_DIR, sprintf("draw_%03d.rds", i))
  done_txt <- file.path(B1_DIR, sprintf("draw_%03d_DONE.txt", i))
  if (file.exists(done_txt)) { log_msg("draw", i, " SKIP (DONE)"); next }
  stopifnot(file.exists(file.path(sdir, "DONE.txt")))

  log_msg("draw", i, " loading matrix...")
  x <- Matrix::readMM(file.path(sdir, "expression_cells_by_genes.mtx"))
  x <- as(x, "dgCMatrix")
  genes <- read.csv(file.path(sdir, "genes.csv"), check.names = FALSE)$gene
  stopifnot(ncol(x) == length(genes))
  nnz <- length(x@x)
  log_msg("  matrix: ", nrow(x), " x ", ncol(x), " nnz=", nnz)

  Ws <- list()          # K -> 5000 x (20*K) 程序权重 (来自 fit$h 转置)
  metrics <- list()     # 每 K: mse, stability
  top50 <- list()       # K -> list(init -> data.frame(program, rank, gene))
  for (k in K_RANGE) {
    W_k <- matrix(0, nrow = ncol(x), ncol = N_INIT * k)   # genes x 20K
    fits <- vector("list", N_INIT)
    t0 <- Sys.time()
    for (j in seq_len(N_INIT)) {
      set.seed(BASE_SEED + i * 100L + k * 1000L + j)   # 用户冻结 seed 规则
      fit <- RcppML::nmf(x, k = k, tol = NMF_TOL, maxit = NMF_MAXIT, verbose = FALSE)
      fits[[j]] <- fit
      W_k[, ((j - 1) * k + 1):(j * k)] <- t(fit$h)   # h: K x genes -> genes x K
    }
    mse <- calc_mse(x, fits[[1]])
    stab <- init_stability(lapply(fits, function(f) t(f$h)))  # 程序权重 genes x K
    metrics[[as.character(k)]] <- data.frame(
      draw = i, K = k, n_init = N_INIT, mse = mse, stability = stab,
      nnz = nnz, n_cells = nrow(x), n_genes = ncol(x),
      elapsed_s = round(as.numeric(Sys.time() - t0, units = "secs"), 1),
      stringsAsFactors = FALSE)
    Ws[[as.character(k)]] <- W_k
    # top50 per init per program
    t50 <- rbindlist(lapply(seq_len(N_INIT), function(j) {
      Wj <- t(fits[[j]]$h)   # genes x K
      do.call(rbind, lapply(seq_len(k), function(p) {
        ord <- order(Wj[, p], decreasing = TRUE)[1:50]
        data.frame(draw = i, K = k, init = j, program = p,
                   seed = BASE_SEED + i * 100L + k * 1000L + j,
                   gene = genes[ord], loading = Wj[ord, p], stringsAsFactors = FALSE)
      }))
    }))
    top50[[as.character(k)]] <- t50
    log_msg("  draw", i, " K=", k, " mse=", round(mse, 5),
            " stab=", round(stab, 4), " (", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s)")
    rm(fits); gc()
  }

  saveRDS(list(
    draw = i, genes = genes,
    Ws = Ws,                       # K -> 5000 x 20K
    metrics = rbindlist(metrics),
    top50 = rbindlist(top50)
  ), out_rds)
  # 原子追加汇总
  met <- rbindlist(metrics)
  kfile <- file.path(B1_DIR, "k_scan_all_draws.csv")
  write.table(met, kfile, sep = ",", row.names = FALSE, col.names = !file.exists(kfile),
              append = file.exists(kfile), quote = TRUE)
  mfile <- file.path(B1_DIR, "raw_program_manifest.csv")
  rp <- rbindlist(top50)
  rp <- rp[, .(n_top50 = .N), by = .(draw, K, init, program, seed)]
  write.table(rp, mfile, sep = ",", row.names = FALSE,
              col.names = !file.exists(mfile), append = file.exists(mfile), quote = TRUE)
  writeLines("DONE", done_txt)
  log_msg("draw", i, " DONE -> ", out_rds)
  rm(x, Ws, metrics, top50); gc()
}
log_msg("== 02b done (draws ", DRAW_START, "-", DRAW_END, ") ==")
