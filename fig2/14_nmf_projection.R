# =============================================================================
# 04_projection.R — 正式 NMF: 03_projection 共识 MP 投影全量细胞
# 方法: RcppML::project(A = t(X), w = consensus_W) -> h (MP x cells)
#       (RcppML 非负最小二乘投影, 每细胞独立求解)
# 验证: 行名=MP, 列名=CellID, CellID 与对象一致, 无重复, 无 NA/Inf/负值,
#       每 MP 有效 usage, source/projection CellID 顺序一致
# 输出: 03_projection/primary_MP_v2_usage_dims16.rds (list(H_norm, dominant))
#       + H_norm_dims16.csv + projection_cell_audit.csv
#       (projection_sha256.txt 由 Python 生成)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

PROJ     <- "D:/data/microglia_new"
NMF_ROOT <- "E:/micro_output/v3_dims16/05_nmf"
OBJ_D15  <- "E:/micro_output/v3_dims16/objects/microglia_atlas_harmony_dims16_v3.rds"
HVG_FILE <- file.path(NMF_ROOT, "objects", "hvg5000.csv")
C2_RDS   <- file.path(NMF_ROOT, "02_consensus", "consensus_programs.rds")
P3_DIR   <- file.path(NMF_ROOT, "03_projection")
dir.create(P3_DIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(OBJ_D15), file.exists(HVG_FILE), file.exists(C2_RDS))

log_msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...), "\n")
log_msg("=== 04 projection start ===")

obj <- readRDS(OBJ_D15)
stopifnot(ncol(obj) == 424222L)
hvg <- readLines(HVG_FILE)
hvg <- hvg[nchar(hvg) > 0][seq_len(5000L)]
cnt <- LayerData(obj, assay = "RNA", layer = "counts")[hvg, ]
X <- as(cnt, "dgCMatrix")                  # 5000 x 384,805 (genes x cells)
cell_ids <- colnames(cnt)
rm(cnt); gc()
log_msg("X: ", nrow(X), " genes x ", ncol(X), " cells")

cons <- readRDS(C2_RDS)
W <- cons$consensus_W                      # 5000 x MP
stopifnot(identical(rownames(W), hvg))
MOD <- colnames(W)
log_msg("consensus W: ", nrow(W), " x ", length(MOD), " MP")

# ---- RcppML project: 固定 W 投影 (每细胞 NNLS) ----
log_msg("RcppML::project ...")
h <- RcppML::project(A = X, w = W, nonneg = TRUE)   # MP x cells
H <- as.matrix(h)
rownames(H) <- MOD
colnames(H) <- cell_ids
log_msg("H: ", nrow(H), " x ", ncol(H))

# ---- H_norm: 每细胞列归一化 (usage 分数 0-1) ----
col_sum <- colSums(H)
col_sum[col_sum < 1e-12] <- 1
H_norm <- sweep(H, 2, col_sum, "/")
stopifnot(!anyNA(H_norm), !any(is.infinite(H_norm)), all(H_norm >= 0),
          all(is.finite(colSums(H_norm))))

dominant <- MOD[apply(H_norm, 2, which.max)]
usage <- list(H_norm = H_norm, dominant = dominant)
saveRDS(usage, file.path(P3_DIR, "primary_MP_v2_usage_dims16.rds"))

# ---- 审计 ----
audit <- data.frame(
  item = c("H_norm_rows", "H_norm_cols", "rownames_are_MP", "colnames_are_CellID",
           "cellid_match_object", "cellid_unique", "no_na", "no_inf", "no_negative",
           "cell_usage_sum", "min_col_sum", "max_col_sum",
           "n_cells_dominant_MP", "all_MP_have_usage"),
  value = c(
    nrow(H_norm), ncol(H_norm),
    all(rownames(H_norm) == MOD),
    all(colnames(H_norm) == cell_ids),
    identical(colnames(H_norm), colnames(obj)),
    !anyDuplicated(colnames(H_norm)),
    !anyNA(H_norm), !any(is.infinite(H_norm)), all(H_norm >= 0),
    "列归一化 (每细胞 sum=1)",
    round(min(colSums(H_norm)), 6), round(max(colSums(H_norm)), 6),
    paste(MOD, table(dominant), sep = "=", collapse = "; "),
    all(rowSums(H_norm) > 0)
  ), stringsAsFactors = FALSE)
write.csv(audit, file.path(P3_DIR, "projection_cell_audit.csv"), row.names = FALSE)
log_msg("projection_cell_audit.csv 已写; dominant: ",
        paste(names(table(dominant)), table(dominant), sep = "=", collapse = ", "))
log_msg("== 04 done ==")
