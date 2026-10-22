# =============================================================================
# 11_nmf_sampling.R — dims16 正式 NMF: 100 draws 分层抽样导出
# 移植自 02a_sampling.R（冻结版）: 仅改输入/输出路径与细胞数断言
# =============================================================================

suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table) })

ATLAS   <- "E:/micro_output/v3_dims16/objects/microglia_atlas_harmony_dims16_v3.rds"
NMF_ROOT <- "E:/micro_output/v3_dims16/05_nmf"
HVG_FILE <- file.path(NMF_ROOT, "objects", "hvg5000.csv")
SAMP_DIR <- file.path(NMF_ROOT, "01_bootstrap", "sampling")
B1_DIR   <- file.path(NMF_ROOT, "01_bootstrap")
dir.create(SAMP_DIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(ATLAS), file.exists(HVG_FILE))

N_DRAWS <- 100L
N_CELLS_PER_DIS <- 1000L
BASE_SEED <- 20260801L

log_msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...), "\n")
log_msg("=== 11 sampling start (100 draws x 1000 cells/disease) ===")

obj <- readRDS(ATLAS)
stopifnot(ncol(obj) == 424222L)
hvg <- readLines(HVG_FILE)
hvg <- hvg[nchar(hvg) > 0]
stopifnot(all(hvg %in% rownames(obj)), length(hvg) >= 5000L)
hvg <- hvg[seq_len(5000L)]
log_msg("HVG: ", length(hvg))

cnt <- LayerData(obj, assay = "RNA", layer = "counts")[hvg, ]
md <- obj[[]]
md$CellID <- colnames(obj)

dis_levels <- unique(md$Disease_Abbrev_Analysis)
dis_levels <- dis_levels[!is.na(dis_levels)]
log_msg("疾病层 (含 Ctrl): ", length(dis_levels))

audit_rows <- list()
for (i in seq_len(N_DRAWS)) {
  ddir <- file.path(SAMP_DIR, sprintf("draw_%03d", i))
  done_txt <- file.path(ddir, "DONE.txt")
  if (file.exists(done_txt)) { log_msg("draw", i, " SKIP (DONE)"); next }
  dir.create(ddir, recursive = TRUE, showWarnings = FALSE)
  set.seed(BASE_SEED + i * 100L)
  samp_cells <- unlist(lapply(dis_levels, function(d) {
    cells_d <- md$CellID[md$Disease_Abbrev_Analysis == d]
    if (length(cells_d) <= N_CELLS_PER_DIS) cells_d else sample(cells_d, N_CELLS_PER_DIS)
  }))
  sub_md <- md[md$CellID %in% samp_cells, ]
  sub_cnt <- cnt[, sub_md$CellID]
  Matrix::writeMM(t(as(sub_cnt, "CsparseMatrix")),
                  file.path(ddir, "expression_cells_by_genes.mtx"))
  write.csv(data.frame(gene = hvg), file.path(ddir, "genes.csv"), row.names = FALSE)
  write.csv(data.frame(CellID = sub_md$CellID), file.path(ddir, "cells.csv"),
            row.names = FALSE)
  sl <- data.frame(draw = i, Disease = dis_levels,
                   n_sampled = vapply(dis_levels, function(d) {
                     sum(sub_md$Disease_Abbrev_Analysis == d)
                   }, integer(1)))
  write.csv(sl, file.path(ddir, "sampling_log.csv"), row.names = FALSE)
  audit_rows[[i]] <- sl
  writeLines("DONE", done_txt)
  if (i %% 10 == 0) log_msg("draw", i, "/", N_DRAWS)
  rm(sub_cnt, sub_md); gc()
}
audit <- rbindlist(audit_rows)
fwrite(audit, file.path(B1_DIR, "sampling_draw_audit.csv"))
log_msg("sampling_draw_audit.csv: ", nrow(audit), " rows")
log_msg("== 11 done ==")
