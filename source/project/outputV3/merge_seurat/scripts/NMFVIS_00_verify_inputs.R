# =============================================================================
# NMFVIS_00_verify_inputs.R — 12 MP 可视化前置：输入总核验 + 细胞映射表导出
# 权威参考: output/manuscript_figures/Figure3_NMF/scripts/00_verify_inputs.R
# 适配: 输入改为 E:/micro_output/v3_dims16/05_nmf 的 12 MP 结果；MP01-MP12；
#       导出 cell_metadata_map.csv 供下游使用（替代旧 cell_sample_map.csv）。
# 本脚本只读 NMF 结果与 atlas meta（不修改 NMF 结果）。
# =============================================================================

suppressPackageStartupMessages({ library(Seurat); library(data.table) })

NV     <- "E:/micro_output/v3_dims16/05_nmf"
ATLAS  <- "E:/micro_output/v3_dims16/objects/microglia_atlas_harmony_dims16_v3.rds"
PC_REF <- "D:/data/microglia_new/authority/03_reference_configs/gencode_v50_protein_coding_genes.csv"
ROOT   <- "E:/micro_output/NMF"
dir.create(file.path(ROOT, "tables"), recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...), "\n")
log_msg("=== NMFVIS_00 verify inputs start ===")
rep <- character(0)
ok <- function(name, val) { rep <<- c(rep, sprintf("%s: %s", name, val)); cat(" ", name, "=", val, "\n") }

# ---- consensus ----
cns <- readRDS(file.path(NV, "02_consensus", "consensus_programs.rds"))
W <- cns$consensus_W
ok("W_dim", paste(dim(W), collapse = "x"))
ok("W_cols", paste(colnames(W), collapse = ","))
stopifnot(identical(colnames(W), sprintf("MP%02d", 1:12)))
ok("W_rownames_HVG", length(rownames(W)))

# ---- usage ----
u <- readRDS(file.path(NV, "03_projection", "primary_MP_v2_usage_dims16.rds"))
H <- u$H_norm
ok("H_norm_dim", paste(dim(H), collapse = "x"))
ok("H_rows", paste(rownames(H), collapse = ","))
stopifnot(identical(rownames(H), colnames(W)))
ok("H_colSums_minmax", paste(round(range(colSums(H)), 6), collapse = "-"))
stopifnot(!anyNA(H), all(is.finite(H)), all(H >= 0))

# ---- assignment ----
a <- fread(file.path(NV, "02_consensus", "consensus_assignment.csv"))
ok("assignment_rows", nrow(a))
ok("assignment MPs", paste(sort(unique(a$cluster_final)), collapse = ","))
ok("assignment draws", uniqueN(a$draw))

# ---- top50 / k_scan / draws ----
t50 <- fread(file.path(NV, "02_consensus", "consensus_programs_top50.csv"))
ok("top50_rows", nrow(t50))
ks <- fread(file.path(NV, "01_bootstrap", "k_scan_all_draws.csv"))
ok("k_scan_rows", nrow(ks))
ok("draw_rds_count", length(list.files(file.path(NV, "01_bootstrap"), "^draw_[0-9]+\\.rds$")))

# ---- atlas meta: 导出 cell map（下游一律用该表对齐）----
cat(sprintf("[%s] loading atlas meta...\n", format(Sys.time(), "%H:%M:%S")))
obj <- readRDS(ATLAS)
md <- as.data.table(obj[[]])
md[, CellID := colnames(obj)]
need <- c("CellID", "DatasetID", "SampleID", "DonorID", "SampleID_Global",
          "Disease_Abbrev_Analysis", "seurat_clusters")
stopifnot(all(need %in% names(md)))
stopifnot(setequal(md$CellID, colnames(H)))
fwrite(md[, ..need], file.path(ROOT, "tables", "cell_metadata_map.csv"))
ok("cell_map_rows", nrow(md))
ok("n_samples", uniqueN(md$SampleID_Global))
ok("n_disease", uniqueN(md$Disease_Abbrev_Analysis))

# ---- 蛋白编码计数核验：不覆盖上游审计（不可变性），修正副本写本目录 ----
pc <- fread(PC_REF)
pc_col <- intersect(c("gene_name", "Gene_Symbol", "gene"), names(pc))[1]
pc_genes <- unique(pc[[pc_col]])
pc_use <- length(intersect(rownames(obj), pc_genes))
audit_path <- file.path(NV, "00_contract", "input_object_audit.csv")
audit <- fread(audit_path)
audit[item == "n_genes_protein_coding",
      value := paste0(as.character(pc_use), " (corrected copy; upstream file untouched)")]
fwrite(audit, file.path(ROOT, "tables", "input_object_audit_corrected_copy.csv"),
       row.names = FALSE)
ok("protein_coding_in_object", pc_use)
rm(obj, md); invisible(gc())

writeLines(c("=== NMFVIS_00 verify report ===", rep),
           file.path(ROOT, "tables", "verify_inputs_report.txt"))
log_msg("== NMFVIS_00 done: all core inputs verified ==")
