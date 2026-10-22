# =============================================================================
# 10_nmf_contract_hvg.R — dims16 正式 NMF: 契约审计 + HVG5000 重算
# 冻结流程移植自 outputV2/03_NMF_dims15（仅改输入/输出路径与断言，方法不变）
# 契约参数: N_BOOTSTRAP=100  MAX_CELLS_PER_DISEASE=1000  K_RANGE=4:15
#           N_INIT=20  BASE_SEED=20260801  RcppML::nmf(x,k,tol=1e-4,maxit=100)
#           NMF 输入 = RNA counts 的 HVG5000(蛋白编码) 子集; 不用 Harmony/UMAP/cluster
# =============================================================================

suppressPackageStartupMessages({ library(Seurat); library(data.table) })

ATLAS  <- "E:/micro_output/v3_dims16/objects/microglia_atlas_harmony_dims16_v3.rds"
NMF_ROOT <- "E:/micro_output/v3_dims16/05_nmf"
PC_REF <- "D:/data/microglia_new/authority/03_reference_configs/gencode_v50_protein_coding_genes.csv"
dir.create(file.path(NMF_ROOT, "00_contract"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(NMF_ROOT, "objects"), recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(ATLAS), file.exists(PC_REF))

log_msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...), "\n")
log_msg("=== 10 contract audit + HVG5000 start ===")

pc <- fread(PC_REF)
pc_genes <- unique(pc$gene_name)
log_msg("protein coding reference genes: ", length(pc_genes))

cat(sprintf("[%s] loading atlas...\n", format(Sys.time(), "%H:%M:%S")))
obj <- readRDS(ATLAS)
md <- obj[[]]
md$CellID <- colnames(obj)
all_genes <- rownames(obj)
pc_use <- intersect(all_genes, pc_genes)
log_msg("atlas: ", ncol(obj), " cells x ", length(all_genes),
        " genes; protein coding in object: ", length(pc_use))

req <- c("SampleID", "DonorID", "DatasetID", "SampleID_Global", "Disease_Abbrev_Analysis")
miss <- setdiff(req, names(md))
stopifnot(length(miss) == 0L)
stopifnot(!anyDuplicated(md$CellID))
n_dis <- uniqueN(md$Disease_Abbrev_Analysis)
log_msg("disease levels: ", n_dis)

cat(sprintf("[%s] vst HVG 5000 on protein-coding subset...\n", format(Sys.time(), "%H:%M:%S")))
obj_pc <- subset(obj, features = pc_use)
obj_pc <- FindVariableFeatures(obj_pc, selection.method = "vst", nfeatures = 5000, verbose = FALSE)
hvg <- VariableFeatures(obj_pc)
stopifnot(length(hvg) == 5000L)
writeLines(hvg, file.path(NMF_ROOT, "objects", "hvg5000.csv"))
log_msg("HVG5000 written")

audit <- data.frame(
  item = c("n_cells", "n_genes_total", "n_genes_protein_coding", "n_hvg_selected",
           "n_datasets", "n_samples", "n_disease_levels", "cellid_unique",
           "required_meta_present", "counts_layer_present", "data_layer_present"),
  value = c(ncol(obj), length(all_genes), length(pc_use), length(hvg),
            uniqueN(md$DatasetID), uniqueN(md$SampleID), n_dis,
            !anyDuplicated(md$CellID), length(miss) == 0L,
            "counts" %in% Layers(obj[["RNA"]]), "data" %in% Layers(obj[["RNA"]])),
  stringsAsFactors = FALSE)
fwrite(audit, file.path(NMF_ROOT, "00_contract", "input_object_audit.csv"), row.names = FALSE)

contract <- c(
  "# NMF_CONTRACT_dims16_v1 (frozen, transplanted from NMF_CONTRACT_dims15_v2)",
  "# 仅改输入/输出路径; 算法/种子/参数与 dims15 冻结版逐项一致",
  "N_BOOTSTRAP=100",
  "MAX_CELLS_PER_DISEASE=1000 (Disease_Abbrev_Analysis 分层, ctrl 作为一个疾病组)",
  "K_RANGE=4:15",
  "N_INIT=20",
  "BASE_SEED=20260801",
  "NMF_TOL=1e-4",
  "NMF_MAXIT=100",
  "抽样 seed = 20260801 + i*100",
  "初始化 seed = 20260801 + i*100 + K*1000 + j",
  "求解器 = RcppML::nmf(x, k, tol=1e-4, maxit=100, verbose=FALSE)",
  "NMF 输入 = RNA counts 的 HVG5000(蛋白编码 vst 重算) 子集",
  "不使用 Harmony/UMAP/cluster/donor 信息",
  "数据驱动选 K (诊断图供裁决), 禁止强行裁成任何预设程序数")
writeLines(contract, file.path(NMF_ROOT, "00_contract", "NMF_CONTRACT_dims16_v1.md"))
log_msg("== 10 done ==")
