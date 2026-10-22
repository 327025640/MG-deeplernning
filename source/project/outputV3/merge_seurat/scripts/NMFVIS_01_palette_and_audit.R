# =============================================================================
# NMFVIS_01_palette_and_audit.R — 冻结 MP01-MP12 调色板 + 输入 SHA256
# 权威参考: 01_palette_and_audit.R（历史 11 色 + 新增第 12 色 #17BECF）
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

ROOT <- "E:/micro_output/NMF"
NV   <- "E:/micro_output/v3_dims16/05_nmf"
ATLAS <- "E:/micro_output/v3_dims16/objects/microglia_atlas_harmony_dims16_v3.rds"
dir.create(file.path(ROOT, "source_data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ROOT, "tables"), recursive = TRUE, showWarnings = FALSE)

# ---- 冻结 MP palette（前 11 色沿用历史 Figure3 冻结版, MP12 新增）----
pal <- data.table(
  MP = sprintf("MP%02d", 1:12),
  Color_Hex = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948",
                "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "#1F3A93", "#17BECF"),
  Status = "frozen_2026-08-30")
fwrite(pal, file.path(ROOT, "source_data", "NMFVIS_MP_palette.csv"))

# ---- SHA256 / MD5 审计 ----
files <- c(
  usage       = file.path(NV, "03_projection", "primary_MP_v2_usage_dims16.rds"),
  W           = file.path(NV, "02_consensus", "consensus_programs.rds"),
  assignment  = file.path(NV, "02_consensus", "consensus_assignment.csv"),
  top50       = file.path(NV, "02_consensus", "consensus_programs_top50.csv"),
  k_scan      = file.path(NV, "01_bootstrap", "k_scan_all_draws.csv"),
  sampling    = file.path(NV, "01_bootstrap", "sampling_draw_audit.csv"),
  hvg5000     = file.path(NV, "objects", "hvg5000.csv"),
  atlas       = ATLAS,
  registry    = "D:/data/microglia_new/authority/03_reference_configs/microglia_function_axis_registry_v1.csv",
  humica      = "D:/data/microglia_new/input/gene_sets/humica_dam_aging_gene_sets_v1.csv",
  disease_pal = "D:/data/microglia_new/authority/03_reference_configs/disease_palette_contract.csv",
  dataset_pal = "D:/data/microglia_new/authority/03_reference_configs/dataset_palette_contract.csv")
stopifnot(all(file.exists(files)))

sha_fun <- function(f) {
  out <- suppressWarnings(system2("sha256sum", shQuote(f), stdout = TRUE))
  if (is.character(out) && length(out) == 1L && grepl("^[0-9a-f]{64}", out[1])) {
    return(sub(" .*", "", out[1]))
  }
  unname(tools::md5sum(f))  # sha256sum 不可用时的回退（列名注明算法）
}
sha <- data.table(role = names(files), path = unname(files))
sha[, checksum := vapply(path, sha_fun, character(1))]
sha[, algo := fifelse(grepl("^[0-9a-f]{64}$", checksum), "sha256", "md5")]
sha[, size_bytes := vapply(path, file.size, numeric(1))]
fwrite(sha, file.path(ROOT, "tables", "NMFVIS_input_checksums.csv"))
cat("palette frozen (12 MP) + checksum audit |", nrow(sha), "inputs\n")
