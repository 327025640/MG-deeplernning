options(stringsAsFactors = FALSE)
set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(SeuratObject)
  library(Matrix)
})

ALL_ROOT <- "E:/micro_output/fig3/MP_FUN_PTG02/pan_disease_allcells_v3"
ALL_DATA <- file.path(ALL_ROOT, "data")
ALL_TABLES <- file.path(ALL_ROOT, "tables")
ALL_FIGURES <- file.path(ALL_ROOT, "figures")
ALL_SOURCE <- file.path(ALL_ROOT, "source_data")
ALL_LOGS <- file.path(ALL_ROOT, "logs")
ALL_SENS <- file.path(ALL_ROOT, "sensitivity")
ALL_SENS_DIRS <- file.path(ALL_SENS, c("data", "tables", "figures", "source_data", "logs"))

ALL_IN <- c(
  merged_obj = "E:/micro_output/v3_dims16/objects/merged_obj.rds",
  cluster_order = "E:/micro_output/fig3/data/05_cluster_row_order.csv",
  cluster_map = "E:/micro_output/fig3/fig4/tables/cluster_name_mapping.csv",
  disease_group = "E:/micro_output/fig3/data/01_Disease17_group_mapping.csv",
  disease_palette = "D:/data/microglia_new/authority/03_reference_configs/disease_palette_contract.csv",
  mp_usage = "E:/micro_output/v3_dims16/05_nmf/03_projection/primary_MP_v2_usage_dims16.rds",
  mp_stability = "E:/micro_output/NMF/tables/MP_gene_stability.csv",
  fun_registry = "E:/micro_output/fig3/shared_function_framework/tables/function_signature_registry.csv",
  fun_palette = "E:/micro_output/fig3/shared_function_framework/tables/function_theme_palette.csv",
  old_fun_scores = "E:/micro_output/fig3/shared_function_framework/data/cell_function_scores.rds",
  old_fun_meta = "E:/micro_output/fig3/shared_function_framework/data/meta_cell_info.csv",
  old_fun_cache_meta = "E:/micro_output/fig3/shared_function_framework/tables/cell_function_scores_cache_meta.csv",
  old_fun_mapping = "E:/micro_output/fig3/shared_function_framework/tables/MP_function_mapping.csv",
  cds = "E:/micro_output/fig3/data/Fig2i_monocle2_cds_default.rds",
  traj_meta = "E:/micro_output/fig3/data/Fig2i_sampled_cells_metadata.csv",
  ptg_membership = "E:/micro_output/fig3/trajectory_S4/data/default/s4h_ptg_membership.csv",
  ptg02_membership_v2 = "E:/micro_output/fig3/MP_FUN_PTG02/pan_disease_v2/P2_branch_repair/tables/PTG02_subprogram_membership_v2.csv",
  old_p2_manifest = "E:/micro_output/fig3/MP_FUN_PTG02/pan_disease_v2/P2_branch_repair/tables/version_change_manifest.csv"
)

PAN14 <- c("C00", "C09", "C01", "C03", "C04", "C07", "C02", "C05", "C06", "C13", "C08", "C19", "C17", "C12")
EXCLUDED_CLUSTERS <- c("C10", "C15", "C20", "C24", "C27", "C28")

all_make_dirs <- function() {
  dirs <- c(ALL_ROOT, ALL_DATA, ALL_TABLES, ALL_FIGURES, ALL_SOURCE, ALL_LOGS, ALL_SENS, ALL_SENS_DIRS)
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

all_check_inputs <- function() {
  required <- c("merged_obj", "cluster_order", "cluster_map", "disease_group", "disease_palette",
                "mp_usage", "fun_registry", "fun_palette", "old_fun_scores", "old_fun_meta",
                "cds", "traj_meta", "ptg_membership", "ptg02_membership_v2")
  missing <- required[!file.exists(ALL_IN[required])]
  if (length(missing)) stop("缺少正式输入：", paste(missing, collapse = ", "))
  invisible(TRUE)
}

all_read_group_map <- function() {
  x <- fread(ALL_IN[["disease_group"]], encoding = "UTF-8")
  if ("Disease_Group3" %in% names(x)) setnames(x, "Disease_Group3", "Group3")
  if (!all(c("Disease", "Group3") %in% names(x))) stop("Disease group mapping 缺少 Disease/Group3")
  x[, Disease := as.character(Disease)]
  x[, Group3 := as.character(Group3)]
  x
}

all_read_cluster_labels <- function() {
  ro <- fread(ALL_IN[["cluster_order"]], encoding = "UTF-8")
  if (!all(c("Cluster", "row_order") %in% names(ro))) stop("cluster row order 缺少 Cluster/row_order")
  cm <- fread(ALL_IN[["cluster_map"]], encoding = "UTF-8")
  if (!all(c("Cluster", "FinalName") %in% names(cm))) stop("cluster name mapping 缺少 Cluster/FinalName")
  x <- merge(ro[Cluster %chin% PAN14, .(Cluster, row_order)], cm, by = "Cluster", all.x = TRUE, sort = FALSE)
  x <- x[order(row_order)]
  if (!identical(x$Cluster, PAN14)) stop("正式 pan14 顺序与冻结 row_order 不一致")
  x[, cluster_label := paste(Cluster, FinalName, sep = " · ")]
  x
}

all_read_disease_palette <- function(diseases) {
  p <- fread(ALL_IN[["disease_palette"]], encoding = "UTF-8")
  if (!all(c("Disease_Abbrev", "Disease_Color_Hex") %in% names(p))) stop("疾病 palette 缺少必要列")
  p <- p[!duplicated(Disease_Abbrev)]
  x <- data.table(Disease = as.character(diseases))
  x <- merge(x, p[, .(Disease = as.character(Disease_Abbrev), Display_Label = if ("Display_Label" %in% names(p)) Display_Label else Disease_Abbrev,
                      Disease_Color_Hex)], by = "Disease", all.x = TRUE, sort = FALSE)
  x[is.na(Disease_Color_Hex), Disease_Color_Hex := "#7F7F7F"]
  x[is.na(Display_Label), Display_Label := Disease]
  x
}

all_load_registry <- function() {
  reg <- fread(ALL_IN[["fun_registry"]], encoding = "UTF-8")
  if (!all(c("FunctionID", "FunctionLabel", "Genes") %in% names(reg))) stop("FUN registry 缺少 FunctionID/FunctionLabel/Genes")
  reg[, FunctionID := as.character(FunctionID)]
  reg[, GeneSetMD5 := if ("GeneSetMD5" %in% names(reg)) as.character(GeneSetMD5) else NA_character_]
  reg
}

all_read_fun_palette <- function() {
  p <- fread(ALL_IN[["fun_palette"]], encoding = "UTF-8")
  if ("FunctionID" %in% names(p)) setnames(p, "FunctionID", "FUN")
  if (!all(c("FUN", "ColorHex") %in% names(p))) stop("FUN palette 缺少 FUN/ColorHex")
  p
}

all_load_mp <- function(cell_ids = NULL) {
  mp_obj <- readRDS(ALL_IN[["mp_usage"]])
  if (!"H_norm" %in% names(mp_obj)) stop("MP usage RDS 缺少 H_norm")
  h <- mp_obj$H_norm
  if (is.null(colnames(h)) || is.null(rownames(h))) stop("H_norm 缺少 CellID 或 MP 名称")
  if (!is.null(cell_ids)) {
    idx <- match(cell_ids, colnames(h))
    if (anyNA(idx)) stop("MP usage 未覆盖全部正式 pan14 CellID，缺失数：", sum(is.na(idx)))
    h <- h[, idx, drop = FALSE]
  }
  dt <- as.data.table(t(h), keep.rownames = "CellID")
  setnames(dt, setdiff(names(dt), "CellID"), rownames(h))
  dt
}

all_load_old_fun <- function() {
  old <- as.data.table(readRDS(ALL_IN[["old_fun_scores"]]))
  if (!"CellID" %in% names(old)) stop("旧 FUN 缓存缺少 CellID")
  old
}

all_load_ptg_sets <- function() {
  p <- fread(ALL_IN[["ptg02_membership_v2"]], encoding = "UTF-8")
  if (!all(c("gene_id", "subprogram") %in% names(p))) stop("P2 PTG02 membership 缺少 gene_id/subprogram")
  p[, gene_id := as.character(gene_id)]
  p[, subprogram := as.character(subprogram)]
  split(p$gene_id, p$subprogram)
}

all_save_session <- function(stage) capture.output(sessionInfo(), file = file.path(ALL_LOGS, paste0(stage, "_sessionInfo.txt")))
all_write_log <- function(stage, lines) writeLines(as.character(lines), file.path(ALL_LOGS, paste0(stage, ".log")), useBytes = TRUE)
all_write_readme <- function(lines) writeLines(lines, file.path(ALL_ROOT, "解读.md"), useBytes = TRUE)

all_summary_long <- function(x, features, by_vars = c("Cluster"), unit_name = "all_cells") {
  y <- melt(x, id.vars = intersect(c("CellID", "SampleID_Global", "Disease", "DatasetID", "DonorID", "Group3", by_vars), names(x)),
            measure.vars = features, variable.name = "Feature", value.name = "value")
  y[, .(mean = mean(value, na.rm = TRUE), median = median(value, na.rm = TRUE),
        q25 = as.numeric(quantile(value, .25, na.rm = TRUE)), q75 = as.numeric(quantile(value, .75, na.rm = TRUE)),
        n_cells = sum(is.finite(value)), n_samples = uniqueN(SampleID_Global), n_donors = uniqueN(DonorID)),
    by = c(by_vars, "Feature")][, unit := unit_name][]
}

all_dynamic_scores <- function(expr, gene_sets) {
  out <- data.table(CellID = colnames(expr))
  for (nm in names(gene_sets)) {
    genes <- intersect(gene_sets[[nm]], rownames(expr))
    out[, (nm) := if (length(genes)) Matrix::colMeans(expr[genes, , drop = FALSE]) else NA_real_]
  }
  out
}
