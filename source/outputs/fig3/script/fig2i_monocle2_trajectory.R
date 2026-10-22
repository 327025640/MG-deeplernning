# Figure 2i: Monocle2 single-cell trajectory of retained microglia subtypes
# Reference: Nie et al., Cell Discovery 2026, doi:10.1038/s41421-026-00888-3
# Paper-aligned core: Monocle2; counts; genes detected in >=10% cells;
# differentialGeneTest P < 0.01; reduceDimension to 2D; orderCells.
# Project adaptations: remove user-adjudicated contaminants; primary scope uses
# Pan-disease retained microglia; deterministic cluster-balanced, sample-balanced sampling.
#
# 2026-09-03 repair notes:
# 1) Root cause of the original failure: the previous run saved cds_default right
#    after reduceDimension WITHOUT a successful orderCells pass (auxOrderingData
#    only contains objective_vals + pr_graph_cell_proj_closest_vertex, no
#    root_cell/pr_graph_cell_proj_tree/branch_points), so pData had no State
#    column and the C00 sensitivity guard stopped.
# 2) The nei() compat patch was wrong for igraph >= 2.1: `igraph::.nei(...)` with
#    a `neimode=` argument breaks [.igraph.vs NSE (namespace-qualified calls are
#    not detected; the argument is named `mode`). v2 rewrites to bare `.nei(` +
#    `mode=`, and resolves `.nei` through a per-function closure environment so
#    no locked namespace is modified.
# 3) A checkpoint of the reduceDimension product (Fig2i_monocle2_cds_default.rds)
#    is reused when it validates, so PCA/Harmony/UMAP/FindNeighbors and the whole
#    Monocle2 pre-processing are NOT recomputed.

options(stringsAsFactors = FALSE)
set.seed(42)

###############################################################################
# 1. Frozen parameters
###############################################################################

INPUT_RDS <- "E:/micro_output/v3_dims16/objects/merged_obj.rds"
TIER_FILE <- "E:/micro_output/fig3/tables/cluster_tier_map.csv"
DISEASE_PALETTE_FILE <- "D:/data/microglia_new/authority/03_reference_configs/disease_palette_contract.csv"
OUT_ROOT <- "E:/micro_output/fig3"

CLUSTER_COL <- "RNA_snn_res.0.6"
TRAJECTORY_SCOPE <- "pan_disease"
CONTAMINANT_CLUSTERS <- sprintf("C%02d", c(10, 15, 20, 24, 27, 28))
PAN_TIER_LABEL <- "Pan-disease type"
HOMEOSTATIC_CLUSTER <- "C00"   # user-adjudicated homeostatic cluster for root sensitivity

MAX_CELLS_PER_CLUSTER <- 1200L
MAX_CELLS_TOTAL <- 30000L
N_CORES <- 4L
SEED <- 42L

MIN_DETECTION_FRACTION <- 0.10
ORDERING_P_CUTOFF <- 0.01
EXPECTED_MONOCLE_VERSION <- "2.28.0"

FIG_DIR <- file.path(OUT_ROOT, "figures")
TABLE_DIR <- file.path(OUT_ROOT, "tables")
SOURCE_DIR <- file.path(OUT_ROOT, "source_data")
DATA_DIR <- file.path(OUT_ROOT, "data")
LOG_DIR <- file.path(OUT_ROOT, "logs")
CHECKPOINT_DEFAULT <- file.path(DATA_DIR, "Fig2i_monocle2_cds_default.rds")

###############################################################################
# 2. Package and input checks (no automatic installation)
###############################################################################

required_pkgs <- c(
  "Seurat", "SeuratObject", "Matrix", "data.table", "ggplot2",
  "monocle", "Biobase", "scales", "viridis", "dplyr", "rlang"
)
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs)) {
  stop(
    "缺少以下R包（脚本不会自动安装）：", paste(missing_pkgs, collapse = ", "),
    "\n参考论文使用 monocle 2.28.0。请安装后重新运行。"
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(monocle)
  library(Biobase)
  library(scales)
  library(viridis)
  library(dplyr)
  library(rlang)
})

# ---- Legacy dplyr patch: Monocle2 still calls group_by_()/select_() ----------
if (utils::packageVersion("dplyr") >= "1.1.0") {
  patch_legacy_dplyr <- function(ns) {
    if (exists("group_by_", envir = ns, inherits = FALSE)) {
      unlockBinding("group_by_", ns)
      assign("group_by_", function(.data, ..., .dots = NULL, add = FALSE) {
        vars <- c(list(...), as.list(.dots))
        if (!length(vars)) return(dplyr::group_by(.data, .add = add))
        dplyr::group_by(.data, !!!vars, .add = add)
      }, envir = ns)
      lockBinding("group_by_", ns)
    }
    if (exists("select_", envir = ns, inherits = FALSE)) {
      unlockBinding("select_", ns)
      assign("select_", function(.data, ..., .dots = NULL) {
        vars <- c(list(...), as.list(.dots))
        if (!length(vars)) return(dplyr::select(.data))
        dplyr::select(.data, !!!vars)
      }, envir = ns)
      lockBinding("select_", ns)
    }
  }
  patch_legacy_dplyr(asNamespace("dplyr"))
  # Monocle resolves these imports through its imports environment.
  patch_legacy_dplyr(parent.env(asNamespace("monocle")))
}

# ---- Legacy igraph patch 1: dfs(neimode=) -> dfs(mode=) ----------------------
# Monocle2 (e.g. extract_ddrtree_ordering via graph.dfs) still passes neimode.
if (utils::packageVersion("igraph") >= "1.3.0") {
  igraph_ns <- asNamespace("igraph")
  dfs_current <- get("dfs", envir = igraph_ns, inherits = FALSE)
  unlockBinding("dfs", igraph_ns)
  assign("dfs", function(graph, root, mode = c("out", "in", "all", "total"),
                          ..., neimode = NULL) {
    if (!is.null(neimode)) mode <- neimode
    dfs_current(graph = graph, root = root, mode = mode, ...)
  }, envir = igraph_ns)
  lockBinding("dfs", igraph_ns)
}

# ---- Legacy igraph patch 2 (v2): nei() -> .nei() -----------------------------
# Monocle2 tree-walking helpers use the old indexing helper
#   V(graph)[suppressWarnings(nei(v, neimode = "all"))]
# igraph >= 2.1 removed nei(); the supported replacement is `.nei()` used inside
# the SAME indexing context, with the argument renamed to `mode`. Two constraints
# discovered by probing igraph 2.3.2:
#   (a) a namespace-qualified call `igraph::::.nei(...)` is NOT recognized by
#       [.igraph.vs non-standard evaluation -> must rewrite to bare `.nei(`;
#   (b) the argument must be renamed neimode -> mode.
# `.nei` is exported by igraph but not imported by monocle; instead of touching
# the locked namespaces, each patched function gets a child environment that
# carries `.nei`, so lexical lookup resolves it first.
igraph_nei_impl <- get(".nei", envir = asNamespace("igraph"), inherits = FALSE)
monocle_ns <- asNamespace("monocle")
nei_patched_ok <- character(0)
nei_patched_fail <- character(0)
for (fun_name in ls(monocle_ns, all.names = TRUE)) {
  fun_obj <- tryCatch(get(fun_name, envir = monocle_ns), error = function(e) NULL)
  if (!is.function(fun_obj)) next
  fun_text <- paste(deparse(body(fun_obj)), collapse = "\n")
  if (!grepl("(?<![[:alnum:]_.])nei\\(", fun_text, perl = TRUE)) next
  new_text <- gsub("(?<![[:alnum:]_.])nei\\(", ".nei(", fun_text, perl = TRUE)
  new_text <- gsub("\\bneimode\\b", "mode", new_text)
  parsed <- tryCatch(parse(text = new_text)[[1]], error = function(e) {
    message("nei-patch parse failure for ", fun_name, ": ", conditionMessage(e))
    NULL
  })
  if (is.null(parsed)) {
    nei_patched_fail <- c(nei_patched_fail, fun_name)
    next
  }
  patched_fun <- fun_obj
  body(patched_fun) <- parsed
  patched_env <- new.env(parent = asNamespace("monocle"))
  assign(".nei", igraph_nei_impl, envir = patched_env)
  environment(patched_fun) <- patched_env
  unlockBinding(fun_name, monocle_ns)
  assign(fun_name, patched_fun, envir = monocle_ns)
  lockBinding(fun_name, monocle_ns)
  nei_patched_ok <- c(nei_patched_ok, fun_name)
}
message("nei-patch v2 applied to: ", paste(nei_patched_ok, collapse = ", "))
if (length(nei_patched_fail)) {
  stop("nei-patch v2 failed on: ", paste(nei_patched_fail, collapse = ", "),
       " —— 请人工检查这些函数，不能带病继续。")
}

for (d in c(FIG_DIR, TABLE_DIR, SOURCE_DIR, DATA_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

monocle_version <- as.character(utils::packageVersion("monocle"))
if (!identical(monocle_version, EXPECTED_MONOCLE_VERSION)) {
  warning(
    "当前 monocle 版本为 ", monocle_version,
    "；参考论文版本为 ", EXPECTED_MONOCLE_VERSION,
    "。脚本继续运行，但必须在方法审计中披露版本差异。"
  )
}

###############################################################################
# 3. Trajectory scope from frozen tier map + C00 label normalization
###############################################################################

if (!file.exists(TIER_FILE)) stop("必要输入文件不存在：", TIER_FILE)
if (!file.exists(DISEASE_PALETTE_FILE)) stop("必要输入文件不存在：", DISEASE_PALETTE_FILE)

tier_dt <- fread(TIER_FILE)
stopifnot(all(c("cluster", "tier", "n_dis_present", "total_cells") %in% names(tier_dt)))
tier_dt[, cluster := as.character(cluster)]

all_clusters <- sprintf("C%02d", 0:29)
retained_microglia <- setdiff(all_clusters, CONTAMINANT_CLUSTERS)
trajectory_clusters <- intersect(
  retained_microglia,
  tier_dt[tier == PAN_TIER_LABEL, cluster]
)
trajectory_clusters <- trajectory_clusters[order(as.integer(sub("C", "", trajectory_clusters)))]
if (length(trajectory_clusters) < 2L) stop("纳入轨迹的cluster少于2个，无法构建轨迹。")

# C00 label normalization: accept equivalent spellings ("0", "C0", "c00", ...)
# and map them to the canonical "C00" used everywhere in this script.
canonical_c00 <- function(candidate_labels) {
  labels <- as.character(candidate_labels)
  normed <- toupper(gsub("[^A-Za-z0-9]", "", labels))
  hit <- labels[!is.na(normed) & normed == "C00"]
  if (length(hit)) hit[1] else NA_character_
}
c00_label <- canonical_c00(c(all_clusters, tier_dt$cluster))
c00_in_primary <- !is.na(c00_label) && c00_label %in% trajectory_clusters
c00_is_contaminant <- !is.na(c00_label) && c00_label %in% CONTAMINANT_CLUSTERS
if (is.na(c00_label)) {
  stop("未在簇标签或 tier map 中找到 C00 的任何等价标签（0/C0/C00），无法执行 C00-root 敏感性分析。")
}
if (c00_is_contaminant) {
  stop("C00 被标记为用户裁决污染簇，与“C00 作为稳态簇”的敏感性分析前提冲突；请先复核污染裁决。")
}

cluster_decision <- data.table(
  Cluster = all_clusters,
  UserDecision = fifelse(
    all_clusters %in% CONTAMINANT_CLUSTERS, "contaminant_excluded", "microglia_retained"
  ),
  Tier = tier_dt$tier[match(all_clusters, tier_dt$cluster)],
  IncludedInPrimaryTrajectory = all_clusters %in% trajectory_clusters,
  ExclusionReason = fifelse(
    all_clusters %in% CONTAMINANT_CLUSTERS,
    "User-adjudicated contaminant",
    fifelse(
      all_clusters %in% trajectory_clusters,
      "Included",
      "Retained microglia but not Pan-disease; excluded from primary trajectory"
    )
  )
)
fwrite(cluster_decision, file.path(TABLE_DIR, "Fig2i_cluster_inclusion_audit.csv"))

eligible_cells <- NULL  # only materialized on the full (non-checkpoint) path

###############################################################################
# 4. Checkpoint reuse: skip resampling/normalization/ordering-gene selection/
#    reduceDimension when the saved post-reduceDimension cds validates
###############################################################################

cds_default <- NULL
resume_checkpoint <- FALSE
checkpoint_note <- "no checkpoint found"

if (file.exists(CHECKPOINT_DEFAULT)) {
  cds_default <- tryCatch(readRDS(CHECKPOINT_DEFAULT), error = function(e) NULL)
}
if (!is.null(cds_default)) {
  ck_ok <- tryCatch({
    stopifnot(inherits(cds_default, "CellDataSet"))
    stopifnot(ncol(cds_default) > 0, nrow(cds_default) > 0)
    pd_ck <- pData(cds_default)
    stopifnot(all(c("CellID", "Cluster", "Disease", "SampleID_Global", "DatasetID") %in% colnames(pd_ck)))
    stopifnot(setequal(levels(pd_ck$Cluster), trajectory_clusters))
    stopifnot(!is.null(minSpanningTree(cds_default)))
    stopifnot(identical(cds_default@dim_reduce_type, "DDRTree"))
    stopifnot(length(reducedDimS(cds_default)) > 0, nrow(reducedDimS(cds_default)) == 2)
    TRUE
  }, error = function(e) {
    message("checkpoint 校验失败：", conditionMessage(e))
    FALSE
  })
  if (ck_ok) {
    resume_checkpoint <- TRUE
    checkpoint_note <- paste0(
      "reused post-reduceDimension checkpoint (", ncol(cds_default),
      " cells; Cluster levels identical to current trajectory scope; ",
      "resampling/normalization/ordering-gene selection/DDRTree NOT recomputed)"
    )
  } else {
    rm(cds_default); invisible(gc())
    checkpoint_note <- "checkpoint invalid -> full pipeline rerun"
  }
}
message("checkpoint: ", checkpoint_note)

if (!resume_checkpoint) {
  if (!file.exists(INPUT_RDS)) stop("必要输入文件不存在：", INPUT_RDS)

  ###########################################################################
  # 4b. Full path: read Seurat object, validate, sample, build cds, reduce
  ###########################################################################

  obj <- readRDS(INPUT_RDS)
  DefaultAssay(obj) <- "RNA"
  md <- obj[[]]

  required_meta <- c(
    CLUSTER_COL, "seurat_clusters", "Disease", "SampleID_Global",
    "DatasetID", "CellID"
  )
  stopifnot(
    ncol(obj) == 424222,
    nrow(obj) == 18687,
    all(required_meta %in% colnames(md)),
    identical(as.character(md[[CLUSTER_COL]]), as.character(md$seurat_clusters)),
    identical(as.character(md$CellID), colnames(obj)),
    all(c("counts", "data") %in% Layers(obj[["RNA"]]))
  )

  cluster_num <- suppressWarnings(as.integer(as.character(md[[CLUSTER_COL]])))
  if (anyNA(cluster_num)) stop(CLUSTER_COL, " 不能完整转换为0-29整数簇号。")
  md$Cluster <- sprintf("C%02d", cluster_num)
  stopifnot(setequal(unique(md$Cluster), sprintf("C%02d", 0:29)))

  eligible_cells <- rownames(md)[md$Cluster %in% trajectory_clusters]
  if (!length(eligible_cells)) stop("没有符合轨迹纳入条件的细胞。")

  cell_info <- data.table(
    CellID = rownames(md),
    Cluster = md$Cluster,
    Disease = as.character(md$Disease),
    SampleID_Global = as.character(md$SampleID_Global),
    DatasetID = as.character(md$DatasetID)
  )[CellID %in% eligible_cells]

  set.seed(SEED)
  sampled_info <- cell_info[, {
    x <- copy(.SD)
    x[, random_key := stats::runif(.N)]
    x[, within_sample_rank := frank(random_key, ties.method = "first"),
      by = SampleID_Global]
    setorder(x, within_sample_rank, SampleID_Global, random_key)
    x <- x[seq_len(min(.N, MAX_CELLS_PER_CLUSTER))]
    x[, c("random_key", "within_sample_rank") := NULL]
    x
  }, by = Cluster]
  if (nrow(sampled_info) > MAX_CELLS_TOTAL) {
    cap_per_cluster <- max(1L, floor(MAX_CELLS_TOTAL / length(trajectory_clusters)))
    sampled_info <- sampled_info[
      , .SD[sample.int(.N, min(.N, cap_per_cluster))],
      by = Cluster
    ]
  }

  setorder(sampled_info, Cluster, Disease, SampleID_Global, CellID)
  sampled_cells <- sampled_info$CellID
  stopifnot(!anyDuplicated(sampled_cells), all(sampled_cells %in% colnames(obj)))

  sampling_audit <- merge(
    cell_info[, .(n_before = .N, n_samples_before = uniqueN(SampleID_Global)), by = .(Cluster, Disease)],
    sampled_info[, .(n_after = .N, n_samples_after = uniqueN(SampleID_Global)), by = .(Cluster, Disease)],
    by = c("Cluster", "Disease"), all = TRUE
  )
  sampling_audit[is.na(n_before), `:=`(n_before = 0L, n_samples_before = 0L)]
  sampling_audit[is.na(n_after), `:=`(n_after = 0L, n_samples_after = 0L)]
  fwrite(sampling_audit, file.path(TABLE_DIR, "Fig2i_sampling_audit_cluster_disease.csv"))
  fwrite(sampled_info, file.path(DATA_DIR, "Fig2i_sampled_cells_metadata.csv"))

  counts_mat <- LayerData(obj, assay = "RNA", layer = "counts")[, sampled_cells, drop = FALSE]
  counts_mat <- as(counts_mat, "dgCMatrix")
  if (any(counts_mat@x < 0)) stop("RNA counts中检测到负值。")

  sampled_md <- md[sampled_cells, , drop = FALSE]
  sampled_md$Cluster <- sprintf(
    "C%02d", as.integer(as.character(sampled_md[[CLUSTER_COL]]))
  )
  sampled_md$Disease <- as.character(sampled_md$Disease)
  sampled_md$SampleID_Global <- as.character(sampled_md$SampleID_Global)
  sampled_md$DatasetID <- as.character(sampled_md$DatasetID)

  pheno <- data.frame(
    CellID = sampled_cells,
    Cluster = factor(sampled_md$Cluster, levels = trajectory_clusters),
    Disease = sampled_md$Disease,
    SampleID_Global = sampled_md$SampleID_Global,
    DatasetID = sampled_md$DatasetID,
    row.names = sampled_cells,
    check.names = FALSE
  )
  feature <- data.frame(
    gene_short_name = rownames(counts_mat),
    row.names = rownames(counts_mat),
    check.names = FALSE
  )

  pd <- new("AnnotatedDataFrame", data = pheno)
  fd <- new("AnnotatedDataFrame", data = feature)
  cds <- newCellDataSet(
    counts_mat,
    phenoData = pd,
    featureData = fd,
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size()
  )

  rm(obj, md, counts_mat)
  invisible(gc())

  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  cds <- detectGenes(cds, min_expr = 0.1)

  min_detected_cells <- ceiling(MIN_DETECTION_FRACTION * ncol(cds))
  eligible_genes <- rownames(subset(
    fData(cds), num_cells_expressed >= min_detected_cells
  ))
  if (length(eligible_genes) < 100L) {
    stop("满足>=10%细胞表达的基因少于100个：", length(eligible_genes))
  }

  ordering_test <- differentialGeneTest(
    cds[eligible_genes, ],
    fullModelFormulaStr = "~Cluster",
    cores = N_CORES
  )
  ordering_test$gene <- rownames(ordering_test)
  ordering_test$detected_fraction <- fData(cds)[ordering_test$gene, "num_cells_expressed"] / ncol(cds)
  ordering_test$ordering_gene <- is.finite(ordering_test$pval) &
    ordering_test$pval < ORDERING_P_CUTOFF &
    ordering_test$detected_fraction >= MIN_DETECTION_FRACTION
  ordering_genes <- ordering_test$gene[ordering_test$ordering_gene]
  if (length(ordering_genes) < 50L) {
    stop("P<0.01且表达比例>=10%的ordering genes少于50个：", length(ordering_genes))
  }

  fwrite(as.data.table(ordering_test), file.path(TABLE_DIR, "Fig2i_ordering_gene_test.csv"))
  cds <- setOrderingFilter(cds, ordering_genes)

  set.seed(SEED)
  cds <- reduceDimension(
    cds,
    max_components = 2,
    reduction_method = "DDRTree",
    norm_method = "log",
    pseudo_expr = 1,
    relative_expr = TRUE,
    auto_param_selection = TRUE,
    verbose = TRUE
  )
  cds_default <- cds
  rm(cds); invisible(gc())
  ordering_gene_count_used <- length(ordering_genes)
  cells_before_sampling_used <- nrow(cell_info)
  cells_after_sampling_used <- nrow(sampled_info)
  checkpoint_note_full <- "full pipeline rerun (no valid checkpoint)"
} else {
  ordering_gene_count_used <- NA_integer_
  if (!is.null(fData(cds_default)$use_for_ordering)) {
    ordering_gene_count_used <- sum(fData(cds_default)$use_for_ordering)
  }
  cells_before_sampling_used <- NA_integer_
  cells_after_sampling_used <- ncol(cds_default)
  checkpoint_note_full <- checkpoint_note
}

# Validate the frozen Disease palette against whatever the trajectory contains.
pal <- fread(DISEASE_PALETTE_FILE)
stopifnot(all(c("Disease_Abbrev", "Disease_Color_Hex") %in% names(pal)))
pal <- pal[!duplicated(Disease_Abbrev, fromLast = TRUE)]
disease_colors <- setNames(pal$Disease_Color_Hex, pal$Disease_Abbrev)
disease_colors[c("Ctrl", "LGG", "GBM")] <- c("#4DAF4A", "#8C1D18", "#B6362D")
missing_disease_colors <- setdiff(unique(pData(cds_default)$Disease), names(disease_colors))
if (length(missing_disease_colors)) {
  stop("Disease固定配色缺失：", paste(missing_disease_colors, collapse = ", "))
}

###############################################################################
# 5. Monocle2 ordering: primary default direction + State guard
###############################################################################

cds_default <- orderCells(cds_default)

# Guard: orderCells must have written State and Pseudotime. If a legacy/patched
# code path returns without them, rebuild both exactly as the DDRTree branch of
# orderCells does, reusing the projection stored in auxOrderingData.
state_guard_note <- "orderCells wrote State and Pseudotime directly"
if (is.null(pData(cds_default)$State) || is.null(pData(cds_default)$Pseudotime)) {
  aux_ddrtree <- cds_default@auxOrderingData[["DDRTree"]]
  proj_tree <- aux_ddrtree$pr_graph_cell_proj_tree
  proj_root <- aux_ddrtree$root_cell
  if (is.null(proj_tree) || is.null(proj_root)) {
    stop("orderCells 未写入 State/Pseudotime，且 auxOrderingData 缺少投影树，无法兜底重建。")
  }
  old_mst <- minSpanningTree(cds_default)
  minSpanningTree(cds_default) <- proj_tree
  cc_cell <- extract_ddrtree_ordering(cds_default, proj_root)
  minSpanningTree(cds_default) <- old_mst
  state_vec <- cc_cell[row.names(pData(cds_default)), ]$cell_state
  ptime_vec <- cc_cell[row.names(pData(cds_default)), ]$pseudo_time
  if (all(is.na(state_vec)) || all(is.na(ptime_vec))) {
    stop("State/Pseudotime 兜底重建失败（全 NA）。")
  }
  pData(cds_default)$State <- factor(state_vec)
  pData(cds_default)$Pseudotime <- ptime_vec
  state_guard_note <- paste0(
    "orderCells returned without State/Pseudotime; rebuilt from projected-tree ",
    "ordering with root cell ", proj_root, " (orderCells DDRTree branch semantics)"
  )
  message("State/Pseudotime 由投影树兜底重建：root=", proj_root)
}

saveRDS(cds_default, file.path(DATA_DIR, "Fig2i_monocle2_cds_default.rds"), compress = FALSE)

###############################################################################
# 6. C00-root sensitivity analysis
###############################################################################

# The default direction is the primary result. The sensitivity analysis only
# re-roots pseudotime at the State carrying the most C00 cells; Monocle2 keeps
# the State partition unchanged under root_state (states are a property of the
# principal graph, not of the chosen root).
state_c00 <- table(pData(cds_default)$State[pData(cds_default)$Cluster == c00_label])
if (!length(state_c00) || sum(state_c00) == 0) {
  stop("默认轨迹中未找到C00细胞（label=", c00_label, "），无法执行C00起点敏感性分析。")
}
c00_root_state <- as.numeric(names(state_c00)[which.max(state_c00)])
c00_cells_in_primary <- sum(state_c00)

cds_c00_root <- orderCells(cds_default, root_state = c00_root_state)
# orderCells(root_state=) intentionally keeps the State partition and only
# rewrites Pseudotime; enforce/preserve that explicitly for the record.
if (is.null(pData(cds_c00_root)$State)) {
  pData(cds_c00_root)$State <- pData(cds_default)$State
}
saveRDS(cds_default, file.path(DATA_DIR, "Fig2i_monocle2_cds.rds"), compress = FALSE)
saveRDS(cds_c00_root, file.path(DATA_DIR, "Fig2i_monocle2_cds_C00_root.rds"), compress = FALSE)

c00_state_note <- paste0(
  "C00 (canonical label ", c00_label, ") is part of the primary pan-disease ",
  "trajectory with ", c00_cells_in_primary, " sampled cells; its majority State ",
  c00_root_state, " was used as the pseudotime origin for the sensitivity run."
)

###############################################################################
# 7. Source data and trajectory summaries
###############################################################################

reduced <- reducedDimS(cds_default)
pd_final <- as.data.table(pData(cds_default), keep.rownames = "CellID")
trajectory_source <- data.table(
  CellID = colnames(reduced),
  DDRTree_1 = as.numeric(reduced[1, ]),
  DDRTree_2 = as.numeric(reduced[2, ])
)
trajectory_source <- merge(
  trajectory_source,
  pd_final[, .(
    CellID, Cluster, Disease, SampleID_Global, DatasetID,
    State, Pseudotime
  )],
  by = "CellID", all.x = TRUE, sort = FALSE
)
stopifnot(nrow(trajectory_source) == ncol(cds_default), !anyNA(trajectory_source$Pseudotime))
fwrite(trajectory_source, file.path(SOURCE_DIR, "Fig2i_monocle2_trajectory_source.csv"))

reduced_c00 <- reducedDimS(cds_c00_root)
pd_c00 <- as.data.table(pData(cds_c00_root), keep.rownames = "CellID")
trajectory_source_c00 <- data.table(
  CellID = colnames(reduced_c00),
  DDRTree_1 = as.numeric(reduced_c00[1, ]),
  DDRTree_2 = as.numeric(reduced_c00[2, ])
)
trajectory_source_c00 <- merge(
  trajectory_source_c00,
  pd_c00[, .(CellID, Cluster, Disease, SampleID_Global, DatasetID, State, Pseudotime)],
  by = "CellID", all.x = TRUE, sort = FALSE
)
fwrite(trajectory_source_c00, file.path(SOURCE_DIR, "Fig2i_monocle2_trajectory_source_C00_root.csv"))

state_cluster <- trajectory_source[
  , .(n_cells = .N), by = .(State, Cluster)
][
  , state_fraction := n_cells / sum(n_cells), by = State
]
fwrite(state_cluster, file.path(TABLE_DIR, "Fig2i_state_cluster_composition.csv"))

cluster_pseudotime <- trajectory_source[
  , .(
    n_cells = .N,
    median_pseudotime = median(Pseudotime),
    mean_pseudotime = mean(Pseudotime),
    q25_pseudotime = quantile(Pseudotime, 0.25),
    q75_pseudotime = quantile(Pseudotime, 0.75)
  ), by = Cluster
][order(median_pseudotime)]
fwrite(cluster_pseudotime, file.path(TABLE_DIR, "Fig2i_cluster_pseudotime_summary.csv"))

c00_root_check <- trajectory_source_c00[
  Cluster == c00_label, .(min_pt = min(Pseudotime), median_pt = median(Pseudotime))
]
fwrite(
  data.table(
    Item = c("C00_canonical_label", "C00_in_primary_trajectory",
             "C00_sampled_cells", "C00_majority_State",
             "C00_root_min_pseudotime", "C00_root_median_pseudotime",
             "State_partition_unchanged_by_root_state"),
    Value = c(c00_label, c00_in_primary, c00_cells_in_primary, c00_root_state,
              c00_root_check$min_pt, c00_root_check$median_pt,
              identical(as.character(pData(cds_default)$State),
                        as.character(pData(cds_c00_root)$State)))
  ),
  file.path(TABLE_DIR, "Fig2i_C00_root_sensitivity_audit.csv")
)

###############################################################################
# 8. Figures: independent PNG/PDF panels
###############################################################################

cluster_palette <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
  "#E6AB02", "#A6761D", "#1F78B4", "#B2DF8A", "#FB9A99",
  "#CAB2D6", "#6A3D9A", "#B15928", "#17BECF"
)
cluster_colors <- setNames(cluster_palette[seq_along(trajectory_clusters)], trajectory_clusters)

base_theme <- theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title = element_text(colour = "black"),
    plot.title = element_text(face = "bold", colour = "black"),
    plot.subtitle = element_text(colour = "black"),
    legend.title = element_text(face = "bold")
  )

p_cluster <- plot_cell_trajectory(
  cds_default, color_by = "Cluster", show_branch_points = TRUE, cell_size = 0.45
) +
  scale_color_manual(values = cluster_colors, drop = FALSE) +
  labs(
    title = "Microglia subtype trajectory",
    subtitle = "Pan-disease subtypes; Monocle2 DDRTree; no root state imposed",
    x = "Component 1", y = "Component 2", colour = "Subcluster"
  ) + base_theme

p_time <- plot_cell_trajectory(
  cds_default, color_by = "Pseudotime", show_branch_points = TRUE, cell_size = 0.45
) +
  scale_color_viridis_c(option = "C", direction = 1) +
  labs(
    title = "Microglia pseudotime",
    subtitle = "Default Monocle2 orientation; no biological root was imposed",
    x = "Component 1", y = "Component 2", colour = "Pseudotime"
  ) + base_theme

p_state <- plot_cell_trajectory(
  cds_default, color_by = "State", show_branch_points = TRUE, cell_size = 0.45
) +
  labs(
    title = "Monocle2 trajectory states",
    subtitle = "Default orderCells state assignment",
    x = "Component 1", y = "Component 2", colour = "State"
  ) + base_theme

p_disease <- plot_cell_trajectory(
  cds_default, color_by = "Disease", show_branch_points = TRUE, cell_size = 0.35
) +
  scale_color_manual(values = disease_colors, drop = FALSE) +
  labs(
    title = "Disease distribution along the trajectory",
    subtitle = "Diagnostic only; disease is not used to construct or orient the trajectory",
    x = "Component 1", y = "Component 2", colour = "Disease"
  ) + base_theme

p_dataset <- plot_cell_trajectory(
  cds_default, color_by = "DatasetID", show_branch_points = TRUE, cell_size = 0.35
) +
  scale_color_viridis_d(option = "D", end = 0.92) +
  labs(
    title = "Dataset distribution along the trajectory",
    subtitle = "Diagnostic only; DatasetID is not used to construct the trajectory",
    x = "Component 1", y = "Component 2", colour = "Dataset"
  ) + base_theme

p_cluster_c00 <- plot_cell_trajectory(
  cds_c00_root, color_by = "Cluster", show_branch_points = TRUE, cell_size = 0.45
) +
  scale_color_manual(values = cluster_colors, drop = FALSE) +
  labs(title = "Microglia subtype trajectory (C00-root sensitivity)",
       subtitle = paste0("C00-associated State ", c00_root_state, " used for orientation"),
       x = "Component 1", y = "Component 2", colour = "Subcluster") + base_theme

p_time_c00 <- plot_cell_trajectory(
  cds_c00_root, color_by = "Pseudotime", show_branch_points = TRUE, cell_size = 0.45
) +
  scale_color_viridis_c(option = "C", direction = 1) +
  labs(title = "Microglia pseudotime (C00-root sensitivity)",
       subtitle = paste0("C00-associated State ", c00_root_state, " defined as pseudotime origin"),
       x = "Component 1", y = "Component 2", colour = "Pseudotime") + base_theme

save_panel <- function(plot_obj, stem, width = 8.2, height = 6.5) {
  ggsave(file.path(FIG_DIR, paste0(stem, ".png")), plot_obj,
         width = width, height = height, dpi = 600, bg = "white")
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), plot_obj,
         width = width, height = height, device = cairo_pdf, bg = "white")
}
save_panel(p_cluster, "Fig2i_monocle2_trajectory_by_cluster")
save_panel(p_time, "Fig2i_monocle2_trajectory_by_pseudotime")
save_panel(p_state, "Fig2i_monocle2_trajectory_by_state")
save_panel(p_disease, "Fig2i_monocle2_trajectory_by_disease", width = 9.2)
save_panel(p_dataset, "Fig2i_monocle2_trajectory_by_dataset", width = 9.2)
save_panel(p_cluster_c00, "Fig2i_monocle2_trajectory_by_cluster_C00_root")
save_panel(p_time_c00, "Fig2i_monocle2_trajectory_by_pseudotime_C00_root")

###############################################################################
# 9. Method audit, interpretation note, and session information
###############################################################################

method_audit <- data.table(
  Item = c(
    "Reference", "Reference_Monocle_version", "Installed_Monocle_version",
    "Input", "User_contaminants_excluded", "Trajectory_scope",
    "Included_clusters", "Sampling", "Cells_before_sampling",
    "Cells_after_sampling", "Per_cluster_cap", "Ordering_gene_detection",
    "Ordering_gene_test", "Ordering_gene_threshold", "Ordering_gene_count",
    "Reduction", "Pseudotime_orientation",
    "Checkpoint_reuse", "Checkpoint_note",
    "C00_actual_label", "C00_label_normalization",
    "C00_in_primary_trajectory", "C00_sampled_cells_in_primary",
    "C00_root_State", "C00_root_State_semantics",
    "Primary_vs_sensitivity", "State_guard", "Nei_patch",
    "Statistical_interpretation"
  ),
  Value = c(
    "Nie et al. Cell Discovery 2026; doi:10.1038/s41421-026-00888-3",
    EXPECTED_MONOCLE_VERSION, monocle_version, INPUT_RDS,
    paste(CONTAMINANT_CLUSTERS, collapse = ","), TRAJECTORY_SCOPE,
    paste(trajectory_clusters, collapse = ","), "cluster-balanced and sample-balanced",
    as.character(cells_before_sampling_used), as.character(cells_after_sampling_used),
    MAX_CELLS_PER_CLUSTER,
    ">=10% sampled cells", "differentialGeneTest(~Cluster)",
    paste0("P<", ORDERING_P_CUTOFF), as.character(ordering_gene_count_used),
    "DDRTree, 2 components", "Primary=orderCells() default; sensitivity=orderCells(root_state)",
    resume_checkpoint, checkpoint_note_full,
    c00_label, "equivalent labels 0/C0/C00 normalized to C00",
    c00_in_primary, c00_cells_in_primary,
    as.character(c00_root_state),
    paste0("State partition is a property of the principal graph and is NOT ",
           "rewritten by orderCells(root_state=); root_state only re-orients ",
           "Pseudotime. cds_C00_root therefore shares the primary State labels."),
    paste0("Primary: default orderCells() direction (no biological root imposed). ",
           "Sensitivity: pseudotime re-rooted at C00-majority State ",
           c00_root_state, "; only reading direction changes, clustering/states ",
           "and trajectory geometry are identical."),
    state_guard_note,
    "v2: bare .nei( + neimode->mode inside [.igraph.vs NSE context; .nei bound ",
    "via per-function closure environment (igraph 2.3.2 compatible)",
    "Descriptive trajectory; not proof of lineage, temporal development, or a universal disease progression axis"
  )
)
fwrite(method_audit, file.path(TABLE_DIR, "Fig2i_method_audit.csv"))

interpretation <- c(
  "# Fig2i Monocle2轨迹解读",
  "",
  "## 分析目的",
  "借鉴Nie等Cell Discovery 2026的Fig.2i，描述跨疾病Pan-disease小胶质亚型在转录连续体中的相对位置和分支结构。",
  "",
  "## 输入数据",
  paste0("正式Seurat对象：", INPUT_RDS),
  paste0("用户裁决污染簇：", paste(CONTAMINANT_CLUSTERS, collapse = ", "), "。"),
  paste0("正式轨迹纳入簇：", paste(trajectory_clusters, collapse = ", "), "。"),
  paste0("本次实际进入轨迹的细胞数：", ncol(cds_default), "。"),
  paste0("检查点复用：", checkpoint_note_full, "。"),
  "",
  "## 方法",
  paste0("Monocle ", monocle_version, "；RNA counts；estimateSizeFactors/estimateDispersions；表达于至少10%细胞且differentialGeneTest(~Cluster) P<0.01的基因作为ordering genes；DDRTree二维降维；orderCells排序。"),
  paste0("兼容性补丁：dplyr group_by_/select_ 重定向；igraph dfs(neimode)→dfs(mode)；monocle nei()→.nei()（v2：裸 .nei 调用 + mode 参数 + 闭包环境绑定，兼容 igraph 2.3.2 的 [.igraph.vs 非标准求值）。"),
  paste0("C00 标签归一化：接受 0/C0/C00 等价写法，统一为 ", c00_label, "。"),
  paste0("主分析：orderCells() 默认方向（不人为指定起点）。"),
  paste0("敏感性分析：orderCells(root_state=", c00_root_state, ")，即 C00 细胞最多的 State。"),
  "两种方向均仅用于描述，不能解释为已证明的稳态到激活态方向。",
  "",
  "## C00 与主轨迹/敏感性分析的关系",
  paste0("C00（", c00_label, "）是 Pan-disease 型、非污染簇，包含在主轨迹中，抽样 ", c00_cells_in_primary, " 个细胞。"),
  paste0("C00 主要位于 State ", c00_root_state, "；敏感性分析以该 State 为拟时序起点。"),
  "State 划分由主图（principal graph）结构决定，与 root 选择无关：orderCells(root_state=) 只重定向 Pseudotime，不改变 State 标签，因此主分析与敏感性分析共享同一 State 划分。",
  "主分析结果（Fig2i_monocle2_cds_default.rds / Fig2i_monocle2_cds.rds）不因敏感性分析被覆盖；C00-root 仅作为方向性稳健性检查（Fig2i_monocle2_cds_C00_root.rds）。",
  "",
  "## 结果读取",
  "主图按subcluster着色；pseudotime图用于读取方向；State图用于分支审计；Disease图仅检查疾病或数据来源是否主导轨迹；C00-root 两图与主图逐细胞对照读取方向差异。",
  "",
  "## 生物学意义与限制",
  paste0("C00 为用户确认的稳态簇；C00-root 敏感性分析检验“默认方向是否把稳态端放在起点”。若默认方向与 C00-root 方向一致，说明方向对起点选择稳健；若相反，主图方向只能视为 Monocle2 的默认取向而非生物学结论。"),
  "该结果是转录相似性排序，不证明真实细胞谱系、不可逆发育过程或所有疾病共享同一时间轴。若轨迹由单一Disease/DatasetID主导，应降低结论强度或进行分疾病/leave-one-dataset-out敏感性分析。"
)
md_path <- file.path(OUT_ROOT, "Fig2i_解读.md")
Encoding(md_path) <- "UTF-8"   # 脚本源为UTF-8字节且C locale下被标记为native,显式纠正
writeLines(interpretation, md_path, useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(TABLE_DIR, "Fig2i_sessionInfo.txt"), useBytes = TRUE)

message(
  "Fig2i完成：", ncol(cds_default), " cells; ", length(trajectory_clusters),
  " clusters; ordering genes=", ordering_gene_count_used,
  "; C00-root State=", c00_root_state, "."
)
