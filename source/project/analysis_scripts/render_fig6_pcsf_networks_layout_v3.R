#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(ggplot2)
  library(ggrepel)
  library(ragg)
  library(digest)
  library(jsonlite)
  library(grid)
})
if (!nzchar(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.65001"))) {
  stop("R无法切换到支持中文路径的Windows locale。")
}

ROOT <- "E:/micro_output/分析审阅中心/待审阅/Fig6_疾病vsCtrl_SHAP_PCSF模块"
NETWORK_DIR <- file.path(ROOT, "networks")
FIG_DIR <- file.path(ROOT, "figures", "network_layout_v3")
SOURCE_DIR <- file.path(FIG_DIR, "source_data")
CONFIG_DIR <- file.path(ROOT, "config", "network_layout_v3")
LOG_DIR <- file.path(ROOT, "logs", "network_layout_v3")
SCRIPT_PATH <- "D:/data/microglia_new/analysis_scripts/render_fig6_pcsf_networks_layout_v3.R"
SCRIPT_COPY <- file.path(ROOT, "scripts", basename(SCRIPT_PATH))
RESUME <- "--resume" %chin% commandArgs(trailingOnly = TRUE)
DISEASES <- c("AD", "ALS", "ALS_FTD", "DLB", "DLBD_PD", "EOAD", "FTD", "FTD_MND",
              "GBM", "HD", "LCBM", "LGG", "MB", "MS", "PD", "PSP")
EXPECTED <- data.table(
  Disease = DISEASES,
  nodes = c(89L, 52L, 0L, 92L, 130L, 0L, 0L, 17L, 30L, 0L, 88L, 153L, 138L, 103L, 81L, 229L),
  forest_edges = c(88L, 51L, 0L, 91L, 129L, 0L, 0L, 16L, 29L, 0L, 87L, 151L, 137L, 102L, 80L, 228L)
)
PALETTE <- c("#4E79A7", "#F28E2B", "#59A14F", "#B07AA1", "#76B7B2", "#E15759",
             "#EDC948", "#FF9DA7", "#9C755F", "#79706E", "#8DA0CB", "#66C2A5")
names(PALETTE) <- sprintf("M%03d", seq_along(PALETTE))
NODE_FILL <- "#D9D9D9"
NODE_STROKE <- "#595959"
CROSS_COLOUR <- "#858585"
SEED <- 42L
LAYOUT_ITERATIONS <- c(4000L, 8000L, 12000L)
PLOT_WIDTH_MM <- 190
DEFAULT_HEIGHT_MM <- 165
DENSE_HEIGHT_MM <- 190
sha256_file <- function(path) as.character(digest::digest(file = path, algo = "sha256", serialize = FALSE))

if (!dir.exists(ROOT) || !dir.exists(NETWORK_DIR)) stop("正式分析根目录或networks目录不存在：", ROOT)
OUTPUT_EXISTS <- any(dir.exists(c(FIG_DIR, SOURCE_DIR, CONFIG_DIR, LOG_DIR)))
if (OUTPUT_EXISTS && !RESUME) {
  stop("v3输出目录已存在；为保护旧结果，本脚本拒绝覆盖：", FIG_DIR)
}
required_packages <- c("data.table", "igraph", "ggplot2", "ggrepel", "ragg", "digest", "jsonlite", "systemfonts")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("缺少绘图依赖；本任务不安装软件：", paste(missing_packages, collapse = ", "))
if (!capabilities("cairo")) stop("当前R不支持Cairo PDF设备；不替换为位图PDF。")
font_match <- systemfonts::match_font("Arial")
if (!file.exists(font_match$path)) stop("Arial字体未解析到本机字体文件。")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CONFIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(SCRIPT_COPY), recursive = TRUE, showWarnings = FALSE)
if (file.exists(SCRIPT_COPY) && !identical(sha256_file(SCRIPT_PATH), copy_hash <- sha256_file(SCRIPT_COPY))) {
  archive_dir <- file.path(dirname(SCRIPT_COPY), "archive_network_layout_v3")
  dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
  archived_script <- file.path(archive_dir, paste0(basename(SCRIPT_COPY), ".", substr(copy_hash, 1L, 10L), ".R"))
  if (!file.exists(archived_script) && !file.copy(SCRIPT_COPY, archived_script, overwrite = FALSE)) {
    stop("无法保留首次渲染脚本副本。")
  }
  if (!file.copy(SCRIPT_PATH, SCRIPT_COPY, overwrite = TRUE)) stop("无法更新渲染脚本副本。")
} else if (!file.exists(SCRIPT_COPY)) {
  if (!file.copy(SCRIPT_PATH, SCRIPT_COPY, overwrite = FALSE)) stop("无法保存渲染脚本副本。")
}

LOG_FILE <- file.path(LOG_DIR, "render_layout_v3.log")
log_msg <- function(...) {
  line <- paste0("[", format(Sys.time(), "%F %T %z"), "] ", paste0(..., collapse = ""), "\n")
  cat(line, file = LOG_FILE, append = TRUE)
  cat(line)
}
read_table <- function(path) data.table::fread(path, encoding = "UTF-8", na.strings = c("", "NA"), data.table = TRUE)
edge_key <- function(a, b) {
  an <- suppressWarnings(as.numeric(as.character(a)))
  bn <- suppressWarnings(as.numeric(as.character(b)))
  left <- ifelse(an <= bn, as.character(a), as.character(b))
  right <- ifelse(an <= bn, as.character(b), as.character(a))
  paste(left, right, sep = "|")
}
module_number <- function(x) as.integer(sub(".*_M", "", as.character(x)))
parse_solver_flag <- function(x, disease) {
  raw <- tolower(trimws(as.character(x)))
  yes <- raw %chin% c("true", "t", "1", "yes")
  no <- raw %chin% c("false", "f", "0", "no")
  if (any(!(yes | no))) stop(disease, " in_solver_forest存在未识别值：", paste(unique(raw[!(yes | no)]), collapse = ","))
  yes
}
manifest_output_hash <- function(manifest, filename) {
  outs <- manifest$outputs
  if (is.null(outs) || is.null(outs$path) || is.null(outs$sha256)) return(NA_character_)
  idx <- which(basename(as.character(outs$path)) == filename)
  if (length(idx) != 1L) return(NA_character_)
  as.character(outs$sha256[[idx]])
}

audit_disease <- function(disease) {
  disease_dir <- file.path(NETWORK_DIR, disease)
  paths <- c(
    nodes = file.path(disease_dir, "selected_nodes.csv"),
    forest = file.path(disease_dir, "selected_forest_edges.csv"),
    induced = file.path(disease_dir, "selected_induced_ppi_edges.csv"),
    membership = file.path(disease_dir, "module_membership.csv"),
    module_summary = file.path(disease_dir, "module_summary.csv"),
    manifest = file.path(disease_dir, "run_manifest.json")
  )
  if (!all(file.exists(paths))) stop(disease, "正式输入缺失：", paste(names(paths)[!file.exists(paths)], collapse = ", "))
  n <- read_table(paths[["nodes"]])
  e <- read_table(paths[["forest"]])
  ind <- read_table(paths[["induced"]])
  mem <- read_table(paths[["membership"]])
  ms <- read_table(paths[["module_summary"]])
  manifest <- jsonlite::fromJSON(paths[["manifest"]], simplifyVector = TRUE)
  expected <- EXPECTED[Disease == disease]
  required <- list(
    nodes = c("Disease", "Entrez", "Symbol", "prize", "mean_abs_SHAP", "mean_abs_SHAP_normalized", "module_id", "rank"),
    forest = c("Disease", "Protein_A_Entrez_ID", "Protein_B_Entrez_ID"),
    induced = c("Disease", "Protein_A_Entrez_ID", "Protein_B_Entrez_ID", "in_solver_forest"),
    membership = c("Disease", "Entrez", "module_id"),
    module_summary = c("Disease", "module_id", "node_count", "edge_count")
  )
  frames <- list(nodes = n, forest = e, induced = ind, membership = mem, module_summary = ms)
  for (nm in names(required)) {
    missing <- setdiff(required[[nm]], names(frames[[nm]]))
    if (length(missing)) stop(disease, " ", nm, "缺字段：", paste(missing, collapse = ", "))
  }
  if (!identical(as.character(manifest$Disease), disease) || !identical(as.character(manifest$counts$status), "complete")) {
    stop(disease, " run_manifest疾病身份或状态不符。")
  }
  if (any(n$Disease != disease) || any(e$Disease != disease) || any(ind$Disease != disease) ||
      any(mem$Disease != disease) || any(ms$Disease != disease)) stop(disease, "输入中混有其他Disease标签。")
  n[, Entrez := as.character(Entrez)]
  mem[, Entrez := as.character(Entrez)]
  if (anyNA(n$Entrez) || anyNA(suppressWarnings(as.numeric(n$Entrez))) || anyDuplicated(n$Entrez)) {
    stop(disease, "正式节点Entrez缺失、非数字或重复。")
  }
  if (anyDuplicated(mem$Entrez) || nrow(n) != nrow(mem) || !setequal(n$Entrez, mem$Entrez)) {
    stop(disease, "selected_nodes与module_membership不是一一对应。")
  }
  mem <- mem[match(n$Entrez, mem$Entrez)]
  if (any(as.character(n$module_id) != as.character(mem$module_id))) stop(disease, "冻结module_id与成员表不一致。")
  n[, `:=`(entrez_num = as.numeric(Entrez), prize = as.numeric(prize), rank_num = suppressWarnings(as.numeric(rank)))]
  if (any(!is.finite(n$prize)) || any(n$prize < 0)) stop(disease, "prize含非有限值或负值。")
  mods <- sort(unique(as.character(n$module_id)), method = "radix")
  mod_n <- module_number(mods)
  if (length(mods) && (anyNA(mod_n) || any(mod_n < 1L | mod_n > length(PALETTE)) || anyDuplicated(mod_n))) {
    stop(disease, "冻结模块编号超出颜色表范围或重复。")
  }
  if (length(mods) != as.integer(manifest$counts$module_count)) stop(disease, "模块数与run_manifest不一致。")
  if (nrow(n) != as.integer(manifest$counts$selected_node_count) ||
      nrow(e) != as.integer(manifest$counts$forest_edge_count) ||
      nrow(ind) != as.integer(manifest$counts$induced_ppi_edge_count)) stop(disease, "节点或边行数与run_manifest不一致。")
  if (nrow(n) != expected$nodes || nrow(e) != expected$forest_edges) {
    stop(disease, "实际正式结果与用户提供的当前核对值不符；实际节点/骨架边=", nrow(n), "/", nrow(e),
         "，参考值=", expected$nodes, "/", expected$forest_edges, "。本病已暂停出图。")
  }
  if (nrow(ms) != length(mods) || (nrow(ms) && sum(as.numeric(ms$node_count), na.rm = TRUE) != nrow(n)) ||
      !setequal(as.character(ms$module_id), mods)) stop(disease, "module_summary未覆盖冻结节点/模块。")

  e[, `:=`(from = as.character(Protein_A_Entrez_ID), to = as.character(Protein_B_Entrez_ID))]
  ind[, `:=`(from = as.character(Protein_A_Entrez_ID), to = as.character(Protein_B_Entrez_ID),
             in_solver_forest_bool = parse_solver_flag(in_solver_forest, disease))]
  if (nrow(e)) {
    e[, `:=`(from_num = as.numeric(from), to_num = as.numeric(to))]
    if (anyNA(e$from_num) || anyNA(e$to_num) || any(e$from_num == e$to_num) ||
        any(!e$from %chin% n$Entrez) || any(!e$to %chin% n$Entrez)) stop(disease, "骨架边有自环或端点不在正式节点中。")
    if (anyDuplicated(edge_key(e$from, e$to))) stop(disease, "正式骨架边包含重复无向边。")
  }
  if (nrow(ind) && (anyNA(suppressWarnings(as.numeric(ind$from))) || anyNA(suppressWarnings(as.numeric(ind$to))) ||
                    any(ind$from == ind$to) || any(!ind$from %chin% n$Entrez) || any(!ind$to %chin% n$Entrez))) {
    stop(disease, "诱导边有自环或端点不在正式节点中。")
  }
  forest_keys <- sort(edge_key(e$from, e$to), method = "radix")
  selected_induced_keys <- sort(unique(edge_key(ind[in_solver_forest_bool == TRUE, from],
                                                   ind[in_solver_forest_bool == TRUE, to])), method = "radix")
  if (!identical(forest_keys, selected_induced_keys)) stop(disease, "骨架边与诱导边in_solver_forest标记不一致。")
  if (nrow(ind) != as.integer(manifest$counts$induced_ppi_edge_count)) stop(disease, "诱导边行数不一致。")
  if (nrow(e)) {
    if (anyDuplicated(edge_key(e$from, e$to))) stop(disease, "标准化无向骨架边重复。")
  }
  for (nm in c("nodes", "forest", "induced", "membership", "module_summary")) {
    expected_hash <- manifest_output_hash(manifest, basename(paths[[nm]]))
    if (is.na(expected_hash) || !identical(tolower(sha256_file(paths[[nm]])), tolower(expected_hash))) {
      stop(disease, "正式输入SHA256与run_manifest不符：", nm)
    }
  }
  if (!nrow(n) && (nrow(e) || length(mods))) stop(disease, "空候选集却带有正式边或模块。")
  if (nrow(n) && nrow(e) != nrow(n) - length(igraph::components(
      igraph::graph_from_data_frame(as.data.frame(e[, .(from, to)]), directed = FALSE,
                                    vertices = data.frame(name = n$Entrez)))$csize)) {
    stop(disease, "正式骨架不是分量内森林（边数不满足V-C）。")
  }
  list(disease = disease, paths = paths, nodes = n, forest = e, induced = ind,
       membership = mem, module_summary = ms, manifest = manifest, modules = mods)
}

orient_pca <- function(xy) {
  xy <- as.matrix(xy)
  if (nrow(xy) <= 2L || !all(is.finite(xy))) return(xy)
  pc <- stats::prcomp(xy, center = TRUE, scale. = FALSE)
  z <- pc$x[, 1:2, drop = FALSE]
  for (j in 1:2) {
    k <- which.max(abs(pc$rotation[, j]))
    if (pc$rotation[k, j] < 0) z[, j] <- -z[, j]
  }
  rownames(z) <- rownames(xy)
  z
}

layout_component <- function(g, ids, iterations) {
  ids <- ids[order(as.numeric(ids), method = "radix")]
  if (length(ids) == 1L) return(matrix(c(0, 0), nrow = 1L, dimnames = list(ids, c("x", "y"))))
  if (length(ids) == 2L) return(matrix(c(-0.5, 0, 0.5, 0), nrow = 2L, byrow = TRUE,
                                       dimnames = list(ids, c("x", "y"))))
  subg <- igraph::induced_subgraph(g, vids = match(ids, igraph::V(g)$name))
  set.seed(SEED)
  xy <- igraph::layout_with_fr(subg, niter = iterations, grid = "nogrid",
                               weights = rep(1, igraph::ecount(subg)))
  rownames(xy) <- igraph::V(subg)$name
  orient_pca(xy[ids, , drop = FALSE])
}

pack_components <- function(g, ids, iterations, gap = 4) {
  comp <- igraph::components(g)
  memb <- setNames(comp$membership, igraph::V(g)$name)
  groups <- split(igraph::V(g)$name, memb)
  min_id <- vapply(groups, function(z) min(as.numeric(z)), numeric(1))
  groups <- groups[order(-lengths(groups), min_id)]
  locals <- lapply(groups, function(z) layout_component(g, z, iterations))
  boxes <- lapply(locals, function(z) c(xmin = min(z[, 1]), xmax = max(z[, 1]),
                                       ymin = min(z[, 2]), ymax = max(z[, 2])))
  centers <- vector("list", length(groups))
  placed <- vector("list", length(groups))
  if (length(groups)) {
    centers[[1L]] <- c(x = 0, y = 0)
    b <- boxes[[1L]]
    placed[[1L]] <- c(xmin = unname(b["xmin"]), xmax = unname(b["xmax"]),
                      ymin = unname(b["ymin"]), ymax = unname(b["ymax"]))
  }
  if (length(groups) > 1L) for (i in 2:length(groups)) {
    b <- boxes[[i]]
    w <- b["xmax"] - b["xmin"]
    h <- b["ymax"] - b["ymin"]
    cand <- matrix(numeric(), ncol = 2L, dimnames = list(NULL, c("x", "y")))
    for (old in placed[seq_len(i - 1L)]) {
      ox <- (old["xmin"] + old["xmax"]) / 2
      oy <- (old["ymin"] + old["ymax"]) / 2
      ow <- old["xmax"] - old["xmin"]
      oh <- old["ymax"] - old["ymin"]
      cand <- rbind(cand, c(ox - ow/2 - gap - w/2, oy), c(ox + ow/2 + gap + w/2, oy),
                    c(ox, oy - oh/2 - gap - h/2), c(ox, oy + oh/2 + gap + h/2),
                    c(ox - ow/2 - gap - w/2, oy - oh/2 - gap - h/2),
                    c(ox + ow/2 + gap + w/2, oy - oh/2 - gap - h/2),
                    c(ox - ow/2 - gap - w/2, oy + oh/2 + gap + h/2),
                    c(ox + ow/2 + gap + w/2, oy + oh/2 + gap + h/2))
    }
    cand <- unique(cand)
    legal <- vector("list", nrow(cand))
    scores <- rep(Inf, nrow(cand))
    for (j in seq_len(nrow(cand))) {
      cx <- unname(cand[j, 1L]); cy <- unname(cand[j, 2L])
      nb <- c(xmin = cx + unname(b["xmin"]), xmax = cx + unname(b["xmax"]),
              ymin = cy + unname(b["ymin"]), ymax = cy + unname(b["ymax"]))
      overlap <- vapply(placed[seq_len(i - 1L)], function(old) {
        !(nb["xmax"] + gap/2 <= old["xmin"] || old["xmax"] + gap/2 <= nb["xmin"] ||
          nb["ymax"] + gap/2 <= old["ymin"] || old["ymax"] + gap/2 <= nb["ymin"])
      }, logical(1))
      if (any(overlap)) next
      allb <- c(placed[seq_len(i - 1L)], list(nb))
      xmin <- min(vapply(allb, `[[`, numeric(1), "xmin")); xmax <- max(vapply(allb, `[[`, numeric(1), "xmax"))
      ymin <- min(vapply(allb, `[[`, numeric(1), "ymin")); ymax <- max(vapply(allb, `[[`, numeric(1), "ymax"))
      area <- max(1e-8, (xmax - xmin) * (ymax - ymin))
      aspect_penalty <- abs(log(max(1e-8, (xmax - xmin) / max(1e-8, ymax - ymin))))
      center_penalty <- (cx^2 + cy^2) * 0.015
      scores[j] <- area + aspect_penalty * area * 0.12 + center_penalty
      legal[[j]] <- nb
    }
    if (!any(is.finite(scores))) stop("Component packer could not place a non-overlapping bounding box.")
    j <- order(scores, cand[, 2L], cand[, 1L])[1L]
    centers[[i]] <- c(x = cand[j, 1L], y = cand[j, 2L])
    placed[[i]] <- legal[[j]]
  }
  xy <- matrix(NA_real_, nrow = length(ids), ncol = 2L, dimnames = list(ids, c("x", "y")))
  comp_id <- setNames(integer(length(ids)), ids)
  for (i in seq_along(groups)) {
    loc <- locals[[i]]
    xy[groups[[i]], ] <- sweep(loc, 2L, centers[[i]], "+")
    comp_id[groups[[i]]] <- i
  }
  xy <- sweep(xy, 2L, colMeans(xy), "-")
  list(xy = xy, comp_id = comp_id, groups = groups, boxes = placed,
       bbox_nonoverlap = TRUE)
}

count_crossings <- function(edges, xy) {
  if (nrow(edges) < 2L) return(0L)
  e <- edges[, .(from, to)]
  coords <- xy[match(e$from, rownames(xy)), , drop = FALSE]
  ends <- xy[match(e$to, rownames(xy)), , drop = FALSE]
  n <- nrow(e); total <- 0L
  for (i in seq_len(n - 1L)) {
    j <- (i + 1L):n
    shared <- e$from[j] == e$from[i] | e$from[j] == e$to[i] | e$to[j] == e$from[i] | e$to[j] == e$to[i]
    j <- j[!shared]
    if (!length(j)) next
    a <- coords[i, ]; b <- ends[i, ]
    c1 <- coords[j, , drop = FALSE]; d1 <- ends[j, , drop = FALSE]
    orient <- function(p, q, r) (q[1] - p[1]) * (r[, 2] - p[2]) - (q[2] - p[2]) * (r[, 1] - p[1])
    o1 <- orient(a, b, c1); o2 <- orient(a, b, d1)
    o3 <- mapply(function(k) (d1[k, 1] - c1[k, 1]) * (a[2] - c1[k, 2]) - (d1[k, 2] - c1[k, 2]) * (a[1] - c1[k, 1]), seq_len(nrow(c1)))
    o4 <- mapply(function(k) (d1[k, 1] - c1[k, 1]) * (b[2] - c1[k, 2]) - (d1[k, 2] - c1[k, 2]) * (b[1] - c1[k, 1]), seq_len(nrow(c1)))
    total <- total + sum((o1 * o2 < -1e-10) & (o3 * o4 < -1e-10))
  }
  as.integer(total)
}

layout_metrics <- function(edges, nodes, xy, fit_width, fit_height) {
  dx <- diff(range(xy[, 1])); dy <- diff(range(xy[, 2]))
  scale <- min(fit_width / max(dx, 1e-8), fit_height / max(dy, 1e-8))
  scaled <- xy * scale
  n <- nrow(nodes)
  overlap_count <- 0L; min_clearance <- Inf
  if (n > 1L) {
    for (i in seq_len(n - 1L)) {
      ddx <- scaled[(i + 1L):n, 1] - scaled[i, 1]
      ddy <- scaled[(i + 1L):n, 2] - scaled[i, 2]
      distance <- sqrt(ddx^2 + ddy^2)
      required <- (nodes$diameter_mm[(i + 1L):n] + nodes$diameter_mm[i]) / 2 + 0.6
      clearance <- distance - required
      min_clearance <- min(min_clearance, clearance)
      overlap_count <- overlap_count + sum(clearance < 0)
    }
  } else min_clearance <- Inf
  crossings <- count_crossings(edges, xy)
  if (nrow(edges)) {
    from <- scaled[match(edges$from, rownames(scaled)), , drop = FALSE]
    to <- scaled[match(edges$to, rownames(scaled)), , drop = FALSE]
    lens <- sqrt(rowSums((from - to)^2))
    med <- stats::median(lens)
    long_edges <- if (is.finite(med) && med > 0) sum(lens > 5 * med) else 0L
    edge_length_ratio <- if (is.finite(med) && med > 0) max(lens) / med else 0
  } else {
    long_edges <- 0L; edge_length_ratio <- 0
  }
  list(overlap_count = as.integer(overlap_count), crossing_count = crossings,
       long_edge_count = as.integer(long_edges), min_node_clearance_mm = min_clearance,
       scale_fit = scale, edge_length_max_to_median = edge_length_ratio)
}

choose_layout <- function(graph, nodes, edges, height_mm) {
  fit_width <- PLOT_WIDTH_MM - 30
  fit_height <- height_mm - 62
  trials <- vector("list", length(LAYOUT_ITERATIONS))
  for (i in seq_along(LAYOUT_ITERATIONS)) {
    niter <- LAYOUT_ITERATIONS[i]
    packed <- pack_components(graph, nodes$Entrez, niter)
    metrics <- layout_metrics(edges, nodes, packed$xy, fit_width, fit_height)
    aspect <- diff(range(packed$xy[, 1])) / max(1e-8, diff(range(packed$xy[, 2])))
    metrics$aspect_mismatch <- abs(log(max(1e-8, aspect) / (fit_width / fit_height)))
    metrics$niter <- niter
    trials[[i]] <- list(packed = packed, metrics = metrics)
  }
  score <- rbindlist(lapply(trials, function(x) as.data.table(x$metrics)), fill = TRUE)
  setorder(score, overlap_count, crossing_count, long_edge_count, aspect_mismatch, niter)
  winner <- trials[[match(score$niter[1L], LAYOUT_ITERATIONS)]]
  list(xy = winner$packed$xy, comp_id = winner$packed$comp_id,
       component_count = length(winner$packed$groups),
       component_bbox_nonoverlap = winner$packed$bbox_nonoverlap,
       selected_niter = winner$metrics$niter, metrics = winner$metrics,
       trials = score)
}

select_labels <- function(nodes, edges) {
  if (!nrow(nodes)) return(list(genes = nodes[0], module = data.table()))
  top <- nodes[order(-prize, rank_num, entrez_num, na.last = TRUE), .SD[1L], by = "module_id"]
  bridge_ids <- character()
  if (nrow(edges)) {
    node_module <- setNames(nodes$module_id, nodes$Entrez)
    a <- data.table(Entrez = c(edges$from, edges$to),
                    other_module = c(unname(node_module[edges$to]), unname(node_module[edges$from])))
    bridge_ids <- a[other_module != unname(node_module[Entrez]), .(n_other_modules = data.table::uniqueN(other_module)), by = Entrez][n_other_modules >= 2L, Entrez]
  }
  gene <- nodes[Entrez %chin% unique(c(top$Entrez, bridge_ids))]
  gene[, `:=`(label_type = "gene", label = fifelse(is.na(Symbol) | Symbol == "", Entrez, Symbol),
              label_reason = vapply(Entrez, function(id) paste(c(
                if (id %chin% top$Entrez) "highest_prize_in_module_rank_then_Entrez_tiebreak",
                if (id %chin% bridge_ids) "adjacent_to_at_least_two_distinct_other_modules"), collapse = ";"), character(1)),
              label_size_pt = 8.2, fontface = "plain", label_colour = "#000000")]
  nodes[, .(Entrez, module_id, layout_component_id, x, y)] -> temp
  module <- temp[, .(n_in_component = .N, min_entrez = min(as.numeric(Entrez)), x = mean(x), y = mean(y)),
                 by = .(module_id, layout_component_id)][order(-n_in_component, min_entrez), .SD[1L], by = "module_id"]
  module[, `:=`(label_type = "module", label = sprintf("M%02d", module_number(module_id)),
                 label_reason = "label_at_largest_real_component_containing_module_nodes",
                 label_size_pt = 7.0, fontface = "bold", label_colour = unname(PALETTE[sprintf("M%03d", module_number(module_id))]))]
  list(genes = gene, module = module)
}

write_plot_safe <- function(plot, png_path, pdf_path, width_mm, height_mm) {
  png_device <- FALSE
  tryCatch({
    ragg::agg_png(png_path, width = width_mm, height = height_mm, units = "mm", res = 600, background = "white")
    png_device <- TRUE
    print(plot)
  }, finally = {
    if (png_device && grDevices::dev.cur() > 1L) grDevices::dev.off()
  })
  pdf_device <- FALSE
  tryCatch({
    grDevices::cairo_pdf(pdf_path, width = width_mm / 25.4, height = height_mm / 25.4, family = "Arial", onefile = TRUE)
    pdf_device <- TRUE
    print(plot)
  }, finally = {
    if (pdf_device && grDevices::dev.cur() > 1L) grDevices::dev.off()
  })
}

render_one <- function(disease) {
  log_msg("START disease=", disease)
  audited <- audit_disease(disease)
  nodes <- copy(audited$nodes)
  edges <- copy(audited$forest)
  induced <- audited$induced
  mods <- audited$modules
  height_mm <- if (disease == "PSP" || nrow(nodes) >= 180L) DENSE_HEIGHT_MM else DEFAULT_HEIGHT_MM
  stem <- paste0(disease, "_PCSF_network_layout_v3")
  png_path <- file.path(FIG_DIR, paste0(stem, ".png"))
  pdf_path <- file.path(FIG_DIR, paste0(stem, ".pdf"))
  if ((file.exists(png_path) || file.exists(pdf_path)) && !RESUME) stop(disease, "新图已存在，拒绝覆盖。")
  colour_table <- data.table(Disease = disease, module_id = mods,
                             module_label = sprintf("M%02d", module_number(mods)),
                             colour = unname(PALETTE[sprintf("M%03d", module_number(mods))]))
  if (!nrow(nodes)) {
    empty_message <- paste("No candidate network returned", "under the frozen PCSF configuration", sep = "\n")
    secondary <- "This does not establish absence of disease-associated genes."
    plot <- ggplot() +
      annotate("text", x = 0, y = 0.2, label = empty_message, family = "Arial", size = 3.3, lineheight = 1.15) +
      annotate("text", x = 0, y = -0.32, label = secondary, family = "Arial", size = 2.35, colour = "#555555") +
      labs(title = paste0(disease, "-associated candidate network"),
           subtitle = "0 genes | 0 PCSF-selected interactions | 0 modules") +
      theme_void(base_family = "Arial") +
      theme(plot.background = element_rect(fill = "white", colour = NA),
            plot.title = element_text(size = 11, face = "bold", hjust = 0, colour = "black"),
            plot.subtitle = element_text(size = 8, hjust = 0, colour = "black"),
            plot.title.position = "plot", plot.margin = margin(7, 8, 8, 8))
    node_source <- data.table(Disease = character(), Entrez = character(), Symbol = character(), prize = numeric(),
                              module_id = character(), node_diameter_mm = numeric(), x = numeric(), y = numeric(),
                              label_status = character(), label_reason = character(), layout_component_id = integer())
    edge_source <- data.table(Disease = character(), from = character(), to = character(), module_from = character(),
                              module_to = character(), edge_class = character(), edge_colour = character(),
                              x = numeric(), y = numeric(), xend = numeric(), yend = numeric())
    label_source <- data.table(Disease = character(), label_type = character(), label = character(),
                               label_reason = character(), anchor_x = numeric(), anchor_y = numeric(),
                               label_colour = character(), fontface = character(), font_size_pt = numeric(),
                               manual_dx_mm = numeric(), manual_dy_mm = numeric(), offset_status = character())
    layout_config <- list(Disease = disease, status = "empty_candidate_set_not_assessed",
                          empty_state_text = c(empty_message, secondary), layout = "not applicable",
                          png_dpi = 600, PDF = "single-page Cairo vector text; Arial")
    write_plot_safe(plot, png_path, pdf_path, PLOT_WIDTH_MM, DEFAULT_HEIGHT_MM)
  } else {
    max_prize <- max(nodes$prize)
    nodes[, node_diameter_mm := if (is.finite(max_prize) && max_prize > 0) 1.2 + 3 * prize / max_prize else 1.2]
    graph_vertices <- as.data.frame(nodes[, .(name = Entrez, module_id, prize)])
    nodes <- nodes[order(entrez_num, method = "radix")]
    edges[, `:=`(from = as.character(Protein_A_Entrez_ID), to = as.character(Protein_B_Entrez_ID))]
    if (nrow(edges)) {
      edges[, `:=`(from_num = as.numeric(from), to_num = as.numeric(to))]
      edges[, `:=`(left_num = pmin(from_num, to_num), right_num = pmax(from_num, to_num))]
      edges <- edges[order(left_num, right_num)]
      edges[, c("from_num", "to_num", "left_num", "right_num") := NULL]
    } else edges[, `:=`(from = character(), to = character())]
    graph <- igraph::graph_from_data_frame(as.data.frame(edges[, .(from, to)]), directed = FALSE,
                                           vertices = graph_vertices)
    if (igraph::vcount(graph) != nrow(nodes) || igraph::ecount(graph) != nrow(edges)) stop(disease, "布局图节点/边数与正式表不符。")
    layout <- choose_layout(graph, nodes, edges, height_mm)
    nodes[, `:=`(x = layout$xy[Entrez, 1L], y = layout$xy[Entrez, 2L],
                 layout_component_id = as.integer(layout$comp_id[Entrez]),
                 module_label = sprintf("M%02d", module_number(module_id)),
                 module_colour = unname(PALETTE[sprintf("M%03d", module_number(module_id))]))]
    module_by_node <- setNames(nodes$module_id, nodes$Entrez)
    if (nrow(edges)) {
      edges[, `:=`(module_from = unname(module_by_node[from]), module_to = unname(module_by_node[to]))]
      edges[, edge_class := ifelse(module_from == module_to, "within-module", "between-module")]
      edges[, edge_colour := ifelse(edge_class == "within-module",
                                    unname(PALETTE[sprintf("M%03d", module_number(module_from))]), CROSS_COLOUR)]
      edges[, `:=`(x = nodes$x[match(from, nodes$Entrez)], y = nodes$y[match(from, nodes$Entrez)],
                   xend = nodes$x[match(to, nodes$Entrez)], yend = nodes$y[match(to, nodes$Entrez)])]
      if (any(!is.finite(edges$x)) || any(!is.finite(edges$y)) || any(!is.finite(edges$xend)) || any(!is.finite(edges$yend))) {
        stop(disease, "绘图边坐标不是有限值。")
      }
    } else {
      edges[, `:=`(module_from = character(), module_to = character(), edge_class = character(), edge_colour = character(),
                   x = numeric(), y = numeric(), xend = numeric(), yend = numeric())]
    }
    selected_labels <- select_labels(nodes, edges)
    genes <- selected_labels$genes
    module_labels <- selected_labels$module
    # One shared repel layer contains every node; unlabelled nodes retain empty labels as obstacles.
    repel_nodes <- nodes[, .(Disease = disease, Entrez, label = "", label_type = "node_obstacle",
                             label_reason = "all_formal_nodes_are_repel_obstacles", x, y,
                             label_colour = "#000000", fontface = "plain", font_size_pt = 0.5)]
    if (nrow(genes)) {
      gene_labels <- genes[, .(Disease = disease, Entrez, label, label_type, label_reason, x, y,
                               label_colour, fontface, font_size_pt = label_size_pt)]
      repel_nodes[match(gene_labels$Entrez, repel_nodes$Entrez), `:=`(
        label = gene_labels$label, label_type = "gene", label_reason = gene_labels$label_reason,
        label_colour = "#000000", fontface = "plain", font_size_pt = gene_labels$font_size_pt)]
    } else gene_labels <- data.table()
    if (nrow(module_labels)) {
      module_label_rows <- module_labels[, .(Disease = disease, Entrez = NA_character_, label, label_type,
                                              label_reason, x, y, label_colour, fontface,
                                              font_size_pt = label_size_pt)]
      repel_nodes <- rbind(repel_nodes, module_label_rows, fill = TRUE)
    } else module_label_rows <- data.table()
    labelled <- rbindlist(list(gene_labels, module_label_rows), use.names = TRUE, fill = TRUE)
    label_source <- labelled[, .(Disease, label_type, label, label_reason, anchor_x = x, anchor_y = y,
                                 label_colour, fontface, font_size_pt, manual_dx_mm = 0,
                                 manual_dy_mm = 0, offset_status = "automatic_joint_ggrepel_seed42")]
    legend_prizes <- if (max_prize > 0) unique(round(stats::quantile(nodes$prize, c(0, .5, 1), names = FALSE, na.rm = TRUE), 4)) else 0
    legend_diameters <- if (max_prize > 0) 1.2 + 3 * legend_prizes / max_prize else rep(1.2, length(legend_prizes))
    xspan_data <- diff(range(nodes$x)); yspan_data <- diff(range(nodes$y))
    fit_width <- PLOT_WIDTH_MM - 30; fit_height <- height_mm - 62
    fitted_scale <- min(fit_width / max(xspan_data, 1e-8), fit_height / max(yspan_data, 1e-8))
    nodes[, `:=`(x = x * fitted_scale, y = y * fitted_scale)]
    edges[, `:=`(x = x * fitted_scale, y = y * fitted_scale, xend = xend * fitted_scale, yend = yend * fitted_scale)]
    repel_nodes[, `:=`(x = x * fitted_scale, y = y * fitted_scale)]
    label_source[, `:=`(anchor_x = anchor_x * fitted_scale, anchor_y = anchor_y * fitted_scale)]
    module_labels[, `:=`(x = x * fitted_scale, y = y * fitted_scale)]
    gene_labels[, `:=`(x = x * fitted_scale, y = y * fitted_scale)]
    # Center the actual network bounds before applying symmetric coordinate limits.
    # Packed layouts are not guaranteed to be centered at (0, 0); without this
    # translation, large connected components can be cropped at a device edge.
    center_x <- mean(range(nodes$x))
    center_y <- mean(range(nodes$y))
    nodes[, `:=`(x = x - center_x, y = y - center_y)]
    edges[, `:=`(x = x - center_x, y = y - center_y,
                 xend = xend - center_x, yend = yend - center_y)]
    repel_nodes[, `:=`(x = x - center_x, y = y - center_y)]
    label_source[, `:=`(anchor_x = anchor_x - center_x, anchor_y = anchor_y - center_y)]
    module_labels[, `:=`(x = x - center_x, y = y - center_y)]
    gene_labels[, `:=`(x = x - center_x, y = y - center_y)]
    # Split plotting layers only after the same scale and translation have been
    # applied to every edge; earlier data.table subsets would retain stale coords.
    edge_within <- edges[edge_class == "within-module"]
    edge_cross <- edges[edge_class == "between-module"]
    xspan <- diff(range(nodes$x)) + 10
    yspan <- diff(range(nodes$y)) + 10
    body_h <- height_mm - 58
    target_aspect <- (PLOT_WIDTH_MM - 18) / body_h
    xspan <- max(xspan, yspan * target_aspect)
    yspan <- xspan / target_aspect
    xlim <- c(-xspan / 2, xspan / 2)
    ylim <- c(-yspan / 2, yspan / 2)
    plot <- ggplot() +
      geom_segment(data = edge_cross, aes(x = x, y = y, xend = xend, yend = yend,
                                          linetype = "Between-module interactions"),
                   colour = CROSS_COLOUR, linewidth = .25, alpha = .65, lineend = "round") +
      geom_segment(data = edge_within, aes(x = x, y = y, xend = xend, yend = yend, colour = edge_colour),
                   linewidth = .30, alpha = .90, lineend = "round", show.legend = TRUE) +
      geom_point(data = nodes, aes(x = x, y = y, size = node_diameter_mm), shape = 21,
                 stroke = .20, fill = NODE_FILL, colour = NODE_STROKE, show.legend = TRUE) +
      ggrepel::geom_text_repel(data = repel_nodes,
        aes(x = x, y = y, label = label, colour = label_colour, fontface = fontface,
            size = font_size_pt / ggplot2::.pt),
        family = "Arial", box.padding = grid::unit(.7, "mm"), point.padding = grid::unit(.8, "mm"),
        max.overlaps = Inf, max.iter = 100000, max.time = 10, min.segment.length = 0,
        segment.colour = "#999999", segment.size = .15, seed = SEED,
        show.legend = FALSE, na.rm = TRUE) +
      scale_colour_identity(name = "Within-module interactions",
        breaks = unname(colour_table$colour), labels = colour_table$module_label, guide = "legend") +
      scale_linetype_manual(name = NULL, values = c("Between-module interactions" = "solid"),
        breaks = "Between-module interactions", guide = "legend") +
      scale_size_identity(name = "Relative SHAP importance", breaks = legend_diameters,
        labels = format(legend_prizes, trim = TRUE, scientific = FALSE), guide = "legend") +
      guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                    override.aes = list(linewidth = .8, alpha = 1, shape = NA)),
             linetype = guide_legend(override.aes = list(colour = CROSS_COLOUR, linewidth = .8, shape = NA)),
             size = guide_legend(override.aes = list(shape = 21, fill = NODE_FILL, colour = NODE_STROKE, stroke = .2))) +
      coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE, clip = "off") +
      labs(title = paste0(disease, "-associated candidate network"),
           subtitle = sprintf("%d genes | %d PCSF-selected interactions | %d modules", nrow(nodes), nrow(edges), length(mods))) +
      theme_void(base_family = "Arial") +
      theme(plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white", colour = NA),
            plot.title = element_text(size = 11, face = "bold", hjust = 0, colour = "black"),
            plot.subtitle = element_text(size = 8, hjust = 0, colour = "black"), plot.title.position = "plot",
            legend.position = "bottom", legend.box = "vertical", legend.title = element_text(size = 7),
            legend.text = element_text(size = 7), legend.key.width = grid::unit(7, "mm"),
            plot.margin = margin(5, 8, 6, 8))
    write_plot_safe(plot, png_path, pdf_path, PLOT_WIDTH_MM, height_mm)
    node_source <- nodes[, .(Disease = disease, Entrez, Symbol, mean_abs_SHAP, mean_abs_SHAP_normalized,
                             prize, module_id, module_label, module_colour, node_diameter_mm,
                             x, y, layout_component_id,
                             label_status = ifelse(Entrez %chin% genes$Entrez, "labelled", "not_selected_by_frozen_label_rule"),
                             label_reason = ifelse(Entrez %chin% genes$Entrez,
                                                  genes$label_reason[match(Entrez, genes$Entrez)],
                                                  "not_selected_by_frozen_label_rule"))]
    edge_source <- edges[, .(Disease = disease, from, to, module_from, module_to, edge_class, edge_colour,
                             x, y, xend, yend)]
    layout_config <- list(
      Disease = disease, status = "rendered", layout_method = "igraph FR on formal PCSF skeleton; connected components only",
      iterations_trials = as.list(layout$trials), selected_niter = layout$selected_niter,
      seed = SEED, pca_orientation = "prcomp rotation for horizontal major axis; display only",
      component_pack_order = "descending node count, then minimum numeric Entrez",
      component_bbox_nonoverlap = layout$component_bbox_nonoverlap,
      node_count = nrow(nodes), forest_edge_count = nrow(edges), induced_edge_count_audit_only = nrow(induced),
      module_count = length(mods), gene_label_count = nrow(genes), module_label_count = nrow(module_labels),
      node_diameter = "1.2 + 3 * prize / max(prize); all-zero disease uses 1.2 mm",
      node_fill = NODE_FILL, node_outline = NODE_STROKE, cross_module_edge = CROSS_COLOUR,
      edge_layout_weights = "unit weights for FR; prize/effective_cost/number not used",
      labels = "one highest-prize gene per module (rank then Entrez tiebreak) plus genes adjacent to >=2 distinct other modules",
      repel = "single joint layer; every formal node included, empty label for unlabelled nodes; seed=42; max.overlaps=Inf; max.iter=100000; max.time=10",
      label_offsets = "no manual offsets; zero offset table saved; automatic joint repel is seed-stable",
      width_mm = PLOT_WIDTH_MM, height_mm = height_mm, png_dpi = 600,
      pdf = "single-page Cairo vector PDF; Arial", font_path = font_match$path,
      metrics = layout$metrics, coord_limits = list(x = xlim, y = ylim))
  }
  input_hashes <- data.table(input = names(audited$paths), path = normalizePath(audited$paths, winslash = "/", mustWork = TRUE),
                             sha256 = vapply(audited$paths, sha256_file, character(1)))
  input_hashes <- rbind(input_hashes,
    data.table(input = c("renderer_script", "renderer_script_copy"),
      path = c(normalizePath(SCRIPT_PATH, winslash = "/", mustWork = TRUE), normalizePath(SCRIPT_COPY, winslash = "/", mustWork = TRUE)),
      sha256 = c(sha256_file(SCRIPT_PATH), sha256_file(SCRIPT_COPY))))
  file_tag <- paste0(disease, "_PCSF_network_layout_v3")
  fwrite(node_source, file.path(SOURCE_DIR, paste0(file_tag, "_nodes.csv")), bom = TRUE, na = "NA")
  fwrite(edge_source, file.path(SOURCE_DIR, paste0(file_tag, "_edges.csv")), bom = TRUE, na = "NA")
  fwrite(label_source, file.path(SOURCE_DIR, paste0(file_tag, "_labels_and_offsets.csv")), bom = TRUE, na = "NA")
  fwrite(colour_table, file.path(SOURCE_DIR, paste0(file_tag, "_module_colours.csv")), bom = TRUE, na = "NA")
  fwrite(input_hashes, file.path(SOURCE_DIR, paste0(file_tag, "_input_hashes.csv")), bom = TRUE, na = "NA")
  config_file <- file.path(CONFIG_DIR, paste0(file_tag, "_config.json"))
  jsonlite::write_json(c(layout_config, list(input_hashes = input_hashes)), config_file,
                       pretty = TRUE, auto_unbox = TRUE, na = "null", dataframe = "rows")
  log_msg("DONE disease=", disease, " nodes=", nrow(nodes), " edges=", nrow(edges),
          " modules=", length(mods), " genes_labeled=", if (nrow(nodes)) sum(node_source$label_status == "labelled") else 0L,
          " height_mm=", if (nrow(nodes)) height_mm else DEFAULT_HEIGHT_MM)
  data.table(Disease = disease, status = if (nrow(nodes)) "rendered" else "empty_candidate_set_not_assessed",
             nodes = nrow(nodes), forest_edges = nrow(edges), induced_edges_audit_only = nrow(induced),
             modules = length(mods), gene_labels = if (nrow(nodes)) sum(node_source$label_status == "labelled") else 0L,
             module_labels = if (nrow(nodes)) nrow(module_labels) else 0L,
             png = basename(png_path), pdf = basename(pdf_path), width_mm = PLOT_WIDTH_MM,
             height_mm = if (nrow(nodes)) height_mm else DEFAULT_HEIGHT_MM,
             input_audit = "passed", output_folder = FIG_DIR)
}

main <- function() {
  log_msg("RUN_START formal inputs only; solver-improvement candidates are excluded")
  summary_rows <- vector("list", length(DISEASES))
  for (i in seq_along(DISEASES)) {
    disease <- DISEASES[i]
    summary_rows[[i]] <- tryCatch(render_one(disease), error = function(err) {
      log_msg("FAILED disease=", disease, " error=", conditionMessage(err))
      trace_file <- file.path(LOG_DIR, paste0(disease, "_error_trace.txt"))
      writeLines(c(conditionMessage(err), capture.output(sys.calls())), trace_file, useBytes = TRUE)
      data.table(Disease = disease, status = "blocked_input_or_render_error", error = conditionMessage(err),
                 nodes = NA_integer_, forest_edges = NA_integer_, modules = NA_integer_,
                 gene_labels = NA_integer_, module_labels = NA_integer_)
    })
  }
  summary <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
  fwrite(summary, file.path(FIG_DIR, "network_layout_v3_render_summary.csv"), bom = TRUE, na = "NA")
  output_paths <- c(list.files(FIG_DIR, recursive = TRUE, full.names = TRUE),
                    file.path(CONFIG_DIR, list.files(CONFIG_DIR, full.names = FALSE)))
  output_paths <- unique(output_paths[file.exists(output_paths) & !dir.exists(output_paths)])
  output_paths <- output_paths[!grepl("render_layout_v3.log$", output_paths)]
  # Avoid self-reference and a stale hash for the run config written below.
  output_paths <- output_paths[!basename(output_paths) %chin%
                                 c("network_layout_v3_output_manifest.csv", "network_layout_v3_run_config.json")]
  output_manifest <- data.table(path = normalizePath(output_paths, winslash = "/", mustWork = TRUE),
                                bytes = as.numeric(file.info(output_paths)$size),
                                sha256 = vapply(output_paths, sha256_file, character(1)))
  fwrite(output_manifest, file.path(SOURCE_DIR, "network_layout_v3_output_manifest.csv"), bom = TRUE)
  run_config <- list(
    status = if (all(summary$status %in% c("rendered", "empty_candidate_set_not_assessed"))) "complete" else "partial",
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), disease_order = DISEASES,
    formal_input_root = NETWORK_DIR, forbidden_input_scope = "logs/pcst_solver_improvement_*",
    output_figure_dir = FIG_DIR, output_source_dir = SOURCE_DIR, output_config_dir = CONFIG_DIR,
    renderer = SCRIPT_COPY, renderer_sha256 = sha256_file(SCRIPT_COPY),
    package_versions = setNames(lapply(required_packages, function(p) as.character(packageVersion(p))), required_packages),
    font = list(family = "Arial", resolved_path = font_match$path),
    palette = as.list(PALETTE), node_fill = NODE_FILL, node_outline = NODE_STROKE,
    cross_module_edge = CROSS_COLOUR, disease_results = summary,
    outputs = output_manifest, no_analysis_recomputed = TRUE, drug_analysis_performed = FALSE)
  jsonlite::write_json(run_config, file.path(CONFIG_DIR, "network_layout_v3_run_config.json"),
                       pretty = TRUE, auto_unbox = TRUE, na = "null", dataframe = "rows")
  log_msg("RUN_COMPLETE status=", run_config$status, " diseases=", nrow(summary),
          " nonempty_png_pdf=", sum(summary$status == "rendered"),
          " empty_png_pdf=", sum(summary$status == "empty_candidate_set_not_assessed"))
  print(summary)
  invisible(summary)
}

main()
