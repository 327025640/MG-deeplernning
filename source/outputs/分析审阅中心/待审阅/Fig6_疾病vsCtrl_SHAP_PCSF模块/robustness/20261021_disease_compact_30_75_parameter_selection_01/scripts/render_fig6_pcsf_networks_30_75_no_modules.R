#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(igraph); library(ggplot2); library(ggrepel)
  library(ragg); library(digest); library(jsonlite); library(grid)
})
if(!nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.65001")))) stop("R无法切换到支持中文路径的Windows locale。")
ROOT <- "E:/micro_output/分析审阅中心/待审阅/Fig6_疾病vsCtrl_SHAP_PCSF模块"
RUN_DIR <- file.path(ROOT, "robustness", "20261021_disease_compact_30_75_parameter_selection_01")
NETWORK_DIR <- file.path(RUN_DIR, "selected_networks")
REGISTRY_PATH <- file.path(RUN_DIR, "summaries", "selected_network_delivery_registry.csv")
MAPPING_PATH <- file.path(ROOT, "data", "orgHs_official_symbol_entrez.tsv")
MAPPING_MANIFEST <- file.path(ROOT, "data", "orgHs_mapping_resource_manifest.json")
INPUT_SHA_MANIFEST <- file.path(RUN_DIR, "inputs", "input_sha256_manifest.csv")
LOCK_PATH <- file.path(RUN_DIR, "config", "selected_parameters.lock.json")
SCRIPT_PATH <- "D:/data/microglia_new/analysis_scripts/render_fig6_pcsf_networks_layout_v3.R"
SCRIPT_COPY <- file.path(RUN_DIR, "scripts", "render_fig6_pcsf_networks_30_75_no_modules.R")
FIG_DIR <- file.path(RUN_DIR, "figures", "networks_no_modules")
SOURCE_DIR <- file.path(RUN_DIR, "figure_source", "networks_no_modules")
CONFIG_DIR <- file.path(RUN_DIR, "config", "networks_no_modules")
LOG_DIR <- file.path(RUN_DIR, "logs", "networks_no_modules")
DISEASE_ORDER <- c("AD", "ALS", "ALS_FTD", "DLB", "DLBD_PD", "EOAD", "FTD", "FTD_MND", "GBM", "HD", "LCBM", "LGG", "MB", "MS", "PD", "PSP")
PLOT_WIDTH_MM <- 195
DEFAULT_HEIGHT_MM <- 195
DENSE_HEIGHT_MM <- 205
NODE_FILL <- "#D9D9D9"
NODE_STROKE <- "#555555"
EDGE_COLOUR <- "#858585"
SEED <- 42L
LAYOUT_ITERATIONS <- c(4000L,8000L,12000L)
LABEL_COUNT <- 10L
sha256_file <- function(path) as.character(digest::digest(file=path,algo="sha256",serialize=FALSE))
read_csv <- function(path) data.table::fread(path,encoding="UTF-8",na.strings=c("","NA"),data.table=TRUE)
edge_key <- function(a,b) { an<-as.numeric(as.character(a)); bn<-as.numeric(as.character(b)); paste(ifelse(an<=bn,as.character(a),as.character(b)),ifelse(an<=bn,as.character(b),as.character(a)),sep="|") }
parse_bool <- function(x,label) { z<-tolower(trimws(as.character(x))); y<-z %chin% c("true","t","1","yes"); n<-z %chin% c("false","f","0","no"); if(any(!(y|n))) stop(label,"含有无法解析的逻辑值"); y }
log_msg <- function(...) { cat(paste0("[",format(Sys.time(),"%F %T %z"),"] ",paste0(...,collapse=""),"\n"),file=file.path(LOG_DIR,"render_no_modules.log"),append=TRUE) }
if (!dir.exists(ROOT) || !dir.exists(RUN_DIR)) stop("分析根目录或当前结果目录不存在。")
required_packages <- c("data.table","igraph","ggplot2","ggrepel","ragg","digest","jsonlite","systemfonts")
missing_packages <- required_packages[!vapply(required_packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing_packages)) stop("本机缺少绘图依赖，不自动安装：",paste(missing_packages,collapse=", "))
if(!capabilities("cairo")) stop("当前R不支持Cairo PDF设备。")
font_match <- systemfonts::match_font("Arial")
if(!file.exists(font_match$path)) stop("本机未解析到Arial字体文件。")
for (p in c(FIG_DIR,SOURCE_DIR,CONFIG_DIR,LOG_DIR,dirname(SCRIPT_COPY))) dir.create(p,recursive=TRUE,showWarnings=FALSE)
if(!file.exists(SCRIPT_PATH) || !file.exists(REGISTRY_PATH) || !file.exists(MAPPING_PATH) || !file.exists(MAPPING_MANIFEST) || !file.exists(INPUT_SHA_MANIFEST) || !file.exists(LOCK_PATH)) stop("绘图脚本、交付登记、参数锁、冻结映射文件或输入SHA清单缺失。")
SCRIPT_SELF <- "D:/data/microglia_new/scripts/render_fig6_pcsf_networks_30_75_no_modules.R"
if(!file.exists(SCRIPT_COPY)) {
  if(!file.copy(SCRIPT_SELF,SCRIPT_COPY,overwrite=FALSE)) stop("无法保存当前无模块脚本副本。")
} else if(!identical(tolower(sha256_file(SCRIPT_COPY)),tolower(sha256_file(SCRIPT_SELF)))) {
  archive_dir<-file.path(LOG_DIR,"script_versions"); dir.create(archive_dir,recursive=TRUE,showWarnings=FALSE)
  archived_script<-file.path(archive_dir,paste0("render_fig6_compact_networks_no_modules.",substr(sha256_file(SCRIPT_COPY),1L,12L),".R"))
  if(!file.exists(archived_script) && !file.copy(SCRIPT_COPY,archived_script,overwrite=FALSE)) stop("无法归档先前run脚本副本。")
  if(!file.copy(SCRIPT_SELF,SCRIPT_COPY,overwrite=TRUE)) stop("无法更新run脚本副本。")
}
registry <- read_csv(REGISTRY_PATH)
if(anyDuplicated(registry$Disease) || any(!registry$Disease %chin% DISEASE_ORDER)) stop("新run网络登记Disease键重复或含未知疾病。")
DISEASES <- DISEASE_ORDER[DISEASE_ORDER %chin% registry$Disease]
if(!length(DISEASES)) stop("新run没有可绘制的已选网络。")
selected_lock <- jsonlite::fromJSON(LOCK_PATH, simplifyVector=FALSE)
if(!identical(as.character(selected_lock$status), "selection_locked_for_evaluation") ||
   !identical(tolower(sha256_file(as.character(selected_lock$selected_parameters_path))),
              tolower(as.character(selected_lock$selected_parameters_sha256)))) stop("参数锁状态或selected_parameters SHA无效。")
locked_parameters <- read_csv(as.character(selected_lock$selected_parameters_path))
if(anyDuplicated(locked_parameters$Disease) || !setequal(locked_parameters$Disease, registry$Disease)) stop("选参锁疾病集合与当前网络登记不一致。")
map_manifest<-jsonlite::fromJSON(MAPPING_MANIFEST,simplifyVector=FALSE)
map_sha_registry<-read_csv(INPUT_SHA_MANIFEST)
norm_path<-function(x) tolower(gsub("\\\\","/",as.character(x)))
mapping_expected<-map_sha_registry[norm_path(path)==norm_path(MAPPING_PATH),sha256]
mapping_manifest_expected<-map_sha_registry[norm_path(path)==norm_path(MAPPING_MANIFEST),sha256]
if(length(mapping_expected)!=1L || !identical(tolower(sha256_file(MAPPING_PATH)),tolower(mapping_expected))) stop("官方Symbol映射SHA未通过冻结输入清单核对。")
if(length(mapping_manifest_expected)!=1L || !identical(tolower(sha256_file(MAPPING_MANIFEST)),tolower(mapping_manifest_expected))) stop("Symbol映射资源manifest SHA未通过冻结输入清单核对。")
map <- read_csv(MAPPING_PATH)
if(nrow(map)!=as.integer(map_manifest$official_mapping_rows) || is.null(map_manifest$version)) stop("冻结映射资源行数或版本元数据与manifest不符。")
if(!all(c("SYMBOL","ENTREZID") %chin% names(map))) stop("冻结official Symbol-Entrez映射字段不符。")
map[,ENTREZID:=as.character(ENTREZID)]; map[,SYMBOL:=as.character(SYMBOL)]
map <- unique(map[!is.na(ENTREZID)&nzchar(ENTREZID)&!is.na(SYMBOL)&nzchar(SYMBOL),.(ENTREZID,SYMBOL)])
map_by_id <- map[,.(Symbol=if(data.table::uniqueN(SYMBOL)==1L) SYMBOL[1L] else NA_character_,symbol_mapping_status=if(data.table::uniqueN(SYMBOL)==1L) "unique_frozen_official_symbol" else "ambiguous_official_symbol"),by=ENTREZID]
sha_in_manifest <- function(m,filename) { o<-m$outputs; if(is.null(o)) return(NA_character_); p<-as.character(vapply(o,function(z) as.character(z$path),character(1))); h<-as.character(vapply(o,function(z) as.character(z$sha256),character(1))); i<-which(basename(p)==filename); if(length(i)!=1L) NA_character_ else h[i] }
audit_one <- function(d) {
  reg<-registry[Disease==d]; if(nrow(reg)!=1L) stop(d,"在交付登记中不是唯一记录。")
  locked<-locked_parameters[Disease==d]; if(nrow(locked)!=1L) stop(d,"参数锁中疾病记录不是唯一。")
  nd<-as.character(reg$selected_nodes_path); ed<-as.character(reg$selected_edges_path); mf<-as.character(reg$manifest_path)
  paths<-c(nodes=nd,edges=ed,manifest=mf); if(!all(file.exists(paths))) stop(d,"冻结网络输入缺失。")
  expected_paths<-c(nodes=file.path(NETWORK_DIR,d,"selected_nodes.csv"),
                    edges=file.path(NETWORK_DIR,d,"selected_forest_edges.csv"),
                    manifest=file.path(NETWORK_DIR,d,"selected_parameter_manifest.json"))
  if(!identical(norm_path(paths),norm_path(expected_paths))) stop(d,"登记路径未指向本轮30–75 selected_networks目录。")
  h<-vapply(paths,sha256_file,character(1))
  if(!identical(tolower(unname(h["nodes"])),tolower(as.character(reg$selected_nodes_sha256))) || !identical(tolower(unname(h["edges"])),tolower(as.character(reg$selected_edges_sha256))) || !identical(tolower(unname(h["manifest"])),tolower(as.character(reg$manifest_sha256)))) stop(d,"输入SHA与新run selected_network_delivery_registry不符。")
  m<-jsonlite::fromJSON(mf,simplifyVector=FALSE)
  if(!identical(as.character(m$Disease),d) || !identical(as.character(m$status),"awaiting_user_review")) stop(d,"parameter manifest疾病或状态不匹配。")
  if(!identical(tolower(as.character(m$selected_parameters_sha256)),tolower(as.character(selected_lock$selected_parameters_sha256)))) stop(d,"疾病manifest未绑定到本轮selected_parameters锁。")
  if(!identical(as.character(m$parameter_id),as.character(locked$parameter_id)) ||
     !isTRUE(all.equal(as.numeric(m$b),as.numeric(locked$b),tolerance=0)) ||
     !isTRUE(all.equal(as.numeric(m$w),as.numeric(locked$w),tolerance=0)) ||
     !isTRUE(all.equal(as.numeric(m$g),as.numeric(locked$g),tolerance=0))) stop(d,"疾病manifest与本轮selected_parameters锁的参数不一致。")
  if(!isTRUE(all.equal(as.numeric(m$b),as.numeric(reg$b),tolerance=0)) || !isTRUE(all.equal(as.numeric(m$w),as.numeric(reg$w),tolerance=0)) || !isTRUE(all.equal(as.numeric(m$g),as.numeric(reg$g),tolerance=0))) stop(d,"parameter manifest与registry的b/w/g不一致。")
  n<-read_csv(nd); e<-read_csv(ed)
  nr<-c("Disease","parameter_id","Entrez","model_SHAP_evidence","mean_abs_SHAP","normalized_prize","solver_prize")
  er<-c("Disease","Entrez_A","Entrez_B")
  if(length(setdiff(nr,names(n)))||length(setdiff(er,names(e)))) stop(d,"节点或边表字段不全。")
  if(any(n$Disease!=d)||any(e$Disease!=d)||any(n$parameter_id!=m$parameter_id)) stop(d,"节点/边疾病或冻结manifest参数身份混杂。")
  if("parameter_id"%chin%names(reg) && any(n$parameter_id!=reg$parameter_id)) stop(d,"节点表参数身份与交付登记不一致。")
  n[,Entrez:=as.character(Entrez)]; n[,entrez_num:=suppressWarnings(as.numeric(Entrez))]
  if(anyNA(n$entrez_num)||anyDuplicated(n$Entrez)) stop(d,"Entrez缺失、非数或重复。")
  n[,model_SHAP_evidence_bool:=parse_bool(model_SHAP_evidence,paste(d,"model_SHAP_evidence"))]
  n[,mean_abs_SHAP:=as.numeric(mean_abs_SHAP)]; n[,normalized_prize:=as.numeric(normalized_prize)]; n[,solver_prize:=as.numeric(solver_prize)]
  if(any(!is.finite(n$normalized_prize)|n$normalized_prize<0)||any(!is.finite(n$solver_prize)|n$solver_prize<0)) stop(d,"奖励含非有限值或负值。")
  if(any(!is.finite(n$mean_abs_SHAP[n$model_SHAP_evidence_bool])|n$mean_abs_SHAP[n$model_SHAP_evidence_bool]<0)) stop(d,"有模型证据节点的mean_abs_SHAP无效。")
  if(any(abs(n$solver_prize-as.numeric(m$b)*n$normalized_prize)>1e-8)) stop(d,"solver_prize与b×normalized_prize不符。")
  if(any(n$normalized_prize[!n$model_SHAP_evidence_bool]!=0)) stop(d,"无模型证据连接节点不是零奖励。")
  e[,from:=as.character(Entrez_A)]; e[,to:=as.character(Entrez_B)]
  e[,from_num:=suppressWarnings(as.numeric(from))]; e[,to_num:=suppressWarnings(as.numeric(to))]
  if(anyNA(e$from_num)|anyNA(e$to_num)|any(e$from==e$to)|any(!e$from%chin%n$Entrez)|any(!e$to%chin%n$Entrez)) stop(d,"边端点无效或包含自环。")
  if(anyDuplicated(edge_key(e$from,e$to))) stop(d,"骨架边存在重复无向边。")
  if(nrow(n)!=as.integer(reg$node_count) || nrow(e)!=as.integer(reg$edge_count) || nrow(e)!=nrow(n)-1L) stop(d,"节点/边数与新run网络登记不符。")
  graph<-igraph::graph_from_data_frame(as.data.frame(e[,.(from,to)]),directed=FALSE,vertices=data.frame(name=n$Entrez))
  if(igraph::vcount(graph)!=nrow(n)||igraph::ecount(graph)!=nrow(e)||igraph::components(graph)$no!=1L) stop(d,"图结构未通过节点覆盖/单分量核对。")
  output_paths<-as.character(vapply(m$outputs,function(z) as.character(z$path),character(1)))
  output_hashes<-as.character(vapply(m$outputs,function(z) as.character(z$sha256),character(1)))
  for(fname in c("selected_nodes.csv","selected_forest_edges.csv")) {
    key<-if(fname=="selected_nodes.csv") "nodes" else "edges"
    idx<-which(norm_path(output_paths)==norm_path(paths[key]))
    if(length(idx)!=1L || !identical(tolower(output_hashes[idx]),tolower(unname(h[key])))) stop(d,"参数manifest输出路径或SHA不匹配：",fname)
  }
  n<-merge(n,map_by_id,by.x="Entrez",by.y="ENTREZID",all.x=TRUE,sort=FALSE)
  n[is.na(symbol_mapping_status),symbol_mapping_status:="not_mapped_in_frozen_official_symbol_table"]
  n[is.na(Symbol)|!nzchar(Symbol),Symbol:=Entrez]
  n[,prize_for_size:=normalized_prize]
  n[,rank_frozen:=if("rank"%chin%names(n)) suppressWarnings(as.numeric(rank)) else NA_real_]
  list(d=d,registry=reg,nodes=n,edges=e,manifest=m,paths=paths,hashes=h,graph=graph)
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

render_one <- function(d, overwrite=FALSE) {
  a<-audit_one(d); n<-copy(a$nodes); e<-copy(a$edges); m<-a$manifest
  n<-n[order(entrez_num,method="radix")]
  max_prize<-max(n$prize_for_size)
  n[,diameter_mm:=if(is.finite(max_prize)&&max_prize>0) 1.2+3*prize_for_size/max_prize else 1.2]
  e[,left_num:=pmin(from_num,to_num)]; e[,right_num:=pmax(from_num,to_num)]; e<-e[order(left_num,right_num)]
  graph<-igraph::graph_from_data_frame(as.data.frame(e[,.(from,to)]),directed=FALSE,vertices=data.frame(name=n$Entrez))
  height_mm<-if(nrow(n)>=180L) DENSE_HEIGHT_MM else DEFAULT_HEIGHT_MM
  layout<-choose_layout(graph,n,e,height_mm)
  n[,`:=`(x=layout$xy[Entrez,1L],y=layout$xy[Entrez,2L],layout_component_id=as.integer(layout$comp_id[Entrez]))]
  genes<-n[order(-prize_for_size,if("rank"%chin%names(n)) suppressWarnings(as.numeric(rank)) else entrez_num,entrez_num,na.last=TRUE)][seq_len(min(LABEL_COUNT,.N))]
  genes[,`:=`(label=Symbol,label_reason=if("rank"%chin%names(n)) "top10_normalized_prize_then_frozen_rank_then_Entrez" else "top10_normalized_prize_then_numeric_Entrez_no_rank_column",label_size_pt=7.5)]
  fit_width<-PLOT_WIDTH_MM-30; fit_height<-height_mm-62
  xspan_data<-diff(range(n$x)); yspan_data<-diff(range(n$y)); fitted_scale<-min(fit_width/max(xspan_data,1e-8),fit_height/max(yspan_data,1e-8))
  n[,`:=`(x=x*fitted_scale,y=y*fitted_scale)]
  e[,`:=`(x=n$x[match(from,n$Entrez)],y=n$y[match(from,n$Entrez)],xend=n$x[match(to,n$Entrez)],yend=n$y[match(to,n$Entrez)])]
  repel_nodes<-n[,.(Entrez,label="",label_colour="#000000",fontface="plain",font_size_pt=.5,x,y)]
  repel_nodes[match(genes$Entrez,repel_nodes$Entrez),`:=`(label=genes$label,font_size_pt=genes$label_size_pt)]
  center_x<-mean(range(n$x)); center_y<-mean(range(n$y)); n[,`:=`(x=x-center_x,y=y-center_y)]
  e[,`:=`(x=x-center_x,y=y-center_y,xend=xend-center_x,yend=yend-center_y)]
  repel_nodes[,`:=`(x=x-center_x,y=y-center_y)]
  genes[,`:=`(x=n$x[match(Entrez,n$Entrez)],y=n$y[match(Entrez,n$Entrez)])]
  repel_nodes[,`:=`(x=n$x[match(Entrez,n$Entrez)],y=n$y[match(Entrez,n$Entrez)])]
  legend_prizes<-if(max_prize>0) unique(round(stats::quantile(n$prize_for_size,c(0,.5,1),names=FALSE,na.rm=TRUE),4)) else 0
  legend_diameters<-if(max_prize>0) 1.2+3*legend_prizes/max_prize else rep(1.2,length(legend_prizes))
  xspan<-diff(range(n$x))+10; yspan<-diff(range(n$y))+10; body_h<-height_mm-58; target_aspect<-(PLOT_WIDTH_MM-18)/body_h
  xspan<-max(xspan,yspan*target_aspect); yspan<-xspan/target_aspect; xlim<-c(-xspan/2,xspan/2); ylim<-c(-yspan/2,yspan/2)
  edge_plot<-e[,.(x,y,xend,yend,edge_label="PCSF-selected backbone interaction")]
  plot<-ggplot() +
    geom_segment(data=edge_plot,aes(x=x,y=y,xend=xend,yend=yend,linetype=edge_label),colour=EDGE_COLOUR,linewidth=.25,alpha=.85,lineend="round") +
    geom_point(data=n,aes(x=x,y=y,size=diameter_mm),shape=21,stroke=.20,fill=NODE_FILL,colour=NODE_STROKE,show.legend=TRUE) +
    ggrepel::geom_text_repel(data=repel_nodes,aes(x=x,y=y,label=label,colour=label_colour,fontface=fontface,size=font_size_pt/ggplot2::.pt),family="Arial",box.padding=grid::unit(.7,"mm"),point.padding=grid::unit(.8,"mm"),max.overlaps=Inf,max.iter=100000,max.time=10,min.segment.length=0,segment.colour="#999999",segment.size=.15,seed=SEED,show.legend=FALSE,na.rm=TRUE) +
    scale_colour_identity() +
    scale_linetype_manual(name=NULL,values=c("PCSF-selected backbone interaction"="solid"),breaks="PCSF-selected backbone interaction",guide="legend") +
    scale_size_identity(name="Relative SHAP importance",breaks=legend_diameters,labels=format(legend_prizes,trim=TRUE,scientific=FALSE),guide="legend") +
    guides(linetype=guide_legend(override.aes=list(colour=EDGE_COLOUR,linewidth=.8)),size=guide_legend(override.aes=list(shape=21,fill=NODE_FILL,colour=NODE_STROKE,stroke=.2))) +
    coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off") +
    labs(title=paste0(d,"-associated candidate network"),subtitle=sprintf("%d genes | %d PCSF-selected interactions",nrow(n),nrow(e))) +
    theme_void(base_family="Arial") + theme(plot.background=element_rect(fill="white",colour=NA),panel.background=element_rect(fill="white",colour=NA),plot.title=element_text(size=11,face="bold",hjust=0,colour="black"),plot.subtitle=element_text(size=8,hjust=0,colour="black"),plot.title.position="plot",legend.position="bottom",legend.box="horizontal",legend.title=element_text(size=7),legend.text=element_text(size=7),legend.key.width=grid::unit(8,"mm"),plot.margin=margin(5,8,6,8))
  stem<-paste0(d,"_PCSF_compact_network_no_modules")
  png_path<-file.path(FIG_DIR,paste0(stem,".png")); pdf_path<-file.path(FIG_DIR,paste0(stem,".pdf"))
  node_path<-file.path(SOURCE_DIR,paste0(d,"_drawn_nodes.csv")); edge_path<-file.path(SOURCE_DIR,paste0(d,"_drawn_edges.csv")); label_path<-file.path(SOURCE_DIR,paste0(d,"_label_selection.csv")); config_path<-file.path(CONFIG_DIR,paste0(d,"_layout_config.json"))
  outs<-c(png_path,pdf_path,node_path,edge_path,label_path,config_path)
  if(any(file.exists(outs))&&!overwrite) stop(d,"已有无模块输出；为保护旧文件拒绝覆盖。")
  write_plot_safe(plot,png_path,pdf_path,PLOT_WIDTH_MM,height_mm)
  node_source<-n[,.(Disease=d,Entrez,Symbol,symbol_mapping_status,mean_abs_SHAP,normalized_prize,solver_prize,prize_for_size,model_SHAP_evidence,model_SHAP_evidence_bool,node_diameter_mm=diameter_mm,x,y,layout_component_id,labelled=Entrez%chin%genes$Entrez,label_reason=ifelse(Entrez%chin%genes$Entrez,genes$label_reason[match(Entrez,genes$Entrez)],"not_selected_top10"))]
  edge_source<-e[,.(Disease=d,Entrez_A=from,Entrez_B=to,edge_key=edge_key(from,to),backbone_identity="selected_forest_edge",draw_colour=EDGE_COLOUR,line_width_mm=.25,x,y,xend,yend)]
  label_source<-genes[,.(Disease=d,Entrez,Symbol,label,normalized_prize,label_reason,rank_rule=if("rank"%chin%names(n)) "frozen_rank_then_Entrez_tiebreak" else "numeric_Entrez_tiebreak_no_frozen_rank")]
  fwrite(node_source,node_path,bom=TRUE,na="NA"); fwrite(edge_source,edge_path,bom=TRUE,na="NA"); fwrite(label_source,label_path,bom=TRUE,na="NA")
  input_hashes<-data.table(input=c("selected_nodes","selected_forest_edges","selected_parameter_manifest","delivery_registry","official_symbol_entrez","official_mapping_manifest","input_sha256_manifest","original_renderer","no_module_renderer"),path=c(a$paths,REGISTRY_PATH,MAPPING_PATH,MAPPING_MANIFEST,INPUT_SHA_MANIFEST,SCRIPT_PATH,SCRIPT_COPY),sha256=c(a$hashes,sha256_file(REGISTRY_PATH),sha256_file(MAPPING_PATH),sha256_file(MAPPING_MANIFEST),sha256_file(INPUT_SHA_MANIFEST),sha256_file(SCRIPT_PATH),sha256_file(SCRIPT_COPY)))
  cfg<-list(Disease=d,status="rendered",show_modules=FALSE,modules_recomputed=FALSE,network_recomputed=FALSE,original_renderer_sha256=sha256_file(SCRIPT_PATH),adapted_script_sha256=sha256_file(SCRIPT_COPY),parameter_id=as.character(m$parameter_id),b=as.numeric(m$b),w=as.numeric(m$w),g=as.numeric(m$g),node_count=nrow(n),backbone_edge_count=nrow(e),layout_method="copied from approved v3: igraph FR, choose_layout geometric scoring",layout_iterations=as.list(LAYOUT_ITERATIONS),layout_seed=SEED,layout_edge_weights="unit weights; prize and costs unused",label_count=nrow(genes),label_rule=if("rank"%chin%names(n)) "normalized_prize descending, frozen rank ascending, numeric Entrez ascending" else "normalized_prize descending, numeric Entrez ascending because frozen rank absent",label_font_pt=7.5,node_size_formula="diameter_mm = 1.2 + 3 * normalized_prize / max(selected normalized_prize); solver_prize excluded",node_size_min_mm=min(n$diameter_mm),node_size_max_mm=max(n$diameter_mm),node_fill=NODE_FILL,node_outline=NODE_STROKE,edge_colour=EDGE_COLOUR,edge_width_mm=.25,figure_width_mm=PLOT_WIDTH_MM,figure_height_mm=height_mm,png_dpi=600,pdf="Cairo vector PDF, Arial",font_path=font_match$path,component_count=igraph::components(graph)$no,layout_metrics=layout$metrics,coordinate_limits=list(x=xlim,y=ylim),official_symbol_mapping=list(resource="org.Hs.eg.db",version=as.character(map_manifest$version),manifest_path=MAPPING_MANIFEST,table_path=MAPPING_PATH,sha256=sha256_file(MAPPING_PATH),unresolved_ids_labeled_by="Entrez"),input_hashes=input_hashes,outputs=list(png=list(path=png_path,sha256=sha256_file(png_path)),pdf=list(path=pdf_path,sha256=sha256_file(pdf_path)),nodes=list(path=node_path,sha256=sha256_file(node_path)),edges=list(path=edge_path,sha256=sha256_file(edge_path)),labels=list(path=label_path,sha256=sha256_file(label_path))))
  jsonlite::write_json(cfg,config_path,pretty=TRUE,auto_unbox=TRUE,na="null",dataframe="rows")
  log_msg("RENDERED ",d," nodes=",nrow(n)," edges=",nrow(e)," labels=",nrow(genes))
  data.table(Disease=d,status="rendered",nodes=nrow(n),edges=nrow(e),labels=nrow(genes),png=png_path,pdf=pdf_path,node_source=node_path,edge_source=edge_path,label_source=label_path,config=config_path)
}
args<-commandArgs(trailingOnly=TRUE)
mode<-if("--only-AD"%chin%args) "AD" else if("--all"%chin%args) "all" else stop("调用参数须为 --only-AD 或 --all。")
expected_renderer_sha<-"5e096ad09b1876c14d5c289c13c87b988ca12b1984206b30874f2d76ae5816d7"
if(!identical(tolower(sha256_file(SCRIPT_PATH)),expected_renderer_sha)) stop("原v3绘图器SHA与已确认版本不符。")
script_self<-SCRIPT_SELF
if(!file.exists(SCRIPT_COPY)) { if(!file.copy(script_self,SCRIPT_COPY,overwrite=FALSE)) stop("无法复制无模块适配脚本至run/scripts。") }
if(!identical(sha256_file(script_self),sha256_file(SCRIPT_COPY))) stop("run脚本副本与审查脚本SHA不一致。")
log_msg("START mode=",mode," renderer_v3_sha=",sha256_file(SCRIPT_PATH)," show_modules=FALSE")
targets<-if(mode=="AD") "AD" else DISEASES
rows<-list()
for(d in targets) {
  png<-file.path(FIG_DIR,paste0(d,"_PCSF_compact_network_no_modules.png")); pdf<-file.path(FIG_DIR,paste0(d,"_PCSF_compact_network_no_modules.pdf")); cfg<-file.path(CONFIG_DIR,paste0(d,"_layout_config.json"))
  if(mode=="all" && all(file.exists(c(png,pdf,cfg)))) {
    old<-jsonlite::fromJSON(cfg,simplifyVector=FALSE)
    audit_one(d)
    if(!identical(as.character(old$Disease),d) || !identical(old$show_modules,FALSE) || !identical(old$network_recomputed,FALSE) || !identical(old$modules_recomputed,FALSE)) stop(d,"现有配置身份不符，拒绝跳过。")
    if(!identical(tolower(as.character(old$adapted_script_sha256)),tolower(sha256_file(script_self))) || !identical(tolower(as.character(old$original_renderer_sha256)),tolower(sha256_file(SCRIPT_PATH)))) stop(d,"已有图件对应的绘图脚本SHA已变化，拒绝跳过。")
    source_files<-c(file.path(SOURCE_DIR,paste0(d,"_drawn_nodes.csv")),file.path(SOURCE_DIR,paste0(d,"_drawn_edges.csv")),file.path(SOURCE_DIR,paste0(d,"_label_selection.csv")))
    required_outputs<-c(png,pdf,source_files)
    if(!all(file.exists(required_outputs)) || length(old$outputs)!=5L) stop(d,"部分已有输出缺失，不能跳过或覆盖。")
    recorded_paths<-vapply(old$outputs,function(z) as.character(z$path),character(1))
    recorded_hashes<-vapply(old$outputs,function(z) as.character(z$sha256),character(1))
    if(!all(file.exists(recorded_paths))) stop(d,"配置记录的已有输出路径缺失。")
    expected_paths<-normalizePath(required_outputs,winslash="/",mustWork=TRUE)
    if(!setequal(norm_path(recorded_paths),norm_path(expected_paths))) stop(d,"已有配置记录的输出集合与本轮要求不一致。")
    actual_hashes<-vapply(recorded_paths,sha256_file,character(1))
    if(!identical(tolower(recorded_hashes),tolower(actual_hashes))) stop(d,"已有图件或source SHA与疾病配置不符。")
    rows[[d]]<-data.table(Disease=d,status="rendered_reused_verified_outputs",nodes=as.integer(old$node_count),edges=as.integer(old$backbone_edge_count),labels=as.integer(old$label_count),png=png,pdf=pdf,node_source=source_files[1],edge_source=source_files[2],label_source=source_files[3],config=cfg)
  } else rows[[d]]<-tryCatch(render_one(d),error=function(err){log_msg("FAILED ",d," : ",conditionMessage(err));data.table(Disease=d,status="failed",error=conditionMessage(err))})
}
summary<-rbindlist(rows,use.names=TRUE,fill=TRUE)
summary_path<-file.path(SOURCE_DIR,"network_no_modules_render_summary.csv"); fwrite(summary,summary_path,bom=TRUE,na="NA")
outputs<-c(summary_path,SCRIPT_PATH,SCRIPT_COPY,list.files(FIG_DIR,full.names=TRUE,recursive=TRUE),list.files(SOURCE_DIR,full.names=TRUE,recursive=TRUE),list.files(CONFIG_DIR,full.names=TRUE,recursive=TRUE))
outputs<-unique(outputs[file.exists(outputs)&!dir.exists(outputs)])
manifest<-list(status=if(mode=="all"&&nrow(summary)==length(DISEASES)&&all(summary$status%chin%c("rendered","rendered_reused_verified_outputs"))) "complete" else "partial",mode=mode,created_at=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),disease_order=DISEASES,unselected_diseases=setdiff(DISEASE_ORDER,DISEASES),show_modules=FALSE,original_renderer=SCRIPT_PATH,original_renderer_sha256=sha256_file(SCRIPT_PATH),adapted_script=SCRIPT_COPY,adapted_script_sha256=sha256_file(SCRIPT_COPY),disease_results=summary,output_manifest=data.table(path=normalizePath(outputs,winslash="/",mustWork=TRUE),sha256=vapply(outputs,sha256_file,character(1))),network_or_SHAP_recomputed=FALSE,dependencies_installed=FALSE)
manifest_path<-file.path(CONFIG_DIR,"network_no_modules_run_manifest.json"); jsonlite::write_json(manifest,manifest_path,pretty=TRUE,auto_unbox=TRUE,na="null",dataframe="rows")
log_msg("END status=",manifest$status," diseases=",nrow(summary))
print(summary)
if(any(summary$status=="failed")) quit(status=1L)

