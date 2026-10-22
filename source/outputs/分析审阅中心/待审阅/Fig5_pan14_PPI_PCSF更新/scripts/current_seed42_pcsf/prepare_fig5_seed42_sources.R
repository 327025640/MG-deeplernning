options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))
root <- "E:/micro_output/分析审阅中心/待审阅/Fig5_pan14_PPI_PCSF更新"
run <- file.path(root, "robustness/20261019_pcsf_parameter_adaptation_01/runs/b10_w3_g3/seed_42")
ass <- file.path(root, "robustness/20261019_pcsf_parameter_adaptation_01/reference_seed_assessment")
out <- file.path(root, "source_data/current_seed42_pcsf")
model_root <- "E:/micro_output/分析审阅中心/待审阅/PTG02_CD_分支疾病与调控候选/ptg14_training_v2"
base_source <- file.path(root, "source_data")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
sha256 <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) stop("Cannot hash missing file: ", path)
  normalized <- normalizePath(path, winslash="\\", mustWork=TRUE)
  out <- suppressWarnings(system2("certutil.exe",
    c("-hashfile", shQuote(normalized), "SHA256"), stdout=TRUE, stderr=TRUE))
  status <- attr(out, "status")
  hash_lines <- grepl("^[A-Fa-f0-9]{64}$", out, perl=TRUE, useBytes=TRUE)
  hashes <- tolower(out[hash_lines])
  if (!is.null(status) && status != 0L) stop("certutil SHA-256 failed for ", path, ": ", paste(out, collapse=" | "))
  if (length(hashes) != 1L) stop("Could not parse streaming SHA-256 for ", path, ": ", paste(out, collapse=" | "))
  hashes[[1L]]
}
as_bool <- function(x) {
  if (is.logical(x)) return(x)
  values <- tolower(trimws(as.character(x)))
  if (any(!is.na(values) & !values %in% c("true", "false", "1", "0"))) stop("Unexpected boolean encoding")
  ifelse(is.na(values), NA, values %in% c("true", "1"))
}
required <- c(file.path(ass,"config/assessment_config.json"), file.path(ass,"assessment_manifest.json"),
  file.path(ass,"inputs/input_sha256_manifest.csv"), file.path(ass,"summaries/seed42_nodes_with_seed_support.csv"),
  file.path(ass,"summaries/seed42_backbone_edges_with_seed_support.csv"), file.path(ass,"summaries/ten_seed_network_summary.csv"),
  file.path(ass,"summaries/all_seed_node_frequency.csv"), file.path(ass,"summaries/seed_representativeness.csv"),
  file.path(ass,"summaries/consensus_set_connectivity.csv"), file.path(run,"selected_nodes.csv"),
  file.path(run,"selected_forest_edges.tsv"), file.path(run,"run_manifest.json"),
  file.path(base_source,"all_seed_full_gene_weights.csv"), file.path(base_source,"all_gene_mapping_audit.csv"),
  file.path(root,"data/ppi_edges_ge2.tsv"), file.path(root,"data/orgHs_mapping_resource_manifest.json"),
  file.path("E:/micro_output/分析审阅中心/待审阅/Fig5_主分析与流程图收尾/source_data/PTG14_seed42_51_majority28_seed_manifest.csv"))
if (!all(file.exists(required))) stop("缺少冻结输入：", paste(required[!file.exists(required)], collapse="; "))
if (!requireNamespace("jsonlite", quietly=TRUE) || !requireNamespace("digest", quietly=TRUE) || !requireNamespace("igraph", quietly=TRUE)) stop("缺少本机已有的jsonlite/digest/igraph；不安装依赖。")
config <- jsonlite::fromJSON(required[1], simplifyVector=FALSE)
manifest <- jsonlite::fromJSON(required[2], simplifyVector=FALSE)
model_contract <- jsonlite::fromJSON(file.path(model_root,"data/contract.json"), simplifyVector=FALSE)
if (!identical(config$parameter_id, "b10_w3_g3") && !identical(as.character(config$parameter_id), "b10_w3_g3")) stop("assessment config parameter identity mismatch")
if (config$reference_seed != 42L || !identical(as.integer(unlist(config$model_seeds)), 42:51)) stop("Reference assessment model-seed identity mismatch")
if (!isTRUE(config$pcsf_rerun == FALSE) || !isTRUE(config$model_or_shap_rerun == FALSE)) stop("Assessment manifest does not establish frozen-only reuse")
if (manifest$status != "complete" || manifest$parameter_id != "b10_w3_g3" || manifest$model_identity$branch != "PTG14 pan14" || manifest$model_identity$features != 4843L) stop("Assessment identity/status mismatch")
rm <- jsonlite::fromJSON(file.path(run,"run_manifest.json"), simplifyVector=FALSE)
if (rm$status != "complete" || rm$model_seed != 42L || rm$parameter_id != "b10_w3_g3" || rm$b != 10 || rm$w != 3 || rm$g != 3 || rm$pruning != "strong" || rm$solution_metrics$real_node_count != 82L || rm$solution_metrics$real_edge_count != 81L || rm$solution_metrics$component_count_after_root_removal != 1L) stop("Seed42 selected-network flat manifest fields do not match the frozen b10/w3/g3 82-node/81-edge result")
nodes <- fread(file.path(run,"selected_nodes.csv"), na.strings=c("","NA"))
edges <- fread(file.path(run,"selected_forest_edges.tsv"), sep="\t", na.strings=c("","NA"))
support <- fread(file.path(ass,"summaries/seed42_nodes_with_seed_support.csv"), na.strings=c("","NA"))
edge_support <- fread(file.path(ass,"summaries/seed42_backbone_edges_with_seed_support.csv"), na.strings=c("","NA"))
weights <- fread(file.path(base_source,"all_seed_full_gene_weights.csv"), na.strings=c("","NA"))
mapping <- fread(file.path(base_source,"all_gene_mapping_audit.csv"), na.strings=c("","NA"))
genes <- fread(file.path(model_root,"data/genes.csv"), na.strings=c("","NA"))
if (!all(c("gene_id","PTG") %in% names(genes)) || nrow(genes)!=4843L || uniqueN(genes$gene_id)!=4843L || anyNA(genes$gene_id) || any(!nzchar(genes$gene_id))) stop("Frozen model genes.csv must contain 4,843 unique official gene_id/PTG rows")
if (!all(c("seed","gene_id","PTG","mean_abs_SHAP","weight") %in% names(weights))) stop("all_seed_full_gene_weights.csv fields do not match the frozen full-SHAP contract")
if (!all(c("seed","gene_id","PTG","mean_abs_SHAP","model_weight_sum_normalized","prize_full_model_max_normalized","mapping_status","chosen_Entrez","mapping_source","in_qualified_PPI") %in% names(mapping))) stop("all_gene_mapping_audit.csv fields do not match the frozen mapping/normalization contract")
if (uniqueN(mapping[,.(seed,gene_id)])!=nrow(mapping) || any(!mapping$gene_id %chin% genes$gene_id) || any(!mapping$seed %in% 42:51)) stop("Mapping audit contains duplicated gene-seed keys or genes outside the frozen model feature space")
if (nrow(nodes)!=82L || uniqueN(nodes$Entrez)!=82L || nrow(edges)!=81L || uniqueN(edges$edge_index)!=81L) stop("Seed42 expected 82 nodes/81 backbone edges not met")
if (sum(nodes$model_evidence %in% TRUE)!=76L || sum(nodes$role=="zero_reward_PPI_connector")!=6L) stop("Seed42 expected 76 evidence/6 connector split not met")
if (uniqueN(support$Entrez)!=82L || !setequal(as.character(nodes$Entrez),as.character(support$Entrez))) stop("Seed42 support audit node identity mismatch")
canon <- function(a,b) paste(pmin(as.numeric(a),as.numeric(b)),pmax(as.numeric(a),as.numeric(b)),sep="-")
forest_keys <- canon(edges$Entrez_A,edges$Entrez_B)
support_keys <- canon(edge_support$Entrez_A,edge_support$Entrez_B)
if (any(edges$Entrez_A==edges$Entrez_B) || uniqueN(forest_keys)!=81L || !all(c(edges$Entrez_A,edges$Entrez_B)%in%nodes$Entrez)) stop("Seed42 backbone edge integrity failure")
if (nrow(edge_support)!=81L || uniqueN(support_keys)!=81L || !setequal(forest_keys,support_keys)) {
  missing_support <- setdiff(forest_keys,support_keys); extra_support <- setdiff(support_keys,forest_keys)
  stop("Canonical edge-support identity mismatch; missing support keys=",paste(missing_support,collapse=","),"; extra support keys=",paste(extra_support,collapse=","))
}
g <- igraph::graph_from_data_frame(data.frame(from=as.character(edges$Entrez_A),to=as.character(edges$Entrez_B)),directed=FALSE,vertices=data.frame(name=as.character(nodes$Entrez)))
if (igraph::components(g)$no!=1L || any(igraph::degree(g)==0L)) stop("Seed42 PCSF backbone is not one connected component")
if (nrow(weights)!=48430L || uniqueN(weights[,.(seed,gene_id)])!=48430L || any(!weights$seed%in%42:51)) stop("Ten-seed full SHAP source is incomplete or duplicated")
if (any(!is.finite(weights$mean_abs_SHAP)) || any(weights$mean_abs_SHAP<0) || any(!is.finite(weights$weight)) || any(weights$weight<0)) stop("Full SHAP weights contain invalid values")
denominator_audit <- vector("list", length(42:51)); names(denominator_audit) <- as.character(42:51)
for (s in 42:51) {
  wm <- weights[seed==s]
  gene_order <- as.character(genes$gene_id)
  if (nrow(wm)!=4843L || uniqueN(wm$gene_id)!=4843L || !setequal(wm$gene_id,gene_order) || !identical(as.character(wm$gene_id),gene_order)) stop("Seed ",s," full-SHAP gene IDs/order differ from frozen official genes.csv")
  gene_ptg <- genes$PTG[match(wm$gene_id,genes$gene_id)]
  if (anyNA(gene_ptg) || any(as.character(wm$PTG)!=as.character(gene_ptg))) stop("Seed ",s," official Symbol/PTG identity differs from genes.csv")
  denom_sum <- sum(wm$mean_abs_SHAP); denom_max <- max(wm$mean_abs_SHAP)
  if (nrow(wm)!=4843L || any(!is.finite(wm$mean_abs_SHAP)) || any(wm$mean_abs_SHAP<0) || !is.finite(denom_sum) || denom_sum<=0 || !is.finite(denom_max) || denom_max<=0 || abs(sum(wm$weight)-1)>1e-8 || max(abs(wm$weight-wm$mean_abs_SHAP/denom_sum))>1e-12) stop("Seed ",s," full 4843-feature SHAP sum-normalization mismatch")
  ma <- mapping[seed==s][match(wm$gene_id,mapping[seed==s]$gene_id)]
  if (nrow(ma)!=4843L || anyNA(ma$gene_id) || any(as.character(ma$PTG)!=as.character(wm$PTG)) || max(abs(as.numeric(ma$mean_abs_SHAP)-wm$mean_abs_SHAP))>1e-12 || max(abs(as.numeric(ma$model_weight_sum_normalized)-wm$weight))>1e-12 || max(abs(as.numeric(ma$prize_full_model_max_normalized)-wm$mean_abs_SHAP/denom_max))>1e-12) stop("Seed ",s," mapping audit symbol/PTG/raw SHAP/full-model denominator or max-prize normalization mismatch")
  denominator_audit[[as.character(s)]] <- data.table(seed=s,n_features=nrow(wm),full_model_mean_abs_SHAP_sum=denom_sum,full_model_mean_abs_SHAP_max=denom_max)
  cp <- file.path(model_root,"runs",paste0("seed_",s))
  complete <- jsonlite::fromJSON(file.path(cp,"complete.json"),simplifyVector=FALSE)
  if (complete$contract != model_contract$id || complete$seed != s || sha256(file.path(cp,"model.pt"))!=complete$outputs$model.pt || sha256(file.path(cp,"weights.csv"))!=complete$outputs$weights.csv) stop("Frozen model/checkpoint hash identity failed for seed ",s)
  wseed <- fread(file.path(cp,"weights.csv"))
  wcomb <- wm[match(wseed$gene_id,wm$gene_id)]
  if (nrow(wseed)!=4843L || anyNA(wcomb$gene_id) || max(abs(wseed$mean_abs_SHAP-wcomb$mean_abs_SHAP))>1e-12 || max(abs(wseed$weight-wcomb$weight))>1e-12) stop("Combined ten-seed SHAP table differs from frozen weights.csv seed ",s)
}
denominator_audit <- rbindlist(denominator_audit)
fwrite(denominator_audit,file.path(out,"seed_full_SHAP_denominator_audit.csv"),na="NA")
# use authoritative Entrez-mapped seed42 network symbols, never gene symbols as node identity
compare_nodes <- merge(nodes[,.(Entrez,gene_symbols,model_evidence,raw_mean_abs_SHAP,normalized_prize,solver_prize,role)],
  support[,.(Entrez,Symbol,seed42_model_evidence,seed42_raw_mean_abs_SHAP,seed42_normalized_prize,seed42_solver_prize,seed42_node_role)],by="Entrez",sort=FALSE)
evidence_match <- compare_nodes$model_evidence == compare_nodes$seed42_model_evidence
evidence_rows <- compare_nodes$model_evidence %in% TRUE
connector_rows <- !evidence_rows
same_raw_shap <- identical(is.na(compare_nodes$raw_mean_abs_SHAP), is.na(compare_nodes$seed42_raw_mean_abs_SHAP)) &&
  all(is.finite(compare_nodes$raw_mean_abs_SHAP[evidence_rows])) &&
  all(abs(compare_nodes$raw_mean_abs_SHAP[evidence_rows]-compare_nodes$seed42_raw_mean_abs_SHAP[evidence_rows]) <= 1e-12)
same_prize <- all(abs(compare_nodes$normalized_prize-compare_nodes$seed42_normalized_prize) <= 1e-12)
same_solver_prize <- all(abs(compare_nodes$solver_prize-compare_nodes$seed42_solver_prize) <= 1e-12)
evidence_symbols_valid <- all(!is.na(compare_nodes$gene_symbols[evidence_rows]) & nzchar(compare_nodes$gene_symbols[evidence_rows]) &
  !is.na(compare_nodes$Symbol[evidence_rows]) & nzchar(compare_nodes$Symbol[evidence_rows]) &
  compare_nodes$gene_symbols[evidence_rows] == compare_nodes$Symbol[evidence_rows])
connectors_valid <- sum(connector_rows) == 6L && all(is.na(compare_nodes$raw_mean_abs_SHAP[connector_rows])) &&
  all(is.na(compare_nodes$seed42_raw_mean_abs_SHAP[connector_rows])) &&
  all(compare_nodes$normalized_prize[connector_rows] == 0) && all(compare_nodes$seed42_normalized_prize[connector_rows] == 0) &&
  all(compare_nodes$solver_prize[connector_rows] == 0) && all(compare_nodes$seed42_solver_prize[connector_rows] == 0) &&
  all(compare_nodes$role[connector_rows] == "zero_reward_PPI_connector") &&
  all(compare_nodes$seed42_node_role[connector_rows] == "zero_reward_PPI_connector")
if (nrow(compare_nodes)!=82L || !evidence_symbols_valid || !connectors_valid ||
  any(is.na(evidence_match) | !evidence_match) || !same_raw_shap || !same_prize ||
  !same_solver_prize ||
  any(as.character(compare_nodes$role)!=as.character(compare_nodes$seed42_node_role)))
  stop("Selected network nodes and support audit fields differ after Entrez-keyed evidence-aware comparison")
node_support <- copy(support)
setnames(node_support,c("Symbol","seed42_model_evidence","seed42_node_role","seed42_raw_mean_abs_SHAP","seed42_normalized_prize","seed42_solver_prize","core_10_of10"),
  c("official_symbols","model_evidence","node_role","raw_mean_abs_SHAP","normalized_prize","solver_prize","core_10_of_10"))
node_support[, model_evidence:=as_bool(model_evidence)]
evidence <- node_support[model_evidence==TRUE]
if (nrow(evidence)!=76L || anyDuplicated(evidence$official_symbols) || any(!evidence$official_symbols %in% genes$gene_id) || any(!evidence$official_symbols %in% weights[seed==42L]$gene_id)) stop("The 76 model evidence nodes do not map one-to-one to frozen official model symbols")
map42_evidence <- mapping[seed==42L & gene_id %chin% evidence$official_symbols]
if (nrow(map42_evidence)!=76L || uniqueN(map42_evidence$gene_id)!=76L || any(!map42_evidence$mapping_status %chin% c("official_unique","alias_unique")) || any(!as.logical(map42_evidence$in_qualified_PPI)) || any(as.character(map42_evidence$chosen_Entrez)!=as.character(evidence$Entrez[match(map42_evidence$gene_id,evidence$official_symbols)]))) stop("Seed42 candidate official-Symbol to Entrez/PPI mapping is not one-to-one or does not match frozen nodes")
if (sum(as.integer(support$selected_count_of_10)==10L)!=15L || sum(as.integer(support$selected_count_of_10)==10L & as_bool(support$seed42_model_evidence))!=13L) stop("Core 10/10 node identity/count differs from frozen reference audit")
ptg <- unique(weights[,.(gene_id,PTG)])
if (anyDuplicated(ptg$gene_id)) stop("PTG annotation is not unique by model gene")
weight_median <- weights[,.(median_weight=median(weight)),by=gene_id]
# seed42 display rank is descriptive ordering only, based on raw SHAP then Entrez tie-break
setorder(evidence,-raw_mean_abs_SHAP,Entrez)
evidence[, selected_rank := seq_len(.N)]
evidence[, `:=`(gene_id=official_symbols, Entrez=as.character(Entrez), selected_runs=as.integer(selected_count_of_10),
  core_10_of_10=as.integer(selected_count_of_10)==10L, seed42_mean_abs_SHAP=as.numeric(raw_mean_abs_SHAP),
  display_rank_definition="seed42 raw mean_abs_SHAP descending; Entrez ascending for ties; display order only")]
evidence <- merge(evidence,ptg,by="gene_id",all.x=TRUE,sort=FALSE)
evidence <- merge(evidence,weight_median,by="gene_id",all.x=TRUE,sort=FALSE)
if (anyNA(evidence$PTG)) stop("A model evidence gene lacks its frozen PTG annotation")
registry <- evidence[,.(gene_id,Entrez,PTG,selected_runs,core_10_of_10,selected_rank,median_weight,seed42_mean_abs_SHAP,selected_seed_list,normalized_prize,display_rank_definition)]
fwrite(registry,file.path(out,"seed42_model_evidence_76_registry.csv"),na="NA")
# every candidate retains all ten original SHAP rows even when not selected in that seed
weights_long <- merge(weights[gene_id%chin%registry$gene_id],registry[,.(gene_id,Entrez,selected_rank,selected_runs,core_10_of_10,selected_seed_list)],by="gene_id",all=TRUE,sort=FALSE)
weights_long[, selected_this_seed := mapply(function(x,s) as.character(s)%chin%strsplit(x,"\\|",fixed=FALSE)[[1]],selected_seed_list,seed)]
if (nrow(weights_long)!=760L || uniqueN(weights_long[,.(gene_id,seed)])!=760L || any(!is.finite(weights_long$weight))) stop("B source must contain 76×10 full-model observations")
setorder(weights_long,selected_rank,seed)
fwrite(weights_long,file.path(out,"seed42_B_weights_76x10.csv"),na="NA")
# full 82-node input with model-evidence and connector roles; size uses median sum-normalized weight across the ten frozen models
cn <- merge(node_support,ptg,by.x="official_symbols",by.y="gene_id",all.x=TRUE,sort=FALSE)
cn <- merge(cn,weight_median,by.x="official_symbols",by.y="gene_id",all.x=TRUE,sort=FALSE)
cn[, `:=`(Entrez=as.numeric(Entrez),gene_id=as.character(official_symbols),
  model_evaluable=as.logical(model_evidence),selected_runs=as.integer(selected_count_of_10),
  connector_status=ifelse(as.logical(model_evidence),"model_SHAP_evidence","zero_reward_PPI_connector_no_SHAP"))]
if (anyNA(cn$model_evaluable)) stop("Seed42 node evidence status contains NA; refusing to filter node rows")
cn[model_evaluable == FALSE,median_weight:=NA_real_]
if (nrow(cn)!=82L || uniqueN(cn$Entrez)!=82L || sum(cn$model_evaluable)!=76L || sum(cn$core_10_of_10)!=15L ||
  sum(cn$core_10_of_10 & cn$model_evaluable)!=13L || sum(cn$core_10_of_10 & !cn$model_evaluable)!=2L) stop("C node table identity/core-role validation failed")
fwrite(cn,file.path(out,"seed42_C_nodes_82.csv"),na="NA")
# attach all-seed support to exactly the 81 actual PCSF backbone edges, canonicalized by Entrez
es <- copy(edges)
es[, edge_key:=canon(Entrez_A,Entrez_B)]
ea <- copy(edge_support); ea[, edge_key:=canon(Entrez_A,Entrez_B)]
if (max(abs(es$effective_cost[match(es$edge_key,ea$edge_key)]-ea$effective_cost[match(es$edge_key,ea$edge_key)]))>1e-12) stop("Reference support edge costs differ from the frozen seed42 forest input")
es <- merge(es,ea[,.(edge_key,backbone_selected_count_of_10,backbone_selected_frequency,selected_seed_list,effective_cost)],by="edge_key",all.x=TRUE,sort=FALSE,suffixes=c("_forest","_support"))
if (nrow(es)!=81L || anyNA(es$backbone_selected_count_of_10) || uniqueN(es$edge_key)!=81L) stop("Backbone-edge support join failed")
nsym <- cn[,.(Entrez,gene_id,PTG,node_role,core_10_of_10)]
es <- merge(es,nsym[,.(Entrez_A=Entrez,Symbol_A=gene_id,PTG_A=PTG,role_A=node_role,core_A=core_10_of_10)],by="Entrez_A",all.x=TRUE,sort=FALSE)
es <- merge(es,nsym[,.(Entrez_B=Entrez,Symbol_B=gene_id,PTG_B=PTG,role_B=node_role,core_B=core_10_of_10)],by="Entrez_B",all.x=TRUE,sort=FALSE)
es[, `:=`(Protein_A_Entrez_ID=as.numeric(Entrez_A),Protein_B_Entrez_ID=as.numeric(Entrez_B),
 selected_in_seed_count=as.integer(backbone_selected_count_of_10),
 is_selected_any_seed=backbone_selected_count_of_10>=1L,
 display_category=ifelse(backbone_selected_count_of_10>=1L,"selected_in_at_least_1_of_10_runs","not_selected_in_other_runs"),
 display_color="#555555",display_linetype="solid")]
if (nrow(es)!=81L || any(es$Protein_A_Entrez_ID==es$Protein_B_Entrez_ID) || anyDuplicated(es[,.(Protein_A_Entrez_ID,Protein_B_Entrez_ID)])) stop("C edge source integrity failure")
fwrite(es,file.path(out,"seed42_C_backbone_edges_81.csv"),na="NA")
# stage only the frozen sources needed by the renderer and hash all source identities
source_paths <- unique(c(required,file.path(model_root,"runs",paste0("seed_",42:51),"complete.json"),file.path(model_root,"runs",paste0("seed_",42:51),"weights.csv"),file.path(model_root,"runs",paste0("seed_",42:51),"model.pt"),file.path(model_root,"data",c("contract.json","genes.csv","cells.csv","expression_cells_x_genes.f32")),file.path(base_source,c("all_seed_full_gene_weights.csv","all_gene_mapping_audit.csv")),file.path(run,c("selected_nodes.csv","selected_forest_edges.tsv","run_manifest.json"))))
hash_table <- data.table(path=source_paths,sha256=vapply(source_paths,sha256,character(1)),size_bytes=file.info(source_paths)$size)
fwrite(hash_table,file.path(out,"input_sha256_manifest.csv"))
summary <- list(status="sources_prepared",model_identity="PTG14 pan14; 4,843 features; Neurological=0/Tumor=1; Ctrl reference only",pcsf=list(b=10,w=3,g=3,pruning="strong",reference_seed=42,seed_provenance="assessment_config.reference_seed and seed42 run path; no separate pcst_fast seed field in the flat selected-network run_manifest"),seed42_network=list(nodes=nrow(nodes),edges=nrow(edges),components=1L,evidence_nodes=sum(nodes$model_evidence),zero_reward_connectors=sum(nodes$role=="zero_reward_PPI_connector")),cross_seed=list(core_nodes_10_of_10=15L,core_model_evidence=13L,core_connectors=2L,seed42_ge6_nodes=sum(as.integer(support$selected_count_of_10)>=6L)),B=list(genes=nrow(registry),observations=nrow(weights_long),weight="mean_abs_SHAP / per-seed sum across all 4,843 features; all seeds included"),C=list(nodes=nrow(cn),backbone_edges=nrow(es),size="median of per-seed sum-normalized mean_abs_SHAP; connectors minimum visible size"),input_manifest_sha256=sha256(file.path(out,"input_sha256_manifest.csv")))
jsonlite::write_json(summary,file.path(out,"source_preparation_manifest.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
cat("SOURCE_PREPARATION_OK evidence=",nrow(registry)," B_rows=",nrow(weights_long)," C_nodes=",nrow(cn)," C_edges=",nrow(es),"\n",sep="")
