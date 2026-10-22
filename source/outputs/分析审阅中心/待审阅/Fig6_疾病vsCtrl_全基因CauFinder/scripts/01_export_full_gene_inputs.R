suppressPackageStartupMessages({library(Seurat); library(Matrix); library(data.table)})
invisible(try(Sys.setlocale("LC_CTYPE", "Chinese_China.utf8"), silent = TRUE))

OBJ <- "E:/micro_output/v3_dims16/objects/merged_obj.rds"
OBJ_AUDIT <- "E:/micro_output/v3_dims16/objects/merged_obj_audit.csv"
OUT <- "E:/micro_output/分析审阅中心/待审阅/Fig6_疾病vsCtrl_全基因CauFinder"
BLOCK_DIR <- file.path(OUT, "input", "expression_csr_blocks")
BLOCK_SIZE <- 1024L
SEED <- 42L
DROP_CLUSTERS <- c("10", "15", "20", "24", "27", "28")
KEEP_GATES <- c("confirmed_microglia", "activated_microglia_candidate")
run_started <- Sys.time()
stopifnot(file.exists(OBJ), file.exists(OBJ_AUDIT), dir.exists(OUT))
dir.create(BLOCK_DIR, recursive=TRUE, showWarnings=FALSE)

get_sha256 <- function(path) {
  out <- system2("C:/WINDOWS/system32/certutil.exe", c("-hashfile", path, "SHA256"),
                 stdout=TRUE, stderr=TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) stop(paste(out, collapse="\n"))
  hex_bytes <- c(as.integer(charToRaw("0123456789abcdefABCDEF")))
  is_hash_line <- vapply(out, function(line) {
    bytes <- as.integer(charToRaw(line))
    length(bytes) == 64L && all(bytes %in% hex_bytes)
  }, logical(1))
  hits <- out[is_hash_line]
  if (length(hits) != 1L) stop("正式对象SHA256解析失败：certutil 输出中未找到唯一64位哈希行")
  hash <- rawToChar(charToRaw(hits[[1L]]))
  Encoding(hash) <- "unknown"
  hash <- tolower(hash)
  if (!grepl("^[0-9a-f]{64}$", hash)) stop("正式对象SHA256解析失败")
  hash
}
sha_before <- get_sha256(OBJ)
audit <- fread(OBJ_AUDIT)
get_audit <- function(key) {
  z <- audit[item == key, value]
  if (length(z) != 1L) stop("对象审计键缺失或重复: ", key)
  as.character(z[[1L]])
}
if (tolower(get_audit("sha256")) != sha_before) stop("对象SHA256与当前权威审计不符")

t_read <- Sys.time()
obj <- readRDS(OBJ)
read_elapsed <- as.numeric(difftime(Sys.time(), t_read, units = "secs"))
stopifnot(inherits(obj, "Seurat"), identical(Assays(obj), "RNA"), "data" %in% Layers(obj[["RNA"]]))
need <- c("BroadLineage","MicrogliaGate","Disease","DatasetID","SampleID_Global","SampleID","DonorID","seurat_clusters")
if (!all(need %in% names(obj@meta.data))) stop("meta.data缺字段: ", paste(setdiff(need,names(obj@meta.data)),collapse=","))
md <- as.data.table(obj@meta.data, keep.rownames="CellID")
genes_all <- rownames(obj[["RNA"]])
stopifnot(identical(colnames(obj),md$CellID), !anyDuplicated(md$CellID), !anyDuplicated(genes_all))
gb04 <- md[DatasetID=="LUNG_BRAIN_MET_2022" & SampleID=="GB04"]
if (nrow(gb04)==0L || !identical(sort(unique(as.character(gb04$Disease))),"GBM")) stop("当前对象GB04 Disease不唯一为GBM")
fwrite(data.table(DatasetID="LUNG_BRAIN_MET_2022",SampleID="GB04",n_cells=nrow(gb04),
                  Disease=unique(as.character(gb04$Disease)),
                  Disease_Abbrev_Analysis=unique(as.character(gb04$Disease_Abbrev_Analysis)),
                  note="仅使用当前对象Disease；遵从本次GB04=GBM指令"),
       file.path(OUT,"tables","GB04_current_Disease_audit.csv"))

cl <- as.character(md$seurat_clusters)
is_mg <- !is.na(md$BroadLineage) & md$BroadLineage=="Microglia" & as.character(md$MicrogliaGate)%in%KEEP_GATES
is_drop <- cl%in%DROP_CLUSTERS
ok <- !is.na(md$Disease)&nzchar(trimws(as.character(md$Disease)))&
      !is.na(md$SampleID_Global)&nzchar(trimws(as.character(md$SampleID_Global)))&
      !is.na(md$DatasetID)&nzchar(trimws(as.character(md$DatasetID)))
selected <- md[is_mg & !is_drop & ok, .(
  CellID, SampleID=as.character(SampleID),SampleID_Global=as.character(SampleID_Global),
  DatasetID=as.character(DatasetID),DonorID=as.character(DonorID),Disease=as.character(Disease),
  Cluster=paste0("C",sprintf("%02d",as.integer(seurat_clusters)))
)]
if (!nrow(selected)||!("Ctrl"%in%selected$Disease)) stop("clean24纳入集为空或没有Ctrl")
selected[,GlobalRow:=.I-1L]
selected[,Split:=NA_character_]

set.seed(SEED)
split_one <- function(i) {
  x <- sample(i,length(i),replace=FALSE); nt <- floor(length(x)*.64); nv <- floor(length(x)*.16)
  s <- rep("final_eval",length(x))
  if(nt>0L)s[seq_len(nt)]<-"train"
  if(nv>0L)s[nt+seq_len(nv)]<-"validation"
  data.table(row=x,Split=s)
}
split_list <- list(split_one(which(selected$Disease=="Ctrl")))
for(d in sort(setdiff(unique(selected$Disease),"Ctrl"))) split_list[[length(split_list)+1L]] <- split_one(which(selected$Disease==d))
sm <- rbindlist(split_list)
selected$Split[sm$row] <- sm$Split
if(anyNA(selected$Split)||anyDuplicated(sm$row)||nrow(sm)!=nrow(selected))stop("split未完整覆盖或有重叠")

sample_labels <- unique(selected[,.(SampleID_Global,SampleID,DatasetID,Disease,DonorID)])
map_audit <- sample_labels[,.(n_DatasetID=uniqueN(DatasetID),n_Disease=uniqueN(Disease),
  DatasetID_values=paste(sort(unique(DatasetID)),collapse=";"),Disease_values=paste(sort(unique(Disease)),collapse=";")),by=SampleID_Global]
if(any(map_audit$n_DatasetID!=1L|map_audit$n_Disease!=1L))stop("SampleID_Global映射DatasetID/Disease冲突")
coverage <- selected[,.(n_cells=.N,n_samples=uniqueN(SampleID_Global),
 n_donors=uniqueN(DonorID[!is.na(DonorID)&nzchar(trimws(DonorID))]),n_datasets=uniqueN(DatasetID)),by=Disease][order(Disease)]
split_counts <- selected[,.(n_cells=.N),by=.(Disease,Split)][order(Disease,Split)]
dataset_cov <- selected[,.(n_cells=.N,n_samples=uniqueN(SampleID_Global),
 n_donors=uniqueN(DonorID[!is.na(DonorID)&nzchar(trimws(DonorID))])),by=.(Disease,DatasetID)][order(Disease,DatasetID)]

x <- LayerData(obj,assay="RNA",layer="data")
stopifnot(inherits(x,"Matrix"),identical(rownames(x),genes_all),identical(colnames(x),colnames(obj)))
if(any(!is.finite(x@x)))stop("RNA data layer含非有限值")
n_global_train <- sum(selected$Split=="train")
sx <- numeric(length(genes_all)); sx2 <- numeric(length(genes_all))
cell_blocks <- split(seq_len(nrow(selected)),ceiling(seq_len(nrow(selected))/BLOCK_SIZE))
t_variance <- Sys.time()
for(b in seq_along(cell_blocks)){
  ii <- cell_blocks[[b]]; ii <- ii[selected$Split[ii]=="train"]; if(!length(ii))next
  z <- Matrix::t(x[,selected$CellID[ii],drop=FALSE])
  sx <- sx+Matrix::colSums(z); sx2 <- sx2+Matrix::colSums(z^2)
  if(b%%40L==0L||b==length(cell_blocks)){cat(sprintf("[variance] %d/%d\n",b,length(cell_blocks)));flush.console()}
}
mu <- sx/n_global_train; variance <- pmax(0,sx2/n_global_train-mu^2)
all_zero <- sx==0; near_zero <- variance<=1e-12; include <- !(all_zero|near_zero)
variance_elapsed <- as.numeric(difftime(Sys.time(), t_variance, units = "secs"))
gene_all <- data.table(gene_id=genes_all,RNA_row_index=seq_along(genes_all),training_mean=mu,
 training_population_variance=variance,all_zero_in_global_training=all_zero,
 near_zero_variance_le_1e_12=near_zero,include_in_model=include,
 removal_reason=fcase(all_zero,"all_zero_in_global_training",near_zero,"population_variance_le_1e-12",default="included"))
gene_model <- genes_all[include]
if(!length(gene_model))stop("零/近零方差过滤后没有模型基因")

fwrite(selected,file.path(OUT,"input","cell_manifest.csv.gz"))
fwrite(unique(selected[,.(SampleID_Global,SampleID,DatasetID,Disease,DonorID)])[order(Disease,DatasetID,SampleID_Global)],
       file.path(OUT,"input","sample_manifest.csv"))
fwrite(map_audit,file.path(OUT,"tables","sample_label_mapping_audit.csv"))
fwrite(coverage,file.path(OUT,"tables","disease_coverage.csv"))
fwrite(split_counts,file.path(OUT,"tables","disease_split_counts.csv"))
fwrite(dataset_cov,file.path(OUT,"tables","disease_dataset_coverage.csv"))
fwrite(gene_all,file.path(OUT,"input","gene_manifest_all_RNA.csv"))
fwrite(gene_all[include],file.path(OUT,"input","gene_manifest_model.csv"))
excluded <- copy(md)
excluded[,Exclusion_Reason:=fcase(!is_mg,"outside_formal_microglia_inclusion_list",is_drop,"adjudicated_contamination_cluster",
 !ok,"missing_Disease_or_SampleID_Global_or_DatasetID",default="included")]
fwrite(excluded[Exclusion_Reason!="included",.(CellID,SampleID,SampleID_Global,DatasetID,Disease,
 Cluster=as.character(seurat_clusters),BroadLineage,MicrogliaGate,Exclusion_Reason)],
 file.path(OUT,"input","excluded_cell_manifest.csv.gz"))
fwrite(excluded[Exclusion_Reason!="included",.(n_cells=.N),by=Exclusion_Reason],
 file.path(OUT,"tables","exclusion_counts.csv"))

t_export <- Sys.time()
block_manifest <- vector("list",length(cell_blocks))
for(b in seq_along(cell_blocks)){
  ii <- cell_blocks[[b]]; ids <- selected$CellID[ii]
  z <- as(Matrix::t(x[gene_model,ids,drop=FALSE]),"RsparseMatrix")
  stopifnot(nrow(z)==length(ii),ncol(z)==length(gene_model),all(is.finite(z@x)),all(z@x>=0))
  stem <- file.path(BLOCK_DIR,sprintf("block_%04d",b))
  writeBin(as.integer(z@p),paste0(stem,".indptr.i32"),size=4L,endian="little")
  writeBin(as.integer(z@j),paste0(stem,".indices.i32"),size=4L,endian="little")
  writeBin(as.numeric(z@x),paste0(stem,".data.f32"),size=4L,endian="little")
  block_manifest[[b]] <- data.table(block_id=b-1L,global_row_start=min(ii)-1L,n_rows=length(ii),
    n_genes=ncol(z),nnz=length(z@x),indptr_file=basename(paste0(stem,".indptr.i32")),
    indices_file=basename(paste0(stem,".indices.i32")),data_file=basename(paste0(stem,".data.f32")))
  rm(z)
  if(b%%20L==0L||b==length(cell_blocks)){cat(sprintf("[export] %d/%d blocks (%d cells)\n",b,length(cell_blocks),min(nrow(selected),b*BLOCK_SIZE)));flush.console()}
}
block_manifest <- rbindlist(block_manifest)
fwrite(block_manifest,file.path(OUT,"input","csr_block_manifest.csv"))
set.seed(SEED); rr<-sample(seq_len(nrow(selected)),min(64L,nrow(selected))); gg<-sample(seq_along(gene_model),min(32L,length(gene_model)))
check <- as.matrix(x[gene_model[gg],selected$CellID[rr],drop=FALSE])
ref <- CJ(cell_order=seq_along(rr),gene_order=seq_along(gg))
ref[,.(GlobalRow=selected$GlobalRow[rr[cell_order]],CellID=selected$CellID[rr[cell_order]],
       gene_id=gene_model[gg[gene_order]],RNA_data_value=as.numeric(check[cbind(gene_order,cell_order)]))] |>
  fwrite(file.path(OUT,"input","csr_export_validation_reference.csv"))

sha_after <- get_sha256(OBJ)
if(sha_before!=sha_after)stop("正式对象在输入导出期间发生变化")
audit_out <- data.table(item=c("object_path","object_sha256","object_audit_sha256","n_cells_object","n_genes_RNA",
 "n_clean24_cells","n_disease_labels","disease_labels","n_global_train_cells","n_genes_removed",
 "n_model_genes","n_csr_blocks","seed","split_rule","RNA_layer","contamination_clusters","GB04_Disease",
 "object_read_seconds","variance_screen_seconds","sparse_export_seconds","total_elapsed_seconds"),
 value=c(OBJ,sha_before,get_audit("sha256"),ncol(obj),length(genes_all),nrow(selected),uniqueN(selected$Disease),
 paste(sort(unique(selected$Disease)),collapse=";"),n_global_train,sum(!include),length(gene_model),nrow(block_manifest),
 SEED,"within Disease; Ctrl split once and reused; 64/16/20","RNA assay data",
 paste0("C",DROP_CLUSTERS,collapse=";"),"GBM",as.character(read_elapsed),
 as.character(variance_elapsed),as.character(as.numeric(difftime(Sys.time(),t_export,units="secs"))),
 as.character(as.numeric(difftime(Sys.time(),run_started,units="secs")))))
fwrite(audit_out,file.path(OUT,"tables","input_and_method_audit.csv"))
fwrite(data.table(stage=c("read_object","variance_filter","CSR_export","total"),
 elapsed_seconds=c(read_elapsed,variance_elapsed,as.numeric(difftime(Sys.time(),t_export,units="secs")),
 as.numeric(difftime(Sys.time(),run_started,units="secs")))),file.path(OUT,"tables","input_preparation_runtime.csv"))
cat(sprintf("[complete] cells=%d disease_labels=%d RNA_genes=%d model_genes=%d blocks=%d\n",
 nrow(selected),uniqueN(selected$Disease),length(genes_all),length(gene_model),nrow(block_manifest)))

