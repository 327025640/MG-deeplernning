# 只读取冻结对象；preflight不导出矩阵，export按细胞分块导出。
args <- commandArgs(trailingOnly=TRUE)
mode <- if (length(args)) args[[1]] else "preflight"
stopifnot(mode %in% c("preflight", "export"))
pkgs <- c("SeuratObject", "Matrix", "jsonlite")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("缺包：", paste(missing, collapse=", "))
root <- "E:/micro_output/分析审阅中心/待审阅/PTG02_CD_分支疾病与调控候选"
# Windows R 4.4.3在当前C.UTF-8 locale下需显式标记中文路径，避免file.exists转换失败。
Encoding(root) <- "UTF-8"
out <- file.path(root, "ptg14_training_v2")
protected <- file.path(out, c("data/expression_cells_x_genes.f32", "data/expression_cells_x_genes.f32.partial", "data/export_complete.json", "data/contract.json"))
if (any(file.exists(protected))) stop("已有导出或锁定产物，停止以保护矩阵与行列清单；请核验后恢复下游步骤。")
dir.create(file.path(out, "data"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out, "tables"), recursive=TRUE, showWarnings=FALSE)
obj_path <- "E:/micro_output/v3_dims16/objects/merged_obj.rds"
ptg_path <- "E:/micro_output/fig3/trajectory_S4/data/default/s4h_ptg_membership.csv"
blocks_path <- file.path(root, "data/pan14_disease/P14D_features_cell_blocks.csv")
stopifnot(all(file.exists(c(obj_path, ptg_path, blocks_path))))
obj <- readRDS(obj_path)
meta <- obj[[]]
stopifnot(ncol(obj)==424222, nrow(obj)==18687,
          identical(as.character(meta$CellID), colnames(obj)),
          identical(as.character(meta$seurat_clusters), as.character(meta$RNA_snn_res.0.6)))
ptg <- read.csv(ptg_path, stringsAsFactors=FALSE)
blocks <- read.csv(blocks_path, stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(ptg$gene_id), !anyDuplicated(blocks$CellID),
          setequal(unique(ptg$PTG), sprintf("PTG%02d",1:4)),
          all(blocks$CellID %in% colnames(obj)), nrow(blocks)==365563)
md <- meta[match(blocks$CellID, colnames(obj)),,drop=FALSE]
stopifnot(identical(as.character(md$SampleID_Global), blocks$SampleID_Global))
group3 <- ifelse(md$Disease=="Ctrl", "Ctrl",
                 ifelse(md$Disease %in% c("LGG","GBM","MB","LCBM"), "Tumor", "Neurological"))
stopifnot(identical(group3, blocks$Group3))
stopifnot(all(blocks$split %in% c("train","val","holdout","ctrl_reference")))
tr <- which(blocks$split=="train")
stopifnot(length(tr)==203080, all(group3[tr] %in% c("Tumor","Neurological")))
sample_split <- split(blocks$split, blocks$SampleID_Global)
stopifnot(all(vapply(sample_split, function(z) length(unique(z))==1L, logical(1))))
x <- SeuratObject::LayerData(obj, assay="RNA", layer="data")
stopifnot(identical(colnames(x), colnames(obj)))
ptg$present <- ptg$gene_id %in% rownames(x)
present <- ptg$gene_id[ptg$present]
xx <- x[present, blocks$CellID[tr], drop=FALSE]
det <- Matrix::rowSums(xx > 0)
mu <- Matrix::rowMeans(xx)
sq <- xx; sq@x <- sq@x^2
variance <- pmax(0, Matrix::rowMeans(sq)-mu^2)
ptg$detected_train_cells <- unname(det[match(ptg$gene_id, present)])
ptg$train_variance <- unname(variance[match(ptg$gene_id, present)])
# 沿用原特征门槛：训练集至少10个细胞检出；去掉完全常量输入。
ptg$keep <- ptg$present & !is.na(ptg$detected_train_cells) &
            ptg$detected_train_cells>=10 & !is.na(ptg$train_variance) & ptg$train_variance>0
ptg$reason <- ifelse(!ptg$present,"absent_RNA",ifelse(ptg$detected_train_cells<10,"detected_lt10",ifelse(ptg$train_variance<=0,"constant_train","retained")))
write.csv(ptg,file.path(out,"tables/gene_audit.csv"),row.names=FALSE)
genes <- ptg$gene_id[ptg$keep]
rm(xx, sq); gc()
stopifnot(length(genes)>0, setequal(ptg$PTG[ptg$keep],sprintf("PTG%02d",1:4)))
cells <- data.frame(CellID=blocks$CellID, split=blocks$split, Group3=group3,
                    SampleID_Global=blocks$SampleID_Global, Disease=md$Disease,
                    DatasetID=md$DatasetID, Cluster=as.character(md$RNA_snn_res.0.6))
write.csv(cells,file.path(out,"data/cells.csv"),row.names=FALSE)
write.csv(ptg[ptg$keep,c("gene_id","PTG")],file.path(out,"data/genes.csv"),row.names=FALSE)
record <- list(mode=mode,n_cells=nrow(cells),n_genes=length(genes),n_train=length(tr),
               object=obj_path,ptg=ptg_path,blocks=blocks_path,
               source_md5=as.list(tools::md5sum(c(obj_path,ptg_path,blocks_path))),
               float_bytes=nrow(cells)*length(genes)*4,layer="RNA data",recompute=FALSE)
jsonlite::write_json(record,file.path(out,"data/preflight.json"),auto_unbox=TRUE,pretty=TRUE)
if(mode=="export") {
  target <- file.path(out,"data/expression_cells_x_genes.f32")
  if(file.exists(target)) stop("导出已存在，请先核验，不覆盖：",target)
  tmp <- paste0(target,".partial")
  if(file.exists(tmp)) stop("有未完成导出，请确认后清理：",tmp)
  con <- file(tmp,"wb")
  tryCatch({
    for(start in seq.int(1L,nrow(cells),by=1024L)) {
      ids <- cells$CellID[start:min(start+1023L,nrow(cells))]
      z <- as.matrix(x[genes,ids,drop=FALSE])
      stopifnot(all(is.finite(z)),all(z>=0))
      writeBin(as.numeric(z),con,size=4,endian="little")
    }
  },finally=close(con))
  stopifnot(file.info(tmp)$size==record$float_bytes,file.rename(tmp,target))
  record$matrix_md5 <- unname(tools::md5sum(target))
  record$cells_md5 <- unname(tools::md5sum(file.path(out,"data/cells.csv")))
  record$genes_md5 <- unname(tools::md5sum(file.path(out,"data/genes.csv")))
  jsonlite::write_json(record,file.path(out,"data/export_complete.json"),auto_unbox=TRUE,pretty=TRUE)
}
cat("完成：",mode,"；细胞",nrow(cells),"；PTG并集保留基因",length(genes),"\n")
