#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(jsonlite);library(digest);library(ragg);library(systemfonts)})
if (!nzchar(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.65001"))) stop("无法切换中文路径locale")
ROOT <- "E:/micro_output/分析审阅中心/待审阅/Fig6_疾病vsCtrl_SHAP_PCSF模块"
RUN <- file.path(ROOT, "robustness", "20261021_disease_compact_30_75_parameter_selection_01")
SRC <- file.path(RUN, "figure_source", "stability_30_75", "S_C_network_size_30_75_source.csv")
PAL <- "D:/data/microglia_new/authority/03_reference_configs/disease_palette_contract.csv"
OUT <- file.path(RUN, "figures", "stability_30_75")
CFG <- file.path(RUN, "config", "final_closeout_30_75", "S_C_render_patch_20261022.json")
DISEASES <- c("AD","ALS","ALS_FTD","DLB","DLBD_PD","EOAD","FTD","FTD_MND","GBM","LCBM","LGG","MB","MS","PD","PSP")
PLOT_SEED <- 42L; WIDTH_MM <- 350; HEIGHT_MM <- 215; TARGET <- c(30L,75L)
sha256 <- function(p) as.character(digest::digest(file=p,algo="sha256",serialize=FALSE))
if (!all(file.exists(c(SRC,PAL))) || !dir.exists(OUT)) stop("只读source、配色或输出目录缺失")
dt <- fread(SRC,encoding="UTF-8",na.strings=c("","NA"))
req <- c("Disease","replicate","random_seed","value","metric","status","control_type","baseline_node_count","disease_color")
if(length(setdiff(req,names(dt)))) stop("S-C source缺字段")
dt[,Disease:=trimws(as.character(Disease))];dt[,replicate:=as.integer(replicate)];dt[,value:=as.numeric(value)]
if(nrow(dt)!=1500L || !setequal(unique(dt$Disease),DISEASES) || anyDuplicated(dt[,.(Disease,replicate)]) || any(dt[, .N,by=Disease]$N!=100L)) stop("S-C source不是完整15病×100次")
if(any(dt$metric!="network_real_node_count") || any(dt$status!="complete") || any(dt$control_type!="edge_cost_perturbation") || any(!is.finite(dt$value) | dt$value<0)) stop("S-C指标或重复状态异常")
pal_dt<-fread(PAL,encoding="UTF-8");if(!all(c("Disease_Abbrev","Disease_Color_Hex")%in%names(pal_dt)))stop("冻结配色缺字段")
pal_dt[,`:=`(Disease=trimws(as.character(Disease_Abbrev)),color=trimws(as.character(Disease_Color_Hex)))]
if(anyDuplicated(pal_dt$Disease)||any(!grepl("^#[0-9A-Fa-f]{6}$",pal_dt$color)))stop("冻结疾病配色无效")
palette<-setNames(pal_dt$color,pal_dt$Disease);if(!"GBM"%in%names(palette)&&"GBM_IDHwt"%in%names(palette))palette[["GBM"]]<-palette[["GBM_IDHwt"]]
if(!all(DISEASES%in%names(palette)))stop("冻结配色未覆盖当前15病")
if(any(as.character(dt$disease_color)!=unname(palette[dt$Disease])))stop("source疾病颜色与冻结配色不一致")
dt[,Disease:=factor(Disease,levels=DISEASES,ordered=TRUE)]
baseline<-unique(dt[,.(Disease,baseline_node_count=as.integer(baseline_node_count))]);if(nrow(baseline)!=15L||anyNA(baseline$baseline_node_count)||anyDuplicated(baseline$Disease))stop("未扰动基线节点数缺失或不唯一")
plot_theme<-theme_classic(base_family="Arial",base_size=9)+theme(plot.title=element_text(face="bold",size=12,colour="black",hjust=0),plot.caption=element_text(size=7.5,colour="black",hjust=0),axis.title=element_text(size=9,colour="black"),axis.text=element_text(size=8,colour="black"),axis.text.x=element_text(angle=45,hjust=1,vjust=1,colour="black"),legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=8,colour="black"),plot.margin=margin(8,10,8,8,unit="mm"))
set.seed(PLOT_SEED)
p<-ggplot(dt,aes(x=Disease,y=value,fill=Disease,colour=Disease))+annotate("rect",xmin=.5,xmax=length(DISEASES)+.5,ymin=TARGET[1],ymax=TARGET[2],fill="grey85",alpha=.45)+geom_boxplot(width=.62,coef=1.5,outlier.shape=NA,alpha=.75,linewidth=.35)+geom_point(aes(shape="Edge-cost repeat"),position=position_jitter(width=.13,height=0,seed=PLOT_SEED),alpha=.28,size=1.05)+geom_point(data=baseline,aes(x=Disease,y=baseline_node_count,shape="Unperturbed original"),inherit.aes=FALSE,colour="black",fill="black",size=2.7,stroke=.35)+scale_fill_manual(values=palette[DISEASES],drop=FALSE)+scale_colour_manual(values=palette[DISEASES],drop=FALSE)+scale_shape_manual(values=c("Edge-cost repeat"=21,"Unperturbed original"=23),breaks="Unperturbed original")+guides(fill="none",colour="none")+scale_x_discrete(drop=FALSE)+scale_y_continuous(expand=expansion(mult=c(.02,.08)))+labs(title="Network size under edge-cost perturbation",x=NULL,y="Number of real network nodes",caption="Shaded band: Baseline candidate range (30–75 nodes); it was not imposed as a constraint during perturbation solves. All observed sizes are shown. Black diamond: unperturbed original network. HD had no eligible baseline network and is not evaluated.",colour=NULL,shape=NULL)+plot_theme
if(!capabilities("cairo")||!requireNamespace("ragg",quietly=TRUE)||!requireNamespace("systemfonts",quietly=TRUE))stop("PNG/PDF绘图设备或字体工具不可用")
font<-systemfonts::match_font("Arial");if(!file.exists(font$path))stop("找不到Arial字体")
png<-file.path(OUT,"S_C_network_size_30_75.png");pdf<-file.path(OUT,"S_C_network_size_30_75.pdf")
ggsave(png,plot=p,device=ragg::agg_png,width=WIDTH_MM,height=HEIGHT_MM,units="mm",dpi=600,bg="white",limitsize=FALSE)
ggsave(pdf,plot=p,device=grDevices::cairo_pdf,width=WIDTH_MM,height=HEIGHT_MM,units="mm",family="Arial",bg="white",limitsize=FALSE)
config<-list(task="S-C-only render; no analysis recomputation",source=normalizePath(SRC,winslash="/"),source_sha256=sha256(SRC),palette=normalizePath(PAL,winslash="/"),palette_sha256=sha256(PAL),main_renderer=normalizePath(file.path(RUN,"scripts","final_closeout_30_75","fig6_30_75_final_closeout.R"),winslash="/"),main_renderer_sha256=sha256(file.path(RUN,"scripts","final_closeout_30_75","fig6_30_75_final_closeout.R")),standalone_renderer=normalizePath(sub("^--file=","",commandArgs(trailingOnly=FALSE)[grep("^--file=",commandArgs(trailingOnly=FALSE))[1]]),winslash="/"),disease_order=DISEASES,target_baseline_range=TARGET,records=nrow(dt),dpi=600,width_mm=WIDTH_MM,height_mm=HEIGHT_MM,output=list(png=list(path=normalizePath(png,winslash="/"),sha256=sha256(png)),pdf=list(path=normalizePath(pdf,winslash="/"),sha256=sha256(pdf))),note="Only current S-C rendering; 1500 frozen repeat rows used unchanged.")
write_json(config,CFG,pretty=TRUE,auto_unbox=TRUE,digits=NA)