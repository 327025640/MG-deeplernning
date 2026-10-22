# Panel F 出版级重排：UP / DOWN 上下堆叠为**一张**组合面板，直接按 Fig.1 预留框尺寸出图。
# 目标 84 x 100 mm（宽高比 0.84），矢量 PDF + 600 dpi PNG，等比例即可落入预留框，无需拉伸。
# 复用冻结 ORA 结果与资格审计表重绘；不改统计量、不改 Top3 选取规则、不改显著性与排序。
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(grid)
})

V3 <- "E:/micro_output/Figure1_Disease_V3_Disease17"
OUT <- "E:/micro_output/分析审阅中心/待审阅/Figure1_图件收尾/figures"
ORA <- file.path(V3, "source_data", "panelF_go_bp_ora_all.csv")
ELIG <- file.path(V3, "tables", "disease_eligibility_audit.csv")
PALETTE <- "D:/data/microglia_new/authority/03_reference_configs/disease_palette_contract.csv"
for (p in c(ORA, ELIG, PALETTE)) stopifnot(file.exists(p))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

WIDTH_MM <- 84
HEIGHT_MM <- 100
DPI <- 600
TOP_K <- 3L
WRAP_CHARS <- 28        # y 轴 GO term 折行宽度（作者建议 25–32）

# 字号按「不叠字、不裁切」反推，不靠试错：
# UP 16 行分到约 39 mm 面板高 → 每行约 2.43 mm；
# 两行折行标签高 = 2 × 1.2 × 字号，故字号 ≤ 2.43/(2×1.2) ≈ 2.9 pt 才不叠字。
# 缩小字号同时把左侧标签列从 34 mm 压到约 17 mm，面板宽度翻倍，
# 分面列宽随之从 2.65 mm 升到约 4.3 mm —— 分面标题才有容身之处。
Y_FONT <- 2.9
X_FONT <- 2.9
STRIP_FONT <- 2.9
LEG_TITLE_FONT <- 2.9
LEG_TEXT_FONT <- 2.6

# ---- 与大脚本一致的派生色表与疾病顺序（不修改权威文件）----
dpal <- fread(PALETTE)
pal <- copy(dpal[Disease_Abbrev != "ctrl"])
if (!"GBM" %in% pal$Disease_Abbrev) {
  gbm_row <- copy(dpal[Disease_Abbrev == "rGBM"][1])
  gbm_row[, `:=`(Disease_Abbrev = "GBM", Disease_Color_Hex = "#B6362D",
                 Display_Label = "GBM")]
  pal <- rbind(pal, gbm_row, fill = TRUE)
}
pal <- pal[!duplicated(Disease_Abbrev)]
dis_color <- setNames(pal$Disease_Color_Hex, pal$Disease_Abbrev)
dis_color["Ctrl"] <- "#4DAF4A"; dis_color["GBM"] <- "#B6362D"; dis_color["LGG"] <- "#8C1D18"
class_of_all <- setNames(pal$Disease_Class, pal$Disease_Abbrev)
class_of_all["Ctrl"] <- "Control"

elig <- fread(ELIG)
formal_diseases <- elig[formal_eligible == TRUE, Disease]
if (!"Ctrl" %in% formal_diseases) stop("Ctrl 未通过既定 Figure 1 资格门")
disease_order <- sort(setdiff(formal_diseases, "Ctrl"))

# ---- 仅供显示的 facet 标签，内部 Disease_Class 一律不改 ----
# Brain_metastasis 与 Neuroinflammatory_immune_disease 各只含 1 个疾病，
# 在 84 mm 宽下分面仅约 3.8 mm，2.9 pt 时最多容纳约 8 字符，故只能取短名。
# 完整名称放 figure legend；用户后续在 Illustrator 中自行改写。
CLASS_SHORT <- c(
  "Control" = "Control",
  "Neurodegenerative_disease" = "Neurodegenerative",
  "Neuroinflammatory_immune_disease" = "Immune",
  "Primary_brain_tumor" = "Brain tumor",
  "Brain_metastasis" = "Mets",
  "Neuropsychiatric_disease" = "Neuropsychiatric"
)

ora_all <- fread(ORA)
stopifnot(all(c("ID", "p.adjust", "Disease", "Direction") %in% names(ora_all)))

readable_term <- function(x) {
  x <- sub("^GOBP_", "", x)
  x <- gsub("_", " ", x)
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}

# ---- 出版级 theme：一切都按最终毫米尺寸设定，不走"先画大再缩小"----
theme_f <- function() {
  theme_classic(base_family = "Arial", base_size = 6) +
    theme(
      plot.margin = margin(1, 1, 1, 1, unit = "mm"),
      panel.spacing.x = unit(1, "mm"),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      axis.ticks.length = unit(0.8, "mm"),
      axis.text.y = element_text(size = Y_FONT, colour = "black", lineheight = 0.9),
      axis.text.x = element_text(size = X_FONT, angle = 45, hjust = 1, vjust = 1,
                                 colour = "black"),
      strip.background = element_rect(fill = "grey95", colour = "grey70",
                                      linewidth = 0.15),
      strip.text.x = element_text(size = STRIP_FONT, face = "bold",
                                  colour = "black", margin = margin(t = 0.4, b = 0.4)),
      plot.title = element_blank(),
      legend.title = element_text(size = LEG_TITLE_FONT),
      legend.text = element_text(size = LEG_TEXT_FONT),
      legend.key.height = unit(1.8, "mm"),
      legend.key.width = unit(1.8, "mm"),
      legend.margin = margin(0, 0, 0, 0)
    )
}

pick_terms <- function(direction) {
  selected <- copy(ora_all[Direction == direction & p.adjust < 0.05])
  setorder(selected, Disease, p.adjust)
  selected <- selected[, head(.SD, TOP_K), by = Disease]
  selected[, Term_Readable := readable_term(ID)]
  selected[, Disease_Class_Short := CLASS_SHORT[class_of_all[Disease]]]
  stopifnot(!anyNA(selected$Disease_Class_Short))

  selected[, Disease_plot := factor(Disease, levels = disease_order)]
  cls_order <- unique(na.omit(class_of_all[disease_order]))
  selected[, Disease_Class_Short := factor(Disease_Class_Short,
      levels = unname(CLASS_SHORT[cls_order]))]

  term_order <- selected[, .(n_diseases = uniqueN(Disease), min_fdr = min(p.adjust)),
                         by = Term_Readable]
  setorder(term_order, -n_diseases, min_fdr)
  selected[, Term_plot := factor(Term_Readable, levels = rev(term_order$Term_Readable))]
  selected[, NegLogFDR := -log10(p.adjust)]
  selected
}

sel_up <- pick_terms("UP")
sel_dn <- pick_terms("DOWN")

# UP 与 DOWN 必须共用同一把 size 标尺：否则两件事同时出问题——
# (1) patchwork 无法合并量程不同的图例，组合图会出现两个 -log10(FDR)；
# (2) 即便强行合并，读者也会把 DOWN 的点大小按 UP 的刻度误读。
size_max <- max(c(sel_up$NegLogFDR, sel_dn$NegLogFDR))
SIZE_LIM <- c(0, ceiling(size_max / 10) * 10)
SIZE_BRK <- pretty(c(0, SIZE_LIM[2]), n = 4)
SIZE_BRK <- SIZE_BRK[SIZE_BRK > 0 & SIZE_BRK <= SIZE_LIM[2]]

make_go <- function(selected) {
  ggplot(selected, aes(Disease_plot, Term_plot,
                       size = NegLogFDR, colour = Disease_plot)) +
    geom_point(alpha = 0.9) +
    scale_colour_manual(values = dis_color[disease_order], guide = "none") +
    scale_size_continuous(limits = SIZE_LIM, breaks = SIZE_BRK, range = c(0.6, 2.4),
                          name = expression(-log[10](FDR))) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = WRAP_CHARS)) +
    facet_grid(. ~ Disease_Class_Short, scales = "free_x", space = "free_x",
               drop = TRUE) +
    theme_f() +
    labs(title = NULL, x = NULL, y = NULL)
}

# 不设标题：UP / DOWN 字样按用户要求删除，方向由上下位置与 figure legend 说明。
p_up <- make_go(sel_up)
p_dn <- make_go(sel_dn)
cat(sprintf("共享 size 标尺: 0-%g, breaks %s\n", SIZE_LIM[2],
            paste(round(SIZE_BRK, 1), collapse = "/")))

# 两图 disease 列已严格一致（同一 disease_order + 同一 facet 结构），
# 因此上幅隐去 x 轴文字以节省纵向空间，横坐标由下幅统一给出。
p_up <- p_up + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# 单独文件：各自保留图例，便于单独查看
ggsave(file.path(OUT, "PanelF_GO_BP_UP_single.pdf"), p_up + theme(legend.position = "right"),
       width = WIDTH_MM, height = HEIGHT_MM / 2, units = "mm", device = cairo_pdf, bg = "white")
ggsave(file.path(OUT, "PanelF_GO_BP_UP_single.png"), p_up + theme(legend.position = "right"),
       width = WIDTH_MM, height = HEIGHT_MM / 2, units = "mm", dpi = DPI, bg = "white")
ggsave(file.path(OUT, "PanelF_GO_BP_DOWN_single.pdf"), p_dn + theme(legend.position = "right"),
       width = WIDTH_MM, height = HEIGHT_MM / 2, units = "mm", device = cairo_pdf, bg = "white")
ggsave(file.path(OUT, "PanelF_GO_BP_DOWN_single.png"), p_dn + theme(legend.position = "right"),
       width = WIDTH_MM, height = HEIGHT_MM / 2, units = "mm", dpi = DPI, bg = "white")

# 组合：UP 上 / DOWN 下，共享一个 -log10(FDR) 图例置于右侧
p_final <- (p_up / p_dn) +
  plot_layout(heights = c(1, 1), guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(OUT, "PanelF_GO_BP_UP_DOWN.pdf"), p_final,
       width = WIDTH_MM, height = HEIGHT_MM, units = "mm",
       device = cairo_pdf, bg = "white")
ggsave(file.path(OUT, "PanelF_GO_BP_UP_DOWN.png"), p_final,
       width = WIDTH_MM, height = HEIGHT_MM, units = "mm", dpi = DPI, bg = "white")

cat(sprintf("COMPLETE: PanelF combined %d x %d mm (ratio %.3f)\n",
            WIDTH_MM, HEIGHT_MM, WIDTH_MM / HEIGHT_MM))
