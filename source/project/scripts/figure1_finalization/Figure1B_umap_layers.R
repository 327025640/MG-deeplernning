# Figure 1B（Panel B）Illustrator 分层导出：点阵层 + 矢量坐标轴层。
# 点云出无图例无文字的 PNG；坐标轴与 UMAP1/UMAP2 文字出矢量 PDF。
# 两层同画布尺寸、同面板矩形，可在 Illustrator 中直接叠准。
# 只读冻结源表，不重算统计、不改变细胞范围或配色。
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pdftools)
})

ROOT <- "E:/micro_output/分析审阅中心/待审阅/Figure1_图件收尾"
OUT <- file.path(ROOT, "figures")
INPUT <- "E:/micro_output/Figure1_Disease_V3_Disease17/source_data/panelB_umap_data.csv"
PALETTE_DISEASE <- "D:/data/microglia_new/authority/03_reference_configs/disease_palette_contract.csv"
stopifnot(file.exists(INPUT), file.exists(PALETTE_DISEASE))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# 与 Figure1_Disease_V3_Disease17_all_in_one.R 保持一致，保证点云叠绘次序不变。
RANDOM_SEED <- 20260901L
WIDTH_IN <- 8.2
HEIGHT_IN <- 6
AXIS_TITLE_SIZE <- 8
AXIS_LINE_WIDTH <- 0.4

cu <- fread(INPUT, select = c("CellID", "UMAP1", "UMAP2", "Disease"))
stopifnot(nrow(cu) == 424222L, !anyDuplicated(cu$CellID),
          !anyNA(cu$UMAP1), !anyNA(cu$UMAP2), !anyNA(cu$Disease))
actual_diseases <- sort(unique(cu$Disease))
stopifnot(length(actual_diseases) == 17L)

dpal <- fread(PALETTE_DISEASE)
stopifnot(all(c("Disease_Abbrev", "Disease_Color_Hex") %in% names(dpal)))
dis_color <- setNames(dpal$Disease_Color_Hex, dpal$Disease_Abbrev)
# 与 Figure1A 一致的三个主疾病映射，不修改共享色表。
dis_color["Ctrl"] <- "#4DAF4A"
dis_color["GBM"] <- "#B6362D"
dis_color["LGG"] <- "#8C1D18"
stopifnot(all(actual_diseases %in% names(dis_color)),
          !anyNA(dis_color[actual_diseases]))

set.seed(RANDOM_SEED)
cu_plot <- cu[sample(.N)]

# 坐标范围直接取自数据，显式写进 coord：两层不依赖 scale 训练，必然同范围。
XLIM <- range(cu$UMAP1)
YLIM <- range(cu$UMAP2)

# 两层共用的面板几何。theme_void 起手，补回轴线与轴标题；
# 刻度线与刻度数字保持 theme_void 的关闭状态（用户裁决：极简轴）。
base_theme <- theme_void(base_family = "Arial") +
  theme(
    plot.margin = margin(6, 6, 6, 6),
    axis.line = element_line(linewidth = AXIS_LINE_WIDTH, colour = "black"),
    axis.title.x = element_text(size = AXIS_TITLE_SIZE, family = "Arial",
                                colour = "black", margin = margin(t = 4)),
    axis.title.y = element_text(size = AXIS_TITLE_SIZE, family = "Arial",
                                colour = "black", angle = 90, margin = margin(r = 4))
  )

p_base <- ggplot(cu_plot, aes(UMAP1, UMAP2, colour = Disease)) +
  scale_colour_manual(values = dis_color[actual_diseases],
                      breaks = actual_diseases, name = "Disease") +
  coord_fixed(xlim = XLIM, ylim = YLIM, expand = FALSE) +
  labs(x = "UMAP1", y = "UMAP2") +
  base_theme +
  theme(legend.position = "none")

# 层1（PNG）：只有点云。轴线与轴标题改为全透明而非 element_blank()，
# 版面空间照旧占住，面板矩形才与层2一致。
# 注意必须逐个改「子元素」axis.line.x/y 与 axis.title.x/y：
# ggplot2 主题继承里子元素优先于父元素，只改父元素 axis.line/axis.title
# 不会生效（子元素已显式指定 colour），文字会照旧画出来。
p_points <- p_base +
  geom_point(size = 0.18, alpha = 0.6, stroke = 0) +
  theme(
    axis.line.x = element_line(linewidth = AXIS_LINE_WIDTH, colour = "transparent"),
    axis.line.y = element_line(linewidth = AXIS_LINE_WIDTH, colour = "transparent"),
    axis.title.x = element_text(size = AXIS_TITLE_SIZE, family = "Arial",
                                colour = "transparent", margin = margin(t = 4)),
    axis.title.y = element_text(size = AXIS_TITLE_SIZE, family = "Arial",
                                colour = "transparent", angle = 90,
                                margin = margin(r = 4))
  )

# 层2（PDF）：只有轴线与 UMAP1/UMAP2 文字。点设 colour = NA 保持数据在场，
# 坐标范围与层1完全相同，同时不落任何墨点。
p_axes <- p_base +
  geom_point(size = 0.18, alpha = 0.6, stroke = 0, colour = NA)

# 两层都必须走 cairo_pdf：agg_png 与 cairo_pdf 对同一张图的内容落点相差约 2 pt
# （字体度量不同导致版心位移），点阵走 agg_png 会与坐标轴 PDF 错开约 0.7 mm。
# 因此点阵层先出 cairo_pdf，再按 600 dpi 栅格化成 PNG，几何与矢量层同源。
pts_pdf <- file.path(OUT, ".Figure1B_umap_points_only.tmp.pdf")
ggsave(pts_pdf, p_points, width = WIDTH_IN, height = HEIGHT_IN,
       device = cairo_pdf, bg = "white")
pdftools::pdf_convert(pts_pdf, format = "png", dpi = 600,
                      filenames = file.path(OUT, "Figure1B_umap_points_only.png"),
                      verbose = FALSE)
unlink(pts_pdf)
ggsave(file.path(OUT, "Figure1B_umap_axes_only.pdf"), p_axes,
       width = WIDTH_IN, height = HEIGHT_IN, device = cairo_pdf, bg = "transparent")

cat("COMPLETE: Figure1B layers; cells=", nrow(cu),
    "; diseases=", length(actual_diseases),
    "; xlim=", paste(sprintf("%.3f", XLIM), collapse = ".."),
    "; ylim=", paste(sprintf("%.3f", YLIM), collapse = ".."), "\n", sep = "")
