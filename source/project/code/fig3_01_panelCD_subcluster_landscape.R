# =====================================================================
# fig3_01_panelCD_subcluster_landscape.R
# 目的: 复刻 TAM 文献 (Cell Discovery 2026, s41421-026-00888-3) Fig 2c/2d
#       Panel C = 亚簇 x 疾病组 存在性热图 (Subcluster exist)
#       Panel D = 亚簇 x 疾病组 占比气泡图 (Subcluster percentage)
# 输入: E:/micro_output/fig3/source_data/cell_meta_merged_obj.csv
#       (由 fig3_00_extract_meta_from_merged_obj.R 生成, 来源 merged_obj.rds,
#        res0.6 = 30 簇)
#       E:/micro_output/v3_dims16/objects/merged_obj_disease_mapping.csv
#       (冻结的 19->17 疾病映射: disease_analyse -> Disease, 即 pGBM+rGBM->GBM,
#        Anaplastic_astrocytoma->LGG; 分析用列 = Disease, 共 17 组)
# 输出: E:/micro_output/fig3/
#   figures/PanelC_subcluster_exist.png + .pdf + source_data/PanelC_source.csv
#   figures/PanelD_subcluster_percentage.png + .pdf + source_data/PanelD_source.csv
#   tables/cluster_tier_map.csv
# =====================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

## ---------------- 0. 冻结参数 ----------------
IN_CSV            <- "E:/micro_output/fig3/source_data/cell_meta_merged_obj.csv"
PATH_DISEASE_MAP  <- "E:/micro_output/v3_dims16/objects/merged_obj_disease_mapping.csv"
OUT_DIR           <- "E:/micro_output/fig3"
N_TOTAL_EXPECT    <- 424222L
MIN_CELLS_PRESENT <- 20L   # "存在"判定: 簇在某疾病组 >= 20 细胞(= 项目 pseudobulk 最低门槛)
PAN_MIN           <- 8L    # Pan-disease type: 出现 >= 8 个疾病组(沿用文献绝对阈值)
SPEC_MAX          <- 2L    # Disease-specific type: 出现 <= 2 个疾病组(沿用文献)
PCT_BREAKS        <- c(1, 10, 25, 50)  # 气泡大小图例刻度(沿用文献)
PENDING_QC        <- character(0)      # res0.6 = 30 簇的 QC 尚未裁决, 暂无标记

FONT        <- "Arial"
COL_PRESENT <- "#D73027"
COL_ABSENT  <- "#F5F5F5"
COL_BLUES   <- c("<1%" = "#C6DBEF", "1-5%" = "#6BAED6",
                 "5-10%" = "#2171B5", ">10%" = "#08306B")

## ---------------- 1. 读取、疾病列覆盖与校验 ----------------
for (p in c(IN_CSV, PATH_DISEASE_MAP)) {
  if (!file.exists(p)) stop("输入不存在: ", p)
}
cells <- fread(IN_CSV)
stopifnot(all(c("CellID", "cluster", "disease_analyse") %in% names(cells)))
stopifnot(nrow(cells) == N_TOTAL_EXPECT, !any(cells$disease_analyse == ""), !anyNA(cells$cluster))

# 疾病列覆盖: disease_analyse -> Disease (17 组), 使用对象冻结映射表, 不重载大对象
dmap <- fread(PATH_DISEASE_MAP)
stopifnot(all(c("Disease_Abbrev_Analysis", "Disease", "n_cells") %in% names(dmap)))
cells[dmap, on = c(disease_analyse = "Disease_Abbrev_Analysis"), Disease := i.Disease]
stopifnot(!anyNA(cells$Disease))
dmap_agg <- dmap[, .(n_cells = sum(n_cells)), by = Disease]  # pGBM/rGBM 等多行同组需先聚合
stopifnot(setequal(unique(cells$Disease), dmap_agg$Disease))  # 双向集合校验, 防任一侧丢组
chk <- cells[, .(n = .N), by = Disease][dmap_agg, on = "Disease"]
stopifnot(chk[, all(n == n_cells)])  # 各组细胞数与冻结映射表完全一致
cat(sprintf("[1] 细胞 %d | 亚簇 %d | 疾病组 %d (Disease, 17 组)\n",
            nrow(cells), uniqueN(cells$cluster), uniqueN(cells$Disease)))

## ---------------- 2. 汇总: 簇 x 疾病 计数 / 占比 / 存在性 ----------------
all_cluster <- sort(unique(cells$cluster))
all_disease <- sort(unique(cells$Disease))
disease_order <- c("Ctrl", setdiff(all_disease, "Ctrl"))  # Ctrl 首列

grid <- as.data.table(CJ(cluster = all_cluster, Disease = disease_order))
grp <- cells[, .(n_cells = .N), by = .(cluster, Disease)]
grp <- grid[grp, on = .(cluster, Disease)]
grp[is.na(n_cells), n_cells := 0L]
tot <- cells[, .(n_tot = .N), by = Disease]
grp <- tot[grp, on = "Disease"]
grp[, pct := 100 * n_cells / n_tot]
grp[, present := as.integer(n_cells >= MIN_CELLS_PRESENT)]
stopifnot(sum(grp$n_cells) == N_TOTAL_EXPECT)

## ---------------- 3. 分层(tier)与全局行序 ----------------
cov <- grp[, .(n_dis_present = sum(present)), by = cluster]
tc  <- cells[, .(total_cells = .N), by = .(cluster)]
cov <- tc[cov, on = "cluster"]
cov[, tier := fcase(
  n_dis_present >= PAN_MIN, "Pan-disease type",
  n_dis_present <= SPEC_MAX, "Disease-specific type",
  default = "Shared type")]

tier_levels <- c("Disease-specific type", "Shared type", "Pan-disease type")
cov[, tier := factor(tier, levels = tier_levels)]
setorder(cov, tier, -n_dis_present, -total_cells)   # 自上而下; 层内覆盖度降序
final_order <- copy(cov)
K <- nrow(final_order)

y_levels <- rev(final_order$cluster)                # ggplot 第一水平在底部
grp[, cluster := factor(cluster, levels = y_levels)]
grp[, Disease := factor(Disease, levels = disease_order)]

final_order[, y_pos := K - .I + 1]
bounds <- final_order[tier != shift(tier), .(y = y_pos + 0.5)]  # 层间界线: 新层首行上缘
bounds <- bounds[!is.na(y)]
tier_mid <- final_order[, .(y_mid = mean(y_pos), y_start = min(y_pos)), by = tier]
setorder(tier_mid, y_start)
tier_mid[, y_start := NULL]

## ---------------- 4. 主题(两 panel 共用) ----------------
theme_landscape <- theme_bw(base_family = FONT, base_size = 9) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5, colour = "black"),
        axis.text.y = element_text(size = 7, colour = "black"),
        plot.title = element_text(face = "bold", size = 10, hjust = 0),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        plot.margin = margin(6, 12, 4, 4))
x_expand <- expansion(add = c(2.2, 0.3))
x_gutter <- -0.75
x_dash_l <- 0.6
x_dash_r <- length(disease_order) + 0.4

## ---------------- 5. Panel C: Subcluster exist ----------------
p_c <- ggplot(grp, aes(Disease, cluster,
                       fill = factor(present, c(1, 0), c("Present", "Absent")))) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_segment(data = bounds, aes(x = x_dash_l, xend = x_dash_r, y = y, yend = y),
               inherit.aes = FALSE, linetype = "dashed", colour = "grey30", linewidth = 0.35) +
  annotate("text", x = x_gutter, y = tier_mid$y_mid, label = tier_mid$tier,
           angle = 90, size = 8 / .pt, colour = "black") +
  scale_fill_manual(values = c(Present = COL_PRESENT, Absent = COL_ABSENT), guide = "none") +
  scale_x_discrete(expand = x_expand) +
  coord_cartesian(clip = "off") +
  labs(title = "Subcluster exist",
       subtitle = sprintf("Red = present (>= %d cells); rows grouped by disease coverage (>= %d / <= %d groups)",
                          MIN_CELLS_PRESENT, PAN_MIN, SPEC_MAX)) +
  theme_landscape
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(OUT_DIR, "figures", "PanelC_subcluster_exist.png"),
       p_c, width = 7.5, height = 8.3, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "figures", "PanelC_subcluster_exist.pdf"),
       p_c, width = 7.5, height = 8.3, device = cairo_pdf, bg = "white")

## ---------------- 6. Panel D: Subcluster percentage ----------------
pd_sub <- grp[present == 1]
pd_sub[, pct_bin := cut(pct, breaks = c(0, 1, 5, 10, Inf),
                        labels = names(COL_BLUES), right = FALSE)]
p_d <- ggplot(pd_sub, aes(Disease, cluster)) +
  geom_point(aes(size = pct, color = pct_bin)) +
  geom_segment(data = bounds, aes(x = x_dash_l, xend = x_dash_r, y = y, yend = y),
               inherit.aes = FALSE, linetype = "dashed", colour = "grey30", linewidth = 0.35) +
  annotate("text", x = x_gutter, y = tier_mid$y_mid, label = tier_mid$tier,
           angle = 90, size = 8 / .pt, colour = "black") +
  scale_size(breaks = PCT_BREAKS, range = c(0.3, 7), name = "Proportion (%)") +
  scale_color_manual(values = COL_BLUES, name = "Proportion", drop = FALSE) +
  scale_x_discrete(expand = x_expand) +
  coord_cartesian(clip = "off") +
  labs(title = "Subcluster percentage",
       subtitle = "Bubble size and color encode cluster proportion within each disease group") +
  theme_landscape +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.28, "cm"),
        legend.box = "vertical")
ggsave(file.path(OUT_DIR, "figures", "PanelD_subcluster_percentage.png"),
       p_d, width = 8.2, height = 8.3, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "figures", "PanelD_subcluster_percentage.pdf"),
       p_d, width = 8.2, height = 8.3, device = cairo_pdf, bg = "white")

## ---------------- 7. source CSV 与簇分层表 ----------------
dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
fwrite(grp[, .(cluster, Disease, n_cells, n_tot_group = n_tot, pct, present)],
       file.path(OUT_DIR, "source_data", "Fig2CD_cellcounts_and_pct.csv"))
fwrite(pd_sub[, .(cluster, Disease, pct, pct_bin)],
       file.path(OUT_DIR, "source_data", "PanelD_source.csv"))
cov[, pending_QC := cluster %in% PENDING_QC]
fwrite(cov[, .(cluster, tier, n_dis_present, total_cells, pending_QC)],
       file.path(OUT_DIR, "tables", "cluster_tier_map.csv"))

## ---------------- 8. 摘要 ----------------
cat("\n==== 分层结果 ====\n")
print(cov[, .(cluster, tier, n_dis_present, total_cells, pending_QC)])
cat("\n==== 各疾病组总细胞数 ====\n")
print(grp[cluster == all_cluster[1], .(Disease, n_tot)])
cat(sprintf("\n[done] 输出目录: %s\n", OUT_DIR))
