# =============================================================================
# ptg02_cd_01_composition_figures.R   （计划 §2 C 图）
# C1: 样本级区段占用比例小提琴（x = shared/arm1/arm2；组内 Ctrl/Neurological/Tumor）
# C2: 区段 100% 堆积柱（细胞级分母）
# 关键约束：
#   - 两种图分母不同，图注不得混写（C1 样本级比例；C2 细胞级组成）
#   - 三层（violin/box/point）必须使用同一 dodge 宽度，严禁虚影错位
#   - 不画超出 0–100% 的密度；不写 P 值/星号
# 输出：待审阅/PTG02_CD_分支疾病与调控候选/{figures,source_data,tables}
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
})

OUT_ROOT <- "E:/micro_output/分析审阅中心/待审阅/PTG02_CD_分支疾病与调控候选"
FIG_DIR <- file.path(OUT_ROOT, "figures")
SRC_DIR <- file.path(OUT_ROOT, "source_data")
TAB_DIR <- file.path(OUT_ROOT, "tables")
LOG_DIR <- file.path(OUT_ROOT, "logs")
for (p in c(FIG_DIR, SRC_DIR, TAB_DIR, LOG_DIR)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOG_DIR, "ptg02_cd_01_figures.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n"); cat(msg, "\n", file = log_file, append = TRUE)
}

TRAJ <- "E:/micro_output/fig3/trajectory_S4"
F_OCC <- file.path(TRAJ, "sample_level_exploration/tables/sample_segment_occupancy_all_audit.csv")
F_MAP <- file.path(OUT_ROOT, "tables/cell_segment_mapping.csv")
for (f in c(F_OCC, F_MAP)) if (!file.exists(f)) stop("Missing input: ", f)

# 固定配色（计划 §2.1）
GROUP_COLS <- c(Ctrl = "#4DAF4A", Neurological = "#4477AA", Tumor = "#B6362D")
GROUP_LEVELS <- c("Ctrl", "Neurological", "Tumor")
SEG_LEVELS <- c("trunk_shared", "arm_1", "arm_2")
SEG_LABELS <- c("shared", "arm1", "arm2")
DODGE_W <- 0.78          # 三层共用，禁止各自不同
MIN_CELLS_TRAJECTORY <- 20L

theme_cd <- function(base_size = 9) {
  theme_classic(base_family = "Arial", base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(size = base_size + 1.5, face = "bold", colour = "black"),
      plot.subtitle = element_text(size = base_size - 1, colour = "black"),
      legend.title = element_text(size = base_size - 0.5, colour = "black"),
      legend.text = element_text(size = base_size - 1, colour = "black"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}
save_pub <- function(plot, stem, width_in, height_in) {
  ggsave(paste0(stem, ".png"), plot, width = width_in, height = height_in,
         units = "in", dpi = 600, device = grDevices::png, bg = "white")
  ggsave(paste0(stem, ".pdf"), plot, width = width_in, height = height_in,
         units = "in", device = grDevices::cairo_pdf, bg = "white")
}

# =============================================================================
# C1：样本级区段占用比例（分母 = 样本在整条轨迹中的细胞数）
# =============================================================================
occ <- fread(F_OCC)
occ[, SampleID_Global := as.character(SampleID_Global)]
occ[, Disease := as.character(Disease)]
occ[, Group3 := as.character(Group3)]
if (anyDuplicated(occ$SampleID_Global)) stop("占用缓存中 SampleID_Global 有重复")

# 纳入规则：整条轨迹总细胞 >= 20（不是每个区段 >= 20，也不是表达评分过滤）
occ[, included_C1 := (eligible_trajectory %in% c(TRUE, "TRUE", 1L)) & n_trajectory >= MIN_CELLS_TRAJECTORY]
log_msg("C1 inclusion: total samples=", nrow(occ), "; included=", sum(occ$included_C1))
print_tab <- occ[, .(n_all = .N, n_included = sum(included_C1)), by = Group3]
log_msg("C1 included by Group3: ",
        paste(print_tab$Group3, print_tab$n_included, sep = "=", collapse = "; "))

d1 <- occ[included_C1 == TRUE]
long <- melt(d1[, .(SampleID_Global, Disease, DatasetID, Group3,
                    trunk_shared = fraction_trunk_shared,
                    arm_1 = fraction_arm_1,
                    arm_2 = fraction_arm_2)],
             id.vars = c("SampleID_Global", "Disease", "DatasetID", "Group3"),
             variable.name = "Segment", value.name = "Fraction")
long[, Segment := factor(as.character(Segment), levels = SEG_LEVELS, labels = SEG_LABELS)]
long[, Group3 := factor(Group3, levels = GROUP_LEVELS)]
# 前置断言：把三层 na.rm = TRUE 的静默丢行变成硬失败
stopifnot(all(is.finite(long$Fraction)), all(long$Fraction >= 0 & long$Fraction <= 1))
setorder(long, Segment, Group3, SampleID_Global)
fwrite(long, file.path(SRC_DIR, "C1_sample_segment_occupancy_source.csv"))

# 比例自洽核对
chk <- d1[, .(s = fraction_trunk_shared + fraction_arm_1 + fraction_arm_2)]
log_msg("C1 fraction row-sum max deviation: ", sprintf("%.3e", max(abs(chk$s - 1))))

p1 <- ggplot(long, aes(x = Segment, y = Fraction, fill = Group3)) +
  # trim = TRUE：密度只在观测范围内估计，不画超出 0-100% 的密度
  geom_violin(position = position_dodge(width = DODGE_W), width = DODGE_W,
              scale = "width", trim = TRUE, colour = "black", linewidth = 0.2,
              na.rm = TRUE) +
  # fill 是固定设定而非映射，会丢掉全局 aes(fill=Group3) → 必须显式重建 x×fill 分组，
  # 否则 dodge 的 n==1 直接返回原值，箱线画在区段正中且统计量是三组合并值（虚影错位）
  geom_boxplot(aes(group = interaction(Segment, Group3)),
               position = position_dodge(width = DODGE_W), width = 0.10,
               outlier.shape = NA, colour = "black", fill = "white",
               linewidth = 0.22, na.rm = TRUE) +
  geom_point(position = position_dodge(width = DODGE_W), size = 0.7,
             colour = "black", alpha = 0.6, stroke = 0, na.rm = TRUE) +
  scale_fill_manual(values = GROUP_COLS, name = NULL, breaks = GROUP_LEVELS) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = label_percent(accuracy = 1),
                     breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Sample-level occupancy of each trajectory segment",
       subtitle = paste0("Frozen trajectory subset; per sample: segment cells / all trajectory cells of that sample; ",
                         "samples with >= ", MIN_CELLS_TRAJECTORY, " trajectory cells (n = ",
                         length(unique(long$SampleID_Global)), ")"),
       x = "Trajectory segment", y = "Sample segment occupancy (%)") +
  theme_cd(base_size = 9) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 9, colour = "black"),
        plot.margin = margin(8, 10, 6, 6))
save_pub(p1, file.path(FIG_DIR, "C1_sample_segment_occupancy_violin"), 9, 5)
log_msg("C1 written")
rm(p1); invisible(gc())

# =============================================================================
# C2：区段 100% 堆积柱（分母 = 该区段全部轨迹细胞，无样本过滤）
# =============================================================================
map <- fread(F_MAP)
setnames(map, "Segment", "SegmentRaw")
map[, SegmentRaw := as.character(SegmentRaw)]
g2 <- melt(map, id.vars = c("SegmentRaw", "TOTAL"),
           measure.vars = GROUP_LEVELS, variable.name = "Group3", value.name = "n_cells")
g2[, Group3 := factor(as.character(Group3), levels = GROUP_LEVELS)]
g2[, Segment := factor(SegmentRaw, levels = SEG_LEVELS, labels = SEG_LABELS)]
g2[, Proportion := n_cells / TOTAL]
g2[, Percentage := 100 * Proportion]
setorder(g2, Segment, Group3)
fwrite(g2[, .(Segment = as.character(Segment), Group3 = as.character(Group3),
              n_cells, segment_total = TOTAL, proportion = Proportion, percentage = Percentage)],
       file.path(SRC_DIR, "C2_segment_group3_composition_source.csv"))

# 逐 Disease 计数表（不额外生成 17 疾病大图）
cells <- fread(file.path(TRAJ, "data/default/cells_with_pseudotime.csv"),
               select = c("CellID", "Disease", "SampleID_Global", "DatasetID", "State"))
cells[, Segment := factor(unname(c("trunk_shared", "arm_1", "arm_2")[as.integer(State)]),
                          levels = SEG_LEVELS, labels = SEG_LABELS)]
if (!all(cells$State %in% 1:3)) stop("cells_with_pseudotime 的 State 超出 1:3")
stopifnot(!anyNA(cells$Segment))
by_dis <- cells[, .(n_cells = .N, n_samples = uniqueN(SampleID_Global),
                    n_datasets = uniqueN(DatasetID)),
                by = .(Segment, Disease)]
seg_tot <- cells[, .(segment_total = .N), by = Segment]
# 补全完整网格：某疾病在某区段无细胞时显式写 0，而不是"没有行"
grid_dis <- CJ(Segment = factor(SEG_LABELS, levels = SEG_LABELS),
               Disease = sort(unique(cells$Disease)), unique = TRUE)
by_dis <- merge(grid_dis, by_dis, by = c("Segment", "Disease"), all.x = TRUE, sort = FALSE)
by_dis[is.na(n_cells), `:=`(n_cells = 0L, n_samples = 0L, n_datasets = 0L)]
by_dis <- merge(by_dis, seg_tot, by = "Segment", all.x = TRUE)
by_dis[, percentage_of_segment := 100 * n_cells / segment_total]
setorder(by_dis, Segment, -n_cells, Disease)
fwrite(by_dis, file.path(TAB_DIR, "C2_segment_disease_counts.csv"))

xlab_map <- setNames(
  sprintf("%s\nN = %s", as.character(seg_tot$Segment),
          format(seg_tot$segment_total, big.mark = ",", scientific = FALSE)),
  as.character(seg_tot$Segment))
log_msg("C2 segment totals: ", paste(names(xlab_map), seg_tot$segment_total, sep = "=", collapse = "; "))

lbl <- g2   # 全部段都标注百分比（最小区段 1.8%，字号 2.2pt 仍可容纳）
p2 <- ggplot(g2, aes(x = Segment, y = Proportion, fill = Group3)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.35) +
  geom_text(data = lbl, aes(label = sprintf("%.1f%%", Percentage)),
            position = position_stack(vjust = 0.5), family = "Arial",
            size = 2.2, colour = "black") +
  scale_fill_manual(values = GROUP_COLS, name = NULL, breaks = GROUP_LEVELS) +
  scale_x_discrete(labels = xlab_map) +
  scale_y_continuous(labels = label_percent(accuracy = 1), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Disease-group composition of each trajectory segment",
       subtitle = "Frozen trajectory subset; cell-level denominator: all trajectory cells in that segment (no sample filter)",
       x = "Trajectory segment", y = "Share of segment cells (%)") +
  theme_cd(base_size = 9) +
  theme(legend.position = "bottom", plot.margin = margin(8, 10, 6, 6))
save_pub(p2, file.path(FIG_DIR, "C2_segment_group3_composition"), 6, 4.5)
log_msg("C2 written")

writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "ptg02_cd_01_sessionInfo.txt"), useBytes = TRUE)
log_msg("C1/C2 DONE")
cat("\n--- C2 百分比 ---\n"); print(g2[, .(Segment, Group3, n_cells, Percentage = round(Percentage, 2))])
cat("\n--- C1 纳入样本 ---\n"); print(print_tab)
