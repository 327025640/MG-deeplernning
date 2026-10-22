#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(digest)
  library(ragg)
  library(systemfonts)
})

if (!nzchar(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.65001"))) {
  stop("无法切换中文路径locale")
}

ROOT <- "E:/micro_output/分析审阅中心/待审阅/Fig6_疾病vsCtrl_SHAP_PCSF模块"
RUN <- file.path(ROOT, "robustness", "20261021_disease_compact_30_75_parameter_selection_01")
SCRIPT <- "D:/data/microglia_new/analysis_scripts/fig6_30_75_final_closeout.R"
PALETTE_PATH <- "D:/data/microglia_new/authority/03_reference_configs/disease_palette_contract.csv"
FIG_DIR <- file.path(RUN, "figures", "stability_30_75")
SOURCE_DIR <- file.path(RUN, "figure_source", "stability_30_75")
CONFIG_DIR <- file.path(RUN, "config", "final_closeout_30_75")
ARCHIVE_DIR <- file.path(RUN, "archive", "20261022_final_closeout_01")
LOG_DIR <- file.path(RUN, "logs", "final_closeout_30_75")
RUN_SCRIPT_DIR <- file.path(RUN, "scripts", "final_closeout_30_75")
DISEASES <- c("AD", "ALS", "ALS_FTD", "DLB", "DLBD_PD", "EOAD", "FTD", "FTD_MND",
              "GBM", "LCBM", "LGG", "MB", "MS", "PD", "PSP")
PLOT_SEED <- 42L
WIDTH_MM <- 350
HEIGHT_JACCARD_MM <- 205
HEIGHT_SIZE_MM <- 215
HEIGHT_COMPARE_MM <- 230
TARGET_SIZE <- c(30L, 75L)
SELECTION_TARGET <- file.path(RUN, "summaries", "selection_candidate_20repeats.csv")

sha256_file <- function(path) {
  as.character(digest::digest(file = path, algo = "sha256", serialize = FALSE))
}
read_csv <- function(path) {
  data.table::fread(path, encoding = "UTF-8", na.strings = c("", "NA"), data.table = TRUE)
}
assert_cols <- function(dt, cols, label) {
  missing <- setdiff(cols, names(dt))
  if (length(missing)) stop(label, "缺少字段：", paste(missing, collapse = ","))
}
write_csv_new <- function(x, path) {
  if (file.exists(path) &&
      !(identical(normalizePath(dirname(path), winslash = "/", mustWork = TRUE),
                  normalizePath(SOURCE_DIR, winslash = "/", mustWork = TRUE)) &&
        !file.exists(file.path(SOURCE_DIR, "final_closeout_artifact_hash_manifest.csv")))) {
    stop("拒绝覆盖既有文件：", path)
  }
  data.table::fwrite(x, path, na = "NA", quote = TRUE)
}
write_json_new <- function(x, path) {
  if (file.exists(path)) stop("拒绝覆盖既有文件：", path)
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA)
}
key_of <- function(dt) paste(dt$Disease, dt$parameter_id, dt$replicate, dt$repeat_record_path, sep = "\034")
num_string <- function(x) formatC(x, format = "g", digits = 17, decimal.mark = ".")

if (!dir.exists(ROOT) || !dir.exists(RUN)) stop("分析根目录或当前run不存在")
required <- c(
  file.path(RUN, "summaries", "selected_parameters.csv"),
  file.path(RUN, "summaries", "evaluation_repeat_metrics.csv"),
  file.path(RUN, "summaries", "degree_null_repeat_metrics.csv"),
  file.path(RUN, "summaries", "final_disease_summary.csv"),
  file.path(RUN, "summaries", "selection_candidate_20repeats.csv"),
  file.path(RUN, "summaries", "parameter_qualification_all_104.csv"),
  file.path(RUN, "summaries", "run_manifest.json"),
  PALETTE_PATH
)
if (!all(file.exists(required))) stop("缺少输入：", paste(required[!file.exists(required)], collapse = "; "))
pkgs <- c("data.table", "ggplot2", "jsonlite", "digest", "ragg", "systemfonts")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) || !capabilities("cairo")) stop("绘图包或Cairo PDF设备不可用；不安装依赖")
font <- systemfonts::match_font("Arial")
if (!file.exists(font$path)) stop("本机找不到Arial字体")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CONFIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ARCHIVE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RUN_SCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)
run_script_copy <- file.path(RUN_SCRIPT_DIR, basename(SCRIPT))
if (!file.exists(run_script_copy) && !file.copy(SCRIPT, run_script_copy, overwrite = FALSE)) {
  stop("无法保存当前run脚本副本")
}
if (!identical(sha256_file(SCRIPT), sha256_file(run_script_copy))) {
  stop("无法保存并核验当前run脚本副本")
}

selected <- read_csv(required[[1L]])
evaluation <- read_csv(required[[2L]])
degree_null <- read_csv(required[[3L]])
final_summary <- read_csv(required[[4L]])
qualification <- read_csv(required[[6L]])
palette_dt <- read_csv(PALETTE_PATH)
assert_cols(selected, c("Disease", "parameter_id", "b", "w", "g"), "已锁定参数表")
assert_cols(evaluation, c("Disease", "parameter_id", "replicate", "random_seed",
  "node_jaccard_to_unperturbed_adopted", "backbone_edge_jaccard_to_unperturbed_adopted",
  "network_real_node_count", "repeat_record_path", "status"), "成本扰动表")
assert_cols(degree_null, c("Disease", "parameter_id", "replicate", "random_seed",
  "node_jaccard_to_unperturbed_adopted", "repeat_record_path", "status"), "度分层奖励置换表")
assert_cols(final_summary, c("Disease", "node_count", "status"), "最终疾病汇总表")
assert_cols(qualification, c("Disease", "parameter_id", "full_model_SHAP_raw_total",
  "mapped_PPI_SHAP_raw_total"), "参数资格表")
assert_cols(palette_dt, c("Disease_Abbrev", "Disease_Color_Hex"), "疾病冻结配色")

palette_dt[, `:=`(Disease = trimws(as.character(Disease_Abbrev)),
                  Disease_Color_Hex = trimws(as.character(Disease_Color_Hex)))]
if (anyDuplicated(palette_dt$Disease) || any(!grepl("^#[0-9A-Fa-f]{6}$", palette_dt$Disease_Color_Hex))) {
  stop("冻结疾病配色存在重复疾病或非法颜色")
}
palette <- setNames(palette_dt$Disease_Color_Hex, palette_dt$Disease)
palette_aliases <- character()
if (!"GBM" %in% names(palette) && "GBM_IDHwt" %in% names(palette)) {
  palette[["GBM"]] <- palette[["GBM_IDHwt"]]
  palette_aliases <- c(GBM = "GBM uses the existing frozen GBM_IDHwt compatibility mapping")
}
if (!all(DISEASES %in% names(palette))) stop("冻结配色未覆盖全部疾病：", paste(setdiff(DISEASES, names(palette)), collapse = ","))

selected[, Disease := trimws(as.character(Disease))]
evaluation[, Disease := trimws(as.character(Disease))]
degree_null[, Disease := trimws(as.character(Disease))]
final_summary[, Disease := trimws(as.character(Disease))]
validate_repeats <- function(dt, label) {
  if (nrow(dt) != 1500L || !setequal(unique(dt$Disease), DISEASES) ||
      anyDuplicated(dt[, .(Disease, replicate)]) || any(dt[, .N, by = Disease]$N != 100L) ||
      any(dt[, .(lo = min(replicate), hi = max(replicate), n = uniqueN(replicate)), by = Disease][,
           lo != 1L | hi != 100L | n != 100L]) || any(dt$status != "complete")) {
    stop(label, "必须含15病×100个唯一完整重复")
  }
  if (any(!is.finite(as.numeric(dt$random_seed)))) stop(label, "随机种子缺失或非有限")
}
validate_repeats(evaluation, "成本扰动表")
validate_repeats(degree_null, "度分层奖励置换表")
if ("HD" %in% evaluation$Disease || "HD" %in% degree_null$Disease) stop("HD不应有本次100次评价数据")

expected_disease_order <- DISEASES
if (!setequal(unique(selected$Disease), DISEASES) || nrow(selected) != 15L || anyDuplicated(selected$Disease)) {
  stop("锁定参数表必须且只能包含15个可评价疾病")
}
if (!setequal(unique(final_summary$Disease[final_summary$status != "no_eligible_candidate"]), DISEASES)) {
  stop("最终汇总的可评价疾病与15病冻结范围不一致")
}
baseline <- final_summary[Disease %in% DISEASES, .(Disease, baseline_node_count = as.integer(node_count))]
if (nrow(baseline) != 15L || anyNA(baseline$baseline_node_count) || anyDuplicated(baseline$Disease)) {
  stop("未扰动网络基线节点数缺失或不唯一")
}

evaluation[, `:=`(
  replicate = as.integer(replicate), b = as.numeric(b), w = as.numeric(w), g = as.numeric(g),
  node_jaccard = as.numeric(node_jaccard_to_unperturbed_adopted),
  edge_jaccard = as.numeric(backbone_edge_jaccard_to_unperturbed_adopted),
  real_node_count = as.integer(network_real_node_count), control_type = "edge_cost_perturbation"
)]
degree_null[, `:=`(
  replicate = as.integer(replicate), b = as.numeric(b), w = as.numeric(w), g = as.numeric(g),
  node_jaccard = as.numeric(node_jaccard_to_unperturbed_adopted),
  control_type = "degree_stratified_reward_shuffle"
)]
if (any(!is.finite(evaluation$node_jaccard) | evaluation$node_jaccard < 0 | evaluation$node_jaccard > 1) ||
    any(!is.finite(evaluation$edge_jaccard) | evaluation$edge_jaccard < 0 | evaluation$edge_jaccard > 1) ||
    any(!is.finite(degree_null$node_jaccard) | degree_null$node_jaccard < 0 | degree_null$node_jaccard > 1) ||
    anyNA(evaluation$real_node_count) || any(evaluation$real_node_count < 0)) stop("重复指标存在NA、越界或非法规模")
if (any(evaluation$parameter_id != selected$parameter_id[match(evaluation$Disease, selected$Disease)]) ||
    any(degree_null$parameter_id != selected$parameter_id[match(degree_null$Disease, selected$Disease)])) {
  stop("重复记录的parameter_id与锁定参数不一致")
}
evaluation[, Disease := factor(Disease, levels = expected_disease_order, ordered = TRUE)]
degree_null[, Disease := factor(Disease, levels = expected_disease_order, ordered = TRUE)]
baseline[, Disease := factor(Disease, levels = expected_disease_order, ordered = TRUE)]

node_source <- evaluation[, .(Disease, parameter_id, b, w, g, replicate, random_seed,
  metric = "node_jaccard_to_unperturbed_adopted", value = node_jaccard,
  repeat_record_path, status, control_type)]
edge_source <- evaluation[, .(Disease, parameter_id, b, w, g, replicate, random_seed,
  metric = "backbone_edge_jaccard_to_unperturbed_adopted", value = edge_jaccard,
  repeat_record_path, status, control_type)]
size_source <- merge(evaluation[, .(Disease, parameter_id, b, w, g, replicate, random_seed,
  metric = "network_real_node_count", value = real_node_count, repeat_record_path, status, control_type)],
  baseline, by = "Disease", all.x = TRUE, sort = FALSE)
size_source[, `:=`(within_30_75 = value >= TARGET_SIZE[1L] & value <= TARGET_SIZE[2L],
                   target_range_lower = TARGET_SIZE[1L], target_range_upper = TARGET_SIZE[2L])]
comparison_source <- rbind(
  node_source[, .(Disease, parameter_id, b, w, g, replicate, random_seed, value,
                  repeat_record_path, status, control_type)],
  degree_null[, .(Disease, parameter_id, b, w, g, replicate, random_seed, value = node_jaccard,
                  repeat_record_path, status, control_type)], use.names = TRUE
)
comparison_source[, analysis_type := factor(
  fifelse(control_type == "edge_cost_perturbation", "Edge-cost perturbation", "Degree-stratified prize shuffle"),
  levels = c("Edge-cost perturbation", "Degree-stratified prize shuffle"), ordered = TRUE)]
if (nrow(node_source) != 1500L || nrow(edge_source) != 1500L || nrow(size_source) != 1500L ||
    nrow(comparison_source) != 3000L || anyDuplicated(node_source[, .(Disease, replicate)]) ||
    anyDuplicated(edge_source[, .(Disease, replicate)]) || anyDuplicated(comparison_source[, .(Disease, control_type, replicate)])) {
  stop("绘图source记录数或唯一性核查失败")
}
node_source[, disease_color := unname(palette[as.character(Disease)])]
edge_source[, disease_color := unname(palette[as.character(Disease)])]
size_source[, disease_color := unname(palette[as.character(Disease)])]
comparison_source[, disease_color := unname(palette[as.character(Disease)])]

plot_theme <- theme_classic(base_family = "Arial", base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 12, colour = "black", hjust = 0),
        plot.caption = element_text(size = 7.5, colour = "black", hjust = 0),
        axis.title = element_text(size = 9, colour = "black"),
        axis.text = element_text(size = 8, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
        legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 8, colour = "black"),
        plot.margin = margin(8, 10, 8, 8, unit = "mm"))
save_plot <- function(plot, stem, height_mm) {
  png <- file.path(FIG_DIR, paste0(stem, ".png"))
  pdf <- file.path(FIG_DIR, paste0(stem, ".pdf"))
  if ((file.exists(png) || file.exists(pdf)) && file.exists(file.path(SOURCE_DIR, "final_closeout_artifact_hash_manifest.csv"))) {
    stop("交付hash manifest已存在，拒绝覆盖已登记图件：", stem)
  }
  ggplot2::ggsave(png, plot = plot, device = ragg::agg_png, width = WIDTH_MM, height = height_mm,
                  units = "mm", dpi = 600, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(pdf, plot = plot, device = grDevices::cairo_pdf, width = WIDTH_MM, height = height_mm,
                  units = "mm", family = "Arial", bg = "white", limitsize = FALSE)
}

set.seed(PLOT_SEED)
p_a <- ggplot(node_source, aes(x = Disease, y = value, fill = Disease, colour = Disease)) +
  geom_boxplot(width = .62, coef = 1.5, outlier.shape = NA, alpha = .52, linewidth = .35) +
  geom_point(position = position_jitter(width = .13, height = 0, seed = PLOT_SEED), alpha = .22, size = .65) +
  scale_fill_manual(values = palette[DISEASES], drop = FALSE) +
  scale_colour_manual(values = palette[DISEASES], drop = FALSE) +
  guides(fill = "none", colour = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), expand = expansion(mult = c(.01, .02))) +
  labs(title = "Node-set stability under edge-cost perturbation", x = NULL,
       y = "Jaccard similarity to unperturbed network",
       caption = "Boxes show Q1–Q3 and median; whiskers extend to the most extreme observation within 1.5×IQR. All 100 cost-perturbation repeats per eligible disease are shown. HD had no eligible baseline network and is not evaluated.") + plot_theme
save_plot(p_a, "S_A_node_Jaccard_30_75", HEIGHT_JACCARD_MM)

p_b <- ggplot(edge_source, aes(x = Disease, y = value, fill = Disease, colour = Disease)) +
  geom_boxplot(width = .62, coef = 1.5, outlier.shape = NA, alpha = .52, linewidth = .35) +
  geom_point(position = position_jitter(width = .13, height = 0, seed = PLOT_SEED), alpha = .22, size = .65) +
  scale_fill_manual(values = palette[DISEASES], drop = FALSE) +
  scale_colour_manual(values = palette[DISEASES], drop = FALSE) +
  guides(fill = "none", colour = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), expand = expansion(mult = c(.01, .02))) +
  labs(title = "Backbone-edge stability under edge-cost perturbation", x = NULL,
       y = "Jaccard similarity to unperturbed backbone",
       caption = "Boxes show Q1–Q3 and median; whiskers extend to the most extreme observation within 1.5×IQR. All 100 cost-perturbation repeats per eligible disease are shown. HD had no eligible baseline network and is not evaluated.") + plot_theme
save_plot(p_b, "S_B_backbone_edge_Jaccard_30_75", HEIGHT_JACCARD_MM)

baseline_points <- baseline[, .(Disease, value = baseline_node_count)]
p_c <- ggplot(size_source, aes(x = Disease, y = value, fill = Disease, colour = Disease)) +
  annotate("rect", xmin = .5, xmax = length(DISEASES) + .5, ymin = TARGET_SIZE[1L], ymax = TARGET_SIZE[2L],
           fill = "grey85", alpha = .45) +
  geom_boxplot(width = .62, coef = 1.5, outlier.shape = NA, alpha = .75, linewidth = .35) +
  geom_point(aes(shape = "Edge-cost repeat"),
             position = position_jitter(width = .13, height = 0, seed = PLOT_SEED), alpha = .28, size = 1.05) +
  geom_point(data = baseline_points, aes(x = Disease, y = value, shape = "Unperturbed original"),
             inherit.aes = FALSE, colour = "black", fill = "black", size = 2.7, stroke = .35) +
  scale_fill_manual(values = palette[DISEASES], drop = FALSE) +
  scale_shape_manual(values = c("Edge-cost repeat" = 21, "Unperturbed original" = 23), breaks = "Unperturbed original") +
  guides(fill = "none", colour = "none") + scale_x_discrete(drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(.02, .08))) +
  annotate("text", x = .62, y = TARGET_SIZE[2L], label = "Baseline candidate range", hjust = 0, vjust = -0.25,
           family = "Arial", size = 3, colour = "#444444") +
  labs(title = "Network size under edge-cost perturbation", x = NULL, y = "Number of real network nodes",
       caption = "Grey band marks the 30–75-node baseline candidate-selection range; it was not imposed as a constraint during perturbation solves. All observed sizes are shown. Black diamond: unperturbed original network. HD had no eligible baseline network and is not evaluated.",
       colour = NULL, shape = NULL) + plot_theme
save_plot(p_c, "S_C_network_size_30_75", HEIGHT_SIZE_MM)

p_d <- ggplot(comparison_source, aes(x = Disease, y = value, fill = analysis_type, colour = analysis_type)) +
  geom_boxplot(position = position_dodge(width = .78), width = .68, coef = 1.5,
               outlier.shape = NA, alpha = .50, linewidth = .30) +
  geom_point(position = position_jitterdodge(jitter.width = .12, jitter.height = 0, dodge.width = .78, seed = PLOT_SEED),
             alpha = .19, size = .58) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), expand = expansion(mult = c(.01, .02))) +
  scale_fill_manual(values = c("Edge-cost perturbation" = "#4E79A7", "Degree-stratified prize shuffle" = "#E15759")) +
  scale_colour_manual(values = c("Edge-cost perturbation" = "#315A7D", "Degree-stratified prize shuffle" = "#A83A3B")) +
  labs(title = "Network overlap under cost perturbation and prize shuffling", x = NULL,
       y = "Node Jaccard similarity to unperturbed network",
       caption = "Both mechanisms compare with the same disease's unperturbed network. Cost perturbation and degree-stratified prize shuffle are distinct descriptive contrasts, not biological causal validation. All 100 repeats from each mechanism are shown; HD was not evaluated.",
       fill = NULL, colour = NULL) + plot_theme
save_plot(p_d, "S_D_cost_perturbation_vs_reward_null_30_75", HEIGHT_COMPARE_MM)

excluded <- data.table(Disease = "HD", status = "not_evaluated", reason = "no eligible baseline network; no 100-repeat distributions")
source_outputs <- list(
  node = list(data = node_source, filename = "S_A_node_Jaccard_30_75_source.csv"),
  edge = list(data = edge_source, filename = "S_B_backbone_edge_Jaccard_30_75_source.csv"),
  size = list(data = size_source, filename = "S_C_network_size_30_75_source.csv"),
  comparison = list(data = comparison_source, filename = "S_D_cost_perturbation_vs_reward_null_30_75_source.csv")
)
for (entry in source_outputs) write_csv_new(entry$data, file.path(SOURCE_DIR, entry$filename))
write_csv_new(excluded, file.path(SOURCE_DIR, "excluded_disease_status.csv"))

# Keyed field backfill. Original bytes are archived before modifying the target table.
selection <- read_csv(SELECTION_TARGET)
assert_cols(selection, c("Disease", "parameter_id", "replicate", "repeat_record_path", "record_origin",
  "full_model_SHAP_raw_total", "mapped_PPI_SHAP_raw_total", "selected_SHAP_raw_sum",
  "retained_SHAP_mass_fraction_full_model", "retained_SHAP_mass_fraction_mapped_PPI"), "选参重复明细")
selection[, `:=`(Disease = trimws(as.character(Disease)), parameter_id = as.character(parameter_id),
                 replicate = as.integer(replicate), repeat_record_path = as.character(repeat_record_path))]
if (nrow(selection) != 1980L || anyDuplicated(selection[, .(Disease, parameter_id, replicate)]) || anyDuplicated(key_of(selection))) {
  stop("选参表行数或唯一键不符合冻结记录")
}
denom_by_disease <- qualification[, .(
  full_unique = uniqueN(full_model_SHAP_raw_total), mapped_unique = uniqueN(mapped_PPI_SHAP_raw_total),
  full_model_SHAP_raw_total = first(full_model_SHAP_raw_total), mapped_PPI_SHAP_raw_total = first(mapped_PPI_SHAP_raw_total)
), by = Disease]
if (nrow(denom_by_disease) != 16L || any(denom_by_disease$full_unique != 1L) ||
    any(denom_by_disease$mapped_unique != 1L) || any(!is.finite(denom_by_disease$full_model_SHAP_raw_total)) ||
    any(!is.finite(denom_by_disease$mapped_PPI_SHAP_raw_total)) || any(denom_by_disease$full_model_SHAP_raw_total <= 0) ||
    any(denom_by_disease$mapped_PPI_SHAP_raw_total <= 0)) stop("qualification表疾病SHAP分母跨参数不一致或无效")

target_cols <- c("full_model_SHAP_raw_total", "mapped_PPI_SHAP_raw_total", "selected_SHAP_raw_sum")
missing_mask <- Reduce(`|`, lapply(target_cols, function(nm) is.na(selection[[nm]])))
if (!all(vapply(target_cols, function(nm) identical(is.na(selection[[nm]]), missing_mask), logical(1))) || sum(missing_mask) != 1200L) {
  stop("仅预期1200条新重复三字段同时缺失；实际模式不同")
}
if (any(selection$record_origin[missing_mask] != "new_30_75_selection") ||
    any(vapply(target_cols, function(nm) anyNA(selection[[nm]][!missing_mask]), logical(1)))) {
  stop("缺失行来源或历史行状态不符")
}
missing_rows <- selection[missing_mask]
if (any(!file.exists(missing_rows$repeat_record_path))) stop("缺失重复记录的JSON路径不存在")
record_rows <- vector("list", nrow(missing_rows))
for (i in seq_len(nrow(missing_rows))) {
  path <- missing_rows$repeat_record_path[[i]]
  rec <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  selected_sum <- as.numeric(rec$network_metrics$selected_SHAP_raw_sum)
  if (length(selected_sum) != 1L || !is.finite(selected_sum) || selected_sum < 0) stop("重复JSON缺少有效network_metrics.selected_SHAP_raw_sum：", path)
  record_rows[[i]] <- data.table(Disease = missing_rows$Disease[[i]], parameter_id = missing_rows$parameter_id[[i]],
    replicate = missing_rows$replicate[[i]], repeat_record_path = path,
    selected_SHAP_raw_sum = selected_sum)
}
json_sums <- rbindlist(record_rows)
denoms <- merge(missing_rows[, .(Disease, parameter_id, replicate, repeat_record_path)],
  denom_by_disease[, .(Disease, full_model_SHAP_raw_total, mapped_PPI_SHAP_raw_total)], by = "Disease", all.x = TRUE)
derived <- merge(json_sums, denoms, by = c("Disease", "parameter_id", "replicate", "repeat_record_path"),
                 all = TRUE, sort = FALSE)
if (nrow(derived) != 1200L || anyNA(derived[, .(full_model_SHAP_raw_total, mapped_PPI_SHAP_raw_total, selected_SHAP_raw_sum)]) ||
    anyDuplicated(key_of(derived))) stop("按复合键关联JSON/分母失败")
derived[, `:=`(retained_full_check = selected_SHAP_raw_sum / full_model_SHAP_raw_total,
               retained_mapped_check = selected_SHAP_raw_sum / mapped_PPI_SHAP_raw_total)]
missing_with_existing <- merge(
  missing_rows[, .(Disease, parameter_id, replicate, repeat_record_path,
    retained_SHAP_mass_fraction_full_model, retained_SHAP_mass_fraction_mapped_PPI)],
  derived[, .(Disease, parameter_id, replicate, repeat_record_path,
    retained_full_check, retained_mapped_check)],
  by = c("Disease", "parameter_id", "replicate", "repeat_record_path"), all.x = TRUE, sort = FALSE
)
if (nrow(missing_with_existing) != 1200L ||
    any(abs(missing_with_existing$retained_full_check - missing_with_existing$retained_SHAP_mass_fraction_full_model) > 1e-10) ||
    any(abs(missing_with_existing$retained_mapped_check - missing_with_existing$retained_SHAP_mass_fraction_mapped_PPI) > 1e-10)) {
  stop("回填原始SHAP总量与已有full/mapped留存比例不一致")
}

archive_copy <- file.path(ARCHIVE_DIR, "selection_candidate_20repeats_pre_closeout.csv")
if (file.exists(archive_copy)) stop("回填前归档已存在，拒绝覆盖")
if (!file.copy(SELECTION_TARGET, archive_copy, overwrite = FALSE)) stop("无法归档回填前选参表")
old_selection_sha <- sha256_file(SELECTION_TARGET)
archive_sha <- sha256_file(archive_copy)
if (!identical(old_selection_sha, archive_sha)) stop("回填前归档SHA不匹配")
run_manifest_path <- file.path(RUN, "summaries", "run_manifest.json")
run_manifest_archive <- file.path(ARCHIVE_DIR, "run_manifest_pre_closeout.json")
if (file.exists(run_manifest_archive)) stop("阶段manifest归档已存在，拒绝覆盖")
if (!file.copy(run_manifest_path, run_manifest_archive, overwrite = FALSE) ||
    !identical(sha256_file(run_manifest_path), sha256_file(run_manifest_archive))) stop("阶段manifest归档失败")
run_manifest_sha_before <- sha256_file(run_manifest_path)

# Build a value map by the full frozen key; textual patching preserves every other CSV field byte-for-byte.
fill_map <- copy(derived[, .(Disease, parameter_id, replicate, repeat_record_path,
  full_model_SHAP_raw_total, mapped_PPI_SHAP_raw_total, selected_SHAP_raw_sum)])
fill_map[, key := key_of(fill_map)]
if (anyDuplicated(fill_map$key)) stop("派生字段复合键重复")
selection[, key := key_of(selection)]
map_idx <- match(selection$key, fill_map$key)
if (sum(missing_mask & !is.na(map_idx)) != 1200L || any(missing_mask & is.na(map_idx))) stop("无法将回填值按复合键定位回原记录")

raw_lines <- readLines(SELECTION_TARGET, warn = FALSE, encoding = "UTF-8")
if (length(raw_lines) != nrow(selection) + 1L || any(grepl('"', raw_lines, fixed = TRUE))) {
  stop("原选参CSV含引号/行数不符，停止文本级回填以保护其他字段")
}
header <- strsplit(raw_lines[[1L]], ",", fixed = TRUE)[[1L]]
col_idx <- match(target_cols, header)
if (anyNA(col_idx) || anyDuplicated(col_idx)) stop("目标字段表头不唯一")
patched_lines <- raw_lines
for (i in which(missing_mask)) {
  fields <- strsplit(raw_lines[[i + 1L]], ",", fixed = TRUE)[[1L]]
  if (length(fields) != length(header)) stop("原表字段数不一致：CSV record ", i)
  k <- map_idx[[i]]
  vals <- c(fill_map$full_model_SHAP_raw_total[[k]], fill_map$mapped_PPI_SHAP_raw_total[[k]], fill_map$selected_SHAP_raw_sum[[k]])
  fields[col_idx] <- vapply(vals, num_string, character(1))
  patched_lines[[i + 1L]] <- paste(fields, collapse = ",")
}
if (file.exists(paste0(SELECTION_TARGET, ".precloseout.tmp"))) stop("发现临时文件，停止避免覆盖")
tmp_selection <- paste0(SELECTION_TARGET, ".precloseout.tmp")
writeLines(enc2utf8(patched_lines), tmp_selection, useBytes = TRUE)
if (!file.rename(tmp_selection, SELECTION_TARGET)) stop("回填文件原子替换失败")
new_selection_sha <- sha256_file(SELECTION_TARGET)
patched <- read_csv(SELECTION_TARGET)
if (nrow(patched) != 1980L || anyDuplicated(patched[, .(Disease, parameter_id, replicate)]) ||
    any(vapply(target_cols, function(nm) anyNA(patched[[nm]]), logical(1)))) stop("回填后表的行数、键或字段完整性失败")
if (!identical(names(patched), names(selection)[names(selection) != "key"])) stop("回填意外改变列名或顺序")
for (nm in setdiff(names(patched), target_cols)) {
  a <- selection[[nm]]; b <- patched[[nm]]
  if (is.numeric(a) || is.integer(a)) {
    if (!isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 0, check.attributes = FALSE))) stop("回填改变了原字段：", nm)
  } else if (!identical(as.character(a), as.character(b))) stop("回填改变了原字段：", nm)
}
for (nm in target_cols) {
  if (any(as.numeric(patched[[nm]][!missing_mask]) != as.numeric(selection[[nm]][!missing_mask]))) stop("原非缺失行被改写：", nm)
}

# Keep the frozen stage flags and every other manifest field unchanged; refresh only the direct file hash entry.
manifest_lines <- readLines(run_manifest_path, warn = FALSE, encoding = "UTF-8")
path_line <- grep("selection_candidate_20repeats.csv", manifest_lines, fixed = TRUE)
if (length(path_line) != 1L || path_line + 2L > length(manifest_lines) ||
    !grepl('"sha256"', manifest_lines[[path_line + 1L]], fixed = TRUE) ||
    !grepl('"size_bytes"', manifest_lines[[path_line + 2L]], fixed = TRUE)) stop("阶段manifest中直接哈希条目结构不唯一")
old_manifest_sha <- sub('.*"sha256"[[:space:]]*:[[:space:]]*"([^"]+)".*', "\\1", manifest_lines[[path_line + 1L]])
old_manifest_size <- as.numeric(sub('.*"size_bytes"[[:space:]]*:[[:space:]]*([0-9]+).*', "\\1", manifest_lines[[path_line + 2L]]))
if (!identical(tolower(old_manifest_sha), tolower(sha256_file(archive_copy))) ||
    !is.finite(old_manifest_size) || old_manifest_size != file.info(archive_copy)$size) {
  stop("阶段manifest旧哈希/大小与回填前归档不一致")
}
manifest_lines[[path_line + 1L]] <- sub('"sha256"[[:space:]]*:[[:space:]]*"[^"]+"',
  paste0('"sha256": "', new_selection_sha, '"'), manifest_lines[[path_line + 1L]])
manifest_lines[[path_line + 2L]] <- sub('"size_bytes"[[:space:]]*:[[:space:]]*[0-9]+',
  paste0('"size_bytes": ', file.info(SELECTION_TARGET)$size), manifest_lines[[path_line + 2L]])
manifest_tmp <- paste0(run_manifest_path, ".closeout.tmp")
if (file.exists(manifest_tmp)) stop("发现阶段manifest临时文件，停止避免覆盖")
writeLines(enc2utf8(manifest_lines), manifest_tmp, useBytes = TRUE)
if (!file.rename(manifest_tmp, run_manifest_path)) stop("更新阶段manifest直接哈希失败")
updated_stage_manifest <- jsonlite::fromJSON(run_manifest_path, simplifyVector = FALSE)
if (!identical(updated_stage_manifest$drug_analysis_run, FALSE) || !identical(updated_stage_manifest$formal_figures_run, FALSE)) {
  stop("阶段manifest范围标记被意外改变")
}
manifest_record <- data.table(
  manifest_path = normalizePath(run_manifest_path, winslash = "/", mustWork = TRUE),
  manifest_sha256_before = run_manifest_sha_before,
  manifest_sha256_after = sha256_file(run_manifest_path),
  selection_table_sha256_before = old_selection_sha,
  selection_table_sha256_after = new_selection_sha,
  stage_flags_preserved = TRUE,
  update_scope = "Only the direct selection_candidate_20repeats.csv output SHA256 and size_bytes were refreshed; drug_analysis_run and formal_figures_run remain false because they describe the frozen selection stage."
)
write_csv_new(manifest_record, file.path(ARCHIVE_DIR, "stage_manifest_hash_refresh_audit.csv"))

archive_manifest <- data.table(
  artifact = c("selection_candidate_20repeats.csv before field backfill", "selection_candidate_20repeats.csv after field backfill"),
  path = c(normalizePath(archive_copy, winslash = "/", mustWork = TRUE), normalizePath(SELECTION_TARGET, winslash = "/", mustWork = TRUE)),
  sha256 = c(archive_sha, new_selection_sha),
  reason = c("Byte-preserved pre-closeout archive", "Only three previously missing raw SHAP total fields backfilled by Disease+parameter_id+replicate+repeat_record_path")
)
write_csv_new(archive_manifest, file.path(ARCHIVE_DIR, "selection_table_backfill_archive_manifest.csv"))
backfill_audit <- copy(derived)
backfill_audit[, `:=`(source_json_sha256 = vapply(repeat_record_path, sha256_file, character(1)),
  denominator_source = normalizePath(required[[6L]], winslash = "/", mustWork = TRUE),
  denominator_source_sha256 = sha256_file(required[[6L]]))]
write_csv_new(backfill_audit, file.path(SOURCE_DIR, "selection_candidate_20repeats_backfill_audit.csv"))

# Produce a current closeout hash manifest; leave the frozen stage run_manifest untouched.
all_outputs <- c(
  required[1:6], required[[8L]], run_manifest_path, archive_copy, run_manifest_archive,
  run_script_copy,
  file.path(ARCHIVE_DIR, "selection_table_backfill_archive_manifest.csv"),
  file.path(ARCHIVE_DIR, "stage_manifest_hash_refresh_audit.csv"),
  file.path(SOURCE_DIR, vapply(source_outputs, `[[`, character(1), "filename")),
  file.path(SOURCE_DIR, "excluded_disease_status.csv"), file.path(SOURCE_DIR, "selection_candidate_20repeats_backfill_audit.csv"),
  file.path(FIG_DIR, c("S_A_node_Jaccard_30_75.png", "S_A_node_Jaccard_30_75.pdf",
    "S_B_backbone_edge_Jaccard_30_75.png", "S_B_backbone_edge_Jaccard_30_75.pdf",
    "S_C_network_size_30_75.png", "S_C_network_size_30_75.pdf",
    "S_D_cost_perturbation_vs_reward_null_30_75.png", "S_D_cost_perturbation_vs_reward_null_30_75.pdf")),
  SCRIPT, PALETTE_PATH
)
all_outputs <- unique(all_outputs[file.exists(all_outputs)])
manifest_dt <- data.table(path = normalizePath(all_outputs, winslash = "/", mustWork = TRUE),
  sha256 = vapply(all_outputs, sha256_file, character(1)), size_bytes = as.numeric(file.info(all_outputs)$size))
write_csv_new(manifest_dt, file.path(SOURCE_DIR, "final_closeout_artifact_hash_manifest.csv"))

message("Closeout figures complete; selection rows=", nrow(patched), "; keyed backfills=", sum(missing_mask),
        "; original SHA=", old_selection_sha, "; updated SHA=", new_selection_sha)
