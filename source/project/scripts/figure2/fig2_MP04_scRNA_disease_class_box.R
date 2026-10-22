# MP04 scRNA disease-class descriptive boxplot.
# Descriptive only: no significance test, no donor-level inference, no NMF/usage recomputation.
# Reuses the verified sample-level source of the frozen 12MP run and cross-checks every
# sample against the frozen pseudobulk usage, metadata and comparability tables before plotting.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales) })

NMF_ROOT  <- "E:/micro_output/NMF"
PANEL_TAG <- Sys.getenv("FIG2_MP04BOX_TAG", unset = "MP04_scRNA_box_preview")
if (!grepl("^[A-Za-z0-9_-]+$", PANEL_TAG)) stop("Unsafe output tag")
# full  = legend + caption, 110 x 90 mm (standalone panel)
# frame = no legend, no caption, canvas set to the reserved black frame in
#         E:/micro_output/主图/fig2/fig2.ai (58.42 x 67.06 mm, w/h = 0.8712)
LAYOUT <- Sys.getenv("FIG2_MP04BOX_LAYOUT", unset = "full")
if (!LAYOUT %in% c("full", "frame")) stop("FIG2_MP04BOX_LAYOUT must be 'full' or 'frame'")
OUT_ROOT <- file.path(NMF_ROOT, paste0("figure2_final_20260929_", PANEL_TAG))
FIG_DIR  <- file.path(OUT_ROOT, "figures")
TAB_DIR  <- file.path(OUT_ROOT, "tables")
LOG_DIR  <- file.path(OUT_ROOT, "logs")

IN_RUN     <- "E:/micro_output/NMF/figure2_final_20260929_run07"
IN_SAMPLE  <- file.path(IN_RUN, "tables/B_MP04_sample_usage_common_source.csv")
IN_USAGE   <- file.path(NMF_ROOT, "source_data/PanelD_pseudobulk_MP_usage_raw_source.csv")
IN_META    <- file.path(NMF_ROOT, "tables/PanelD_pseudobulk_metadata.csv")
IN_COMP    <- file.path(NMF_ROOT, "stage2_comparability/sample_comparability.csv")
IN_DIS_PAL <- file.path(NMF_ROOT, "tables/Figure2_12MP_Disease17_palette.csv")
INPUTS <- c(IN_SAMPLE, IN_USAGE, IN_META, IN_COMP, IN_DIS_PAL)
if (!all(file.exists(INPUTS))) stop("Required frozen input missing: ", paste(INPUTS[!file.exists(INPUTS)], collapse = "; "))
if (dir.exists(OUT_ROOT) && length(list.files(OUT_ROOT, recursive = TRUE, all.files = TRUE, no.. = TRUE))) stop("Output directory is not empty; refusing to overwrite: ", OUT_ROOT)

# Deterministic string collation: this environment cannot set LC_COLLATE=C.UTF-8, and the
# frozen run assigned different region codes under different collations before this was pinned.
suppressWarnings(invisible(Sys.setlocale("LC_COLLATE", "C")))

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
sink(file.path(LOG_DIR, "fig2_MP04_scRNA_disease_class_box.log"), split = TRUE)
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

NEURO_DISEASES <- c("AD", "ALS", "ALS_FTD", "DLBD_PD", "EOAD", "FTD", "MS", "PSP")
TUMOR_DISEASES <- c("GBM", "LGG", "LCBM", "MB")
CLASS_LEVELS   <- c("Non-neoplastic neurological disease", "Tumor")
CTRL_REASON    <- "Ctrl belongs to neither disease class in this contract; the scRNA platform contributes only one Ctrl sample, which cannot form a control distribution."

run <- fread(IN_SAMPLE); usage <- fread(IN_USAGE); meta <- fread(IN_META)
comp <- fread(IN_COMP);  dis_pal <- fread(IN_DIS_PAL)

need_run <- c("SampleID_Global", "DatasetID", "platform", "Disease", "platform_source",
              "Brain_Region_Authority", "region_exact", "n_cells", "usage")
if (!all(need_run %in% names(run))) stop("Sample source lacks: ", paste(setdiff(need_run, names(run)), collapse = ", "))
if (nrow(run) != 432L || anyDuplicated(run$SampleID_Global)) stop("Sample source is not the frozen 432-sample universe")
if (!all(c("Disease", "Color_Hex") %in% names(dis_pal)) || anyDuplicated(dis_pal$Disease)) stop("Frozen Disease palette schema mismatch")
if (!all(c("PseudobulkID", "MP04") %in% names(usage))) stop("Frozen usage table lacks PseudobulkID or MP04")
if (!all(c("PseudobulkID", "SampleID_Global", "DatasetID", "Disease", "n_cells", "included_ge20") %in% names(meta))) stop("Frozen metadata schema mismatch")
if (!all(c("SampleID_Global", "sc_or_sn", "Brain_Region_Authority", "included_ge20") %in% names(comp))) stop("Frozen comparability schema mismatch")

# ---- per-sample cross-check against the frozen sources ----------------------
chk <- merge(run[, .(SampleID_Global, DatasetID_run = DatasetID, Disease_run = Disease,
                     n_cells_run = n_cells, usage, region_exact,
                     Brain_Region_Authority_run = Brain_Region_Authority)],
             meta[included_ge20 == TRUE, .(SampleID_Global, PseudobulkID,
                                           DatasetID_meta = DatasetID, Disease_meta = Disease,
                                           n_cells_meta = n_cells)],
             by = "SampleID_Global", all.x = TRUE)
chk <- merge(chk, usage[, .(PseudobulkID, MP04_frozen = MP04)], by = "PseudobulkID", all.x = TRUE)
chk <- merge(chk, comp[included_ge20 == TRUE, .(SampleID_Global, sc_or_sn, Brain_Region_Authority = Brain_Region_Authority)],
             by = "SampleID_Global", all.x = TRUE)
if (nrow(chk) != 432L || anyNA(chk$MP04_frozen) || anyNA(chk$sc_or_sn)) stop("Cross-check join lost samples or platform")

norm_region <- function(x) {
  x <- trimws(tolower(as.character(x)))
  x <- gsub("[[:space:]]+", "_", x)
  x[is.na(x) | x %in% c("", "not_reported", "unknown", "na", "n/a")] <- NA_character_
  x
}
# fwrite writes NA as an empty field, so read the round-tripped columns back as NA too
reg_from_run  <- norm_region(chk$region_exact)
reg_from_comp <- norm_region(chk$Brain_Region_Authority)
same_region <- (is.na(reg_from_comp) & is.na(reg_from_run)) |
  (!is.na(reg_from_comp) & !is.na(reg_from_run) & reg_from_comp == reg_from_run)
if (!all(same_region)) {
  bad <- which(!same_region)
  cat("Region mismatch examples (SampleID | from source | from comparability):\n")
  for (i in head(bad, 5)) cat("   ", chk$SampleID_Global[i], "|", chk$region_exact[i], "|", chk$Brain_Region_Authority[i], "\n")
  stop("Exact brain region disagrees between plotted source and frozen comparability for ", length(bad), " sample(s)")
}
if (sum(is.na(reg_from_run)) != 41L) warning("Expected 41 samples without an exact authority brain region; found ", sum(is.na(reg_from_run)))

ok <- c(
  dataset = identical(as.character(chk$DatasetID_run), as.character(chk$DatasetID_meta)),
  disease = identical(as.character(chk$Disease_run), as.character(chk$Disease_meta)),
  ncells  = identical(as.integer(chk$n_cells_run), as.integer(chk$n_cells_meta)),
  usage   = isTRUE(all.equal(chk$usage, chk$MP04_frozen, tolerance = 1e-12)),
  region  = all(same_region))
if (!all(ok)) stop("Frozen cross-check failed: ", paste(names(ok)[!ok], collapse = ", "))
cat("Cross-check vs frozen pseudobulk usage/metadata/comparability: 432/432 samples matched\n")
cat("  max |usage deviation| =", max(abs(chk$usage - chk$MP04_frozen)), "\n")
cat("  exact brain region agrees after the frozen normalisation for all 432 samples\n")

# ---- confirmed scRNA samples only -------------------------------------------
plat <- chk$sc_or_sn
plat[plat == "scRNA-seq"] <- "scRNA"
cat("Platform counts (all 432):", paste(names(table(plat)), as.integer(table(plat)), collapse = "; "), "\n")
if (any(!plat %in% c("scRNA", "snRNA"))) stop("Unresolved platform present; this panel is restricted to confirmed scRNA samples")

d <- chk[plat == "scRNA"]
cat("scRNA samples:", nrow(d), "| datasets:", uniqueN(d$DatasetID_run), "\n")
unmapped <- setdiff(unique(d$Disease_run), c(NEURO_DISEASES, TUMOR_DISEASES, "Ctrl"))
if (length(unmapped)) stop("scRNA disease label not mapped to a class; adjudication required: ", paste(unmapped, collapse = ", "))

d[, disease_class := NA_character_]
d[Disease_run %in% NEURO_DISEASES, disease_class := CLASS_LEVELS[1]]
d[Disease_run %in% TUMOR_DISEASES, disease_class := CLASS_LEVELS[2]]
d[, included := !is.na(disease_class)]
excluded <- d[included == FALSE]
if (nrow(excluded) && any(excluded$Disease_run != "Ctrl")) stop("A non-Ctrl sample was excluded; contract does not allow silent exclusion")
excluded[, exclusion_reason := CTRL_REASON]

cat("included:", sum(d$included), "| excluded:", nrow(excluded), "\n")
pl <- d[included == TRUE, .(n_samples = .N, n_diseases = uniqueN(Disease_run), n_datasets = uniqueN(DatasetID_run)),
        by = .(disease_class = as.character(disease_class))]
setorder(pl, disease_class)
print(pl)
if (sum(d$included) != 158L) stop("Included scRNA sample count is ", sum(d$included), ", expected 158 (66 non-neoplastic + 92 tumour)")

# ---- colours (project-frozen Disease palette) -------------------------------
# Legend follows the two classes rather than the alphabet, so the 8 non-neoplastic
# colours and the 4 tumour colours stay contiguous in the key.
present <- unique(d[included == TRUE, Disease_run])
used_diseases <- c(NEURO_DISEASES[NEURO_DISEASES %in% present], TUMOR_DISEASES[TUMOR_DISEASES %in% present])
if (!setequal(used_diseases, present)) stop("Disease class lists do not cover every plotted disease")
dis_colors <- setNames(dis_pal$Color_Hex, dis_pal$Disease)
miss <- setdiff(used_diseases, names(dis_colors))
if (length(miss)) stop("Disease lacks a frozen palette entry: ", paste(miss, collapse = ", "))
color_key <- dis_pal[Disease %in% used_diseases][order(match(Disease, used_diseases))][, in_this_panel := TRUE]
color_key_out <- rbind(color_key[, .(Disease, Color_Hex, in_this_panel)],
                       dis_pal[!Disease %in% used_diseases][, .(Disease, Color_Hex, in_this_panel = FALSE)])
setorder(color_key_out, -in_this_panel, Disease)
fwrite(color_key_out, file.path(TAB_DIR, "MP04_scRNA_disease_color_key.csv"))
cat("\nFrozen Disease colours used in this panel:\n"); print(color_key_out[in_this_panel == TRUE])

# ---- plot -------------------------------------------------------------------
pd <- copy(d[included == TRUE])
pd[, disease_class := factor(disease_class, levels = CLASS_LEVELS)]
pd[, Disease_run := factor(Disease_run, levels = used_diseases)]

ymax_raw <- max(pd$usage)
ymax <- signif(ymax_raw * 1.12, 2)
if (ymax <= ymax_raw * 1.05) ymax <- signif(ymax_raw * 1.10, 2)
ybreaks <- pretty(c(0, ymax), n = 5)
ybreaks <- ybreaks[ybreaks >= 0 & ybreaks <= ymax]

grp <- pd[, .(n_samples = .N, n_diseases = uniqueN(Disease_run), n_datasets = uniqueN(DatasetID_run),
              mean_usage = mean(usage), median_usage = median(usage),
              q1 = quantile(usage, .25, names = FALSE), q3 = quantile(usage, .75, names = FALSE)),
          by = disease_class]
grp[, iqr := q3 - q1]
setorder(grp, disease_class)
grp[, class_label := if (LAYOUT == "full")
      sprintf("%s\n(n = %d; %d diseases; %d dataset%s)",
              as.character(disease_class), n_samples, n_diseases, n_datasets,
              ifelse(n_datasets == 1L, "", "s"))
    else
      sprintf("%s\n(n = %d)", gsub(" ", "\n", as.character(disease_class)), n_samples)]

wrap <- function(x, width) paste(strwrap(x, width = width), collapse = "\n")
cap <- if (LAYOUT == "full") paste(
  wrap("Each point is one SampleID_Global (mean MP04 usage across its cells); no sample is removed. Box = Q1-Q3, centre line = median, whiskers = the most extreme value within 1.5 x IQR, black diamond = group mean.", 118),
  wrap("Descriptive sample-level comparison within one scRNA platform; disease class is confounded with dataset, so disease and dataset effects cannot be separated; usage is a relative programme contribution and cannot be read as absolute activity. Samples may be correlated within a donor; no donor-level inference is made.", 118),
  wrap(sprintf("All %d non-neoplastic neurological disease samples come from TUDDENHAM_2024, so that group reflects a single cohort.", sum(pd$disease_class == CLASS_LEVELS[1])), 118),
  sep = "\n") else NULL

p <- ggplot(pd, aes(x = disease_class, y = usage)) +
  # Box is drawn first so that every sample point stays visible on top; the boxplot's own
  # outlier symbol is suppressed because all observations are already plotted.
  geom_boxplot(width = .42, fill = "white", colour = "black", linewidth = .38, outlier.shape = NA) +
  geom_jitter(aes(colour = Disease_run), position = position_jitter(width = .21, height = 0, seed = 20260929),
              size = 1.25, alpha = .82) +
  # black diamond = group mean; a hairline white edge keeps it legible where it sits on
  # top of the black median line
  geom_point(data = grp, aes(x = disease_class, y = mean_usage), shape = 23, fill = "black",
             colour = "white", stroke = .45, size = 2.5, inherit.aes = FALSE) +
  scale_colour_manual(values = dis_colors, breaks = used_diseases, name = NULL,
                      guide = guide_legend(ncol = 6, byrow = FALSE)) +
  scale_x_discrete(labels = setNames(grp$class_label, as.character(grp$disease_class))) +
  scale_y_continuous(limits = c(0, ymax), breaks = ybreaks,
                     labels = scales::label_number(accuracy = .01, trim = TRUE), expand = c(0, 0)) +
  labs(x = NULL, y = "Relative MP04 usage", caption = cap) +
  theme_minimal(base_family = "Arial", base_size = 7) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.line = element_line(colour = "black", linewidth = .35),
        axis.ticks = element_line(colour = "black", linewidth = .35),
        axis.text.x = element_text(colour = "black", size = 6.4, lineheight = 1.05),
        axis.text.y = element_text(colour = "black", size = 6.2),
        axis.title.y = element_text(colour = "black", size = 7),
        plot.caption = element_text(colour = "black", size = 4.4, hjust = 0, lineheight = 1.15),
        legend.position = if (LAYOUT == "full") "bottom" else "none",
        legend.key.size = unit(3.0, "mm"),
        legend.text = element_text(size = 5.6), legend.margin = margin(t = -2),
        plot.margin = if (LAYOUT == "full") margin(4, 6, 4, 4) else margin(2, 3, 2, 2))

# frame layout: canvas equals the reserved black frame in fig2.ai (58.42 x 67.06 mm, w/h 0.8712)
FIG_W <- if (LAYOUT == "full") 110 else 58.42
FIG_H <- if (LAYOUT == "full") 90 else 67.06
base <- file.path(FIG_DIR, "MP04_scRNA_disease_class_box")
for (dev in list(c("png", NA), c("pdf", "cairo"))) {
  args <- if (is.na(dev[2])) list(filename = paste0(base, ".png"), plot = p, width = FIG_W, height = FIG_H, units = "mm", dpi = 600, bg = "white")
          else list(filename = paste0(base, ".pdf"), plot = p, width = FIG_W, height = FIG_H, units = "mm", device = cairo_pdf, bg = "white")
  withCallingHandlers(do.call(ggsave, args),
                      warning = function(w) { cat("RENDER WARNING:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
}

# ---- source and summary tables ----------------------------------------------
mk <- function(dt) data.table(
  SampleID_Global = dt$SampleID_Global, DatasetID = dt$DatasetID_run,
  platform = "scRNA", Disease = as.character(dt$Disease_run),
  Brain_Region_Authority = dt$Brain_Region_Authority_run, region_exact = dt$region_exact,
  n_cells = dt$n_cells_run, MP04_usage = dt$usage,
  disease_class = dt$disease_class, included = dt$included)
inc <- mk(d[included == TRUE]); inc[, exclusion_reason := NA_character_]
exc <- mk(excluded);           exc[, exclusion_reason := CTRL_REASON]
all_s <- rbind(inc, exc, fill = TRUE)
setorder(all_s, -included, disease_class, Disease, DatasetID, SampleID_Global)
fwrite(all_s, file.path(TAB_DIR, "MP04_scRNA_disease_class_source_narrow.csv"))

inc[, `:=`(plot_x_order = match(disease_class, CLASS_LEVELS), plot_y_axis_min = 0, plot_y_axis_max = ymax,
           point_colour = unname(dis_colors[Disease]),
           box_rule = "Q1-Q3 box; median centre line; whiskers to the most extreme value within 1.5 x IQR; outlier symbols suppressed because all points are plotted",
           mean_symbol = "black diamond = group mean", jitter_width = .21)]
setorder(inc, disease_class, Disease, DatasetID, SampleID_Global)
fwrite(inc, file.path(TAB_DIR, "MP04_scRNA_disease_class_source.csv"))

grp[, `:=`(plot_y_axis_max = ymax,
           box_rule = "Q1-Q3; median; whisker = most extreme value within 1.5 x IQR",
           mean_symbol = "black diamond = group mean")]
fwrite(grp[, .(disease_class, n_samples, n_diseases, n_datasets, mean_usage, median_usage, q1, q3, iqr,
               plot_y_axis_max, box_rule, mean_symbol)],
       file.path(TAB_DIR, "MP04_scRNA_disease_class_summary.csv"))

cnt <- d[included == TRUE, .(n_samples = .N),
         by = .(disease_class = as.character(disease_class), Disease = as.character(Disease_run), DatasetID = DatasetID_run)]
setorder(cnt, disease_class, Disease, DatasetID)
cnt[, pct_within_class := round(100 * n_samples / sum(n_samples), 2), by = disease_class]
fwrite(cnt, file.path(TAB_DIR, "MP04_scRNA_disease_by_dataset_counts.csv"))

cat("\nGroup summary:\n"); print(grp[, .(disease_class, n_samples, n_diseases, n_datasets, mean_usage, median_usage, q1, q3, iqr)])
cat("\nExcluded samples:\n"); print(all_s[included == FALSE, .(SampleID_Global, Disease, DatasetID)])
cat("\nY axis: 0 to", ymax, "; breaks:", paste(ybreaks, collapse = ", "), "; observed max:", ymax_raw, "\n")
cat("Figure size:", FIG_W, "x", FIG_H, "mm at 600 dpi\n")
cat("Outputs:", OUT_ROOT, "\n")
graphics.off()
cat("Graphics devices after explicit close:", length(dev.list()), "\n")
cat("End:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
sink()
