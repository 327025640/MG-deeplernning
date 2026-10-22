from __future__ import annotations

import copy
import gc
import hashlib
import json
import os
import random
import shutil
import sys
import time
import traceback
from pathlib import Path
from types import SimpleNamespace

SCRATCH = Path(r"E:\micro_output\分析审阅中心\待审阅\Fig6_疾病vsCtrl_全基因CauFinder\scratch")
NUMBA_CACHE = SCRATCH / "numba_cache"
PY_TEMP = SCRATCH / "python_temp"
NUMBA_CACHE.mkdir(parents=True, exist_ok=True)
PY_TEMP.mkdir(parents=True, exist_ok=True)
os.environ["NUMBA_CACHE_DIR"] = str(NUMBA_CACHE)
os.environ["TEMP"] = str(PY_TEMP)
os.environ["TMP"] = str(PY_TEMP)

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import shap
import torch
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    precision_score,
    recall_score,
    roc_auc_score,
)

ROOT = Path(r"E:\micro_output\分析审阅中心\待审阅\Fig6_疾病vsCtrl_全基因CauFinder")
INPUT = ROOT / "input"
BLOCK_DIR = INPUT / "expression_csr_blocks"
ADAPTER = ROOT / "adapter"
TABLES = ROOT / "tables"
MODELS = ROOT / "models"
LOGS = ROOT / "logs"
TMP = ROOT / "scratch"
for path in (TABLES, MODELS, LOGS, TMP):
    path.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(ADAPTER))
from CauFinder import caufinder_main  # noqa: E402
from CauFinder.caufinder_main import CausalFinder, ShapModel  # noqa: E402
from CauFinder.models import DCDVAE  # noqa: E402
from CauFinder.utils import set_seed  # noqa: E402

SEED = 42
BATCH_SIZE = 128
SHAP_BATCH_SIZE = 8
INFERENCE_BATCH_SIZE = 512
GRAD_BATCH_SIZE = 128
GRAD_GENE_BLOCK = 128
N_ALPHA = 5
N_BETA = 10
MODEL_KWARGS = {
    "n_latent": 25,
    "n_causal": 5,
    "n_hidden": 128,
    "n_layers_encoder": 0,
    "n_layers_decoder": 0,
    "n_layers_dpd": 0,
    "dropout_rate_encoder": 0.0,
    "dropout_rate_decoder": 0.0,
    "dropout_rate_dpd": 0.0,
    "use_batch_norm": "none",
    "use_batch_norm_dpd": True,
    "pdp_linear": True,
}


def log(event: str, **values) -> None:
    record = {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "event": event, **values}
    line = json.dumps(record, ensure_ascii=False, default=str)
    with (LOGS / "run_events.jsonl").open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")
        handle.flush()
    print(line, flush=True)


def available_memory_gb() -> float | None:
    try:
        import psutil

        return round(psutil.virtual_memory().available / (1024**3), 2)
    except Exception:
        try:
            import ctypes

            class MemoryStatus(ctypes.Structure):
                _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                            ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                            ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                            ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                            ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]

            status = MemoryStatus()
            status.dwLength = ctypes.sizeof(status)
            if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
                return round(status.ullAvailPhys / (1024**3), 2)
        except Exception:
            pass
        return None


def process_rss_gb() -> float | None:
    try:
        import psutil

        return round(psutil.Process().memory_info().rss / (1024**3), 2)
    except Exception:
        try:
            import ctypes

            class ProcessMemoryCounters(ctypes.Structure):
                _fields_ = [("cb", ctypes.c_ulong), ("PageFaultCount", ctypes.c_ulong),
                            ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                            ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                            ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                            ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]

            counters = ProcessMemoryCounters()
            counters.cb = ctypes.sizeof(counters)
            kernel32 = ctypes.windll.kernel32
            kernel32.GetCurrentProcess.restype = ctypes.c_void_p
            psapi = ctypes.windll.psapi
            psapi.GetProcessMemoryInfo.argtypes = [ctypes.c_void_p, ctypes.POINTER(ProcessMemoryCounters), ctypes.c_ulong]
            psapi.GetProcessMemoryInfo.restype = ctypes.c_int
            if psapi.GetProcessMemoryInfo(kernel32.GetCurrentProcess(), ctypes.byref(counters), counters.cb):
                return round(counters.WorkingSetSize / (1024**3), 2)
        except Exception:
            pass
        return None


def csv_write(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, encoding="utf-8-sig")


class DiskCSRStore:
    """Read-only memmapped CSR blocks; only requested minibatches become dense."""

    def __init__(self, block_manifest: pd.DataFrame, block_size: int, n_genes: int):
        self.block_size = int(block_size)
        self.n_genes = int(n_genes)
        self.blocks = {}
        for row in block_manifest.itertuples(index=False):
            bid = int(row.block_id)
            n_rows = int(row.n_rows)
            nnz = int(row.nnz)
            self.blocks[bid] = {
                "start": int(row.global_row_start),
                "n_rows": n_rows,
                "indptr": np.memmap(row.indptr_file, dtype="<i4", mode="r", shape=(n_rows + 1,)),
                "indices": (np.memmap(row.indices_file, dtype="<i4", mode="r", shape=(nnz,))
                            if nnz else np.empty(0, dtype="<i4")),
                "data": (np.memmap(row.data_file, dtype="<f4", mode="r", shape=(nnz,))
                         if nnz else np.empty(0, dtype="<f4")),
            }

    def dense(self, global_rows: np.ndarray) -> np.ndarray:
        rows = np.asarray(global_rows, dtype=np.int64)
        out = np.zeros((len(rows), self.n_genes), dtype=np.float32)
        block_ids = rows // self.block_size
        for bid in np.unique(block_ids):
            positions = np.flatnonzero(block_ids == bid)
            block = self.blocks[int(bid)]
            local_rows = rows[positions] - block["start"]
            for out_row, local_row in zip(positions, local_rows):
                a = int(block["indptr"][local_row])
                b = int(block["indptr"][local_row + 1])
                if b > a:
                    out[out_row, block["indices"][a:b]] = block["data"][a:b]
        return out


class MemoryBoundBatch:
    def __init__(self, X: np.ndarray, obs: pd.DataFrame):
        self.X = X
        self.obs = obs
        self.shape = X.shape


input_started = time.perf_counter()
manifest = pd.read_csv(INPUT / "cell_manifest.csv.gz")
gene_manifest_all = pd.read_csv(INPUT / "gene_manifest_all_RNA.csv")
gene_annotation = pd.read_csv(INPUT / "gene_manifest_model.csv")
genes = gene_annotation["gene_id"].astype(str).tolist()
block_manifest = pd.read_csv(INPUT / "csr_block_manifest.csv")
for block_column in ("indptr_file", "indices_file", "data_file"):
    block_manifest[block_column] = block_manifest[block_column].map(
        lambda name: str(BLOCK_DIR / Path(str(name)).name)
    )
coverage_table = pd.read_csv(TABLES / "disease_coverage.csv")
if manifest.empty or not genes or manifest["GlobalRow"].duplicated().any():
    raise RuntimeError("输入清单为空或 GlobalRow 不唯一。")
if len(genes) == 0 or manifest["GlobalRow"].min() != 0 or manifest["GlobalRow"].max() != len(manifest) - 1:
    raise RuntimeError("输入空间为空或 GlobalRow 不连续；请复核正式输入清单。")
if not manifest["Disease"].eq("Ctrl").any() or manifest["Disease"].nunique() < 2:
    raise RuntimeError("正式训练输入缺少 Ctrl 或目标疾病。")

cell_count = len(manifest)
store = DiskCSRStore(block_manifest, int(block_manifest.n_rows.max()), len(genes))
for bid, block in store.blocks.items():
    ptr = block["indptr"]
    if int(ptr[0]) != 0 or int(ptr[-1]) != len(block["indices"]) or np.any(np.diff(ptr) < 0):
        raise RuntimeError(f"CSR block {bid} 的 row pointer 损坏。")
    if len(block["indices"]) and (int(block["indices"].min()) < 0 or int(block["indices"].max()) >= len(genes)):
        raise RuntimeError(f"CSR block {bid} 的 feature index 越界。")
reference = pd.read_csv(INPUT / "csr_export_validation_reference.csv")
check_rows = reference.GlobalRow.drop_duplicates().to_numpy(np.int64)
check_genes = reference.gene_id.drop_duplicates().astype(str).tolist()
check_X = store.dense(check_rows)
row_pos = {int(row): i for i, row in enumerate(check_rows)}
gene_pos = {gene: i for i, gene in enumerate(genes)}
expected = reference.apply(lambda r: np.float32(r.RNA_data_value), axis=1).to_numpy(np.float32)
observed = np.array([check_X[row_pos[int(r.GlobalRow)], gene_pos[str(r.gene_id)]]
                     for r in reference.itertuples(index=False)], dtype=np.float32)
if not np.array_equal(expected, observed):
    raise RuntimeError("CSR 缓存与 RNA data layer 抽样核对不一致。")
input_validation = {
    "n_cells": int(len(manifest)),
    "n_diseases": int(manifest.Disease.nunique()),
    "n_all_RNA_genes": int(len(gene_manifest_all)),
    "n_model_genes_after_global_zero_variance_filter": int(len(genes)),
    "n_removed_all_zero_or_near_zero_variance": int(len(gene_manifest_all) - len(genes)),
    "n_csr_blocks": int(len(block_manifest)),
    "csr_block_pointers_valid": True,
    "RNA_data_layer_roundtrip_rows": int(len(reference)),
    "RNA_data_layer_roundtrip_float32_exact": True,
    "input_read_validation_elapsed_s": time.perf_counter() - input_started,
}
(INPUT / "CSR_export_validation.json").write_text(json.dumps(input_validation, indent=2), encoding="utf-8")
log("input_read_validation_complete", **input_validation,
    available_memory_gb=available_memory_gb(), process_rss_gb=process_rss_gb())
diseases = sorted(set(manifest.Disease.astype(str).unique()) - {"Ctrl"})
if not diseases:
    raise RuntimeError("输入清单没有非Ctrl疾病。")

# The trainable metadata-only AnnData carries a tiny all-zero CSR placeholder.
# Actual expression is supplied by the on-disk CSR batch sampler below.
def no_cell_resplit(adata, train_size=1.0, validation_size=None, use_gpu=None, **kwargs):
    return adata, None


def disk_batch_sampler(adata, batch_size, shuffle=True, drop_last=False):
    n_obs = adata.n_obs
    order = np.random.permutation(n_obs) if shuffle else np.arange(n_obs)
    stop = n_obs - (n_obs % batch_size) if drop_last else n_obs
    row_ids = adata.obs["global_row"].to_numpy(np.int64)
    for start in range(0, stop, batch_size):
        batch_ix = order[start : min(start + batch_size, stop)]
        if drop_last and len(batch_ix) < batch_size:
            continue
        yield MemoryBoundBatch(store.dense(row_ids[batch_ix]), adata.obs.iloc[batch_ix])


caufinder_main.data_splitter = no_cell_resplit
caufinder_main.batch_sampler = disk_batch_sampler


def predict_prob(module: torch.nn.Module, rows: np.ndarray, batch_size: int) -> np.ndarray:
    module.eval()
    device = next(module.parameters()).device
    probs = np.empty(len(rows), dtype=np.float32)
    with torch.no_grad():
        for start in range(0, len(rows), batch_size):
            batch_rows = rows[start : start + batch_size]
            x = torch.from_numpy(store.dense(batch_rows)).to(device)
            out = module(x, use_mean=True)
            probs[start : start + len(batch_rows)] = out["alpha_dpd"]["prob"].detach().cpu().numpy().reshape(-1)
    return probs


def auc_metrics(y: np.ndarray, p: np.ndarray) -> tuple[float, float]:
    if np.unique(y).size < 2:
        return np.nan, np.nan
    return float(roc_auc_score(y, p)), float(average_precision_score(y, p))


def train_adata_for(disease: str) -> ad.AnnData:
    sub = manifest.loc[(manifest.Split == "train") & manifest.Disease.isin([disease, "Ctrl"]),
                       ["GlobalRow", "CellID", "Disease"]].copy()
    if set(sub.Disease.astype(str).unique()) - {disease, "Ctrl"}:
        raise RuntimeError(f"{disease}: other diseases entered the binary training AnnData")
    sub["global_row"] = sub.GlobalRow.astype(np.int64)
    sub["labels"] = (sub.Disease == disease).astype(np.float32)
    sub.index = sub.CellID.astype(str)
    var = pd.DataFrame(index=pd.Index(genes, name="gene_id"))
    zero_placeholder = sp.csr_matrix((len(sub), len(genes)), dtype=np.float32)
    return ad.AnnData(X=zero_placeholder, obs=sub[["GlobalRow", "global_row", "Disease", "labels"]], var=var)


def evaluate_split(disease: str, split: str, rows: np.ndarray, probabilities: np.ndarray):
    part = manifest.set_index("GlobalRow").loc[rows].reset_index()
    y = (part.Disease.to_numpy() == disease).astype(np.int8)
    predicted = (probabilities >= 0.5).astype(np.int8)
    auroc, auprc = auc_metrics(y, probabilities)
    tn, fp, fn, tp = confusion_matrix(y, predicted, labels=[0, 1]).ravel()
    metrics = {
        "Disease": disease,
        "split": split,
        "n_cells": len(y),
        "n_positive": int(y.sum()),
        "positive_fraction": float(y.mean()),
        "AUROC": auroc,
        "AUPRC": auprc,
        "AUPRC_positive_baseline": float(y.mean()),
        "recall": float(recall_score(y, predicted, zero_division=0)),
        "precision": float(precision_score(y, predicted, zero_division=0)),
        "balanced_accuracy": float(balanced_accuracy_score(y, predicted)),
        "threshold": 0.5,
        "TN": int(tn), "FP": int(fp), "FN": int(fn), "TP": int(tp),
    }
    pred = part[["CellID", "Disease", "DatasetID", "SampleID_Global", "Cluster", "Split"]].copy()
    pred["target_disease"] = disease
    pred["y_true"] = y
    pred["p_target"] = probabilities
    pred["y_pred_threshold_0_5"] = predicted
    ref_rows = []
    for ref_disease, grp in pred.groupby("Disease", sort=True):
        gprob = grp.p_target.to_numpy()
        gpred = grp.y_pred_threshold_0_5.to_numpy()
        is_target = ref_disease == disease
        ref_rows.append({
            "target_disease": disease,
            "split": split,
            "reference_disease": ref_disease,
            "actual_class": int(is_target),
            "n_cells": len(grp),
            "mean_p_target": float(gprob.mean()),
            "median_p_target": float(np.median(gprob)),
            "predicted_positive_fraction": float(gpred.mean()),
            "false_negative_cells": int((gpred == 0).sum()) if is_target else 0,
            "false_positive_cells": int(gpred.sum()) if not is_target else 0,
        })
    return metrics, pred, ref_rows


def shap_summaries(disease: str, module: torch.nn.Module, train_rows: np.ndarray,
                   y_train: np.ndarray, train_diseases: np.ndarray, out_dir: Path):
    rng = np.random.default_rng(SEED)
    pos = train_rows[y_train == 1]
    neg = train_rows[y_train == 0]
    n_pos_bg = min(50, len(pos))
    n_neg_bg = min(50, len(neg))
    bg_rows = np.concatenate([
        rng.choice(pos, n_pos_bg, replace=False) if n_pos_bg else np.array([], dtype=np.int64),
        rng.choice(neg, n_neg_bg, replace=False) if n_neg_bg else np.array([], dtype=np.int64),
    ])
    if len(bg_rows) == 0:
        raise RuntimeError("SHAP 背景细胞为空。")
    csv_write(manifest.set_index("GlobalRow").loc[bg_rows].reset_index(),
              out_dir / "shap_background_cells.csv")
    csv_write(manifest.set_index("GlobalRow").loc[train_rows].reset_index(),
              out_dir / "shap_explanation_cells.csv.gz")
    log("shap_background", disease=disease, requested=100, actual=len(bg_rows),
        positive=n_pos_bg, negative=n_neg_bg, available_memory_gb=available_memory_gb())
    module.eval()
    device = next(module.parameters()).device
    background = torch.from_numpy(store.dense(bg_rows)).to(device)
    explainer = shap.DeepExplainer(ShapModel(module, "prob").to(device), background)
    group_names = [f"disease:{disease}", "disease:Ctrl"]
    abs_sum = {name: np.zeros(len(genes), dtype=np.float64) for name in group_names}
    signed_sum = {name: np.zeros(len(genes), dtype=np.float64) for name in group_names}
    group_n = {name: 0 for name in group_names}
    started = time.perf_counter()
    for start in range(0, len(train_rows), SHAP_BATCH_SIZE):
        rows = train_rows[start : start + SHAP_BATCH_SIZE]
        X = store.dense(rows)
        shap_raw = explainer.shap_values(torch.from_numpy(X).to(device), check_additivity=False)
        if isinstance(shap_raw, list):
            shap_raw = shap_raw[0]
        values = np.asarray(shap_raw, dtype=np.float32)
        if values.ndim == 3:
            values = values[..., 0]
        if values.shape != X.shape:
            raise RuntimeError(f"SHAP 维度异常: {values.shape} vs {X.shape}")
        local_diseases = train_diseases[start : start + len(rows)]
        local_y = y_train[start : start + len(rows)]
        masks = {f"disease:{disease}": local_diseases == disease,
                 "disease:Ctrl": local_diseases == "Ctrl"}
        for name, mask in masks.items():
            if np.any(mask):
                block = values[mask]
                abs_sum[name] += np.abs(block).sum(axis=0, dtype=np.float64)
                signed_sum[name] += block.sum(axis=0, dtype=np.float64)
                group_n[name] += int(mask.sum())
        if (start // SHAP_BATCH_SIZE + 1) % 500 == 0:
            log("shap_progress", disease=disease, completed=min(start + len(rows), len(train_rows)),
                total=len(train_rows), elapsed_s=round(time.perf_counter() - started, 1),
                available_memory_gb=available_memory_gb(), process_rss_gb=process_rss_gb())
        del X, values, shap_raw
    shap_compute_s = time.perf_counter() - started
    summary_rows = []
    for name in group_names:
        n = group_n[name]
        if n == 0:
            continue
        mean_abs = abs_sum[name] / n
        total = mean_abs.sum()
        normalized = mean_abs / total if total > 0 else np.zeros_like(mean_abs)
        for gene, raw, norm_value, signed in zip(genes, mean_abs, normalized, signed_sum[name] / n):
            summary_rows.append({
                "target_disease": disease,
                "contribution_group": name,
                "n_cells": n,
                "gene_id": gene,
                "mean_abs_SHAP": float(raw),
                "mean_abs_SHAP_normalized": float(norm_value),
                "mean_signed_SHAP_descriptive_only": float(signed),
            })
    write_started = time.perf_counter()
    csv_write(pd.DataFrame(summary_rows), out_dir / "shap_importance_all_groups.csv")
    shap_write_s = time.perf_counter() - write_started
    return group_n, shap_compute_s, shap_write_s


def gradient_summary(disease: str, module: torch.nn.Module, rows: np.ndarray, out_dir: Path):
    cache_file = TMP / f"{disease}_target_cells_gradients.f32"
    gradient_map = np.memmap(cache_file, dtype="<f4", mode="w+", shape=(len(rows), len(genes)))
    module.eval()
    device = next(module.parameters()).device
    sum_grad = np.zeros(len(genes), dtype=np.float64)
    pos_count = np.zeros(len(genes), dtype=np.int64)
    neg_count = np.zeros(len(genes), dtype=np.int64)
    zero_count = np.zeros(len(genes), dtype=np.int64)
    started = time.perf_counter()
    for start in range(0, len(rows), GRAD_BATCH_SIZE):
        batch_rows = rows[start : start + GRAD_BATCH_SIZE]
        x_np = store.dense(batch_rows)
        x = torch.from_numpy(x_np).to(device).requires_grad_(True)
        outputs = module(x, use_mean=True)
        probability = outputs["alpha_dpd"]["prob"]
        gradients = torch.autograd.grad(probability.sum(), x, retain_graph=False, create_graph=False)[0]
        g = gradients.detach().cpu().numpy().astype(np.float32, copy=False)
        gradient_map[start : start + len(g)] = g
        sum_grad += g.sum(axis=0, dtype=np.float64)
        pos_count += (g > 0).sum(axis=0)
        neg_count += (g < 0).sum(axis=0)
        zero_count += (g == 0).sum(axis=0)
        if (start // GRAD_BATCH_SIZE + 1) % 500 == 0:
            log("gradient_progress", disease=disease, completed=min(start + len(g), len(rows)),
                total=len(rows), elapsed_s=round(time.perf_counter() - started, 1),
                available_memory_gb=available_memory_gb(), process_rss_gb=process_rss_gb())
        del x_np, x, outputs, probability, gradients, g
    gradient_map.flush()
    medians = np.empty(len(genes), dtype=np.float64)
    for g_start in range(0, len(genes), GRAD_GENE_BLOCK):
        g_end = min(g_start + GRAD_GENE_BLOCK, len(genes))
        medians[g_start:g_end] = np.median(np.asarray(gradient_map[:, g_start:g_end]), axis=0)
    n = len(rows)
    gradient_compute_s = time.perf_counter() - started
    frame = pd.DataFrame({
        "target_disease": disease,
        "gene_id": genes,
        "n_target_disease_cells": n,
        "gradient_mean_dP_dRNAdata": sum_grad / n,
        "gradient_median_dP_dRNAdata": medians,
        "gradient_positive_fraction": pos_count / n,
        "gradient_negative_fraction": neg_count / n,
        "gradient_zero_fraction": zero_count / n,
        "direction_interpretation": np.where(
            (pos_count > 0) & (neg_count > 0), "mixed_sign",
            np.where(pos_count > 0, "positive_local_gradient",
                     np.where(neg_count > 0, "negative_local_gradient", "zero"))),
    })
    write_started = time.perf_counter()
    csv_write(frame, out_dir / "causal_branch_gradient_direction.csv")
    gradient_write_s = time.perf_counter() - write_started
    gradient_map.flush()
    del gradient_map
    gc.collect()
    cache_file.unlink(missing_ok=True)
    return gradient_compute_s, gradient_write_s


def stage_timing(history: dict, disease: str, out_dir: Path):
    write_started = time.perf_counter()
    frame = pd.DataFrame(history)
    csv_write(frame.assign(Disease=disease, epoch=np.arange(1, len(frame) + 1)),
              out_dir / "training_log.csv")
    weight_cols = [c for c in frame.columns if c.startswith("weight_")]
    timing_rows = []
    if weight_cols:
        keys = frame[weight_cols].round(10).astype(str).agg("|".join, axis=1)
        group_id = keys.ne(keys.shift()).cumsum()
        for _, part in frame.groupby(group_id, sort=True):
            timing_rows.append({
                "Disease": disease,
                "epoch_start": int(part.index.min() + 1),
                "epoch_end": int(part.index.max() + 1),
                "n_epochs": len(part),
                "training_elapsed_s": float(part.epoch_duration_s.sum()),
                "loss_weight_schedule": json.dumps(part.iloc[0][weight_cols].to_dict()),
            })
    csv_write(pd.DataFrame(timing_rows), out_dir / "training_stage_timing.csv")
    return frame, time.perf_counter() - write_started


def run_one(disease: str, core38: set[str], ptg_sets: dict[str, set[str]]) -> dict:
    out_dir = MODELS / disease
    out_dir.mkdir(parents=True, exist_ok=True)
    started_model = time.perf_counter()
    train_frame = manifest.loc[manifest.Split == "train"]
    train_rows = train_frame.GlobalRow.to_numpy(np.int64)
    train_diseases = train_frame.Disease.to_numpy(str)
    y_train = (train_diseases == disease).astype(np.int8)
    n_positive = int(y_train.sum())
    n_negative = int(len(y_train) - n_positive)
    if n_positive == 0 or n_negative == 0:
        raise RuntimeError(f"{disease}: 训练集缺少一类，无法拟合。")
    n_train = len(train_rows)
    w_positive = n_train / (2.0 * n_positive)
    w_negative = n_train / (2.0 * n_negative)
    class_weights = [w_negative, w_positive]
    auto_epochs = min(round((20000 / n_train) * 400), 400)
    val_rows = manifest.loc[manifest.Split == "validation", "GlobalRow"].to_numpy(np.int64)
    val_y = (manifest.set_index("GlobalRow").loc[val_rows].Disease.to_numpy() == disease).astype(np.int8)
    best = {"auroc": -np.inf, "auprc": -np.inf, "epoch": None}
    checkpoint_write_s = 0.0
    validation_elapsed_s = 0.0
    ckpt_file = out_dir / "best_validation_checkpoint.pt"
    epoch_log = LOGS / f"{disease}_epoch_training.jsonl"
    if epoch_log.exists():
        epoch_log.unlink()

    adata_train = train_adata_for(disease)
    set_seed(SEED)
    model = CausalFinder(adata=adata_train, **MODEL_KWARGS)

    config = {
        "Disease": disease,
        "seed": SEED,
        "split_unit": "cell",
        "split_proportions_requested": {"train": 0.64, "validation": 0.16, "final_eval": 0.20},
        "n_train_cells": n_train,
        "n_positive_train": n_positive,
        "n_negative_train": n_negative,
        "w_positive": w_positive,
        "w_negative": w_negative,
        "class_weight_application": "fixed per-cell BCE multiplier; one application; batch mean",
        "class_weights_order": ["negative", "positive"],
        "features": len(genes),
        "model_architecture": MODEL_KWARGS,
        "training": {"max_epochs": auto_epochs, "auto_rule": "round((20000/n_train_cells)*400), capped at 400",
                     "stage_training": True, "n_alpha": N_ALPHA, "n_beta": N_BETA,
                     "batch_size": BATCH_SIZE, "optimizer": "Adam", "lr": 0.0005},
        "checkpoint_selection": "maximum validation AUROC; validation AUPRC as tie-break; fixed before fitting",
        "classification_threshold": 0.5,
        "input_scale": "frozen Seurat RNA data layer; float32 storage; no recomputation",
        "gpu_available": torch.cuda.is_available(),
    }
    (out_dir / "model_config.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")

    def validation_callback(module, epoch, epoch_losses):
        nonlocal best, checkpoint_write_s, validation_elapsed_s
        validation_started = time.perf_counter()
        p = predict_prob(module, val_rows, INFERENCE_BATCH_SIZE)
        validation_elapsed_s += time.perf_counter() - validation_started
        auroc, auprc = auc_metrics(val_y, p)
        if np.isnan(auroc):
            raise RuntimeError(f"{disease}: validation AUROC undefined; validation lacks one class.")
        is_better = (auroc > best["auroc"]) or (np.isclose(auroc, best["auroc"]) and auprc > best["auprc"])
        if is_better:
            best = {"auroc": auroc, "auprc": auprc, "epoch": int(epoch + 1)}
            state = {k: v.detach().cpu().clone() for k, v in module.state_dict().items()}
            tmp_path = ckpt_file.with_suffix(".tmp")
            checkpoint_started = time.perf_counter()
            torch.save({"state_dict": state, "best_validation": best, "config": config}, tmp_path)
            os.replace(tmp_path, ckpt_file)
            checkpoint_write_s += time.perf_counter() - checkpoint_started
        row = {"Disease": disease, "epoch": epoch + 1,
               "val_AUROC": auroc, "val_AUPRC": auprc,
               "is_best_checkpoint": bool(is_better),
               "best_epoch_so_far": best["epoch"],
               "best_val_AUROC_so_far": best["auroc"],
               "best_val_AUPRC_so_far": best["auprc"]}
        for key, values in epoch_losses.items():
            if values and key not in {"val_AUROC", "val_AUPRC"}:
                row[key] = values[-1]
        with epoch_log.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")
            handle.flush()
        log("validation_checkpoint", **row, available_memory_gb=available_memory_gb(), process_rss_gb=process_rss_gb())
        return {"val_AUROC": auroc, "val_AUPRC": auprc}

    train_started = time.perf_counter()
    log("training_start", disease=disease, n_train=n_train, n_positive=n_positive,
        n_negative=n_negative, w_positive=w_positive, w_negative=w_negative,
        expected_auto_epochs=auto_epochs, available_memory_gb=available_memory_gb(),
        process_rss_gb=process_rss_gb())
    model.train(max_epochs=None, stage_training=True, train_size=1.0, validation_size=None,
                batch_size=BATCH_SIZE, drop_last=False, n_alpha=N_ALPHA, n_beta=N_BETA,
                class_weights=class_weights, epoch_callback=validation_callback)
    training_s = time.perf_counter() - train_started
    if not ckpt_file.exists():
        raise RuntimeError(f"{disease}: validation 未产生可选检查点。")
    checkpoint = torch.load(ckpt_file, map_location="cpu", weights_only=False)
    model.module.load_state_dict(checkpoint["state_dict"])
    model.module.eval()
    history, training_log_write_s = stage_timing(model.history, disease, out_dir)
    log("training_complete", disease=disease, elapsed_s=round(training_s, 2),
        best_epoch=best["epoch"], best_val_AUROC=best["auroc"], best_val_AUPRC=best["auprc"],
        n_epochs=len(history), available_memory_gb=available_memory_gb(), process_rss_gb=process_rss_gb())

    # Save train, validation, and untouched final-evaluation predictions.
    predictions, metrics, by_reference = [], [], []
    split_rows = {"train": all_train_rows, "validation": all_val_rows, "final_eval": all_eval_rows}
    prediction_started = time.perf_counter()
    for split, rows in split_rows.items():
        p = predict_prob(model.module, rows, INFERENCE_BATCH_SIZE)
        metric, pred, refs = evaluate_split(disease, split, rows, p)
        predictions.append(pred)
        metrics.append(metric)
        by_reference.extend(refs)
    prediction_s = time.perf_counter() - prediction_started
    pred_all = pd.concat(predictions, ignore_index=True)
    prediction_write_started = time.perf_counter()
    csv_write(pred_all, TABLES / f"{disease}_predictions_cells.csv.gz")
    csv_write(pd.DataFrame(metrics), TABLES / f"{disease}_metrics.csv")
    csv_write(pd.DataFrame(by_reference), TABLES / f"{disease}_by_reference_disease.csv")
    prediction_write_s = time.perf_counter() - prediction_write_started

    # Full training-set Deep SHAP, accumulated in memory by group and gene.
    shap_started = time.perf_counter()
    shap_groups, shap_s, shap_write_s = shap_summaries(disease, model.module, train_rows, y_train,
                                                       train_diseases, out_dir)
    shap_frame = pd.read_csv(out_dir / "shap_importance_all_groups.csv")
    target_shap = shap_frame.loc[shap_frame.contribution_group == f"disease:{disease}"].copy()
    target_shap = target_shap.sort_values("mean_abs_SHAP_normalized", ascending=False).reset_index(drop=True)
    target_shap["rank"] = np.arange(1, len(target_shap) + 1)
    target_shap = target_shap.merge(gene_annotation, on="gene_id", how="left", validate="one_to_one")
    target_shap["is_original_38_gene"] = target_shap.gene_id.isin(core38)
    for ptg_name, ptg in ptg_sets.items():
        target_shap[f"in_{ptg_name}"] = target_shap.gene_id.isin(ptg)
    ranking_write_started = time.perf_counter()
    csv_write(target_shap, TABLES / f"{disease}_complete_gene_ranking.csv")
    for cutoff in (20, 50, 100):
        csv_write(target_shap.head(cutoff), TABLES / f"{disease}_top{cutoff}_preview.csv")
    ranking_write_s = time.perf_counter() - ranking_write_started

    # Local prediction-probability gradients for all cells of the target disease.
    target_rows = manifest.loc[manifest.Disease == disease, "GlobalRow"].to_numpy(np.int64)
    gradient_s, gradient_write_s = gradient_summary(disease, model.module, target_rows, out_dir)
    grad = pd.read_csv(out_dir / "causal_branch_gradient_direction.csv")
    grad = grad.merge(target_shap[["gene_id", "rank", "is_original_38_gene"] +
                                  [f"in_{name}" for name in ptg_sets]], on="gene_id", how="left")
    relation_write_started = time.perf_counter()
    csv_write(grad, TABLES / f"{disease}_candidate_gradient_relation.csv")
    relation_write_s = time.perf_counter() - relation_write_started

    eval_row = next(row for row in metrics if row["split"] == "final_eval")
    low_coverage = int(manifest.Disease.eq(disease).sum()) < 1000
    no_discrimination = (not np.isfinite(eval_row["AUROC"]) or eval_row["AUROC"] <= 0.5 or
                         eval_row["AUPRC"] <= eval_row["AUPRC_positive_baseline"])
    if low_coverage:
        status = "待解释/不支持正式候选：目标疾病细胞数低于1000，估计不稳定"
    elif no_discrimination:
        status = "待解释/不支持正式候选：最终评估未优于基线"
    else:
        status = "完成；模型有区分信号，候选仍需下游审阅"
    write_started = time.perf_counter()
    torch.save({"state_dict": {k: v.detach().cpu() for k, v in model.module.state_dict().items()},
                "config": config, "best_validation": best}, out_dir / "selected_model.pt")
    config["best_validation"] = best
    config["status"] = status
    (out_dir / "model_config.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")
    writing_s = (time.perf_counter() - write_started + checkpoint_write_s + training_log_write_s +
                 prediction_write_s + shap_write_s + ranking_write_s + gradient_write_s + relation_write_s)
    elapsed = time.perf_counter() - started_model
    runtime = {
        "Disease": disease,
        "n_cells_total": int(manifest.Disease.eq(disease).sum()),
        "n_samples": int(coverage_table.loc[coverage_table.Disease == disease, "n_samples"].iloc[0]),
        "n_datasets": int(coverage_table.loc[coverage_table.Disease == disease, "n_datasets"].iloc[0]),
        "n_Fig1_DEGs_in_union": int(fig1_deg_source.loc[fig1_deg_source.Disease == disease, "gene"].nunique()),
        "n_train": n_train,
        "n_train_positive": n_positive,
        "auto_epochs": len(history),
        "best_epoch": best["epoch"],
        "training_s": training_s,
        "validation_prediction_s": validation_elapsed_s,
        "checkpoint_write_s": checkpoint_write_s,
        "prediction_s": prediction_s,
        "shap_s": shap_s,
        "shap_write_s": shap_write_s,
        "gradient_s": gradient_s,
        "gradient_write_s": gradient_write_s,
        "write_s": writing_s,
        "total_model_s": elapsed,
        "available_memory_gb_after": available_memory_gb(),
        "process_rss_gb_after": process_rss_gb(),
        "status": status,
    }
    log("model_complete", **runtime, shap_group_cells=shap_groups)
    (out_dir / "runtime.json").write_text(json.dumps(runtime, ensure_ascii=False, indent=2), encoding="utf-8")
    del model, adata_train, checkpoint, history, shap_frame, target_shap, grad
    gc.collect()
    return {**runtime, "best_val_AUROC": best["auroc"], "best_val_AUPRC": best["auprc"],
            "holdout_AUROC": eval_row["AUROC"], "holdout_AUPRC": eval_row["AUPRC"],
            "holdout_AUPRC_baseline": eval_row["AUPRC_positive_baseline"]}


def main():
    core38_file = Path(r"E:\micro_output\分析审阅中心\待审阅\Fig4_38核心网络收尾\fig4C_network\source_data\Fig4C_38candidate_network_nodes_source.csv")
    core38_df = pd.read_csv(core38_file)
    core38 = set(core38_df.gene_id.astype(str))
    if len(core38) != 38:
        raise RuntimeError(f"冻结 38 基因网络节点数异常: {len(core38)}")
    genes_in_space = set(gene_manifest_all.gene_id.astype(str))
    genes_in_training = set(genes)
    csv_write(pd.DataFrame({
        "gene_id": sorted(core38),
        "is_in_Fig1_union_training_space": [gene in genes_in_space for gene in sorted(core38)],
        "is_in_model_input_after_variance_filter": [gene in genes_in_training for gene in sorted(core38)],
    }), INPUT / "original_38_gene_coverage.csv")
    ptg_sets = {}
    ptg_root = Path(r"E:\micro_output\fig3\trajectory_S4\tables\default\module_gene_lists")
    for i in range(1, 5):
        name = f"PTG0{i}"
        frame = pd.read_csv(ptg_root / f"{name}_genes.csv")
        col = "gene_id" if "gene_id" in frame else "gene_short_name"
        ptg_sets[name] = set(frame[col].dropna().astype(str))
    source_log = pd.read_csv(INPUT / "input_sources.csv")
    source_log = pd.concat([source_log, pd.DataFrame([
        {"role": "frozen_38_gene_set", "path": str(core38_file),
         "record": "38 frozen nodes read from the existing network source; no model/network/drug results modified"},
        {"role": "frozen_PTG_programs", "path": str(ptg_root),
         "record": "PTG01-PTG04 membership read from existing trajectory gene lists"},
    ])], ignore_index=True)
    csv_write(source_log, INPUT / "input_sources.csv")
    genes_file = INPUT / "gene_manifest.csv"
    gene_frame = pd.read_csv(genes_file)
    gene_frame["is_original_38_gene"] = gene_frame.gene_id.isin(core38)
    for name, members in ptg_sets.items():
        gene_frame[f"in_{name}"] = gene_frame.gene_id.isin(members)
    csv_write(gene_frame, INPUT / "gene_manifest_with_frozen_program_links.csv")
    source_rows = []
    for _, row in gene_frame.iterrows():
        source_rows.append({"gene_id": row.gene_id,
                            "is_original_38_gene": row.is_original_38_gene,
                            **{f"in_{name}": bool(row[f"in_{name}"]) for name in ptg_sets}})
    links = pd.DataFrame(source_rows)
    status_path = TABLES / "disease_model_status.csv"
    if status_path.exists():
        old_status = pd.read_csv(status_path)
        results = old_status.to_dict(orient="records")
        completed = {str(row["Disease"]) for row in results
                     if str(row.get("status", "")).startswith(("完成", "待解释/不支持正式候选"))}
    else:
        results = []
        completed = set()
    requested = os.environ.get("OVR_DISEASES")
    if requested:
        diseases_to_run = [d.strip() for d in requested.split(",") if d.strip()]
        unknown = sorted(set(diseases_to_run) - set(diseases))
        if unknown:
            raise RuntimeError(f"请求的疾病不在输入清单中: {unknown}")
    else:
        diseases_to_run = diseases
    for disease in diseases_to_run:
        if disease in completed:
            log("model_already_completed_skip", disease=disease)
            continue
        try:
            results.append(run_one(disease, core38, ptg_sets))
        except Exception as exc:
            failure = {
                "Disease": disease,
                "n_cells_total": int(manifest.Disease.eq(disease).sum()),
                "status": "失败",
                "error": str(exc),
                "traceback": traceback.format_exc(),
                "available_memory_gb_at_failure": available_memory_gb(),
                "process_rss_gb_at_failure": process_rss_gb(),
            }
            failure_path = MODELS / disease / "failure.json"
            failure_path.parent.mkdir(parents=True, exist_ok=True)
            failure_path.write_text(json.dumps(failure, ensure_ascii=False, indent=2), encoding="utf-8")
            log("model_failed", **failure)
            results.append(failure)
        csv_write(pd.DataFrame(results), TABLES / "disease_model_status.csv")
        # Update aggregate candidate-to-frozen-program relation after each completed model.
        if results:
            ranking_parts = []
            for done in results:
                path = TABLES / f"{done['Disease']}_complete_gene_ranking.csv"
                if path.exists():
                    ranking_parts.append(pd.read_csv(path))
            if ranking_parts:
                csv_write(pd.concat(ranking_parts, ignore_index=True),
                          TABLES / "candidate_to_38_and_PTG_relationships.csv")
    csv_write(pd.DataFrame(results), TABLES / "disease_model_status.csv")
    log("all_models_complete", n_models=len(results), available_memory_gb=available_memory_gb())




def atomic_save(payload: dict, path: Path) -> None:
    temp_path = path.with_suffix(path.suffix + ".tmp")
    torch.save(payload, temp_path)
    os.replace(temp_path, path)


def _binary_rows(disease: str, split: str) -> pd.DataFrame:
    frame = manifest.loc[(manifest.Split == split) & manifest.Disease.isin([disease, "Ctrl"])]
    frame = frame.sort_values("GlobalRow", kind="stable").reset_index(drop=True)
    if frame.empty or set(frame.Disease.astype(str).unique()) != {disease, "Ctrl"}:
        raise RuntimeError(f"{disease}/{split}: 目标病或 Ctrl 缺失")
    if frame.Disease.isin(set(diseases) - {disease}).any():
        raise RuntimeError(f"{disease}/{split}: 其他疾病错误进入模型")
    return frame


def _loss_schedule(module: torch.nn.Module, n_epochs: int) -> tuple[pd.DataFrame, list[dict]]:
    rows = []
    for epoch0 in range(n_epochs):
        weights = module.update_loss_weights(epoch0, n_epochs, scheme=None)
        rows.append({"epoch": epoch0 + 1, "epoch_index_zero_based": epoch0,
                     **{f"weight_{k}": float(v) for k, v in weights.items()}})
    frame = pd.DataFrame(rows)
    weight_cols = [c for c in frame.columns if c.startswith("weight_")]
    group = frame[weight_cols].round(12).astype(str).agg("|".join, axis=1).ne(
        frame[weight_cols].round(12).astype(str).agg("|".join, axis=1).shift()).cumsum()
    stages = []
    for stage_no, (_, part) in enumerate(frame.groupby(group, sort=True), start=1):
        stages.append({
            "stage": stage_no, "epoch_start_1based": int(part.epoch.min()),
            "epoch_end_1based": int(part.epoch.max()), "n_epochs": int(len(part)),
            "loss_weights": {c.removeprefix("weight_"): float(part.iloc[0][c]) for c in weight_cols},
        })
    return frame, stages


def run_one_fig6(disease: str) -> dict:
    if disease not in diseases:
        raise ValueError(f"Disease不在当前对象: {disease}")
    out = MODELS / disease
    out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    split_frames = {s: _binary_rows(disease, s) for s in ("train", "validation", "final_eval")}
    train_frame = split_frames["train"]
    train_rows = train_frame.GlobalRow.to_numpy(np.int64)
    train_diseases = train_frame.Disease.astype(str).to_numpy()
    y_train = (train_diseases == disease).astype(np.int8)
    n_train, n_pos = int(len(train_rows)), int(y_train.sum())
    n_neg = n_train - n_pos
    if n_pos == 0 or n_neg == 0:
        raise RuntimeError(f"{disease}: train缺少一类")
    weights = [n_train / (2*n_neg), n_train / (2*n_pos)]
    n_epochs = min(round((20000 / n_train) * 400), 400)
    if n_epochs < 1:
        raise RuntimeError(f"{disease}: 自动轮数小于1")
    val_rows = split_frames["validation"].GlobalRow.to_numpy(np.int64)
    val_y = (split_frames["validation"].Disease.to_numpy(str) == disease).astype(np.int8)

    for split, frame in split_frames.items():
        csv_write(frame, out / f"model_{split}_cells.csv.gz")
    set_seed(SEED)
    adata = train_adata_for(disease)
    model = CausalFinder(adata=adata, **MODEL_KWARGS)
    stage_table, stages = _loss_schedule(model.module, n_epochs)
    csv_write(stage_table, out / "training_stage_schedule.csv")
    coverage = coverage_table.set_index("Disease")
    audit = pd.read_csv(TABLES / "input_and_method_audit.csv").set_index("item")["value"].to_dict()
    config = {
        "Disease": disease, "positive_class": disease, "negative_class": "Ctrl",
        "other_diseases_enter_model": False, "input_object_sha256": audit["object_sha256"],
        "seed": SEED, "split_unit": "cell; not independent-patient validation",
        "split": {"train": 0.64, "validation": 0.16, "final_eval": 0.20, "shared_ctrl_split": True},
        "split_counts": {s: {"n": len(f), "positive_target": int(f.Disease.eq(disease).sum()),
                             "negative_Ctrl": int(f.Disease.eq("Ctrl").sum())} for s, f in split_frames.items()},
        "coverage_target": coverage.loc[disease].to_dict(), "coverage_Ctrl": coverage.loc["Ctrl"].to_dict(),
        "n_train": n_train, "n_positive_train": n_pos, "n_negative_train": n_neg,
        "class_weights": {"negative_Ctrl": weights[0], "positive_target": weights[1]},
        "class_weight_formula": "N_train/(2*N_class_train)",
        "class_weight_use": "one multiplier in per-cell BCE, then batch mean",
        "n_model_genes": len(genes), "n_all_RNA_genes": len(gene_manifest_all),
        "model_kwargs": MODEL_KWARGS,
        "training": {"max_epochs": n_epochs,
          "epoch_rule": "min(round((20000/n_train_cells)*400),400), Python round",
          "learning_rate": 0.0005, "weight_decay": 0.000001, "optimizer": "Adam",
          "batch_size": BATCH_SIZE, "n_x": 5, "n_alpha": N_ALPHA, "n_beta": N_BETA,
          "stage_training": True, "spurious_fold": 2.0, "stage_schedule": stages},
        "checkpoint_selection": "validation AUROC; AUPRC tie-break",
        "final_eval_used_for_selection": False, "threshold": 0.5,
        "checkpoint_files": ["best_validation_checkpoint.pt", "end_of_training_checkpoint.pt",
                             "latest_training_checkpoint.pt"],
        "attribution": {"method": "Deep SHAP; original CauFinder forward; check_additivity=False as audited",
          "background_target_max": 50, "background_Ctrl_max": 50, "background_sampling": "without replacement",
          "seed": SEED, "batch_size": SHAP_BATCH_SIZE,
          "explained_cells": "all target-disease and Ctrl training cells; saved in exact explanation order"},
        "resumable_fields": ["model_state_dict", "optimizer_state_dict", "epoch_completed",
          "python_rng_state", "numpy_rng_state", "torch_rng_state", "cuda_rng_state_all", "epoch_losses"],
        "gpu_available": torch.cuda.is_available(),
    }
    config_fingerprint = hashlib.sha256(json.dumps(config, sort_keys=True, ensure_ascii=False,
                                                    separators=(",", ":")).encode("utf-8")).hexdigest()
    config["configuration_fingerprint"] = config_fingerprint
    (out / "model_config.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")
    best_path, latest_path = out/"best_validation_checkpoint.pt", out/"latest_training_checkpoint.pt"
    epoch_log = out/"epoch_training.jsonl"
    best = {"auroc": -np.inf, "auprc": -np.inf, "epoch": None}
    resume_state = None
    if latest_path.exists():
        resume_state = torch.load(latest_path, map_location="cpu", weights_only=False)
        if resume_state.get("configuration_fingerprint") != config_fingerprint:
            raise RuntimeError(f"{disease}: latest checkpoint与冻结配置不一致")
        model.module.load_state_dict(resume_state["model_state_dict"])
        best = resume_state["best_validation"]
    elif epoch_log.exists() and epoch_log.stat().st_size:
        raise RuntimeError(f"{disease}: 存在epoch日志但没有可恢复checkpoint")
    val_seconds, checkpoint_seconds = 0.0, 0.0

    def validation_callback(module, epoch, epoch_losses):
        nonlocal best, val_seconds, checkpoint_seconds
        t0 = time.perf_counter()
        p = predict_prob(module, val_rows, INFERENCE_BATCH_SIZE)
        val_seconds += time.perf_counter() - t0
        auroc, auprc = auc_metrics(val_y, p)
        if not np.isfinite(auroc) or not np.isfinite(auprc):
            raise RuntimeError(f"{disease}: validation缺少正负类")
        improved = (auroc > best["auroc"]) or (np.isclose(auroc, best["auroc"]) and auprc > best["auprc"])
        if improved:
            best = {"auroc": float(auroc), "auprc": float(auprc), "epoch": int(epoch + 1)}
            checkpoint_started = time.perf_counter()
            atomic_save({"state_dict": {k:v.detach().cpu().clone() for k,v in module.state_dict().items()},
                         "best_validation": best, "configuration_fingerprint": config_fingerprint,
                         "config": config}, best_path)
            checkpoint_seconds += time.perf_counter() - checkpoint_started
        logrow = {"Disease": disease, "epoch": epoch+1, "val_AUROC": float(auroc),
                  "val_AUPRC": float(auprc), "is_best_checkpoint": bool(improved),
                  "best_epoch_so_far": best["epoch"], "best_val_AUROC_so_far": best["auroc"],
                  "best_val_AUPRC_so_far": best["auprc"]}
        for key, vals in epoch_losses.items():
            if vals and key not in logrow: logrow[key] = vals[-1]
        with epoch_log.open("a", encoding="utf-8") as f:
            f.write(json.dumps(logrow, ensure_ascii=False, default=float)+"\n")
            f.flush()
        log("validation_epoch", **logrow, available_memory_gb=available_memory_gb(),
            process_rss_gb=process_rss_gb())
        return {"val_AUROC": float(auroc), "val_AUPRC": float(auprc)}

    def state_callback(module, optimizer, epoch_completed, max_epochs_arg, history):
        payload = {
            "Disease": disease, "configuration_fingerprint": config_fingerprint,
            "max_epochs": int(max_epochs_arg), "epoch_completed": int(epoch_completed),
            "model_state_dict": {k:v.detach().cpu().clone() for k,v in module.state_dict().items()},
            "optimizer_state_dict": optimizer.state_dict(),
            "python_rng_state": random.getstate(), "numpy_rng_state": np.random.get_state(),
            "torch_rng_state": torch.get_rng_state(),
            "cuda_rng_state_all": torch.cuda.get_rng_state_all() if torch.cuda.is_available() else None,
            "epoch_losses": copy.deepcopy(history), "best_validation": copy.deepcopy(best),
        }
        atomic_save(payload, latest_path)

    train_started = time.perf_counter()
    completed = int(resume_state["epoch_completed"]) if resume_state is not None else 0
    if completed < n_epochs:
        model.train(max_epochs=n_epochs, lr=5e-4, weight_decay=1e-6, n_x=5, causal_weight=1.0,
          spurious_fold=2.0, stage_training=True, train_size=1.0, validation_size=None,
          batch_size=BATCH_SIZE, drop_last=False, n_alpha=N_ALPHA, n_beta=N_BETA,
          class_weights=weights, epoch_callback=validation_callback,
          training_state_callback=state_callback, resume_state=resume_state)
    else:
        model.history = resume_state["epoch_losses"]
    train_seconds = time.perf_counter()-train_started
    if not best_path.is_file() or not latest_path.is_file():
        raise RuntimeError(f"{disease}: best/end checkpoint缺失")
    end_path = out/"end_of_training_checkpoint.pt"
    if end_path.exists(): end_path.unlink()
    try:
        os.link(latest_path, end_path)
    except OSError:
        shutil.copy2(latest_path, end_path)
    history, history_write_s = stage_timing(model.history, disease, out)
    best_ckpt = torch.load(best_path, map_location="cpu", weights_only=False)
    model.module.load_state_dict(best_ckpt["state_dict"])
    model.module.eval()

    preds, metrics = [], []
    t_pred = time.perf_counter()
    for split, frame in split_frames.items():
        rows = frame.GlobalRow.to_numpy(np.int64)
        p = predict_prob(model.module, rows, INFERENCE_BATCH_SIZE)
        metric, pred, _ = evaluate_split(disease, split, rows, p)
        preds.append(pred); metrics.append(metric)
    prediction_seconds = time.perf_counter()-t_pred
    prediction_frame = pd.concat(preds, ignore_index=True)
    csv_write(prediction_frame, TABLES/f"{disease}_cell_predictions.csv.gz")
    csv_write(pd.DataFrame(metrics), TABLES/f"{disease}_evaluation_metrics.csv")

    raw_gate = model.module.feature_mapper.weight.detach().cpu().numpy().reshape(-1)
    gate_threshold = float(model.module.feature_mapper.thresh.detach().cpu())
    gate_probability = 1.0 / (1.0 + np.exp(-raw_gate))
    effective_gate = gate_probability > gate_threshold
    prediction_ranges = prediction_frame.groupby("Split", sort=False).p_target.agg(
        n="size", minimum="min", maximum="max", mean="mean", median="median").reset_index()
    csv_write(prediction_ranges, out/"alpha_prediction_range_by_split.csv")
    gate_summary = {
        "Disease": disease, "n_genes": int(len(raw_gate)),
        "effective_gate_count": int(effective_gate.sum()),
        "effective_gate_fraction": float(effective_gate.mean()),
        "gate_threshold": gate_threshold, "comparison": "sigmoid(weight) > threshold",
        "sigmoid_gate_min": float(gate_probability.min()),
        "sigmoid_gate_median": float(np.median(gate_probability)),
        "sigmoid_gate_max": float(gate_probability.max()),
        "prediction_min": float(prediction_frame.p_target.min()),
        "prediction_max": float(prediction_frame.p_target.max()),
        "prediction_split_ranges": prediction_ranges.to_dict(orient="records"),
        "prediction_branch": "alpha_dpd probability from original CauFinder forward",
    }
    (out/"gate_alpha_prediction_summary.json").write_text(
        json.dumps(gate_summary, ensure_ascii=False, indent=2), encoding="utf-8")

    shap_counts, shap_seconds, shap_write_seconds = shap_summaries(
        disease, model.module, train_rows, y_train, train_diseases, out
    )
    shap_all = pd.read_csv(out/"shap_importance_all_groups.csv")
    for group, label, name in (
        (f"disease:{disease}", "target_disease_train", f"{disease}_SHAP_target_full_ranking.csv"),
        ("disease:Ctrl", "Ctrl_train", f"{disease}_SHAP_Ctrl_full_ranking.csv"),
    ):
        ranking = shap_all.loc[shap_all.contribution_group == group].copy()
        ranking = ranking.sort_values(["mean_abs_SHAP","gene_id"],ascending=[False,True]).reset_index(drop=True)
        if len(ranking)!=len(genes) or set(ranking.gene_id.astype(str))!=set(genes):
            raise RuntimeError(f"{disease}/{label}: SHAP未完整覆盖全部模型基因")
        ranking["rank"] = np.arange(1,len(ranking)+1)
        ranking["analysis_group"] = label
        ranking = ranking.merge(gene_annotation,on="gene_id",how="left",validate="one_to_one",suffixes=("","_manifest"))
        csv_write(ranking,TABLES/name)

    final = next(row for row in metrics if row["split"]=="final_eval")
    n_target = int(coverage.loc[disease,"n_cells"])
    low = n_pos < 200 or n_target < 1000
    weak = (not np.isfinite(final["AUROC"]) or not np.isfinite(final["AUPRC"]) or
            final["AUROC"]<=0.5 or final["AUPRC"]<=final["AUPRC_positive_baseline"])
    status = ("calculated; interpretation limited by low target-cell count" if low else
              "calculated; final evaluation shows no clear discrimination" if weak else
              "calculated; discrimination signal present; candidate interpretation remains descriptive")
    best_stage = next((s["stage"] for s in stages if s["epoch_start_1based"]<=best["epoch"]<=s["epoch_end_1based"]),None)
    runtime = {
        "Disease":disease,"n_cells_target_total":n_target,
        "n_cells_Ctrl_total":int(coverage.loc["Ctrl","n_cells"]),
        "n_train":n_train,"n_positive_train":n_pos,"n_negative_train":n_neg,
        "n_features":len(genes),"auto_epochs":n_epochs,"best_epoch":best["epoch"],"best_stage":best_stage,
        "best_val_AUROC":best["auroc"],"best_val_AUPRC":best["auprc"],
        "effective_causal_gate_count":gate_summary["effective_gate_count"],
        "effective_causal_gate_fraction":gate_summary["effective_gate_fraction"],
        "alpha_prediction_min":gate_summary["prediction_min"],
        "alpha_prediction_max":gate_summary["prediction_max"],
        "final_eval_AUROC":final["AUROC"],"final_eval_AUPRC":final["AUPRC"],
        "final_eval_AUPRC_positive_baseline":final["AUPRC_positive_baseline"],
        "n_target_train_SHAP_cells":int(shap_counts[f"disease:{disease}"]),
        "n_Ctrl_train_SHAP_cells":int(shap_counts["disease:Ctrl"]),
        "training_seconds":train_seconds,"validation_seconds":val_seconds,
        "checkpoint_seconds":checkpoint_seconds,"prediction_seconds":prediction_seconds,
        "SHAP_seconds":shap_seconds,"SHAP_write_seconds":shap_write_seconds,
        "history_write_seconds":history_write_s,"total_model_seconds":time.perf_counter()-started,
        "available_memory_gb_after":available_memory_gb(),"process_rss_gb_after":process_rss_gb(),
        "status":status}
    (out/"runtime.json").write_text(json.dumps(runtime,ensure_ascii=False,indent=2),encoding="utf-8")
    log("model_SHAP_complete",**runtime)
    return runtime



def write_fig6_batch_plan() -> None:
    config_dir = ROOT / "config"
    config_dir.mkdir(parents=True, exist_ok=True)
    audit = pd.read_csv(TABLES / "input_and_method_audit.csv").set_index("item")["value"].to_dict()
    dummy = ad.AnnData(X=sp.csr_matrix((1, 1), dtype=np.float32),
                       obs=pd.DataFrame({"labels": [0.0]}, index=["config_only"]),
                       var=pd.DataFrame(index=["config_only"]))
    set_seed(SEED)
    schedule_model = CausalFinder(adata=dummy, **MODEL_KWARGS).module
    counts = manifest.groupby(["Disease", "Split"]).size().to_dict()
    plan_rows, schedule_rows = [], []
    plan = {
        "input_object_sha256": audit["object_sha256"],
        "seed": SEED, "current_diseases": diseases,
        "ctrl_split_generated_once_and_reused": True,
        "split": {"train": .64, "validation": .16, "final_eval": .20, "unit": "cell"},
        "model_kwargs": MODEL_KWARGS, "batch_size": BATCH_SIZE, "learning_rate": .0005,
        "optimizer": "Adam", "weight_decay": .000001, "n_x": 5,
        "n_alpha": N_ALPHA, "n_beta": N_BETA, "stage_training": True,
        "epoch_rule": "min(round((20000/n_train_cells)*400),400), Python round; per disease+Ctrl",
        "class_weight_rule": "N_train/(2*N_class_train), applied once to per-cell BCE before batch mean",
        "disease_configurations": {},
    }
    for disease in diseases:
        n_pos = int(counts[(disease, "train")])
        n_neg = int(counts[("Ctrl", "train")])
        n_train = n_pos + n_neg
        n_epochs = min(round((20000 / n_train) * 400), 400)
        _, stages = _loss_schedule(schedule_model, n_epochs)
        record = {
            "Disease": disease, "n_train": n_train, "n_positive_train": n_pos,
            "n_negative_ctrl_train": n_neg,
            "w_positive": n_train/(2*n_pos), "w_negative": n_train/(2*n_neg),
            "max_epochs": n_epochs, "batch_size": BATCH_SIZE, "learning_rate": .0005,
            "optimizer": "Adam", "weight_decay": .000001,
            "n_x": 5, "n_alpha": N_ALPHA, "n_beta": N_BETA,
            "stage_schedule": stages,
        }
        plan["disease_configurations"][disease] = record
        plan_rows.append({k: v for k, v in record.items() if k != "stage_schedule"} |
                         {"stage_schedule_json": json.dumps(stages, ensure_ascii=False)})
        for stage in stages:
            schedule_rows.append({"Disease": disease, **stage,
                                  "loss_weights_json": json.dumps(stage["loss_weights"], ensure_ascii=False)})
    (config_dir / "all_disease_training_configurations.json").write_text(
        json.dumps(plan, ensure_ascii=False, indent=2), encoding="utf-8")
    csv_write(pd.DataFrame(plan_rows), config_dir / "all_disease_training_configurations.csv")
    csv_write(pd.DataFrame(schedule_rows), config_dir / "all_disease_stage_boundaries.csv")

def main_fig6() -> None:
    if len(sys.argv)==2 and sys.argv[1]=="--plan-only":
        write_fig6_batch_plan()
        log("all_disease_training_plan_written", n_diseases=len(diseases))
        return
    if len(sys.argv)!=3 or sys.argv[1]!="--disease":
        raise SystemExit("usage: Fig6_train_disease_vs_ctrl.py --plan-only | --disease <Disease>")
    disease=sys.argv[2]
    try:
        write_fig6_batch_plan()
        result=run_one_fig6(disease)
    except Exception as exc:
        failure={"Disease":disease,"status":"failed","error":str(exc),
                 "traceback":traceback.format_exc(),"available_memory_gb":available_memory_gb(),
                 "process_rss_gb":process_rss_gb()}
        path=MODELS/disease/"failure.json"; path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(failure,ensure_ascii=False,indent=2),encoding="utf-8")
        log("model_failed",**failure)
        raise
    status_path=TABLES/"disease_model_status.csv"
    prior=pd.read_csv(status_path) if status_path.exists() else pd.DataFrame()
    if "Disease" in prior.columns:
        prior=prior.loc[prior["Disease"].astype(str)!=disease]
    elif len(prior):
        raise RuntimeError(f"{status_path}存在记录但缺少Disease列")
    csv_write(pd.concat([prior,pd.DataFrame([result])],ignore_index=True),status_path)

if __name__ == "__main__":
    try:
        main_fig6()
    except Exception as exc:
        log("fatal_error", error=str(exc), traceback=traceback.format_exc(),
            available_memory_gb=available_memory_gb(), process_rss_gb=process_rss_gb())
        raise
