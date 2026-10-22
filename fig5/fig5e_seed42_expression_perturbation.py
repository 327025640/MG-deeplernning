#!/usr/bin/env python3
"""Replay frozen PTG14 holdout perturbations for the 76 seed42 PCSF model-evidence genes."""
from __future__ import annotations

import ctypes
import hashlib
import json
import os
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

MODEL_ROOT = Path(r"E:\micro_output\分析审阅中心\待审阅\PTG02_CD_分支疾病与调控候选\ptg14_training_v2")
DATA_DIR = MODEL_ROOT / "data"
RUNS_DIR = MODEL_ROOT / "runs"
ADAPTER_DIR = Path(r"E:\micro_output\分析审阅中心\待审阅\PTG02_CD_分支疾病与调控候选\data\pan14_disease\caufinder_adapter")
OUT_ROOT = Path(r"E:\micro_output\分析审阅中心\待审阅\Fig5_pan14_PPI_PCSF更新")
OLD_ROOT = Path(r"E:\micro_output\分析审阅中心\待审阅\Fig5_主分析与流程图收尾")
SOURCE_DIR = OUT_ROOT / "source_data" / "current_seed42_pcsf"
FIGURE_DIR = OUT_ROOT / "figures" / "current_seed42_pcsf"
LOG_DIR = OUT_ROOT / "logs" / "current_seed42_pcsf"
SEEDS = tuple(range(42, 52))
GROUPS = ("Neurological", "Tumor")
SCENARIOS = (("down50", 0.5), ("up50", 1.5))
BATCH_SIZE = 128
EXPECTED_FEATURES = 4843
EXPECTED_CELLS = 365563
EXPECTED_HOLDOUT = 61355
CONTRACT_ID = "a02f5b0437435e2029a8b1a4ab2f3896f275208f97804983c3fbdeff8d20a597"
PTG_COLORS = {"PTG01": "#2166AC", "PTG02": "#B2182B", "PTG03": "#1B9E77", "PTG04": "#7570B3"}
CORE_RED = "#D73027"
BLUE = "#2166AC"
WHITE = "#F7F7F7"

OUT_DELTAS = SOURCE_DIR / "Fig5E_seed42_PCSF_paired_deltas_by_seed.csv"
OUT_BASELINE = SOURCE_DIR / "Fig5E_seed42_PCSF_baseline_by_seed.csv"
OUT_HEATMAP = SOURCE_DIR / "Fig5E_seed42_PCSF_heatmap_values.csv"
OUT_PARAMS = SOURCE_DIR / "Fig5E_seed42_PCSF_display_parameters.csv"
OUT_PROVENANCE = SOURCE_DIR / "Fig5E_seed42_PCSF_provenance.json"
OUT_FAILURE = LOG_DIR / "fig5e_failure.json"
OUT_PNG = FIGURE_DIR / "Fig5E_seed42_PCSF_expression_perturbation.png"
OUT_PDF = FIGURE_DIR / "Fig5E_seed42_PCSF_expression_perturbation.pdf"
OUT_LOG = LOG_DIR / "fig5e_pan14_expression_perturbation.log"
PROGRESS_DIR = LOG_DIR / "fig5e_progress"
PROGRESS_STATE = PROGRESS_DIR / "identity.json"
PROGRESS_BASELINE = PROGRESS_DIR / "baseline_by_seed.csv"
PROGRESS_DELTAS = PROGRESS_DIR / "paired_deltas_by_seed.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_model_code_fingerprints(contract: dict) -> dict[str, str]:
    fingerprints = contract.get("fingerprints", {})
    code_entries = {str(raw_path): str(expected) for raw_path, expected in fingerprints.items()
                    if str(raw_path).lower().endswith(".py")}
    if not code_entries:
        raise ValueError("PTG14 contract has no Python source fingerprints for frozen model code")
    observed = {}
    for raw_path, expected in code_entries.items():
        path = Path(raw_path)
        if not path.is_file():
            raise FileNotFoundError("Frozen model code listed by contract is unavailable: {}".format(path))
        actual = sha256(path)
        if actual != expected:
            raise ValueError("Frozen model code SHA differs from contract: {}".format(path))
        observed[str(path)] = actual
    return observed


def log(message: str) -> None:
    line = "[{}] {}".format(datetime.now().astimezone().isoformat(timespec="seconds"), message)
    print(line, flush=True)
    OUT_LOG.parent.mkdir(parents=True, exist_ok=True)
    with OUT_LOG.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    frame.to_csv(temporary, index=False)
    os.replace(temporary, path)


def atomic_text(text: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def atomic_figure_save(fig, path: Path, file_format: str, dpi: int | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.stem + ".partial" + path.suffix)
    kwargs = {"format": file_format, "facecolor": "white"}
    if dpi is not None:
        kwargs["dpi"] = dpi
    fig.savefig(temporary, **kwargs)
    os.replace(temporary, path)


def read_csv_compat(path: Path, **kwargs) -> pd.DataFrame:
    """Read Unicode Windows paths through a Python handle for the frozen pandas build."""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return pd.read_csv(handle, **kwargs)


def load_progress_csv(path: Path, keys: list[str], required: set[str]) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    frame = read_csv_compat(path)
    if not required.issubset(frame.columns) or frame.duplicated(keys).any():
        raise ValueError("Persisted Fig5E progress file has invalid columns or duplicate keys: {}".format(path))
    return frame


def memory_status() -> dict:
    class MEMORYSTATUSEX(ctypes.Structure):
        _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
    status = MEMORYSTATUSEX(); status.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
        raise OSError("GlobalMemoryStatusEx failed")
    return {"available_ram_gib": status.ullAvailPhys / 2**30,
            "total_ram_gib": status.ullTotalPhys / 2**30,
            "available_pagefile_gib": status.ullAvailPageFile / 2**30}


def read_inputs():
    registry_path = SOURCE_DIR / "seed42_model_evidence_76_registry.csv"
    nodes_path = SOURCE_DIR / "seed42_C_nodes_82.csv"
    if not registry_path.is_file() or not nodes_path.is_file():
        raise FileNotFoundError("Required prepared seed42 evidence/node source is missing: {} or {}".format(registry_path, nodes_path))
    registry = read_csv_compat(registry_path, dtype={"gene_id": str, "PTG": str})
    nodes = read_csv_compat(nodes_path, dtype={"gene_id": str, "Entrez": str})
    required = {"gene_id", "Entrez", "PTG", "selected_runs", "core_10_of_10", "selected_rank"}
    if not required.issubset(registry.columns) or len(registry) != 76 or registry.gene_id.duplicated().any() or registry.Entrez.astype(str).duplicated().any():
        raise ValueError("Seed42 model-evidence registry fields/uniqueness do not satisfy the frozen PCSF contract")
    if not {"gene_id", "Entrez", "model_evaluable"}.issubset(nodes.columns) or len(nodes) != 82 or nodes.Entrez.astype(str).duplicated().any():
        raise ValueError("Seed42 selected-node source does not match the frozen 82-node network identity")
    evidence_nodes = nodes[nodes.model_evaluable.astype(str).str.lower().isin(["true", "1"])].copy()
    evidence_identity = registry[["gene_id", "Entrez"]].astype(str).sort_values(["gene_id", "Entrez"]).reset_index(drop=True)
    node_identity = evidence_nodes[["gene_id", "Entrez"]].astype(str).sort_values(["gene_id", "Entrez"]).reset_index(drop=True)
    if len(evidence_nodes) != 76 or not evidence_identity.equals(node_identity):
        raise ValueError("E perturbation candidates are not exactly the 76 model-evaluable nodes from seed42 PCSF")
    cells = read_csv_compat(DATA_DIR / "cells.csv", dtype=str)
    genes = read_csv_compat(DATA_DIR / "genes.csv", dtype=str)
    contract = json.loads((DATA_DIR / "contract.json").read_text(encoding="utf-8"))
    export = json.loads((DATA_DIR / "export_complete.json").read_text(encoding="utf-8"))
    fingerprints = contract.get("fingerprints", {})
    for path in (DATA_DIR / "cells.csv", DATA_DIR / "genes.csv", DATA_DIR / "export_complete.json"):
        expected = fingerprints.get(str(path))
        if not expected:
            raise ValueError("Frozen contract lacks a SHA fingerprint for {}".format(path))
        if sha256(path) != expected:
            raise ValueError("Frozen contract input SHA mismatch: {}".format(path))
    code_hashes = verify_model_code_fingerprints(contract)
    if contract.get("id") != CONTRACT_ID or export.get("n_cells") != EXPECTED_CELLS or export.get("n_genes") != EXPECTED_FEATURES:
        raise ValueError("Frozen PTG14 contract/export dimensions mismatch")
    if len(cells) != EXPECTED_CELLS or len(genes) != EXPECTED_FEATURES or genes.gene_id.duplicated().any():
        raise ValueError("Frozen cell/gene manifest dimensions or uniqueness mismatch")
    if not set(cells.Group3.unique()) >= {"Neurological", "Tumor", "Ctrl"}:
        raise ValueError("Frozen manifest labels differ from Neurological/Tumor/Ctrl reference contract")
    if set(cells.loc[cells.split == "holdout", "Group3"].unique()) != set(GROUPS):
        raise ValueError("Frozen holdout must contain only Neurological and Tumor cells")
    holdout = np.flatnonzero(np.asarray(cells.split.eq("holdout"), dtype=bool))
    if len(holdout) != EXPECTED_HOLDOUT or not np.array_equal(holdout, np.sort(holdout)):
        raise ValueError("Frozen holdout count or input row order mismatch")
    if not set(registry.gene_id.astype(str)).issubset(set(genes.gene_id.astype(str))):
        raise ValueError("A seed42 PCSF model-evidence candidate is absent from the frozen PTG14 feature set")
    feature_index = genes.set_index("gene_id").index.get_indexer(registry.gene_id.astype(str))
    if np.any(feature_index < 0):
        raise ValueError("A majority candidate cannot be indexed in the frozen feature order")
    return registry.sort_values("selected_rank").reset_index(drop=True), cells, genes, holdout, feature_index, contract, export, code_hashes


def verify_seed_runs():
    runs = {}
    for seed in SEEDS:
        run_dir = RUNS_DIR / "seed_{}".format(seed)
        complete_path = run_dir / "complete.json"
        checkpoint_path = run_dir / "model.pt"
        weight_path = run_dir / "weights.csv"
        if not all(p.is_file() for p in (complete_path, checkpoint_path, weight_path)):
            raise FileNotFoundError("seed {} model/SHAP output is incomplete".format(seed))
        complete = json.loads(complete_path.read_text(encoding="utf-8"))
        if complete.get("seed") != seed or complete.get("contract") != CONTRACT_ID:
            raise ValueError("seed {} checkpoint complete.json identity mismatch".format(seed))
        observed = sha256(checkpoint_path)
        if observed != complete.get("outputs", {}).get("model.pt"):
            raise ValueError("seed {} model.pt SHA mismatch".format(seed))
        weight_hash = sha256(weight_path)
        if weight_hash != complete.get("outputs", {}).get("weights.csv"):
            raise ValueError("seed {} weights.csv SHA mismatch".format(seed))
        runs[seed] = {"checkpoint": checkpoint_path, "complete": complete, "checkpoint_sha256": observed,
                      "weights_sha256": weight_hash}
    return runs


def load_frozen_model(seed: int, checkpoint_path: Path, genes: pd.DataFrame):
    sys.path.insert(0, str(ADAPTER_DIR))
    from CauFinder.caufinder_main import CausalFinder, ShapModel
    from CauFinder.utils import set_seed
    set_seed(seed)
    tiny = ad.AnnData(X=np.zeros((1, len(genes)), dtype=np.float32),
        obs=pd.DataFrame({"labels": np.array([0.0], dtype=np.float32)}),
        var=pd.DataFrame(index=genes.gene_id.astype(str)))
    model = CausalFinder(adata=tiny, n_latent=25, n_causal=5, n_hidden=128,
        n_layers_encoder=0, n_layers_decoder=0, n_layers_dpd=0,
        dropout_rate_encoder=0.0, dropout_rate_decoder=0.0, dropout_rate_dpd=0.0,
        use_batch_norm="none", use_batch_norm_dpd=True, pdp_linear=True)
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    if checkpoint.get("contract") != CONTRACT_ID or checkpoint.get("seed") != seed:
        raise ValueError("seed {} checkpoint identity mismatch".format(seed))
    model.module.load_state_dict(checkpoint["state_dict"], strict=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.module.to(device).eval()
    wrapped = ShapModel(model.module, "prob").to(device).eval()
    return wrapped, device


def predict_gene(wrapped, device, matrix, holdout_rows, gene_index, multiplier):
    result = np.empty(len(holdout_rows), dtype=np.float32)
    with torch.no_grad():
        for start in range(0, len(holdout_rows), BATCH_SIZE):
            stop = min(start + BATCH_SIZE, len(holdout_rows))
            row_ids = holdout_rows[start:stop]
            batch = np.array(matrix[row_ids, :], dtype=np.float32, copy=True)
            original_log = batch[:, gene_index].astype(np.float64)
            if not np.isfinite(original_log).all() or np.any(original_log < 0):
                raise ValueError("Frozen RNA input is not finite nonnegative log1p data")
            perturbed_log = np.log1p(np.expm1(original_log) * multiplier)
            if not np.isfinite(perturbed_log).all():
                raise ValueError("Linear-scale perturbation produced non-finite values")
            batch[:, gene_index] = perturbed_log.astype(np.float32)
            pred = wrapped(torch.as_tensor(batch, dtype=torch.float32, device=device)).detach().cpu().numpy().reshape(-1)
            if not np.isfinite(pred).all() or np.any(pred < 0) or np.any(pred > 1):
                raise ValueError("Frozen model returned invalid Tumor probabilities")
            result[start:stop] = pred
    return result


def predict_baseline(wrapped, device, matrix, holdout_rows):
    result = np.empty(len(holdout_rows), dtype=np.float32)
    with torch.no_grad():
        for start in range(0, len(holdout_rows), BATCH_SIZE):
            stop = min(start + BATCH_SIZE, len(holdout_rows))
            batch = np.array(matrix[holdout_rows[start:stop], :], dtype=np.float32, copy=True)
            pred = wrapped(torch.as_tensor(batch, dtype=torch.float32, device=device)).detach().cpu().numpy().reshape(-1)
            if not np.isfinite(pred).all() or np.any(pred < 0) or np.any(pred > 1):
                raise ValueError("Frozen model returned invalid baseline Tumor probabilities")
            result[start:stop] = pred
    return result


def verify_cache(runs, cells, genes, holdout_rows):
    old_source = OLD_ROOT / "source_data"
    prov_path = old_source / "Fig5E_majority28_provenance.json"
    delta_path = old_source / "Fig5E_majority28_paired_deltas_by_seed.csv"
    baseline_path = old_source / "Fig5E_majority28_baseline_by_seed.csv"
    required = (prov_path, delta_path, baseline_path)
    if not all(p.is_file() for p in required):
        return None, None, "legacy perturbation cache incomplete"
    prov = json.loads(prov_path.read_text(encoding="utf-8"))
    expected_cell_hash = hashlib.sha256("\n".join(cells.iloc[holdout_rows].CellID.astype(str)).encode()).hexdigest()
    expected_feature_hash = hashlib.sha256("\n".join(genes.gene_id.astype(str)).encode()).hexdigest()
    checks = [prov.get("status") == "complete", prov.get("contract") == CONTRACT_ID,
        prov.get("feature_count") == EXPECTED_FEATURES, prov.get("seeds") == list(SEEDS),
        prov.get("n_holdout") == EXPECTED_HOLDOUT, prov.get("holdout_cellid_order_sha256") == expected_cell_hash,
        prov.get("input_feature_order_sha256") == expected_feature_hash,
        prov.get("input_transform") == "x_model=log1p(linear_normalized_expression); perturbation uses expm1(x_model)*{0.5,1.5}, then log1p",
        prov.get("batch_size") == BATCH_SIZE,
        prov.get("batch_semantics") == "original holdout sequence, all Group3 labels mixed; aggregate after full-batch predictions",
        prov.get("model_wrapper") == "CauFinder.caufinder_main.ShapModel(module, 'prob')"]
    if not all(checks):
        return None, None, "legacy perturbation provenance does not exactly match current frozen contract/order/parameters"
    current_checkpoints = {str(s): runs[s]["checkpoint_sha256"] for s in SEEDS}
    if prov.get("seed_checkpoints") != current_checkpoints:
        return None, None, "legacy cache model checkpoint SHA differs"
    old_hashes = prov.get("output_sha256", {})
    for path in (delta_path, baseline_path):
        if old_hashes.get(str(path)) != sha256(path):
            return None, None, "legacy perturbation output SHA mismatch"
    delta = read_csv_compat(delta_path)
    baseline = read_csv_compat(baseline_path)
    required_delta = {"gene_id", "seed", "Group3", "scenario", "linear_normalized_expression_multiplier", "n_cells",
                      "baseline_mean_tumor_probability", "perturbed_mean_tumor_probability", "mean_paired_delta_tumor_probability"}
    if not required_delta.issubset(delta.columns) or delta.duplicated(["gene_id", "seed", "Group3", "scenario"]).any():
        return None, None, "legacy cache key/columns invalid"
    if baseline.duplicated(["seed", "Group3"]).any() or len(baseline) != len(SEEDS) * len(GROUPS):
        return None, None, "legacy baseline cache incomplete"
    required_baseline = {"seed", "Group3", "n_cells", "baseline_mean_tumor_probability", "checkpoint_sha256"}
    if not required_baseline.issubset(baseline.columns):
        return None, None, "legacy baseline cache fields incomplete"
    expected_counts = {group: int(np.sum(cells.iloc[holdout_rows].Group3.astype(str).to_numpy() == group)) for group in GROUPS}
    baseline_key_set = set()
    for row in baseline.to_dict("records"):
        seed, group = int(row["seed"]), str(row["Group3"])
        if seed not in SEEDS or group not in GROUPS or int(row["n_cells"]) != expected_counts[group]:
            return None, None, "legacy baseline cache has invalid seed/group cell count"
        if row["checkpoint_sha256"] != runs[seed]["checkpoint_sha256"] or not np.isfinite(float(row["baseline_mean_tumor_probability"])):
            return None, None, "legacy baseline cache checkpoint/value invalid"
        baseline_key_set.add((seed, group))
    if baseline_key_set != {(seed, group) for seed in SEEDS for group in GROUPS}:
        return None, None, "legacy baseline cache does not cover every seed/group key"
    baseline_lookup = {(int(row.seed), str(row.Group3)): float(row.baseline_mean_tumor_probability) for row in baseline.itertuples()}
    expected_scenarios = {name: value for name, value in SCENARIOS}
    valid_gene_set = set()
    for gene in sorted(set(delta.gene_id.astype(str))):
        gene_rows = delta[delta.gene_id.astype(str) == gene]
        valid = len(gene_rows) == len(SEEDS) * len(GROUPS) * len(SCENARIOS)
        if valid:
            valid = not gene_rows.duplicated(["seed", "Group3", "scenario"]).any()
        if valid:
            observed_keys = {(int(row.seed), str(row.Group3), str(row.scenario)) for row in gene_rows.itertuples()}
            expected_keys = {(seed, group, scenario) for seed in SEEDS for group in GROUPS for scenario in expected_scenarios}
            valid = observed_keys == expected_keys
        if valid:
            for row in gene_rows.itertuples():
                seed, group, scenario = int(row.seed), str(row.Group3), str(row.scenario)
                values = (float(row.linear_normalized_expression_multiplier), float(row.n_cells),
                          float(row.baseline_mean_tumor_probability), float(row.perturbed_mean_tumor_probability),
                          float(row.mean_paired_delta_tumor_probability))
                expected_multiplier = expected_scenarios.get(scenario)
                if (not np.isfinite(values).all() or expected_multiplier is None or values[0] != expected_multiplier or int(values[1]) != expected_counts[group] or
                    abs(values[2] - baseline_lookup[(seed, group)]) > 1e-10 or
                    abs(values[4] - (values[3] - values[2])) > 1e-7):
                    valid = False
                    break
        if valid:
            valid_gene_set.add(gene)
    invalid_genes = sorted(set(delta.gene_id.astype(str)) - valid_gene_set)
    delta = delta[delta.gene_id.astype(str).isin(valid_gene_set)].copy()
    note = "cache identity/hashes verified; reusable genes={}; invalid/incomplete cache genes={}".format(len(valid_gene_set), len(invalid_genes))
    return delta, baseline, note


def verify_existing_final_outputs(runs, registry, input_hashes):
    final_files = [OUT_DELTAS, OUT_BASELINE, OUT_HEATMAP, OUT_PARAMS, OUT_PNG, OUT_PDF]
    existing = [path for path in final_files + [OUT_PROVENANCE] if path.exists()]
    if not existing:
        return False
    if not OUT_PROVENANCE.is_file():
        raise FileExistsError("Partial Fig5E final outputs exist without provenance; refusing to treat them as complete or overwrite: {}".format(", ".join(map(str, existing))))
    provenance = json.loads(OUT_PROVENANCE.read_text(encoding="utf-8"))
    if provenance.get("status") != "complete" or any(not path.is_file() for path in final_files):
        raise FileExistsError("Fig5E final output set is partial or provenance is not complete; refusing overwrite")
    if provenance.get("candidate_count") != len(registry) or provenance.get("contract") != CONTRACT_ID or provenance.get("seeds") != list(SEEDS):
        raise ValueError("Existing Fig5E provenance belongs to a different candidate set, contract, or model seeds")
    expected_checkpoints = {str(seed): runs[seed]["checkpoint_sha256"] for seed in SEEDS}
    if provenance.get("seed_checkpoints") != expected_checkpoints:
        raise ValueError("Existing Fig5E provenance checkpoint identities differ from current frozen models")
    recorded_inputs = provenance.get("input_files_sha256", {})
    if any(recorded_inputs.get(path) != digest for path, digest in input_hashes.items()):
        raise ValueError("Existing Fig5E provenance input hashes differ from the current frozen inputs")
    recorded_outputs = provenance.get("output_sha256", {})
    for path in final_files:
        if recorded_outputs.get(str(path)) != sha256(path):
            raise ValueError("Existing Fig5E output hash does not match its provenance: {}".format(path))
    return True


def main():
    started = time.time()
    final_files = (OUT_DELTAS, OUT_BASELINE, OUT_HEATMAP, OUT_PARAMS, OUT_PNG, OUT_PDF)
    existing = [p for p in final_files + (OUT_PROVENANCE,) if p.exists()]
    if existing and not OUT_PROVENANCE.is_file():
        raise FileExistsError("Partial Fig5E final output exists without provenance; refusing overwrite: {}".format(", ".join(map(str, existing))))
    if existing and any(not p.is_file() for p in final_files):
        raise FileExistsError("Fig5E output set is incomplete despite an existing provenance; refusing overwrite")
    registry, cells, genes, holdout, feature_index, contract, export, model_code_hashes = read_inputs()
    runs = verify_seed_runs()
    cache_delta, cache_baseline, cache_note = verify_cache(runs, cells, genes, holdout)
    log("preflight candidates={}, seeds=10, holdout={}, cache={}; RAM={}".format(len(registry), len(holdout), cache_note, memory_status()))
    matrix_path = DATA_DIR / "expression_cells_x_genes.f32"
    expected_bytes = EXPECTED_CELLS * EXPECTED_FEATURES * np.dtype("<f4").itemsize
    if matrix_path.stat().st_size != expected_bytes:
        raise ValueError("Frozen expression data file length differs from the contract")
    # Hash the multi-GB frozen layer in streaming mode; never densify it in memory.
    matrix_hash = sha256(matrix_path)
    expected_matrix_hash = contract.get("fingerprints", {}).get(str(matrix_path))
    if expected_matrix_hash and expected_matrix_hash != matrix_hash:
        raise ValueError("Frozen expression matrix SHA differs from model contract")
    current_input_hashes = {str(p): sha256(p) for p in [DATA_DIR / "cells.csv", DATA_DIR / "genes.csv", DATA_DIR / "contract.json",
        DATA_DIR / "export_complete.json", matrix_path, SOURCE_DIR / "seed42_model_evidence_76_registry.csv"]}
    current_input_hashes.update(model_code_hashes)
    for seed in SEEDS:
        current_input_hashes[str(runs[seed]["checkpoint"])] = runs[seed]["checkpoint_sha256"]
        current_input_hashes[str(RUNS_DIR / "seed_{}".format(seed) / "weights.csv")] = runs[seed]["weights_sha256"]
        current_input_hashes[str(RUNS_DIR / "seed_{}".format(seed) / "complete.json")] = sha256(RUNS_DIR / "seed_{}".format(seed) / "complete.json")
    if verify_existing_final_outputs(runs, registry, current_input_hashes):
        log("Existing Fig5E final outputs and provenance verified exactly; no recomputation required")
        return
    prior_input_hashes = {}
    if cache_delta is not None:
        prior_input_hashes = json.loads((OLD_ROOT / "source_data" / "Fig5E_majority28_provenance.json").read_text(encoding="utf-8")).get("input_files_sha256", {})
        for path, observed in ((matrix_path, matrix_hash), (DATA_DIR / "cells.csv", sha256(DATA_DIR / "cells.csv")),
                               (DATA_DIR / "genes.csv", sha256(DATA_DIR / "genes.csv")),
                               (DATA_DIR / "export_complete.json", sha256(DATA_DIR / "export_complete.json"))):
            if prior_input_hashes.get(str(path)) != observed:
                cache_delta, cache_baseline = None, None
                cache_note = "legacy cache input hash changed: {}".format(path)
                break
    matrix = np.memmap(matrix_path, dtype="<f4", mode="r", shape=(EXPECTED_CELLS, EXPECTED_FEATURES))
    probe_rows = np.unique(np.linspace(holdout[0], holdout[-1], 257).astype(int))
    probe = np.asarray(matrix[probe_rows, :], dtype=np.float32)
    if not np.isfinite(probe).all() or np.any(probe < 0):
        raise ValueError("RNA data probe is not finite nonnegative log1p expression")
    del probe
    labels = cells.iloc[holdout].Group3.astype(str).values
    groups = {group: np.flatnonzero(labels == group) for group in GROUPS}
    cell_ids = cells.iloc[holdout].CellID.astype(str).values
    if len(np.unique(cell_ids)) != EXPECTED_HOLDOUT or any(len(groups[g]) == 0 for g in GROUPS):
        raise ValueError("Holdout CellID or group partitions are invalid")
    progress_identity = {
        "contract": CONTRACT_ID,
        "script_sha256": sha256(Path(__file__).resolve()),
        "candidate_registry_sha256": sha256(SOURCE_DIR / "seed42_model_evidence_76_registry.csv"),
        "expression_matrix_sha256": matrix_hash,
        "cells_sha256": sha256(DATA_DIR / "cells.csv"),
        "genes_sha256": sha256(DATA_DIR / "genes.csv"),
        "holdout_cellid_order_sha256": hashlib.sha256("\n".join(cell_ids).encode()).hexdigest(),
        "feature_order_sha256": hashlib.sha256("\n".join(genes.gene_id.astype(str)).encode()).hexdigest(),
        "seed_checkpoints": {str(s): runs[s]["checkpoint_sha256"] for s in SEEDS},
        "seed_model_weights_sha256": {str(s): runs[s]["weights_sha256"] for s in SEEDS},
        "seed_complete_manifest_sha256": {str(s): current_input_hashes[str(RUNS_DIR / "seed_{}".format(s) / "complete.json")] for s in SEEDS},
        "frozen_model_code_sha256": model_code_hashes,
        "seeds": list(SEEDS), "groups": list(GROUPS),
        "scenarios": [{"name": name, "multiplier": value} for name, value in SCENARIOS],
        "batch_size": BATCH_SIZE,
    }
    PROGRESS_DIR.mkdir(parents=True, exist_ok=True)
    if PROGRESS_STATE.exists():
        prior_identity = json.loads(PROGRESS_STATE.read_text(encoding="utf-8"))
        if prior_identity != progress_identity:
            raise ValueError("Existing Fig5E progress belongs to different frozen inputs/configuration; refusing reuse")
    else:
        PROGRESS_STATE.write_text(json.dumps(progress_identity, ensure_ascii=False, indent=2), encoding="utf-8")
    # Reuse cached baseline means only when their provenance was fully validated.
    if cache_baseline is not None:
        for seed in SEEDS:
            for group in GROUPS:
                row = cache_baseline[(cache_baseline.seed.astype(int) == seed) & (cache_baseline.Group3 == group)]
                if len(row) != 1 or int(row.iloc[0].n_cells) != len(groups[group]) or row.iloc[0].checkpoint_sha256 != runs[seed]["checkpoint_sha256"]:
                    cache_delta, cache_baseline = None, None
                    cache_note = "legacy baseline identity/count check failed"
                    break
            if cache_baseline is None:
                break
    if len(registry) == 0:
        baseline = pd.DataFrame.from_records([], columns=["seed", "Group3", "n_cells", "baseline_mean_tumor_probability", "checkpoint_sha256", "cache_status"])
        log("no model-evaluable seed42 candidates; baseline/perturbation inference is not applicable")
    elif cache_baseline is None:
        baseline_rows = []
        baseline_required = {"seed", "Group3", "n_cells", "baseline_mean_tumor_probability", "checkpoint_sha256"}
        baseline_progress = load_progress_csv(PROGRESS_BASELINE, ["seed", "Group3"], baseline_required)
        if len(baseline_progress):
            for _, saved in baseline_progress.iterrows():
                seed = int(saved["seed"])
                group = str(saved["Group3"])
                if seed not in runs or group not in GROUPS or saved["checkpoint_sha256"] != runs[seed]["checkpoint_sha256"] or int(saved["n_cells"]) != len(groups[group]) or not np.isfinite(float(saved["baseline_mean_tumor_probability"])):
                    raise ValueError("Persisted Fig5E baseline checkpoint identity/count mismatch")
            baseline_rows = baseline_progress.to_dict("records")
        for seed in SEEDS:
            saved_groups = {str(r["Group3"]) for r in baseline_rows if int(r["seed"]) == seed}
            missing_groups = [group for group in GROUPS if group not in saved_groups]
            if not missing_groups:
                log("seed {} baseline restored from durable progress".format(seed))
                continue
            wrapped, device = load_frozen_model(seed, runs[seed]["checkpoint"], genes)
            baseline_pred = predict_baseline(wrapped, device, matrix, holdout)
            for group in missing_groups:
                ids = groups[group]
                baseline_rows.append({"seed": seed, "Group3": group, "n_cells": int(len(ids)),
                    "baseline_mean_tumor_probability": float(baseline_pred[ids].mean(dtype=np.float64)),
                    "checkpoint_sha256": runs[seed]["checkpoint_sha256"],
                    "cache_status": "recomputed_from_frozen_checkpoint_and_holdout"})
            del wrapped, baseline_pred
            if torch.cuda.is_available(): torch.cuda.empty_cache()
            atomic_csv(pd.DataFrame(baseline_rows), PROGRESS_BASELINE)
            log("seed {} frozen baseline predictions recomputed".format(seed))
        baseline = pd.DataFrame(baseline_rows)
    else:
        baseline = cache_baseline.copy()
    if len(baseline) != len(SEEDS) * len(GROUPS) or baseline.duplicated(["seed", "Group3"]).any():
        raise ValueError("Baseline source must contain exactly one row for every frozen seed×Group3 key")
    output_rows = []
    missing_candidates = []
    old_gene_set = set(cache_delta.gene_id.astype(str)) if cache_delta is not None else set()
    for row in registry.to_dict("records"):
        gene = str(row["gene_id"])
        if cache_delta is not None and gene in old_gene_set:
            cached = cache_delta[cache_delta.gene_id.astype(str) == gene].copy()
            if len(cached) != len(SEEDS) * len(GROUPS) * len(SCENARIOS):
                raise ValueError("Cached candidate {} does not have complete seed/group/scenario rows".format(gene))
            for seed in SEEDS:
                subset = cached[cached.seed.astype(int) == seed]
                if len(subset) != len(GROUPS) * len(SCENARIOS):
                    raise ValueError("Cached candidate {} seed {} rows are incomplete".format(gene, seed))
            cached["PTG"] = row.get("PTG", "")
            cached["selected_rank"] = int(row["selected_rank"])
            cached["selected_runs"] = int(row["selected_runs"])
            cached["core_10_of_10"] = bool(row["core_10_of_10"])
            cached["cache_status"] = "reused_exact_verified_legacy_cache"
            output_rows.extend(cached.to_dict("records"))
            log("candidate {} reused verified matched cache".format(gene))
            continue
        missing_candidates.append(gene)
    if missing_candidates:
        feature_index_by_gene = dict(zip(registry.gene_id.astype(str), feature_index.astype(int)))
        candidate_rows = registry.set_index("gene_id").to_dict("index")
        delta_required = {"gene_id", "PTG", "selected_rank", "selected_runs", "core_10_of_10", "seed", "Group3", "scenario",
            "linear_normalized_expression_multiplier", "n_cells", "baseline_mean_tumor_probability", "perturbed_mean_tumor_probability",
            "mean_paired_delta_tumor_probability", "delta_definition", "batch_inference", "cache_status"}
        delta_progress = load_progress_csv(PROGRESS_DELTAS, ["gene_id", "seed", "Group3", "scenario"], delta_required)
        candidate_set = set(missing_candidates)
        if len(delta_progress) and (not set(delta_progress.gene_id.astype(str)).issubset(candidate_set) or
            not set(delta_progress.seed.astype(int)).issubset(set(SEEDS)) or
            not set(delta_progress.Group3.astype(str)).issubset(set(GROUPS)) or
            not set(delta_progress.scenario.astype(str)).issubset({x[0] for x in SCENARIOS})):
            raise ValueError("Persisted Fig5E perturbation progress contains rows outside current frozen candidate keys")
        if len(delta_progress):
            multipliers = {name: value for name, value in SCENARIOS}
            for saved in delta_progress.itertuples():
                gene, seed, group, scenario = str(saved.gene_id), int(saved.seed), str(saved.Group3), str(saved.scenario)
                values = (float(saved.linear_normalized_expression_multiplier), float(saved.n_cells),
                          float(saved.baseline_mean_tumor_probability), float(saved.perturbed_mean_tumor_probability),
                          float(saved.mean_paired_delta_tumor_probability))
                expected_row = candidate_rows[gene]
                expected_baseline = float(baseline[(baseline.seed.astype(int) == seed) & (baseline.Group3.astype(str) == group)].iloc[0].baseline_mean_tumor_probability)
                if (not np.isfinite(values).all() or values[0] != multipliers[scenario] or int(values[1]) != len(groups[group]) or
                    abs(values[2] - expected_baseline) > 1e-10 or abs(values[4] - (values[3] - values[2])) > 1e-7 or
                    str(saved.PTG) != str(expected_row.get("PTG", "")) or int(saved.selected_rank) != int(expected_row["selected_rank"]) or
                    int(saved.selected_runs) != int(expected_row["selected_runs"])):
                    raise ValueError("Persisted Fig5E perturbation progress has a mismatched or invalid completed key: {} / {} / {} / {}".format(gene, seed, group, scenario))
        for seed in SEEDS:
            base_seed = baseline[baseline.seed.astype(int) == seed].set_index("Group3")
            for gene in missing_candidates:
                gene_progress = delta_progress[(delta_progress.gene_id.astype(str) == gene) & (delta_progress.seed.astype(int) == seed)] if len(delta_progress) else pd.DataFrame()
                observed_keys = {(str(row.Group3), str(row.scenario)) for row in gene_progress.itertuples()}
                expected_keys = {(group, scenario) for group in GROUPS for scenario, _ in SCENARIOS}
                missing_keys = expected_keys - observed_keys
                if not missing_keys:
                    continue
                row = candidate_rows[gene]
                gix = feature_index_by_gene[gene]
                gene_rows = []
                if "wrapped" not in locals() or wrapped_seed != seed:
                    if "wrapped" in locals():
                        del wrapped
                        if torch.cuda.is_available(): torch.cuda.empty_cache()
                    wrapped, device = load_frozen_model(seed, runs[seed]["checkpoint"], genes)
                    wrapped_seed = seed
                for scenario, multiplier in SCENARIOS:
                    groups_missing = [group for group in GROUPS if (group, scenario) in missing_keys]
                    if not groups_missing:
                        continue
                    perturbed = predict_gene(wrapped, device, matrix, holdout, gix, multiplier)
                    for group in groups_missing:
                        ids = groups[group]
                        perturbed_mean = float(perturbed[ids].mean(dtype=np.float64))
                        base_mean = float(base_seed.loc[group, "baseline_mean_tumor_probability"])
                        gene_rows.append({"gene_id": gene, "PTG": row.get("PTG", ""), "selected_rank": int(row["selected_rank"]),
                            "selected_runs": int(row["selected_runs"]), "core_10_of_10": bool(row["core_10_of_10"]),
                            "seed": seed, "Group3": group, "scenario": scenario,
                            "linear_normalized_expression_multiplier": multiplier, "n_cells": int(len(ids)),
                            "baseline_mean_tumor_probability": base_mean, "perturbed_mean_tumor_probability": perturbed_mean,
                            "mean_paired_delta_tumor_probability": perturbed_mean - base_mean,
                            "delta_definition": "mean cell-paired probability change; calculated as difference of means over identical frozen holdout cells",
                            "batch_inference": "all holdout rows in frozen order; mixed Group3 batches of 128",
                            "cache_status": "new_inference_frozen_checkpoint"})
                if gene_rows:
                    delta_progress = pd.concat([delta_progress, pd.DataFrame(gene_rows)], ignore_index=True)
                    if delta_progress.duplicated(["gene_id", "seed", "Group3", "scenario"]).any():
                        raise ValueError("Resume logic attempted to append an already complete gene×seed×Group3×scenario key")
                    atomic_csv(delta_progress, PROGRESS_DELTAS)
                if len(delta_progress[(delta_progress.gene_id.astype(str) == gene) & (delta_progress.seed.astype(int) == seed)]) == len(GROUPS) * len(SCENARIOS):
                    log("seed {} candidate {} complete".format(seed, gene))
            if "wrapped" in locals():
                del wrapped
                del wrapped_seed
                if torch.cuda.is_available(): torch.cuda.empty_cache()
            log("seed {} candidate progress checkpointed".format(seed))
        output_rows.extend(delta_progress.to_dict("records"))
    delta_columns = ["gene_id", "PTG", "selected_rank", "selected_runs", "core_10_of_10", "seed", "Group3", "scenario",
        "linear_normalized_expression_multiplier", "n_cells", "baseline_mean_tumor_probability", "perturbed_mean_tumor_probability",
        "mean_paired_delta_tumor_probability", "delta_definition", "batch_inference", "cache_status"]
    delta_df = pd.DataFrame.from_records(output_rows, columns=delta_columns)
    expected_rows = len(registry) * len(SEEDS) * len(GROUPS) * len(SCENARIOS)
    if len(delta_df) != expected_rows or delta_df.duplicated(["gene_id", "seed", "Group3", "scenario"]).any():
        raise ValueError("Perturbation source rows are incomplete or duplicated")
    column_order = ["Neurological_down50", "Neurological_up50", "Tumor_down50", "Tumor_up50"]
    if len(registry):
        heatmap = delta_df.groupby(["gene_id", "Group3", "scenario"]).apply(lambda group: pd.Series({
            "median_paired_delta_tumor_probability": group["mean_paired_delta_tumor_probability"].median(),
            "mean_paired_delta_tumor_probability": group["mean_paired_delta_tumor_probability"].mean(),
            "n_seeds": group["seed"].nunique()})).reset_index()
        heatmap["column"] = heatmap.Group3 + "_" + heatmap.scenario
        heatmap = heatmap.pivot(index="gene_id", columns="column", values="median_paired_delta_tumor_probability")
        heatmap = heatmap.reindex(index=registry.gene_id.astype(str), columns=column_order).reset_index()
    else:
        heatmap = pd.DataFrame.from_records([], columns=["gene_id"] + column_order)
    if heatmap.shape != (len(registry), 5) or (len(registry) and heatmap[column_order].isna().any().any()):
        raise ValueError("Final heatmap values do not match model-evaluable majority candidates")
    plot_heatmap(registry, heatmap, column_order)
    atomic_csv(delta_df, OUT_DELTAS)
    atomic_csv(baseline, OUT_BASELINE)
    atomic_csv(heatmap, OUT_HEATMAP)
    delta_limit = float(np.abs(np.asarray(heatmap[column_order], dtype=float)).max()) if len(registry) else 0.0
    params = pd.DataFrame([
        ("candidate_count", len(registry)), ("seeds", ";".join(map(str, SEEDS))),
        ("contract", CONTRACT_ID), ("holdout_count", len(holdout)),
        ("holdout_groups", "Neurological;Tumor;Ctrl remains reference-only and is not in holdout"),
        ("scenarios", "down50=0.5;up50=1.5 on linear normalized RNA@data scale"),
        ("input_transform", "x_model=log1p(linear_normalized_expression); x_perturbed=log1p(expm1(x_model)*multiplier)"),
        ("batch_size", BATCH_SIZE), ("batch_semantics", "frozen holdout order with original mixed labels"),
        ("delta", "mean paired cell-level Tumor probability change per seed/group/scenario, then median over ten frozen models"),
        ("reused_legacy_gene_cache", ";".join(sorted(set(registry.gene_id.astype(str)) & old_gene_set))),
        ("newly_inferred_genes", ";".join(missing_candidates)), ("color_limits", "-{} to +{}".format(delta_limit, delta_limit)),
        ("png_dpi", 600), ("pdf_fonttype", 42), ("font", "Arial")], columns=["parameter", "value"])
    atomic_csv(params, OUT_PARAMS)
    input_hashes = current_input_hashes
    prov = {"status": "complete", "analysis": "PTG14 pan14 PCSF candidate fixed-model expression perturbation",
        "candidate_count": len(registry), "seeds": list(SEEDS), "features": EXPECTED_FEATURES,
        "contract": CONTRACT_ID, "n_holdout": len(holdout),
        "holdout_group_counts": {g: len(ix) for g, ix in groups.items()},
        "holdout_cellid_order_sha256": hashlib.sha256("\n".join(cell_ids).encode()).hexdigest(),
        "feature_order_sha256": hashlib.sha256("\n".join(genes.gene_id.astype(str)).encode()).hexdigest(),
        "input_transform": "x_model=log1p(linear_normalized_expression); perturbation uses expm1(x_model)*{0.5,1.5}, then log1p",
        "batch_size": BATCH_SIZE, "batch_semantics": "original holdout sequence, all Group3 labels mixed; aggregate after full-batch predictions",
        "model_wrapper": "CauFinder.caufinder_main.ShapModel(module, 'prob')",
        "seed_checkpoints": {str(s): runs[s]["checkpoint_sha256"] for s in SEEDS},
        "cache_reuse_status": cache_note, "new_inference_genes": missing_candidates,
        "input_files_sha256": input_hashes,
        "output_sha256": {str(p): sha256(p) for p in [OUT_DELTAS, OUT_BASELINE, OUT_HEATMAP, OUT_PARAMS, OUT_PNG, OUT_PDF]},
        "elapsed_seconds": round(time.time() - started, 2),
        "interpretation_limit": "descriptive frozen-model perturbation sensitivity; seed models are not biological replicates and perturbation is not experimental validation"}
    atomic_text(json.dumps(prov, ensure_ascii=False, indent=2), OUT_PROVENANCE)
    log("COMPLETE: candidates={}, reused={}, new inference={}; {}".format(len(registry), len(set(registry.gene_id.astype(str)) & old_gene_set), len(missing_candidates), OUT_PNG))


def plot_heatmap(registry, heatmap, columns):
    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42, "axes.unicode_minus": False})
    if len(registry) == 0:
        fig, ax = plt.subplots(figsize=(7.4, 2.5)); ax.axis("off")
        ax.text(0.5, 0.5, "No model-evaluable gene met the ≥6/10 PCSF selection rule", ha="center", va="center", family="Arial", fontsize=10)
        atomic_figure_save(fig, OUT_PNG, "png", dpi=600)
        atomic_figure_save(fig, OUT_PDF, "pdf"); plt.close(fig); return
    values = np.asarray(heatmap[columns], dtype=float)
    vmax = float(np.nanmax(np.abs(values)))
    if not np.isfinite(vmax) or vmax <= 0: vmax = 1e-6
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
    from matplotlib.patches import Patch
    cmap = LinearSegmentedColormap.from_list("blue_white_red", [BLUE, WHITE, "#D73027"])
    height = max(4.5, 0.24 * len(registry) + 1.8)
    fig, ax = plt.subplots(figsize=(7.4, height))
    fig.subplots_adjust(left=0.24, right=0.88, top=0.90, bottom=0.13)
    image = ax.imshow(values, cmap=cmap, norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax), aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(4)); ax.set_xticklabels(["Neurological\n−50%", "Neurological\n+50%", "Tumor\n−50%", "Tumor\n+50%"], fontsize=8, color="black")
    ax.xaxis.tick_top(); ax.tick_params(axis="x", length=0, pad=7)
    ax.set_yticks(np.arange(len(registry))); ax.set_yticklabels(registry.gene_id.astype(str), fontsize=8, color="black")
    for tick, core in zip(ax.get_yticklabels(), registry.core_10_of_10.astype(str).str.lower().isin(["true", "1"])):
        tick.set_color(CORE_RED if core else "black")
    for i, ptg in enumerate(registry.PTG.astype(str)):
        ax.add_patch(plt.Rectangle((-0.016, i - 0.48), 0.012, 0.96, transform=ax.get_yaxis_transform(),
            color=PTG_COLORS.get(ptg, "#999999"), clip_on=False, linewidth=0))
    ax.set_xticks(np.arange(-0.5, 4, 1), minor=True); ax.set_yticks(np.arange(-0.5, len(registry), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8); ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values(): spine.set_visible(False)
    cbar = fig.colorbar(image, ax=ax, fraction=0.04, pad=0.025); cbar.set_label("Median paired Δ Tumor probability", fontsize=8, color="black")
    cbar.ax.tick_params(labelsize=7, colors="black")
    fig.text(0.24, 0.055, "Red labels: selected in 10/10 seeds", va="center", ha="left", fontsize=7.5, color="black", family="Arial")
    ptg_handles = [Patch(facecolor=color, edgecolor="none", label=ptg) for ptg, color in PTG_COLORS.items()]
    fig.legend(handles=ptg_handles, title="Row color: PTG module", loc="lower center", bbox_to_anchor=(0.63, 0.025), ncol=4, frameon=False, fontsize=7, title_fontsize=7)
    atomic_figure_save(fig, OUT_PNG, "png", dpi=600); atomic_figure_save(fig, OUT_PDF, "pdf"); plt.close(fig)


if __name__ == "__main__":
    started_at = time.time()
    try:
        main()
    except BaseException as exc:
        failure = {"status": "failed", "error_type": type(exc).__name__, "error": str(exc),
                   "traceback": traceback.format_exc(), "failed_at": datetime.now(timezone.utc).isoformat()}
        try:
            OUT_FAILURE.parent.mkdir(parents=True, exist_ok=True)
            OUT_FAILURE.write_text(json.dumps(failure, ensure_ascii=False, indent=2), encoding="utf-8")
            log("FAILED: {}: {}".format(type(exc).__name__, exc))
        except Exception:
            pass
        raise

