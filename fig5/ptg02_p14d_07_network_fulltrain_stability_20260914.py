"""Full-training CauFinder network rebuild and 80% subsampling stability.

This script reads the frozen CSR expression cache and adapter source. It does
not train a model, recompute SHAP, or alter old outputs. `benchmark` must pass
before `run` is allowed. Correlations are Spearman correlations from average
ranks including implicit CSR zeros; the rank matrix is an ephemeral memmap.
"""

from __future__ import annotations

import argparse
import ast
import ctypes
import gc
import hashlib
import json
import os
import platform
import shutil
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path

import cvxpy as cp
import networkx as nx
import numpy as np
import pandas as pd
import scipy
import scipy.sparse as sp
from scipy.stats import rankdata
from threadpoolctl import threadpool_info, threadpool_limits


ROOT = Path(r"E:/micro_output/分析审阅中心/待审阅/PTG02_CD_分支疾病与调控候选")
D14 = ROOT / "data/pan14_disease"
ADAPTER = D14 / "caufinder_adapter/CauFinder"
SCRIPT_DIR = Path(r"E:/micro_output/fig3/script")
OUT = ROOT / "network_fulltrain_stability"
DATA = D14
CSR_DATA = DATA / "P14D_data_features_csr_data.f32"
CSR_INDICES = DATA / "P14D_data_features_csr_indices.i32"
CSR_INDPTR = DATA / "P14D_data_features_csr_indptr.i32"
CELL_PATH = DATA / "P14D_counts_features_cells.csv"
GENE_PATH = DATA / "P14D_counts_features_genes.txt"
BLOCK_PATH = DATA / "P14D_features_cell_blocks.csv"
MANIFEST_PATH = DATA / "P14D_cell_manifest.csv.gz"
PRIOR_PATH = DATA / "resources/NicheNet_human.csv"
RANK_PATH = ROOT / "tables/P14D_candidate_ranking.csv"
OLD_CELL_PATH = ROOT / "network_visualization/data/P14D_network_training_cells.csv"
OLD_NODE_PATH = ROOT / "network_visualization/tables/P14D_network_nodes_all.csv"
OLD_EDGE_PATH = ROOT / "network_visualization/tables/P14D_network_edges_all.csv"
OLD_SCORE_PATH = ROOT / "network_visualization/tables/P14D_network_score_reconstruction_comparison.csv"
OLD_CANDIDATE_PATH = ROOT / "network_visualization/tables/P14D_CauFVS_candidates_all31.csv"
PTG_PATH = Path(r"E:/micro_output/fig3/trajectory_S4/data/default/s4h_ptg_membership.csv")
UTILS_PATH = ADAPTER / "utils.py"
DRIVER_PATH = ADAPTER / "driver_regulators.py"
MAIN_PATH = ADAPTER / "caufinder_main.py"
PERTURB_SCRIPT = SCRIPT_DIR / "ptg02_p14d_04_perturb.py"
EXPORT_SCRIPT = SCRIPT_DIR / "ptg02_p14d_05_network_export.py"

PREPARE_PARAMS = {"add_edges_pct": 0.01, "corr_cutoff": 0.7, "uppercase": False}
NETWORK_PARAMS = {"keep_self_loops": True, "average_degree": 10, "weight_degree": True}
SCORE_PARAMS = {"lam": 1.0, "power": 1.0}
OPTIMIZER_PARAMS = {"topK": 50, "driver_union": True, "solver": "SCIP", "lam": 0.5}
OLD_SELECTION_SEED = 42
REPEAT_SEEDS = tuple(range(42, 52))
REPEAT_FRACTION = 0.8
SCIP_SEED = 42
BLAS_THREADS = min(8, os.cpu_count() or 1)


def utc_stamp() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_fingerprint(paths: list[Path]) -> pd.DataFrame:
    rows = []
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
        rows.append({"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    return pd.DataFrame(rows)


def memory_status() -> dict:
    class MEMORYSTATUSEX(ctypes.Structure):
        _fields_ = [
            ("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
            ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
            ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
            ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
            ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
        ]

    status = MEMORYSTATUSEX()
    status.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
        return {"memory_query": "GlobalMemoryStatusEx_failed"}
    return {
        "total_ram_gib": round(status.ullTotalPhys / 2**30, 3),
        "available_ram_gib": round(status.ullAvailPhys / 2**30, 3),
        "available_pagefile_gib": round(status.ullAvailPageFile / 2**30, 3),
    }


def require_inputs() -> None:
    required = [
        CSR_DATA, CSR_INDICES, CSR_INDPTR, CELL_PATH, GENE_PATH, BLOCK_PATH,
        MANIFEST_PATH, PRIOR_PATH, RANK_PATH, OLD_CELL_PATH, OLD_NODE_PATH,
        OLD_EDGE_PATH, OLD_SCORE_PATH, OLD_CANDIDATE_PATH, PTG_PATH, UTILS_PATH, DRIVER_PATH,
        MAIN_PATH, PERTURB_SCRIPT, EXPORT_SCRIPT,
    ]
    absent = [str(path) for path in required if not path.is_file()]
    if absent:
        raise FileNotFoundError("Required frozen inputs are missing:\n" + "\n".join(absent))


def load_frozen_inputs() -> dict:
    require_inputs()
    genes = pd.Index(Path(GENE_PATH).read_text(encoding="utf-8").splitlines(), name="gene_id")
    cells = pd.read_csv(CELL_PATH, dtype=str).iloc[:, 0].astype(str)
    blocks = pd.read_csv(BLOCK_PATH, dtype={"CellID": str})
    meta = pd.read_csv(MANIFEST_PATH, dtype={"CellID": str})
    if blocks["CellID"].astype(str).reset_index(drop=True).equals(cells.reset_index(drop=True)) is False:
        raise ValueError("CSR row order and P14D_features_cell_blocks CellID order differ")
    if meta["CellID"].astype(str).reset_index(drop=True).equals(cells.reset_index(drop=True)) is False:
        raise ValueError("CSR row order and P14D_cell_manifest CellID order differ")
    if genes.has_duplicates or cells.duplicated().any():
        raise ValueError("Frozen gene or CellID order contains duplicates")
    data = np.memmap(CSR_DATA, dtype=np.float32, mode="r")
    indices = np.memmap(CSR_INDICES, dtype=np.int32, mode="r")
    indptr = np.memmap(CSR_INDPTR, dtype=np.int32, mode="r")
    if len(data) != len(indices) or len(indptr) != len(cells) + 1:
        raise ValueError("Frozen CSR binary lengths do not match metadata")
    if np.isnan(data).any() or not np.isfinite(data).all() or float(data.min()) < 0:
        raise ValueError("Frozen RNA data contains non-finite or negative values")
    if int(indptr[-1]) != len(data) or (indices.size and (int(indices.min()) < 0 or int(indices.max()) >= len(genes))):
        raise ValueError("Frozen CSR pointers or feature indices are invalid")
    matrix = sp.csr_matrix((data, indices, indptr), shape=(len(cells), len(genes)), copy=False)
    if matrix.nnz != len(data):
        raise ValueError("Unexpected explicit-zero or CSR canonicalization change")

    train_mask = blocks["split"].astype(str).to_numpy() == "train"
    train_rows = np.flatnonzero(train_mask)
    if len(train_rows) != 203080:
        raise ValueError(f"Expected 203080 frozen train rows, observed {len(train_rows)}")
    if not set(meta.loc[train_rows, "Group3"].astype(str)).issubset({"Neurological", "Tumor"}):
        raise ValueError("Train split contains a group outside frozen Neurological/Tumor labels")
    if not (blocks.loc[train_rows, "Group3"].astype(str).to_numpy() == meta.loc[train_rows, "Group3"].astype(str).to_numpy()).all():
        raise ValueError("Cell-block Group3 and manifest Group3 disagree for train cells")
    if not (blocks.loc[train_rows, "SampleID_Global"].astype(str).to_numpy() == meta.loc[train_rows, "SampleID_Global"].astype(str).to_numpy()).all():
        raise ValueError("Cell-block sample labels and manifest sample labels disagree")

    prior = pd.read_csv(PRIOR_PATH, dtype=str)
    if not {"from", "to"}.issubset(prior.columns):
        raise ValueError("Frozen prior network must contain from/to columns")
    rank = pd.read_csv(RANK_PATH)
    if not {"gene_id", "shap_total"}.issubset(rank.columns) or rank["gene_id"].duplicated().any():
        raise ValueError("Frozen SHAP ranking table is missing a unique gene_id/shap_total mapping")
    rank = rank.set_index("gene_id")
    if set(genes) - set(rank.index):
        raise ValueError("SHAP weights do not cover all frozen feature genes")
    if not np.isfinite(rank.loc[genes, "shap_total"].to_numpy(dtype=float)).all():
        raise ValueError("Frozen shap_total has non-finite values")
    ptg = pd.read_csv(PTG_PATH, dtype=str)
    if not {"gene_id", "PTG"}.issubset(ptg.columns):
        raise ValueError("Frozen full PTG membership table needs gene_id and PTG")
    if ptg["gene_id"].duplicated().any():
        duplicated = ptg.loc[ptg["gene_id"].duplicated(keep=False), "gene_id"].unique().tolist()
        raise ValueError(f"Frozen PTG table has multi-assigned genes; verify before continuing: {duplicated[:30]}")
    if not set(ptg["PTG"].dropna().astype(str)).issubset({"PTG01", "PTG02", "PTG03", "PTG04"}):
        raise ValueError("Frozen membership has unexpected PTG labels")
    return {
        "genes": genes, "cells": cells, "blocks": blocks, "meta": meta,
        "matrix": matrix, "train_rows": train_rows, "prior": prior,
        "rank": rank, "ptg": ptg,
    }


def old_stage_a_cell_indices(blocks: pd.DataFrame, meta: pd.DataFrame) -> np.ndarray:
    train_rows = np.flatnonzero(blocks["split"].astype(str).to_numpy() == "train")
    label = blocks["Group3"].map({"Neurological": 0.0, "Tumor": 1.0}).to_numpy(dtype=float)
    disease = meta["Disease"].astype(str).to_numpy()
    rng = np.random.RandomState(OLD_SELECTION_SEED)
    picked = []
    for group in (0.0, 1.0):
        group_rows = train_rows[label[train_rows] == group]
        disease_counts = pd.Series(disease[group_rows]).value_counts()
        for disease_name, count in disease_counts.items():
            take = min(int(count), max(1, int(round(5000 / 2 / disease_counts.size))))
            strata_rows = group_rows[disease[group_rows] == disease_name]
            picked.extend(rng.choice(strata_rows, size=take, replace=False).tolist())
    return np.array(sorted(set(picked)), dtype=np.int64)[:10000]


def extract_influence_function() -> callable:
    source = UTILS_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(UTILS_PATH))
    function = next((node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "get_influence_score"), None)
    if function is None:
        raise RuntimeError("Adapted get_influence_score function not found")
    namespace = {"np": np, "pd": pd, "nx": nx}
    module = ast.Module(body=[function], type_ignores=[])
    exec(compile(module, str(UTILS_PATH), "exec"), namespace)
    return namespace["get_influence_score"]


def make_prior_graph(prior: pd.DataFrame, genes: pd.Index) -> tuple[nx.DiGraph, pd.DataFrame, list[str]]:
    filtered = prior.loc[prior["from"].isin(genes) & prior["to"].isin(genes), ["from", "to"]].drop_duplicates(subset=["from", "to"])
    graph = nx.from_pandas_edgelist(filtered, "from", "to", create_using=nx.DiGraph)
    ordered_nodes = list(graph.nodes())
    if not ordered_nodes:
        raise ValueError("No prior network nodes overlap frozen model features")
    return graph, filtered, ordered_nodes


def exact_spearman_by_rank_matrix(x: np.ndarray, feature_names: list[str], rank_work_path: Path) -> tuple[np.ndarray, dict]:
    n_cells, n_genes = x.shape
    expected_bytes = n_cells * n_genes * np.dtype(np.float64).itemsize
    if rank_work_path.exists():
        raise FileExistsError(f"Refusing to overwrite rank work file: {rank_work_path}")
    ranks = np.memmap(rank_work_path, dtype=np.float64, mode="w+", shape=(n_genes, n_cells))
    t_rank = time.perf_counter()
    center = (n_cells + 1.0) / 2.0
    for j in range(n_genes):
        ranks[j, :] = rankdata(x[:, j], method="average") - center
    ranks.flush()
    rank_seconds = time.perf_counter() - t_rank
    t_corr = time.perf_counter()
    covariance = np.empty((n_genes, n_genes), dtype=np.float64)
    with threadpool_limits(limits=BLAS_THREADS):
        np.matmul(ranks, ranks.T, out=covariance)
    variance = np.einsum("ij,ij->i", ranks, ranks, optimize=True)
    denom = np.sqrt(variance)
    denominator = denom[:, None] * denom[None, :]
    corr = np.full_like(covariance, np.nan)
    np.divide(covariance, denominator, out=corr, where=denominator > 0)
    corr_seconds = time.perf_counter() - t_corr
    zero_variance = [feature_names[i] for i in np.flatnonzero(denom == 0)]
    del covariance, denominator, denom, variance, ranks
    gc.collect()
    rank_work_path.unlink()
    return corr, {
        "n_cells": int(n_cells), "n_genes": int(n_genes),
        "rank_work_bytes": int(expected_bytes), "rank_seconds": rank_seconds,
        "correlation_seconds": corr_seconds, "zero_variance_genes": zero_variance,
    }


def reconstruct_from_correlation(
    prior: pd.DataFrame,
    model_genes: pd.Index,
    ordered_nodes: list[str],
    spearman: np.ndarray,
) -> tuple[nx.DiGraph, pd.DataFrame, dict]:
    node_index = {gene: i for i, gene in enumerate(ordered_nodes)}
    corr_abs = np.abs(spearman)
    prepare_corr = corr_abs.copy()
    # The source prepare_network zeros the diagonal only for added-edge search;
    # get_network independently recomputes self-loop rho (=1 unless constant).
    np.fill_diagonal(prepare_corr, 0.0)
    corr_df = pd.DataFrame(prepare_corr, index=ordered_nodes, columns=ordered_nodes)

    graph, filtered_prior, checked_order = make_prior_graph(prior, model_genes)
    if checked_order != ordered_nodes:
        raise RuntimeError("Prior graph node order changed during network reconstruction")
    top_k = int(len(graph.nodes()) * (len(graph.nodes()) - 1) / 2 * PREPARE_PARAMS["add_edges_pct"])
    flattened = corr_df.stack().reset_index()
    flattened.columns = ["from", "to", "weight"]
    additional = flattened.loc[flattened["weight"] > PREPARE_PARAMS["corr_cutoff"]]
    additional = additional.nlargest(top_k, "weight")[["from", "to"]]
    existing_edges = set(graph.edges())
    new_edges = [
        (row["from"], row["to"])
        for _, row in additional.iterrows()
        if (row["from"], row["to"]) not in existing_edges
    ]
    graph.add_edges_from(new_edges)

    for source, target in graph.edges():
        rho = corr_abs[node_index[source], node_index[target]]
        if NETWORK_PARAMS["weight_degree"]:
            rho *= (graph.degree(source) + graph.degree(target)) / 2.0
        graph[source][target]["weight"] = float(rho)
    if not NETWORK_PARAMS["keep_self_loops"]:
        graph.remove_edges_from(nx.selfloop_edges(graph))
    target_edges = int(NETWORK_PARAMS["average_degree"] * len(graph.nodes()))
    sorted_edges = sorted(graph.edges(data=True), key=lambda edge: edge[2]["weight"], reverse=True)
    final = nx.DiGraph()
    for source, target, attrs in sorted_edges[:target_edges]:
        final.add_edge(source, target, **attrs)
    largest_component = max(nx.weakly_connected_components(final), key=len)
    final = final.subgraph(largest_component).copy()

    prior_pairs = set(map(tuple, filtered_prior[["from", "to"]].itertuples(index=False, name=None)))
    edge_rows = [
        {"source": source, "target": target, "weight": attrs.get("weight", np.nan),
         "edge_origin": "prior_network" if (source, target) in prior_pairs else "expression_correlation_added"}
        for source, target, attrs in final.edges(data=True)
    ]
    return final, pd.DataFrame(edge_rows, columns=["source", "target", "weight", "edge_origin"]), {
        "prior_filtered_unique_directed_edges": int(len(filtered_prior)),
        "nodes_after_prepare_network": int(len(graph.nodes())),
        "edges_after_prepare_network": int(len(graph.edges())),
        "correlation_edges_selected_before_prior_overlap_filter": int(len(additional)),
        "new_expression_edges_added": int(len(new_edges)),
        "get_network_edge_target": int(target_edges),
        "final_nodes": int(len(final.nodes())),
        "final_edges": int(len(final.edges())),
        "final_weak_components": int(nx.number_weakly_connected_components(final)),
        "largest_strong_component_nodes": int(max((len(c) for c in nx.strongly_connected_components(final)), default=0)),
        "final_edges_with_undefined_weight": int((~np.isfinite(pd.to_numeric(pd.Series([e[2].get("weight", np.nan) for e in final.edges(data=True)]), errors="coerce"))).sum()),
    }


def compare_old_benchmark(graph: nx.DiGraph, edge_frame: pd.DataFrame, score_frame: pd.DataFrame) -> dict:
    old_nodes = pd.read_csv(OLD_NODE_PATH, keep_default_na=True)
    old_edges = pd.read_csv(OLD_EDGE_PATH, keep_default_na=True)
    old_score = pd.read_csv(OLD_SCORE_PATH, keep_default_na=True)
    observed_nodes = set(map(str, graph.nodes()))
    expected_nodes = set(old_nodes["gene_id"].astype(str))
    observed_edges = {(str(r.source), str(r.target)): float(r.weight) if pd.notna(r.weight) else np.nan for r in edge_frame.itertuples(index=False)}
    expected_edges = {(str(r.source), str(r.target)): float(r.weight) if pd.notna(r.weight) else np.nan for r in old_edges.itertuples(index=False)}
    node_equal = observed_nodes == expected_nodes
    observed_edge_keys, expected_edge_keys = set(observed_edges), set(expected_edges)
    edge_equal = observed_edge_keys == expected_edge_keys
    common_edge_keys = observed_edge_keys & expected_edge_keys
    nan_pattern_equal = all(np.isnan(observed_edges[key]) == np.isnan(expected_edges[key]) for key in common_edge_keys)
    finite_diffs = [abs(observed_edges[k] - expected_edges[k]) for k in common_edge_keys if np.isfinite(observed_edges[k]) and np.isfinite(expected_edges[k])]
    old_score_map = old_nodes.set_index("gene_id")["recomputed_network_score"]
    score_map = score_frame["score"]
    common = old_score_map.index.intersection(score_map.index)
    old_vals = pd.to_numeric(old_score_map.reindex(common), errors="coerce").to_numpy(dtype=float)
    new_vals = pd.to_numeric(score_map.reindex(common), errors="coerce").to_numpy(dtype=float)
    score_nan_equal = bool(np.array_equal(np.isnan(old_vals), np.isnan(new_vals)))
    score_diffs = np.abs(old_vals - new_vals)
    finite_score_diffs = score_diffs[np.isfinite(score_diffs)]
    return {
        "old_final_node_count": int(len(expected_nodes)),
        "reconstructed_final_node_count": int(len(observed_nodes)),
        "node_set_equal": bool(node_equal),
        "old_final_edge_count": int(len(expected_edges)),
        "reconstructed_final_edge_count": int(len(observed_edges)),
        "directed_edge_set_equal": bool(edge_equal),
        "directed_edges_missing_vs_old": int(len(expected_edge_keys - observed_edge_keys)),
        "directed_edges_extra_vs_old": int(len(observed_edge_keys - expected_edge_keys)),
        "common_directed_edges": int(len(common_edge_keys)),
        "edge_nan_weight_pattern_equal": bool(nan_pattern_equal),
        "max_abs_edge_weight_difference": float(max(finite_diffs, default=np.nan)),
        "old_score_rows": int(len(old_score)),
        "score_common_node_count": int(len(common)),
        "score_undefined_pattern_equal": score_nan_equal,
        "max_abs_influence_score_difference": float(finite_score_diffs.max()) if finite_score_diffs.size else np.nan,
        "old_expression_edges": int((old_edges["edge_origin"] == "expression_correlation_added").sum()),
        "reconstructed_expression_edges": int((edge_frame["edge_origin"] == "expression_correlation_added").sum()),
    }


def run_benchmark() -> dict:
    inputs = load_frozen_inputs()
    root = OUT
    for folder in ("data", "tables", "logs", "source_data", "figures"):
        (root / folder).mkdir(parents=True, exist_ok=True)
    attempt_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    benchmark_dir = root / f"data/benchmark_old_4852/attempts/{attempt_id}"
    protected_outputs = [
        benchmark_dir / "old_4852_training_cells_reproduced.csv",
        benchmark_dir / "reconstructed_edges.csv",
        root / f"tables/benchmark_old_4852_equivalence_{attempt_id}.csv",
        root / f"tables/benchmark_old_4852_zero_variance_genes_{attempt_id}.csv",
        root / f"tables/resource_timing_estimate_{attempt_id}.csv",
        root / f"data/benchmark_old_4852_environment_{attempt_id}.json",
        root / f"logs/benchmark_old_4852_{attempt_id}.log",
    ]
    existing_outputs = [str(path) for path in protected_outputs if path.exists()]
    if existing_outputs:
        raise FileExistsError("Refusing to overwrite existing benchmark artifacts:\n" + "\n".join(existing_outputs))
    benchmark_dir.mkdir(parents=True, exist_ok=True)
    log_path = root / f"logs/benchmark_old_4852_{attempt_id}.log"

    def log(message: str) -> None:
        line = f"[{utc_stamp()}] {message}"
        print(line, flush=True)
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write(line + "\n")

    all_gene_names = inputs["genes"]
    base_graph, _, network_nodes = make_prior_graph(inputs["prior"], all_gene_names)
    del base_graph
    gene_to_index = {gene: i for i, gene in enumerate(all_gene_names)}
    network_col_idx = np.array([gene_to_index[g] for g in network_nodes], dtype=np.int64)
    old_rows = old_stage_a_cell_indices(inputs["blocks"], inputs["meta"])
    stored_cells = pd.read_csv(OLD_CELL_PATH, dtype=str).iloc[:, 0].astype(str)
    reproduced_cells = inputs["cells"].iloc[old_rows].reset_index(drop=True)
    selection_exact = reproduced_cells.equals(stored_cells.reset_index(drop=True))
    if not selection_exact or len(old_rows) != 4852:
        raise RuntimeError("Original 4,852-cell Stage A selection does not reproduce the frozen cell manifest exactly")
    pd.DataFrame({"CellID": reproduced_cells}).to_csv(benchmark_dir / "old_4852_training_cells_reproduced.csv", index=False)

    log(f"benchmark inputs validated; old sample={len(old_rows)}, network genes={len(network_nodes)}, train={len(inputs['train_rows'])}")
    t = time.perf_counter()
    x_old = inputs["matrix"][old_rows, :][:, network_col_idx].toarray()
    dense_seconds = time.perf_counter() - t
    log(f"old-sample CSR to dense: {dense_seconds:.2f}s, shape={x_old.shape}, dtype={x_old.dtype}")

    t = time.perf_counter()
    reference = pd.DataFrame(x_old, columns=network_nodes).corr(method="spearman").to_numpy(dtype=np.float64)
    reference_seconds = time.perf_counter() - t
    log(f"pandas Spearman reference: {reference_seconds:.2f}s")
    rank_path = benchmark_dir / "rank_working.dat"
    observed, rank_stats = exact_spearman_by_rank_matrix(x_old, network_nodes, rank_path)
    finite_mask = np.isfinite(reference) & np.isfinite(observed)
    corr_nan_equal = bool(np.array_equal(np.isnan(reference), np.isnan(observed)))
    corr_max_diff = float(np.max(np.abs(reference[finite_mask] - observed[finite_mask]))) if finite_mask.any() else np.nan
    corr_match = corr_nan_equal and corr_max_diff <= 1e-12
    log(f"rank-dot vs pandas: mask_equal={corr_nan_equal}, max_abs_diff={corr_max_diff:.3g}")

    final_graph, edge_frame, network_stats = reconstruct_from_correlation(
        inputs["prior"], all_gene_names, network_nodes, observed,
    )
    influence_fn = extract_influence_function()
    gene_info = inputs["rank"].loc[all_gene_names, ["shap_total"]].rename(columns={"shap_total": "weight"})
    influence = influence_fn(final_graph, gene_info, lam=SCORE_PARAMS["lam"], power=SCORE_PARAMS["power"])
    comparison = compare_old_benchmark(final_graph, edge_frame, influence)
    comparison.update({
        "old_cell_selection_order_equal": bool(selection_exact),
        "old_cell_count": int(len(old_rows)),
        "network_feature_count": int(len(network_nodes)),
        "pandas_reference_seconds": float(reference_seconds),
        "rank_seconds": float(rank_stats["rank_seconds"]),
        "rank_dot_seconds": float(rank_stats["correlation_seconds"]),
        "dense_extract_seconds": float(dense_seconds),
        "rank_dot_vs_pandas_nan_mask_equal": corr_nan_equal,
        "rank_dot_vs_pandas_max_abs_difference": corr_max_diff,
        "rank_dot_matches_pandas": bool(corr_match),
        **network_stats,
    })
    comparison["equivalence_pass"] = bool(
        corr_match and comparison["node_set_equal"] and comparison["directed_edge_set_equal"]
        and comparison["edge_nan_weight_pattern_equal"]
        and comparison["max_abs_edge_weight_difference"] <= 1e-10
        and comparison["score_undefined_pattern_equal"]
        and comparison["max_abs_influence_score_difference"] <= 1e-10
    )
    old_edges = pd.read_csv(OLD_EDGE_PATH, keep_default_na=True)
    observed_edge_set = set(zip(edge_frame["source"].astype(str), edge_frame["target"].astype(str)))
    expected_edge_set = set(zip(old_edges["source"].astype(str), old_edges["target"].astype(str)))
    edge_diff_rows = []
    for source, target in sorted(expected_edge_set - observed_edge_set):
        edge_diff_rows.append({"source": source, "target": target, "status": "missing_vs_old"})
    for source, target in sorted(observed_edge_set - expected_edge_set):
        edge_diff_rows.append({"source": source, "target": target, "status": "extra_vs_old"})
    pd.DataFrame(edge_diff_rows, columns=["source", "target", "status"]).to_csv(
        benchmark_dir / "directed_edge_set_differences.csv", index=False,
    )
    pd.DataFrame([comparison]).to_csv(root / f"tables/benchmark_old_4852_equivalence_{attempt_id}.csv", index=False)
    pd.DataFrame([{"gene_id": gene, "zero_variance_in_old_sample": True} for gene in rank_stats["zero_variance_genes"]]).to_csv(
        root / f"tables/benchmark_old_4852_zero_variance_genes_{attempt_id}.csv", index=False,
    )

    # A separate 20,000-train-cell timing probe estimates dense rank/correlation
    # cost at a larger scale; it is not a biological analysis or a repeat.
    probe_rows = inputs["train_rows"][:20000]
    t = time.perf_counter()
    x_probe = inputs["matrix"][probe_rows, :][:, network_col_idx].toarray()
    probe_dense_seconds = time.perf_counter() - t
    probe_corr, probe_stats = exact_spearman_by_rank_matrix(
        x_probe, network_nodes, benchmark_dir / "rank_probe_working.dat",
    )
    del probe_corr, x_probe
    gc.collect()
    projected_cell_ratio = len(inputs["train_rows"]) / len(probe_rows)
    projected_full_corr_seconds = (probe_stats["rank_seconds"] + probe_stats["correlation_seconds"]) * projected_cell_ratio
    projected_main_plus_10_repeats_seconds = projected_full_corr_seconds * (1 + 10 * REPEAT_FRACTION)
    row_nnz = np.diff(inputs["matrix"].indptr)
    n_train_nnz = int(row_nnz[inputs["train_rows"]].sum())
    # Conservative working-set estimate: dense expression + float64 ranks +
    # one sparse row subset + correlation + 1 GiB for graph/stack/Python overhead.
    projected_x_bytes = len(inputs["train_rows"]) * len(network_nodes) * np.dtype(np.float32).itemsize
    projected_rank_bytes = len(inputs["train_rows"]) * len(network_nodes) * np.dtype(np.float64).itemsize
    projected_sparse_copy_bytes = n_train_nnz * (np.dtype(np.float32).itemsize + np.dtype(np.int32).itemsize)
    projected_corr_bytes = len(network_nodes) ** 2 * np.dtype(np.float64).itemsize
    projected_working_gib = (projected_x_bytes + projected_rank_bytes + projected_sparse_copy_bytes + projected_corr_bytes) / 2**30 + 1.0
    available_gib = memory_status().get("available_ram_gib", np.nan)
    disk_free_gib = shutil.disk_usage(OUT.drive + "\\").free / 2**30
    resource_gate = bool(
        np.isfinite(available_gib) and projected_working_gib <= available_gib - 2.0
        and projected_rank_bytes / 2**30 <= disk_free_gib - 2.0
    )
    pd.DataFrame([{
        "timing_probe_cells": len(probe_rows), "timing_probe_dense_extract_seconds": probe_dense_seconds,
        "timing_probe_rank_seconds": probe_stats["rank_seconds"],
        "timing_probe_rank_matrix_multiply_seconds": probe_stats["correlation_seconds"],
        "estimated_full_train_rank_plus_correlation_seconds": projected_full_corr_seconds,
        "estimated_main_plus_10_repeats_rank_plus_correlation_seconds": projected_main_plus_10_repeats_seconds,
        "estimated_peak_working_memory_gib": projected_working_gib,
        "currently_available_ram_gib": available_gib, "e_drive_free_gib": disk_free_gib,
        "largest_transient_rank_file_gib": projected_rank_bytes / 2**30,
        "resource_gate_pass": resource_gate,
    }]).to_csv(root / f"tables/resource_timing_estimate_{attempt_id}.csv", index=False)
    pd.DataFrame([{"source": s, "target": t, "weight": w, "edge_origin": o}
                  for (s, t), (w, o) in _edge_origin_map(edge_frame).items()]).to_csv(
        benchmark_dir / "reconstructed_edges.csv", index=False,
    )
    environment = {
        "timestamp": utc_stamp(), "phase": "benchmark_old_4852",
        "python": sys.version, "platform": platform.platform(), "cpu_count": os.cpu_count(),
        **memory_status(), "e_drive_free_gib": round(shutil.disk_usage(OUT.drive + "\\").free / 2**30, 3),
        "versions": {"numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
                     "networkx": nx.__version__, "cvxpy": cp.__version__},
        "BLAS_threads_cap_for_full_run": BLAS_THREADS,
        "threadpools": threadpool_info(),
        "rank_matrix_strategy": "one ephemeral float64 memmap, feature-wise scipy.rankdata(method=average), centered ranks, full matrix multiplication; implicit sparse zeros are materialized in the selected expression matrix and tied normally",
        "timings_seconds": {"dense_extract_old_4852": dense_seconds, "pandas_reference_old_4852": reference_seconds,
                            "rank_old_4852": rank_stats["rank_seconds"], "rank_dot_old_4852": rank_stats["correlation_seconds"],
                            "dense_extract_probe_20000": probe_dense_seconds,
                            "rank_probe_20000": probe_stats["rank_seconds"],
                            "rank_dot_probe_20000": probe_stats["correlation_seconds"]},
        "benchmark": comparison,
        "resource_gate_pass": resource_gate,
        "projected_full_train_rank_plus_correlation_seconds": projected_full_corr_seconds,
        "projected_full_train_plus_10_repeats_rank_plus_correlation_seconds": projected_main_plus_10_repeats_seconds,
        "projected_peak_working_memory_gib": projected_working_gib,
    }
    (root / f"data/benchmark_old_4852_environment_{attempt_id}.json").write_text(json.dumps(environment, indent=2, ensure_ascii=False), encoding="utf-8")
    log(f"old 4,852-cell benchmark equivalence_pass={comparison['equivalence_pass']}; missing_edges={comparison['directed_edges_missing_vs_old']}; extra_edges={comparison['directed_edges_extra_vs_old']}; max_weight_diff={comparison['max_abs_edge_weight_difference']}; max_score_diff={comparison['max_abs_influence_score_difference']}")
    log(f"timing/resource gate pass={resource_gate}; estimated all correlation seconds including 10 repeats={projected_main_plus_10_repeats_seconds:.1f}; estimated peak working memory GiB={projected_working_gib:.2f}; available RAM GiB={available_gib:.2f}; E free GiB={disk_free_gib:.2f}")
    return comparison


def _edge_origin_map(edge_frame: pd.DataFrame) -> dict:
    return {(str(row.source), str(row.target)): (row.weight, str(row.edge_origin))
            for row in edge_frame.itertuples(index=False)}


def write_csv_new(frame: pd.DataFrame, path: Path, *, compression=None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    frame.to_csv(path, index=False, compression=compression, na_rep="NA", encoding="utf-8-sig")


def write_json_new(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str), encoding="utf-8")


def make_disease_stratified_sample(train_rows: np.ndarray, disease: np.ndarray, seed: int) -> tuple[np.ndarray, pd.DataFrame]:
    rng = np.random.RandomState(seed)
    selected = []
    records = []
    for disease_name in sorted(pd.unique(disease[train_rows]).astype(str)):
        stratum = train_rows[disease[train_rows].astype(str) == disease_name]
        take = int(np.floor(REPEAT_FRACTION * len(stratum)))
        selected.extend(rng.choice(stratum, size=take, replace=False).tolist())
        records.append({"Disease": disease_name, "n_train": int(len(stratum)), "n_selected": take,
                        "fraction_selected": take / len(stratum) if len(stratum) else np.nan,
                        "rounding_rule": "floor(0.8*n_disease)"})
    return np.array(sorted(selected), dtype=np.int64), pd.DataFrame(records)


def load_driver_function():
    import importlib.util

    spec = importlib.util.spec_from_file_location("caufinder_adapted_driver_fulltrain", DRIVER_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import adapted driver regulator module: {DRIVER_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.driver_regulators


def optimize_candidates(driver_function, graph: nx.DiGraph, scores: pd.DataFrame) -> tuple[pd.DataFrame, list[dict], str | None]:
    solve_events: list[dict] = []
    original_solve = cp.Problem.solve

    def solve_with_fixed_scip_seed(problem, *args, **kwargs):
        solver = kwargs.get("solver", args[0] if args else None)
        if solver != cp.SCIP:
            return original_solve(problem, *args, **kwargs)
        solver_options = dict(kwargs.get("scip_params", {}))
        solver_options.setdefault("randomization/randomseedshift", SCIP_SEED)
        solver_options.setdefault("parallel/maxnthreads", 1)
        kwargs["scip_params"] = solver_options
        started = time.perf_counter()
        event = {"solver": "SCIP", "random_seed_shift": SCIP_SEED, "threads": 1}
        try:
            answer = original_solve(problem, *args, **kwargs)
            event["cvxpy_status"] = problem.status
            event["objective_value"] = float(problem.value) if problem.value is not None else np.nan
            extra = getattr(problem.solver_stats, "extra_stats", None)
            model = extra.get("model") if isinstance(extra, dict) else None
            if model is not None:
                event.update({
                    "scip_status": model.getStatus(),
                    "scip_gap": _safe_scip_value(model, "getGap"),
                    "scip_primal_objective": _safe_scip_value(model, "getObjVal"),
                    "scip_dual_bound": _safe_scip_value(model, "getDualbound"),
                    "scip_solving_seconds": _safe_scip_value(model, "getSolvingTime"),
                    "scip_nodes": _safe_scip_value(model, "getNNodes"),
                })
            else:
                event["scip_status"] = "solver_model_not_exposed"
            return answer
        except Exception as exc:
            event["error"] = repr(exc)
            event["scip_status"] = "exception"
            raise
        finally:
            event["wall_seconds"] = time.perf_counter() - started
            solve_events.append(event)

    cp.Problem.solve = solve_with_fixed_scip_seed
    error = None
    try:
        drivers = driver_function(
            graph,
            scores,
            topK=OPTIMIZER_PARAMS["topK"],
            driver_union=OPTIMIZER_PARAMS["driver_union"],
            solver=OPTIMIZER_PARAMS["solver"],
            lam=OPTIMIZER_PARAMS["lam"],
        )
    except Exception as exc:
        drivers = pd.DataFrame(columns=["score", "is_MFVS_driver", "is_CauFVS_driver", "is_CauW_driver"])
        error = f"{type(exc).__name__}: {exc}"
    finally:
        cp.Problem.solve = original_solve
    return drivers, solve_events, error


def _safe_scip_value(model, method_name: str):
    try:
        value = getattr(model, method_name)()
        if isinstance(value, (np.integer, int)):
            return int(value)
        if isinstance(value, (np.floating, float)):
            return float(value)
        return str(value)
    except Exception:
        return np.nan


def build_run_frames(
    graph: nx.DiGraph,
    edges: pd.DataFrame,
    scores: pd.DataFrame,
    drivers: pd.DataFrame,
    optimizer_valid: bool,
    all_features: pd.Index,
    shap_weights: pd.Series,
    ptg_membership: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ordered_genes = sorted(map(str, graph.nodes()))
    score_rows = scores.copy()
    score_rows.index = score_rows.index.astype(str)
    score_rows.index.name = "gene_id"
    score_rows["gene_id"] = score_rows.index
    score_rows["network_rank"] = score_rows["score"].rank(ascending=False, method="min")
    score_rows["shap_total"] = pd.Series(shap_weights).reindex(score_rows.index)
    score_rows["PTG"] = pd.Series(ptg_membership).reindex(score_rows.index).fillna("not_in_PTGs")
    score_rows["in_degree"] = pd.Series(dict(graph.in_degree())).reindex(score_rows.index)
    score_rows["out_degree"] = pd.Series(dict(graph.out_degree())).reindex(score_rows.index)
    score_rows = score_rows.reset_index(drop=True)

    raw = drivers.copy()
    if len(raw):
        raw.index = raw.index.astype(str)
        raw.index.name = "gene_id"
        raw["gene_id"] = raw.index
        raw["shap_total"] = pd.Series(shap_weights).reindex(raw.index).to_numpy()
        raw["PTG"] = pd.Series(ptg_membership).reindex(raw.index).fillna("not_in_PTGs").to_numpy()
        raw = raw.reset_index(drop=True)
    else:
        raw = pd.DataFrame(columns=["gene_id", "score", "is_MFVS_driver", "is_CauFVS_driver", "is_CauW_driver", "shap_total", "PTG"])

    flag_cols = ["is_MFVS_driver", "is_CauFVS_driver", "is_CauW_driver"]
    flag_table = pd.DataFrame({"gene_id": ordered_genes})
    score_lookup = scores["score"].copy()
    score_lookup.index = score_lookup.index.astype(str)
    flag_table["network_score"] = flag_table["gene_id"].map(score_lookup)
    flag_table["shap_total"] = flag_table["gene_id"].map(shap_weights)
    flag_table["PTG"] = flag_table["gene_id"].map(ptg_membership).fillna("not_in_PTGs")
    finite_score = np.isfinite(pd.to_numeric(flag_table["network_score"], errors="coerce"))
    for col in flag_cols:
        if optimizer_valid:
            flag_table[col] = False
            if col in raw.columns and len(raw):
                returned = raw.set_index("gene_id")[col]
                flag_table.loc[flag_table["gene_id"].isin(returned.index), col] = flag_table.loc[
                    flag_table["gene_id"].isin(returned.index), "gene_id"
                ].map(returned).to_numpy()
            flag_table.loc[~finite_score, col] = pd.NA
        else:
            flag_table[col] = pd.NA
    flag_table["optimizer_valid"] = bool(optimizer_valid)
    flag_table["score_status"] = np.where(finite_score, "finite", "undefined_network_score")
    flag_table["in_degree"] = flag_table["gene_id"].map(dict(graph.in_degree()))
    flag_table["out_degree"] = flag_table["gene_id"].map(dict(graph.out_degree()))
    if optimizer_valid:
        is_undefined = ~finite_score
        is_selected = flag_table["is_CauFVS_driver"].fillna(False).astype(bool)
        flag_table["candidate_status"] = np.select(
            [is_undefined, is_selected],
            ["undefined_network_score_selection_not_interpreted", "selected_CauFVS"],
            default="not_selected_by_CauFVS",
        )
    else:
        flag_table["candidate_status"] = "optimizer_invalid_selection_not_interpreted"
    return score_rows, raw, flag_table


def get_passed_benchmark() -> tuple[Path, dict]:
    candidates = sorted((OUT / "data").glob("benchmark_old_4852_environment_*.json"), key=lambda p: p.stat().st_mtime, reverse=True)
    for path in candidates:
        payload = json.loads(path.read_text(encoding="utf-8"))
        benchmark = payload.get("benchmark", {})
        if benchmark.get("equivalence_pass") and payload.get("resource_gate_pass"):
            return path, payload
    raise RuntimeError("No old 4,852-cell equivalence and resource benchmark has passed; full run is blocked")


def run_one(
    run_id: str,
    rows: np.ndarray,
    selection_seed: int | None,
    stratum_table: pd.DataFrame | None,
    inputs: dict,
    network_nodes: list[str],
    network_col_idx: np.ndarray,
    gene_to_index: dict,
    driver_function,
    input_fingerprint_path: Path,
    source_fingerprint_path: Path,
    job_config_hash: str,
) -> dict:
    attempt = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_data_dir = OUT / f"data/runs/{run_id}/{attempt}"
    run_table_dir = OUT / f"tables/runs/{run_id}/{attempt}"
    run_source_dir = OUT / f"source_data/runs/{run_id}/{attempt}"
    log_path = OUT / f"logs/{run_id}_{attempt}.log"
    for path in (run_data_dir, run_table_dir, run_source_dir):
        path.mkdir(parents=True, exist_ok=False)

    def log(message: str) -> None:
        line = f"[{utc_stamp()}] [{run_id}] {message}"
        print(line, flush=True)
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write(line + "\n")

    row_hash = hashlib.sha256(np.asarray(rows, dtype="<i8").tobytes()).hexdigest()
    run_config = {
        "job_config_sha256": job_config_hash, "run_id": run_id,
        "rows_sha256": row_hash, "cell_count": int(len(rows)),
        "selection_seed": selection_seed,
        "selection_definition": "all split=train rows" if selection_seed is None else "within-Disease floor(0.8*n), without replacement, RandomState(seed), rows restored to frozen order",
        "prepare_network": PREPARE_PARAMS, "get_network": NETWORK_PARAMS,
        "get_influence_score": SCORE_PARAMS, "driver_regulators": OPTIMIZER_PARAMS,
        "scip_random_seed_shift": SCIP_SEED, "scip_threads": 1, "blas_threads": BLAS_THREADS,
    }
    run_config_hash = hashlib.sha256(json.dumps(run_config, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
    run_manifest = {
        "run_id": run_id, "attempt": attempt, "status": "running", "job_config_sha256": job_config_hash,
        "run_config": run_config,
        "run_config_sha256": run_config_hash, "started_at": utc_stamp(), "network_valid": False,
        "optimizer_valid": False, "input_fingerprints": str(input_fingerprint_path),
        "source_fingerprints": str(source_fingerprint_path),
        "script_path": str(Path(__file__).resolve()),
    }
    manifest_path = run_data_dir / "run_manifest.json"
    write_json_new(run_manifest, manifest_path)
    log(f"starting run with {len(rows)} cells; row hash={row_hash}")

    cell_meta = inputs["meta"].iloc[rows][["CellID", "Disease", "Group3", "SampleID_Global", "DatasetID"]].copy()
    cell_meta["split"] = "train"
    cell_path = run_source_dir / "training_cells.csv.gz"
    write_csv_new(cell_meta, cell_path, compression="gzip")
    coverage = (cell_meta.groupby("Disease", dropna=False)
                .agg(n_cells=("CellID", "size"), n_samples=("SampleID_Global", "nunique"),
                     n_datasets=("DatasetID", "nunique"), Group3=("Group3", "first"))
                .reset_index())
    coverage_path = run_source_dir / "cell_coverage_by_disease.csv"
    write_csv_new(coverage, coverage_path)
    if stratum_table is not None:
        write_csv_new(stratum_table, run_source_dir / "sampling_strata_counts.csv")

    t = time.perf_counter()
    x = inputs["matrix"][rows, :][:, network_col_idx].toarray()
    dense_seconds = time.perf_counter() - t
    rank_path = run_data_dir / "rank_working.dat"
    spearman, rank_stats = exact_spearman_by_rank_matrix(x, network_nodes, rank_path)
    log(f"expression materialized in {dense_seconds:.1f}s; average ranks {rank_stats['rank_seconds']:.1f}s; full Spearman matrix {rank_stats['correlation_seconds']:.1f}s")
    del x
    gc.collect()

    t = time.perf_counter()
    graph, edges, network_stats = reconstruct_from_correlation(
        inputs["prior"], inputs["genes"], network_nodes, spearman,
    )
    graph_seconds = time.perf_counter() - t
    del spearman
    gc.collect()
    influence_fn = extract_influence_function()
    gene_info = inputs["rank"].loc[inputs["genes"], ["shap_total"]].rename(columns={"shap_total": "weight"})
    t = time.perf_counter()
    scores = influence_fn(graph, gene_info, lam=SCORE_PARAMS["lam"], power=SCORE_PARAMS["power"])
    score_seconds = time.perf_counter() - t
    scores.index = scores.index.astype(str)
    scores.index.name = "gene_id"

    node_path = run_table_dir / "network_nodes.csv"
    preoptimizer_node_path = run_table_dir / "network_nodes_pre_optimizer.csv"
    edge_path = run_table_dir / "network_edges.csv"
    score_path = run_table_dir / "network_scores.csv"
    undefined_score_path = run_table_dir / "undefined_network_scores.csv"
    score_table, _, initial_flags = build_run_frames(
        graph, edges, scores, pd.DataFrame(), False, inputs["genes"],
        inputs["rank"]["shap_total"], inputs["ptg"].set_index("gene_id")["PTG"],
    )
    write_csv_new(score_table, score_path)
    write_csv_new(score_table, preoptimizer_node_path)
    write_csv_new(edges, edge_path)
    undefined_scores = score_table.loc[~np.isfinite(pd.to_numeric(score_table["score"], errors="coerce")),
                                      ["gene_id", "score_out", "score_in", "score", "network_rank", "in_degree", "out_degree"]]
    # Persist diagnostics before the optimizer so a solver failure cannot erase
    # the evidence needed for a minimal-method decision.
    write_csv_new(undefined_scores, undefined_score_path)

    solve_calls_path = run_table_dir / "solver_calls.csv"
    optimizer_table_path = run_table_dir / "optimizer_driver_output.csv"
    candidates_path = run_table_dir / "CauFVS_candidates_full_set.csv"
    flags_path = run_table_dir / "all_final_node_candidate_flags.csv"
    candidate_universe = pd.DataFrame(columns=["gene_id", "score", "is_MFVS_driver", "is_CauFVS_driver", "is_CauW_driver", "shap_total", "PTG"])
    drivers = pd.DataFrame(columns=["score", "is_MFVS_driver", "is_CauFVS_driver", "is_CauW_driver"])
    solver_events: list[dict] = []
    optimizer_error = None
    optimizer_valid = False
    try:
        log(f"network built: {network_stats}; NA scores={len(undefined_scores)}; starting adapted SCIP candidate optimizer")
        drivers, solver_events, optimizer_error = optimize_candidates(driver_function, graph, scores)
        raw_statuses = [str(event.get("scip_status", "")) for event in solver_events]
        solver_calls_valid = (not solver_events) or all(
            status == "optimal" and str(event.get("cvxpy_status", "")) == "optimal"
            for status, event in zip(raw_statuses, solver_events)
        )
        optimizer_valid = optimizer_error is None and solver_calls_valid
    except Exception as exc:
        optimizer_error = f"{type(exc).__name__}: {exc}"
    finally:
        write_csv_new(pd.DataFrame(solver_events), solve_calls_path)

    score_rows, raw_driver_rows, candidate_universe = build_run_frames(
        graph, edges, scores, drivers, optimizer_valid, inputs["genes"],
        inputs["rank"]["shap_total"], inputs["ptg"].set_index("gene_id")["PTG"],
    )
    node_flags = candidate_universe[["gene_id", "is_MFVS_driver", "is_CauFVS_driver", "is_CauW_driver",
                                     "optimizer_valid", "score_status", "candidate_status"]]
    node_table = score_rows.merge(node_flags, on="gene_id", how="left", validate="one_to_one")
    write_csv_new(node_table, node_path)
    write_csv_new(raw_driver_rows, optimizer_table_path)
    write_csv_new(candidate_universe.loc[candidate_universe["is_CauFVS_driver"].eq(True)].copy(), candidates_path)
    write_csv_new(candidate_universe, flags_path)
    selected_genes = candidate_universe.loc[candidate_universe["is_CauFVS_driver"].eq(True), "gene_id"].astype(str).tolist() if optimizer_valid else []
    nonfinite_scores = int(len(undefined_scores))
    optimizer_status = (
        "optimal" if optimizer_valid and solver_events
        else "completed_by_graph_reduction_no_ILP" if optimizer_valid
        else "failed_or_inaccurate"
    )
    network_valid = bool(len(graph) > 0 and len(edges) > 0)
    runtime = {
        "dense_extract_seconds": dense_seconds,
        "rank_seconds": rank_stats["rank_seconds"],
        "rank_dot_seconds": rank_stats["correlation_seconds"],
        "graph_prepare_weight_prune_seconds": graph_seconds,
        "influence_score_seconds": score_seconds,
        "total_seconds_before_optimizer": dense_seconds + rank_stats["rank_seconds"] + rank_stats["correlation_seconds"] + graph_seconds + score_seconds,
        "optimizer_wall_seconds": float(sum(float(row.get("wall_seconds", 0) or 0) for row in solver_events)),
        "wall_seconds": time.perf_counter() - t,
    }
    run_manifest.update({
        "status": "complete" if network_valid and optimizer_valid else "partial_invalid_optimizer" if network_valid else "failed_network",
        "network_valid": network_valid, "optimizer_valid": optimizer_valid,
        "optimizer_status": optimizer_status, "optimizer_error": optimizer_error,
        "solver_call_count": len(solver_events), "selected_CauFVS_candidate_count": len(selected_genes),
        "undefined_network_score_count": nonfinite_scores,
        "network_stats": network_stats, "runtime_seconds": runtime,
        "paths": {"nodes": str(node_path), "edges": str(edge_path), "scores": str(score_path),
                  "nodes_pre_optimizer": str(preoptimizer_node_path),
                  "undefined_scores": str(undefined_score_path), "optimizer_output": str(optimizer_table_path),
                  "candidates": str(candidates_path), "candidate_flags": str(flags_path),
                  "solver_calls": str(solve_calls_path), "cells": str(cell_path), "coverage": str(coverage_path)},
        "finished_at": utc_stamp(),
    })
    if stratum_table is not None:
        run_manifest["sampling_strata_count"] = int(len(stratum_table))
        run_manifest["sampling_selected_cell_count"] = int(len(rows))
    write_json_overwrite_own(run_manifest, manifest_path)
    log(f"run completed: network_valid={network_valid}; optimizer_valid={optimizer_valid}; status={optimizer_status}; selected={len(selected_genes)}; NA_scores={nonfinite_scores}; runtime={runtime}")
    return {
        "run_id": run_id, "manifest_path": str(manifest_path), "manifest": run_manifest,
        "nodes": node_table, "edges": edges, "candidate_flags": candidate_universe,
        "candidates": set(selected_genes), "network_valid": network_valid, "optimizer_valid": optimizer_valid,
    }


def write_json_overwrite_own(payload: dict, path: Path) -> None:
    # This manifest was created by this run attempt immediately before compute.
    if not path.is_file():
        raise FileNotFoundError(f"Run-owned status manifest vanished: {path}")
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str), encoding="utf-8")


def load_completed_run(run_id: str, expected_job_hash: str) -> dict | None:
    manifests = sorted((OUT / f"data/runs/{run_id}").glob("*/run_manifest.json"), key=lambda p: p.stat().st_mtime, reverse=True)
    for path in manifests:
        manifest = json.loads(path.read_text(encoding="utf-8"))
        if manifest.get("status") != "complete" or manifest.get("job_config_sha256") != expected_job_hash:
            continue
        paths = manifest.get("paths", {})
        if not all(Path(paths.get(key, "")).is_file() for key in ("nodes", "edges", "candidate_flags", "candidates")):
            continue
        nodes = pd.read_csv(paths["nodes"], keep_default_na=True)
        edges = pd.read_csv(paths["edges"], keep_default_na=True)
        flags = pd.read_csv(paths["candidate_flags"], keep_default_na=True)
        candidate_frame = pd.read_csv(paths["candidates"], keep_default_na=True)
        return {"run_id": run_id, "manifest_path": str(path), "manifest": manifest,
                "nodes": nodes, "edges": edges, "candidate_flags": flags,
                "candidates": set(candidate_frame["gene_id"].astype(str)),
                "network_valid": True, "optimizer_valid": True, "resumed": True}
    return None


def make_job_config(inputs: dict, input_hashes: pd.DataFrame, source_hashes: pd.DataFrame) -> tuple[dict, str]:
    payload = {
        "input_sha256": dict(zip(input_hashes["path"], input_hashes["sha256"])),
        "source_sha256": dict(zip(source_hashes["path"], source_hashes["sha256"])),
        "prepare_network": PREPARE_PARAMS, "get_network": NETWORK_PARAMS,
        "get_influence_score": SCORE_PARAMS, "driver_regulators": OPTIMIZER_PARAMS,
        "repeat_seeds": REPEAT_SEEDS, "repeat_fraction": REPEAT_FRACTION,
        "repeat_stratification": "Disease, floor(0.8*n), without replacement",
        "scip_seed": SCIP_SEED, "scip_threads": 1, "blas_threads": BLAS_THREADS,
        "input_cells": int(len(inputs["cells"])), "input_genes": int(len(inputs["genes"])),
        "training_cells": int(len(inputs["train_rows"])),
    }
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
    return payload, digest


def write_run_summaries(results: list[dict], inputs: dict, job_config_hash: str) -> dict:
    job_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    table_dir = OUT / f"tables/summary_{job_id}"
    table_dir.mkdir(parents=True, exist_ok=False)
    source_dir = OUT / f"source_data/summary_{job_id}"
    source_dir.mkdir(parents=True, exist_ok=False)
    available_network_runs = [r for r in results if r["network_valid"]]
    available_optimizer_runs = [r for r in results if r["optimizer_valid"]]
    full = next((r for r in results if r["run_id"] == "full_train"), None)
    full_candidates = full["candidates"] if full and full["optimizer_valid"] else set()
    repeat_results = [r for r in results if r["run_id"].startswith("repeat_")]
    valid_repeat_optimizer = [r for r in repeat_results if r["optimizer_valid"]]
    valid_repeat_networks = [r for r in repeat_results if r["network_valid"]]
    old_candidates = set(pd.read_csv(OLD_CANDIDATE_PATH)["gene_id"].astype(str))
    ptg_map = inputs["ptg"].set_index("gene_id")["PTG"].astype(str)
    candidate_union = set(full_candidates)
    for result in valid_repeat_optimizer:
        candidate_union.update(result["candidates"])
    rows = []
    for gene in sorted(candidate_union | old_candidates):
        candidate_repeat_presence = [gene in r["candidates"] for r in valid_repeat_optimizer]
        network_presence = [gene in set(r["nodes"]["gene_id"].astype(str)) for r in valid_repeat_networks]
        all_node_observations = []
        for r in available_network_runs:
            node = r["nodes"].loc[r["nodes"]["gene_id"].astype(str) == gene]
            if len(node):
                all_node_observations.append(node.iloc[0])
        score_values = pd.to_numeric(pd.Series([n.get("score") for n in all_node_observations]), errors="coerce").dropna().to_numpy(dtype=float)
        rank_values = pd.to_numeric(pd.Series([n.get("network_rank") for n in all_node_observations]), errors="coerce").dropna().to_numpy(dtype=float)
        rows.append({
            "gene_id": gene, "PTG": ptg_map.get(gene, "not_in_PTGs"),
            "is_old_31_candidate": gene in old_candidates,
            "selected_full_train": (gene in full_candidates) if full and full["optimizer_valid"] else pd.NA,
            "repeat_selected_count_valid_runs": int(sum(candidate_repeat_presence)),
            "repeat_optimizer_valid_denominator": int(len(valid_repeat_optimizer)),
            "repeat_selection_frequency_valid_runs": (float(np.mean(candidate_repeat_presence)) if valid_repeat_optimizer else np.nan),
            "formal_repeat_frequency_denominator_10": 10 if len(valid_repeat_optimizer) == 10 else pd.NA,
            "formal_repeat_frequency_of_10": (float(np.mean(candidate_repeat_presence)) if len(valid_repeat_optimizer) == 10 else np.nan),
            "network_appearance_count_valid_runs": int(sum(network_presence)),
            "network_valid_repeat_denominator": int(len(valid_repeat_networks)),
            "score_median_all_network_valid_runs": float(np.median(score_values)) if score_values.size else np.nan,
            "score_q25_all_network_valid_runs": float(np.quantile(score_values, 0.25)) if score_values.size else np.nan,
            "score_q75_all_network_valid_runs": float(np.quantile(score_values, 0.75)) if score_values.size else np.nan,
            "score_min_all_network_valid_runs": float(np.min(score_values)) if score_values.size else np.nan,
            "score_max_all_network_valid_runs": float(np.max(score_values)) if score_values.size else np.nan,
            "rank_median_all_network_valid_runs": float(np.median(rank_values)) if rank_values.size else np.nan,
            "rank_q25_all_network_valid_runs": float(np.quantile(rank_values, 0.25)) if rank_values.size else np.nan,
            "rank_q75_all_network_valid_runs": float(np.quantile(rank_values, 0.75)) if rank_values.size else np.nan,
            "old_fulltrain_comparison": "retained" if gene in old_candidates and gene in full_candidates else "new_fulltrain_candidate" if gene not in old_candidates and gene in full_candidates else "old_not_selected_fulltrain" if gene in old_candidates and full and full["optimizer_valid"] else "not_evaluated_fulltrain" if not full or not full["optimizer_valid"] else "repeat_candidate_only",
        })
    stability = pd.DataFrame(rows)
    # Keep current candidate union in the main table while retaining old-only
    # comparison candidates explicitly.
    stability["is_in_full_or_repeat_candidate_union"] = stability["gene_id"].isin(candidate_union)
    write_csv_new(stability, table_dir / "candidate_stability_union_full_and_repeats.csv")

    long_rows = []
    for gene in sorted(candidate_union | old_candidates):
        for result in results:
            node = result["nodes"].loc[result["nodes"]["gene_id"].astype(str) == gene]
            if len(node):
                record = node.iloc[0]
                score = record.get("score", np.nan)
                rank = record.get("network_rank", np.nan)
                node_present = True
                absence_reason = ""
            else:
                score, rank = np.nan, np.nan
                node_present = False
                absence_reason = "not_in_final_network" if result["network_valid"] else "network_run_invalid_or_incomplete"
            if result["optimizer_valid"]:
                selected = gene in result["candidates"]
                selection_status = "selected" if selected else "not_selected"
            else:
                selected = pd.NA
                selection_status = "optimizer_invalid_or_incomplete"
            long_rows.append({
                "gene_id": gene, "PTG": ptg_map.get(gene, "not_in_PTGs"), "run_id": result["run_id"],
                "network_valid": result["network_valid"], "optimizer_valid": result["optimizer_valid"],
                "node_present_in_final_network": node_present, "absence_reason": absence_reason,
                "network_score": score, "network_rank": rank,
                "CauFVS_selected": selected, "selection_status": selection_status,
            })
    write_csv_new(pd.DataFrame(long_rows), table_dir / "candidate_by_run_long.csv")

    jaccard_rows = []
    for result in repeat_results:
        if result["optimizer_valid"] and full and full["optimizer_valid"]:
            union = full_candidates | result["candidates"]
            jac = len(full_candidates & result["candidates"]) / len(union) if union else 1.0
        else:
            jac = np.nan
        jaccard_rows.append({"repeat": result["run_id"], "optimizer_valid": result["optimizer_valid"],
                             "full_train_optimizer_valid": bool(full and full["optimizer_valid"]),
                             "jaccard_vs_full_train": jac,
                             "n_full_train_candidates": len(full_candidates) if full and full["optimizer_valid"] else pd.NA,
                             "n_repeat_candidates": len(result["candidates"]) if result["optimizer_valid"] else pd.NA})
    jaccard = pd.DataFrame(jaccard_rows)
    write_csv_new(jaccard, table_dir / "repeat_candidate_jaccard_vs_full_train.csv")

    full_edges = full["edges"] if full and full["network_valid"] else pd.DataFrame(columns=["source", "target", "edge_origin"])
    all_edge_frames = [("full_train", full_edges)] + [(r["run_id"], r["edges"]) for r in valid_repeat_networks]
    edge_keys = set()
    for _, frame in all_edge_frames:
        if len(frame):
            edge_keys.update((str(x.source), str(x.target)) for x in frame.itertuples(index=False))
    valid_repeat_edge_sets = {
        result["run_id"]: set(zip(result["edges"]["source"].astype(str), result["edges"]["target"].astype(str)))
        for result in valid_repeat_networks
    }
    valid_repeat_node_sets = {
        result["run_id"]: set(result["nodes"]["gene_id"].astype(str))
        for result in valid_repeat_networks
    }
    full_edge_set = set(zip(full_edges["source"].astype(str), full_edges["target"].astype(str))) if len(full_edges) else set()
    candidate_related_genes = candidate_union | old_candidates
    edge_origin_values: dict[tuple[str, str], set[str]] = {}
    for _, frame in all_edge_frames:
        if len(frame):
            for row in frame.itertuples(index=False):
                key = (str(row.source), str(row.target))
                edge_origin_values.setdefault(key, set()).add(str(row.edge_origin))
    edge_rows = []
    for source, target in sorted(edge_keys):
        repeat_contains = []
        endpoint_coobserved = 0
        for result in valid_repeat_networks:
            node_set = valid_repeat_node_sets[result["run_id"]]
            if source in node_set and target in node_set:
                endpoint_coobserved += 1
            repeat_contains.append((source, target) in valid_repeat_edge_sets[result["run_id"]])
        origins = sorted(edge_origin_values.get((source, target), set()))
        edge_rows.append({
            "source": source, "target": target,
            "edge_origin": "|".join(origins) if origins else "not_observed",
            "edge_origin_consistent_across_runs": len(origins) <= 1,
            "present_in_full_train": (source, target) in full_edge_set,
            "candidate_related_edge": source in candidate_related_genes or target in candidate_related_genes,
            "repeat_edge_count_valid_network_runs": int(sum(repeat_contains)),
            "valid_network_repeat_denominator": int(len(valid_repeat_networks)),
            "repeat_edge_frequency_all_valid_network_runs": float(np.mean(repeat_contains)) if valid_repeat_networks else np.nan,
            "repeat_endpoint_coobserved_count": int(endpoint_coobserved),
            "repeat_edge_frequency_when_both_nodes_present": float(sum(repeat_contains) / endpoint_coobserved) if endpoint_coobserved else np.nan,
        })
    edge_stability = pd.DataFrame(edge_rows, columns=[
        "source", "target", "edge_origin", "present_in_full_train", "candidate_related_edge",
        "edge_origin_consistent_across_runs", "repeat_edge_count_valid_network_runs",
        "valid_network_repeat_denominator", "repeat_edge_frequency_all_valid_network_runs",
        "repeat_endpoint_coobserved_count", "repeat_edge_frequency_when_both_nodes_present",
    ])
    write_csv_new(edge_stability, table_dir / "directed_edge_stability_all_network_edges.csv")

    run_rows = []
    for result in results:
        manifest = result["manifest"]
        run_rows.append({"run_id": result["run_id"], "manifest_path": result["manifest_path"],
                         "network_valid": result["network_valid"], "optimizer_valid": result["optimizer_valid"],
                         "optimizer_status": manifest.get("optimizer_status"),
                         "cell_count": manifest["run_config"]["cell_count"],
                         "candidate_count": len(result["candidates"]) if result["optimizer_valid"] else pd.NA,
                         "solver_call_count": manifest.get("solver_call_count"),
                         "undefined_network_score_count": manifest.get("undefined_network_score_count"),
                         "runtime_seconds": manifest.get("runtime_seconds", {}).get("wall_seconds")})
    run_table = pd.DataFrame(run_rows)
    write_csv_new(run_table, table_dir / "run_validity_and_coverage.csv")
    summary_manifest = {
        "summary_job_id": job_id, "job_config_sha256": job_config_hash,
        "created_at": utc_stamp(), "full_train_valid": bool(full and full["optimizer_valid"]),
        "repeat_count_requested": 10, "repeat_optimizer_valid_count": len(valid_repeat_optimizer),
        "repeat_network_valid_count": len(valid_repeat_networks),
        "formal_repeat_frequency_10_valid": bool(len(valid_repeat_optimizer) == 10),
        "old_candidate_count": len(old_candidates),
        "full_train_candidate_count": len(full_candidates) if full and full["optimizer_valid"] else None,
        "full_train_retained_old_candidates": len(old_candidates & full_candidates) if full and full["optimizer_valid"] else None,
        "full_train_new_candidates_vs_old31": len(full_candidates - old_candidates) if full and full["optimizer_valid"] else None,
        "old_candidates_not_selected_full_train": sorted(old_candidates - full_candidates) if full and full["optimizer_valid"] else None,
        "paths": {"candidate_stability": str(table_dir / "candidate_stability_union_full_and_repeats.csv"),
                  "candidate_by_run_long": str(table_dir / "candidate_by_run_long.csv"),
                  "repeat_jaccard": str(table_dir / "repeat_candidate_jaccard_vs_full_train.csv"),
                  "directed_edge_stability": str(table_dir / "directed_edge_stability_all_network_edges.csv"),
                  "run_validity": str(table_dir / "run_validity_and_coverage.csv")},
    }
    write_json_new(summary_manifest, table_dir / "summary_manifest.json")
    return summary_manifest


def run_full_and_repeats() -> dict:
    benchmark_path, benchmark = get_passed_benchmark()
    inputs = load_frozen_inputs()
    root = OUT
    for folder in ("data", "tables", "logs", "source_data", "figures"):
        (root / folder).mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    hash_inputs = [CSR_DATA, CSR_INDICES, CSR_INDPTR, CELL_PATH, GENE_PATH, BLOCK_PATH,
                   MANIFEST_PATH, PRIOR_PATH, RANK_PATH, PTG_PATH]
    hash_sources = [UTILS_PATH, DRIVER_PATH, MAIN_PATH, PERTURB_SCRIPT, EXPORT_SCRIPT, Path(__file__).resolve()]
    input_hashes = file_fingerprint(hash_inputs)
    source_hashes = file_fingerprint(hash_sources)
    fingerprint_path = root / f"source_data/input_fingerprints_{timestamp}.csv"
    write_csv_new(input_hashes.assign(fingerprint_type="frozen_input"), fingerprint_path)
    source_fingerprint_path = root / f"source_data/source_code_fingerprints_{timestamp}.csv"
    write_csv_new(source_hashes.assign(fingerprint_type="source_code"), source_fingerprint_path)
    benchmark_copy = root / f"data/benchmark_used_for_fullrun_{timestamp}.json"
    write_json_new({"benchmark_file": str(benchmark_path), "benchmark": benchmark}, benchmark_copy)

    config_payload, job_hash = make_job_config(inputs, input_hashes, source_hashes)
    config_path = root / f"data/job_configuration_{timestamp}.json"
    write_json_new({"job_config": config_payload, "job_config_sha256": job_hash,
                    "benchmark": str(benchmark_path), "created_at": utc_stamp()}, config_path)
    log_path = root / f"logs/fulltrain_stability_job_{timestamp}.log"

    def log(message: str) -> None:
        line = f"[{utc_stamp()}] {message}"
        print(line, flush=True)
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write(line + "\n")

    genes = inputs["genes"]
    base_graph, _, network_nodes = make_prior_graph(inputs["prior"], genes)
    del base_graph
    gene_to_index = {gene: i for i, gene in enumerate(genes)}
    network_col_idx = np.array([gene_to_index[gene] for gene in network_nodes], dtype=np.int64)
    if len(network_nodes) != 4325:
        raise ValueError(f"Adapted prior-network overlap changed; observed {len(network_nodes)} genes, expected 4325")
    driver_function = load_driver_function()
    log(f"verified benchmark={benchmark_path}; full run starts: train_cells={len(inputs['train_rows'])}; features={len(genes)}; network_nodes_before_pruning={len(network_nodes)}")
    log(f"resource estimate: rank-only total for full plus repeats ~{benchmark['projected_full_train_plus_10_repeats_rank_plus_correlation_seconds']:.1f}s; projected_peak={benchmark['projected_peak_working_memory_gib']:.2f}GiB; current_memory={memory_status()}; threads={BLAS_THREADS}")

    results = []
    planned = [("full_train", inputs["train_rows"], None, None)]
    diseases = inputs["meta"]["Disease"].astype(str).to_numpy()
    for seed in REPEAT_SEEDS:
        repeat_rows, stratum = make_disease_stratified_sample(inputs["train_rows"], diseases, seed)
        planned.append((f"repeat_{seed}", repeat_rows, seed, stratum))

    # Resume only exact, completed, optimizer-valid runs with matching immutable inputs/code.
    for run_id, rows, seed, strata in planned:
        resumed = load_completed_run(run_id, job_hash)
        if resumed is not None:
            log(f"reusing completed matching run {run_id}: {resumed['manifest_path']}")
            results.append(resumed)
            continue
        try:
            result = run_one(run_id, rows, seed, strata, inputs, network_nodes,
                             network_col_idx, gene_to_index, driver_function,
                             fingerprint_path, source_fingerprint_path, job_hash)
            results.append(result)
            if not result["network_valid"]:
                log(f"STOP: {run_id} did not produce a valid network; saved diagnostics; no downstream repeat was launched")
                break
            if not result["optimizer_valid"]:
                undefined_count = int(result["manifest"].get("undefined_network_score_count", 0))
                if undefined_count > 0:
                    log(f"STOP: {run_id} optimizer invalid with {undefined_count} undefined influence scores; no values were zero-filled; save minimum-method decision and stop")
                    break
                log(f"WARNING: {run_id} optimizer status is invalid/inaccurate; continuing other independent scheduled runs but excluding this run from selection frequencies")
        except Exception as exc:
            log(f"FATAL RUN ERROR in {run_id}: {type(exc).__name__}: {exc}\n{traceback.format_exc()}")
            raise

    summary = write_run_summaries(results, inputs, job_hash)
    write_json_new(summary, root / f"data/stability_summary_manifest_{timestamp}.json")
    log(f"summary complete: full_train_valid={summary['full_train_valid']}; repeat_optimizer_valid={summary['repeat_optimizer_valid_count']}/10; frequency_denominator_10_valid={summary['formal_repeat_frequency_10_valid']}")
    return {"summary": summary, "job_hash": job_hash, "results": results,
            "timestamp": timestamp, "job_log": str(log_path)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["benchmark", "run"], default="benchmark")
    args = parser.parse_args()
    if args.phase == "benchmark":
        result = run_benchmark()
        if not result.get("equivalence_pass"):
            raise SystemExit("Benchmark mismatch: full training and stability runs are blocked pending diagnosis.")
        print("Benchmark passed; implementation may proceed to full-train and repeat runs.")
    elif args.phase == "run":
        result = run_full_and_repeats()
        print(json.dumps({"summary": result["summary"], "job_hash": result["job_hash"],
                          "job_log": result["job_log"]}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
