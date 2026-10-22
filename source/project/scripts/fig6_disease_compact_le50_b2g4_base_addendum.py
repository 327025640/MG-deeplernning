#!/usr/bin/env python3
"""Exploratory shared-default PCSF parameter assessment for Fig6 Disease-vs-Ctrl.

This script reads only the frozen Fig6 SHAP-to-PPI tables and invokes the
installed OmicsIntegrator2/pcst-fast implementation. It writes exclusively to
the new robustness run directory and never modifies formal networks.
"""
from __future__ import print_function

import argparse
import csv
import hashlib
import io
import itertools
import json
import math
import os
import platform
import sys
import time
import traceback
from datetime import datetime

import networkx as nx
import numpy as np
import pandas as pd
import pkg_resources
from OmicsIntegrator import graph as oi_graph


ANALYSIS_ROOT = r"E:\micro_output\分析审阅中心\待审阅\Fig6_疾病vsCtrl_SHAP_PCSF模块"
MODEL_ROOT = r"E:\micro_output\分析审阅中心\待审阅\Fig6_疾病vsCtrl_全基因CauFinder"
RUN_ID = "20261021_disease_compact_le50_parameter_selection_01_b2g4_addendum"
REUSE_ROOT = os.path.join(ANALYSIS_ROOT, "robustness", "20261020_disease_parameter_adaptation_08")
RUN = os.path.join(ANALYSIS_ROOT, "robustness", "20261021_disease_compact_le50_parameter_selection_01", "base_scan_b2g4_addendum")
PPI_TABLE = os.path.join(ANALYSIS_ROOT, "tables", "ppi_edges_ge2.tsv")
PPI_SOURCE = r"D:\data\microglia_new\input\drug_repurposing_38core\mentor_20260923\Network-1.txt"
PPI_NODE_TABLE = os.path.join(ANALYSIS_ROOT, "tables", "ppi_node_audit.csv")
INPUT_SUMMARY = os.path.join(ANALYSIS_ROOT, "tables", "input_disease_audit_summary.csv")
BASE_CONFIG = os.path.join(ANALYSIS_ROOT, "config", "run_config.json")
DISEASES = ["AD", "ALS", "ALS_FTD", "DLB", "DLBD_PD", "EOAD", "FTD", "FTD_MND",
            "GBM", "HD", "LCBM", "LGG", "MB", "MS", "PD", "PSP"]
ACTIVE_DISEASES = list(DISEASES)
B_VALUES = [2.0]
W_VALUES = [3, 6]
G_VALUES = [4]
REPEATS = 0
SIGMA = 0.1
PRUNING = "strong"
BASE_COST = 1.0
CTRL_TYPES = []
PRIOR_INPUT_HASHES = {}
INPUT_HASH_PROVENANCE = {}


def stamp():
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _input_sha256(path):
    absolute = os.path.abspath(path)
    key = os.path.normcase(absolute)
    previous = PRIOR_INPUT_HASHES.get(key)
    size = os.path.getsize(path)
    modified = datetime.fromtimestamp(os.path.getmtime(path)).astimezone().isoformat(timespec="seconds")
    if previous is not None:
        if int(previous["size_bytes"]) != size or previous["last_modified"] != modified:
            raise RuntimeError("Frozen input metadata changed since the accepted input SHA audit: {}".format(path))
        if size > 100 * 1024 * 1024:
            value = previous["sha256"]
            method = "reused_prior_SHA_after_size_and_timestamp_match"
        else:
            value = sha256(path)
            if value != previous["sha256"]:
                raise RuntimeError("Frozen input SHA changed since the accepted input audit: {}".format(path))
            method = "recomputed_and_matched_prior_SHA"
    else:
        value = sha256(path)
        method = "fresh_SHA256"
    INPUT_HASH_PROVENANCE[key] = {"path": absolute, "sha256": value, "size_bytes": size,
        "last_modified": modified, "method": method,
        "prior_audit_path": os.path.join(REUSE_ROOT, "inputs", "input_sha256_manifest.csv") if previous is not None else ""}
    return value


def _load_prior_input_hashes():
    global PRIOR_INPUT_HASHES
    prior_path = os.path.join(REUSE_ROOT, "inputs", "input_sha256_manifest.csv")
    if not os.path.isfile(prior_path):
        raise FileNotFoundError("Accepted prior input hash manifest missing: {}".format(prior_path))
    prior_df = read_table(prior_path)
    required = {"path", "sha256", "size_bytes", "last_modified"}
    if not required.issubset(prior_df.columns):
        raise RuntimeError("Prior input hash manifest has unexpected columns")
    PRIOR_INPUT_HASHES = {}
    for row in prior_df.to_dict("records"):
        key = os.path.normcase(os.path.abspath(row["path"]))
        if key in PRIOR_INPUT_HASHES:
            raise RuntimeError("Duplicate path in prior input SHA manifest: {}".format(row["path"]))
        PRIOR_INPUT_HASHES[key] = row


def write_json(path, value):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True,
                  allow_nan=False)
        handle.write("\n")
    os.replace(tmp, path)


def write_rows(path, rows, fields=None, delimiter=","):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if fields is None:
        fields = list(rows[0].keys()) if rows else ["status"]
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter,
                                extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            clean = {}
            for key in fields:
                value = row.get(key, "")
                if value is None or (isinstance(value, float) and not math.isfinite(value)):
                    value = ""
                clean[key] = value
            writer.writerow(clean)
    os.replace(tmp, path)


def read_table(path, sep=None):
    actual_sep = sep if sep is not None else ("\t" if path.lower().endswith((".tsv", ".txt")) else ",")
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        return pd.read_csv(handle, sep=actual_sep, dtype=str, keep_default_na=False)


def log(message):
    line = "[{}] {}".format(stamp(), message)
    os.makedirs(os.path.join(RUN, "logs"), exist_ok=True)
    with open(os.path.join(RUN, "logs", "execution.log"), "a", encoding="utf-8") as handle:
        handle.write(line + "\n")
    print(line, flush=True)


def param_id(b, w, g):
    b_value = float(b)
    if not math.isfinite(b_value):
        raise ValueError("Non-finite b value")
    b_text = format(b_value, ".15g").replace("-", "m").replace(".", "p")
    return "b{}_w{}_g{}".format(b_text, int(w), int(g))


def stable_seed(key):
    digest = hashlib.sha256(key.encode("utf-8")).digest()
    return int.from_bytes(digest[:4], byteorder="little", signed=False)


def jaccard(a, b):
    a, b = set(a), set(b)
    if not a and not b:
        return None
    return float(len(a & b)) / float(len(a | b)) if a | b else None


def degree_bins(indices, degrees, min_size=100):
    by_degree = {}
    for ix in indices:
        by_degree.setdefault(int(degrees[int(ix)]), []).append(int(ix))
    pending, bins = [], []
    low = high = None
    for deg, members in sorted(by_degree.items()):
        if low is None:
            low = deg
        high = deg
        pending.extend(sorted(members))
        if len(pending) >= min_size:
            bins.append({"degree_min": low, "degree_max": high, "indices": pending})
            pending, low, high = [], None, None
    if pending:
        if bins:
            bins[-1]["degree_max"] = high
            bins[-1]["indices"].extend(pending)
        else:
            bins.append({"degree_min": low, "degree_max": high, "indices": pending})
    if len(indices) >= min_size and any(len(x["indices"]) < min_size for x in bins):
        raise RuntimeError("Degree-bin minimum-size rule failed")
    return bins


def initialize_layout():
    for folder in ["config", "inputs", "solver_audit", "grid", "perturbations",
                   "randomizations", "summaries", "figure_source", "logs", "scripts"]:
        os.makedirs(os.path.join(RUN, folder), exist_ok=True)


def write_run_config():
    import shutil
    source_script = os.path.abspath(__file__)
    script_hash = sha256(source_script)
    copy_path = os.path.join(RUN, "scripts", os.path.basename(source_script))
    if os.path.abspath(source_script).lower() != os.path.abspath(copy_path).lower():
        if os.path.isfile(copy_path) and sha256(copy_path) != script_hash:
            archive_script = os.path.join(RUN, "scripts", "preflight_attempt_01_" + sha256(copy_path)[:12] + ".py")
            if os.path.isfile(archive_script):
                raise RuntimeError("Preflight script archive already exists; refusing to overwrite it")
            os.replace(copy_path, archive_script)
        if not os.path.isfile(copy_path):
            shutil.copy2(source_script, copy_path)
    config = {"run_id": RUN_ID, "analysis": "Fig6 Disease-vs-Ctrl PCSF compact le50 b=2,g=4 missing-cell base-grid addendum",
        "status": "computes only the missing in-grid b=2,w=3/6,g=4 base cells; formal networks unchanged",
        "diseases": DISEASES, "model_root": MODEL_ROOT, "analysis_root": ANALYSIS_ROOT,
        "model_identity": "frozen disease-specific Disease-vs-Ctrl; model training seed 42; no retraining or SHAP recomputation",
        "model_features": 17645,
        "reward": "target training-cell mean_abs_SHAP divided by disease-wide maximum over all model genes before PPI mapping; mapped Entrez max rule reused; b applied once; no PPI renormalization",
        "ppi": {"source": r"D:\data\microglia_new\input\drug_repurposing_38core\mentor_20260923\Network-1.txt",
            "filter": "number>=2", "undirected": True, "remove_self_loops": True,
            "deduplicate_unordered_edges": True, "retain_all_components": True,
            "expected_nodes": 12512, "expected_edges": 83065, "base_cost": 1.0},
        "solver": {"implementation": "OmicsIntegrator2 2.3.10 + pcst-fast 1.0.7",
            "pruning": PRUNING, "dummy_mode": "terminals", "seed": 42, "edge_noise": 0,
            "hub_penalty": "10^g*deg(u)*deg(v)/((N-deg(u)-1)*(N-deg(v)-1)+deg(u)*deg(v))",
            "objective": "sum(selected real/root edge costs)+sum(unselected real-node scaled prizes); OI's w cost is represented by root edges exactly once"},
        "grid": {"b": B_VALUES, "w": W_VALUES, "g": G_VALUES,
            "expected_disease_parameter_runs": len(DISEASES)*len(B_VALUES)*len(W_VALUES)*len(G_VALUES),
            "fractional_b_path_encoding": "decimal point replaced with p; no integer truncation",
            "predefined_comparison_settings": [param_id(5,3,5),param_id(10,3,3)]},
        "reuse": {"source_run_id": "20261020_disease_parameter_adaptation_08", "source_root": REUSE_ROOT, "exact_b2_g4_cache_found": False, "external_exact_cache_assembly": "A manifest scan found no exact b=2,g=4 result; this addendum computes only b=2,w=3/6,g=4 for the 16 frozen diseases. Other grid cells remain separately source-verified."},
        "perturbations": {"run": False, "planned_repeats": 0, "types": [], "reason": "base-network scan only", "lognormal_sigma": SIGMA,
            "edge_cost": "effective real-edge costs multiplied by exp(N(0,sigma^2)-sigma^2/2); virtual-root edges unchanged",
            "degree_bins": "keep equal-degree strata intact; accumulate ordered degree groups to >=100 nodes; merge trailing remainder into prior bin; retain a single small bin"},
        "comparison_witnesses": ["raw_pcst_fast", "root_only_empty", "best_single_node", "best_single_real_PPI_edge", "root_star_prize_gt_w"],
        "execution_script_path": source_script, "execution_script_sha256": script_hash,
        "selection_principle": "per-disease metrics first; across diseases equal-weight summaries; multi-metric Pareto and neighborhood continuity; no weighted composite score",
        "created_at": stamp()}
    cfg_path = os.path.join(RUN, "config", "run_config.json")
    if os.path.isfile(cfg_path):
        with open(cfg_path, "r", encoding="utf-8") as handle:
            old = json.load(handle)
        if old.get("run_id") != RUN_ID:
            raise RuntimeError("Existing run id differs; refusing to mix outputs")
        if old.get("execution_script_sha256") != script_hash:
            has_solver_outputs = any(
                os.path.isfile(os.path.join(root, "run_manifest.json"))
                for root, _, files in os.walk(os.path.join(RUN, "grid"))
            )
            if has_solver_outputs:
                raise RuntimeError("Script changed after solver outputs exist; refusing to mix outputs")
            archived_config = os.path.join(RUN, "config", "preflight_attempt_01_" + str(old.get("execution_script_sha256", "unknown"))[:12] + ".json")
            if os.path.isfile(archived_config):
                raise RuntimeError("Preflight config archive already exists; refusing to overwrite it")
            os.replace(cfg_path, archived_config)
        else:
            return old
    write_json(cfg_path, config)
    return config


def load_inputs():
    _load_prior_input_hashes()
    for p in [PPI_TABLE, PPI_NODE_TABLE, INPUT_SUMMARY, BASE_CONFIG]:
        if not os.path.isfile(p):
            raise FileNotFoundError(p)
    ppi_df = read_table(PPI_TABLE, "\t")
    node_df = read_table(PPI_NODE_TABLE)
    disease_audit = read_table(INPUT_SUMMARY)
    with open(BASE_CONFIG, "r", encoding="utf-8-sig") as handle:
        frozen_config = json.load(handle)
    if (set(frozen_config.get("diseases", [])) != set(DISEASES) or
            frozen_config.get("expected_model_genes") != 17645 or
            frozen_config.get("pcsf", {}).get("implementation") != "OmicsIntegrator2 author implementation + pcst_fast"):
        raise RuntimeError("Frozen Fig6 configuration identity does not match expected Disease-vs-Ctrl inputs")
    required_edges = {"Protein_A_Entrez_ID", "Protein_B_Entrez_ID", "number_values"}
    if not required_edges.issubset(set(ppi_df.columns)):
        raise RuntimeError("PPI字段不完整: {}".format(sorted(required_edges - set(ppi_df.columns))))
    if not {"Entrez", "background_degree"}.issubset(set(node_df.columns)):
        raise RuntimeError("PPI节点审计缺字段")
    ppi_df["Protein_A_Entrez_ID"] = ppi_df["Protein_A_Entrez_ID"].astype(str)
    ppi_df["Protein_B_Entrez_ID"] = ppi_df["Protein_B_Entrez_ID"].astype(str)
    ppi_number = pd.to_numeric(ppi_df["number_values"], errors="raise")
    if not np.all(np.isfinite(ppi_number.values)) or np.any(ppi_number.values < 2):
        raise RuntimeError("冻结ppi_edges_ge2.tsv包含number<2或非有限number记录")
    ppi_df["number_numeric_audit"] = ppi_number.astype(float)
    if (ppi_df["Protein_A_Entrez_ID"] == ppi_df["Protein_B_Entrez_ID"]).any():
        raise RuntimeError("冻结PPI表中仍存在自环")
    endpoint_keys = [tuple(sorted(x)) for x in zip(ppi_df["Protein_A_Entrez_ID"], ppi_df["Protein_B_Entrez_ID"])]
    if len(endpoint_keys) != len(set(endpoint_keys)):
        raise RuntimeError("冻结PPI表中存在重复无向端点")
    ppi_ids = set(ppi_df["Protein_A_Entrez_ID"]) | set(ppi_df["Protein_B_Entrez_ID"])
    node_ids = set(node_df["Entrez"].astype(str))
    if len(ppi_df) != 83065 or len(ppi_ids) != 12512 or ppi_ids != node_ids:
        raise RuntimeError("PPI输入规模/节点对应与冻结清单不一致: edges={}, nodes={}, audit_nodes={}".format(len(ppi_df), len(ppi_ids), len(node_ids)))
    if disease_audit["Disease"].astype(str).duplicated().any() or set(DISEASES) != set(disease_audit["Disease"].astype(str)):
        raise RuntimeError("16病冻结疾病审计清单不一致")
    input_files = [PPI_SOURCE, PPI_TABLE, PPI_NODE_TABLE, INPUT_SUMMARY, BASE_CONFIG,
                   os.path.join(ANALYSIS_ROOT, "data", "orgHs_mapping_resource_manifest.json"),
                   os.path.join(ANALYSIS_ROOT, "data", "orgHs_official_symbol_entrez.tsv"),
                   os.path.join(ANALYSIS_ROOT, "data", "orgHs_alias_entrez.tsv"),
                   os.path.join(MODEL_ROOT,"input","gene_manifest_model.csv"),
                   os.path.join(MODEL_ROOT,"input","gene_manifest_all_RNA.csv")]
    for path in input_files:
        if not os.path.isfile(path): raise FileNotFoundError(path)
    feature_manifest=read_table(os.path.join(MODEL_ROOT,"input","gene_manifest_model.csv"))
    feature_manifest=feature_manifest[feature_manifest["include_in_model"].astype(str).str.lower()=="true"].copy()
    feature_manifest["_rna_row_numeric_order"]=pd.to_numeric(feature_manifest["RNA_row_index"],errors="raise")
    feature_manifest=feature_manifest.sort_values("_rna_row_numeric_order").drop(columns=["_rna_row_numeric_order"])
    if len(feature_manifest)!=17645 or feature_manifest["gene_id"].astype(str).nunique()!=17645:
        raise RuntimeError("冻结模型特征manifest未能提供17,645个唯一特征")
    frozen_feature_order=feature_manifest["gene_id"].astype(str).tolist()
    frozen_feature_index=dict(zip(frozen_feature_order,pd.to_numeric(feature_manifest["RNA_row_index"],errors="raise").astype(int)))
    feature_order_sha256=hashlib.sha256("\n".join(frozen_feature_order).encode("utf-8")).hexdigest()
    input_rows = []
    model_rows = []
    disease_data = {}
    blocked_disease_rows = []
    for disease in DISEASES:
        try:
            prize_path = os.path.join(ANALYSIS_ROOT, "tables", disease + "_gene_prizes.csv")
            audit_path = os.path.join(ANALYSIS_ROOT, "tables", disease + "_all_gene_selection_audit.csv")
            model_dir = os.path.join(MODEL_ROOT, "models", disease)
            model_config_path = os.path.join(model_dir, "model_config.json")
            runtime_path = os.path.join(model_dir, "runtime.json")
            epoch_training_path = os.path.join(model_dir, "epoch_training.jsonl")
            target_shap_path = os.path.join(MODEL_ROOT, "tables", disease + "_SHAP_target_full_ranking.csv")
            ctrl_shap_path = os.path.join(MODEL_ROOT, "tables", disease + "_SHAP_Ctrl_full_ranking.csv")
            training_log_path=os.path.join(model_dir,"training_log.csv")
            network_manifest_path = os.path.join(ANALYSIS_ROOT, "networks", disease, "run_manifest.json")
            formal_nodes_path = os.path.join(ANALYSIS_ROOT, "networks", disease, "selected_nodes.csv")
            formal_edges_path = os.path.join(ANALYSIS_ROOT, "networks", disease, "selected_forest_edges.csv")
            if not os.path.isfile(network_manifest_path):
                raise FileNotFoundError("{} frozen network provenance manifest missing".format(disease))
            for p in [prize_path, audit_path, model_config_path, runtime_path, target_shap_path, ctrl_shap_path,
                  formal_nodes_path, formal_edges_path,training_log_path,epoch_training_path]:
                if not os.path.isfile(p):
                    raise FileNotFoundError("{} input missing for {}".format(p, disease))
                input_files.append(p)
            with open(model_config_path, "r", encoding="utf-8-sig") as handle:
                model_cfg = json.load(handle)
            with open(runtime_path, "r", encoding="utf-8-sig") as handle:
                runtime = json.load(handle)
            with open(network_manifest_path, "r", encoding="utf-8-sig") as handle:
                network_manifest = json.load(handle)
            formal_nodes_df=read_table(formal_nodes_path)
            formal_edges_df=read_table(formal_edges_path)
            if {"Entrez_A","Entrez_B"}.issubset(set(formal_edges_df.columns)):
                formal_edge_a,formal_edge_b="Entrez_A","Entrez_B"
            elif {"Protein_A_Entrez_ID","Protein_B_Entrez_ID"}.issubset(set(formal_edges_df.columns)):
                formal_edge_a,formal_edge_b="Protein_A_Entrez_ID","Protein_B_Entrez_ID"
            else:
                formal_edge_a,formal_edge_b=None,None
            if "Entrez" not in formal_nodes_df.columns or formal_edge_a is None:
                raise RuntimeError("{} formal b5w3g5 network files lack expected node/edge fields".format(disease))
            if (model_cfg.get("Disease") != disease or model_cfg.get("positive_class") != disease or
                    model_cfg.get("negative_class") != "Ctrl" or model_cfg.get("other_diseases_enter_model") is not False or
                    int(model_cfg.get("n_model_genes", -1)) != 17645 or int(model_cfg.get("seed", -1)) != 42):
                raise RuntimeError("{} model identity does not match frozen Disease-vs-Ctrl contract".format(disease))
            shap_df = read_table(target_shap_path)
            required_shap = {"target_disease", "contribution_group", "n_cells", "gene_id", "mean_abs_SHAP", "rank"}
            if not required_shap.issubset(set(shap_df.columns)) or len(shap_df) != 17645 or shap_df["gene_id"].nunique() != 17645:
                raise RuntimeError("{} target SHAP table has wrong schema/feature coverage".format(disease))
            shap_genes = shap_df["gene_id"].astype(str)
            if "RNA_row_index" not in shap_df.columns:
                raise RuntimeError("{} target SHAP lacks RNA_row_index needed to verify frozen feature order".format(disease))
            shap_indices = pd.to_numeric(shap_df["RNA_row_index"], errors="raise").astype(int)
            if (shap_genes.duplicated().any() or set(shap_genes) != set(frozen_feature_order) or
                    not np.array_equal(shap_indices.values, np.asarray([frozen_feature_index[x] for x in shap_genes], dtype=int))):
                raise RuntimeError("{} target SHAP feature IDs/order do not match gene_manifest_model.csv".format(disease))
            shap_rank = pd.to_numeric(shap_df["rank"], errors="raise").astype(int).values
            if sorted(shap_rank.tolist()) != list(range(1, 17646)):
                raise RuntimeError("{} target SHAP rank is incomplete or duplicated".format(disease))
            if set(shap_df["target_disease"].astype(str)) != {disease}:
                raise RuntimeError("{} SHAP disease identity mismatch".format(disease))
            if (network_manifest.get("target_shap_sha256", "").lower() != _input_sha256(target_shap_path).lower() or
                    network_manifest.get("model_config_sha256", "").lower() != _input_sha256(model_config_path).lower() or
                    network_manifest.get("pcsf", {}).get("b") != 5 or network_manifest.get("pcsf", {}).get("w") != 3 or
                    network_manifest.get("pcsf", {}).get("g") != 5):
                raise RuntimeError("{} frozen network provenance differs from the current model/SHAP or original parameter record".format(disease))
            raw_shap = pd.to_numeric(shap_df["mean_abs_SHAP"], errors="raise").astype(float).values
            expected_shap_cells = int(runtime.get("n_target_train_SHAP_cells", -1))
            observed_shap_cells = pd.to_numeric(shap_df["n_cells"], errors="raise").astype(int).unique().tolist()
            if observed_shap_cells != [expected_shap_cells]:
                raise RuntimeError("{} SHAP explained-cell count differs from runtime: {} vs {}".format(disease, observed_shap_cells, expected_shap_cells))
            if not np.all(np.isfinite(raw_shap)) or np.any(raw_shap < 0):
                raise RuntimeError("{} SHAP contains invalid mean_abs_SHAP".format(disease))
            contribution_groups = sorted(set(shap_df["contribution_group"].astype(str)))
            if contribution_groups != ["disease:" + disease]:
                raise RuntimeError("{} SHAP contribution group mismatch: {}".format(disease, contribution_groups))
            train_cfg = model_cfg.get("training", {})
            expected_epochs = int(runtime.get("auto_epochs", -1))
            best_epoch = int(runtime.get("best_epoch", -1))
            final_auroc = pd.to_numeric(pd.Series([runtime.get("final_eval_AUROC")]), errors="coerce").iloc[0]
            final_auprc = pd.to_numeric(pd.Series([runtime.get("final_eval_AUPRC")]), errors="coerce").iloc[0]
            training_log = read_table(training_log_path)
            if "epoch" not in training_log.columns or training_log.empty:
                raise RuntimeError("{} training log lacks epoch records".format(disease))
            logged_epochs = pd.to_numeric(training_log["epoch"], errors="raise").astype(int)
            if (int(train_cfg.get("max_epochs", -1)) != expected_epochs or
                    sorted(logged_epochs.tolist()) != list(range(1, expected_epochs + 1)) or
                    best_epoch < 1 or best_epoch > expected_epochs or not np.isfinite(final_auroc) or
                    not np.isfinite(final_auprc) or not runtime.get("best_stage")):
                raise RuntimeError("{} runtime/config/training-log/checkpoint epoch or evaluation metadata is inconsistent".format(disease))
            best_log_row = training_log[pd.to_numeric(training_log["epoch"], errors="raise").astype(int) == best_epoch]
            if best_log_row.empty or not {"val_AUROC", "val_AUPRC"}.issubset(set(training_log.columns)):
                raise RuntimeError("{} training log has no validation metrics at the recorded best epoch".format(disease))
            best_log_auroc = float(pd.to_numeric(best_log_row["val_AUROC"], errors="raise").iloc[0])
            best_log_auprc = float(pd.to_numeric(best_log_row["val_AUPRC"], errors="raise").iloc[0])
            if (not math.isclose(best_log_auroc, float(runtime.get("best_val_AUROC")), rel_tol=1e-9, abs_tol=1e-12) or
                    not math.isclose(best_log_auprc, float(runtime.get("best_val_AUPRC")), rel_tol=1e-9, abs_tol=1e-12)):
                raise RuntimeError("{} recorded best validation metrics differ from training log".format(disease))
            scheduled_stage = next((int(item["stage"]) for item in train_cfg.get("stage_schedule", [])
                                    if int(item["epoch_start_1based"]) <= best_epoch <= int(item["epoch_end_1based"])), None)
            if scheduled_stage is None or int(runtime.get("best_stage", -1)) != scheduled_stage:
                raise RuntimeError("{} best stage does not match the frozen epoch schedule".format(disease))
            runtime_status = str(runtime.get("status", ""))
            if not runtime_status.startswith("calculated;"):
                raise RuntimeError("{} runtime does not report a completed model calculation: {}".format(disease, runtime_status))
            with open(epoch_training_path, "r", encoding="utf-8-sig") as handle:
                epoch_records = [json.loads(line) for line in handle if line.strip()]
            epoch_record_numbers = [int(item.get("epoch", -1)) for item in epoch_records]
            epoch_groups = {}
            for item in epoch_records:
                epoch_groups.setdefault(int(item.get("epoch", -1)), []).append(item)
            unique_epoch_records = []
            for epoch_number in sorted(epoch_groups):
                group = epoch_groups[epoch_number]
                signatures = set()
                for item in group:
                    signatures.add((float(item.get("val_AUROC")), float(item.get("val_AUPRC")),
                                    bool(item.get("is_best_checkpoint")), int(item.get("best_epoch_so_far", -1)),
                                    float(item.get("best_val_AUROC_so_far")), float(item.get("best_val_AUPRC_so_far"))))
                if len(signatures) != 1:
                    raise RuntimeError("{} duplicate epoch {} log records disagree on validation/checkpoint fields".format(disease, epoch_number))
                unique_epoch_records.append(group[-1])
            best_records = [item for item in unique_epoch_records if item.get("is_best_checkpoint") is True]
            best_epoch_records = [item for item in best_records if int(item.get("epoch", -1)) == best_epoch]
            if (sorted(set(epoch_record_numbers)) != list(range(1, expected_epochs + 1)) or
                    len(best_epoch_records) != 1 or int(unique_epoch_records[-1].get("best_epoch_so_far", -1)) != best_epoch):
                raise RuntimeError("{} epoch_training.jsonl does not uniquely confirm the recorded best checkpoint".format(disease))
            best_event = best_epoch_records[0]
            if (not math.isclose(float(best_event.get("best_val_AUROC_so_far")), float(runtime.get("best_val_AUROC")), rel_tol=1e-9, abs_tol=1e-12) or
                    not math.isclose(float(best_event.get("best_val_AUPRC_so_far")), float(runtime.get("best_val_AUPRC")), rel_tol=1e-9, abs_tol=1e-12)):
                raise RuntimeError("{} best-checkpoint event differs from runtime validation selection".format(disease))
            ctrl_df = read_table(ctrl_shap_path)
            required_ctrl = {"target_disease", "contribution_group", "n_cells", "gene_id", "mean_abs_SHAP", "rank", "RNA_row_index"}
            if (not required_ctrl.issubset(set(ctrl_df.columns)) or len(ctrl_df) != 17645 or
                    ctrl_df["gene_id"].astype(str).nunique() != 17645 or
                    set(ctrl_df["gene_id"].astype(str)) != set(frozen_feature_order)):
                raise RuntimeError("{} Ctrl SHAP table does not cover the frozen model feature set".format(disease))
            ctrl_indices = pd.to_numeric(ctrl_df["RNA_row_index"], errors="raise").astype(int)
            ctrl_genes = ctrl_df["gene_id"].astype(str)
            if not np.array_equal(ctrl_indices.values, np.asarray([frozen_feature_index[x] for x in ctrl_genes], dtype=int)):
                raise RuntimeError("{} Ctrl SHAP feature order differs from gene_manifest_model.csv".format(disease))
            if set(ctrl_df["target_disease"].astype(str)) != {disease} or set(ctrl_df["contribution_group"].astype(str)) != {"disease:Ctrl"}:
                raise RuntimeError("{} Ctrl SHAP disease/group identity mismatch".format(disease))
            ctrl_n = pd.to_numeric(ctrl_df["n_cells"], errors="raise").astype(int).unique().tolist()
            if ctrl_n != [int(runtime.get("n_Ctrl_train_SHAP_cells", -1))]:
                raise RuntimeError("{} Ctrl SHAP cell count differs from runtime".format(disease))
            ctrl_rank = pd.to_numeric(ctrl_df["rank"], errors="raise").astype(int).values
            if sorted(ctrl_rank.tolist()) != list(range(1, 17646)):
                raise RuntimeError("{} Ctrl SHAP rank is incomplete or duplicated".format(disease))
            prize_df = read_table(prize_path)
            if len(prize_df) != 12512 or prize_df["Entrez"].astype(str).nunique() != 12512:
                raise RuntimeError("{} gene prize table does not cover the fixed PPI node set".format(disease))
            prize_df = prize_df.set_index(prize_df["Entrez"].astype(str), drop=False)
            if set(prize_df.index) != ppi_ids:
                raise RuntimeError("{} prize table node IDs differ from frozen PPI".format(disease))
            gene_audit = read_table(audit_path)
            required_audit = {"gene_id", "mean_abs_SHAP", "Entrez", "in_qualified_PPI"}
            if (not required_audit.issubset(set(gene_audit.columns)) or len(gene_audit) != 17645 or
                    gene_audit["gene_id"].astype(str).nunique() != 17645):
                raise RuntimeError("{} full gene audit does not cover the unique frozen model feature space".format(disease))
            if set(gene_audit["gene_id"].astype(str)) != set(shap_genes):
                raise RuntimeError("{} SHAP and mapping audit gene sets differ".format(disease))
            audit_shap = pd.Series(pd.to_numeric(gene_audit["mean_abs_SHAP"], errors="raise").astype(float).values,
                                   index=gene_audit["gene_id"].astype(str))
            shap_by_gene = pd.Series(pd.to_numeric(shap_df["mean_abs_SHAP"], errors="raise").astype(float).values,
                                     index=shap_genes)
            if not np.allclose(audit_shap.loc[shap_by_gene.index].values, shap_by_gene.values,
                               rtol=1e-10, atol=1e-14):
                raise RuntimeError("{} full gene audit SHAP values differ from frozen target ranking".format(disease))
            mapped_rows = gene_audit[(gene_audit["Entrez"].astype(str) != "") &
                                     (gene_audit["in_qualified_PPI"].astype(str).str.lower() == "true")].copy()
            mapped_rows["Entrez"] = mapped_rows["Entrez"].astype(str)
            mapped_rows["mean_abs_SHAP"] = pd.to_numeric(mapped_rows["mean_abs_SHAP"], errors="raise").astype(float)
            mapped_max = mapped_rows.groupby("Entrez")["mean_abs_SHAP"].max()
            shap_max = float(np.max(raw_shap))
            if len(mapped_max) != 11968 or shap_max <= 0:
                raise RuntimeError("{} mapping count or full-model reward denominator is invalid".format(disease))
            expected_evidence_ids = set(mapped_max.index.astype(str))
            observed_evidence_ids = set(prize_df.index[pd.to_numeric(prize_df["mean_abs_SHAP"], errors="coerce").notna()])
            if expected_evidence_ids != observed_evidence_ids:
                raise RuntimeError("{} mapped Entrez evidence set differs between audit and prize table".format(disease))
            for eid, expected_raw in mapped_max.items():
                actual_raw = float(prize_df.loc[eid, "mean_abs_SHAP"])
                actual_prize = float(prize_df.loc[eid, "prize"])
                if not math.isclose(actual_raw, float(expected_raw), rel_tol=1e-10, abs_tol=1e-14):
                    raise RuntimeError("{} raw SHAP aggregation differs at Entrez {}".format(disease, eid))
                if not math.isclose(actual_prize, float(expected_raw)/shap_max, rel_tol=1e-9, abs_tol=1e-13):
                    raise RuntimeError("{} prize is not full-model-max normalized exactly once at Entrez {}".format(disease, eid))
            raw = pd.to_numeric(prize_df["mean_abs_SHAP"], errors="coerce").values.astype(float)
            rewards = pd.to_numeric(prize_df["prize"], errors="raise").values.astype(float)
            evidence = np.isfinite(raw)
            if not np.all(np.isfinite(rewards)) or np.any(rewards < 0):
                raise RuntimeError("{} prizes invalid".format(disease))
            if np.any((~evidence) & (rewards != 0)):
                raise RuntimeError("{} no-evidence connectors have nonzero reward".format(disease))
            if np.any(evidence & (rewards == 0) & (raw != 0)):
                raise RuntimeError("{} measured nonzero SHAP unexpectedly has zero prize".format(disease))
            if int(evidence.sum()) != 11968:
                raise RuntimeError("{} evidence node count is not 11968".format(disease))
            checkpoints = model_cfg.get("checkpoint_files", [])
            best_name = "best_validation_checkpoint.pt"
            if best_name not in checkpoints:
                raise RuntimeError("{} config has no declared best checkpoint".format(disease))
            best_checkpoint = os.path.join(model_dir, best_name)
            if not os.path.isfile(best_checkpoint) or os.path.getsize(best_checkpoint) == 0:
                raise RuntimeError("{} best checkpoint missing/empty".format(disease))
            best_checkpoint_sha256 = _input_sha256(best_checkpoint)
            disease_data[disease] = {"prize_df": prize_df, "raw": raw, "rewards": rewards,
                                     "evidence": evidence, "runtime": runtime, "model_cfg": model_cfg,
                                     "target_shap": target_shap_path, "audit": gene_audit,
                                     "network_manifest": network_manifest,
                                     "feature_order_sha256": feature_order_sha256,
                                     "formal_nodes": set(formal_nodes_df["Entrez"].astype(str)),
                                     "formal_edges": {"{}|{}".format(*sorted((str(a),str(b)))) for a,b in zip(formal_edges_df[formal_edge_a],formal_edges_df[formal_edge_b])},
                                     "full_shap_total": float(np.sum(raw_shap))}
            input_files.append(network_manifest_path)
            input_files.append(best_checkpoint)
            model_rows.append({"Disease": disease, "model_config_path": model_config_path,
                "model_config_sha256": _input_sha256(model_config_path), "best_checkpoint_path": best_checkpoint,
                "best_checkpoint_sha256": best_checkpoint_sha256, "training_seed": 42,
                "positive_class": disease, "negative_class": "Ctrl",
                "other_diseases_enter_model": False, "model_features": 17645,
                "target_SHAP_path": target_shap_path, "target_SHAP_sha256": _input_sha256(target_shap_path),
                "frozen_feature_order_sha256": feature_order_sha256,
                "training_log_path": training_log_path, "training_log_sha256": _input_sha256(training_log_path),
            "epoch_training_path": epoch_training_path, "epoch_training_sha256": _input_sha256(epoch_training_path),
                "runtime_path": runtime_path, "runtime_sha256": _input_sha256(runtime_path),
                "training_epochs_expected_and_logged": expected_epochs, "best_epoch": best_epoch,
                "ctrl_SHAP_path": ctrl_shap_path, "ctrl_SHAP_sha256": _input_sha256(ctrl_shap_path),
                "SHAP_aggregation": "mean_abs_SHAP across all target-disease training cells; full 17,645-gene ranking",
                "SHAP_target_cells": runtime.get("n_target_train_SHAP_cells", ""),
                "SHAP_ctrl_cells": runtime.get("n_Ctrl_train_SHAP_cells", ""),
                "n_mapped_model_genes_in_PPI": int(evidence.sum()),
                "mapped_model_feature_fraction": float(evidence.sum())/17645.0,
                "model_features_not_mapped_to_PPI": int(17645-evidence.sum()),
                "n_positive_reward_nodes": int(np.sum(rewards > 0)),
                "mapped_SHAP_mass_fraction": float(disease_audit.loc[disease_audit.Disease == disease, "retained_SHAP_mass_fraction"].iloc[0]),
                "model_final_eval_AUROC": runtime.get("final_eval_AUROC", ""),
                "model_final_eval_AUPRC": runtime.get("final_eval_AUPRC", ""),
                "model_status_note": runtime.get("status", "")})
        except Exception as exc:
            disease_data.pop(disease, None)
            model_rows[:] = [row for row in model_rows if row.get("Disease") != disease]
            blocked_disease_rows.append({"Disease": disease, "status": "blocked_input_identity", "error": repr(exc)})
    ACTIVE_DISEASES[:] = [d for d in DISEASES if d in disease_data]
    write_rows(os.path.join(RUN, "inputs", "blocked_disease_inputs.csv"), blocked_disease_rows, ["Disease", "status", "error"])
    for path in input_files:
        input_rows.append({"path": path, "sha256": _input_sha256(path), "size_bytes": os.path.getsize(path),
                           "last_modified": datetime.fromtimestamp(os.path.getmtime(path)).astimezone().isoformat(timespec="seconds")})
    write_rows(os.path.join(RUN, "inputs", "input_sha256_manifest.csv"), input_rows,
               ["path", "sha256", "size_bytes", "last_modified"])
    write_rows(os.path.join(RUN, "inputs", "input_sha256_provenance.csv"), list(INPUT_HASH_PROVENANCE.values()))
    write_rows(os.path.join(RUN, "inputs", "disease_model_identity_audit.csv"), model_rows)
    return ppi_df, node_df, disease_data


def build_graph(disease, g, ppi_df, disease_data):
    interactome = pd.DataFrame({"protein1": "ENTREZ:" + ppi_df["Protein_A_Entrez_ID"].astype(str),
                               "protein2": "ENTREZ:" + ppi_df["Protein_B_Entrez_ID"].astype(str),
                               "cost": np.ones(len(ppi_df), dtype=float)})
    params = {"b": 1, "w": 3, "g": int(g), "edge_noise": 0, "dummy_mode": "terminals",
              "seed": 42, "skip_checks": False}
    graph = oi_graph.Graph(io.StringIO(interactome.to_csv(sep="\t", index=False)), params)
    nodes = [str(x).split(":", 1)[1] for x in graph.nodes]
    if len(nodes) != 12512 or len(set(nodes)) != 12512:
        raise RuntimeError("{} OI node order/size changed".format(disease))
    ids = set(disease_data[disease]["prize_df"].index)
    if set(nodes) != ids:
        raise RuntimeError("{} OI node identifiers differ from SHAP prize table".format(disease))
    index_by_id = {eid: ix for ix, eid in enumerate(nodes)}
    prize_df = disease_data[disease]["prize_df"]
    bare = np.asarray([float(prize_df.loc[eid, "prize"]) for eid in nodes], dtype=np.float64)
    raw = np.asarray([float(prize_df.loc[eid, "mean_abs_SHAP"]) if str(prize_df.loc[eid, "mean_abs_SHAP"]) not in ("", "nan", "NaN") else np.nan for eid in nodes], dtype=np.float64)
    evidence = np.isfinite(raw)
    if len(graph.edges) != 83065 or not np.array_equal(np.asarray(graph.edges, dtype=np.int64), _expected_edge_indices(ppi_df, nodes)):
        raise RuntimeError("{} OI edge ordering differs from frozen PPI edge table".format(disease))
    degree = np.asarray(graph.node_degrees, dtype=np.int64)
    top_threshold = float(np.quantile(degree, 0.99))
    uv=np.asarray(graph.edges,dtype=np.int64); du=degree[uv[:,0]].astype(float);dv=degree[uv[:,1]].astype(float)
    denom=(len(nodes)-du-1.0)*(len(nodes)-dv-1.0)+du*dv
    penalties=(10.0**int(g))*du*dv/denom
    if not np.allclose(np.asarray(graph.costs,dtype=float),BASE_COST+penalties,rtol=1e-10,atol=1e-12):
        raise RuntimeError("{} OI initial edge costs differ from the frozen hub penalty at g={}".format(disease,g))
    return {"graph": graph, "nodes": nodes, "index_by_id": index_by_id, "bare": bare,
            "raw": raw, "evidence": evidence, "degree": degree, "top1pct_threshold": top_threshold,
            "base_costs": np.ones(len(graph.edges), dtype=np.float64),"edge_penalties":penalties}


def _expected_edge_indices(ppi_df, nodes):
    node_index = {eid: i for i, eid in enumerate(nodes)}
    return np.asarray([(node_index[str(a)], node_index[str(b)]) for a, b in zip(
        ppi_df["Protein_A_Entrez_ID"], ppi_df["Protein_B_Entrez_ID"])], dtype=np.int64)


def set_graph_state(ctx, b, w, g, prizes=None, costs=None):
    graph = ctx["graph"]
    graph._reset_hyperparameters(params={"b": float(b), "w": int(w), "g": int(g), "edge_noise": 0,
        "dummy_mode": "terminals", "seed": 42, "skip_checks": False})
    bare = ctx["bare"] if prizes is None else np.asarray(prizes, dtype=np.float64)
    graph.bare_prizes = bare.copy()
    graph.prizes = graph.bare_prizes * float(b)
    graph.terminals = np.flatnonzero(graph.prizes > 0).astype(np.int64)
    if costs is not None:
        graph.costs = np.asarray(costs, dtype=np.float64).copy()
    else:
        uv = np.asarray(graph.edges, dtype=np.int64)
        deg = np.asarray(graph.node_degrees, dtype=np.float64)
        n_nodes = len(graph.nodes)
        du, dv = deg[uv[:, 0]], deg[uv[:, 1]]
        denom = (n_nodes - du - 1.0) * (n_nodes - dv - 1.0) + du * dv
        penalty = (10.0 ** int(g)) * du * dv / denom
        expected_cost = BASE_COST + penalty
        if not np.all(np.isfinite(expected_cost)) or not np.allclose(
                np.asarray(graph.costs, dtype=np.float64), expected_cost, rtol=1e-10, atol=1e-12):
            raise RuntimeError("OI effective edge costs do not match the frozen hub-penalty formula for g={}".format(g))
        graph.edge_penalties = penalty
        ctx["edge_penalties"] = penalty.copy()
    return graph


def call_solver(graph, pruning=PRUNING):
    captured = []
    original = oi_graph.pcst_fast
    def capture(*args, **kwargs):
        result = original(*args, **kwargs)
        captured.append((args, kwargs, result))
        return result
    oi_graph.pcst_fast = capture
    try:
        stripped_vertices, stripped_edges = graph.pcsf(pruning=pruning, verbosity_level=0)
    finally:
        oi_graph.pcst_fast = original
    if len(captured) != 1:
        raise RuntimeError("Expected one pcst_fast call, captured {}".format(len(captured)))
    args, kwargs, raw = captured[0]
    actual_uv = np.asarray(args[0], dtype=np.int64)
    actual_prizes = np.asarray(args[1], dtype=np.float64)
    actual_costs = np.asarray(args[2], dtype=np.float64)
    root = int(args[3])
    n_real = len(graph.nodes)
    expected_uv = np.vstack([np.asarray(graph.edges, dtype=np.int64),
                              np.column_stack([np.full(len(graph.terminals), n_real, dtype=np.int64),
                                               np.asarray(graph.terminals, dtype=np.int64)])])
    expected_prizes = np.concatenate([np.asarray(graph.prizes, dtype=np.float64), np.asarray([0.0])])
    expected_costs = np.concatenate([np.asarray(graph.costs, dtype=np.float64),
                                     np.full(len(graph.terminals), float(graph.params.w), dtype=np.float64)])
    if (root != n_real or int(args[4]) != 1 or not np.array_equal(actual_uv, expected_uv) or
            not np.allclose(actual_prizes, expected_prizes, rtol=0, atol=0) or
            not np.allclose(actual_costs, expected_costs, rtol=1e-12, atol=1e-12)):
        raise RuntimeError("Captured pcst_fast arrays/root/cluster count differ from independently reconstructed frozen input")
    raw_vertices = np.asarray(raw[0], dtype=np.int64)
    raw_edges = np.asarray(raw[1], dtype=np.int64)
    expected_vertices = np.asarray([v for v in raw_vertices if int(v) != root], dtype=np.int64)
    expected_edges = np.asarray([e for e in raw_edges if int(e) < len(graph.interactome_dataframe)], dtype=np.int64)
    if not np.array_equal(np.sort(expected_vertices), np.sort(np.asarray(stripped_vertices, dtype=np.int64))) or not np.array_equal(np.sort(expected_edges), np.sort(np.asarray(stripped_edges, dtype=np.int64))):
        raise RuntimeError("OI root stripping differs from captured raw solver output")
    return {"raw_vertices": raw_vertices, "raw_edges": raw_edges,
            "stripped_vertices": np.asarray(stripped_vertices, dtype=np.int64),
            "stripped_edges": np.asarray(stripped_edges, dtype=np.int64),
            "call": {"n_nodes": len(args[1]), "n_edges": len(args[0]), "root": int(args[3]),
                     "num_clusters": int(args[4]), "pruning": str(args[5]), "verbosity": int(args[6]),
                     "vertices_dtype": str(np.asarray(args[0]).dtype), "prizes_dtype": str(np.asarray(args[1]).dtype),
                     "costs_dtype": str(np.asarray(args[2]).dtype)},
            "solver_edges": actual_uv, "solver_prizes": actual_prizes,
            "solver_costs": actual_costs, "solver_root": root,
            "solver_num_clusters": int(args[4]), "solver_kwargs": kwargs}


def candidate_audit(graph, ctx, raw_vertices, raw_edges):
    n = len(graph.nodes)
    m = len(graph.edges)
    root = n
    terminal = np.asarray(graph.terminals, dtype=np.int64)
    prizes = np.asarray(graph.prizes, dtype=np.float64)
    costs = np.asarray(graph.costs, dtype=np.float64)
    root_costs = np.full(len(terminal), float(graph.params.w), dtype=np.float64)
    aug_uv = np.vstack([np.asarray(graph.edges, dtype=np.int64),
                        np.column_stack([np.full(len(terminal), root, dtype=np.int64), terminal])])
    aug_costs = np.concatenate([costs, root_costs])
    terminal_pos = {int(node): i for i, node in enumerate(terminal)}
    candidates = [{"solution_type": "raw_pcst_fast", "vertices": list(map(int, raw_vertices)), "edges": list(map(int, raw_edges))}]
    candidates.append({"solution_type": "root_only_empty", "vertices": [root], "edges": []})
    if len(terminal):
        # Every legal singleton uses exactly one root edge at cost w; its objective
        # is minimized by the terminal with the largest current prize.
        ix = int(terminal[np.argmax(prizes[terminal])])
        candidates.append({"solution_type": "best_single_node", "vertices": [root, ix],
                           "edges": [m + terminal_pos[ix]]})
    else:
        candidates.append({"solution_type": "best_single_node", "vertices": [root], "edges": []})
    uv = np.asarray(graph.edges, dtype=np.int64)
    terminal_mask = np.zeros(n, dtype=bool); terminal_mask[terminal] = True
    eligible = terminal_mask[uv[:, 0]] | terminal_mask[uv[:, 1]]
    eligible_ix = np.flatnonzero(eligible)
    if len(eligible_ix):
        possible = uv[eligible_ix]
        candidate_values = (costs[eligible_ix] + float(graph.params.w) + float(np.sum(prizes))
                            - prizes[possible[:, 0]] - prizes[possible[:, 1]])
        local_best = int(np.argmin(candidate_values))
        edge_ix = int(eligible_ix[local_best]); u, v = map(int, uv[edge_ix])
        attach = u if terminal_mask[u] else v
        best_edge = {"solution_type": "best_single_real_PPI_edge", "vertices": [root, u, v],
                     "edges": [edge_ix, m + terminal_pos[attach]]}
        candidates.append(best_edge)
    else:
        candidates.append({"solution_type": "best_single_real_PPI_edge", "vertices": [root], "edges": []})
    star_ix = [int(ix) for ix in terminal if prizes[int(ix)] > float(graph.params.w)]
    candidates.append({"solution_type": "root_star_prize_gt_w", "vertices": [root] + star_ix,
                       "edges": [m + terminal_pos[ix] for ix in star_ix]})
    rows = []
    for cand in candidates:
        vtx = [int(x) for x in cand["vertices"]]
        edg = [int(x) for x in cand["edges"]]
        valid, why = _validate_rooted(vtx, edg, aug_uv, root)
        obj = _objective(vtx, edg, prizes, aug_costs, m, n) if valid else None
        chosen_real = [ix for ix in vtx if ix < n]
        chosen_cost = float(np.sum(aug_costs[np.asarray(edg, dtype=np.int64)])) if valid and edg else 0.0
        chosen_reward = float(np.sum(prizes[chosen_real])) if chosen_real else 0.0
        rows.append(dict(cand, feasible=valid, feasibility_reason=why, objective=obj,
                         selected_reward=chosen_reward, unselected_reward=float(np.sum(prizes))-chosen_reward,
                         selected_edge_cost=chosen_cost,
                         raw_vertices="|".join(map(str, vtx)), raw_edges="|".join(map(str, edg))))
    legal = [r for r in rows if r["feasible"]]
    best_obj = min(r["objective"] for r in legal)
    best = next((r for r in legal if r["solution_type"] == "raw_pcst_fast" and abs(r["objective"] - best_obj) <= 1e-12), None)
    if best is None:
        best = min(legal, key=lambda r: (r["objective"], r["solution_type"]))
    return rows, best, aug_uv, aug_costs


def _objective(vertices, edges, prizes, costs, n_real_edges, n_real_nodes):
    chosen = set(map(int, vertices))
    selected_nodes = [ix for ix in chosen if ix < n_real_nodes]
    selected_prize = float(np.sum(prizes[selected_nodes])) if selected_nodes else 0.0
    total_prize = float(np.sum(prizes))
    selected_cost = float(np.sum(costs[np.asarray(edges, dtype=np.int64)])) if edges else 0.0
    return selected_cost + total_prize - selected_prize


def _validate_rooted(vertices, edges, aug_uv, root):
    vset = set(map(int, vertices)); eset = list(map(int, edges))
    if root not in vset: return False, "root_missing"
    if len(vset) != len(vertices): return False, "duplicate_vertex"
    if len(eset) != len(set(eset)): return False, "duplicate_edge"
    if len(eset) != max(0, len(vset) - 1): return False, "not_tree_edge_count"
    adj = {v: [] for v in vset}
    for ei in eset:
        if ei < 0 or ei >= len(aug_uv): return False, "edge_out_of_bounds"
        u, v = map(int, aug_uv[ei])
        if u not in vset or v not in vset: return False, "edge_endpoint_missing"
        adj[u].append(v); adj[v].append(u)
    seen = set(); stack = [root]
    while stack:
        v = stack.pop()
        if v not in seen:
            seen.add(v); stack.extend(adj[v])
    return (True, "") if seen == vset else (False, "not_root_connected")


def solution_metrics(vertices, edges, graph, ctx):
    n = len(graph.nodes); m = len(graph.edges)
    real_nodes = sorted(int(x) for x in vertices if int(x) < n)
    real_edges = sorted(int(x) for x in edges if int(x) < m)
    root_edges = [int(x) - m for x in edges if int(x) >= m]
    adj = {ix: set() for ix in real_nodes}
    for ei in real_edges:
        u, v = map(int, graph.edges[ei]); adj[u].add(v); adj[v].add(u)
    comps = []
    remaining = set(real_nodes)
    while remaining:
        first = min(remaining); stack = [first]; seen = set()
        while stack:
            u = stack.pop()
            if u not in seen:
                seen.add(u); stack.extend(adj[u] - seen)
        remaining -= seen; comps.append(seen)
    comp_sizes = [len(x) for x in comps]
    evidence = [ix for ix in real_nodes if ctx["evidence"][ix]]
    measured_zero = [ix for ix in evidence if ctx["raw"][ix] == 0]
    connectors = [ix for ix in real_nodes if not ctx["evidence"][ix]]
    raw_den = float(np.nansum(ctx["raw"][ctx["evidence"]]))
    raw_ret = float(np.nansum(ctx["raw"][evidence]))
    degrees = ctx["degree"][real_nodes] if real_nodes else np.asarray([], dtype=float)
    top_nodes = set(np.flatnonzero(ctx["degree"] >= ctx["top1pct_threshold"]).tolist())
    top_fraction = float(sum(ix in top_nodes for ix in real_nodes)) / len(real_nodes) if real_nodes else None
    incident_edges = sum(1 for ei in real_edges if any(int(u) in top_nodes for u in graph.edges[ei]))
    penalties = np.asarray(getattr(graph,"edge_penalties",np.zeros(len(graph.edges))),dtype=float)
    effective_costs = np.asarray(graph.costs,dtype=float)
    all_degree = np.asarray(ctx["degree"],dtype=float)
    return {"real_node_count": len(real_nodes), "model_evidence_node_count": len(evidence),
        "measured_zero_SHAP_evidence_count": len(measured_zero), "zero_reward_no_evidence_connector_count": len(connectors),
        "connector_fraction": float(len(connectors))/len(real_nodes) if real_nodes else None,
        "real_skeleton_edge_count": len(real_edges), "component_count": len(comps),
        "isolated_node_count": sum(1 for x in comp_sizes if x == 1),
        "largest_component_node_count": max(comp_sizes) if comp_sizes else 0,
        "largest_component_fraction": max(comp_sizes)/float(len(real_nodes)) if real_nodes else None,
        "retained_SHAP_mass_fraction_mapped_PPI": raw_ret/raw_den if raw_den else None,
        "retained_SHAP_mass_fraction_full_model": raw_ret/ctx["full_shap_total"] if ctx.get("full_shap_total") else None,
        "selected_SHAP_raw_sum": raw_ret, "mapped_SHAP_raw_sum": raw_den,
        "selected_degree_mean": float(np.mean(degrees)) if len(degrees) else None,
        "selected_degree_median": float(np.median(degrees)) if len(degrees) else None,
        "selected_degree_q90": float(np.quantile(degrees, .9)) if len(degrees) else None,
        "selected_degree_max": int(np.max(degrees)) if len(degrees) else None,
        "top1pct_degree_threshold": ctx["top1pct_threshold"],
        "selected_node_fraction_top1pct_degree": top_fraction,
        "skeleton_edge_fraction_incident_top1pct_degree": incident_edges/float(len(real_edges)) if real_edges else None,
        "background_degree_mean":float(np.mean(all_degree)),"background_degree_median":float(np.median(all_degree)),
        "background_degree_q90":float(np.quantile(all_degree,.9)),"background_degree_max":int(np.max(all_degree)),
        "background_effective_edge_cost_mean":float(np.mean(effective_costs)),
        "background_effective_edge_cost_median":float(np.median(effective_costs)),
        "background_effective_edge_cost_q90":float(np.quantile(effective_costs,.9)),
        "background_effective_edge_cost_max":float(np.max(effective_costs)),
        "background_hub_penalty_mean":float(np.mean(penalties)),
        "background_hub_penalty_median":float(np.median(penalties)),
        "background_hub_penalty_q90":float(np.quantile(penalties,.9)),
        "background_hub_penalty_max":float(np.max(penalties)),
        "component_sizes": "|".join(map(str, sorted(comp_sizes, reverse=True))) }


def _validated_prior_base(disease, b, w, g, current_input_manifest_hash, expected_terminal_count):
    if not (math.isclose(float(b), 5.0, rel_tol=0.0, abs_tol=1e-12) and int(w) in (3, 6) and int(g) in (2, 3, 5)):
        return None
    pid = param_id(b, w, g)
    source_dir = os.path.join(REUSE_ROOT, "grid", disease, pid)
    manifest_path = os.path.join(source_dir, "run_manifest.json")
    if not os.path.isfile(manifest_path):
        raise RuntimeError("Expected reusable prior base result is missing: {}".format(manifest_path))
    with open(manifest_path, "r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    if (manifest.get("status") != "complete" or manifest.get("config_id") != "Fig6_disease_pcsf_grid_v1" or
            manifest.get("Disease") != disease or
            not math.isclose(float(manifest.get("b")), float(b), rel_tol=0.0, abs_tol=1e-12) or
            int(manifest.get("w")) != int(w) or int(manifest.get("g")) != int(g) or
            manifest.get("input_manifest_sha256") != current_input_manifest_hash or
            manifest.get("pruning") != PRUNING or manifest.get("dummy_mode") != "terminals"):
        raise RuntimeError("Prior base result identity does not match frozen inputs/config: {}".format(manifest_path))
    prior_config_path = os.path.join(REUSE_ROOT, "config", "run_config.json")
    with open(prior_config_path, "r", encoding="utf-8") as handle:
        prior_config = json.load(handle)
    solver_cfg = prior_config.get("solver", {})
    grid_cfg = prior_config.get("grid", {})
    if (sha256(prior_config_path) != manifest.get("run_config_sha256") or
            prior_config.get("execution_script_sha256") != manifest.get("execution_script_sha256") or
            sha256(prior_config.get("execution_script_path", "")) != manifest.get("execution_script_sha256") or
            solver_cfg.get("implementation") != "OmicsIntegrator2 2.3.10 + pcst-fast 1.0.7" or
            solver_cfg.get("pruning") != PRUNING or solver_cfg.get("dummy_mode") != "terminals" or
            int(solver_cfg.get("seed", -1)) != 42 or float(solver_cfg.get("edge_noise", -1)) != 0.0 or
            solver_cfg.get("hub_penalty") != "10^g*deg(u)*deg(v)/((N-deg(u)-1)*(N-deg(v)-1)+deg(u)*deg(v))" or
            solver_cfg.get("objective") != "sum(selected real/root edge costs)+sum(unselected real-node scaled prizes); OI's w cost is represented by root edges exactly once" or
            prior_config.get("reward") != "target training-cell mean_abs_SHAP divided by disease-wide maximum over all model genes before PPI mapping; mapped Entrez max rule reused; b applied once; no PPI renormalization" or
            prior_config.get("ppi", {}).get("filter") != "number>=2" or
            float(prior_config.get("ppi", {}).get("base_cost")) != BASE_COST or
            grid_cfg.get("b") != [2, 5, 10] or grid_cfg.get("w") != [1, 3, 6] or grid_cfg.get("g") != [2, 3, 5, 6]):
        raise RuntimeError("Prior full solver/reward/PPI configuration differs: {}".format(prior_config_path))
    call = manifest.get("solver_call", {})
    # OI adds one synthetic root edge per positive-reward terminal. The count
    # therefore varies by disease; validate against the frozen input rather
    # than an incorrect global edge-count constant.
    expected_solver_edges = 83065 + int(expected_terminal_count)
    if (call.get("pruning") != PRUNING or int(call.get("root", -1)) != 12512 or
            int(call.get("num_clusters", -1)) != 1 or int(call.get("n_nodes", -1)) != 12513 or
            int(call.get("n_edges", -1)) != expected_solver_edges or
            call.get("vertices_dtype") != "int64" or call.get("prizes_dtype") != "float64" or
            call.get("costs_dtype") != "float64"):
        raise RuntimeError("Prior per-run pcst_fast call semantics differ: {}".format(manifest_path))
    expected_outputs = {os.path.basename(x["path"]): x["sha256"] for x in manifest.get("outputs", [])}
    required = {"adopted_nodes.csv", "adopted_forest_edges.tsv", "adopted_component_membership.csv",
                "adopted_virtual_root_edges.tsv", "raw_solver_vertices.csv", "raw_solver_edges.csv",
                "candidate_objective_comparison.csv", "solver_arrays.npz"}
    if not required.issubset(expected_outputs):
        raise RuntimeError("Prior reuse manifest lacks required outputs: {}".format(manifest_path))
    for name, expected_hash in expected_outputs.items():
        path = os.path.join(source_dir, name)
        if not os.path.isfile(path) or sha256(path) != expected_hash:
            raise RuntimeError("Prior reused output hash mismatch: {}".format(path))
    metrics = manifest.get("metrics")
    if not isinstance(metrics, dict) or metrics.get("Disease") != disease or metrics.get("parameter_id") != pid:
        raise RuntimeError("Prior metrics payload is incomplete: {}".format(manifest_path))
    comparison_path = os.path.join(source_dir, "candidate_objective_comparison.csv")
    candidates = read_table(comparison_path).to_dict("records")
    reuse = {"Disease": disease, "parameter_id": pid, "b": float(b), "w": int(w), "g": int(g),
        "reuse_status": "verified_reuse", "source_run_id": "20261020_disease_parameter_adaptation_08",
        "source_directory": source_dir, "source_manifest_path": manifest_path,
        "source_manifest_sha256": sha256(manifest_path),
        "source_nodes_sha256": sha256(os.path.join(source_dir, "adopted_nodes.csv")),
        "source_edges_sha256": sha256(os.path.join(source_dir, "adopted_forest_edges.tsv")),
        "source_arrays_sha256": sha256(os.path.join(source_dir, "solver_arrays.npz"))}
    return metrics, candidates, reuse


def run_base(ppi_df, node_df, disease_data):
    grid_rows = []
    comparison_rows = []
    status_rows = []
    reuse_rows = []
    started = time.time()
    config_hash = sha256(os.path.join(RUN, "config", "run_config.json"))
    input_manifest_path = os.path.join(RUN, "inputs", "input_sha256_manifest.csv")
    input_manifest_hash = sha256(input_manifest_path)
    script_hash = sha256(os.path.abspath(__file__))
    for disease in ACTIVE_DISEASES:
        for g in G_VALUES:
            ctx = build_graph(disease, g, ppi_df, disease_data)
            for b, w in itertools.product(B_VALUES, W_VALUES):
                pid = param_id(b, w, g)
                out_dir = os.path.join(RUN, "grid", disease, pid)
                os.makedirs(out_dir, exist_ok=True)
                manifest_path = os.path.join(out_dir, "run_manifest.json")
                if not os.path.isfile(manifest_path):
                    prior = _validated_prior_base(
                        disease, b, w, g, input_manifest_hash,
                        expected_terminal_count=int(np.sum(ctx["bare"] > 0)))
                    if prior is not None:
                        prior_metrics, prior_candidates, reuse_record = prior
                        grid_rows.append(prior_metrics)
                        status_rows.append({"Disease": disease, "parameter_id": pid, "status": "complete_reused_prior_base",
                                            "runtime_seconds": prior_metrics.get("solver_runtime_seconds", ""),
                                            "source_directory": reuse_record["source_directory"]})
                        comparison_rows.extend([{**row, "Disease": disease, "parameter_id": pid} for row in prior_candidates])
                        reuse_rows.append(reuse_record)
                        write_rows(os.path.join(RUN, "summaries", "base_run_status.csv"), status_rows)
                        write_rows(os.path.join(RUN, "summaries", "reuse_registry.csv"), reuse_rows)
                        continue
                if os.path.isfile(manifest_path):
                    with open(manifest_path, "r", encoding="utf-8") as handle:
                        old = json.load(handle)
                    expected_outputs = {os.path.basename(x["path"]): x["sha256"] for x in old.get("outputs", [])}
                    outputs_valid = bool(expected_outputs)
                    for output_name, expected_hash in expected_outputs.items():
                        output_path = os.path.join(out_dir, output_name)
                        if not os.path.isfile(output_path) or sha256(output_path) != expected_hash:
                            outputs_valid = False
                            break
                    if (old.get("status") == "complete" and old.get("config_id") == "Fig6_disease_pcsf_grid_v1" and
                            old.get("execution_script_sha256") == script_hash and
                            old.get("run_config_sha256") == config_hash and
                            old.get("input_manifest_sha256") == input_manifest_hash and outputs_valid):
                        grid_rows.append(old["metrics"]); status_rows.append({"Disease":disease,"parameter_id":pid,"status":"complete_reused","runtime_seconds":old.get("runtime_seconds","")})
                        comparison_path=os.path.join(out_dir,"candidate_objective_comparison.csv")
                        if os.path.isfile(comparison_path):
                            comparison_rows.extend([{**r,"Disease":disease,"parameter_id":pid} for r in read_table(comparison_path).to_dict("records")])
                        continue
                    raise RuntimeError("Existing grid output identity or artifact hashes differ; refusing to reuse/mix it: {}".format(manifest_path))
                t0 = time.time()
                graph = set_graph_state(ctx, b, w, g)
                solver_t0 = time.time()
                raw = call_solver(graph)
                solver_runtime = time.time() - solver_t0
                candidates, best, aug_uv, aug_costs = candidate_audit(graph, ctx, raw["raw_vertices"], raw["raw_edges"])
                n = len(graph.nodes); m = len(graph.edges); prizes = np.asarray(graph.prizes, dtype=np.float64)
                metrics = solution_metrics(best["vertices"], best["edges"], graph,
                    dict(ctx, full_shap_total=disease_data[disease]["full_shap_total"]))
                raw_objective = next(x["objective"] for x in candidates if x["solution_type"]=="raw_pcst_fast")
                metrics.update({"Disease": disease, "parameter_id": pid, "b":b, "w":w, "g":g,
                    "raw_solver_real_node_count":len(raw["stripped_vertices"]), "raw_solver_real_edge_count":len(raw["stripped_edges"]),
                    "solver_runtime_seconds":solver_runtime,
                    "raw_objective":raw_objective,
                    "adopted_solution_source":best["solution_type"], "adopted_objective":best["objective"],
                    "empty_objective":float(np.sum(prizes)),
                    "improvement_over_raw":(raw_objective-best["objective"] if raw_objective is not None else None)})
                if (b,w,g)==(5,3,5):
                    raw_node_ids={ctx["nodes"][int(ix)] for ix in raw["stripped_vertices"]}
                    raw_edge_keys={"{}|{}".format(*sorted((ctx["nodes"][int(graph.edges[int(ei)][0])],
                                                            ctx["nodes"][int(graph.edges[int(ei)][1])])))
                                   for ei in raw["stripped_edges"]}
                    metrics.update({"formal_b5w3g5_node_jaccard":jaccard(raw_node_ids,disease_data[disease]["formal_nodes"]),
                        "formal_b5w3g5_backbone_edge_jaccard":jaccard(raw_edge_keys,disease_data[disease]["formal_edges"]),
                        "formal_b5w3g5_exact_node_match":raw_node_ids==disease_data[disease]["formal_nodes"],
                        "formal_b5w3g5_exact_edge_match":raw_edge_keys==disease_data[disease]["formal_edges"]})
                _save_network(out_dir, disease, pid, best, raw, graph, ctx, ppi_df, metrics, candidates)
                elapsed = time.time() - t0
                manifest = {"Disease":disease,"parameter_id":pid,"b":b,"w":w,"g":g,"seed":42,
                    "implementation":"OmicsIntegrator2 2.3.10 / pcst-fast 1.0.7", "pruning":PRUNING,
                    "dummy_mode":"terminals", "status":"complete", "config_id":"Fig6_disease_pcsf_grid_v1",
                    "started_at":stamp(),"runtime_seconds":elapsed,"solver_runtime_seconds":solver_runtime,
                    "solver_call":raw["call"],"metrics":metrics,
                    "execution_script_sha256":script_hash,"run_config_sha256":config_hash,
                    "input_manifest_sha256":input_manifest_hash,
                    "comparison_candidates":candidates,"outputs":_output_hashes(out_dir)}
                write_json(manifest_path, manifest)
                grid_rows.append(metrics)
                status_rows.append({"Disease":disease,"parameter_id":pid,"status":"complete","runtime_seconds":elapsed})
                comparison_rows.extend([{**row,"Disease":disease,"parameter_id":pid} for row in candidates])
                write_rows(os.path.join(RUN,"summaries","base_run_status.csv"),status_rows)
                log("base complete disease={} param={} nodes={} edges={} adopted={} J={:.6f} sec={:.2f}".format(disease,pid,metrics["real_node_count"],metrics["real_skeleton_edge_count"],best["solution_type"],best["objective"],elapsed))
    write_rows(os.path.join(RUN,"summaries","parameter_network_metrics.csv"),grid_rows)
    write_rows(os.path.join(RUN,"summaries","solver_objective_audit.csv"),comparison_rows)
    write_rows(os.path.join(RUN,"summaries","base_run_status.csv"),status_rows)
    write_rows(os.path.join(RUN,"summaries","reuse_registry.csv"),reuse_rows)
    # Neighbor similarities are assembled after all results are loaded through the provenance registry.
    log("base grid complete n_runs={} elapsed_seconds={:.1f}".format(len(grid_rows),time.time()-started))


def _json_clean(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.ndarray):
        return [_json_clean(x) for x in value.tolist()]
    if isinstance(value, dict):
        return {str(k): _json_clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_clean(x) for x in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _control_run_path(disease, pid, control_type, replicate):
    return os.path.join(RUN, "perturbations", disease, pid, control_type,
                        "replicate_{:02d}.json".format(int(replicate)))


def _record_digest(record):
    payload=dict(record);payload.pop("record_payload_sha256",None)
    canonical=json.dumps(_json_clean(payload),ensure_ascii=False,sort_keys=True,
                         separators=(",",":"),allow_nan=False).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def run_one_control(disease, b, w, g, control_type, replicate, ppi_df, gctx,
                    disease_data, degree_bin_map, config_hash, input_manifest_hash,
                    script_hash):
    pid = param_id(b, w, g)
    result_path = _control_run_path(disease, pid, control_type, replicate)
    if os.path.isfile(result_path):
        with open(result_path, "r", encoding="utf-8") as handle:
            old = json.load(handle)
        if (old.get("status") == "complete" and old.get("Disease") == disease and
                old.get("parameter_id") == pid and old.get("control_type") == control_type and
                int(old.get("replicate", -1)) == int(replicate) and
                old.get("run_config_sha256") == config_hash and
                old.get("input_manifest_sha256") == input_manifest_hash and
                old.get("execution_script_sha256") == script_hash and
                old.get("record_payload_sha256") == _record_digest(old)):
            return old
        raise RuntimeError("Perturbation checkpoint identity mismatch; refusing reuse: {}".format(result_path))

    bseed = stable_seed("{}|{}|{}|{}|rep{}".format(RUN_ID, disease, pid, control_type, replicate))
    rng = np.random.RandomState(bseed)
    base_bare = gctx["bare"].copy()
    cost_z = None
    permutation_indices = None
    if control_type == "edge_cost":
        cost_z = rng.normal(0.0, SIGMA, size=len(gctx["base_costs"]))
        multiplier = np.exp(cost_z - (SIGMA ** 2) / 2.0)
        new_costs = gctx["base_costs"] * multiplier + gctx["edge_penalties"]
        altered_bare = base_bare
    elif control_type == "reward_shuffle":
        evidence_ix = np.flatnonzero(gctx["evidence"])
        permutation_indices = rng.permutation(len(evidence_ix)).astype(np.int32)
        altered_bare = base_bare.copy()
        altered_bare[evidence_ix] = base_bare[evidence_ix[permutation_indices]]
        new_costs = gctx["base_costs"] + gctx["edge_penalties"]
    elif control_type == "degree_stratified_reward_shuffle":
        evidence_ix = np.flatnonzero(gctx["evidence"])
        permutation_indices = np.arange(len(evidence_ix), dtype=np.int32)
        for bin_row in degree_bin_map[disease]:
            members = np.asarray(bin_row["indices"], dtype=np.int64)
            positions = np.searchsorted(evidence_ix, members)
            perm = rng.permutation(len(members))
            permutation_indices[positions] = positions[perm]
        altered_bare = base_bare.copy()
        altered_bare[evidence_ix] = base_bare[evidence_ix[permutation_indices]]
        new_costs = gctx["base_costs"] + gctx["edge_penalties"]
    else:
        raise ValueError("Unknown control type: {}".format(control_type))

    graph = set_graph_state(gctx, b, w, g, prizes=altered_bare, costs=new_costs)
    solver_started = time.time()
    raw = call_solver(graph)
    solver_runtime_seconds = time.time() - solver_started
    n, m = len(graph.nodes), len(graph.edges)
    raw_vertices = list(map(int, raw["raw_vertices"]))
    raw_edges = list(map(int, raw["raw_edges"]))
    candidate_rows, adopted, aug_uv, aug_costs = candidate_audit(graph, gctx, raw_vertices, raw_edges)
    if not adopted.get("feasible"):
        raise RuntimeError("Control solver has no legal feasible candidate {}: {}".format(disease, control_type))
    adopted_vertices = list(map(int, adopted["vertices"]))
    adopted_edges = list(map(int, adopted["edges"]))
    valid, why = _validate_rooted(adopted_vertices, adopted_edges, aug_uv, n)
    if not valid:
        raise RuntimeError("Adopted control witness is infeasible {}: {} {}".format(disease, control_type, why))
    objective = float(adopted["objective"])
    metrics = solution_metrics(adopted_vertices, adopted_edges, graph, gctx)
    real_node_ix = sorted(ix for ix in set(adopted_vertices) if ix < n)
    real_edge_ix = sorted(ix for ix in set(adopted_edges) if ix < m)
    raw_real_node_ix = sorted(ix for ix in set(raw_vertices) if ix < n)
    raw_real_edge_ix = sorted(ix for ix in set(raw_edges) if ix < m)
    ref_path = os.path.join(RUN, "grid", disease, pid, "adopted_nodes.csv")
    ref_node_ids = set(read_table(ref_path)["Entrez"].astype(str))
    ref_node_ix = {gctx["index_by_id"][eid] for eid in ref_node_ids}
    ref_edge_path = os.path.join(RUN, "grid", disease, pid, "adopted_forest_edges.tsv")
    ref_edges_df = read_table(ref_edge_path, "\t")
    ref_edge_keys = {"{}|{}".format(*sorted((str(a), str(b))))
                     for a, b in zip(ref_edges_df["Entrez_A"], ref_edges_df["Entrez_B"])}
    observed_node_ids = {gctx["nodes"][ix] for ix in real_node_ix}
    observed_edge_keys = {"{}|{}".format(*sorted((gctx["nodes"][int(gctx["graph"].edges[ix][0])],
                                                    gctx["nodes"][int(gctx["graph"].edges[ix][1])])) )
                          for ix in real_edge_ix}
    selected_cost = float(np.sum(aug_costs[np.asarray(adopted_edges, dtype=np.int64)])) if adopted_edges else 0.0
    total_reward = float(np.sum(np.asarray(raw["solver_prizes"], dtype=np.float64)[:n]))
    selected_reward = float(np.sum(np.asarray(graph.prizes, dtype=np.float64)[real_node_ix])) if real_node_ix else 0.0
    root_terminal_ix = [int(graph.terminals[ix - m]) for ix in adopted_edges if ix >= m]
    components = []
    component_adj = {ix: set() for ix in real_node_ix}
    for ei in real_edge_ix:
        u, v = map(int, graph.edges[ei]); component_adj[u].add(v); component_adj[v].add(u)
    remaining = set(real_node_ix)
    while remaining:
        first = min(remaining); stack = [first]; members = set()
        while stack:
            ix = stack.pop()
            if ix not in members:
                members.add(ix); stack.extend(component_adj[ix] - members)
        remaining -= members; components.append(sorted(members))
    components.sort(key=lambda xs: xs[0])
    component_id = {ix: ci for ci, xs in enumerate(components, 1) for ix in xs}
    if control_type == "degree_stratified_reward_shuffle":
        perm_record = "within frozen degree bins; map stored as evidence-index permutation"
    elif control_type == "reward_shuffle":
        perm_record = "uniform permutation across all mapped model-evidence Entrez nodes, including measured zero SHAP"
    else:
        perm_record = "not applicable"
    record = {"Disease": disease, "parameter_id": pid, "b": float(b), "w": int(w), "g": int(g),
        "control_type": control_type, "replicate": int(replicate), "rng_algorithm": "NumPy RandomState MT19937",
        "numpy_version": np.__version__, "random_seed": int(bseed), "status": "complete",
        "termination_status": "pcst_fast_returned", "is_feasible": True, "feasibility_reason": "",
        "solver_name": "pcst_fast", "solver_version": "1.0.7", "pruning": PRUNING,
        "virtual_root_index": n, "virtual_root_selected": n in set(raw_vertices),
        "raw_solver_vertex_indices": raw_vertices, "raw_solver_edge_indices": raw_edges,
        "selected_raw_vertex_indices": adopted_vertices, "selected_augmented_edge_indices": adopted_edges,
        "adopted_solution_source": adopted["solution_type"],
        "candidate_objective_comparison": candidate_rows,
        "raw_solver_is_feasible": bool(_validate_rooted(raw_vertices, raw_edges,
            np.asarray(raw["solver_edges"], dtype=np.int64), int(raw["solver_root"]))[0]),
        "raw_solver_real_node_indices": raw_real_node_ix, "raw_solver_real_edge_indices": raw_real_edge_ix,
        "selected_real_node_indices": real_node_ix, "selected_real_edge_indices": real_edge_ix,
        "selected_virtual_root_terminal_indices": root_terminal_ix,
        "component_membership": [{"Entrez": gctx["nodes"][ix], "component_id": component_id[ix]} for ix in real_node_ix],
        "real_component_count": len(components), "selected_reward": selected_reward,
        "unselected_reward": total_reward-selected_reward, "selected_augmented_edge_cost": selected_cost,
        "objective_recomputed": objective, "objective_terms": {"selected_edge_cost": selected_cost,
            "unselected_node_reward": total_reward-selected_reward}, "best_bound": None,
        "absolute_gap": None, "relative_gap": None, "solver_runtime_seconds": solver_runtime_seconds,
        "network_metrics": metrics, "node_jaccard_to_unperturbed_adopted": jaccard(real_node_ix, ref_node_ix),
        "backbone_edge_jaccard_to_unperturbed_adopted": jaccard(observed_edge_keys, ref_edge_keys),
        "selected_node_Entrez": [gctx["nodes"][ix] for ix in real_node_ix],
        "selected_edge_keys": sorted(observed_edge_keys), "reward_transform": perm_record,
        "reward_evidence_node_order_sha256": hashlib.sha256("\n".join([gctx["nodes"][ix] for ix in np.flatnonzero(gctx["evidence"])]).encode("utf-8")).hexdigest(),
        "edge_order_sha256": hashlib.sha256("\n".join(["{}|{}".format(*sorted((gctx["nodes"][int(u)],gctx["nodes"][int(v)]))) for u,v in graph.edges]).encode("utf-8")).hexdigest(),
        "cost_noise_seed_replay": ("RandomState(seed).normal(0,0.1,n_edges), multiplier=exp(z-0.1^2/2), applied to base cost before fixed hub penalty" if control_type == "edge_cost" else None),
        "prize_permutation_seed_replay": ("RandomState(seed) permutation over ordered evidence-node array; no-evidence connector prizes fixed at zero" if control_type != "edge_cost" else None),
        "solver_call": raw["call"], "run_config_sha256": config_hash,
        "input_manifest_sha256": input_manifest_hash, "execution_script_sha256": script_hash,
        "completed_at": stamp()}
    record["record_payload_sha256"]=_record_digest(record)
    write_json(result_path, _json_clean(record))
    return record


def run_controls(ppi_df, disease_data):
    config_hash = sha256(os.path.join(RUN, "config", "run_config.json"))
    input_manifest_path = os.path.join(RUN, "inputs", "input_sha256_manifest.csv")
    input_manifest_hash = sha256(input_manifest_path)
    script_hash = sha256(os.path.abspath(__file__))
    bins_by_disease = {}
    bin_rows = []
    if not ACTIVE_DISEASES:
        raise RuntimeError("No disease has passed frozen model/SHAP/PPI input identity checks")
    # The PPI topology and g-specific hub costs are shared by all diseases.
    # Keep one graph template per g and attach each disease's frozen rewards in
    # the verified OI node order; this avoids rebuilding the 83,065-edge graph
    # thousands of times during the 20-repeat rotation.
    graph_templates = {g: build_graph(ACTIVE_DISEASES[0], g, ppi_df, disease_data)
                       for g in G_VALUES}
    for disease in ACTIVE_DISEASES:
        ctx = dict(graph_templates[G_VALUES[0]])
        prize_df = disease_data[disease]["prize_df"]
        ctx["bare"] = np.asarray([float(prize_df.loc[eid, "prize"]) for eid in ctx["nodes"]], dtype=np.float64)
        ctx["raw"] = np.asarray([float(prize_df.loc[eid, "mean_abs_SHAP"]) if str(prize_df.loc[eid, "mean_abs_SHAP"]) not in ("", "nan", "NaN") else np.nan for eid in ctx["nodes"]], dtype=np.float64)
        ctx["evidence"] = np.isfinite(ctx["raw"])
        evidence_ix = np.flatnonzero(ctx["evidence"])
        bins = degree_bins(evidence_ix, ctx["degree"], min_size=100)
        bins_by_disease[disease] = bins
        for bin_id, item in enumerate(bins, 1):
            for ix in item["indices"]:
                bin_rows.append({"Disease": disease, "degree_bin": bin_id,
                    "degree_min": item["degree_min"], "degree_max": item["degree_max"],
                    "Entrez": ctx["nodes"][int(ix)], "background_degree": int(ctx["degree"][int(ix)]),
                    "bin_size": len(item["indices"]), "tail_handling": "merge trailing remainder into prior bin"})
    write_rows(os.path.join(RUN,"config","degree_stratified_reward_shuffle_bins.csv"), bin_rows,
        ["Disease","degree_bin","degree_min","degree_max","Entrez","background_degree","bin_size","tail_handling"])
    total = len(ACTIVE_DISEASES)*len(B_VALUES)*len(W_VALUES)*len(G_VALUES)*len(CTRL_TYPES)*REPEATS
    completed = 0
    stage_start = time.time()
    # Disease coverage is advanced in lock-step: each replicate/type visits all 16 diseases before the next round.
    for replicate in range(1, REPEATS+1):
        for control_type in CTRL_TYPES:
            for disease in ACTIVE_DISEASES:
                for g in G_VALUES:
                    ctx = dict(graph_templates[g])
                    prize_df = disease_data[disease]["prize_df"]
                    ctx["bare"] = np.asarray([float(prize_df.loc[eid, "prize"]) for eid in ctx["nodes"]], dtype=np.float64)
                    ctx["raw"] = np.asarray([float(prize_df.loc[eid, "mean_abs_SHAP"]) if str(prize_df.loc[eid, "mean_abs_SHAP"]) not in ("", "nan", "NaN") else np.nan for eid in ctx["nodes"]], dtype=np.float64)
                    ctx["evidence"] = np.isfinite(ctx["raw"])
                    for b, w in itertools.product(B_VALUES, W_VALUES):
                        record = run_one_control(disease, b, w, g, control_type, replicate,
                            ppi_df, ctx, disease_data, bins_by_disease, config_hash, input_manifest_hash, script_hash)
                        completed += 1
                # persist per-disease checkpoint/status summary after each disease turn
                write_json(os.path.join(RUN,"perturbations",disease,"disease_control_progress.json"),
                    {"Disease":disease,"completed_records_in_current_run":_count_disease_controls(disease),
                     "planned_records":len(B_VALUES)*len(W_VALUES)*len(G_VALUES)*len(CTRL_TYPES)*REPEATS,
                     "updated_at":stamp(),"run_config_sha256":config_hash,"input_manifest_sha256":input_manifest_hash})
            log("control rotation complete replicate={} type={} records_completed={} planned={} elapsed_seconds={:.1f}".format(
                replicate, control_type, completed, total, time.time()-stage_start))
    summarize_controls()


def _count_disease_controls(disease):
    root = os.path.join(RUN, "perturbations", disease)
    count = 0
    if not os.path.isdir(root): return count
    for base, _, files in os.walk(root):
        count += sum(1 for f in files if f.startswith("replicate_") and f.endswith(".json"))
    return count


def summarize_controls():
    rows=[]
    root=os.path.join(RUN,"perturbations")
    for disease in ACTIVE_DISEASES:
        for b,w,g in itertools.product(B_VALUES,W_VALUES,G_VALUES):
            pid=param_id(b,w,g)
            for control_type in CTRL_TYPES:
                for replicate in range(1,REPEATS+1):
                    path=_control_run_path(disease,pid,control_type,replicate)
                    if not os.path.isfile(path):
                        continue
                    with open(path,"r",encoding="utf-8") as handle: row=json.load(handle)
                    row_summary={k:row.get(k) for k in ["Disease","parameter_id","b","w","g","control_type","replicate","status","termination_status","is_feasible","real_component_count","objective_recomputed","selected_reward","unselected_reward","selected_augmented_edge_cost","node_jaccard_to_unperturbed_adopted","backbone_edge_jaccard_to_unperturbed_adopted","random_seed","solver_runtime_seconds","termination_status"]}
                    metrics=row.get("network_metrics",{})
                    for key,value in metrics.items(): row_summary["network_"+key]=value
                    row_summary["selected_node_Entrez"]="|".join(row.get("selected_node_Entrez",[]))
                    row_summary["selected_edge_keys"]=";".join(row.get("selected_edge_keys",[]))
                    row_summary["selected_virtual_root_terminal_indices"]="|".join(map(str,row.get("selected_virtual_root_terminal_indices",[])))
                    rows.append(row_summary)
    fields=sorted({k for row in rows for k in row}) if rows else ["Disease","parameter_id","control_type","replicate","status"]
    write_rows(os.path.join(RUN,"summaries","replicate_network_metrics.csv"),rows,fields)
    # Per disease/parameter/control type node and edge selection frequencies, retaining actual repeat denominator.
    node_freq=[]; edge_freq=[]; controls_summary=[]
    grouped={}
    for row in rows: grouped.setdefault((row["Disease"],row["parameter_id"],row["control_type"]),[]).append(row)
    for key,group in grouped.items():
        disease,pid,control_type=key; observed=len(group)
        node_counts={}; edge_counts={}
        for row in group:
            for eid in str(row.get("selected_node_Entrez","")).split("|"):
                if eid: node_counts[eid]=node_counts.get(eid,0)+1
            for edge in str(row.get("selected_edge_keys","")).split(";"):
                if edge: edge_counts[edge]=edge_counts.get(edge,0)+1
        for eid,count in node_counts.items(): node_freq.append({"Disease":disease,"parameter_id":pid,"control_type":control_type,"Entrez":eid,"selected_count":count,"completed_repeats":observed,"selection_frequency":count/float(observed)})
        for edge,count in edge_counts.items(): edge_freq.append({"Disease":disease,"parameter_id":pid,"control_type":control_type,"edge_key":edge,"selected_count":count,"completed_repeats":observed,"selection_frequency":count/float(observed)})
        vals=[r.get("node_jaccard_to_unperturbed_adopted") for r in group]
        vals=[float(x) for x in vals if x is not None]
        edge_vals=[r.get("backbone_edge_jaccard_to_unperturbed_adopted") for r in group]
        edge_vals=[float(x) for x in edge_vals if x is not None]
        controls_summary.append({"Disease":disease,"parameter_id":pid,"control_type":control_type,"completed_repeats":observed,
            "median_node_jaccard":float(np.median(vals)) if vals else None,"q25_node_jaccard":float(np.quantile(vals,.25)) if vals else None,
            "median_edge_jaccard":float(np.median(edge_vals)) if edge_vals else None,"q25_edge_jaccard":float(np.quantile(edge_vals,.25)) if edge_vals else None})
    write_rows(os.path.join(RUN,"summaries","control_node_selection_frequency.csv"),node_freq,
        ["Disease","parameter_id","control_type","Entrez","selected_count","completed_repeats","selection_frequency"])
    write_rows(os.path.join(RUN,"summaries","control_edge_selection_frequency.csv"),edge_freq,
        ["Disease","parameter_id","control_type","edge_key","selected_count","completed_repeats","selection_frequency"])
    write_rows(os.path.join(RUN,"summaries","control_stability_by_disease_parameter.csv"),controls_summary)
    log("control summaries written records={} groups={}".format(len(rows),len(grouped)))


def summarize_shared_defaults():
    metrics_path=os.path.join(RUN,"summaries","parameter_network_metrics.csv")
    neighbors_path=os.path.join(RUN,"summaries","parameter_neighbor_similarity.csv")
    if not os.path.isfile(metrics_path) or not os.path.isfile(neighbors_path):
        raise RuntimeError("Base-grid summaries are missing; cannot compare shared defaults")
    metrics=read_table(metrics_path)
    neighbors=read_table(neighbors_path)
    for col in ["b","w","g","real_node_count","real_skeleton_edge_count","component_count",
                "isolated_node_count","largest_component_fraction","retained_SHAP_mass_fraction_mapped_PPI",
                "connector_fraction","selected_node_fraction_top1pct_degree",
                "skeleton_edge_fraction_incident_top1pct_degree","adopted_objective","raw_objective"]:
        if col in metrics.columns: metrics[col]=pd.to_numeric(metrics[col],errors="coerce")
    for col in ["node_jaccard","backbone_edge_jaccard"]:
        neighbors[col]=pd.to_numeric(neighbors[col],errors="coerce")
    neighbor_values={}
    for _,row in neighbors.iterrows():
        key=(row["Disease"],row["parameter_from"])
        neighbor_values.setdefault(key,[]).append((row.get("node_jaccard"),row.get("backbone_edge_jaccard")))
        key=(row["Disease"],row["parameter_to"])
        neighbor_values.setdefault(key,[]).append((row.get("node_jaccard"),row.get("backbone_edge_jaccard")))
    for index,row in metrics.iterrows():
        key=(row["Disease"],row["parameter_id"]); vals=neighbor_values.get(key,[])
        node=[float(x[0]) for x in vals if x[0] is not None and np.isfinite(float(x[0]))]
        edge=[float(x[1]) for x in vals if x[1] is not None and np.isfinite(float(x[1]))]
        metrics.loc[index,"neighbor_node_jaccard_mean"] = float(np.mean(node)) if node else np.nan
        metrics.loc[index,"neighbor_edge_jaccard_mean"] = float(np.mean(edge)) if edge else np.nan
    controls_path=os.path.join(RUN,"summaries","control_stability_by_disease_parameter.csv")
    controls=read_table(controls_path) if os.path.isfile(controls_path) else pd.DataFrame()
    if not controls.empty:
        for col in ["median_node_jaccard","q25_node_jaccard","median_edge_jaccard","q25_edge_jaccard"]:
            controls[col]=pd.to_numeric(controls[col],errors="coerce")
    control_wide=[]
    if not controls.empty:
        for (disease,pid),group in controls.groupby(["Disease","parameter_id"],sort=False):
            row={"Disease":disease,"parameter_id":pid}
            for _,item in group.iterrows():
                prefix=str(item["control_type"])
                for col in ["completed_repeats","median_node_jaccard","q25_node_jaccard","median_edge_jaccard","q25_edge_jaccard"]:
                    row[prefix+"_"+col]=item.get(col)
            control_wide.append(row)
    per_disease=metrics.merge(pd.DataFrame(control_wide),on=["Disease","parameter_id"],how="left") if control_wide else metrics.copy()
    write_rows(os.path.join(RUN,"summaries","shared_default_per_disease_metrics.csv"),per_disease.to_dict("records"))

    summary=[]
    expected_all_diseases=set(DISEASES)
    active_all_diseases=set(ACTIVE_DISEASES)
    base_status_path=os.path.join(RUN,"summaries","base_run_status.csv")
    base_status=read_table(base_status_path) if os.path.isfile(base_status_path) else pd.DataFrame()
    replicate_path=os.path.join(RUN,"summaries","replicate_network_metrics.csv")
    replicate_rows=read_table(replicate_path) if os.path.isfile(replicate_path) else pd.DataFrame()
    for pid,group in metrics.groupby("parameter_id",sort=False):
        row={"parameter_id":pid,"b":float(group["b"].iloc[0]),"w":int(group["w"].iloc[0]),"g":int(group["g"].iloc[0]),
            "valid_disease_count":int(group["Disease"].nunique()),
            "connected_network_disease_count":int((group["real_skeleton_edge_count"]>0).sum()),
            "empty_solution_diseases":"|".join(group.loc[group["real_node_count"]==0,"Disease"].astype(str)),
            "isolated_only_diseases":"|".join(group.loc[(group["real_node_count"]>0)&(group["real_skeleton_edge_count"]==0),"Disease"].astype(str)),
            "other_non_tree_or_multicomponent_diseases":"|".join(group.loc[(group["component_count"]>1)|(group["isolated_node_count"]>0),"Disease"].astype(str))}
        base_group=base_status[base_status["parameter_id"].astype(str)==str(pid)] if not base_status.empty else pd.DataFrame()
        base_ok=(active_all_diseases==expected_all_diseases and set(group["Disease"].astype(str))==expected_all_diseases and
                 not base_group.empty and set(base_group["Disease"].astype(str))==expected_all_diseases and
                 set(base_group["status"].astype(str)).issubset({"complete","complete_reused"}))
        controls_ok=bool(base_ok and not replicate_rows.empty)
        if controls_ok:
            for disease in DISEASES:
                for control_type in CTRL_TYPES:
                    subset=replicate_rows[(replicate_rows["Disease"].astype(str)==disease)&
                        (replicate_rows["parameter_id"].astype(str)==str(pid))&
                        (replicate_rows["control_type"].astype(str)==control_type)&
                        (replicate_rows["status"].astype(str)=="complete")]
                    reps=set(pd.to_numeric(subset["replicate"],errors="coerce").dropna().astype(int))
                    if reps!=set(range(1,REPEATS+1)):
                        controls_ok=False;break
                if not controls_ok: break
        row["all_16_frozen_disease_inputs_valid"]=active_all_diseases==expected_all_diseases
        row["base_grid_complete_for_all_16_diseases"]=bool(base_ok)
        row["all_three_controls_complete_20_repeats_per_disease"]=bool(controls_ok)
        row["common_default_evidence_complete"]=bool(base_ok and controls_ok)
        for col in ["real_node_count","model_evidence_node_count","zero_reward_no_evidence_connector_count",
                    "real_skeleton_edge_count","component_count","isolated_node_count","largest_component_fraction",
                    "retained_SHAP_mass_fraction_mapped_PPI","connector_fraction","selected_node_fraction_top1pct_degree",
                    "skeleton_edge_fraction_incident_top1pct_degree","neighbor_node_jaccard_mean","neighbor_edge_jaccard_mean"]:
            vals=pd.to_numeric(group[col],errors="coerce").dropna().values
            row[col+"_median_across_diseases"]=float(np.median(vals)) if len(vals) else None
            row[col+"_q25_across_diseases"]=float(np.quantile(vals,.25)) if len(vals) else None
        if not controls.empty:
            for control_type in CTRL_TYPES:
                cg=controls[(controls["parameter_id"]==pid)&(controls["control_type"]==control_type)]
                if cg.empty: continue
                for col in ["median_node_jaccard","q25_node_jaccard","median_edge_jaccard","q25_edge_jaccard"]:
                    vals=pd.to_numeric(cg[col],errors="coerce").dropna().values
                    row[control_type+"_"+col+"_across_diseases_median"]=float(np.median(vals)) if len(vals) else None
                    row[control_type+"_"+col+"_across_diseases_q25"]=float(np.quantile(vals,.25)) if len(vals) else None
                row[control_type+"_diseases_with_20_repeats"]=int((pd.to_numeric(cg["completed_repeats"],errors="coerce")==REPEATS).sum())
        summary.append(row)
    summary_df=pd.DataFrame(summary)
    # The parameter frontier is transparent and multi-metric; no weighted composite score is constructed.
    maximize=["edge_cost_q25_node","edge_cost_q25_edge","neighbor_node","neighbor_edge","retention"]
    minimize=["connector","hub_node","hub_edge"]
    for source,target in [("edge_cost_q25_node","edge_cost_q25_node"),("edge_cost_q25_edge","edge_cost_q25_edge"),
                          ("neighbor_node","neighbor_node"),("neighbor_edge","neighbor_edge"),
                          ("retention","retention"),("connector","connector"),
                          ("hub_node","hub_node"),("hub_edge","hub_edge")]:
        if source.startswith("edge_cost"):
            col="edge_cost_"+source.replace("edge_cost_q25_", "q25_")+"_jaccard_across_diseases_median"
        elif source.startswith("neighbor"):
            col=source.replace("neighbor_node","neighbor_node_jaccard_mean").replace("neighbor_edge","neighbor_edge_jaccard_mean")+"_median_across_diseases"
        elif source=="retention": col="retained_SHAP_mass_fraction_mapped_PPI_median_across_diseases"
        elif source=="connector": col="connector_fraction_median_across_diseases"
        elif source=="hub_node": col="selected_node_fraction_top1pct_degree_median_across_diseases"
        else: col="skeleton_edge_fraction_incident_top1pct_degree_median_across_diseases"
        if col not in summary_df.columns:
            summary_df[col]=np.nan
        target_col="frontier_"+target
        summary_df[target_col]=pd.to_numeric(summary_df[col],errors="coerce")
    direction={"frontier_edge_cost_q25_node":1,"frontier_edge_cost_q25_edge":1,
        "frontier_neighbor_node":1,"frontier_neighbor_edge":1,"frontier_retention":1,
        "frontier_connector":-1,"frontier_hub_node":-1,"frontier_hub_edge":-1}
    dominated=[]; comparison_complete=[]
    frontier_cols=list(direction.keys())
    for i, a in summary_df.iterrows():
        a_values=pd.to_numeric(a[frontier_cols],errors="coerce").values.astype(float)
        a_complete=bool(np.all(np.isfinite(a_values)) and summary_df.loc[i,"common_default_evidence_complete"])
        comparison_complete.append(a_complete)
        is_dom=False
        if a_complete:
            a_vals={c:direction[c]*float(a[c]) for c in frontier_cols}
            for j,brow in summary_df.iterrows():
                if i == j: continue
                if not bool(brow.get("common_default_evidence_complete",False)):
                    continue
                b_values=pd.to_numeric(brow[frontier_cols],errors="coerce").values.astype(float)
                if not np.all(np.isfinite(b_values)):
                    continue
                b_vals={c:direction[c]*float(brow[c]) for c in frontier_cols}
                if all(b_vals[c]>=a_vals[c] for c in frontier_cols) and any(b_vals[c]>a_vals[c] for c in frontier_cols):
                    is_dom=True;break
        dominated.append(not is_dom)
    summary_df["multi_metric_nondominated_frontier"]=dominated
    summary_df["pareto_comparison_complete_same_metric_set"]=comparison_complete
    # Keep the three candidates traceable with lexicographic tie-breaks, not a composite score.
    frontier=summary_df[summary_df["multi_metric_nondominated_frontier"] & summary_df["pareto_comparison_complete_same_metric_set"]].copy()
    candidate_rows=[]
    if not frontier.empty:
        # Offer up to three transparent tradeoff representatives, without a
        # composite score: strongest cost robustness, highest SHAP retention,
        # and lowest combined connector/hub-dependence indicators.
        bases=[("cost_robustness","frontier_edge_cost_q25_node",False),
               ("SHAP_retention","frontier_retention",False),
               ("lower_connector_and_hub_dependence","frontier_connector",True)]
        chosen=set()
        for basis,col,ascending in bases:
            available=frontier[pd.to_numeric(frontier[col],errors="coerce").notna()].copy()
            if basis=="lower_connector_and_hub_dependence":
                dep_cols=["frontier_connector","frontier_hub_node","frontier_hub_edge"]
                available=available[np.all(np.isfinite(available[dep_cols].apply(pd.to_numeric,errors="coerce").values.astype(float)),axis=1)]
                if available.empty: continue
                # Minimize each listed dependence measure jointly by Pareto
                # eligibility; use their rank tuple only to break exact ties.
                for metric in dep_cols:
                    available[metric+"_rank"]=available[metric].rank(method="min",ascending=True)
                available["_tie_rank_sum"]=available[[x+"_rank" for x in dep_cols]].sum(axis=1)
                selected=available.sort_values(["_tie_rank_sum","parameter_id"],ascending=[True,True]).iloc[0]
            else:
                selected=available.sort_values([col,"parameter_id"],ascending=[ascending,True]).iloc[0]
            pid0=str(selected["parameter_id"])
            if pid0 not in chosen:
                chosen.add(pid0)
                candidate_rows.append({**selected.to_dict(),"candidate_basis":basis,
                    "candidate_note":"Pareto tradeoff candidate; no automatic formal parameter selection"})
        # If one complete parameter dominates every other complete parameter on
        # the same metric set, record that fact; otherwise leave the decision open.
        complete=summary_df[summary_df["pareto_comparison_complete_same_metric_set"]]
        if len(complete)>1:
            for idx,row in complete.iterrows():
                others=complete[complete.index!=idx]
                av={c:direction[c]*float(row[c]) for c in frontier_cols}
                if all(all(direction[c]*float(other[c])<=av[c] for c in frontier_cols) and
                       any(direction[c]*float(other[c])<av[c] for c in frontier_cols)
                       for _,other in others.iterrows()):
                    for cand in candidate_rows:
                        if str(cand["parameter_id"])==str(row["parameter_id"]):
                            cand["candidate_note"]="唯一在统一完整指标集上支配其余参数；仍需用户确认后登记正式参数"
                            cand["uniquely_dominates_complete_candidates"]=True
                    break
    if not candidate_rows:
        candidate_rows=[{"status":"no_complete_parameter_evidence_set",
            "reason":"需要全部16病输入、完整基础网格及每病每参数三类扰动各20次，并具备同一完整指标集"}]
    write_rows(os.path.join(RUN,"summaries","common_parameter_comparison.csv"),summary_df.to_dict("records"))
    write_rows(os.path.join(RUN,"summaries","common_default_candidates_top3_review_order.csv"),candidate_rows)

    # Parameter-local frontiers and selected-versus-neighbor evidence, one disease at a time.
    disease_frontier=[]
    for disease,group in per_disease.groupby("Disease",sort=False):
        cols=["edge_cost_median_node_jaccard","edge_cost_median_edge_jaccard","retained_SHAP_mass_fraction_mapped_PPI",
              "connector_fraction","selected_node_fraction_top1pct_degree","skeleton_edge_fraction_incident_top1pct_degree"]
        local_metric_directions={cols[0]:1,cols[1]:1,cols[2]:1,cols[3]:-1,cols[4]:-1,cols[5]:-1}
        for idx,row in group.iterrows():
            dominated_row=False
            row_values=pd.to_numeric(row[cols],errors="coerce").values.astype(float)
            controls20=all(pd.notna(row.get(control_type+"_completed_repeats")) and
                           int(float(row.get(control_type+"_completed_repeats")))==REPEATS
                           for control_type in CTRL_TYPES)
            row_complete=bool(np.all(np.isfinite(row_values)) and controls20 and
                              active_all_diseases==expected_all_diseases)
            for other_idx,other in group.iterrows():
                if other["parameter_id"]==row["parameter_id"]: continue
                other_controls_complete=all(pd.notna(other.get(control_type+"_completed_repeats")) and
                    int(float(other.get(control_type+"_completed_repeats")))==REPEATS for control_type in CTRL_TYPES)
                other_values=pd.to_numeric(other[cols],errors="coerce").values.astype(float)
                if (not row_complete or not other_controls_complete or
                        active_all_diseases!=expected_all_diseases or not np.all(np.isfinite(other_values))):
                    continue
                a={c:local_metric_directions[c]*float(row[c]) for c in cols}
                bvals={c:local_metric_directions[c]*float(other[c]) for c in cols}
                if all(bvals[c]>=a[c] for c in cols) and any(bvals[c]>a[c] for c in cols): dominated_row=True;break
            disease_frontier.append({"Disease":disease,"parameter_id":row["parameter_id"],
                "per_disease_nondominated":bool(row_complete and not dominated_row),"pareto_comparison_complete_same_metric_set":row_complete,"real_node_count":row["real_node_count"],
                "real_skeleton_edge_count":row["real_skeleton_edge_count"],
                "edge_cost_median_node_jaccard":row.get("edge_cost_median_node_jaccard"),
                "edge_cost_median_edge_jaccard":row.get("edge_cost_median_edge_jaccard"),
                "retained_SHAP_mass_fraction_mapped_PPI":row.get("retained_SHAP_mass_fraction_mapped_PPI"),
                "connector_fraction":row.get("connector_fraction"),
                "hub_node_fraction":row.get("selected_node_fraction_top1pct_degree"),
                "neighbor_node_jaccard_mean":row.get("neighbor_node_jaccard_mean"),
                "neighbor_edge_jaccard_mean":row.get("neighbor_edge_jaccard_mean")})
    write_rows(os.path.join(RUN,"summaries","per_disease_parameter_frontier.csv"),disease_frontier)
    common_vs_local=[]
    top_pids=[str(row["parameter_id"]) for row in candidate_rows if row.get("parameter_id")]
    for common_pid in top_pids:
        for disease in ACTIVE_DISEASES:
            common_row=per_disease[(per_disease["Disease"]==disease)&(per_disease["parameter_id"]==common_pid)]
            local_rows=[x for x in disease_frontier if x["Disease"]==disease and
                        x["per_disease_nondominated"] and x["pareto_comparison_complete_same_metric_set"]]
            if common_row.empty: continue
            c=common_row.iloc[0]
            local=pd.DataFrame(local_rows)
            out={"common_parameter_id":common_pid,"Disease":disease,
                 "common_real_node_count":c.get("real_node_count"),"common_real_skeleton_edge_count":c.get("real_skeleton_edge_count"),
                 "common_edge_cost_node_jaccard_median":c.get("edge_cost_median_node_jaccard"),
                 "common_edge_cost_edge_jaccard_median":c.get("edge_cost_median_edge_jaccard"),
                 "common_neighbor_node_jaccard_mean":c.get("neighbor_node_jaccard_mean"),
                 "common_neighbor_edge_jaccard_mean":c.get("neighbor_edge_jaccard_mean"),
                 "common_retained_SHAP_mass_fraction":c.get("retained_SHAP_mass_fraction_mapped_PPI"),
                 "common_connector_fraction":c.get("connector_fraction"),
                 "common_hub_node_fraction":c.get("selected_node_fraction_top1pct_degree"),
                 "local_nondominated_parameter_count":len(local)}
            compare_fields=[("edge_cost_node_jaccard_median","edge_cost_median_node_jaccard","max"),
                ("edge_cost_edge_jaccard_median","edge_cost_median_edge_jaccard","max"),
                ("retained_SHAP_mass_fraction","retained_SHAP_mass_fraction_mapped_PPI","max"),
                ("connector_fraction","connector_fraction","min"),
                ("hub_node_fraction","hub_node_fraction","min"),
                ("neighbor_node_jaccard_mean","neighbor_node_jaccard_mean","max"),
                ("neighbor_edge_jaccard_mean","neighbor_edge_jaccard_mean","max")]
            local_key_map={"hub_node_fraction":"hub_node_fraction"}
            for output_name,field,direction_name in compare_fields:
                common_value=out.get("common_"+output_name)
                if field=="edge_cost_median_node_jaccard": source="edge_cost_median_node_jaccard"
                elif field=="edge_cost_median_edge_jaccard": source="edge_cost_median_edge_jaccard"
                elif field=="hub_node_fraction": source="hub_node_fraction"
                else: source=field
                values=pd.to_numeric(local[source],errors="coerce").dropna().values if (not local.empty and source in local.columns) else np.asarray([])
                reference=(float(np.max(values)) if direction_name=="max" else float(np.min(values))) if len(values) else None
                out["local_frontier_best_"+output_name]=reference
                out["common_minus_local_frontier_"+output_name]=(float(common_value)-reference if common_value is not None and reference is not None and pd.notna(common_value) else None)
            common_vs_local.append(out)
    write_rows(os.path.join(RUN,"summaries","common_candidates_vs_per_disease_frontier.csv"),common_vs_local)

    # Selection frequencies across the single unperturbed base grid are explicitly parameter frequencies.
    node_freq=[];edge_freq=[]
    for disease in ACTIVE_DISEASES:
        node_counts={};edge_counts={}
        for b,w,g in itertools.product(B_VALUES,W_VALUES,G_VALUES):
            pid=param_id(b,w,g); directory=os.path.join(RUN,"grid",disease,pid)
            npth=os.path.join(directory,"adopted_nodes.csv");epth=os.path.join(directory,"adopted_forest_edges.tsv")
            if not os.path.isfile(npth) or not os.path.isfile(epth): continue
            for eid in read_table(npth)["Entrez"].astype(str): node_counts[eid]=node_counts.get(eid,0)+1
            edf=read_table(epth,"\t")
            for a,bid in zip(edf.get("Entrez_A",[]),edf.get("Entrez_B",[])):
                ek="{}|{}".format(*sorted((str(a),str(bid))));edge_counts[ek]=edge_counts.get(ek,0)+1
        for eid,count in node_counts.items(): node_freq.append({"Disease":disease,"Entrez":eid,"selected_parameter_count":count,"grid_denominator":36,"parameter_frequency":count/36.0})
        for edge,count in edge_counts.items(): edge_freq.append({"Disease":disease,"edge_key":edge,"selected_parameter_count":count,"grid_denominator":36,"parameter_frequency":count/36.0})
    write_rows(os.path.join(RUN,"summaries","base_grid_node_selection_frequency.csv"),node_freq)
    write_rows(os.path.join(RUN,"summaries","base_grid_edge_selection_frequency.csv"),edge_freq)
    data_dictionary=[
        {"file":"parameter_network_metrics.csv","meaning":"Disease x parameter base network metrics and objective comparisons"},
        {"file":"solver_objective_audit.csv","meaning":"raw solver and legal same-objective witness objective terms"},
        {"file":"parameter_neighbor_similarity.csv","meaning":"within-disease Jaccard for one-axis adjacent grid settings"},
        {"file":"replicate_network_metrics.csv","meaning":"all completed edge-cost and reward-shuffle replicate metrics and selected ID references"},
        {"file":"control_node_selection_frequency.csv","meaning":"per Disease/parameter/control type node selection count divided by actual completed replicates"},
        {"file":"control_edge_selection_frequency.csv","meaning":"per Disease/parameter/control type backbone edge selection count divided by actual completed replicates"},
        {"file":"common_parameter_comparison.csv","meaning":"cross-disease unweighted medians/q25 and multi-metric nondominated flag; no weighted score"},
        {"file":"common_default_candidates_top3_review_order.csv","meaning":"transparent lexicographic review order within nondominated combinations, not automated approval"},
        {"file":"per_disease_parameter_frontier.csv","meaning":"per-disease sensitivity frontier, not disease-specific formal parameter selection"},
        {"file":"base_grid_node_selection_frequency.csv","meaning":"frequency across 36 parameter settings, not replicate stability"},
        {"file":"base_grid_edge_selection_frequency.csv","meaning":"frequency across 36 parameter settings, not replicate stability"}]
    write_rows(os.path.join(RUN,"figure_source","data_dictionary.csv"),data_dictionary)
    log("cross-disease summaries complete; no composite parameter score was calculated")


def _save_network(out_dir,disease,pid,best,raw,graph,ctx,ppi_df,metrics,candidates):
    node_rows=[]; edge_rows=[]; component_rows=[]
    n=len(graph.nodes); m=len(graph.edges)
    chosen=set(int(x) for x in best["vertices"] if int(x)<n)
    raw_chosen=set(int(x) for x in raw["stripped_vertices"])
    for ix in sorted(chosen):
        node_rows.append({"Disease":disease,"parameter_id":pid,"Entrez":ctx["nodes"][ix],"model_SHAP_evidence":bool(ctx["evidence"][ix]),
            "mean_abs_SHAP":ctx["raw"][ix],"normalized_prize":ctx["bare"][ix],"solver_prize":float(graph.prizes[ix]),
            "background_degree":int(ctx["degree"][ix]),"adopted_solution_source":best["solution_type"],"in_raw_pcst_fast":ix in raw_chosen})
    for ei in sorted(int(x) for x in best["edges"] if int(x)<m):
        u,v=map(int,graph.edges[ei]); edge_rows.append({"Disease":disease,"parameter_id":pid,"edge_index":ei,
            "Entrez_A":ctx["nodes"][u],"Entrez_B":ctx["nodes"][v],"effective_cost":float(graph.costs[ei]),"base_cost":1.0,"hub_penalty":float(graph.edge_penalties[ei])})
    component_adj = {ix: set() for ix in chosen}
    for ei in sorted(int(x) for x in best["edges"] if int(x) < m):
        u, v = map(int, graph.edges[ei]); component_adj[u].add(v); component_adj[v].add(u)
    components = []
    remaining = set(chosen)
    while remaining:
        root_ix = min(remaining)
        stack = [root_ix]; members = set()
        while stack:
            ix = stack.pop()
            if ix not in members:
                members.add(ix); stack.extend(component_adj[ix] - members)
        remaining -= members; components.append(members)
    components.sort(key=lambda members: min(members))
    for component_id, members in enumerate(components, 1):
        for ix in sorted(members):
            component_rows.append({"Disease": disease, "parameter_id": pid,
                "component_id": component_id, "Entrez": ctx["nodes"][ix],
                "model_SHAP_evidence": bool(ctx["evidence"][ix])})
    raw_nodes=[{"Disease":disease,"parameter_id":pid,"solver_vertex_index":int(ix),"Entrez":ctx["nodes"][int(ix)] if int(ix)<n else "VIRTUAL_ROOT","is_virtual_root":int(ix)==n} for ix in raw["raw_vertices"]]
    raw_edges=[]
    terminal=list(map(int,graph.terminals))
    for ei in raw["raw_edges"]:
        ei=int(ei)
        if ei<m:
            u,v=map(int,graph.edges[ei]); a=ctx["nodes"][u]; b=ctx["nodes"][v]; cost=float(graph.costs[ei]); role="real_PPI"
        else:
            pos=ei-m; node_ix=terminal[pos]; a="VIRTUAL_ROOT"; b=ctx["nodes"][node_ix]; cost=float(graph.params.w); role="virtual_root_connection"
        raw_edges.append({"Disease":disease,"parameter_id":pid,"solver_edge_index":ei,"Entrez_A":a,"Entrez_B":b,"cost":cost,"edge_role":role})
    adopted_root=[]
    for ei in sorted(int(x) for x in best["edges"] if int(x)>=m):
        pos=ei-m
        adopted_root.append({"Disease":disease,"parameter_id":pid,"virtual_root":"VIRTUAL_ROOT",
            "Entrez":ctx["nodes"][terminal[pos]],"virtual_root_cost":float(graph.params.w),"terminal_index":terminal[pos]})
    write_rows(os.path.join(out_dir,"adopted_nodes.csv"),node_rows,
        ["Disease","parameter_id","Entrez","model_SHAP_evidence","mean_abs_SHAP","normalized_prize","solver_prize","background_degree","adopted_solution_source","in_raw_pcst_fast"])
    write_rows(os.path.join(out_dir,"adopted_forest_edges.tsv"),edge_rows,
        ["Disease","parameter_id","edge_index","Entrez_A","Entrez_B","effective_cost","base_cost","hub_penalty"],delimiter="\t")
    write_rows(os.path.join(out_dir,"adopted_component_membership.csv"),component_rows,
        ["Disease","parameter_id","component_id","Entrez","model_SHAP_evidence"])
    write_rows(os.path.join(out_dir,"adopted_virtual_root_edges.tsv"),adopted_root,
        ["Disease","parameter_id","virtual_root","Entrez","virtual_root_cost","terminal_index"],delimiter="\t")
    write_rows(os.path.join(out_dir,"raw_solver_vertices.csv"),raw_nodes,
        ["Disease","parameter_id","solver_vertex_index","Entrez","is_virtual_root"])
    write_rows(os.path.join(out_dir,"raw_solver_edges.csv"),raw_edges,
        ["Disease","parameter_id","solver_edge_index","Entrez_A","Entrez_B","cost","edge_role"])
    write_rows(os.path.join(out_dir,"candidate_objective_comparison.csv"),candidates,
        ["solution_type","vertices","edges","feasible","feasibility_reason","objective","selected_reward","unselected_reward","selected_edge_cost","raw_vertices","raw_edges"])
    np.savez_compressed(os.path.join(out_dir,"solver_arrays.npz"), raw_vertices=raw["raw_vertices"], raw_edges=raw["raw_edges"],
        adopted_vertices=np.asarray(best["vertices"],dtype=np.int64), adopted_edges=np.asarray(best["edges"],dtype=np.int64),
        prizes=np.asarray(graph.prizes,dtype=np.float64), real_costs=np.asarray(graph.costs,dtype=np.float64),
        root_terminals=np.asarray(graph.terminals,dtype=np.int64), edge_uv=np.asarray(graph.edges,dtype=np.int64),
        solver_edges=np.asarray(raw["solver_edges"],dtype=np.int64),
        solver_prizes=np.asarray(raw["solver_prizes"],dtype=np.float64),
        solver_costs=np.asarray(raw["solver_costs"],dtype=np.float64),
        solver_root=np.asarray([raw["solver_root"]],dtype=np.int64),
        solver_num_clusters=np.asarray([raw["solver_num_clusters"]],dtype=np.int64))


def _output_hashes(out_dir):
    rows=[]
    for filename in ["adopted_nodes.csv","adopted_forest_edges.tsv","adopted_component_membership.csv","adopted_virtual_root_edges.tsv","raw_solver_vertices.csv","raw_solver_edges.csv","candidate_objective_comparison.csv","solver_arrays.npz"]:
        path=os.path.join(out_dir,filename)
        rows.append({"path":path,"sha256":sha256(path),"size_bytes":os.path.getsize(path)})
    return rows


def build_neighbor_similarity(metrics, ppi_df):
    by = {(r["Disease"],float(r["b"]),int(r["w"]),int(r["g"])):r for r in metrics}
    grid_edges=[]
    for disease in ACTIVE_DISEASES:
        for b,w,g in itertools.product(B_VALUES,W_VALUES,G_VALUES):
            here=by.get((disease,b,w,g))
            if here is None: continue
            for axis,vals in [("b",B_VALUES),("w",W_VALUES),("g",G_VALUES)]:
                ix=vals.index({"b":b,"w":w,"g":g}[axis])
                if ix+1>=len(vals): continue
                nxt={"b":b,"w":w,"g":g}; nxt[axis]=vals[ix+1]
                other=by.get((disease,nxt["b"],nxt["w"],nxt["g"]))
                if other is None: continue
                node_a=_read_ids(disease,b,w,g,"adopted_nodes.csv")
                node_b=_read_ids(disease,nxt["b"],nxt["w"],nxt["g"],"adopted_nodes.csv")
                edge_a=_read_edge_keys(disease,b,w,g)
                edge_b=_read_edge_keys(disease,nxt["b"],nxt["w"],nxt["g"])
                grid_edges.append({"Disease":disease,"parameter_from":param_id(b,w,g),"parameter_to":param_id(nxt["b"],nxt["w"],nxt["g"]),
                    "changed_parameter":axis,"node_jaccard":jaccard(node_a,node_b),"backbone_edge_jaccard":jaccard(edge_a,edge_b),
                    "from_nodes":len(node_a),"to_nodes":len(node_b),"from_edges":len(edge_a),"to_edges":len(edge_b)})
    return grid_edges


def _read_ids(disease,b,w,g,filename):
    path=os.path.join(RUN,"grid",disease,param_id(b,w,g),filename)
    if not os.path.isfile(path): return set()
    df=read_table(path)
    if "Entrez" not in df.columns: raise RuntimeError("Missing Entrez column in {}".format(path))
    return set(df["Entrez"].astype(str))


def _read_edge_keys(disease,b,w,g):
    path=os.path.join(RUN,"grid",disease,param_id(b,w,g),"adopted_forest_edges.tsv")
    if not os.path.isfile(path): return set()
    df=read_table(path,"\t")
    if not {"Entrez_A","Entrez_B"}.issubset(set(df.columns)):
        raise RuntimeError("Missing edge endpoint columns in {}".format(path))
    return {"{}|{}".format(*sorted((str(a),str(b)))) for a,b in zip(df["Entrez_A"],df["Entrez_B"])}


def audit_only(ppi_df,node_df,disease_data):
    degrees=pd.to_numeric(node_df["background_degree"],errors="raise").astype(float).values
    ppi_rows=[{"ppi_node_count":len(node_df),"ppi_edge_count":len(ppi_df),
        "ppi_component_count":int(node_df["component_id"].nunique()) if "component_id" in node_df.columns else None,
        "number_min":float(pd.to_numeric(ppi_df["number_numeric_audit"],errors="raise").min()),
        "number_median":float(pd.to_numeric(ppi_df["number_numeric_audit"],errors="raise").median()),
        "number_max":float(pd.to_numeric(ppi_df["number_numeric_audit"],errors="raise").max()),
        "degree_mean":float(np.mean(degrees)),"degree_median":float(np.median(degrees)),
        "degree_q90":float(np.quantile(degrees,.9)),"degree_max":float(np.max(degrees))}]
    write_rows(os.path.join(RUN,"inputs","ppi_background_summary.csv"),ppi_rows)
    rows=[]
    for disease in ACTIVE_DISEASES:
        d=disease_data[disease]; cfg=d["model_cfg"]; runtime=d["runtime"]
        rows.append({"Disease":disease,"training_seed":cfg["seed"],"positive_class":cfg["positive_class"],"negative_class":cfg["negative_class"],
            "other_diseases_enter_model":cfg["other_diseases_enter_model"],"model_features":cfg["n_model_genes"],
            "target_SHAP_train_cells":runtime.get("n_target_train_SHAP_cells"),"ctrl_SHAP_train_cells":runtime.get("n_Ctrl_train_SHAP_cells"),
            "target_SHAP_raw_total":float(np.nansum(d["raw"][d["evidence"]])),"model_mapped_nodes":int(d["evidence"].sum()),
            "positive_reward_nodes":int(np.sum(d["rewards"]>0)),"measured_zero_SHAP_nodes":int(np.sum(d["evidence"] & (d["raw"]==0))),
            "no_SHAP_evidence_PPI_connectors":int(np.sum(~d["evidence"])),"model_final_eval_AUROC":runtime.get("final_eval_AUROC"),
            "model_final_eval_AUPRC":runtime.get("final_eval_AUPRC"),"model_status":runtime.get("status")})
    write_rows(os.path.join(RUN,"summaries","input_identity_and_reward_audit.csv"),rows)
    log("input audit complete diseases={} ppi_nodes={} ppi_edges={}".format(len(rows),len(node_df),len(ppi_df)))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["audit", "base"], default="audit")
    args=parser.parse_args()
    initialize_layout()
    write_run_config()
    ppi_df,node_df,disease_data=load_inputs()
    if len(disease_data) != len(DISEASES) or set(disease_data) != set(DISEASES):
        raise RuntimeError("All 16 frozen disease inputs must pass identity checks before this scan")
    audit_only(ppi_df,node_df,disease_data)
    if args.phase == "base":
        run_base(ppi_df,node_df,disease_data)
    log("phase {} finished".format(args.phase))


if __name__=="__main__":
    try:
        main()
    except Exception:
        os.makedirs(os.path.join(RUN,"logs"),exist_ok=True)
        with open(os.path.join(RUN,"logs","fatal_error.log"),"a",encoding="utf-8") as handle:
            handle.write("[{}]\n{}\n".format(stamp(),traceback.format_exc()))
        raise

