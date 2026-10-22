#!/usr/bin/env python3
"""Exploratory pan14 PCSF grid and preliminary robustness workflow.

Uses only frozen PTG14 SHAP-to-PPI arrays and the installed OI2/pcst_fast
implementation. It never writes to formal networks/ or Fig5 panels.
"""
from __future__ import print_function

import argparse
import ast
import csv
import hashlib
import itertools
import json
import math
import os
import platform
import sys
import time
import traceback
from datetime import datetime

import numpy as np
import pandas as pd
import pkg_resources
from OmicsIntegrator import graph as oi_graph

ROOT = r"E:\micro_output\分析审阅中心\待审阅\Fig5_pan14_PPI_PCSF更新"
AUDIT = os.path.join(ROOT, "robustness", "20261019_pcsf_quality_01")
RUN_ID = "20261019_pcsf_parameter_adaptation_01"
RUN = os.path.join(ROOT, "robustness", RUN_ID)
SEEDS = list(range(42, 52))
EXPECTED = {"nodes": 12512, "edges": 83065, "components": 34,
            "evidence_nodes": 4175, "features": 4843, "root_edges": 4175}
B_VALUES = [2, 5, 10]
W_VALUES = [1, 3, 6]
G_VALUES = [2, 3, 5, 6]
REPEATS = 20
SIGMA = 0.1
PRUNING = "strong"
NUM_CLUSTERS = 1
BASE_COST = 1.0
CTRL_TYPES = ["edge_cost", "reward_shuffle", "degree_stratified_reward_shuffle"]


def stamp():
    return datetime.now().astimezone().isoformat(timespec="seconds")


def log(message):
    line = "[{}] {}".format(stamp(), message)
    print(line, flush=True)
    with open(os.path.join(RUN, "logs", "execution.log"), "a", encoding="utf-8") as f:
        f.write(line + "\n")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sha_array(arr):
    a = np.ascontiguousarray(arr)
    return hashlib.sha256(a.tobytes(order="C")).hexdigest()


def write_json(path, value):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(value, f, ensure_ascii=False, indent=2, sort_keys=True,
                  allow_nan=False)
        f.write("\n")
    os.replace(tmp, path)


def write_rows(path, rows, fields, delimiter=","):
    tmp = path + ".tmp"
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(tmp, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter=delimiter,
                                extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            cooked = {}
            for field in fields:
                value = row.get(field, "")
                if value is None or (isinstance(value, float) and not math.isfinite(value)):
                    value = ""
                cooked[field] = value
            writer.writerow(cooked)
    os.replace(tmp, path)


def pandas_read_path(path, **kwargs):
    # The installed pandas C parser cannot initialize some paths containing
    # non-ASCII characters on this Windows setup. Opening the file first keeps
    # the fast parser while preserving the same bytes and parsing options.
    with open(path, "rb") as handle:
        return pd.read_csv(handle, **kwargs)


def read_tsv(path):
    return pandas_read_path(path, sep="\t", dtype=str, keep_default_na=False)


def read_csv(path):
    return pandas_read_path(path, dtype=str, keep_default_na=False)


def param_id(b, w, g):
    return "b{}_w{}_g{}".format(b, w, g)


def ensure_layout(phase="grid"):
    # Never mix this execution with an unrelated or stale run directory. An
    # empty pre-created directory is fine; any existing files require an
    # explicit matching run identity and script hash before we touch outputs.
    if os.path.isdir(RUN):
        existing_files=[]
        for dirpath,_,filenames in os.walk(RUN):
            existing_files.extend(os.path.join(dirpath,name) for name in filenames)
        if existing_files:
            cfgpath=os.path.join(RUN,"config","run_config.json")
            if not os.path.isfile(cfgpath):
                raise RuntimeError("Existing run directory contains files but no run_config.json; refusing to modify it")
            with open(cfgpath,"r",encoding="utf-8") as f:
                cfg=json.load(f)
            current_script_hash=sha256(os.path.abspath(__file__))
            frozen_script_hash=cfg.get("execution_script_sha256")
            if cfg.get("run_id")!=RUN_ID:
                raise RuntimeError("Existing run directory identity mismatch; refusing to resume")
            if frozen_script_hash!=current_script_hash:
                patch_path=os.path.join(RUN,"config","compatibility_patch_registry.json")
                if not os.path.isfile(patch_path):
                    raise RuntimeError("Existing run directory script hash mismatch without a reviewed compatibility patch")
                patch=json.load(open(patch_path,encoding="utf-8"))
                if (patch.get("parent_script_sha256")!=frozen_script_hash or
                        patch.get("patched_script_sha256")!=current_script_hash or
                        patch.get("scope")!="control metric propagation and summary loop arity only; no input, objective, solver, or parameter changes" or
                        phase not in ("controls","summarize")):
                    raise RuntimeError("Compatibility patch registry does not authorize this phase/script/run")
    for folder in ["config", "inputs", "scripts", "runs", "perturbations",
                   "nulls", "summaries", "figure_source", "logs"]:
        os.makedirs(os.path.join(RUN, folder), exist_ok=True)
    if not os.path.isdir(RUN):
        raise RuntimeError("Run directory could not be created")


def copy_frozen_inputs():
    import shutil
    dst = os.path.join(RUN, "inputs", "frozen")
    os.makedirs(dst, exist_ok=True)
    files = [
        "inputs/oi_base_node_order.tsv", "inputs/oi_base_edge_array.tsv",
        "inputs/oi_input_nodes.csv", "inputs/oi_input_edges.tsv",
        "inputs/ppi_nodes_ge2.csv", "inputs/ppi_edges_ge2.tsv",
        "inputs/full_gene_mapping_audit.csv", "inputs/full_seed_gene_weights.csv",
        "inputs/seed_source_chain_verification.json",
        "inputs/source_input_sha256_manifest.csv", "config/solver_quality_audit_manifest.json",
        "solver_audit/solver_objective_audit.csv",
    ]
    accepted_manifest_path = os.path.join(AUDIT, "config", "artifact_sha256_manifest.csv")
    if not os.path.isfile(accepted_manifest_path):
        raise FileNotFoundError(accepted_manifest_path)
    accepted = {str(r["path"]).replace("\\", "/"): str(r["sha256"]).lower()
                for r in read_csv(accepted_manifest_path).to_dict("records")}
    rows = []
    for rel in files:
        source = os.path.join(AUDIT, *rel.split("/"))
        if not os.path.isfile(source):
            raise FileNotFoundError(source)
        source_hash = sha256(source)
        if accepted.get(rel.replace("\\", "/")) != source_hash:
            raise RuntimeError("Frozen input differs from accepted artifact manifest: {}".format(rel))
        target = os.path.join(dst, os.path.basename(source))
        if not os.path.isfile(target):
            shutil.copy2(source, target)
        rows.append({"source_path": source, "copied_path": target,
                     "source_sha256": source_hash, "copied_sha256": sha256(target),
                     "size_bytes": os.path.getsize(source)})
    # The accepted chain's manifest_sha256 identifies the per-seed pan14
    # network run_manifest.json, whereas complete.json is a model-training
    # artifact checked by source_input_sha256_manifest.csv. Preserve and verify
    # these two manifest types separately instead of equating their hashes.
    chain_path=os.path.join(AUDIT,"inputs","seed_source_chain_verification.json")
    chain=json.load(open(chain_path,"r",encoding="utf-8"))
    manifest_copy_dir=os.path.join(dst,"network_run_manifests")
    os.makedirs(manifest_copy_dir,exist_ok=True)
    for seed in SEEDS:
        source=os.path.join(ROOT,"networks","seed_{}".format(seed),"run_manifest.json")
        expected=str(chain.get("seeds",{}).get(str(seed),{}).get("manifest_sha256","")).lower()
        if not os.path.isfile(source) or not expected or sha256(source).lower()!=expected:
            raise RuntimeError("Seed {} pan14 network run_manifest does not match accepted source-chain hash".format(seed))
        target=os.path.join(manifest_copy_dir,"seed_{}_run_manifest.json".format(seed))
        if not os.path.isfile(target):
            shutil.copy2(source,target)
        if sha256(target).lower()!=expected:
            raise RuntimeError("Seed {} copied pan14 network run_manifest hash mismatch".format(seed))
        rows.append({"source_path":source,"copied_path":target,"source_sha256":sha256(source),
                     "copied_sha256":sha256(target),"size_bytes":os.path.getsize(source)})
    if any(r["source_sha256"] != r["copied_sha256"] for r in rows):
        raise RuntimeError("Frozen input copy hash mismatch")
    write_rows(os.path.join(RUN, "inputs", "input_manifest.csv"), rows,
               ["source_path", "copied_path", "source_sha256", "copied_sha256", "size_bytes"])
    verify_source_chain(dst)
    return {os.path.basename(r["source_path"]): r["source_sha256"] for r in rows}


def verify_source_chain(frozen_dir):
    chain = json.load(open(os.path.join(frozen_dir, "seed_source_chain_verification.json"), "r", encoding="utf-8"))
    quality_manifest = json.load(open(os.path.join(frozen_dir, "solver_quality_audit_manifest.json"), "r", encoding="utf-8"))
    if chain.get("status") != "verified" or set(map(int, chain.get("seeds", {}).keys())) != set(SEEDS):
        raise RuntimeError("Accepted PTG14 source chain is not verified for all ten seeds")
    source_rows = read_csv(os.path.join(frozen_dir, "source_input_sha256_manifest.csv")).to_dict("records")
    source_by_path = {os.path.normcase(os.path.normpath(r["path"])): r for r in source_rows}
    for seed in SEEDS:
        item = chain["seeds"][str(seed)]
        if item.get("status") != "verified" or int(item.get("input_feature_count", -1)) != EXPECTED["features"] or int(item.get("full_shap_rows", -1)) != EXPECTED["features"]:
            raise RuntimeError("Seed {} source-chain feature/SHAP verification failed".format(seed))
        verified = item.get("verified_input_hashes", {})
        expected_suffixes = {"weights.csv": r"\runs\seed_{}\weights.csv".format(seed),
                             "candidates.csv": r"\runs\seed_{}\candidates.csv".format(seed),
                             "model.pt": r"\runs\seed_{}\model.pt".format(seed),
                             "complete.json": r"\runs\seed_{}\complete.json".format(seed)}
        for key, suffix in expected_suffixes.items():
            row = next((r for r in source_rows if r["path"].lower().endswith(suffix.lower())), None)
            expected_hash = verified.get(key)
            if row is None or (key!="complete.json" and
                               str(expected_hash or "").lower() != str(row["sha256"]).lower()):
                raise RuntimeError("Seed {} {} hash differs from accepted source manifest".format(seed, key))
            if key=="complete.json":
                with open(row["path"],"r",encoding="utf-8-sig") as f:
                    completion=json.load(f)
                completion_outputs=completion.get("outputs",{})
                for output_key in ["model.pt","weights.csv","candidates.csv"]:
                    if str(completion_outputs.get(output_key,"")).lower()!=str(verified.get(output_key,"")).lower():
                        raise RuntimeError("Seed {} complete.json output hash disagrees with verified {} hash".format(seed,output_key))
        common_suffixes = {
            "mentor_PPI": "Network-1.txt",
            "official_map": "orgHs_official_symbol_entrez.tsv",
            "alias_map": "orgHs_alias_entrez.tsv",
            "model_contract_json": r"\data\contract.json",
            "model_genes_csv": r"\data\genes.csv",
            "model_cells_csv": r"\data\cells.csv",
        }
        for key, suffix in common_suffixes.items():
            row = next((r for r in source_rows if r["path"].lower().endswith(suffix.lower())), None)
            expected_hash = item.get("model_contract_json_sha256") if key == "model_contract_json" else verified.get(key)
            if row is None or str(expected_hash or "").lower() != str(row["sha256"]).lower():
                raise RuntimeError("Seed {} {} hash differs from accepted source manifest".format(seed, key))
        if quality_manifest.get("model_contract_id_by_seed", {}).get(str(seed)) != item.get("contract_id"):
            raise RuntimeError("Seed {} model contract ID differs from solver quality audit".format(seed))
        network_manifest_live=os.path.join(ROOT,"networks","seed_{}".format(seed),"run_manifest.json")
        network_manifest_copy=os.path.join(frozen_dir,"network_run_manifests","seed_{}_run_manifest.json".format(seed))
        expected_manifest_hash=str(item.get("manifest_sha256","")).lower()
        if (not expected_manifest_hash or not os.path.isfile(network_manifest_live) or
                sha256(network_manifest_live).lower()!=expected_manifest_hash or
                not os.path.isfile(network_manifest_copy) or sha256(network_manifest_copy).lower()!=expected_manifest_hash):
            raise RuntimeError("Seed {} pan14 network run_manifest provenance hash mismatch".format(seed))
    # Live source hashes are checked once against the previously accepted chain.
    live = []
    for row in source_rows:
        path = row["path"]
        if not os.path.isfile(path):
            raise FileNotFoundError("Accepted source file is missing: {}".format(path))
        actual = sha256(path)
        if actual.lower() != str(row["sha256"]).lower():
            raise RuntimeError("Accepted source file hash changed: {}".format(path))
        live.append({"path": path, "sha256": actual, "status": "verified"})
    write_rows(os.path.join(RUN, "inputs", "source_chain_live_hash_check.csv"), live,
               ["path", "sha256", "status"])


def load_frozen():
    inp = os.path.join(RUN, "inputs", "frozen")
    node_df = read_tsv(os.path.join(inp, "oi_base_node_order.tsv"))
    edge_df = read_tsv(os.path.join(inp, "oi_base_edge_array.tsv"))
    evidence_df = read_csv(os.path.join(inp, "full_gene_mapping_audit.csv"))
    weights_df = read_csv(os.path.join(inp, "full_seed_gene_weights.csv"))
    ppi_nodes_df = read_csv(os.path.join(inp, "ppi_nodes_ge2.csv"))
    if len(node_df) != EXPECTED["nodes"] or node_df["internal_id"].nunique() != len(node_df):
        raise RuntimeError("Frozen PPI node count/order differs from accepted audit")
    if len(edge_df) != EXPECTED["edges"] or edge_df["edge_index"].nunique() != len(edge_df):
        raise RuntimeError("Frozen PPI edge count/order differs from accepted audit")
    if not np.array_equal(pd.to_numeric(node_df["node_index"], errors="raise").astype(np.int64).values,
                          np.arange(len(node_df), dtype=np.int64)):
        raise RuntimeError("Frozen node indices are not contiguous in accepted order")
    if not np.array_equal(pd.to_numeric(edge_df["edge_index"], errors="raise").astype(np.int64).values,
                          np.arange(len(edge_df), dtype=np.int64)):
        raise RuntimeError("Frozen edge indices are not contiguous in accepted order")
    if node_df["Entrez"].astype(str).duplicated().any() or not np.array_equal(
            node_df["internal_id"].astype(str).values,
            np.asarray(["ENTREZ:" + str(x) for x in node_df["Entrez"]], dtype=object)):
        raise RuntimeError("Frozen node index/internal ID/Entrez mapping is inconsistent")
    if len(ppi_nodes_df) != EXPECTED["nodes"]:
        raise RuntimeError("Frozen PPI node metadata count mismatch")
    if len(evidence_df) != EXPECTED["features"] * len(SEEDS):
        raise RuntimeError("Full gene mapping audit row count mismatch")
    if len(weights_df) != EXPECTED["features"] * len(SEEDS):
        raise RuntimeError("Full seed SHAP table row count mismatch")
    if set(pd.to_numeric(weights_df["seed"], errors="raise").astype(int)) != set(SEEDS):
        raise RuntimeError("Frozen SHAP table seed coverage mismatch")
    nodes = node_df["Entrez"].astype(str).tolist()
    node_ix = {eid: i for i, eid in enumerate(nodes)}
    edge_uv = edge_df[["node_u_index", "node_v_index"]].astype(np.int32).values
    expected_u = node_df.iloc[edge_uv[:, 0]]["internal_id"].astype(str).values
    expected_v = node_df.iloc[edge_uv[:, 1]]["internal_id"].astype(str).values
    if not np.array_equal(expected_u, edge_df["internal_id_u"].astype(str).values) or not np.array_equal(
            expected_v, edge_df["internal_id_v"].astype(str).values):
        raise RuntimeError("Frozen edge array endpoints do not map to the accepted node order")
    oi_nodes = read_csv(os.path.join(inp, "oi_input_nodes.csv"))
    oi_edges = read_tsv(os.path.join(inp, "oi_input_edges.tsv"))
    if oi_nodes["internal_id"].astype(str).tolist() != node_df["internal_id"].astype(str).tolist():
        raise RuntimeError("OI input node order differs from the accepted array node order")
    expected_oi_edges = pd.DataFrame({"protein1": expected_u, "protein2": expected_v,
                                      "cost": pd.to_numeric(edge_df["base_cost"], errors="raise").astype(float)})
    if not oi_edges[["protein1", "protein2"]].astype(str).equals(expected_oi_edges[["protein1", "protein2"]].astype(str)) or not np.allclose(
            pd.to_numeric(oi_edges["cost"], errors="raise").astype(float).values,
            expected_oi_edges["cost"].astype(float).values, rtol=0, atol=1e-12):
        raise RuntimeError("OI text input edge order/endpoints/cost differs from frozen indexed arrays")
    base_cost = pd.to_numeric(edge_df["base_cost"], errors="raise").astype(float).values
    frozen_g5 = pd.to_numeric(edge_df["effective_cost"], errors="raise").astype(float).values
    if edge_uv.min() < 0 or edge_uv.max() >= len(nodes):
        raise RuntimeError("Frozen PPI edge index out of bounds")
    if not np.all(edge_uv[:, 0] < edge_uv[:, 1]):
        raise RuntimeError("Frozen PPI edges are not normalized undirected endpoints")
    if not np.allclose(base_cost, BASE_COST, rtol=0, atol=1e-12):
        raise RuntimeError("Frozen PPI base costs are not all 1")
    degree = np.bincount(edge_uv.ravel(), minlength=len(nodes)).astype(np.int64)
    adjacency = [[] for _ in nodes]
    for u, v in edge_uv:
        adjacency[int(u)].append(int(v)); adjacency[int(v)].append(int(u))
    remaining = set(range(len(nodes))); observed_components = 0
    while remaining:
        observed_components += 1
        stack = [min(remaining)]; seen = set()
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x); stack.extend(y for y in adjacency[x] if y not in seen)
        remaining -= seen
    if observed_components != EXPECTED["components"]:
        raise RuntimeError("Independent connected-component count differs from accepted PPI")
    n_nodes = len(nodes)
    denom = ((n_nodes - degree[edge_uv[:, 0]] - 1) *
             (n_nodes - degree[edge_uv[:, 1]] - 1) +
             degree[edge_uv[:, 0]] * degree[edge_uv[:, 1]])
    penalty_g5 = 10.0 ** 5 * degree[edge_uv[:, 0]] * degree[edge_uv[:, 1]] / denom
    if not np.allclose(base_cost + penalty_g5, frozen_g5, rtol=1e-11, atol=1e-12):
        raise RuntimeError("Frozen effective costs do not match OI g=5 hub penalty")
    comp_count = ppi_nodes_df["component_id"].nunique()
    if comp_count != EXPECTED["components"]:
        raise RuntimeError("Frozen PPI component count mismatch")
    ppi_degree = dict(zip(ppi_nodes_df["Entrez"].astype(str),
                          pd.to_numeric(ppi_nodes_df["background_degree"], errors="raise").astype(int)))
    ppi_component = dict(zip(ppi_nodes_df["Entrez"].astype(str), ppi_nodes_df["component_id"].astype(str)))
    evidence = {}
    prizes = {}
    raw_shap = {}
    symbols = {}
    for seed in SEEDS:
        wseed = weights_df[pd.to_numeric(weights_df["seed"], errors="raise").astype(int) == seed].copy()
        e = evidence_df[pd.to_numeric(evidence_df["seed"], errors="raise").astype(int) == seed]
        if len(wseed) != EXPECTED["features"] or wseed["gene_id"].astype(str).nunique() != EXPECTED["features"]:
            raise RuntimeError("Seed {} SHAP table does not contain 4,843 unique model genes".format(seed))
        if len(e) != EXPECTED["features"] or e["gene_id"].astype(str).nunique() != EXPECTED["features"]:
            raise RuntimeError("Seed {} mapping audit does not contain 4,843 unique model genes".format(seed))
        if set(wseed["gene_id"].astype(str)) != set(e["gene_id"].astype(str)):
            raise RuntimeError("Seed {} SHAP and mapping gene keys differ".format(seed))
        feature_order_hash = hashlib.sha256("\n".join(wseed["gene_id"].astype(str).tolist()).encode("utf-8")).hexdigest()
        chain = json.load(open(os.path.join(inp, "seed_source_chain_verification.json"), "r", encoding="utf-8"))
        if feature_order_hash != chain["seeds"][str(seed)].get("input_feature_order_sha256"):
            raise RuntimeError("Seed {} full SHAP gene order differs from accepted model feature contract".format(seed))
        wseed["mean_abs_SHAP_num"] = pd.to_numeric(wseed["mean_abs_SHAP"], errors="raise").astype(float)
        if not np.all(np.isfinite(wseed["mean_abs_SHAP_num"].values)) or np.any(wseed["mean_abs_SHAP_num"].values < 0):
            raise RuntimeError("Seed {} SHAP contains invalid or negative weights".format(seed))
        wmap = wseed.set_index(wseed["gene_id"].astype(str))
        emap = e.set_index(e["gene_id"].astype(str))
        for gene_id in wmap.index:
            if not np.isclose(float(wmap.loc[gene_id, "mean_abs_SHAP_num"]),
                              float(pd.to_numeric(emap.loc[gene_id, "mean_abs_SHAP"], errors="raise")),
                              rtol=1e-10, atol=1e-14):
                raise RuntimeError("Seed {} SHAP value mismatch for {}".format(seed, gene_id))
        full_max = float(wseed["mean_abs_SHAP_num"].max())
        if full_max <= 0:
            raise RuntimeError("Seed {} full-feature SHAP maximum is not positive".format(seed))
        expected_gene_prize = wmap["mean_abs_SHAP_num"].astype(float).values / full_max
        observed_gene_prize = pd.to_numeric(emap.loc[wmap.index, "prize_full_model_max_normalized"], errors="raise").astype(float).values
        if not np.all(np.isfinite(observed_gene_prize)) or np.any(observed_gene_prize < 0) or not np.allclose(
                expected_gene_prize, observed_gene_prize, rtol=1e-10, atol=1e-13):
            raise RuntimeError("Seed {} full-feature reward normalization differs from accepted source".format(seed))
        valid = e[e["chosen_Entrez"].astype(str).isin(nodes)].copy()
        valid["mean_abs_SHAP_num"] = pd.to_numeric(valid["mean_abs_SHAP"], errors="raise").astype(float)
        valid["prize_num"] = pd.to_numeric(valid["prize_full_model_max_normalized"], errors="raise").astype(float)
        # Collapsed node evidence retains maximum gene-level SHAP per unique Entrez.
        collapsed = valid.sort_values(["mean_abs_SHAP_num", "gene_id"], ascending=[False, True]).drop_duplicates("chosen_Entrez")
        evidence[seed] = set(collapsed["chosen_Entrez"].astype(str))
        if len(evidence[seed]) != EXPECTED["evidence_nodes"]:
            raise RuntimeError("Seed {} unique mapped evidence nodes differ: {}".format(seed, len(evidence[seed])))
        raw_shap[seed] = dict(zip(collapsed["chosen_Entrez"].astype(str), collapsed["mean_abs_SHAP_num"].astype(float)))
        symbols[seed] = {}
        for eid, group in valid.groupby("chosen_Entrez", sort=False):
            symbols[seed][str(eid)] = "|".join(sorted(set(group["gene_id"].astype(str))))
        p = np.zeros(len(nodes), dtype=np.float64)
        for row in collapsed.itertuples(index=False):
            eid = str(getattr(row, "chosen_Entrez"))
            p[node_ix[eid]] = float(getattr(row, "prize_num"))
        # Independent cross-check: audited node-prize input has the same order and normalized prizes.
        nd = read_tsv(os.path.join(AUDIT, "solver_audit", "seed_{}".format(seed), "node_prize_array.tsv"))
        if nd["Entrez"].astype(str).tolist() != nodes:
            raise RuntimeError("Seed {} reward node order differs from frozen PPI".format(seed))
        p_audit = pd.to_numeric(nd["prize_unscaled"], errors="raise").astype(float).values
        if not np.allclose(p, p_audit, rtol=1e-11, atol=1e-14):
            raise RuntimeError("Seed {} reconstructed normalized prize mismatch".format(seed))
        if not np.isclose(np.max(p), 1.0, rtol=1e-12, atol=1e-14):
            raise RuntimeError("Seed {} full-model max-normalized prize denominator mismatch".format(seed))
        if int(np.count_nonzero(p > 0)) != EXPECTED["root_edges"]:
            raise RuntimeError("Seed {} positive terminal/root-edge count mismatch".format(seed))
        prizes[seed] = p
    return {"nodes": nodes, "node_ix": node_ix, "edge_uv": edge_uv,
            "edge_ids": [tuple(map(str, edge_df.iloc[i][["internal_id_u", "internal_id_v"]].tolist())) for i in range(len(edge_df))],
            "base_cost": base_cost, "degree": degree, "ppi_degree": ppi_degree,
            "ppi_component": ppi_component,
            "components": comp_count, "evidence": evidence, "prizes": prizes,
            "raw_shap": raw_shap, "symbols": symbols, "weights": weights_df,
            "edge_df": edge_df}


def make_adjacency(edge_uv, n_nodes):
    adj = [[] for _ in range(n_nodes)]
    for ix, (u, v) in enumerate(edge_uv):
        adj[int(u)].append((int(v), ix))
        adj[int(v)].append((int(u), ix))
    return adj


def finite_summary(values):
    a=np.asarray(values,dtype=np.float64);a=a[np.isfinite(a)]
    if not len(a):return {"n":0,"sum":None,"mean":None,"sd":None,"min":None,"p25":None,"median":None,"p75":None,"p90":None,"p95":None,"p99":None,"max":None}
    q=np.quantile(a,[.25,.5,.75,.9,.95,.99])
    return {"n":int(len(a)),"sum":float(np.sum(a)),"mean":float(np.mean(a)),"sd":float(np.std(a,ddof=1)) if len(a)>1 else 0.0,
        "min":float(np.min(a)),"p25":float(q[0]),"median":float(q[1]),"p75":float(q[2]),
        "p90":float(q[3]),"p95":float(q[4]),"p99":float(q[5]),"max":float(np.max(a))}


def save_input_audits(data):
    mapping_rows=[];reward_rows=[];component_rows=[];degree_rows=[];cost_rows=[]
    for seed in SEEDS:
        raw_total=float(pd.to_numeric(data["weights"].loc[pd.to_numeric(data["weights"]["seed"]).astype(int)==seed,"mean_abs_SHAP"],errors="raise").sum())
        mapped_total=float(sum(data["raw_shap"][seed].values()))
        p=data["prizes"][seed]
        maprows=data["evidence"][seed]
        mapping_rows.append({"model_seed":seed,"model_features":EXPECTED["features"],
            "unique_model_evidence_nodes_in_PPI":len(maprows),"PPI_nodes":len(data["nodes"]),
            "mapped_unique_raw_SHAP_sum":mapped_total,"full_model_raw_SHAP_sum":raw_total,
            "unique_node_SHAP_mass_fraction_of_full_features":mapped_total/raw_total if raw_total else None,
            "positive_reward_nodes":int(np.count_nonzero(p>0)),"zero_reward_mapped_evidence_nodes":int(sum(1 for eid in maprows if data["prizes"][seed][data["node_ix"][eid]]==0)),
            "PPI_only_zero_reward_connectors":len(data["nodes"])-len(maprows)})
        s=finite_summary(p[ p>0 ])
        reward_rows.append({"model_seed":seed,"reward_scale":"full 4,843-feature max-normalized prize before b",
            "positive_reward_node_count":int(np.count_nonzero(p>0)),**s,
            "top_1pct_reward_share":float(np.sort(p)[-max(1,int(math.ceil(len(p)*.01))):].sum()/p.sum()),
            "top_5pct_reward_share":float(np.sort(p)[-max(1,int(math.ceil(len(p)*.05))):].sum()/p.sum()),
            "top_10pct_reward_share":float(np.sort(p)[-max(1,int(math.ceil(len(p)*.10))):].sum()/p.sum())})
        by_comp={}
        for i,eid in enumerate(data["nodes"]):
            c=data["ppi_component"][eid]
            item=by_comp.setdefault(c,{"nodes":0,"evidence":0,"reward":0.0})
            item["nodes"]+=1
            if eid in maprows:item["evidence"]+=1;item["reward"]+=float(p[i])
        for c,x in sorted(by_comp.items()):
            component_rows.append({"model_seed":seed,"component_id":c,"PPI_nodes":x["nodes"],
                "model_evidence_nodes":x["evidence"],"normalized_reward_sum":x["reward"],
                "reward_fraction":x["reward"]/float(p.sum()) if p.sum() else None})
        evidence_idx=np.asarray([data["node_ix"][eid] for eid in maprows],dtype=np.int64)
        top=evidence_idx[np.argsort(p[evidence_idx])[::-1][:max(1,int(math.ceil(len(evidence_idx)*.01)))]]
        de=np.asarray([data["degree"][i] for i in evidence_idx],dtype=float)
        pr=p[evidence_idx]
        rank_corr=pd.Series(pr).corr(pd.Series(de),method="spearman") if len(pr)>1 else np.nan
        degree_rows.append({"model_seed":seed,"model_evidence_nodes":len(evidence_idx),
            "mean_evidence_background_degree":float(np.mean(de)),"top_1pct_evidence_mean_degree":float(np.mean(data["degree"][top])),
            "reward_degree_spearman":float(rank_corr) if np.isfinite(rank_corr) else None,
            "top_20_gene_entrez_prize_degree":"|".join("{}:{:.8g}:{}".format(data["nodes"][int(i)],p[int(i)],data["degree"][int(i)]) for i in evidence_idx[np.argsort(pr)[::-1][:20]])})
    for g in G_VALUES:
        vals=effective_costs(data,g)
        cost_rows.append({"g":g,"base_cost":BASE_COST,"hub_penalty_formula":"10^g*d(u)*d(v)/((N-d(u)-1)*(N-d(v)-1)+d(u)*d(v))",
            **finite_summary(vals)})
    write_rows(os.path.join(RUN,"summaries","mapping_coverage_by_seed.csv"),mapping_rows,list(mapping_rows[0]))
    write_rows(os.path.join(RUN,"summaries","reward_distribution_by_seed.csv"),reward_rows,list(reward_rows[0]))
    write_rows(os.path.join(RUN,"summaries","reward_by_ppi_component.csv"),component_rows,list(component_rows[0]))
    write_rows(os.path.join(RUN,"summaries","reward_degree_profile.csv"),degree_rows,list(degree_rows[0]))
    write_rows(os.path.join(RUN,"summaries","edge_cost_distribution_by_g.csv"),cost_rows,list(cost_rows[0]))


def components_from_edges(selected_nodes, selected_edges, edge_uv):
    selected = set(int(x) for x in selected_nodes)
    adj = {v: [] for v in selected}
    for ei in selected_edges:
        u, v = map(int, edge_uv[int(ei)])
        if u not in selected or v not in selected:
            raise RuntimeError("Selected edge endpoint is absent from selected nodes")
        adj[u].append(v); adj[v].append(u)
    comps = []
    remaining = set(selected)
    while remaining:
        root = min(remaining)
        stack = [root]; seen = set()
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            stack.extend(y for y in adj[x] if y not in seen)
        remaining -= seen
        comps.append(sorted(seen))
    comps.sort(key=lambda c: (-len(c), c[0] if c else -1))
    return comps


def connected_rooted(root, nodes, edges, aug_uv):
    node_set = set(int(x) for x in nodes)
    if root not in node_set:
        return False, "root_missing"
    edge_list = [int(x) for x in edges]
    if len(edge_list) != len(set(edge_list)):
        return False, "duplicate_edges"
    adj = {v: [] for v in node_set}
    for ei in edge_list:
        if ei < 0 or ei >= len(aug_uv):
            return False, "edge_index_out_of_bounds"
        u, v = map(int, aug_uv[ei])
        if u not in node_set or v not in node_set:
            return False, "edge_endpoint_missing"
        adj[u].append(v); adj[v].append(u)
    seen = set(); stack = [root]
    while stack:
        v = stack.pop()
        if v not in seen:
            seen.add(v); stack.extend(adj[v])
    if seen != node_set:
        return False, "not_root_connected"
    if len(edge_list) != max(0, len(node_set) - 1):
        return False, "not_tree_edge_node_identity"
    return True, ""


def objective(nodes, edges, prizes, costs, root, n_real):
    chosen = set(int(x) for x in nodes)
    selected_cost = float(np.sum(costs[np.asarray(edges, dtype=np.int64)])) if len(edges) else 0.0
    total_prize = float(np.sum(prizes[:n_real]))
    selected_ids = np.asarray([i for i in chosen if i < n_real], dtype=np.int64)
    selected_prize = float(np.sum(prizes[selected_ids])) if len(selected_ids) else 0.0
    lost_prize = total_prize - selected_prize
    return {"objective": selected_cost + lost_prize, "selected_cost": selected_cost,
            "unselected_prize": lost_prize, "selected_prize": selected_prize,
            "net_gain_vs_empty": selected_prize - selected_cost}


def evaluate_solution(name, nodes, edges, real_prizes, real_costs, root_costs,
                      edge_uv, root, n_real, terminals, raw_preferred=False):
    n_real_edges = len(real_costs)
    terminal_list = list(map(int, terminals))
    aug_uv = []
    root_positions = []
    endpoint_error = ""
    for local_ix, ei in enumerate(map(int, edges)):
        if 0 <= ei < n_real_edges:
            aug_uv.append(edge_uv[ei])
        elif n_real_edges <= ei < n_real_edges + len(terminal_list):
            pos = ei - n_real_edges
            aug_uv.append((root, terminal_list[pos]))
            root_positions.append(pos)
        else:
            endpoint_error = "edge_index_out_of_bounds"
            break
    aug_uv = np.asarray(aug_uv, dtype=np.int64).reshape((-1, 2))
    node_list = list(map(int, nodes))
    legal, reason = connected_rooted(root, node_list, list(range(len(aug_uv))), aug_uv)
    if endpoint_error:
        legal, reason = False, endpoint_error
    if len(node_list) != len(set(node_list)):
        legal, reason = False, "duplicate_vertices"
    if legal:
        for v in node_list:
            if int(v) < 0 or int(v) > root:
                legal, reason = False, "vertex_index_out_of_bounds"; break
        selected_set = set(node_list)
        if legal and any(int(u) not in selected_set or int(v) not in selected_set for u, v in aug_uv):
            legal, reason = False, "unselected_edge_endpoint"
    terms = None
    if legal:
        selected_real_nodes = [i for i in node_list if i < n_real]
        selected_real_edges = [ei for ei in map(int, edges) if ei < n_real_edges]
        selected_cost = float(np.sum(real_costs[np.asarray(selected_real_edges, dtype=np.int64)])) if selected_real_edges else 0.0
        selected_cost += float(np.sum(root_costs[np.asarray(root_positions, dtype=np.int64)])) if root_positions else 0.0
        selected_prize = float(np.sum(real_prizes[np.asarray(selected_real_nodes, dtype=np.int64)])) if selected_real_nodes else 0.0
        total_prize = float(np.sum(real_prizes))
        terms = {"objective": selected_cost + total_prize - selected_prize,
                 "selected_cost": selected_cost, "selected_prize": selected_prize,
                 "unselected_prize": total_prize - selected_prize,
                 "net_gain_vs_empty": selected_prize - selected_cost}
    return {"candidate": name, "nodes": sorted(set(map(int, nodes))),
            "edges": sorted(set(map(int, edges))), "feasible": bool(legal),
            "reason": reason, "objective": terms["objective"] if terms else None,
            "selected_cost": terms["selected_cost"] if terms else None,
            "selected_prize": terms["selected_prize"] if terms else None,
            "unselected_prize": terms["unselected_prize"] if terms else None,
            "net_gain_vs_empty": terms["net_gain_vs_empty"] if terms else None,
            "raw_preferred": bool(raw_preferred)}


def candidate_solutions(raw_vertices, raw_edges, prizes, real_costs, root_costs,
                        edge_uv, terminals, w):
    n_real = len(prizes); n_real_edges = len(edge_uv); root = n_real
    aug_uv = np.vstack([edge_uv, np.column_stack([np.full(len(terminals), root, dtype=np.int32), terminals])])
    aug_costs = np.concatenate([real_costs, root_costs])
    all_prizes = np.concatenate([prizes, [0.0]])
    out = []
    # Raw arrays are exact return from pcst_fast before OI root stripping.
    out.append(evaluate_solution("raw_pcst_fast", raw_vertices, raw_edges, prizes,
                real_costs, root_costs, edge_uv, root, n_real, terminals, True))
    out.append(evaluate_solution("root_only_empty", [root], [], prizes,
                real_costs, root_costs, edge_uv, root, n_real, terminals))
    positive_terms = [int(i) for i in terminals]
    if positive_terms:
        best_single = min(positive_terms, key=lambda i: (objective([root, i], [n_real_edges + positive_terms.index(i)], all_prizes, aug_costs, root, n_real)["objective"], i))
        best_single_edge = n_real_edges + positive_terms.index(best_single)
        out.append(evaluate_solution("best_singleton", [root, best_single], [best_single_edge], prizes,
                    real_costs, root_costs, edge_uv, root, n_real, terminals))
    else:
        out.append(evaluate_solution("best_singleton", [root], [], prizes,
                    real_costs, root_costs, edge_uv, root, n_real, terminals))
    best_pair = None
    terminal_position = {int(node): i for i, node in enumerate(terminals)}
    total_prize = float(np.sum(prizes))
    for ei, (u, v) in enumerate(edge_uv):
        eligible = [int(x) for x in (u, v) if int(x) in terminal_position]
        for terminal in eligible:
            selected_prize = float(prizes[int(u)] + prizes[int(v)])
            value = float(total_prize - selected_prize + real_costs[int(ei)] + w)
            if best_pair is None or value < best_pair["objective"] - 1e-12:
                root_ei = n_real_edges + terminal_position[terminal]
                best_pair = {"candidate": "best_single_real_PPI_edge",
                    "nodes": [root, int(u), int(v)], "edges": [int(ei), root_ei],
                    "feasible": True, "reason": "", "objective": value,
                    "selected_cost": float(real_costs[int(ei)] + w),
                    "selected_prize": selected_prize,
                    "unselected_prize": float(total_prize - selected_prize),
                    "net_gain_vs_empty": float(selected_prize-real_costs[int(ei)]-w),
                    "raw_preferred": False}
    if best_pair:
        # Run the same independent rooted-tree legality check for the winning witness.
        best_pair = evaluate_solution("best_single_real_PPI_edge", best_pair["nodes"], best_pair["edges"],
            prizes, real_costs, root_costs, edge_uv, root, n_real, terminals)
        out.append(best_pair)
    else:
        out.append(evaluate_solution("best_single_real_PPI_edge", [root], [], prizes,
                    real_costs, root_costs, edge_uv, root, n_real, terminals))
    star_nodes = [root] + [int(i) for i in terminals if prizes[int(i)] > w]
    star_edges = [n_real_edges + terminal_position[int(i)] for i in terminals if prizes[int(i)] > w]
    out.append(evaluate_solution("root_star_prize_gt_w", star_nodes, star_edges, prizes,
                real_costs, root_costs, edge_uv, root, n_real, terminals))
    return out


def summarize_network(nodes, edge_ids, data, seed):
    edge_uv = data["edge_uv"]
    comps = components_from_edges(nodes, edge_ids, edge_uv)
    selected = set(int(x) for x in nodes)
    evidence_ids = data["evidence"][seed]
    evidence_selected = [i for i in selected if data["nodes"][i] in evidence_ids]
    connector_selected = [i for i in selected if data["nodes"][i] not in evidence_ids]
    singleton_count = sum(len(c) == 1 for c in comps)
    backbone_comps = [c for c in comps if len(c) >= 2]
    all_evidence_raw = float(sum(data["raw_shap"][seed].values()))
    retained = float(sum(data["raw_shap"][seed].get(data["nodes"][i], 0.0) for i in evidence_selected))
    evidence_fraction = retained / all_evidence_raw if all_evidence_raw > 0 else None
    node_fraction = len(evidence_selected) / len(evidence_ids) if evidence_ids else None
    enrich = evidence_fraction / node_fraction if evidence_fraction is not None and node_fraction else None
    deg_evidence = [int(data["degree"][i]) for i in evidence_selected]
    deg_conn = [int(data["degree"][i]) for i in connector_selected]
    ks = ks_statistic(deg_evidence, deg_conn)
    return {"real_node_count": len(selected), "real_skeleton_edge_count": len(edge_ids),
            "component_count_after_root_removal": len(comps), "isolated_node_count": singleton_count,
            "isolated_node_fraction": singleton_count / len(selected) if selected else None,
            "non_singleton_component_count": len(backbone_comps),
            "largest_component_nodes": max([len(c) for c in comps] or [0]),
            "largest_component_fraction": (max([len(c) for c in comps] or [0]) / len(selected)) if selected else None,
            "model_evidence_node_count": len(evidence_selected),
            "zero_evidence_PPI_connector_count": len(connector_selected),
            "connector_fraction": len(connector_selected) / len(selected) if selected else None,
            "retained_mapped_SHAP_mass_fraction": evidence_fraction,
            "SHAP_mass_enrichment_vs_node_fraction": enrich,
            "selected_evidence_mean_background_degree": float(np.mean(deg_evidence)) if deg_evidence else None,
            "selected_connector_mean_background_degree": float(np.mean(deg_conn)) if deg_conn else None,
            "evidence_vs_connector_degree_KS": ks,
            "all_mapped_evidence_raw_SHAP_sum": all_evidence_raw,
            "selected_evidence_raw_SHAP_sum": retained,
            "components": comps}


def ks_statistic(a, b):
    if not a or not b:
        return None
    x = np.sort(np.asarray(a)); y = np.sort(np.asarray(b))
    vals = np.unique(np.concatenate([x, y]))
    return float(max(abs(np.searchsorted(x, v, side="right") / float(len(x)) -
                         np.searchsorted(y, v, side="right") / float(len(y))) for v in vals))


def induced_edges(selected_nodes, adj):
    selected = set(int(x) for x in selected_nodes)
    found = set()
    for u in selected:
        for v, ei in adj[u]:
            if v in selected:
                found.add(int(ei))
    return sorted(found)


def solution_metrics(candidate, prizes, costs, root_costs, edge_uv, data, seed):
    if not candidate["feasible"]:
        return {"objective": None, "validity": "invalid"}
    real_nodes = [i for i in candidate["nodes"] if i < len(prizes)]
    real_edges = [i for i in candidate["edges"] if i < len(costs)]
    root_edges = [i - len(costs) for i in candidate["edges"] if i >= len(costs)]
    obj = objective(candidate["nodes"], candidate["edges"],
                    np.concatenate([prizes, [0.0]]),
                    np.concatenate([costs, root_costs]), len(prizes), len(prizes))
    net = summarize_network(real_nodes, real_edges, data, seed)
    return {"objective": obj["objective"], "selected_prize": obj["selected_prize"],
            "unselected_prize": obj["unselected_prize"],
            "real_edge_cost": float(np.sum(costs[np.asarray(real_edges, dtype=np.int64)])) if real_edges else 0.0,
            "root_edge_cost": float(np.sum(root_costs[np.asarray(root_edges, dtype=np.int64)])) if root_edges else 0.0,
            "component_count_after_root_removal": net["component_count_after_root_removal"],
            "real_node_count": net["real_node_count"], "real_edge_count": net["real_skeleton_edge_count"],
            "zero_evidence_PPI_connector_count": net["zero_evidence_PPI_connector_count"],
            "retained_mapped_SHAP_mass_fraction": net["retained_mapped_SHAP_mass_fraction"],
            "isolated_node_count": net["isolated_node_count"], "validity": "feasible"}


def deterministic_seed(key):
    digest = hashlib.sha256(key.encode("utf-8")).digest()
    return int.from_bytes(digest[:4], byteorder="little", signed=False)


def stable_key(param, seed, kind, rep):
    return "{}|{}|{}|{}".format(param, seed, kind, rep)


def degree_bins(evidence_indices, degrees, min_size=100):
    # Preserve equal-degree groups; merge neighboring degree groups greedily.
    by_degree = {}
    for i in evidence_indices:
        by_degree.setdefault(int(degrees[int(i)]), []).append(int(i))
    groups = [(d, sorted(v)) for d, v in sorted(by_degree.items())]
    bins = []; pending = []; dmin = None; dmax = None
    for d, members in groups:
        if dmin is None:
            dmin = d
        dmax = d; pending.extend(members)
        if len(pending) >= min_size:
            bins.append({"degree_min": dmin, "degree_max": dmax, "indices": pending})
            pending = []; dmin = None; dmax = None
    if pending:
        if bins:
            bins[-1]["degree_max"] = dmax
            bins[-1]["indices"].extend(pending)
        else:
            bins.append({"degree_min": dmin, "degree_max": dmax, "indices": pending})
    if any(len(b["indices"]) < min_size for b in bins) and len(evidence_indices) >= min_size:
        raise RuntimeError("Degree bin merge rule failed minimum bin size")
    return bins


def run_manifest_path(param, seed, kind=None, rep=None):
    if kind is None:
        return os.path.join(RUN, "runs", param, "seed_{}".format(seed), "run_manifest.json")
    base = "perturbations" if kind == "edge_cost" else "nulls"
    return os.path.join(RUN, base, param, "seed_{}".format(seed), kind,
                        "rep_{:02d}".format(rep), "run_manifest.json")


def is_valid_done(manifest_path, config_hash):
    if not os.path.isfile(manifest_path):
        return False
    try:
        with open(manifest_path, "r", encoding="utf-8") as f:
            m = json.load(f)
        if m.get("status") != "complete" or m.get("config_sha256") != config_hash:
            return False
        for path, expected_hash in m.get("outputs", {}).items():
            # best_solution.npz is a small derivative of the already persisted
            # selected tables and can be rebuilt without another solver call.
            if os.path.basename(path) == "best_solution.npz":
                continue
            if not os.path.isfile(path) or sha256(path) != expected_hash:
                return False
        return True
    except Exception:
        return False


def read_prior_g5_costs(data, seed):
    edge = read_tsv(os.path.join(AUDIT, "solver_audit", "seed_{}".format(seed), "pcst_fast_call.json"))
    return None


def run_pcst_call(graph, prizes, costs, w, rng=None, cost_noise=False):
    # Invoke the same pcst_fast function and rooted augmentation used by OI2 Graph.pcsf.
    all_nodes = list(range(len(graph.nodes)))
    terminals = np.asarray(graph.terminals, dtype=np.int64)
    root = len(graph.nodes)
    dummy_edges = np.column_stack([np.full(len(terminals), root, dtype=np.int64), terminals]).astype(np.int64)
    augmented_edges = np.concatenate([np.asarray(graph.edges, dtype=np.int64), dummy_edges], axis=0)
    root_costs = np.full(len(terminals), float(w), dtype=np.float64)
    augmented_costs = np.concatenate([np.asarray(costs, dtype=np.float64), root_costs])
    augmented_prizes = np.concatenate([np.asarray(prizes, dtype=np.float64), [0.0]])
    # pcst_fast is called via OI's imported module function; no solver substitution.
    result = oi_graph.pcst_fast(augmented_edges, augmented_prizes, augmented_costs,
                                root, NUM_CLUSTERS, PRUNING, 0)
    return result, augmented_edges, augmented_costs, augmented_prizes, root, terminals, root_costs


def run_one(data, seed, b, w, g, config_hash, adj, graph_cache):
    pid = param_id(b, w, g)
    out_dir = os.path.join(RUN, "runs", pid, "seed_{}".format(seed))
    manifest_path = os.path.join(out_dir, "run_manifest.json")
    if is_valid_done(manifest_path, config_hash):
        return json.load(open(manifest_path, encoding="utf-8"))
    os.makedirs(out_dir, exist_ok=True)
    graph = graph_cache[seed]
    params = {"b": b, "w": w, "g": g, "edge_noise": 0, "dummy_mode": "terminals",
              "seed": 42, "skip_checks": False}
    graph._reset_hyperparameters(params=params)
    prize_rows = pd.DataFrame({"name": ["ENTREZ:" + eid for eid in data["nodes"] if data["prizes"][seed][data["node_ix"][eid]] > 0],
                               "prize": [float(data["prizes"][seed][data["node_ix"][eid]]) for eid in data["nodes"] if data["prizes"][seed][data["node_ix"][eid]] > 0]})
    graph._prepare_prizes(prize_rows)
    prizes = np.asarray(graph.prizes, dtype=np.float64)
    costs = np.asarray(graph.costs, dtype=np.float64)
    edges = np.asarray(graph.edges, dtype=np.int64)
    if len(graph.nodes) != EXPECTED["nodes"] or len(edges) != EXPECTED["edges"]:
        raise RuntimeError("OI graph dimensions changed during grid run")
    # Original input edge and node ordering are frozen from the accepted OI arrays.
    if [str(x) for x in graph.nodes] != ["ENTREZ:" + eid for eid in data["nodes"]]:
        raise RuntimeError("OI node order differs from frozen audit")
    if not np.array_equal(edges, data["edge_uv"]):
        raise RuntimeError("OI edge order/endpoints differ from frozen audit")
    if not np.allclose(prizes, data["prizes"][seed] * b, rtol=1e-11, atol=1e-13):
        raise RuntimeError("OI scaled prize vector differs from frozen seed SHAP")
    denom = ((len(graph.nodes) - graph.node_degrees[edges[:, 0]] - 1) *
             (len(graph.nodes) - graph.node_degrees[edges[:, 1]] - 1) +
             graph.node_degrees[edges[:, 0]] * graph.node_degrees[edges[:, 1]])
    expected_penalty = (10.0 ** g) * (graph.node_degrees[edges[:, 0]] * graph.node_degrees[edges[:, 1]] / denom)
    if not np.allclose(graph.edge_penalties, expected_penalty, rtol=1e-11, atol=1e-12):
        raise RuntimeError("OI hub penalty does not match the audited formula")
    start = time.time()
    raw_capture = []
    original = oi_graph.pcst_fast
    def capture(*args, **kwargs):
        result = original(*args, **kwargs)
        raw_capture.append({"args": args, "kwargs": kwargs, "result": result})
        return result
    oi_graph.pcst_fast = capture
    try:
        vertices_oi, edges_oi = graph.pcsf(pruning=PRUNING, verbosity_level=0)
    finally:
        oi_graph.pcst_fast = original
    elapsed = time.time() - start
    if len(raw_capture) != 1:
        raise RuntimeError("Expected one OI -> pcst_fast call")
    raw_v = np.asarray(raw_capture[0]["result"][0], dtype=np.int64)
    raw_e = np.asarray(raw_capture[0]["result"][1], dtype=np.int64)
    root = len(data["nodes"])
    terminals = np.asarray(graph.terminals, dtype=np.int64)
    root_costs = np.full(len(terminals), float(w), dtype=np.float64)
    comp = candidate_solutions(raw_v, raw_e, prizes, costs, root_costs,
                               data["edge_uv"], terminals, w)
    feasible = [c for c in comp if c["feasible"]]
    min_obj = min(c["objective"] for c in feasible)
    best = next(c for c in feasible if c["candidate"] == "raw_pcst_fast" and abs(c["objective"]-min_obj) <= 1e-12) if any(
        c["candidate"] == "raw_pcst_fast" and abs(c["objective"]-min_obj) <= 1e-12 for c in feasible) else min(
        feasible, key=lambda c: (c["objective"], c["candidate"]))
    # Ensure raw return is preferred in exact ties, per frozen comparator rule.
    raw_obj = next((c["objective"] for c in feasible if c["candidate"] == "raw_pcst_fast"), None)
    if raw_obj is not None and abs(raw_obj - min_obj) <= 1e-12:
        best = next(c for c in feasible if c["candidate"] == "raw_pcst_fast")
    metrics = solution_metrics(best, prizes, costs, root_costs, data["edge_uv"], data, seed)
    real_nodes = [i for i in best["nodes"] if i < len(data["nodes"])]
    real_edges = [i for i in best["edges"] if i < len(costs)]
    root_edges = [i-len(costs) for i in best["edges"] if i >= len(costs)]
    induced = induced_edges(real_nodes, adj)
    comps = components_from_edges(real_nodes, real_edges, data["edge_uv"])
    # Save the exact raw augmented return even if it contains an invalid index.
    raw_vertex_rows=[]
    for value in raw_v:
        ix=int(value)
        label="VIRTUAL_ROOT" if ix==root else ("ENTREZ:"+data["nodes"][ix] if 0<=ix<len(data["nodes"]) else "OUT_OF_RANGE")
        raw_vertex_rows.append({"index":ix,"internal_id":label,"is_virtual_root":ix==root,
                                "index_valid":0<=ix<=root})
    write_rows(os.path.join(out_dir,"pcst_fast_raw_vertices.csv"),raw_vertex_rows,
        ["index","internal_id","is_virtual_root","index_valid"])
    raw_edge_rows=[]
    solver_aug_edges=np.asarray(raw_capture[0]["args"][0],dtype=np.int64)
    solver_aug_costs=np.asarray(raw_capture[0]["args"][2],dtype=np.float64)
    for value in raw_e:
        ei=int(value)
        if 0<=ei<len(solver_aug_edges):
            u,v=map(int,solver_aug_edges[ei]); c=float(solver_aug_costs[ei])
            etype="virtual_root" if ei>=len(costs) else "real_PPI"
            edge_valid=True
        else:
            u=v=None;c=None;etype="invalid_index";edge_valid=False
        raw_edge_rows.append({"index":ei,"edge_type":etype,"node_u_index":u,"node_v_index":v,
                              "cost":c,"index_valid":edge_valid})
    write_rows(os.path.join(out_dir,"pcst_fast_raw_edges.csv"),raw_edge_rows,
        ["index","edge_type","node_u_index","node_v_index","cost","index_valid"])
    write_rows(os.path.join(out_dir, "feasible_comparison.csv"), comp,
        ["candidate", "nodes", "edges", "feasible", "reason", "objective", "selected_cost",
         "selected_prize", "unselected_prize", "net_gain_vs_empty", "raw_preferred"])
    write_rows(os.path.join(out_dir, "selected_nodes.csv"), [
        {"node_index": i, "Entrez": data["nodes"][i],
         "gene_symbols": data["symbols"][seed].get(data["nodes"][i], ""),
         "model_evidence": data["nodes"][i] in data["evidence"][seed],
         "raw_mean_abs_SHAP": data["raw_shap"][seed].get(data["nodes"][i]),
         "normalized_prize": data["prizes"][seed][i], "solver_prize": prizes[i],
         "background_degree": int(data["degree"][i]),
         "role": "model_evidence" if data["nodes"][i] in data["evidence"][seed] else "zero_reward_PPI_connector"}
        for i in real_nodes], ["node_index", "Entrez", "gene_symbols", "model_evidence", "raw_mean_abs_SHAP",
                              "normalized_prize", "solver_prize", "background_degree", "role"])
    write_rows(os.path.join(out_dir, "selected_forest_edges.tsv"), [
        {"edge_index": i, "Entrez_A": data["nodes"][int(data["edge_uv"][i, 0])],
         "Entrez_B": data["nodes"][int(data["edge_uv"][i, 1])], "effective_cost": costs[i]}
        for i in real_edges], ["edge_index", "Entrez_A", "Entrez_B", "effective_cost"], delimiter="\t")
    write_rows(os.path.join(out_dir, "selected_virtual_root_edges.tsv"), [
        {"augmented_edge_index": len(costs) + int(i), "root_index": root,
         "terminal_node_index": int(terminals[int(i)]),
         "terminal_Entrez": data["nodes"][int(terminals[int(i)])], "cost": float(root_costs[int(i)])}
        for i in root_edges],
        ["augmented_edge_index", "root_index", "terminal_node_index", "terminal_Entrez", "cost"], delimiter="\t")
    write_rows(os.path.join(out_dir, "selected_induced_ppi_edges.tsv"), [
        {"edge_index": i, "Entrez_A": data["nodes"][int(data["edge_uv"][i, 0])],
         "Entrez_B": data["nodes"][int(data["edge_uv"][i, 1])], "effective_cost": costs[i],
         "selected_in_forest": i in set(real_edges)} for i in induced],
        ["edge_index", "Entrez_A", "Entrez_B", "effective_cost", "selected_in_forest"], delimiter="\t")
    membership = []
    summary = []
    for j, group in enumerate(comps, 1):
        component_id = "C{:04d}".format(j)
        group_edges = [i for i in real_edges if int(data["edge_uv"][i, 0]) in set(group) and int(data["edge_uv"][i, 1]) in set(group)]
        is_candidate = len(group) > 1 and bool(group_edges)
        for i in group:
            membership.append({"node_index": i, "Entrez": data["nodes"][i], "component_id": component_id,
                               "is_candidate_connected_subnetwork": is_candidate})
        ev = [i for i in group if data["nodes"][i] in data["evidence"][seed]]
        conn = [i for i in group if data["nodes"][i] not in data["evidence"][seed]]
        sum_mass = sum(data["raw_shap"][seed].get(data["nodes"][i], 0.0) for i in ev)
        total_mass = float(sum(data["raw_shap"][seed].values()))
        fraction = sum_mass / total_mass if total_mass else None
        node_fraction = len(ev) / float(len(data["evidence"][seed])) if data["evidence"][seed] else None
        summary.append({"component_id": component_id, "node_count": len(group), "skeleton_edge_count": len(group_edges),
                        "is_candidate_connected_subnetwork": is_candidate, "evidence_node_count": len(ev),
                        "connector_node_count": len(conn), "selected_evidence_SHAP_sum": sum_mass,
                        "fraction_of_all_mapped_SHAP_mass": fraction,
                        "SHAP_mass_enrichment_vs_evidence_node_fraction": fraction / node_fraction if fraction is not None and node_fraction else None,
                        "mean_evidence_background_degree": float(np.mean([data["degree"][i] for i in ev])) if ev else None,
                        "mean_connector_background_degree": float(np.mean([data["degree"][i] for i in conn])) if conn else None})
    write_rows(os.path.join(out_dir, "component_membership.csv"), membership,
               ["node_index", "Entrez", "component_id", "is_candidate_connected_subnetwork"])
    write_rows(os.path.join(out_dir, "component_summary.csv"), summary,
               ["component_id", "node_count", "skeleton_edge_count", "is_candidate_connected_subnetwork",
                "evidence_node_count", "connector_node_count", "selected_evidence_SHAP_sum",
                "fraction_of_all_mapped_SHAP_mass", "SHAP_mass_enrichment_vs_evidence_node_fraction",
                "mean_evidence_background_degree", "mean_connector_background_degree"])
    exact_arrays_path = os.path.join(out_dir, "solver_input_arrays.npz")
    np.savez_compressed(exact_arrays_path,
        real_edges=np.asarray(graph.edges, dtype=np.int64), real_prizes=prizes,
        real_costs=costs, augmented_edges=np.asarray(raw_capture[0]["args"][0], dtype=np.int64),
        augmented_prizes=np.asarray(raw_capture[0]["args"][1], dtype=np.float64),
        augmented_costs=np.asarray(raw_capture[0]["args"][2], dtype=np.float64),
        terminals=terminals, root_index=np.asarray([root], dtype=np.int64),
        root_edge_costs=root_costs)
    best_delta = (raw_obj - best["objective"]) if raw_obj is not None else None
    run_info = {"run_id": RUN_ID, "parameter_id": pid, "b": b, "w": w, "g": g,
        "model_seed": seed, "status": "complete", "raw_solver_result": "pcst_fast_raw_vertices.csv + pcst_fast_raw_edges.csv",
        "best_feasible_result_source": best["candidate"], "raw_objective": raw_obj,
        "best_feasible_objective": best["objective"], "improvement_over_raw": best_delta,
        "best_feasible_solution_improved_by_comparator": bool(best_delta is not None and best_delta > 1e-12),
        "solution_metrics": metrics, "runtime_seconds": elapsed,
        "solver": "OmicsIntegrator2 2.3.10 Graph.pcsf -> pcst_fast 1.0.7",
        "pruning": PRUNING, "num_clusters": NUM_CLUSTERS, "dummy_mode": "terminals",
        "root_index": root, "root_edge_count": len(terminals), "root_edge_cost": w,
        "array_sha256": {"real_edges": sha_array(edges), "scaled_prizes": sha_array(prizes),
                         "effective_costs": sha_array(costs), "augmented_edges": sha_array(raw_capture[0]["args"][0]),
                         "augmented_costs": sha_array(raw_capture[0]["args"][2])},
        "outputs": {}, "config_sha256": config_hash,
        "execution_script_sha256": sha256(os.path.abspath(__file__))}
    for filename in ["pcst_fast_raw_vertices.csv", "pcst_fast_raw_edges.csv", "feasible_comparison.csv",
                     "selected_nodes.csv", "selected_forest_edges.tsv", "selected_virtual_root_edges.tsv", "selected_induced_ppi_edges.tsv",
                     "component_membership.csv", "component_summary.csv"]:
        p = os.path.join(out_dir, filename); run_info["outputs"][p] = sha256(p)
    run_info["outputs"][exact_arrays_path] = sha256(exact_arrays_path)
    write_json(manifest_path, run_info)
    return run_info


def save_grid_config(data, input_hashes, script_hash):
    grid = [{"parameter_id": param_id(b,w,g), "b": b, "w": w, "g": g,
             "pruning": PRUNING, "dummy_mode": "terminals", "pcsf_seed": 42}
            for b,w,g in itertools.product(B_VALUES, W_VALUES, G_VALUES)]
    write_rows(os.path.join(RUN, "config", "parameter_grid.csv"), grid,
               ["parameter_id", "b", "w", "g", "pruning", "dummy_mode", "pcsf_seed"])
    cfg = {"run_id": RUN_ID, "created_local": stamp(), "status": "running",
        "objective": "J=sum(selected real+root edge costs)+sum(unselected real-node prizes)",
        "model": {"branch": "PTG14 pan14", "features": EXPECTED["features"], "labels": {"Neurological": 0, "Tumor": 1},
                  "Ctrl": "reference only", "model_seeds": SEEDS, "reward": "mean_abs_SHAP / max across all 4,843 features"},
        "ppi": {"source_audit": AUDIT, "filter": "number>=2", "nodes": EXPECTED["nodes"],
                "edges": EXPECTED["edges"], "components": EXPECTED["components"], "base_cost": BASE_COST,
                "all_components_retained": True},
        "solver": {"implementation": "OmicsIntegrator2 Graph.pcsf + installed pcst_fast", "OI_version": "2.3.10",
                   "pcst_fast_version": "1.0.7", "pruning": PRUNING, "dummy_mode": "terminals",
                   "num_clusters": NUM_CLUSTERS, "seed": 42, "edge_noise": 0,
                   "pcst_fast_low_level_last_argument": "verbosity_level=0; installed API has no per-call RNG seed argument",
                   "hub_penalty": "10^g * degree(u)*degree(v) / ((N-degree(u)-1)*(N-degree(v)-1)+degree(u)*degree(v))",
                   "root_prize": 0, "root_edge_cost": "w", "no_extra_component_penalty": True,
                   "OI_input_path_compatibility": "byte-identical ASCII-path cache for installed pandas 0.23 C-parser Unicode-path limitation; source remains RUN/inputs/frozen/oi_input_edges.tsv"},
        "grid": {"b": B_VALUES, "w": W_VALUES, "g": G_VALUES, "n_parameter_sets": len(grid), "n_base_calls": len(grid)*len(SEEDS)},
        "controls": {"repeats": REPEATS, "types": CTRL_TYPES, "edge_cost_sigma": SIGMA,
                     "edge_cost_noise_applied_to": "original effective real-edge cost (base cost + frozen OI hub penalty); virtual-root edges remain unchanged",
                     "eligible_parameter_rule": "at least one model seed with >=1 selected real PCSF backbone edge"},
        "feasible_comparator": ["raw pcst_fast", "root_only_empty", "best_singleton", "best_single_real_PPI_edge", "root_star_prize_gt_w"],
        "input_hashes": input_hashes, "data_dictionary_sha256": sha256(os.path.join(RUN,"config","data_dictionary.csv")),
        "execution_script_sha256": script_hash,
        "software": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__,
                     "pandas": pd.__version__, "omicsintegrator": pkg_resources.get_distribution("OmicsIntegrator").version,
                     "pcst_fast": pkg_resources.get_distribution("pcst-fast").version}}
    write_json(os.path.join(RUN, "config", "run_config.json"), cfg)
    return sha256(os.path.join(RUN, "config", "run_config.json"))


def make_graph(data, seed):
    # Hold only one OI graph in memory; the networkX-backed object is large.
    graph_path = os.path.join(RUN, "inputs", "frozen", "oi_input_edges.tsv")
    # OI2's pandas 0.23 C parser fails when given this frozen file through the
    # non-ASCII E: path. Use a byte-identical ASCII-path cache; this changes
    # only path handling, and the content hash is checked before every load.
    import shutil
    ascii_graph_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  "_pan14_pcsf_oi_input_edges.tsv")
    if not os.path.isfile(ascii_graph_path) or sha256(ascii_graph_path)!=sha256(graph_path):
        shutil.copy2(graph_path,ascii_graph_path)
    if sha256(ascii_graph_path)!=sha256(graph_path):
        raise RuntimeError("ASCII-path OI interactome copy differs from frozen graph input")
    graph = oi_graph.Graph(ascii_graph_path, {"w": 3, "b": 5, "g": 5, "edge_noise": 0,
                "dummy_mode": "terminals", "seed": 42, "skip_checks": False})
    if [str(x) for x in graph.nodes] != ["ENTREZ:" + e for e in data["nodes"]]:
        raise RuntimeError("Seed {} OI node order differs from accepted order".format(seed))
    if not np.array_equal(np.asarray(graph.edges, dtype=np.int64), data["edge_uv"]):
        raise RuntimeError("Seed {} OI edge order differs from accepted order".format(seed))
    return graph


def run_grid(data, config_hash):
    adj = make_adjacency(data["edge_uv"], len(data["nodes"]))
    grid = list(itertools.product(B_VALUES, W_VALUES, G_VALUES))
    status_path=os.path.join(RUN, "summaries", "base_run_status.csv")
    previous=read_csv(status_path).to_dict("records") if os.path.isfile(status_path) else []
    status_map={(r.get("parameter_id"),int(r.get("model_seed",-1))):r for r in previous}
    for seed in SEEDS:
        graph = make_graph(data, seed)
        for b,w,g in grid:
            pid = param_id(b,w,g)
            started = stamp()
            try:
                m = run_one(data, seed, b,w,g, config_hash, adj, {seed:graph})
                status_map[(pid,seed)]={"parameter_id": pid, "model_seed": seed, "status": "complete",
                                    "started_at": started, "finished_at": stamp(),
                                    "best_feasible_source": m.get("best_feasible_result_source", ""),
                                    "best_objective": m.get("best_feasible_objective", ""),
                                    "error": ""}
                log("base complete {} seed{} nodes={} edges={} source={} J={:.6f}".format(
                    pid, seed, m["solution_metrics"]["real_node_count"], m["solution_metrics"]["real_edge_count"],
                    m["best_feasible_result_source"], m["best_feasible_objective"]))
            except Exception as exc:
                status_map[(pid,seed)]={"parameter_id": pid, "model_seed": seed, "status": "failed",
                                    "started_at": started, "finished_at": stamp(), "error": repr(exc)}
                log("base FAILED {} seed{}: {}".format(pid, seed, traceback.format_exc()))
            write_rows(os.path.join(RUN, "summaries", "base_run_status.csv"), list(status_map.values()),
                       ["parameter_id", "model_seed", "status", "started_at", "finished_at",
                        "best_feasible_source", "best_objective", "error"])
        del graph
    return list(status_map.values())


def load_base_manifests(data):
    manifests = []
    for b,w,g in itertools.product(B_VALUES,W_VALUES,G_VALUES):
        pid = param_id(b,w,g)
        for seed in SEEDS:
            p = run_manifest_path(pid, seed)
            if os.path.isfile(p):
                with open(p, "r", encoding="utf-8") as f:
                    m = json.load(f)
                manifests.append(m)
    return manifests


def run_control(data, base_manifest, pid, seed, kind, rep, config_hash, adj, graph):
    if base_manifest is not None:
        b,w,g = int(base_manifest["b"]), int(base_manifest["w"]), int(base_manifest["g"])
    else:
        bits=pid.split("_")
        b,w,g=int(bits[0][1:]),int(bits[1][1:]),int(bits[2][1:])
    out_manifest = run_manifest_path(pid, seed, kind, rep)
    if is_valid_done(out_manifest, config_hash):
        with open(out_manifest, "r", encoding="utf-8") as f: return json.load(f)
    out_dir = os.path.dirname(out_manifest); os.makedirs(out_dir, exist_ok=True)
    key = stable_key(pid, seed, kind, rep)
    rng_seed = deterministic_seed(RUN_ID + "|" + key)
    rng = np.random.RandomState(rng_seed)
    graph._reset_hyperparameters(params={"b": b, "w": w, "g": g, "edge_noise": 0,
                  "dummy_mode": "terminals", "seed": 42, "skip_checks": False})
    # Original full-model prizes become scaled rewards exactly once.
    prizes = data["prizes"][seed].copy() * float(b)
    root = len(data["nodes"])
    terms = np.flatnonzero(prizes > 0).astype(np.int64)
    roots = np.full(len(terms), float(w), dtype=np.float64)
    base_real_costs = np.asarray(graph.costs, dtype=np.float64).copy()
    expected_costs = effective_costs(data, g)
    if not np.allclose(base_real_costs, expected_costs, rtol=1e-11, atol=1e-12):
        raise RuntimeError("Control graph effective costs differ from the frozen OI hub-penalty formula")
    if kind == "edge_cost":
        noise = rng.normal(0.0, SIGMA, size=len(base_real_costs))
        real_costs = base_real_costs * np.exp(noise - (SIGMA ** 2) / 2.0)
        prizes_for_call = prizes
        control_vector = np.exp(noise - (SIGMA ** 2) / 2.0).astype(np.float64)
    elif kind == "reward_shuffle":
        evidence_idx = np.asarray(sorted(data["node_ix"][eid] for eid in data["evidence"][seed]), dtype=np.int64)
        perm = rng.permutation(len(evidence_idx)).astype(np.int32)
        prizes_for_call = prizes.copy()
        prizes_for_call[evidence_idx] = prizes[evidence_idx[perm]]
        real_costs = base_real_costs
        control_vector = perm
    elif kind == "degree_stratified_reward_shuffle":
        evidence_idx = np.asarray(sorted(data["node_ix"][eid] for eid in data["evidence"][seed]), dtype=np.int64)
        bins = degree_bins(evidence_idx, data["degree"])
        prizes_for_call = prizes.copy()
        control_vector = np.empty(len(evidence_idx), dtype=np.int32)
        pos = {int(ix): j for j, ix in enumerate(evidence_idx)}
        bin_rows = []
        for bi, item in enumerate(bins, 1):
            idx = np.asarray(item["indices"], dtype=np.int64)
            perm = rng.permutation(len(idx)).astype(np.int32)
            prizes_for_call[idx] = prizes[idx[perm]]
            for source_local, target_local in enumerate(perm):
                control_vector[pos[int(idx[source_local])]] = int(idx[target_local])
            bin_rows.append({"bin_id": bi, "degree_min": item["degree_min"], "degree_max": item["degree_max"],
                             "n_evidence_nodes": len(idx), "members_sha256": sha_array(idx),
                             "member_node_indices": ";".join(map(str, idx.tolist())),
                             "member_Entrez": ";".join(data["nodes"][int(i)] for i in idx)})
        write_rows(os.path.join(out_dir, "degree_bins.csv"), bin_rows,
                   ["bin_id", "degree_min", "degree_max", "n_evidence_nodes", "members_sha256",
                    "member_node_indices", "member_Entrez"])
        real_costs = base_real_costs
    else:
        raise ValueError("Unknown control kind: {}".format(kind))
    # Root terminals follow the positive rewards after the permitted permutation.
    terms = np.flatnonzero(prizes_for_call > 0).astype(np.int64)
    roots = np.full(len(terms), float(w), dtype=np.float64)
    # Explicitly rebuild root augmentation exactly as Graph.pcsf; only costs or
    # evidence prizes are altered for the declared control.
    aug_edges = np.concatenate([np.asarray(graph.edges, dtype=np.int64),
             np.column_stack([np.full(len(terms), root, dtype=np.int64), terms])], axis=0)
    aug_costs = np.concatenate([real_costs, roots])
    aug_prizes = np.concatenate([prizes_for_call, [0.0]])
    start = time.time()
    raw_vertices, raw_edges = oi_graph.pcst_fast(aug_edges, aug_prizes, aug_costs, root,
                                                  NUM_CLUSTERS, PRUNING, 0)
    elapsed = time.time() - start
    raw_vertices = np.asarray(raw_vertices, dtype=np.int64)
    raw_edges = np.asarray(raw_edges, dtype=np.int64)
    comps = candidate_solutions(raw_vertices, raw_edges, prizes_for_call, real_costs,
                                roots, data["edge_uv"], terms, w)
    legal = [x for x in comps if x["feasible"]]
    min_obj = min(x["objective"] for x in legal)
    raw_c = next((x for x in legal if x["candidate"] == "raw_pcst_fast"), None)
    best = raw_c if raw_c is not None and abs(raw_c["objective"] - min_obj) <= 1e-12 else min(legal, key=lambda x:(x["objective"],x["candidate"]))
    mets = solution_metrics(best, prizes_for_call, real_costs, roots, data["edge_uv"], data, seed)
    real_nodes = np.asarray([x for x in best["nodes"] if x < len(data["nodes"])], dtype=np.int32)
    real_edges = np.asarray([x for x in best["edges"] if x < len(real_costs)], dtype=np.int32)
    root_edge_ids = np.asarray([x-len(real_costs) for x in best["edges"] if x >= len(real_costs)], dtype=np.int32)
    induced = np.asarray(induced_edges(real_nodes, adj), dtype=np.int32)
    # Compact binary preserves the exact random vector/permutation and selected indices.
    np.savez_compressed(os.path.join(out_dir, "replicate_arrays.npz"),
        control_vector=control_vector, raw_vertices=raw_vertices.astype(np.int32),
        raw_edges=raw_edges.astype(np.int32), selected_nodes=real_nodes,
        selected_real_edges=real_edges, selected_root_edge_positions=root_edge_ids,
        selected_induced_edges=induced, augmented_edges=aug_edges,
        augmented_prizes=aug_prizes, augmented_costs=aug_costs,
        terminals=terms.astype(np.int32), root_index=np.asarray([root],dtype=np.int32))
    base_solution_path=os.path.join(RUN,"runs",pid,"seed_{}".format(seed),"best_solution.npz")
    if base_manifest is not None and os.path.isfile(base_solution_path):
        base_npz = np.load(base_solution_path)
        base_nodes = set(int(x) for x in base_npz["nodes"])
        base_edges = set(int(x) for x in base_npz["edges"])
        node_j = jaccard(base_nodes, set(real_nodes.tolist()))
        edge_j = jaccard(base_edges, set(real_edges.tolist()))
        base_comparison_status="available"
    else:
        node_j = None; edge_j = None; base_comparison_status="base_run_unavailable"
    result = {"run_id": RUN_ID, "parameter_id": pid, "model_seed": seed, "control_type": kind,
        "replicate": rep, "random_key": key, "rng_algorithm": "numpy.random.RandomState(MT19937)",
        "rng_seed": rng_seed, "status": "complete", "termination": "pcst_fast returned",
        "runtime_seconds": elapsed, "best_feasible_source": best["candidate"],
        "objective": best["objective"], "raw_objective": raw_c["objective"] if raw_c else None,
        "objective_improvement_over_raw": raw_c["objective"]-best["objective"] if raw_c else None,
        "node_jaccard_vs_base": node_j, "skeleton_edge_jaccard_vs_base": edge_j,
        "base_comparison_status": base_comparison_status,
        "real_node_count": mets["real_node_count"], "real_edge_count": mets["real_edge_count"],
        "component_count": mets["component_count_after_root_removal"],
        "isolated_node_count": mets["isolated_node_count"],
        "connector_count": mets["zero_evidence_PPI_connector_count"],
        "retained_SHAP_fraction": mets["retained_mapped_SHAP_mass_fraction"],
        "control_vector_sha256": sha_array(control_vector), "node_order_sha256": sha_array(np.arange(len(data["nodes"]),dtype=np.int32)),
        "evidence_index_order_sha256": sha_array(np.asarray(sorted(data["node_ix"][eid] for eid in data["evidence"][seed]),dtype=np.int32)),
        "edge_order_sha256": sha_array(data["edge_uv"]),
        "augmented_edges_sha256": sha_array(aug_edges), "augmented_prizes_sha256": sha_array(aug_prizes),
        "augmented_costs_sha256": sha_array(aug_costs), "virtual_root_index": root,
        "virtual_root_edges": len(terms), "virtual_root_edge_cost": w,
        "outputs": {}, "config_sha256": config_hash,
        "execution_script_sha256": sha256(os.path.abspath(__file__))}
    write_rows(os.path.join(out_dir, "feasible_comparison.csv"), comps,
        ["candidate", "nodes", "edges", "feasible", "reason", "objective", "selected_cost",
         "selected_prize", "unselected_prize", "net_gain_vs_empty", "raw_preferred"])
    output_names = ["replicate_arrays.npz", "feasible_comparison.csv"]
    if kind == "degree_stratified_reward_shuffle":
        output_names.append("degree_bins.csv")
    for name in output_names:
        p = os.path.join(out_dir, name)
        if not os.path.isfile(p):
            raise RuntimeError("Required control output is missing: {}".format(p))
        result["outputs"][p] = sha256(p)
    write_json(out_manifest, result)
    return result


def jaccard(a, b):
    a = set(a); b = set(b)
    if not a and not b:
        return None
    return len(a & b) / float(len(a | b)) if a | b else None


def summarize_grid(data):
    run_rows=[]; metric_rows=[]; node_selections={}; edge_selections={}; component_members=[];component_summaries=[];quality_rows=[]
    base_status_path=os.path.join(RUN,"summaries","base_run_status.csv")
    previous_status=read_csv(base_status_path).to_dict("records") if os.path.isfile(base_status_path) else []
    failed_lookup={(r.get("parameter_id"),int(r.get("model_seed",-1))):r for r in previous_status if r.get("status")=="failed"}
    for b,w,g in itertools.product(B_VALUES,W_VALUES,G_VALUES):
        pid=param_id(b,w,g)
        for seed in SEEDS:
            mp=run_manifest_path(pid,seed)
            if not os.path.isfile(mp):
                fail=failed_lookup.get((pid,seed),{})
                run_rows.append({"parameter_id":pid,"model_seed":seed,"status":fail.get("status","missing"),"error":fail.get("error","")}); continue
            with open(mp,"r",encoding="utf-8") as f:m=json.load(f)
            if not is_valid_done(mp, m.get("config_sha256")):
                raise RuntimeError("Completed base run has a missing or hash-mismatched required artifact: {}".format(mp))
            nodes=pandas_read_path(os.path.join(RUN,"runs",pid,"seed_{}".format(seed),"selected_nodes.csv"),dtype={"Entrez":str})
            edges=read_tsv(os.path.join(RUN,"runs",pid,"seed_{}".format(seed),"selected_forest_edges.tsv"))
            node_ixs=set(pd.to_numeric(nodes.get("node_index",pd.Series(dtype=int)),errors="coerce").dropna().astype(int))
            edge_ixs=set(pd.to_numeric(edges.get("edge_index",pd.Series(dtype=int)),errors="coerce").dropna().astype(int))
            # Reconstruct metrics using frozen prizes and selected forest arrays.
            metrics=summarize_network(sorted(node_ixs),sorted(edge_ixs),data,seed)
            run_dir=os.path.join(RUN,"runs",pid,"seed_{}".format(seed))
            exact=np.load(os.path.join(run_dir,"solver_input_arrays.npz"))
            best_sol=np.load(os.path.join(run_dir,"best_solution.npz"))
            obj_terms=objective(best_sol["rooted_nodes"],best_sol["rooted_edges"],
                exact["augmented_prizes"],exact["augmented_costs"],len(data["nodes"]),len(data["nodes"]))
            if not np.isclose(obj_terms["objective"], float(m["best_feasible_objective"]), rtol=1e-10, atol=1e-8):
                raise RuntimeError("Seed {} {} objective does not reproduce from saved augmented arrays".format(seed,pid))
            best=m["best_feasible_result_source"]
            quality_rows.append({"parameter_id":pid,"model_seed":seed,"b":b,"w":w,"g":g,
                "raw_objective":m.get("raw_objective"),"best_feasible_objective":m.get("best_feasible_objective"),
                "best_feasible_source":best,"improved_raw":m.get("best_feasible_solution_improved_by_comparator"),
                "objective_delta_raw_minus_best":m.get("improvement_over_raw"),"runtime_seconds":m.get("runtime_seconds"),
                "raw_validity":next((x.get("feasible") for x in read_rows_csv(os.path.join(run_dir,"feasible_comparison.csv")) if x.get("candidate")=="raw_pcst_fast"),None),
                "comparison_rule":"minimum J among raw and four legal witnesses; raw wins ties"})
            for c in read_rows_csv(os.path.join(run_dir,"component_membership.csv")):
                component_members.append({"parameter_id":pid,"model_seed":seed,**c})
            for c in read_rows_csv(os.path.join(run_dir,"component_summary.csv")):
                component_summaries.append({"parameter_id":pid,"model_seed":seed,**c})
            row={"parameter_id":pid,"model_seed":seed,"b":b,"w":w,"g":g,"status":"complete",
                 "raw_objective":m.get("raw_objective"),"best_feasible_objective":m.get("best_feasible_objective"),
                 "best_solution_source":best,"improved_raw":m.get("best_feasible_solution_improved_by_comparator"),
                 "n_real_nodes":metrics["real_node_count"],"n_real_skeleton_edges":metrics["real_skeleton_edge_count"],
                 "n_induced_ppi_edges":len(pandas_read_path(os.path.join(RUN,"runs",pid,"seed_{}".format(seed),"selected_induced_ppi_edges.tsv"),sep="\t")),
                 "objective_recomputed_from_exact_augmented_arrays":obj_terms["objective"],
                 "retained_mapped_SHAP_mass_fraction":metrics["retained_mapped_SHAP_mass_fraction"],
                 "SHAP_mass_enrichment":metrics["SHAP_mass_enrichment_vs_node_fraction"],
                 "runtime_seconds":m.get("runtime_seconds"),"n_components":metrics["component_count_after_root_removal"],
                 "n_non_singleton_components":metrics["non_singleton_component_count"],
                 "n_isolated_nodes":metrics["isolated_node_count"],"isolated_fraction":metrics["isolated_node_fraction"],
                 "largest_component_nodes":metrics["largest_component_nodes"],
                 "largest_component_fraction":metrics["largest_component_fraction"],
                 "n_model_evidence_nodes":metrics["model_evidence_node_count"],
                 "n_zero_reward_PPI_connectors":metrics["zero_evidence_PPI_connector_count"],
                 "connector_fraction":metrics["connector_fraction"],
                 "selected_evidence_mean_degree":metrics["selected_evidence_mean_background_degree"],
                 "connector_mean_degree":metrics["selected_connector_mean_background_degree"],
                 "degree_KS_evidence_vs_connector":metrics["evidence_vs_connector_degree_KS"]}
            run_rows.append(row); metric_rows.append(row)
            key=(pid,seed)
            node_selections[key]=node_ixs; edge_selections[key]=edge_ixs
    fields=sorted(set(k for r in run_rows for k in r))
    write_rows(os.path.join(RUN,"summaries","base_run_summary_status.csv"),run_rows,fields)
    write_rows(os.path.join(RUN,"summaries","network_metrics_by_run.csv"),metric_rows,fields)
    write_rows(os.path.join(RUN,"summaries","solver_quality_by_run.csv"),quality_rows,
        ["parameter_id","model_seed","b","w","g","raw_objective","best_feasible_objective","best_feasible_source",
         "improved_raw","objective_delta_raw_minus_best","runtime_seconds","raw_validity","comparison_rule"])
    write_rows(os.path.join(RUN,"summaries","component_membership.csv"),component_members,
        ["parameter_id","model_seed","node_index","Entrez","component_id","is_candidate_connected_subnetwork"])
    write_rows(os.path.join(RUN,"summaries","component_summary.csv"),component_summaries,
        ["parameter_id","model_seed","component_id","node_count","skeleton_edge_count","is_candidate_connected_subnetwork",
         "evidence_node_count","connector_node_count","selected_evidence_SHAP_sum","fraction_of_all_mapped_SHAP_mass",
         "SHAP_mass_enrichment_vs_evidence_node_fraction","mean_evidence_background_degree","mean_connector_background_degree"])
    return node_selections,edge_selections,run_rows


def read_rows_csv(path):
    return read_csv(path).to_dict("records")


def effective_costs(data,g):
    uv=data["edge_uv"];deg=data["degree"];N=len(data["nodes"])
    d=((N-deg[uv[:,0]]-1)*(N-deg[uv[:,1]]-1)+deg[uv[:,0]]*deg[uv[:,1]])
    return data["base_cost"]+(10.0**g)*(deg[uv[:,0]]*deg[uv[:,1]]/d)


def write_best_npz_for_base():
    for b,w,g in itertools.product(B_VALUES,W_VALUES,G_VALUES):
        pid=param_id(b,w,g)
        for seed in SEEDS:
            d=os.path.join(RUN,"runs",pid,"seed_{}".format(seed));mfile=os.path.join(d,"run_manifest.json")
            if not os.path.isfile(mfile):continue
            m=json.load(open(mfile,encoding="utf-8"))
            best_path = os.path.join(d, "best_solution.npz")
            recorded_hash = m.get("outputs", {}).get(best_path)
            if os.path.isfile(best_path) and recorded_hash and sha256(best_path) == recorded_hash:
                continue
            best_nodes=pandas_read_path(os.path.join(d,"selected_nodes.csv"),dtype={"Entrez":str})["node_index"].astype(int).values
            best_edges=read_tsv(os.path.join(d,"selected_forest_edges.tsv"))["edge_index"].astype(int).values
            best_roots=read_tsv(os.path.join(d,"selected_virtual_root_edges.tsv"))["augmented_edge_index"].astype(int).values
            comp=read_rows_csv(os.path.join(d,"feasible_comparison.csv"))
            choice=next(x for x in comp if x["candidate"]==m["best_feasible_result_source"])
            rooted_nodes=np.asarray(ast.literal_eval(choice["nodes"]),dtype=np.int32)
            rooted_edges=np.asarray(ast.literal_eval(choice["edges"]),dtype=np.int32)
            np.savez_compressed(best_path,nodes=best_nodes.astype(np.int32),
                edges=best_edges.astype(np.int32),root_edge_indices=best_roots.astype(np.int32),
                rooted_nodes=rooted_nodes,rooted_edges=rooted_edges)
            m.setdefault("outputs",{})[best_path]=sha256(best_path)
            write_json(mfile,m)


def run_controls(data, config_hash):
    # User-frozen exploratory rule: evaluate a parameter setting only when at
    # least one completed base seed produced a real PCSF backbone edge. Once a
    # setting is eligible, run controls for all ten model seeds; a missing base
    # run only makes its base-Jaccards unavailable, not its controls ineligible.
    adj=make_adjacency(data["edge_uv"],len(data["nodes"]))
    base=load_base_manifests(data)
    base_by_key={}
    eligible=set()
    for m in base:
        pid=m.get("parameter_id");seed=int(m.get("model_seed",-1))
        if not pid or seed not in SEEDS:
            continue
        mpath=run_manifest_path(pid,seed)
        if m.get("status")=="complete" and is_valid_done(mpath,config_hash):
            base_by_key[(pid,seed)]=m
            if int(m.get("solution_metrics",{}).get("real_edge_count",0))>0:
                eligible.add(pid)
    all_params=set(param_id(b,w,g) for b,w,g in itertools.product(B_VALUES,W_VALUES,G_VALUES))
    write_json(os.path.join(RUN,"config","control_eligible_parameters.json"),
        {"rule":"Run controls for a parameter setting iff at least one valid completed base seed has >=1 real PCSF backbone edge; if eligible, run all ten model seeds. This is the user-frozen exploratory scope rule, not a suitability threshold.",
         "eligible_parameter_ids":sorted(eligible),"ineligible_parameter_ids":sorted(all_params-eligible),
         "eligible_control_seed_count_per_parameter":len(SEEDS),"control_types":CTRL_TYPES,"repeats_per_type_seed":REPEATS})
    status_path=os.path.join(RUN,"summaries","control_run_status.csv")
    previous=read_csv(status_path).to_dict("records") if os.path.isfile(status_path) else []
    status_map={(r.get("parameter_id"),int(r.get("model_seed",-1)),r.get("control_type"),int(r.get("replicate",-1))):r for r in previous}
    total=len(eligible)*len(SEEDS)*len(CTRL_TYPES)*REPEATS
    done=0
    for seed in SEEDS:
        graph=make_graph(data,seed)
        for pid in sorted(eligible):
            base_manifest=base_by_key.get((pid,seed))
            for kind in CTRL_TYPES:
                for rep in range(1,REPEATS+1):
                    started=stamp()
                    try:
                        r=run_control(data,base_manifest,pid,seed,kind,rep,config_hash,adj,graph)
                        status_map[(pid,seed,kind,rep)]={"parameter_id":pid,"model_seed":seed,"control_type":kind,
                            "replicate":rep,"status":"complete","started_at":started,"finished_at":stamp(),
                            "runtime_seconds":r["runtime_seconds"],"error":""}
                    except Exception:
                        error=traceback.format_exc()
                        status_map[(pid,seed,kind,rep)]={"parameter_id":pid,"model_seed":seed,"control_type":kind,
                            "replicate":rep,"status":"failed","started_at":started,"finished_at":stamp(),"error":error}
                        log("control FAILED {}/{}/{}/{}: {}".format(pid,seed,kind,rep,error))
                    done+=1
                    if done%25==0 or done==total:
                        write_rows(status_path,list(status_map.values()),
                            ["parameter_id","model_seed","control_type","replicate","status","started_at","finished_at","runtime_seconds","error"])
                        log("controls checkpoint {}/{}".format(done,total))
        del graph
    write_rows(status_path,list(status_map.values()),
        ["parameter_id","model_seed","control_type","replicate","status","started_at","finished_at","runtime_seconds","error"])
def write_data_dictionary():
    rows=[
      {"table":"input_manifest.csv","field":"source_path / copied_path","definition":"Frozen accepted-audit input path and byte-identical run-local copy."},
      {"table":"mapping_coverage_by_seed.csv","field":"unique_model_evidence_nodes_in_PPI","definition":"Unique Entrez-mapped model genes represented in the frozen PPI node set."},
      {"table":"reward_distribution_by_seed.csv","field":"positive_reward_node_count / quantiles","definition":"Full-feature max-normalized prizes before multiplication by b; distribution among positive mapped model nodes."},
      {"table":"reward_by_ppi_component.csv","field":"reward_fraction","definition":"Seed-specific normalized reward mass in each of the 34 frozen PPI components divided by total mapped reward mass."},
      {"table":"base_run_status.csv","field":"status","definition":"One row per parameter_id x model_seed; complete/failed attempt status and error text."},
      {"table":"solver_quality_by_run.csv","field":"raw_objective / best_feasible_objective","definition":"Objective comparison of raw PCST output and the lowest-J legal raw-or-witness solution; not a global optimality claim."},
      {"table":"network_metrics_by_run.csv","field":"n_real_nodes / n_real_skeleton_edges / components","definition":"Structure and model-evidence metrics for the selected legal rooted forest after removing virtual root edges."},
      {"table":"component_membership.csv","field":"component_id","definition":"Connected components of the real selected PCSF backbone; only components with a real edge qualify as connected candidate subnetworks."},
      {"table":"seed_pairwise_similarity.csv","field":"node_jaccard / evidence_node_jaccard / skeleton_edge_jaccard","definition":"Pairwise comparisons within a parameter setting; empty-empty is NA, failed runs are not treated as empty."},
      {"table":"parameter_neighbor_similarity.csv","field":"node_jaccard / evidence_node_jaccard / skeleton_edge_jaccard","definition":"Same-seed comparisons across adjacent listed values of exactly one parameter."},
      {"table":"node_selection_frequency.csv","field":"selected_count_observed / selection_frequency_of_10","definition":"Sparse node table: nodes never selected are omitted and are zero only when all ten seed runs completed; frequencies are NA for incomplete ten-seed settings."},
      {"table":"node_selection_frequency_by_seed.csv","field":"selected","definition":"Per-seed status for every node selected in at least one completed run; null on failed/missing runs."},
      {"table":"edge_selection_frequency.csv","field":"selected_count_observed / frequency_of_10","definition":"Sparse edge table; frequencies are only reported when all ten seed runs completed."},
      {"table":"control_coverage_by_parameter.csv","field":"control_status / missing_replicates","definition":"Coverage of the user-frozen control scope; noneligible settings are explicitly not applicable, not zero-result controls."},
      {"table":"perturbation_summary.csv","field":"mean_node_jaccard_vs_base / network metrics","definition":"Exploratory summaries across successfully recorded runs; not p-values, FDR, or independent model replicates."},
      {"table":"null_replicate_denominators.csv","field":"successful_replicates / failed_replicates / missing_replicates","definition":"Per setting, model seed, and reward-null type; omitted nodes in sparse selection summaries have count zero among successful repeats."},
      {"table":"cost_perturbation_node_selection_frequency.csv","field":"selected_count / replicate_denominator / selection_frequency","definition":"Sparse per-node inclusion frequency under edge-cost perturbations, using only hash-verified completed repeats; omitted nodes have count zero among successful repeats."},
      {"table":"parameter_recommendation.csv","field":"pareto_scope_status / pareto_nondominated","definition":"Multi-criteria, non-composite comparison; only complete ten-seed and complete pre-specified control settings are Pareto-comparable."},
      {"table":"run_status.csv","field":"run_kind / status / error","definition":"Unified base and control attempt ledger; missing attempts are counted separately in phase and coverage summaries."},
      {"table":"run_manifest.json","field":"outputs / config_sha256 / array_sha256","definition":"Per-run provenance, objective, selected solution, exact input-array hashes, software and output hashes."},
      {"table":"replicate_arrays.npz","field":"control_vector / augmented_* / selected_*","definition":"Exact per-control perturbation or permutation, solver arrays, raw return, selected nodes/edges, and root-edge positions."}
    ]
    write_rows(os.path.join(RUN,"config","data_dictionary.csv"),rows,["table","field","definition"])


def inspect_phase_status(phase):
    details={}
    cfgpath=os.path.join(RUN,"config","run_config.json")
    cfg=json.load(open(cfgpath,encoding="utf-8"))
    config_hash=hashlib.sha256((sha256(cfgpath)+str(cfg["execution_script_sha256"])).encode("ascii")).hexdigest()
    base_path=os.path.join(RUN,"summaries","base_run_status.csv")
    base_rows=read_csv(base_path).to_dict("records") if os.path.isfile(base_path) else []
    base_map={(r.get("parameter_id"),int(r.get("model_seed",-1))):r.get("status") for r in base_rows}
    expected_base={(param_id(b,w,g),seed) for b,w,g in itertools.product(B_VALUES,W_VALUES,G_VALUES) for seed in SEEDS}
    base_complete=base_failed=base_invalid=0
    for pid,seed in expected_base:
        status=base_map.get((pid,seed))
        if status=="complete":
            if is_valid_done(run_manifest_path(pid,seed),config_hash):base_complete+=1
            else:base_invalid+=1
        elif status=="failed":base_failed+=1
    base_missing=len(expected_base)-base_complete-base_failed-base_invalid
    details["base"]={"expected":len(expected_base),"complete":base_complete,"failed":base_failed,
                      "invalid_artifacts":base_invalid,"missing":base_missing,
                      "status":"complete" if base_complete==len(expected_base) else "partial"}
    eligible_path=os.path.join(RUN,"config","control_eligible_parameters.json")
    eligible=json.load(open(eligible_path,encoding="utf-8")).get("eligible_parameter_ids",[]) if os.path.isfile(eligible_path) else []
    control_expected={(pid,seed,kind,rep) for pid in eligible for seed in SEEDS for kind in CTRL_TYPES for rep in range(1,REPEATS+1)}
    ctrl_path=os.path.join(RUN,"summaries","control_run_status.csv")
    ctrl_rows=read_csv(ctrl_path).to_dict("records") if os.path.isfile(ctrl_path) else []
    ctrl_map={(r.get("parameter_id"),int(r.get("model_seed",-1)),r.get("control_type"),int(r.get("replicate",-1))):r.get("status") for r in ctrl_rows}
    ctrl_complete=ctrl_failed=ctrl_invalid=0
    for pid,seed,kind,rep in control_expected:
        status=ctrl_map.get((pid,seed,kind,rep))
        if status=="complete":
            if is_valid_done(run_manifest_path(pid,seed,kind,rep),config_hash):ctrl_complete+=1
            else:ctrl_invalid+=1
        elif status=="failed":ctrl_failed+=1
    ctrl_missing=len(control_expected)-ctrl_complete-ctrl_failed-ctrl_invalid
    control_status=("not_started" if not os.path.isfile(eligible_path) or not os.path.isfile(ctrl_path) else
                    ("complete" if ctrl_complete==len(control_expected) else "partial"))
    details["controls"]={"eligible_parameters":len(eligible),"expected":len(control_expected),"complete":ctrl_complete,
                         "failed":ctrl_failed,"invalid_artifacts":ctrl_invalid,"missing":ctrl_missing,
                         "status":control_status}
    required=[]
    # Every phase depends on a complete, hash-verified base grid. A request to
    # inspect the grid itself must therefore validate it as well.
    required.append(details["base"]["status"]=="complete")
    if phase in ("controls","all","summarize"):required.append(details["controls"]["status"]=="complete")
    if phase in ("summarize","all"):
        summary_names=["parameter_grid_summary.csv","network_metrics_by_run.csv","solver_quality_by_run.csv",
                       "node_selection_frequency.csv","seed_pairwise_similarity.csv","parameter_neighbor_similarity.csv",
                       "perturbation_summary.csv","null_selection_summary.csv","null_replicate_denominators.csv",
                       "cost_perturbation_node_selection_frequency.csv","control_coverage_by_parameter.csv",
                       "parameter_recommendation.csv"]
        missing_summary=[name for name in summary_names if not os.path.isfile(os.path.join(RUN,"summaries",name)) or os.path.getsize(os.path.join(RUN,"summaries",name))==0]
        details["summaries"]={"required":summary_names,"missing_or_empty":missing_summary,
                              "status":"complete" if not missing_summary else "partial"}
        required.append(not missing_summary)
    phase_ok=all(required) if required else False
    return ("complete" if phase_ok else "partial"),details

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--phase",choices=["grid","controls","summarize","all"],default="grid")
    args=ap.parse_args()
    ensure_layout(args.phase)
    input_hashes=copy_frozen_inputs()
    data=load_frozen()
    save_input_audits(data)
    write_data_dictionary()
    audit_names = ["mapping_coverage_by_seed.csv", "reward_distribution_by_seed.csv",
                   "reward_by_ppi_component.csv", "reward_degree_profile.csv",
                   "edge_cost_distribution_by_g.csv"]
    audit_paths = [os.path.join(RUN, "summaries", name) for name in audit_names]
    if not all(os.path.isfile(path) and os.path.getsize(path) > 0 for path in audit_paths):
        raise RuntimeError("Pre-grid input audit did not produce its required outputs")
    write_json(os.path.join(RUN, "config", "pre_grid_audit_gate.json"),
               {"status": "passed", "completed_at": stamp(),
                "required_outputs": {path: sha256(path) for path in audit_paths},
                "statement": "Reward, mapping, component and effective-cost audits completed before PCSF grid"})
    script_path=os.path.abspath(__file__)
    shutil_copy=__import__("shutil").copy2
    shutil_copy(script_path,os.path.join(RUN,"scripts",os.path.basename(script_path)))
    script_hash=sha256(script_path)
    cfgpath=os.path.join(RUN,"config","run_config.json")
    if os.path.isfile(cfgpath):
        existing=json.load(open(cfgpath,encoding="utf-8"))
        frozen_script_hash = existing.get("execution_script_sha256")
        if frozen_script_hash != script_hash:
            patch_path=os.path.join(RUN,"config","compatibility_patch_registry.json")
            if not os.path.isfile(patch_path):
                raise RuntimeError("Execution script hash differs from frozen run config and no reviewed compatibility patch is registered")
            patch=json.load(open(patch_path,encoding="utf-8"))
            if (patch.get("parent_script_sha256") != frozen_script_hash or
                    patch.get("patched_script_sha256") != script_hash or
                    patch.get("scope") != "control metric propagation and summary loop arity only; no input, objective, solver, or parameter changes" or
                    args.phase not in ("controls", "summarize")):
                raise RuntimeError("Registered compatibility patch does not match this script, frozen parent, or permitted phase")
        # Keep the original analysis config hash stable so completed base-grid
        # manifests remain verifiable after this output-field compatibility fix.
        config_hash=hashlib.sha256((sha256(cfgpath)+str(frozen_script_hash)).encode("ascii")).hexdigest()
    else:
        cfg_sha=save_grid_config(data,input_hashes,script_hash)
        config_hash=hashlib.sha256((cfg_sha+script_hash).encode("ascii")).hexdigest()
    adj=make_adjacency(data["edge_uv"],len(data["nodes"]))
    log("preflight complete: {} nodes, {} edges, {} PPI components; 10 seeds; phase={}".format(
        len(data["nodes"]),len(data["edge_uv"]),data["components"],args.phase))
    if args.phase in ("grid","all"):
        run_grid(data,config_hash)
    if args.phase in ("grid","controls","summarize","all"):
        write_best_npz_for_base()
    if args.phase in ("controls","all"):
        run_controls(data,config_hash)
    if args.phase in ("summarize","all"):
        summarize_all(data)
    phase_status,phase_details=inspect_phase_status(args.phase)
    state_path=os.path.join(RUN,"config","run_state.json")
    old_state=json.load(open(state_path,encoding="utf-8")) if os.path.isfile(state_path) else {}
    new_state={"status":phase_status,"last_attempted_phase":args.phase,
               "last_phase_checked_at":stamp(),"phase_details":phase_details,
               "config_sha256":sha256(cfgpath),"execution_script_sha256":script_hash}
    if phase_status=="complete":
        new_state["last_completed_phase"]=args.phase
        new_state["last_phase_completed_at"]=stamp()
    elif old_state.get("last_completed_phase"):
        new_state["last_completed_phase"]=old_state["last_completed_phase"]
        new_state["last_phase_completed_at"]=old_state.get("last_phase_completed_at")
    write_json(state_path,new_state)
    log("phase {}: {}".format(args.phase,phase_status))


def summarize_all(data):
    node_sets,edge_sets,run_rows=summarize_grid(data)
    # Main eligible parameter table and node-level frequency denominator fixed at 10.
    param_rows=[];freq_rows=[];freq_seed_rows=[];seed_pair=[];neighbor_rows=[];edge_frequency_rows=[]
    run_status_lookup={(r.get("parameter_id"),int(r.get("model_seed",-1))):r.get("status") for r in run_rows}
    for b,w,g in itertools.product(B_VALUES,W_VALUES,G_VALUES):
        pid=param_id(b,w,g)
        vals=[next((r for r in run_rows if r.get("parameter_id")==pid and int(r.get("model_seed",-1))==s),{}) for s in SEEDS]
        valid_vals=[v for v in vals if v.get("status")=="complete"]
        n_conn=sum(int(v.get("n_real_skeleton_edges") or 0)>0 for v in valid_vals)
        n_nodes=[int(v.get("n_real_nodes") or 0) for v in valid_vals]
        param_rows.append({"parameter_id":pid,"b":b,"w":w,"g":g,"completed_seeds":sum(v.get("status")=="complete" for v in vals),
            "failed_or_missing_seeds":sum(v.get("status")!="complete" for v in vals),
            "seeds_with_real_backbone_edges":n_conn,"mean_nodes":mean_or_na(n_nodes),"median_nodes":median_or_na(n_nodes),
            "mean_skeleton_edges":mean_or_na([float(v.get("n_real_skeleton_edges")) for v in valid_vals]),
            "mean_SHAP_mass_retained":mean_or_na([v.get("retained_mapped_SHAP_mass_fraction") for v in valid_vals]),
            "mean_connector_fraction":mean_or_na([v.get("connector_fraction") for v in valid_vals]),
            "mean_isolated_fraction":mean_or_na([v.get("isolated_fraction") for v in valid_vals]),
            "mean_largest_component_fraction":mean_or_na([v.get("largest_component_fraction") for v in valid_vals])})
        all_nodes=set().union(*[node_sets.get((pid,s),set()) for s in SEEDS])
        for ix in sorted(all_nodes):
            eid=data["nodes"][ix]
            completed_seeds=[s for s in SEEDS if run_status_lookup.get((pid,s))=="complete"]
            selected=[ix in node_sets.get((pid,s),set()) for s in completed_seeds]
            backbone=[]
            for s in completed_seeds:
                selected_edges=edge_sets.get((pid,s),set())
                nonisolated=any(ix in data["edge_uv"][int(ei)] for ei in selected_edges)
                backbone.append(ix in node_sets.get((pid,s),set()) and nonisolated)
            selected_count=sum(selected)
            backbone_count=sum(backbone)
            frequency_complete=(len(completed_seeds)==len(SEEDS))
            shap_values=[data["raw_shap"][s][eid] for s in SEEDS if eid in data["raw_shap"][s]]
            evidence_seeds=[s for s in SEEDS if eid in data["evidence"][s]]
            connector_seeds=[s for s in SEEDS if eid not in data["evidence"][s]]
            freq_rows.append({"parameter_id":pid,"Entrez":eid,"selected_count_observed":selected_count,
                "completed_seed_n":len(completed_seeds),"frequency_denominator":10,
                "selection_frequency_of_10":selected_count/10.0 if frequency_complete else None,
                "backbone_node_count_observed":backbone_count,
                "backbone_frequency_of_10":backbone_count/10.0 if frequency_complete else None,
                "frequency_status":"complete_10_seed_denominator" if frequency_complete else "incomplete_runs_frequency_not_estimated",
                "role_any_seed": "model_evidence_in_any_seed" if evidence_seeds else "zero_reward_connector_in_all_seeds",
                "evidence_seed_list": "|".join(map(str,evidence_seeds)),
                "connector_seed_list": "|".join(map(str,connector_seeds)),
                "raw_mean_abs_SHAP_median_across_model_seeds":float(np.median(shap_values)) if shap_values else None,
                "background_degree":int(data["degree"][ix]),
                "majority_ge6_of10":(selected_count>=6) if frequency_complete else None,
                "core_10_of10":(selected_count==10) if frequency_complete else None})
            for seed in SEEDS:
                seed_status=run_status_lookup.get((pid,seed),"missing")
                freq_seed_rows.append({"parameter_id":pid,"model_seed":seed,"Entrez":eid,
                    "run_status":seed_status,
                    "selected":bool(ix in node_sets.get((pid,seed),set())) if seed_status=="complete" else None,
                    "model_evidence_for_seed":bool(eid in data["evidence"][seed]),
                    "role_for_seed":"model_evidence" if eid in data["evidence"][seed] else "zero_reward_PPI_connector",
                    "raw_mean_abs_SHAP":data["raw_shap"][seed].get(eid),
                    "background_degree":int(data["degree"][ix])})
        selected_union_edges=set().union(*[edge_sets.get((pid,s),set()) for s in SEEDS])
        for ei in sorted(selected_union_edges):
            completed_seed_n=sum(run_status_lookup.get((pid,s))=="complete" for s in SEEDS)
            observed_count=sum(ei in edge_sets.get((pid,s),set()) for s in SEEDS if run_status_lookup.get((pid,s))=="complete")
            edge_frequency_rows.append({"parameter_id":pid,"edge_index":ei,
                "Entrez_A":data["nodes"][int(data["edge_uv"][ei,0])],
                "Entrez_B":data["nodes"][int(data["edge_uv"][ei,1])],
                "selected_count_observed":observed_count,
                "frequency_denominator":10,"completed_seed_n":completed_seed_n,
                "frequency_of_10":observed_count/10.0 if completed_seed_n==len(SEEDS) else None,
                "frequency_status":"complete_10_seed_denominator" if completed_seed_n==len(SEEDS) else "incomplete_runs_frequency_not_estimated"})
        for i,s1 in enumerate(SEEDS):
            for s2 in SEEDS[i+1:]:
                n1=node_sets.get((pid,s1),set());n2=node_sets.get((pid,s2),set())
                e1={x for x in n1 if data["nodes"][x] in data["evidence"][s1]};e2={x for x in n2 if data["nodes"][x] in data["evidence"][s2]}
                both_complete=run_status_lookup.get((pid,s1))=="complete" and run_status_lookup.get((pid,s2))=="complete"
                seed_pair.append({"parameter_id":pid,"seed_a":s1,"seed_b":s2,
                    "comparison_status":"complete" if both_complete else "failed_or_missing_run",
                    "node_jaccard":jaccard(n1,n2) if both_complete else None,
                    "evidence_node_jaccard":jaccard(e1,e2) if both_complete else None,
                    "skeleton_edge_jaccard":jaccard(edge_sets.get((pid,s1),set()),edge_sets.get((pid,s2),set())) if both_complete else None})
    for pid,seed in node_sets:
        b,w,g=map(int,[pid.split("_")[0][1:],pid.split("_")[1][1:],pid.split("_")[2][1:]])
        for param,values in [("b",B_VALUES),("w",W_VALUES),("g",G_VALUES)]:
            value={"b":b,"w":w,"g":g}[param];ix=values.index(value)
            if ix==0:continue
            neighbor_val=values[ix-1]
            other={"b":b,"w":w,"g":g};other[param]=neighbor_val
            n_pid=param_id(other["b"],other["w"],other["g"])
            if n_pid not in {k[0] for k in node_sets}:continue
            a=node_sets.get((pid,seed),set());bb=node_sets.get((n_pid,seed),set())
            ae={x for x in a if data["nodes"][x] in data["evidence"][seed]};be={x for x in bb if data["nodes"][x] in data["evidence"][seed]}
            both_complete=run_status_lookup.get((pid,seed))=="complete" and run_status_lookup.get((n_pid,seed))=="complete"
            neighbor_rows.append({"model_seed":seed,"parameter_from":n_pid,"parameter_to":pid,
                "changed_parameter":param,"comparison_status":"complete" if both_complete else "failed_or_missing_run",
                "node_jaccard":jaccard(a,bb) if both_complete else None,"evidence_node_jaccard":jaccard(ae,be) if both_complete else None,
                "skeleton_edge_jaccard":jaccard(edge_sets.get((n_pid,seed),set()),edge_sets.get((pid,seed),set())) if both_complete else None})
    write_rows(os.path.join(RUN,"summaries","parameter_grid_summary.csv"),param_rows,list(param_rows[0].keys()))
    write_rows(os.path.join(RUN,"summaries","node_selection_frequency.csv"),freq_rows,
        ["parameter_id","Entrez","selected_count_observed","completed_seed_n","frequency_denominator","selection_frequency_of_10",
         "backbone_node_count_observed","backbone_frequency_of_10","frequency_status","role_any_seed","evidence_seed_list","connector_seed_list",
         "raw_mean_abs_SHAP_median_across_model_seeds","background_degree","majority_ge6_of10","core_10_of10"])
    write_rows(os.path.join(RUN,"summaries","node_selection_frequency_by_seed.csv"),freq_seed_rows,
        ["parameter_id","model_seed","Entrez","run_status","selected","model_evidence_for_seed","role_for_seed","raw_mean_abs_SHAP","background_degree"])
    write_rows(os.path.join(RUN,"summaries","edge_selection_frequency.csv"),edge_frequency_rows,
        ["parameter_id","edge_index","Entrez_A","Entrez_B","selected_count_observed","frequency_denominator","completed_seed_n","frequency_of_10","frequency_status"])
    write_rows(os.path.join(RUN,"summaries","seed_pairwise_similarity.csv"),seed_pair,
        ["parameter_id","seed_a","seed_b","comparison_status","node_jaccard","evidence_node_jaccard","skeleton_edge_jaccard"])
    write_rows(os.path.join(RUN,"summaries","parameter_neighbor_similarity.csv"),neighbor_rows,
        ["model_seed","parameter_from","parameter_to","changed_parameter","comparison_status","node_jaccard","evidence_node_jaccard","skeleton_edge_jaccard"])
    summarize_controls(data)
    write_combined_run_status()
    recommend(data)


def write_combined_run_status():
    rows=[]
    base_path=os.path.join(RUN,"summaries","base_run_status.csv")
    if os.path.isfile(base_path):
        for item in read_csv(base_path).to_dict("records"):
            rows.append({"run_kind":"base","parameter_id":item.get("parameter_id"),
                "model_seed":item.get("model_seed"),"control_type":"","replicate":"",
                "status":item.get("status"),"started_at":item.get("started_at"),
                "finished_at":item.get("finished_at"),"runtime_seconds":"",
                "error":item.get("error","")})
    control_path=os.path.join(RUN,"summaries","control_run_status.csv")
    if os.path.isfile(control_path):
        for item in read_csv(control_path).to_dict("records"):
            rows.append({"run_kind":"control","parameter_id":item.get("parameter_id"),
                "model_seed":item.get("model_seed"),"control_type":item.get("control_type"),
                "replicate":item.get("replicate"),"status":item.get("status"),
                "started_at":item.get("started_at"),"finished_at":item.get("finished_at"),
                "runtime_seconds":item.get("runtime_seconds"),"error":item.get("error","")})
    write_rows(os.path.join(RUN,"summaries","run_status.csv"),rows,
        ["run_kind","parameter_id","model_seed","control_type","replicate","status",
         "started_at","finished_at","runtime_seconds","error"])


def mean_or_na(values):
    clean=[]
    for v in values:
        try:
            x=float(v)
            if math.isfinite(x):clean.append(x)
        except Exception:pass
    return float(np.mean(clean)) if clean else None


def degree_of_node(ix,edges,edge_uv):
    return sum(1 for ei in edges if ix in edge_uv[int(ei)])


def summarize_control_coverage():
    cfgpath=os.path.join(RUN,"config","run_config.json")
    cfg=json.load(open(cfgpath,encoding="utf-8"))
    config_hash=hashlib.sha256((sha256(cfgpath)+str(cfg["execution_script_sha256"])).encode("ascii")).hexdigest()
    eligible_path=os.path.join(RUN,"config","control_eligible_parameters.json")
    eligible=set(json.load(open(eligible_path,encoding="utf-8")).get("eligible_parameter_ids",[])) if os.path.isfile(eligible_path) else set()
    status_path=os.path.join(RUN,"summaries","control_run_status.csv")
    rows=read_csv(status_path).to_dict("records") if os.path.isfile(status_path) else []
    status={(r.get("parameter_id"),int(r.get("model_seed",-1)),r.get("control_type"),int(r.get("replicate",-1))):r.get("status") for r in rows}
    out=[]
    for b,w,g in itertools.product(B_VALUES,W_VALUES,G_VALUES):
        pid=param_id(b,w,g)
        if pid not in eligible:
            out.append({"parameter_id":pid,"control_eligible_by_frozen_rule":False,"expected_replicates":0,
                        "completed_replicates":0,"failed_replicates":0,"invalid_artifacts":0,"missing_replicates":0,
                        "control_status":"not_applicable_by_prespecified_rule"})
            continue
        expected={(pid,seed,kind,rep) for seed in SEEDS for kind in CTRL_TYPES for rep in range(1,REPEATS+1)}
        complete=failed=invalid=0
        for key in expected:
            state=status.get(key)
            if state=="complete":
                if is_valid_done(run_manifest_path(key[0],key[1],key[2],key[3]),config_hash):complete+=1
                else:invalid+=1
            elif state=="failed":failed+=1
        missing=len(expected)-complete-failed-invalid
        out.append({"parameter_id":pid,"control_eligible_by_frozen_rule":True,"expected_replicates":len(expected),
                    "completed_replicates":complete,"failed_replicates":failed,"invalid_artifacts":invalid,"missing_replicates":missing,
                    "control_status":"complete" if complete==len(expected) else "partial"})
    write_rows(os.path.join(RUN,"summaries","control_coverage_by_parameter.csv"),out,
        ["parameter_id","control_eligible_by_frozen_rule","expected_replicates","completed_replicates",
         "failed_replicates","invalid_artifacts","missing_replicates","control_status"])
    return out

def summarize_controls(data):
    summarize_control_coverage()
    rows=[]; verified_manifests=[]; null_node_rows=[]; cost_node_rows=[]
    cfgpath=os.path.join(RUN,"config","run_config.json")
    cfg=json.load(open(cfgpath,encoding="utf-8"))
    config_hash=hashlib.sha256((sha256(cfgpath)+str(cfg["execution_script_sha256"])).encode("ascii")).hexdigest()
    status_path=os.path.join(RUN,"summaries","control_run_status.csv")
    statuses=read_csv(status_path) if os.path.isfile(status_path) else pd.DataFrame(
        columns=["parameter_id","model_seed","control_type","replicate","status"])
    denominator_rows=[]
    eligible_path=os.path.join(RUN,"config","control_eligible_parameters.json")
    eligible=json.load(open(eligible_path,encoding="utf-8")).get("eligible_parameter_ids",[]) if os.path.isfile(eligible_path) else []
    for pid in eligible:
        for seed in SEEDS:
            for kind in CTRL_TYPES:
                grp=statuses[(statuses.parameter_id==pid)&(statuses.model_seed.astype(str)==str(seed))&
                             (statuses.control_type==kind)] if not statuses.empty else pd.DataFrame()
                failed=int((grp.status=="failed").sum()) if not grp.empty else 0
                invalid=0; successful=0
                if not grp.empty:
                    for rep_row in grp.to_dict("records"):
                        if rep_row.get("status")!="complete":continue
                        mpath=run_manifest_path(pid,seed,kind,int(rep_row["replicate"]))
                        if is_valid_done(mpath,config_hash):
                            successful+=1
                        else:
                            invalid+=1
                missing=max(0,REPEATS-successful-failed-invalid)
                denominator_rows.append({"parameter_id":pid,"model_seed":seed,"control_type":kind,
                    "successful_replicates":successful,"expected_replicates":REPEATS,
                    "failed_replicates":failed,"invalid_artifacts":invalid,"missing_replicates":missing,
                    "zero_frequency_encoding":"PPI nodes omitted from sparse per-node selection table have count zero among hash-verified successful replicates"})
    write_rows(os.path.join(RUN,"summaries","null_replicate_denominators.csv"),denominator_rows,
        ["parameter_id","model_seed","control_type","successful_replicates","expected_replicates",
         "failed_replicates","invalid_artifacts","missing_replicates","zero_frequency_encoding"])
    # A complete status row is not sufficient: every summarized result must
    # still match the configuration hash and all output hashes in its manifest.
    for row in statuses.to_dict("records"):
        if row.get("status")!="complete":continue
        mpath=run_manifest_path(row["parameter_id"],int(row["model_seed"]),row["control_type"],int(row["replicate"]))
        if not is_valid_done(mpath,config_hash):continue
        m=json.load(open(mpath,encoding="utf-8"))
        rows.append(m)
        verified_manifests.append(m)
    group_keys=["parameter_id","control_type"]
    out=[]
    if rows:
        df=pd.DataFrame(rows)
        for keys,grp in df.groupby(group_keys):
            out.append({"parameter_id":keys[0],"control_type":keys[1],"n_runs":len(grp),
                "n_model_seeds":grp.model_seed.nunique(),"mean_node_jaccard_vs_base":mean_or_na(grp.node_jaccard_vs_base),
                "median_node_jaccard_vs_base":median_or_na(grp.node_jaccard_vs_base),
                "mean_edge_jaccard_vs_base":mean_or_na(grp.skeleton_edge_jaccard_vs_base),
                "mean_real_nodes":mean_or_na(grp.real_node_count),"mean_real_edges":mean_or_na(grp.real_edge_count),
                "mean_connector_count":mean_or_na(grp.connector_count),"mean_retained_SHAP_fraction":mean_or_na(grp.retained_SHAP_fraction),
                "interpretation":"exploratory frequencies, not p-values/FDR or independent model replicates"})
    write_rows(os.path.join(RUN,"summaries","perturbation_summary.csv"),out,
       ["parameter_id","control_type","n_runs","n_model_seeds","mean_node_jaccard_vs_base","median_node_jaccard_vs_base",
        "mean_edge_jaccard_vs_base","mean_real_nodes","mean_real_edges","mean_connector_count","mean_retained_SHAP_fraction","interpretation"])
    # Per-node inclusion counts are sparse; omitted PPI nodes have count zero
    # among the hash-verified successful repeats in that denominator group.
    for m in verified_manifests:
        arrpath=os.path.join(os.path.dirname(run_manifest_path(m["parameter_id"],int(m["model_seed"]),m["control_type"],int(m["replicate"]))), "replicate_arrays.npz")
        with np.load(arrpath) as z:
            chosen=set(map(int,z["selected_nodes"]))
        for ix in chosen:
            target=cost_node_rows if m["control_type"]=="edge_cost" else null_node_rows
            target.append({"parameter_id":m["parameter_id"],"model_seed":m["model_seed"],
                "control_type":m["control_type"],"replicate":m["replicate"],"node_index":ix})
    denominator_lookup={(r["parameter_id"],str(r["model_seed"]),r["control_type"]):r for r in denominator_rows}
    def write_node_frequency(rows_in,path):
        fields=["parameter_id","model_seed","control_type","Entrez","selected_count","replicate_denominator",
                "expected_replicates","failed_replicates","invalid_artifacts","missing_replicates",
                "selection_frequency","background_degree","role","zero_frequency_encoding"]
        if not rows_in:
            write_rows(path,[],fields)
            return
        nd=pd.DataFrame(rows_in)
        grouped=nd.groupby(["parameter_id","model_seed","control_type","node_index"]).size().reset_index(name="selected_count")
        out_rows=[]
        for rec in grouped.to_dict("records"):
            key=(rec["parameter_id"],str(rec["model_seed"]),rec["control_type"])
            den=denominator_lookup.get(key,{})
            n_success=int(den.get("successful_replicates",0))
            ix=int(rec["node_index"])
            out_rows.append({"parameter_id":rec["parameter_id"],"model_seed":rec["model_seed"],
                "control_type":rec["control_type"],"Entrez":data["nodes"][ix],
                "selected_count":int(rec["selected_count"]),"replicate_denominator":n_success,
                "expected_replicates":int(den.get("expected_replicates",REPEATS)),
                "failed_replicates":int(den.get("failed_replicates",0)),
                "invalid_artifacts":int(den.get("invalid_artifacts",0)),
                "missing_replicates":int(den.get("missing_replicates",REPEATS-n_success)),
                "selection_frequency":int(rec["selected_count"])/float(n_success) if n_success else None,
                "background_degree":int(data["degree"][ix]),
                "role":"model_evidence" if data["nodes"][ix] in data["evidence"][int(rec["model_seed"])] else "zero_reward_PPI_connector",
                "zero_frequency_encoding":"PPI nodes omitted from this sparse table have count zero among hash-verified successful repeats"})
        write_rows(path,out_rows,fields)
    write_node_frequency(null_node_rows,os.path.join(RUN,"summaries","null_selection_summary.csv"))
    write_node_frequency(cost_node_rows,os.path.join(RUN,"summaries","cost_perturbation_node_selection_frequency.csv"))


def median_or_na(values):
    a=[]
    for v in values:
        try:
            x=float(v)
            if math.isfinite(x):a.append(x)
        except Exception:pass
    return float(np.median(a)) if a else None


def recommend(data):
    p=pandas_read_path(os.path.join(RUN,"summaries","parameter_grid_summary.csv"))
    if p.empty:
        write_rows(os.path.join(RUN,"summaries","parameter_recommendation.csv"),[],["parameter_id","recommendation_status"]);return
    # Multi-criteria comparison; network scale is reported but never used as a
    # direction or composite score. Null similarity is minimized to describe
    # dependence on the real reward locations, not as a significance statistic.
    controls=pandas_read_path(os.path.join(RUN,"summaries","perturbation_summary.csv")) if os.path.isfile(os.path.join(RUN,"summaries","perturbation_summary.csv")) else pd.DataFrame()
    edge_stability=controls[controls.control_type=="edge_cost"][["parameter_id","mean_node_jaccard_vs_base","mean_edge_jaccard_vs_base"]].rename(
        columns={"mean_node_jaccard_vs_base":"cost_perturb_node_jaccard","mean_edge_jaccard_vs_base":"cost_perturb_edge_jaccard"}) if not controls.empty else pd.DataFrame()
    nulls=controls[controls.control_type!="edge_cost"][["parameter_id","control_type","mean_node_jaccard_vs_base"]].pivot(
        index="parameter_id",columns="control_type",values="mean_node_jaccard_vs_base").reset_index()
    nulls=nulls.rename(columns={"reward_shuffle":"ordinary_reward_shuffle_node_jaccard",
                                "degree_stratified_reward_shuffle":"degree_stratified_shuffle_node_jaccard"})
    if "ordinary_reward_shuffle_node_jaccard" not in nulls.columns:
        nulls["ordinary_reward_shuffle_node_jaccard"]=np.nan
    if "degree_stratified_shuffle_node_jaccard" not in nulls.columns:
        nulls["degree_stratified_shuffle_node_jaccard"]=np.nan
    seed_stability=pandas_read_path(os.path.join(RUN,"summaries","seed_pairwise_similarity.csv"))
    seed_stability=seed_stability[seed_stability.comparison_status=="complete"].groupby("parameter_id",as_index=False)[
        ["node_jaccard","evidence_node_jaccard","skeleton_edge_jaccard"]].mean().rename(columns={
            "node_jaccard":"mean_model_seed_node_jaccard",
            "evidence_node_jaccard":"mean_model_seed_evidence_jaccard",
            "skeleton_edge_jaccard":"mean_model_seed_edge_jaccard"})
    coverage_path=os.path.join(RUN,"summaries","control_coverage_by_parameter.csv")
    coverage=pandas_read_path(coverage_path) if os.path.isfile(coverage_path) else pd.DataFrame()
    x=p.merge(edge_stability,on="parameter_id",how="left").merge(nulls,on="parameter_id",how="left").merge(seed_stability,on="parameter_id",how="left")
    if not coverage.empty:
        x=x.merge(coverage,on="parameter_id",how="left")
    x["technical_status"]=np.where(x.completed_seeds==len(SEEDS),"all_10_seeds_complete","incomplete")
    x["connection_status"]=np.where(x.seeds_with_real_backbone_edges>0,"has_real_backbone_in_at_least_one_seed","no_real_backbone_observed")
    if "control_eligible_by_frozen_rule" not in x.columns:
        x["control_eligible_by_frozen_rule"]=False
        x["control_status"]="missing_coverage_record"
    x["pareto_scope_status"]=[("eligible_complete" if bool(e) and cs=="complete" and ts=="all_10_seeds_complete" else
                                "eligible_incomplete" if bool(e) else
                                "excluded_no_backbone_by_prespecified_rule")
                               for e,cs,ts in zip(x.control_eligible_by_frozen_rule,x.control_status,x.technical_status)]
    pareto_directions={
        "seeds_with_real_backbone_edges":1,
        "mean_SHAP_mass_retained":1,
        "mean_model_seed_node_jaccard":1,
        "mean_model_seed_evidence_jaccard":1,
        "cost_perturb_node_jaccard":1,
        "cost_perturb_edge_jaccard":1,
        "ordinary_reward_shuffle_node_jaccard":-1,
        "degree_stratified_shuffle_node_jaccard":-1,
        "mean_connector_fraction":-1,
        "mean_isolated_fraction":-1,
    }
    usable=[]
    for row in x.to_dict("records"):
        if row.get("pareto_scope_status")!="eligible_complete":
            usable.append(None)
            continue
        vals=[]
        for metric,direction in pareto_directions.items():
            value=row.get(metric)
            try:value=float(value)
            except Exception:value=np.nan
            if not np.isfinite(value):break
            vals.append(value*direction)
        usable.append(vals if len(vals)==len(pareto_directions) else None)
    x["pareto_dimensions_used"]="|".join(pareto_directions.keys())
    x["pareto_comparable"]=[v is not None for v in usable]
    pareto=[]
    for i,vi in enumerate(usable):
        if vi is None:
            pareto.append(None)
            continue
        dominated=False
        for j,vj in enumerate(usable):
            if i==j or vj is None:continue
            if all(a>=b-1e-12 for a,b in zip(vj,vi)) and any(a>b+1e-12 for a,b in zip(vj,vi)):
                dominated=True;break
        pareto.append(not dominated)
    x["pareto_nondominated"] = pareto
    x["composite_score"]="not_calculated"
    write_rows(os.path.join(RUN,"summaries","parameter_recommendation.csv"),x.to_dict("records"),list(x.columns))
    # Figure-source long table retains all profile dimensions without a composite score.
    write_rows(os.path.join(RUN,"figure_source","parameter_comparison_long.csv"),
        [{"parameter_id":r["parameter_id"],"metric":k,"value":r.get(k)} for r in x.to_dict("records")
         for k in x.columns if k not in ("parameter_id","b","w","g","technical_status","connection_status")],
        ["parameter_id","metric","value"])


if __name__ == "__main__":
    main()

























