#!/usr/bin/env python3
"""Build frozen PTG14 pan14 PPI+PCSF runs; no training or SHAP recomputation."""
from __future__ import print_function

import csv
import hashlib
import json
import math
import os
import sys
import time
import traceback
from collections import defaultdict
from datetime import datetime

import numpy as np
import pandas as pd
import pkg_resources
from OmicsIntegrator import graph as oi_graph

OUT_ROOT = r"E:\micro_output\分析审阅中心\待审阅\Fig5_pan14_PPI_PCSF更新"
FIG5_ROOT = r"E:\micro_output\分析审阅中心\待审阅\Fig5_主分析与流程图收尾"
MODEL_ROOT = r"E:\micro_output\分析审阅中心\待审阅\PTG02_CD_分支疾病与调控候选\ptg14_training_v2"
MANIFEST = os.path.join(FIG5_ROOT, "source_data", "PTG14_seed42_51_majority28_seed_manifest.csv")
PPI_PATH = r"D:\data\microglia_new\input\drug_repurposing_38core\mentor_20260923\Network-1.txt"
FIG6_ROOT = r"E:\micro_output\分析审阅中心\待审阅\Fig6_疾病vsCtrl_SHAP_PCSF模块"
OFFICIAL_MAP = os.path.join(FIG6_ROOT, "data", "orgHs_official_symbol_entrez.tsv")
ALIAS_MAP = os.path.join(FIG6_ROOT, "data", "orgHs_alias_entrez.tsv")
MAP_META = os.path.join(FIG6_ROOT, "data", "orgHs_mapping_resource_manifest.json")
SEEDS = list(range(42, 52))
EXPECTED_CONTRACT = "a02f5b0437435e2029a8b1a4ab2f3896f275208f97804983c3fbdeff8d20a597"
EXPECTED_FEATURES = 4843
OI_PARAMS = {"w": 3, "b": 5, "g": 5, "edge_noise": 0, "dummy_mode": "terminals", "seed": 42, "skip_checks": False}
PRUNING = "strong"
SCRIPT_PATH = os.path.abspath(__file__)

for sub in ("config", "data", "networks", "source_data", "figures", "scripts", "logs"):
    os.makedirs(os.path.join(OUT_ROOT, sub), exist_ok=True)
LOG = os.path.join(OUT_ROOT, "logs", "run_pan14_pcsf.log")


def log(message):
    line = "[{}] {}".format(datetime.now().astimezone().isoformat(timespec="seconds"), message)
    print(line, flush=True)
    with open(LOG, "a", encoding="utf-8") as handle:
        handle.write(line + "\n")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def json_write(path, value):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    os.replace(tmp, path)


def csv_write(path, rows, fields, delimiter=","):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter, extrasaction="ignore")
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


def pandas_read_csv(path, **kwargs):
    """Read Unicode Windows paths through a Python handle for the frozen pandas build."""
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        return pd.read_csv(handle, **kwargs)


def read_mapping(path, key_name):
    result = defaultdict(set)
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row.get(key_name) or "").strip()
            entrez = (row.get("ENTREZID") or "").strip()
            if key and entrez:
                result[key].add(entrez)
    return result


def sort_id(value):
    return (0, int(value)) if str(value).isdigit() else (1, str(value))


def map_symbol(symbol, official, alias):
    ids = official.get(symbol, set())
    if len(ids) == 1:
        return "official_unique", next(iter(ids)), "official_symbol", "|".join(sorted(ids))
    if len(ids) > 1:
        return "official_ambiguous", "", "official_symbol", "|".join(sorted(ids))
    ids = alias.get(symbol, set())
    if len(ids) == 1:
        return "alias_unique", next(iter(ids)), "alias", "|".join(sorted(ids))
    if len(ids) > 1:
        return "alias_ambiguous", "", "alias", "|".join(sorted(ids))
    return "unmapped", "", "none", ""


def load_ppi():
    raw = []
    mismatches = 0
    with open(PPI_PATH, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Protein_A_Entrez_ID", "Protein_B_Entrez_ID", "data_source.s.", "number"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Mentor PPI columns differ: {}".format(reader.fieldnames))
        for record_id, row in enumerate(reader, 1):
            a = (row.get("Protein_A_Entrez_ID") or "").strip()
            b = (row.get("Protein_B_Entrez_ID") or "").strip()
            source = row.get("data_source.s.") or ""
            source_items = [x.strip() for x in source.split(",") if x.strip()]
            try:
                number = int(float((row.get("number") or "").strip()))
            except Exception:
                number = None
            source_match = "" if number is None else str(len(set(source_items)) == number).lower()
            if number is not None and len(set(source_items)) != number:
                mismatches += 1
            if number is None:
                status = "invalid_number"
            elif number < 2:
                status = "excluded_number_lt_2"
            elif not a or not b:
                status = "invalid_endpoint"
            elif a == b:
                status = "excluded_self_loop"
            else:
                status = "eligible_for_dedup"
            raw.append({"record_id": record_id, "A": a, "B": b, "number_raw": row.get("number", ""),
                        "number": number, "source_text": source, "unique_source_count": len(set(source_items)),
                        "source_count_matches_number": source_match, "disposition": status})
    grouped = defaultdict(list)
    for row in raw:
        if row["disposition"] != "eligible_for_dedup":
            continue
        a, b = sorted((row["A"], row["B"]), key=sort_id)
        grouped[(a, b)].append(row)
    edges = []
    degree = defaultdict(int)
    nodes = set()
    for (a, b), records in sorted(grouped.items(), key=lambda item: (sort_id(item[0][0]), sort_id(item[0][1]))):
        sources = sorted({s.strip() for r in records for s in r["source_text"].split(",") if s.strip()})
        numbers = sorted({int(r["number"]) for r in records})
        edges.append({"Protein_A_Entrez_ID": a, "Protein_B_Entrez_ID": b, "base_cost": 1.0,
                      "number_values": ";".join(map(str, numbers)), "data_source_union": ";".join(sources),
                      "raw_source_texts": " || ".join(sorted({r["source_text"] for r in records})),
                      "raw_record_count": len(records), "raw_record_ids": ";".join(str(r["record_id"]) for r in records)})
        nodes.update((a, b)); degree[a] += 1; degree[b] += 1
    # Deterministic component IDs over the complete eligible PPI; no component is dropped.
    adjacency = defaultdict(set)
    for edge in edges:
        a, b = edge["Protein_A_Entrez_ID"], edge["Protein_B_Entrez_ID"]
        adjacency[a].add(b); adjacency[b].add(a)
    components = []
    remaining = set(nodes)
    while remaining:
        root = min(remaining, key=sort_id)
        stack = [root]; seen = set()
        while stack:
            current = stack.pop()
            if current in seen:
                continue
            seen.add(current)
            stack.extend(adjacency[current] - seen)
        remaining -= seen
        components.append(sorted(seen, key=sort_id))
    components.sort(key=lambda values: (-len(values), sort_id(values[0])))
    comp_id = {node: "PPI_C{:05d}".format(i) for i, group in enumerate(components, 1) for node in group}
    raw_fields = ["record_id", "A", "B", "number_raw", "number", "source_text", "unique_source_count",
                  "source_count_matches_number", "disposition"]
    csv_write(os.path.join(OUT_ROOT, "data", "ppi_raw_record_audit.csv"), raw, raw_fields)
    edge_fields = ["Protein_A_Entrez_ID", "Protein_B_Entrez_ID", "base_cost", "number_values", "data_source_union",
                   "raw_source_texts", "raw_record_count", "raw_record_ids"]
    csv_write(os.path.join(OUT_ROOT, "data", "ppi_edges_ge2.tsv"), edges, edge_fields, delimiter="\t")
    node_rows = [{"Entrez": node, "background_degree": degree[node], "component_id": comp_id[node]}
                 for node in sorted(nodes, key=sort_id)]
    csv_write(os.path.join(OUT_ROOT, "data", "ppi_nodes_ge2.csv"), node_rows,
              ["Entrez", "background_degree", "component_id"])
    audit_rows = []
    for key, value in [("raw_record_count", len(raw)),
                       ("number_ge2_record_count", sum(r["number"] is not None and r["number"] >= 2 for r in raw)),
                       ("source_count_mismatch_count", mismatches),
                       ("invalid_number_count", sum(r["number"] is None for r in raw)),
                       ("self_loop_ge2_record_count", sum(r["disposition"] == "excluded_self_loop" for r in raw)),
                       ("eligible_nonself_raw_record_count", sum(r["disposition"] == "eligible_for_dedup" for r in raw)),
                       ("deduplicated_edge_count", len(edges)), ("node_count", len(nodes)),
                       ("component_count", len(components))]:
        audit_rows.append({"metric": key, "value": value})
    csv_write(os.path.join(OUT_ROOT, "data", "ppi_filter_summary.csv"), audit_rows, ["metric", "value"])
    return edges, nodes, degree, comp_id, raw


def safe_outputs(seed):
    folder = os.path.join(OUT_ROOT, "networks", "seed_{}".format(seed))
    names = ["gene_prizes.csv", "all_gene_mapping_audit.csv", "selected_nodes.csv",
             "selected_forest_edges.csv", "selected_induced_ppi_edges.csv", "oi_input_nodes.csv",
             "oi_input_edges.tsv", "oi_positive_prizes.tsv", "solver_call_audit.json", "run_manifest.json"]
    return folder, [os.path.join(folder, name) for name in names]


def main():
    os.environ.setdefault("OMP_NUM_THREADS", "8")
    oi_version = pkg_resources.get_distribution("OmicsIntegrator").version
    pcst_version = pkg_resources.get_distribution("pcst_fast").version
    if oi_version != "2.3.10" or pcst_version != "1.0.7":
        raise RuntimeError("Installed solver versions differ from frozen Fig6 standard OI2.3.10/pcst_fast1.0.7: {} / {}".format(oi_version, pcst_version))
    official = read_mapping(OFFICIAL_MAP, "SYMBOL")
    alias = read_mapping(ALIAS_MAP, "ALIAS")
    official_by_eid = defaultdict(set)
    for symbol, ids in official.items():
        for entrez in ids:
            official_by_eid[entrez].add(symbol)
    edges, ppi_nodes, degree, components, raw_ppi = load_ppi()
    ppi_edge_index = {(e["Protein_A_Entrez_ID"], e["Protein_B_Entrez_ID"]): i for i, e in enumerate(edges)}
    genes = pandas_read_csv(os.path.join(MODEL_ROOT, "data", "genes.csv"), dtype=str)
    if genes.columns.tolist()[:2] != ["gene_id", "PTG"] or len(genes) != EXPECTED_FEATURES or genes.gene_id.nunique() != EXPECTED_FEATURES:
        raise ValueError("Frozen PTG14 data/genes.csv identity/order/features validation failed")
    genes_order = genes.gene_id.astype(str).tolist()
    cell_manifest = pandas_read_csv(os.path.join(MODEL_ROOT, "data", "cells.csv"), dtype=str)
    contract = json.load(open(os.path.join(MODEL_ROOT, "data", "contract.json"), encoding="utf-8"))
    export = json.load(open(os.path.join(MODEL_ROOT, "data", "export_complete.json"), encoding="utf-8"))
    if contract.get("id") != EXPECTED_CONTRACT or export.get("n_genes") != EXPECTED_FEATURES:
        raise ValueError("Frozen PTG14 contract/export metadata mismatch")
    if not set(cell_manifest.get("Group3", [])) >= {"Neurological", "Tumor"}:
        raise ValueError("Frozen cells manifest lacks Neurological/Tumor identities")
    if set(cell_manifest.loc[cell_manifest.split == "ctrl_reference", "Group3"].astype(str)) != {"Ctrl"}:
        raise ValueError("Frozen Ctrl reference cells are not isolated in the ctrl_reference split")
    model_splits = cell_manifest[cell_manifest.split.isin(["train", "val", "holdout"])]
    if not set(model_splits.Group3.astype(str).unique()).issubset({"Neurological", "Tumor"}):
        raise ValueError("A Ctrl or other class is included in the frozen model train/validation/holdout splits")
    seed_manifest = pandas_read_csv(MANIFEST, dtype=str)
    if set(seed_manifest.seed.astype(int)) != set(SEEDS) or any(seed_manifest.contract != EXPECTED_CONTRACT) or any(seed_manifest.n_features.astype(int) != EXPECTED_FEATURES):
        raise ValueError("Original Fig5 seed manifest does not freeze all ten PTG14 runs")
    seed_info = {}
    for seed in SEEDS:
        run_dir = os.path.join(MODEL_ROOT, "runs", "seed_{}".format(seed))
        paths = {name: os.path.join(run_dir, name) for name in ("weights.csv", "candidates.csv", "complete.json", "model.pt")}
        if not all(os.path.isfile(p) for p in paths.values()):
            raise FileNotFoundError("Seed {} missing actual full-run inputs".format(seed))
        complete = json.load(open(paths["complete.json"], encoding="utf-8"))
        registry_row = seed_manifest[seed_manifest.seed.astype(int) == seed].iloc[0]
        weight_hash = sha256(paths["weights.csv"])
        candidate_hash = sha256(paths["candidates.csv"])
        model_hash = sha256(paths["model.pt"])
        if complete.get("seed") != seed or complete.get("contract") != EXPECTED_CONTRACT:
            raise ValueError("Seed {} complete.json identity/contract mismatch".format(seed))
        outputs = complete.get("outputs", {})
        if weight_hash.lower() != str(outputs.get("weights.csv", "")).lower() or candidate_hash.lower() != str(outputs.get("candidates.csv", "")).lower() or model_hash.lower() != str(outputs.get("model.pt", "")).lower():
            raise ValueError("Seed {} actual weights/candidates/checkpoint hash mismatch".format(seed))
        if weight_hash.lower() != str(registry_row.weights_sha256).lower() or candidate_hash.lower() != str(registry_row.candidates_sha256).lower():
            raise ValueError("Seed {} actual full data hashes mismatch original seed registry".format(seed))
        weights = pandas_read_csv(paths["weights.csv"], dtype={"gene_id": str, "PTG": str})
        if weights.gene_id.tolist() != genes_order or weights.PTG.tolist() != genes.PTG.tolist():
            raise ValueError("Seed {} gene order/PTG differs from frozen full gene manifest".format(seed))
        if not {"mean_abs_SHAP", "weight"}.issubset(weights.columns):
            raise ValueError("Seed {} weights missing full SHAP fields".format(seed))
        shap = np.asarray(pd.to_numeric(weights.mean_abs_SHAP, errors="coerce"), dtype=float)
        norm = np.asarray(pd.to_numeric(weights.weight, errors="coerce"), dtype=float)
        if not np.isfinite(shap).all() or np.any(shap < 0) or not np.isfinite(norm).all():
            raise ValueError("Seed {} has invalid frozen SHAP values".format(seed))
        expected_norm = shap / shap.sum()
        if not np.allclose(norm, expected_norm, rtol=1e-12, atol=1e-14):
            raise ValueError("Seed {} original normalized weight formula mismatch".format(seed))
        seed_info[seed] = {"run_dir": run_dir, "paths": paths, "complete": complete,
                           "weights": weights, "weights_sha256": weight_hash,
                           "candidates_sha256": candidate_hash, "model_sha256": model_hash,
                           "shap_sum": float(shap.sum()), "shap_max": float(shap.max())}

    # Save immutable input/config provenance before any solver call.
    input_paths = [MANIFEST, PPI_PATH, OFFICIAL_MAP, ALIAS_MAP, MAP_META,
                   os.path.join(MODEL_ROOT, "data", "contract.json"),
                   os.path.join(MODEL_ROOT, "data", "export_complete.json"),
                   os.path.join(MODEL_ROOT, "data", "genes.csv"),
                   os.path.join(MODEL_ROOT, "data", "cells.csv")]
    input_rows = [{"path": p, "sha256": sha256(p), "size_bytes": os.path.getsize(p),
                   "modified_local": datetime.fromtimestamp(os.path.getmtime(p)).astimezone().isoformat(timespec="seconds")}
                  for p in input_paths]
    for seed in SEEDS:
        for name in ("weights.csv", "candidates.csv", "complete.json", "model.pt"):
            p = seed_info[seed]["paths"][name]
            input_rows.append({"path": p, "sha256": sha256(p), "size_bytes": os.path.getsize(p),
                               "modified_local": datetime.fromtimestamp(os.path.getmtime(p)).astimezone().isoformat(timespec="seconds")})
    csv_write(os.path.join(OUT_ROOT, "source_data", "input_sha256_manifest.csv"), input_rows,
              ["path", "sha256", "size_bytes", "modified_local"])
    config = {"analysis": "Fig5 pan14 PTG14 ten-seed mentor PPI + PCSF update", "status": "running",
              "model_identity": {"contract": EXPECTED_CONTRACT, "features": EXPECTED_FEATURES,
                                 "group_labels": ["Neurological", "Tumor"], "control_role": "reference only; no Ctrl class in this pan14 classifier",
                                 "seeds": SEEDS, "input_transform": "RNA data layer log1p normalized expression; frozen PTG14 seed model input",
                                 "frozen_contract_config": contract.get("config"),
                                 "cell_split_group_counts": {str(k): int(v) for k, v in cell_manifest.groupby(["split", "Group3"]).size().items()},
                                 "feature_order_source": "ptg14_training_v2/data/genes.csv"},
              "seed_manifest": MANIFEST, "source_hashes": {r["path"]: r["sha256"] for r in input_rows},
              "ppi": {"source": PPI_PATH, "filter": "number >= 2", "graph": "undirected; self loops removed; unordered endpoints deduplicated; all components retained", "base_cost": 1.0,
                      "raw_record_count": len(raw_ppi), "eligible_edge_count": len(edges), "node_count": len(ppi_nodes),
                      "component_count": len(set(components.values())), "source_count_mismatches": sum(r["source_count_matches_number"] == "false" for r in raw_ppi)},
              "mapping": {"resource": MAP_META, "priority": "unique official Symbol; only if absent, unique alias", "ambiguous": "unresolved", "many_symbols_one_entrez": "retain source mappings; node reward is maximum; no duplicate summation"},
              "reward": {"raw": "seed-specific mean_abs_SHAP over all 4,843 features", "prize": "mean_abs_SHAP / max over all 4,843 features before mapping/PPI filtering", "outside_model_nodes": "prize=0; SHAP fields missing"},
              "pcsf": {"implementation": "OmicsIntegrator2 Graph + installed pcst_fast", "b": 5, "w": 3, "g": 5,
                       "dummy_mode": "terminals", "pruning": PRUNING, "seed": 42, "edge_noise": 0,
                       "base_edge_cost": 1.0, "solve_count_per_model_seed": 1,
                       "hub_penalty_formula": "10^g * deg(u)*deg(v) / ((N-deg(u)-1)*(N-deg(v)-1)+deg(u)*deg(v)); effective_cost=base_cost+penalty"},
              "cross_seed_selection": {"majority_threshold": "selected in at least 6 of 10 runs", "stable_core": "selected in 10 of 10 runs", "denominator": 10},
              "software": {"python": sys.version, "omicsintegrator": oi_version, "pcst_fast": pcst_version},
              "notes": ["No retraining, SHAP recomputation, Fig6 model reuse, solver improvement heuristic, parameter tuning, MFVS, Louvain or drug analysis.",
                        "PTG14 model weights are re-used from complete full 4,843-feature seed runs; legacy 28-gene candidates are not used as PCSF input."]}
    json_write(os.path.join(OUT_ROOT, "config", "run_config.json"), config)
    # Persist read-only mapping provenance and resource copies in the new analysis branch.
    for source, target in [(OFFICIAL_MAP, "orgHs_official_symbol_entrez.tsv"),
                           (ALIAS_MAP, "orgHs_alias_entrez.tsv"),
                           (MAP_META, "orgHs_mapping_resource_manifest.json")]:
        dest = os.path.join(OUT_ROOT, "data", target)
        if not os.path.exists(dest):
            with open(source, "rb") as src, open(dest, "wb") as dst:
                dst.write(src.read())

    ppi_node_rows = pandas_read_csv(os.path.join(OUT_ROOT, "data", "ppi_nodes_ge2.csv"), dtype=str).to_dict("records")
    ppi_meta = {row["Entrez"]: row for row in ppi_node_rows}
    all_seed_selected_symbols = []
    all_seed_selected_nodes = []
    run_summaries = []
    mapping_audit_all = []
    mapping_audit_by_seed = {}

    # OI environment is imported above; the exact author solver is kept unchanged.
    for position, seed in enumerate(SEEDS, 1):
        run_dir, output_paths = safe_outputs(seed)
        os.makedirs(run_dir, exist_ok=True)
        existing_manifest_path = os.path.join(run_dir, "run_manifest.json")
        if os.path.isfile(existing_manifest_path):
            prior = json.load(open(existing_manifest_path, encoding="utf-8"))
            current_seed_inputs = {"weights.csv": seed_info[seed]["weights_sha256"],
                "candidates.csv": seed_info[seed]["candidates_sha256"], "model.pt": seed_info[seed]["model_sha256"],
                "mentor_PPI": sha256(PPI_PATH), "official_map": sha256(OFFICIAL_MAP), "alias_map": sha256(ALIAS_MAP),
                "model_contract_json": sha256(os.path.join(MODEL_ROOT, "data", "contract.json")),
                "model_genes_csv": sha256(os.path.join(MODEL_ROOT, "data", "genes.csv")),
                "model_cells_csv": sha256(os.path.join(MODEL_ROOT, "data", "cells.csv"))}
            expected_output_names = {os.path.basename(p) for p in output_paths if os.path.basename(p) != "run_manifest.json"}
            prior_output_hashes = prior.get("output_hashes", {})
            outputs_intact = (set(prior_output_hashes) == expected_output_names and
                all(os.path.basename(name) == name and os.path.isfile(os.path.join(run_dir, name)) and
                    sha256(os.path.join(run_dir, name)) == expected_hash
                    for name, expected_hash in prior_output_hashes.items()))
            current_runner_sha = sha256(SCRIPT_PATH)
            solver_identity_matches = (prior.get("status") == "complete" and prior.get("seed") == seed and
                prior.get("contract") == EXPECTED_CONTRACT and prior.get("input_hashes") == current_seed_inputs and
                prior.get("OI_params") == OI_PARAMS and prior.get("pruning") == PRUNING)
            if solver_identity_matches and outputs_intact:
                if prior.get("runner_script_sha256") != current_runner_sha:
                    log("seed {} solver inputs/parameters and all output hashes verified; postprocessing script hash changed, reusing frozen solver outputs".format(seed))
                log("seed {} existing verified artifacts reused".format(seed))
                node_rows = pandas_read_csv(os.path.join(run_dir, "selected_nodes.csv"), dtype=str).to_dict("records")
                selected_ids = {r["Entrez"] for r in node_rows}
                seed_symbols = {s for r in node_rows for s in (r.get("model_gene_symbols", "") or "").split(";") if s}
                all_seed_selected_symbols.extend([{"gene_id": sym, "seed": seed} for sym in sorted(seed_symbols)])
                all_seed_selected_nodes.extend([{"Entrez": r["Entrez"], "seed": seed,
                    "model_gene_symbols": r.get("model_gene_symbols", ""),
                    "model_evidence_status": r.get("model_evidence_status", "no_SHAP_evidence_PPI_connector")} for r in node_rows])
                mapping_audit_by_seed[seed] = pandas_read_csv(os.path.join(run_dir, "all_gene_mapping_audit.csv"), dtype=str).to_dict("records")
                mapping_audit_all.extend(mapping_audit_by_seed[seed])
                run_summaries.append(prior["summary"])
                continue
            raise RuntimeError("Seed {} prior partial/incompatible artifacts exist; refusing overwrite".format(seed))

        started = time.time()
        weights = seed_info[seed]["weights"].copy()
        weights["gene_id"] = weights.gene_id.astype(str)
        weights["mean_abs_SHAP"] = pd.to_numeric(weights.mean_abs_SHAP, errors="coerce")
        max_shap = float(weights.mean_abs_SHAP.max())
        if not math.isfinite(max_shap) or max_shap <= 0:
            raise ValueError("Seed {} complete model SHAP has no positive maximum".format(seed))
        weights["prize_full_model_max_normalized"] = weights.mean_abs_SHAP / max_shap
        by_entrez = defaultdict(list)
        audit_rows = []
        for row in weights.to_dict("records"):
            symbol = str(row["gene_id"])
            status, entrez, source, candidate_ids = map_symbol(symbol, official, alias)
            in_ppi = entrez in ppi_nodes if entrez else False
            audit = {"seed": seed, "gene_id": symbol, "PTG": row.get("PTG", ""),
                     "mean_abs_SHAP": float(row["mean_abs_SHAP"]),
                     "model_weight_sum_normalized": float(row["weight"]),
                     "prize_full_model_max_normalized": float(row["prize_full_model_max_normalized"]),
                     "mapping_status": status, "chosen_Entrez": entrez,
                     "mapping_source": source, "candidate_Entrez_ids": candidate_ids,
                     "in_qualified_PPI": bool(in_ppi), "PPI_background_degree": ppi_meta[entrez]["background_degree"] if in_ppi else "",
                     "PPI_component_id": ppi_meta[entrez]["component_id"] if in_ppi else ""}
            audit_rows.append(audit)
            if in_ppi:
                by_entrez[entrez].append({"gene_id": symbol, "PTG": row.get("PTG", ""),
                    "mean_abs_SHAP": float(row["mean_abs_SHAP"]), "weight": float(row["weight"]),
                    "prize_normalized": float(row["prize_full_model_max_normalized"]),
                    "mapping_status": status, "mapping_source": source})
        mapping_audit_by_seed[seed] = audit_rows
        mapping_audit_all.extend(audit_rows)
        node_gene_records = {}
        for entrez, records in by_entrez.items():
            winner = sorted(records, key=lambda r: (-r["mean_abs_SHAP"], r["gene_id"]))[0]
            node_gene_records[entrez] = {"model_gene_symbols": ";".join(sorted({r["gene_id"] for r in records})),
                "PTG_groups": ";".join(sorted({r["PTG"] for r in records if r["PTG"]})),
                "mean_abs_SHAP": winner["mean_abs_SHAP"], "model_weight_sum_normalized": winner["weight"],
                "prize": winner["prize_normalized"], "prize_source_symbol": winner["gene_id"],
                "mapping_statuses": ";".join(sorted({r["mapping_status"] for r in records})),
                "model_evidence_status": "positive_model_reward" if winner["mean_abs_SHAP"] > 0 else "measured_zero_reward"}
        # All nodes of the approved PPI are available as zero-reward connectors; only positive rewards are terminals.
        node_order = sorted(ppi_nodes, key=sort_id)
        node_ix = {node: i for i, node in enumerate(node_order)}
        edge_frame = pd.DataFrame(edges)
        encoded_edges = [{"protein1": "ENTREZ:" + e["Protein_A_Entrez_ID"],
                          "protein2": "ENTREZ:" + e["Protein_B_Entrez_ID"], "cost": 1.0} for e in edges]
        solver_edges_path = os.path.join(run_dir, "oi_input_edges.tsv")
        pd.DataFrame(encoded_edges, columns=["protein1", "protein2", "cost"]).to_csv(solver_edges_path, sep="\t", index=False)
        positive_prizes = [{"name": "ENTREZ:" + eid, "prize": float(node_gene_records[eid]["prize"])}
                           for eid in sorted(node_gene_records, key=sort_id) if node_gene_records[eid]["prize"] > 0]
        prize_path = os.path.join(run_dir, "oi_positive_prizes.tsv")
        pd.DataFrame(positive_prizes, columns=["name", "prize"]).to_csv(prize_path, sep="\t", index=False)
        oi_nodes_path = os.path.join(run_dir, "oi_input_nodes.csv")
        csv_write(oi_nodes_path, [{"node_index": node_ix[n], "internal_id": "ENTREZ:" + n,
                    "Entrez": n, "prize_unscaled": node_gene_records.get(n, {}).get("prize", 0.0),
                    "model_evidence_status": node_gene_records.get(n, {}).get("model_evidence_status", "no_SHAP_evidence_PPI_connector")}
                    for n in node_order], ["node_index", "internal_id", "Entrez", "prize_unscaled", "model_evidence_status"])
        # Each seed has a standalone network input, while the shared, audited PPI and OI configuration remain frozen.
        with open(solver_edges_path, "r", encoding="utf-8-sig", newline="") as edge_handle:
            graph = oi_graph.Graph(edge_handle, OI_PARAMS)
        with open(prize_path, "r", encoding="utf-8-sig", newline="") as prize_handle:
            graph.prepare_prizes(prize_handle)
        if len(graph.nodes) != len(ppi_nodes) or len(graph.interactome_dataframe) != len(edges):
            raise RuntimeError("Seed {} OI network arrays differ from approved PPI".format(seed))
        if len(graph.terminals) != len(positive_prizes):
            raise RuntimeError("Seed {} OI terminal count differs from positive prize count".format(seed))
        expected_prizes = {"ENTREZ:" + n: float(node_gene_records.get(n, {}).get("prize", 0.0)) for n in node_order}
        for i, internal_id in enumerate(graph.nodes):
            raw_prize = expected_prizes.get(str(internal_id), 0.0)
            if not math.isclose(float(graph.bare_prizes[i]), raw_prize, rel_tol=1e-12, abs_tol=1e-15):
                raise RuntimeError("Seed {} OI bare prize mismatch at {}".format(seed, internal_id))
            if not math.isclose(float(graph.prizes[i]), OI_PARAMS["b"] * raw_prize, rel_tol=1e-12, abs_tol=1e-14):
                raise RuntimeError("Seed {} OI scaled prize mismatch at {}".format(seed, internal_id))
        # Record OI's actual factorized vertex order, not the pre-sorted export order.
        csv_write(oi_nodes_path, [{"node_index": i, "internal_id": str(node),
                    "Entrez": str(node).split(":", 1)[1], "prize_unscaled": float(graph.bare_prizes[i]),
                    "model_evidence_status": node_gene_records.get(str(node).split(":", 1)[1], {}).get("model_evidence_status", "no_SHAP_evidence_PPI_connector")}
                    for i, node in enumerate(graph.nodes)],
                  ["node_index", "internal_id", "Entrez", "prize_unscaled", "model_evidence_status"])
        degree_array = graph.node_degrees.astype(np.float64)
        u, v = graph.edges[:, 0], graph.edges[:, 1]
        denom = ((len(graph.nodes) - degree_array[u] - 1) * (len(graph.nodes) - degree_array[v] - 1)
                 + degree_array[u] * degree_array[v])
        penalty = (10 ** OI_PARAMS["g"]) * (degree_array[u] * degree_array[v] / denom)
        if not np.allclose(graph.edge_penalties, penalty, rtol=1e-12, atol=1e-12):
            raise RuntimeError("Seed {} OI hub penalty formula mismatch".format(seed))
        if not np.allclose(graph.costs, 1.0 + penalty, rtol=1e-12, atol=1e-12):
            raise RuntimeError("Seed {} effective costs mismatch base+hub penalty".format(seed))
        vertices, edge_indices = graph.pcsf(pruning=PRUNING, verbosity_level=0)
        if any(int(i) < 0 or int(i) >= len(graph.nodes) for i in vertices):
            raise RuntimeError("Seed {} solver returned invalid/dummy vertex index".format(seed))
        if any(int(i) < 0 or int(i) >= len(edges) for i in edge_indices):
            raise RuntimeError("Seed {} solver returned invalid/dummy edge index".format(seed))
        selected_internal = {str(graph.nodes[int(i)]) for i in vertices}
        selected_ids = {value.split(":", 1)[1] for value in selected_internal if value.startswith("ENTREZ:")}
        if len(selected_ids) != len(selected_internal) or not selected_ids.issubset(ppi_nodes):
            raise RuntimeError("Seed {} solver output is not a subset of approved PPI nodes".format(seed))
        forest_rows = []
        for ix in edge_indices:
            edge = edges[int(ix)]
            forest_rows.append({"seed": seed, **edge, "hub_penalty": float(graph.edge_penalties[int(ix)]),
                                "effective_cost": float(graph.costs[int(ix)]), "in_solver_forest": True})
        forest_keys = {(r["Protein_A_Entrez_ID"], r["Protein_B_Entrez_ID"]) for r in forest_rows}
        induced_rows = []
        for edge in edges:
            a, b = edge["Protein_A_Entrez_ID"], edge["Protein_B_Entrez_ID"]
            if a in selected_ids and b in selected_ids:
                ix = ppi_edge_index[(a, b)]
                induced_rows.append({"seed": seed, **edge, "hub_penalty_OI_annotation_only": float(graph.edge_penalties[ix]),
                    "effective_cost_OI_annotation_only": float(graph.costs[ix]),
                    "in_solver_forest": (a, b) in forest_keys,
                    "edge_role": "solver_forest_edge" if (a, b) in forest_keys else "induced_only_edge"})
        node_rows = []
        for eid in sorted(selected_ids, key=sort_id):
            evidence = node_gene_records.get(eid)
            background = ppi_meta[eid]
            if evidence:
                record = dict(evidence)
                mean_abs = record["mean_abs_SHAP"]
                weight = record["model_weight_sum_normalized"]
                prize = record["prize"]
                symbols = record["model_gene_symbols"]
                ptgs = record["PTG_groups"]
                status = record["model_evidence_status"]
            else:
                mean_abs, weight, prize, symbols, ptgs, status = float("nan"), float("nan"), 0.0, "", "", "no_SHAP_evidence_PPI_connector"
            node_rows.append({"seed": seed, "Entrez": eid, "official_symbols": "|".join(sorted(official_by_eid.get(eid, set()))),
                "model_gene_symbols": symbols, "PTG_groups": ptgs, "mean_abs_SHAP": mean_abs,
                "model_weight_sum_normalized": weight, "prize": prize, "prize_source_symbol": evidence.get("prize_source_symbol", "") if evidence else "",
                "model_evidence_status": status, "background_degree": int(background["background_degree"]),
                "component_id": background["component_id"], "skeleton_degree": 0, "induced_ppi_degree": 0})
        node_map = {r["Entrez"]: r for r in node_rows}
        for row in forest_rows:
            node_map[row["Protein_A_Entrez_ID"]]["skeleton_degree"] += 1
            node_map[row["Protein_B_Entrez_ID"]]["skeleton_degree"] += 1
        for row in induced_rows:
            node_map[row["Protein_A_Entrez_ID"]]["induced_ppi_degree"] += 1
            node_map[row["Protein_B_Entrez_ID"]]["induced_ppi_degree"] += 1
        summary = {"seed": seed, "status": "complete", "n_model_features": EXPECTED_FEATURES,
            "n_ppi_background_nodes": len(ppi_nodes), "n_ppi_background_edges": len(edges),
            "n_positive_reward_entrez": len(positive_prizes), "n_selected_nodes": len(node_rows),
            "n_model_evidence_nodes": sum(r["model_evidence_status"] != "no_SHAP_evidence_PPI_connector" for r in node_rows),
            "n_connector_nodes": sum(r["model_evidence_status"] == "no_SHAP_evidence_PPI_connector" for r in node_rows),
            "n_true_zero_reward_nodes": sum(r["model_evidence_status"] == "measured_zero_reward" for r in node_rows),
            "n_skeleton_edges": len(forest_rows), "n_induced_ppi_edges": len(induced_rows),
            "n_ppi_components_in_background": len(set(components.values())), "elapsed_seconds": round(time.time()-started, 3)}
        call_audit = {"solver": "OmicsIntegrator2 Graph.pcsf -> installed pcst_fast",
            "OI_params": OI_PARAMS, "pruning": PRUNING, "verbosity_level": 0,
            "graph_nodes": len(graph.nodes), "graph_edges": len(graph.edges),
            "graph_edge_array_shape": list(graph.edges.shape), "graph_edge_array_dtype": str(graph.edges.dtype),
            "graph_cost_array_shape": list(graph.costs.shape), "graph_cost_array_dtype": str(graph.costs.dtype),
            "graph_prize_array_shape": list(graph.prizes.shape), "graph_prize_array_dtype": str(graph.prizes.dtype),
            "positive_terminal_count": len(graph.terminals), "virtual_root_index": len(graph.nodes),
            "virtual_root_mode": "terminals", "virtual_root_edge_count": len(graph.terminals),
            "virtual_root_edge_cost": float(OI_PARAMS["w"]), "num_clusters_inside_OI": 1,
            "dummy_root_removed_inside_OI_before_export": True,
            "selected_vertex_indices": [int(i) for i in vertices],
            "selected_real_edge_indices": [int(i) for i in edge_indices],
            "objective_terms": {"sum_unselected_scaled_prize": float(sum(float(graph.prizes[i]) for i in range(len(graph.nodes)) if i not in set(int(x) for x in vertices))),
                "sum_selected_real_edge_cost": float(sum(float(graph.costs[int(i)]) for i in edge_indices)),
                "virtual_root_edge_cost_sum": float(len(set(int(x) for x in edge_indices)) * 0.0),
                "component_count_after_root_removal": None,
                "OI_public_objective_includes_w_times_components": True},
            "input_edge_order_matches_interactome_dataframe": True,
            "base_cost_all_edges": 1.0,
            "finite_costs_and_prizes": bool(np.isfinite(graph.costs).all() and np.isfinite(graph.prizes).all()),
            "nonnegative_prizes_and_costs": bool(np.all(graph.prizes >= 0) and np.all(graph.costs > 0))}
        # Root-edge count equals the number of real forest components after removing the dummy root.
        forest_adjacency = defaultdict(set)
        for row in forest_rows:
            forest_adjacency[row["Protein_A_Entrez_ID"]].add(row["Protein_B_Entrez_ID"])
            forest_adjacency[row["Protein_B_Entrez_ID"]].add(row["Protein_A_Entrez_ID"])
        remaining_vertices = set(selected_ids)
        component_count = 0
        while remaining_vertices:
            component_count += 1
            stack = [remaining_vertices.pop()]
            while stack:
                current = stack.pop()
                neighbors = forest_adjacency[current] & remaining_vertices
                remaining_vertices.difference_update(neighbors)
                stack.extend(neighbors)
        call_audit["objective_terms"]["component_count_after_root_removal"] = component_count
        call_audit["objective_terms"]["virtual_root_edge_cost_sum"] = float(component_count * OI_PARAMS["w"])
        call_audit["objective_terms"]["objective_value"] = float(
            call_audit["objective_terms"]["sum_unselected_scaled_prize"] +
            call_audit["objective_terms"]["sum_selected_real_edge_cost"] +
            call_audit["objective_terms"]["virtual_root_edge_cost_sum"])
        node_fields = list(node_rows[0].keys()) if node_rows else ["seed", "Entrez", "model_gene_symbols", "model_evidence_status"]
        forest_fields = ["seed", "Protein_A_Entrez_ID", "Protein_B_Entrez_ID", "base_cost", "number_values", "data_source_union", "raw_source_texts", "raw_record_count", "raw_record_ids", "hub_penalty", "effective_cost", "in_solver_forest"]
        induced_fields = forest_fields[:9] + ["hub_penalty_OI_annotation_only", "effective_cost_OI_annotation_only", "in_solver_forest", "edge_role"]
        csv_write(os.path.join(run_dir, "gene_prizes.csv"), [dict(seed=seed, gene_id=r["gene_id"], PTG=r.get("PTG", ""), mean_abs_SHAP=float(r["mean_abs_SHAP"]), model_weight_sum_normalized=float(r["weight"]), prize_full_model_max_normalized=float(r["prize_full_model_max_normalized"])) for r in weights.to_dict("records")], ["seed", "gene_id", "PTG", "mean_abs_SHAP", "model_weight_sum_normalized", "prize_full_model_max_normalized"])
        csv_write(os.path.join(run_dir, "all_gene_mapping_audit.csv"), audit_rows,
                  ["seed", "gene_id", "PTG", "mean_abs_SHAP", "model_weight_sum_normalized", "prize_full_model_max_normalized", "mapping_status", "chosen_Entrez", "mapping_source", "candidate_Entrez_ids", "in_qualified_PPI", "PPI_background_degree", "PPI_component_id"])
        csv_write(os.path.join(run_dir, "selected_nodes.csv"), node_rows, node_fields)
        csv_write(os.path.join(run_dir, "selected_forest_edges.csv"), forest_rows, forest_fields)
        csv_write(os.path.join(run_dir, "selected_induced_ppi_edges.csv"), induced_rows, induced_fields)
        json_write(os.path.join(run_dir, "solver_call_audit.json"), call_audit)
        summary["status"] = "complete"
        outputs = [p for p in output_paths if os.path.basename(p) != "run_manifest.json"]
        outputs.append(os.path.join(run_dir, "solver_call_audit.json"))
        outputs.append(solver_edges_path); outputs.append(prize_path)
        output_hashes = {os.path.basename(p): sha256(p) for p in outputs}
        manifest = {"status": "complete", "seed": seed, "contract": EXPECTED_CONTRACT,
            "input_hashes": {"weights.csv": seed_info[seed]["weights_sha256"], "candidates.csv": seed_info[seed]["candidates_sha256"],
                "model.pt": seed_info[seed]["model_sha256"], "mentor_PPI": sha256(PPI_PATH), "official_map": sha256(OFFICIAL_MAP), "alias_map": sha256(ALIAS_MAP),
                "model_contract_json": sha256(os.path.join(MODEL_ROOT, "data", "contract.json")),
                "model_genes_csv": sha256(os.path.join(MODEL_ROOT, "data", "genes.csv")),
                "model_cells_csv": sha256(os.path.join(MODEL_ROOT, "data", "cells.csv"))},
            "runner_script_path": SCRIPT_PATH, "runner_script_sha256": sha256(SCRIPT_PATH),
            "input_feature_count": EXPECTED_FEATURES, "input_feature_order_sha256": hashlib.sha256("\n".join(genes_order).encode("utf-8")).hexdigest(),
            "input_transform": "RNA data layer log1p normalized expression; frozen PTG14 seed model input",
            "split": {str(k): int(v) for k, v in cell_manifest.groupby(["split", "Group3"]).size().items()},
            "model_seed": seed, "pcsf_seed": 42, "OI_params": OI_PARAMS, "pruning": PRUNING,
            "hub_penalty_formula": config["pcsf"]["hub_penalty_formula"], "summary": summary,
            "output_hashes": output_hashes, "started_at_local": datetime.fromtimestamp(started).astimezone().isoformat(timespec="seconds"),
            "finished_at_local": datetime.now().astimezone().isoformat(timespec="seconds")}
        # Core required outputs are hashed into manifest; manifest itself is the final atomic completion marker.
        json_write(existing_manifest_path, manifest)
        run_summaries.append(summary)
        symbols = {s for r in node_rows for s in (r.get("model_gene_symbols", "") or "").split(";") if s}
        all_seed_selected_symbols.extend([{"gene_id": sym, "seed": seed} for sym in sorted(symbols)])
        all_seed_selected_nodes.extend([{"Entrez": r["Entrez"], "seed": seed,
            "model_gene_symbols": r.get("model_gene_symbols", ""),
            "model_evidence_status": r.get("model_evidence_status", "no_SHAP_evidence_PPI_connector")} for r in node_rows])
        log("seed {} complete: selected nodes={}, model evidence={}, connectors={}, forest edges={}, induced edges={}".format(
            seed, summary["n_selected_nodes"], summary["n_model_evidence_nodes"], summary["n_connector_nodes"], summary["n_skeleton_edges"], summary["n_induced_ppi_edges"]))
        if seed == 42:
            # Technical gate: verify the first run has all required outputs and identity before any further seed.
            for path in outputs:
                if not os.path.isfile(path) or sha256(path) != output_hashes[os.path.basename(path)]:
                    raise RuntimeError("Seed42 technical gate output/hash check failed: {}".format(path))
            log("seed42 technical gate passed; continuing serially with frozen seeds 43-51")

    # Recurrence is counted over actual selected Entrez nodes, including model-external connectors.
    selected_nodes_df = pd.DataFrame(all_seed_selected_nodes)
    if not selected_nodes_df.empty and selected_nodes_df.duplicated(["Entrez", "seed"]).any():
        selected_nodes_df = selected_nodes_df.drop_duplicates(["Entrez", "seed"])
    feature_ptg = genes.drop_duplicates("gene_id").set_index("gene_id").PTG.to_dict()
    frequency = selected_nodes_df.groupby("Entrez").seed.nunique().rename("selected_runs").reset_index() if len(selected_nodes_df) else pd.DataFrame.from_records([], columns=["Entrez", "selected_runs"])
    frequency["selected_runs"] = frequency.selected_runs.astype(int)
    frequency["selection_frequency"] = frequency.selected_runs / 10.0
    frequency["core_10_of_10"] = frequency.selected_runs == 10
    node_symbols = {eid: sorted({s for value in selected_nodes_df.loc[selected_nodes_df.Entrez.astype(str).eq(str(eid)), "model_gene_symbols"].fillna("") for s in str(value).split(";") if s}) for eid in frequency.Entrez.astype(str)}
    frequency["model_gene_symbols"] = frequency.Entrez.astype(str).map(lambda eid: ";".join(node_symbols.get(eid, [])))
    frequency["model_evaluable"] = frequency.model_gene_symbols.ne("")
    frequency["connector_node"] = ~frequency.model_evaluable
    frequency["PTG"] = frequency.model_gene_symbols.map(lambda text: ";".join(sorted({feature_ptg[s] for s in text.split(";") if s in feature_ptg})))
    frequency = frequency.sort_values(["selected_runs", "Entrez"], ascending=[False, True], kind="mergesort")
    frequency["selection_rank"] = np.arange(1, len(frequency)+1)
    majority = frequency[frequency.selected_runs >= 6].copy()
    model_frequency = []
    for eid in majority.Entrez.astype(str):
        ent_row = majority[majority.Entrez.astype(str).eq(eid)].iloc[0]
        for symbol in node_symbols.get(eid, []):
            model_frequency.append({"gene_id": symbol, "Entrez": eid, "PTG": feature_ptg.get(symbol, ""),
                "selected_runs": int(ent_row.selected_runs), "selection_frequency": float(ent_row.selection_frequency),
                "core_10_of_10": bool(ent_row.core_10_of_10), "model_evaluable": True})
    model_majority = pd.DataFrame(model_frequency)
    if len(model_majority):
        model_majority["median_weight"] = model_majority.gene_id.map({g: float(np.median([float(seed_info[s]["weights"].loc[seed_info[s]["weights"].gene_id.astype(str).eq(g), "weight"].iloc[0]) for s in SEEDS])) for g in model_majority.gene_id})
        model_majority = model_majority.sort_values(["selected_runs", "gene_id"], ascending=[False, True], kind="mergesort")
        model_majority["selected_rank"] = np.arange(1, len(model_majority)+1)
    else:
        model_majority = pd.DataFrame.from_records([], columns=["gene_id", "Entrez", "PTG", "selected_runs", "selection_frequency", "core_10_of_10", "model_evaluable", "median_weight", "selected_rank"])
    freq_fields = ["Entrez", "model_gene_symbols", "PTG", "selected_runs", "selection_frequency", "core_10_of_10", "model_evaluable", "connector_node", "selection_rank"]
    csv_write(os.path.join(OUT_ROOT, "source_data", "all_seed_selected_node_frequency.csv"), frequency.to_dict("records"), freq_fields)
    csv_write(os.path.join(OUT_ROOT, "source_data", "majority_selected_nodes.csv"), majority.to_dict("records"), freq_fields)
    csv_write(os.path.join(OUT_ROOT, "source_data", "majority_selected_model_genes.csv"), model_majority.to_dict("records"),
              ["gene_id", "Entrez", "PTG", "selected_runs", "selection_frequency", "core_10_of_10", "median_weight", "model_evaluable", "selected_rank"])
    long_weights = []
    selected_by_symbol_seed = {}
    for row in model_majority.to_dict("records"):
        selected_by_symbol_seed[row["gene_id"]] = {seed for seed in SEEDS
            if not selected_nodes_df[(selected_nodes_df.seed.astype(int) == seed) & (selected_nodes_df.Entrez.astype(str) == str(row["Entrez"]))].empty}
    for row in model_majority.to_dict("records"):
        for seed in SEEDS:
            wrow = seed_info[seed]["weights"].loc[seed_info[seed]["weights"].gene_id.astype(str).eq(row["gene_id"])].iloc[0]
            is_selected = seed in selected_by_symbol_seed[row["gene_id"]]
            long_weights.append({"gene_id": row["gene_id"], "PTG": row["PTG"], "selected_runs": int(row["selected_runs"]),
                "core_10_of_10": bool(row["core_10_of_10"]), "seed": seed,
                "mean_abs_SHAP": float(wrow["mean_abs_SHAP"]), "weight": float(wrow["weight"]),
                "selected_this_seed": bool(is_selected)})
    csv_write(os.path.join(OUT_ROOT, "source_data", "majority_selected_weights_long.csv"), long_weights,
              ["gene_id", "PTG", "selected_runs", "core_10_of_10", "seed", "mean_abs_SHAP", "weight", "selected_this_seed"])

    # C uses only nodes selected in at least 6/10 PCSF runs, including model-external connectors.
    # No shortest-path filler nodes or edges are added.
    display_ids = set(majority.Entrez.astype(str))
    node_output = []
    median_weights = model_majority.set_index("gene_id").median_weight.to_dict() if len(model_majority) else {}
    for eid in sorted(display_ids, key=sort_id):
        symbols = node_symbols.get(eid, [])
        model_gene = bool(symbols)
        p = ppi_meta[eid]
        frequency_row = majority[majority.Entrez.astype(str).eq(eid)].iloc[0]
        node_seed_weights = []
        for seed in SEEDS:
            seed_w = seed_info[seed]["weights"]
            seed_rows = seed_w[seed_w.gene_id.astype(str).isin(symbols)]
            if not seed_rows.empty:
                node_seed_weights.append(float(seed_rows.weight.max()))
        node_output.append({"Entrez": eid, "gene_id": ";".join(symbols), "official_symbols": "|".join(sorted(official_by_eid.get(eid, set()))),
            "PTG": ";".join(sorted({feature_ptg.get(s, "") for s in symbols if feature_ptg.get(s)})),
            "selected_runs": int(frequency_row.selected_runs),
            "core_10_of_10": bool(frequency_row.core_10_of_10),
            "median_weight": float(np.median(node_seed_weights)) if node_seed_weights else float("nan"),
            "model_evaluable": model_gene, "node_role": "majority_model_gene" if model_gene else "model_outside_connector",
            "background_degree": int(p["background_degree"]), "component_id": p["component_id"]})
    node_set = {r["Entrez"] for r in node_output}
    edge_frequency = defaultdict(int)
    for seed in SEEDS:
        run_edges = pandas_read_csv(os.path.join(OUT_ROOT, "networks", "seed_{}".format(seed), "selected_forest_edges.csv"), dtype=str)
        for row in run_edges.to_dict("records"):
            edge_frequency[(row["Protein_A_Entrez_ID"], row["Protein_B_Entrez_ID"])] += 1
    display_edges = []
    for edge in edges:
        a, b = edge["Protein_A_Entrez_ID"], edge["Protein_B_Entrez_ID"]
        if a in node_set and b in node_set:
            count = edge_frequency.get((a, b), 0)
            display_edges.append({"Protein_A_Entrez_ID": a, "Protein_B_Entrez_ID": b, **edge,
                                  "selected_in_seed_count": count,
                                  "display_category": "selected_in_at_least_1_of_10_runs" if count >= 1 else "PPI_only_not_selected_in_any_run"})
    csv_write(os.path.join(OUT_ROOT, "source_data", "Fig5C_PPI_PCSF_candidate_nodes.csv"), node_output,
              ["Entrez", "gene_id", "official_symbols", "PTG", "selected_runs", "core_10_of_10", "median_weight", "model_evaluable", "node_role", "background_degree", "component_id"])
    csv_write(os.path.join(OUT_ROOT, "source_data", "Fig5C_PPI_PCSF_candidate_edges.csv"), display_edges,
              ["Protein_A_Entrez_ID", "Protein_B_Entrez_ID", "base_cost", "number_values", "data_source_union", "raw_source_texts", "raw_record_count", "raw_record_ids", "selected_in_seed_count", "display_category"])
    csv_write(os.path.join(OUT_ROOT, "source_data", "all_gene_mapping_audit.csv"), mapping_audit_all,
              ["seed", "gene_id", "PTG", "mean_abs_SHAP", "model_weight_sum_normalized", "prize_full_model_max_normalized", "mapping_status", "chosen_Entrez", "mapping_source", "candidate_Entrez_ids", "in_qualified_PPI", "PPI_background_degree", "PPI_component_id"])
    all_weights = []
    for seed in SEEDS:
        for row in seed_info[seed]["weights"].to_dict("records"):
            all_weights.append({"seed": seed, "gene_id": row["gene_id"], "PTG": row["PTG"],
                "mean_abs_SHAP": float(row["mean_abs_SHAP"]), "weight": float(row["weight"])})
    csv_write(os.path.join(OUT_ROOT, "source_data", "all_seed_full_gene_weights.csv"), all_weights,
              ["seed", "gene_id", "PTG", "mean_abs_SHAP", "weight"])
    csv_write(os.path.join(OUT_ROOT, "source_data", "run_factsheet.csv"), run_summaries,
              list(run_summaries[0].keys()))
    config["status"] = "complete"
    config["cross_seed_results"] = {"selection_frequency_union_size": len(frequency), "majority_count_ge6": len(majority),
                                     "model_evaluable_majority_count": int(majority.model_evaluable.sum()),
                                     "core_10_of_10_count": int(majority.core_10_of_10.sum()),
                                     "majority_gene_connector_nodes_for_C": sum(r["node_role"] == "model_outside_connector" for r in node_output),
                                     "majority_gene_network_nodes_for_C": sum(r["node_role"] == "majority_model_gene" for r in node_output)}
    json_write(os.path.join(OUT_ROOT, "config", "run_config.json"), config)
    status = {"status": "complete", "seed_runs": run_summaries, "cross_seed_results": config["cross_seed_results"],
              "completed_at": datetime.now().astimezone().isoformat(timespec="seconds")}
    json_write(os.path.join(OUT_ROOT, "logs", "run_status.json"), status)
    log("COMPLETE: majority>=6={}, model-evaluable={}, core10/10={}".format(
        len(majority), int(majority.model_evaluable.sum()), int(majority.core_10_of_10.sum())))


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = {"status": "failed", "error_type": type(exc).__name__, "error": str(exc),
                   "traceback": traceback.format_exc(), "failed_at": datetime.now().astimezone().isoformat(timespec="seconds")}
        try:
            json_write(os.path.join(OUT_ROOT, "logs", "failure.json"), failure)
            log("FAILED: {}: {}".format(type(exc).__name__, exc))
        except Exception:
            pass
        raise
