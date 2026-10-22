#!/usr/bin/env python3
"""Per-disease compact Fig6 PCSF parameter selection and independent stability.

The solver runtime is copied from the accepted size-target scan and locally
patched only to use the already-corrected full effective-edge-cost perturbation
and to resolve comparison networks from the frozen base scan. This driver never
writes to formal networks and never retrains models or recomputes SHAP.
"""
from __future__ import print_function

import argparse
import csv
import hashlib
import importlib.util
import itertools
import json
import math
import os
import shutil
import sys
import time
import traceback

import numpy as np
import pandas as pd


ANALYSIS_ROOT = r"E:\micro_output\分析审阅中心\待审阅\Fig6_疾病vsCtrl_SHAP_PCSF模块"
RUN_ID = "20261021_disease_compact_le50_parameter_selection_01"
SCAN_ID = "20261021_disease_compact_le50_parameter_selection_01_base_scan"
RUN_ROOT = os.path.join(ANALYSIS_ROOT, "robustness", RUN_ID)
SCAN_ROOT = os.path.join(RUN_ROOT, "base_scan")
DISEASES = ["AD", "ALS", "ALS_FTD", "DLB", "DLBD_PD", "EOAD", "FTD", "FTD_MND",
            "GBM", "HD", "LCBM", "LGG", "MB", "MS", "PD", "PSP"]
B_VALUES = [0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
W_VALUES = [3, 6]
G_VALUES = [2, 3, 4, 5]
MASTER_SEED = 42
STAGE_SEED_LABELS = {"selection": "compact_le50_selection_v1",
                     "evaluation": "compact_le50_evaluation_v1",
                     "degree_null": "compact_le50_degree_stratified_shuffle_v1"}
SIGMA = 0.1
SELECTION_REPEATS = 20
EVALUATION_REPEATS = 100
DEGREE_NULL_REPEATS = 100
MIN_NODES = 30
MAX_NODES = 50
BASE_RUNTIME_SOURCE = r"D:\data\microglia_new\scripts\fig6_disease_pcsf_size_target_scan.py"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bytes_sha256(value):
    return hashlib.sha256(value).hexdigest()


def float_or_nan(value):
    """Convert nullable scalar metrics from frozen CSV/JSON audits safely."""
    if value is None or pd.isna(value):
        return float("nan")
    return float(value)


def write_json(path, value):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True,
                  allow_nan=False)
        handle.write("\n")
    os.replace(tmp, path)


def log(message):
    line = "[{}] {}".format(time.strftime("%Y-%m-%dT%H:%M:%S%z"), message)
    path = os.path.join(RUN_ROOT, "logs", "execution.log")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a", encoding="utf-8") as handle:
        handle.write(line + "\n")
    print(line, flush=True)


def atomic_csv(path, frame):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    # pandas 0.23 in the frozen Fig6 Python runtime uses the legacy parameter name.
    frame.to_csv(tmp, index=False, encoding="utf-8-sig", line_terminator="\n")
    os.replace(tmp, path)


def read_csv_unicode(path, **kwargs):
    """Read CSV through an open handle for the frozen pandas/Windows runtime."""
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        return pd.read_csv(handle, **kwargs)


def runtime_copy_and_patch():
    if not os.path.isfile(BASE_RUNTIME_SOURCE):
        raise FileNotFoundError(BASE_RUNTIME_SOURCE)
    runtime_dir = os.path.join(RUN_ROOT, "scripts")
    os.makedirs(runtime_dir, exist_ok=True)
    runtime_path = os.path.join(runtime_dir, "pcsf_runtime_patched.py")
    source_hash = sha256(BASE_RUNTIME_SOURCE)
    if not os.path.isfile(runtime_path):
        with open(BASE_RUNTIME_SOURCE, "r", encoding="utf-8-sig") as handle:
            source = handle.read()
        replacements = [
            ('new_costs = gctx["base_costs"] * multiplier + gctx["edge_penalties"]',
             'frozen_effective = np.asarray(gctx["frozen_effective_costs"], dtype=np.float64)\n'
             '        independently_rebuilt = np.asarray(gctx["base_costs"], dtype=np.float64) + np.asarray(gctx["edge_penalties"], dtype=np.float64)\n'
             '        if not np.array_equal(frozen_effective, independently_rebuilt):\n'
             '            raise RuntimeError("Frozen full effective costs differ from base plus hub penalty")\n'
             '        new_costs = frozen_effective * multiplier'),
            ('"cost_noise_seed_replay": ("RandomState(seed).normal(0,0.1,n_edges), multiplier=exp(z-0.1^2/2), applied to base cost before fixed hub penalty" if control_type == "edge_cost" else None),',
             '"cost_noise_seed_replay": ("RandomState(seed).normal(0,sigma,n_real_edges) in frozen PPI edge order; exp(z-sigma^2/2) multiplies full effective cost (base+hub); virtual-root edges remain w" if control_type == "edge_cost" else None),'),
            ('ref_path = os.path.join(RUN, "grid", disease, pid, "adopted_nodes.csv")',
             'ref_path = frozen_base_result_path(disease, pid, "adopted_nodes.csv")'),
            ('ref_edge_path = os.path.join(RUN, "grid", disease, pid, "adopted_forest_edges.tsv")',
             'ref_edge_path = frozen_base_result_path(disease, pid, "adopted_forest_edges.tsv")'),
            ('record["record_payload_sha256"]=_record_digest(record)\n    write_json(result_path, _json_clean(record))\n    return record',
             'record["record_payload_sha256"]=_record_digest(record)\n    return record'),
        ]
        for old, new in replacements:
            count = source.count(old)
            if count != 1:
                raise RuntimeError("Runtime patch expected exactly one occurrence; got {} for {}".format(count, old[:80]))
            source = source.replace(old, new)
        # Verify solver-call costs before objective auditing and persistence.
        call_old = 'raw = call_solver(graph)\n    solver_runtime_seconds = time.time() - solver_started'
        call_new = ('raw = call_solver(graph)\n'
                    '    n_real_edges = len(graph.edges)\n'
                    '    actual_costs = np.asarray(raw["solver_costs"], dtype=np.float64)\n'
                    '    expected_roots = np.full(len(graph.terminals), float(w), dtype=np.float64)\n'
                    '    if not np.array_equal(actual_costs[:n_real_edges], np.asarray(new_costs, dtype=np.float64)):\n'
                    '        raise RuntimeError("Solver real-edge costs differ from the requested effective costs")\n'
                    '    if not np.array_equal(actual_costs[n_real_edges:], expected_roots):\n'
                    '        raise RuntimeError("Solver virtual-root costs differ from frozen w")\n'
                    '    solver_runtime_seconds = time.time() - solver_started')
        if source.count(call_old) != 1:
            raise RuntimeError("Could not uniquely add solver cost-array validation")
        source = source.replace(call_old, call_new)
        with open(runtime_path, "w", encoding="utf-8", newline="") as handle:
            handle.write(source)
    runtime_hash = sha256(runtime_path)
    return runtime_path, source_hash, runtime_hash


def load_runtime(runtime_path):
    name = "fig6_stability_solver_runtime"
    spec = importlib.util.spec_from_file_location(name, runtime_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parameter_id(mod, b, w, g):
    return mod.param_id(float(b), int(w), int(g))


def expected_parameter_keys():
    return set((d, float(b), int(w), int(g))
               for d in DISEASES for b, w, g in itertools.product(B_VALUES, W_VALUES, G_VALUES))


def verify_scan_summary(mod):
    metrics_path = os.path.join(SCAN_ROOT, "summaries", "per_disease_parameter_network_metrics.csv")
    status_path = os.path.join(SCAN_ROOT, "summaries", "base_run_status.csv")
    if not os.path.isfile(metrics_path) or not os.path.isfile(status_path):
        raise FileNotFoundError("The hash-assembled 104-point base-grid tables are missing")
    frame = mod.read_table(metrics_path)
    status = mod.read_table(status_path)
    expected = set((d, float(b), int(w), int(g))
                   for d in DISEASES for b, w, g in itertools.product(B_VALUES, W_VALUES, G_VALUES))
    observed=[]
    for _, row in frame.iterrows():
        observed.append((str(row["Disease"]), float(row["b"]), int(row["w"]), int(row["g"])))
    if len(frame) != 1664 or len(set(observed)) != 1664 or set(observed) != expected:
        raise RuntimeError("Assembled base grid is not the unique frozen 16x104 parameter set")
    if len(status) != 1664 or status.duplicated(["Disease", "parameter_id"]).any():
        raise RuntimeError("Assembled base status registry is incomplete or duplicated")
    joined=frame[["Disease","parameter_id","result_origin","base_run_manifest_path","base_run_manifest_sha256"]].merge(
        status[["Disease","parameter_id","result_origin","base_run_manifest_path","base_run_manifest_sha256","status"]],
        on=["Disease","parameter_id"],how="outer",validate="one_to_one",suffixes=("_metrics","_status"),indicator=True)
    if len(joined)!=1664 or not joined["_merge"].eq("both").all():
        raise RuntimeError("Assembled base metrics and status keys differ")
    for field in ["result_origin","base_run_manifest_path","base_run_manifest_sha256"]:
        if not joined[field+"_metrics"].astype(str).equals(joined[field+"_status"].astype(str)):
            raise RuntimeError("Assembled base {} values differ between metrics and status registry".format(field))
    allowed={"complete","complete_reused_prior_base"}
    if not set(status["status"]).issubset(allowed):
        raise RuntimeError("Assembled grid has non-complete source rows: {}".format(sorted(set(status["status"])-allowed)))
    return frame


def qualify_base_networks(mod, metrics, ppi_df, disease_data):
    current_input_manifest = os.path.join(RUN_ROOT, "inputs", "input_sha256_manifest.csv")
    if not os.path.isfile(current_input_manifest):
        raise FileNotFoundError(current_input_manifest)
    current_input_hash = sha256(current_input_manifest)
    rows, candidate_files = [], {}
    for _, source in metrics.iterrows():
        disease=str(source["Disease"]); b=float(source["b"]); w=int(source["w"]); g=int(source["g"])
        pid=str(source["parameter_id"])
        expected_pid=parameter_id(mod,b,w,g)
        if pid != expected_pid:
            raise RuntimeError("Fractional-b-safe parameter id mismatch: {} vs {}".format(pid,expected_pid))
        manifest_path=str(source["base_run_manifest_path"])
        config_path=str(source["source_run_config_path"])
        if not os.path.isfile(manifest_path) or not os.path.isfile(config_path):
            raise FileNotFoundError("Source provenance file missing for {} {}".format(disease,pid))
        if sha256(manifest_path)!=str(source["base_run_manifest_sha256"]):
            raise RuntimeError("Source base manifest SHA mismatch: {}".format(manifest_path))
        config=mod.read_json(config_path) if hasattr(mod,"read_json") else json.load(open(config_path,"r",encoding="utf-8-sig"))
        config_sha=sha256(config_path)
        with open(manifest_path,"r",encoding="utf-8-sig") as handle: manifest=json.load(handle)
        if (manifest.get("status")!="complete" or manifest.get("Disease")!=disease or
                manifest.get("parameter_id")!=pid or not math.isclose(float(manifest.get("b")),b,rel_tol=0.0,abs_tol=1e-12) or
                int(manifest.get("w"))!=w or int(manifest.get("g"))!=g):
            raise RuntimeError("Source manifest identity/status mismatch: {}".format(manifest_path))
        if (manifest.get("run_config_sha256")!=config_sha or
                config_sha!=str(source["source_run_config_sha256"]) or
                manifest.get("execution_script_sha256")!=config.get("execution_script_sha256") or
                manifest.get("input_manifest_sha256")!=current_input_hash or
                manifest.get("input_manifest_sha256")!=str(source["input_manifest_sha256"])):
            raise RuntimeError("Source input/config/script identity mismatch: {}".format(manifest_path))
        impl=str(manifest.get("implementation",""))
        if "OmicsIntegrator2 2.3.10" not in impl or "pcst-fast 1.0.7" not in impl:
            raise RuntimeError("Unexpected solver identity in {}: {}".format(manifest_path,impl))
        if manifest.get("pruning")!="strong" or manifest.get("dummy_mode")!="terminals":
            raise RuntimeError("Frozen solver call settings differ: {}".format(manifest_path))
        output_by_name={os.path.basename(item["path"]):item for item in manifest.get("outputs",[])}
        required=["adopted_nodes.csv","adopted_forest_edges.tsv","solver_arrays.npz",
                  "candidate_objective_comparison.csv","raw_solver_vertices.csv","raw_solver_edges.csv"]
        if not set(required).issubset(output_by_name):
            raise RuntimeError("Required base outputs missing from manifest: {}".format(manifest_path))
        for name,item in output_by_name.items():
            path=str(item["path"])
            if not os.path.isfile(path) or sha256(path)!=item.get("sha256") or os.path.getsize(path)!=int(item.get("size_bytes",-1)):
                raise RuntimeError("Base output SHA/size mismatch: {}".format(path))
        node_path=output_by_name["adopted_nodes.csv"]["path"]
        edge_path=output_by_name["adopted_forest_edges.tsv"]["path"]
        nodes=mod.read_table(node_path)
        edges=mod.read_table(edge_path,"\t")
        if "Entrez" not in nodes.columns or not {"Entrez_A","Entrez_B"}.issubset(edges.columns):
            raise RuntimeError("Frozen selected node/edge schema mismatch: {}".format(manifest_path))
        node_ids=nodes["Entrez"].astype(str).tolist()
        if len(node_ids)!=len(set(node_ids)):
            raise RuntimeError("Duplicate Entrez in selected nodes: {}".format(node_path))
        edge_keys=["{}|{}".format(*sorted((str(a),str(b)))) for a,b in zip(edges["Entrez_A"],edges["Entrez_B"])]
        if len(edge_keys)!=len(set(edge_keys)) or any(a==b for a,b in [key.split("|") for key in edge_keys]):
            raise RuntimeError("Duplicate undirected edge or self-loop in selected backbone: {}".format(edge_path))
        node_set=set(node_ids)
        if any(a not in node_set or b not in node_set for a,b in [key.split("|") for key in edge_keys]):
            raise RuntimeError("Selected backbone endpoint is absent from selected nodes: {}".format(edge_path))
        graph=mod.nx.Graph();graph.add_nodes_from(node_ids);graph.add_edges_from([tuple(k.split("|")) for k in edge_keys])
        n=len(node_ids);m=len(edge_keys);components=mod.nx.number_connected_components(graph) if n else 0
        is_tree=bool(n>0 and components==1 and m==n-1 and mod.nx.is_tree(graph))
        recorded=manifest.get("metrics",{})
        if (int(recorded.get("real_node_count",-1))!=n or int(recorded.get("real_skeleton_edge_count",-1))!=m or
                int(recorded.get("component_count",-1))!=components or int(float(source["real_node_count"]))!=n or
                int(float(source["real_skeleton_edge_count"]))!=m or int(float(source["component_count"]))!=components):
            raise RuntimeError("Base manifest/summary/node-edge counts disagree: {}".format(manifest_path))
        adopted=str(recorded.get("adopted_solution_source",""))
        witnesses=[x for x in manifest.get("comparison_candidates",[]) if x.get("solution_type")==adopted]
        if not witnesses or not any(x.get("feasible") is True for x in witnesses):
            raise RuntimeError("Adopted solution lacks a feasible witness: {}".format(manifest_path))
        size_ok=MIN_NODES<=n<=MAX_NODES
        eligible=bool(size_ok and is_tree and manifest.get("status")=="complete")
        reasons=[]
        if not size_ok: reasons.append("node_count_outside_30_50")
        if not is_tree: reasons.append("not_one_connected_true_backbone_tree")
        row={"Disease":disease,"b":b,"w":w,"g":g,"parameter_id":pid,
             "result_origin":str(source["result_origin"]),"source_run_id":str(source["source_run_id"]),
             "real_node_count":n,"model_evidence_node_count":int(recorded.get("model_evidence_node_count",-1)),
             "measured_zero_SHAP_evidence_count":int(recorded.get("measured_zero_SHAP_evidence_count",-1)),
             "zero_reward_no_evidence_connector_count":int(recorded.get("zero_reward_no_evidence_connector_count",-1)),
             "real_skeleton_edge_count":m,"component_count":components,"isolated_node_count":int(recorded.get("isolated_node_count",-1)),
             "largest_component_fraction":float_or_nan(recorded.get("largest_component_fraction")),
             "retained_SHAP_mass_fraction_full_model":float_or_nan(recorded.get("retained_SHAP_mass_fraction_full_model")),
             "retained_SHAP_mass_fraction_mapped_PPI":float_or_nan(recorded.get("retained_SHAP_mass_fraction_mapped_PPI")),
             "adopted_solution_source":adopted,"adopted_objective":float_or_nan(recorded.get("adopted_objective")),
             "result_complete":True,"adopted_solution_feasible":True,"traceable_node_edge_tables_hash_verified":True,
             "one_component_tree_backbone":is_tree,"within_30_50_nodes":size_ok,"eligible":eligible,
             "rejection_reason":"|".join(reasons),"base_run_manifest_path":manifest_path,
             "base_run_manifest_sha256":sha256(manifest_path),"source_run_config_path":config_path,
             "source_run_config_sha256":config_sha,"source_execution_script_sha256":manifest.get("execution_script_sha256"),
             "selected_nodes_path":node_path,"selected_nodes_sha256":output_by_name["adopted_nodes.csv"]["sha256"],
             "selected_edges_path":edge_path,"selected_edges_sha256":output_by_name["adopted_forest_edges.tsv"]["sha256"],
             "solver_arrays_path":output_by_name["solver_arrays.npz"]["path"],"solver_arrays_sha256":output_by_name["solver_arrays.npz"]["sha256"]}
        rows.append(row);candidate_files[(disease,b,w,g)]=row
    out=pd.DataFrame(rows).sort_values(["Disease","b","w","g"]).reset_index(drop=True)
    if len(out)!=1664 or out.duplicated(["Disease","b","w","g"]).any():
        raise RuntimeError("Eligibility audit is not the unique 16x104 complete base grid")
    atomic_csv(os.path.join(RUN_ROOT,"summaries","parameter_qualification_all_104.csv"),out)
    no_eligible=[]
    for disease,group in out.groupby("Disease",sort=False):
        eg=group[group["eligible"]]
        if eg.empty:
            counts=group["rejection_reason"].value_counts().to_dict()
            no_eligible.append({"Disease":disease,"status":"no_eligible_candidate","candidate_count":0,"grid_rows":len(group),
                                "rejection_reason_counts":json.dumps(counts,ensure_ascii=False,sort_keys=True)})
    atomic_csv(os.path.join(RUN_ROOT,"summaries","no_eligible_candidates.csv"),pd.DataFrame(no_eligible,
        columns=["Disease","status","candidate_count","grid_rows","rejection_reason_counts"]))
    return out,candidate_files


def audit_frozen_solver_arrays(mod, ppi_df, disease_data, eligibility):
    rows = []
    contexts = {}
    for disease in DISEASES:
        for g in G_VALUES:
            contexts[(disease, g)] = mod.build_graph(disease, g, ppi_df, disease_data)
    for _, row in eligibility.iterrows():
        disease = str(row["Disease"])
        b, w, g = float(row["b"]), int(row["w"]), int(row["g"])
        pid = str(row["parameter_id"])
        path = str(row["solver_arrays_path"])
        if not os.path.isfile(path) or sha256(path) != str(row["solver_arrays_sha256"]):
            raise RuntimeError("Frozen solver-array file/hash mismatch: {}".format(path))
        ctx = contexts[(disease, g)]
        graph = ctx["graph"]
        expected_prizes = np.asarray(ctx["bare"], dtype=np.float64) * b
        expected_terminals = np.flatnonzero(expected_prizes > 0).astype(np.int64)
        expected_real_edges = np.asarray(graph.edges, dtype=np.int64)
        expected_real_costs = np.asarray(graph.costs, dtype=np.float64)
        root = len(ctx["nodes"])
        expected_solver_edges = np.vstack([
            expected_real_edges,
            np.column_stack([np.full(len(expected_terminals), root, dtype=np.int64), expected_terminals])
        ])
        expected_solver_prizes = np.concatenate([expected_prizes, np.asarray([0.0], dtype=np.float64)])
        expected_solver_costs = np.concatenate([expected_real_costs, np.full(len(expected_terminals), float(w))])
        with np.load(path, allow_pickle=False) as saved:
            checks = {
                "edge_order": np.array_equal(saved["edge_uv"], expected_real_edges),
                "prize_scale_applied_once": np.array_equal(saved["prizes"], expected_prizes),
                "effective_costs": np.allclose(saved["real_costs"], expected_real_costs, rtol=1e-12, atol=1e-12),
                "root_terminals": np.array_equal(saved["root_terminals"], expected_terminals),
                "augmented_edges": np.array_equal(saved["solver_edges"], expected_solver_edges),
                "augmented_prizes": np.array_equal(saved["solver_prizes"], expected_solver_prizes),
                "augmented_costs": np.array_equal(saved["solver_costs"], expected_solver_costs),
                "root_index": int(np.asarray(saved["solver_root"]).reshape(-1)[0]) == root,
                "num_clusters": int(np.asarray(saved["solver_num_clusters"]).reshape(-1)[0]) == 1,
            }
        if not all(checks.values()):
            raise RuntimeError("Frozen solver-array reconstruction mismatch for {} {}: {}".format(
                disease, pid, checks))
        rows.append({"Disease": disease, "parameter_id": pid, "b": b, "w": w, "g": g,
                     "solver_arrays_path": path, "solver_arrays_sha256": sha256(path),
                     "status": "verified", "checks": "|".join(checks.keys()),
                     "effective_cost_formula_verified": True, "prize_b_applied_once_verified": True})
    out = pd.DataFrame(rows).sort_values(["Disease", "b", "w", "g"]).reset_index(drop=True)
    if len(out) != 1664 or not out["status"].eq("verified").all():
        raise RuntimeError("Frozen base solver-array audit did not cover all 1664 results")
    atomic_csv(os.path.join(RUN_ROOT, "summaries", "frozen_solver_array_audit.csv"), out)
    return out


def make_source_manifest(mod, runtime_source_hash, runtime_hash, driver_hash):
    paths = [
        os.path.join(SCAN_ROOT, "config", "run_config.json"),
        os.path.join(SCAN_ROOT, "summaries", "size_target_summary_manifest.json"),
        os.path.join(SCAN_ROOT, "summaries", "per_disease_parameter_network_metrics.csv"),
        os.path.join(SCAN_ROOT, "summaries", "base_run_status.csv"),
        os.path.join(SCAN_ROOT, "inputs", "input_sha256_manifest.csv"),
        os.path.join(ANALYSIS_ROOT, "config", "run_config.json"),
        os.path.join(ANALYSIS_ROOT, "tables", "ppi_edges_ge2.tsv"),
        os.path.join(ANALYSIS_ROOT, "tables", "ppi_node_audit.csv"),
    ]
    rows = []
    for path in paths:
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
        rows.append({"path": path, "sha256": sha256(path), "size_bytes": os.path.getsize(path)})
    rows.extend([
        {"path": BASE_RUNTIME_SOURCE, "sha256": runtime_source_hash, "size_bytes": os.path.getsize(BASE_RUNTIME_SOURCE)},
        {"path": os.path.join(RUN_ROOT, "scripts", "pcsf_runtime_patched.py"), "sha256": runtime_hash,
         "size_bytes": os.path.getsize(os.path.join(RUN_ROOT, "scripts", "pcsf_runtime_patched.py"))},
        {"path": os.path.abspath(__file__), "sha256": driver_hash, "size_bytes": os.path.getsize(__file__)},
    ])
    atomic_csv(os.path.join(RUN_ROOT, "inputs", "source_artifact_sha256_manifest.csv"), pd.DataFrame(rows))
    return rows


def expected_config(runtime_source_hash, runtime_hash, driver_hash, source_manifest):
    return {
        "run_id": RUN_ID,
        "analysis": "Per-disease compact PCSF 30-50 node parameter selection and independent stability evaluation",
        "status": "compact_le50_selection_and_stability; formal_networks_unchanged",
        "analysis_root": ANALYSIS_ROOT,
        "scan_root": SCAN_ROOT,
        "model_root": r"E:\micro_output\分析审阅中心\待审阅\Fig6_疾病vsCtrl_全基因CauFinder",
        "diseases": DISEASES,
        "input_identity": {"features_per_disease": 17645, "mapped_model_entrez_per_disease": 11968,
                           "ppi_nodes": 12512, "ppi_edges": 83065, "ppi_filter": "number>=2",
                           "retain_all_background_components": True},
        "base_grid": {"b": B_VALUES, "w": W_VALUES, "g": G_VALUES,
                      "expected_rows": 1664, "source_run": SCAN_ID,
                      "recompute_base_grid": False,
                      "reuse_policy": "hash-verified b=2 and accepted b=3..6 base outputs; only missing b=0.25,0.5,1,1.5,2.5 solved"},
        "qualification": {"node_count_inclusive": [MIN_NODES, MAX_NODES],
                          "full_selected_skeleton_one_component": True,
                          "require_tree_backbone": True,
                          "all_true_nodes_included_no_component_extraction": True,
                          "no_candidate_status": "no_eligible_candidate",
                          "size_is_design_target_not_biological_truth": True},
        "selection": {"repeats": SELECTION_REPEATS, "noise_sigma": SIGMA,
                      "formula": "c'_e=(base_cost+hub_penalty)*exp(epsilon-sigma^2/2); real PPI edges only; virtual-root costs=w",
                      "paired_across_candidate_parameters_within_disease": True,
                      "seed_key": "42|Disease|compact_le50_selection_v1|replicate; SHA256 first 4 bytes little-endian; parameter_id excluded",
                      "rule": ["node_Jaccard_Q25_desc", "node_Jaccard_median_desc",
                               "baseline_full_model_SHAP_retention_desc", "b_asc", "w_asc", "g_asc"],
                      "quantile": "float64; NumPy linear interpolation; unrounded values"},
        "evaluation": {"repeats": EVALUATION_REPEATS, "stage_seed_key": "42|Disease|compact_le50_evaluation_v1|replicate",
                       "independent_from_selection": True, "selection_parameters_sha_locked": True},
        "degree_stratified_null": {"repeats": DEGREE_NULL_REPEATS,
                       "stage_seed_key": "42|Disease|compact_le50_degree_stratified_shuffle_v1|replicate",
                       "degree_bin_method": "accepted equal-degree groups accumulated to >=100; trailing remainder merged into prior bin; single small bin retained",
                       "only_mapped_model_evidence_nodes_permuted": True,
                       "measured_zero_reward_nodes_included": True,
                       "no_evidence_connectors_remain_zero": True,
                       "root_terminals_recomputed_after_shuffle": True},
        "solver": {"implementation": "OmicsIntegrator2 2.3.10 + pcst-fast 1.0.7",
                   "pruning": "strong", "dummy_mode": "terminals", "seed": 42,
                   "base_cost": 1.0, "hub_penalty_formula": "10^g*deg(u)*deg(v)/((N-deg(u)-1)*(N-deg(v)-1)+deg(u)*deg(v))",
                   "objective": "sum(selected augmented-edge costs)+sum(unselected scaled node prizes); root-edge w charged once",
                   "legal_witness_comparison": ["raw_pcst_fast", "root_only_empty", "best_single_node",
                                                "best_single_real_PPI_edge", "root_star_prize_gt_w"]},
        "software": {"python": "3.7.9", "rng": "NumPy RandomState MT19937", "numpy_version": np.__version__,
                     "pandas_version": pd.__version__},
        "script_hashes": {"driver_sha256": driver_hash, "source_runtime_sha256": runtime_source_hash,
                          "patched_runtime_sha256": runtime_hash},
        "source_artifacts": source_manifest,
        "prohibited_in_this_run": ["model_training", "SHAP_recompute", "formal_network_replacement",
                                   "drug_analysis", "GO", "cross_disease_overlap", "formal_figures"]
    }


def setup_module(mod, runtime_path):
    mod.ANALYSIS_ROOT = ANALYSIS_ROOT
    mod.MODEL_ROOT = r"E:\micro_output\分析审阅中心\待审阅\Fig6_疾病vsCtrl_全基因CauFinder"
    mod.BASE_RUN_ID = SCAN_ID
    mod.BASE_RUN = SCAN_ROOT
    mod.RUN_ID = RUN_ID
    mod.RUN = RUN_ROOT
    mod.DISEASES = list(DISEASES)
    mod.ACTIVE_DISEASES = list(DISEASES)
    mod.B_VALUES = list(B_VALUES)
    mod.W_VALUES = list(W_VALUES)
    mod.G_VALUES = list(G_VALUES)
    mod.REPEATS = SELECTION_REPEATS
    mod.SIGMA = SIGMA
    mod.PRUNE = "strong"
    mod.PRUNING = "strong"

    def stage_seed(_ignored_key):
        return int(mod.CURRENT_RANDOM_SEED)

    def staged_path(disease, pid, control_type, replicate):
        phase = str(mod.CURRENT_PHASE)
        root = "randomizations" if phase == "degree_null" else "perturbations"
        return os.path.join(RUN_ROOT, root, phase, disease, pid, control_type,
                            "replicate_{:03d}.json".format(int(replicate)))

    mod.stable_seed = stage_seed
    mod._control_run_path = staged_path

    def frozen_base_result_path(disease, pid, filename):
        key = (str(disease), str(pid))
        if key not in mod.FROZEN_BASE_RESULT_PATHS:
            raise RuntimeError("No verified frozen base-network path registered for {} {}".format(*key))
        path = os.path.join(mod.FROZEN_BASE_RESULT_PATHS[key], str(filename))
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
        return path

    mod.frozen_base_result_path = frozen_base_result_path
    return mod


def digest_json_record(mod, record):
    rec = dict(record)
    rec.pop("record_payload_sha256", None)
    canonical = json.dumps(mod._json_clean(rec), ensure_ascii=False, sort_keys=True,
                           separators=(",", ":"), allow_nan=False).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def prepare_seed_registry(eligible, diseases):
    rows, keys, seeds = [], set(), {}
    stage_repeats = [("selection", SELECTION_REPEATS),
                     ("evaluation", EVALUATION_REPEATS),
                     ("degree_null", DEGREE_NULL_REPEATS)]
    for disease in diseases:
        candidate_ids = "|".join(sorted(eligible.loc[eligible["Disease"] == disease, "parameter_id"].astype(str).unique()))
        for stage, count in stage_repeats:
            for rep in range(1, count + 1):
                key = "{}|{}|{}|{}".format(MASTER_SEED, disease, STAGE_SEED_LABELS[stage], rep)
                seed = int.from_bytes(hashlib.sha256(key.encode("utf-8")).digest()[:4], "little", signed=False)
                if key in keys:
                    raise RuntimeError("Duplicate seed key: {}".format(key))
                if seed in seeds and seeds[seed] != key:
                    raise RuntimeError("Stable 32-bit seed collision between {} and {}".format(key, seeds[seed]))
                keys.add(key); seeds[seed] = key
                rows.append({"Disease": disease, "stage": stage, "replicate": rep,
                             "seed_key": key, "seed_stage_label": STAGE_SEED_LABELS[stage], "random_seed": seed,
                             "candidate_parameter_ids_sharing_noise": candidate_ids if stage == "selection" else "",
                             "rng_algorithm": "NumPy RandomState MT19937",
                             "derivation": "SHA256(key)[0:4] interpreted little-endian unsigned"})
    out = pd.DataFrame(rows)
    if out["random_seed"].duplicated().any():
        raise RuntimeError("Random seed registry contains duplicate seeds")
    historical_paths = [
        os.path.join(ANALYSIS_ROOT, "robustness", "20261020_disease_per_disease_compact_stability_01", "config", "seed_registry.csv"),
        os.path.join(ANALYSIS_ROOT, "robustness", "20261020_disease_parameter_adaptation_08", "summaries", "run_registry.csv")]
    historical_seeds=set()
    for path in historical_paths:
        if os.path.isfile(path):
            old=read_csv_unicode(path,dtype=str,keep_default_na=False)
            col=next((c for c in ["random_seed","seed"] if c in old.columns),None)
            if col is not None:
                historical_seeds.update(int(x) for x in old[col] if str(x).strip().isdigit())
    collision=set(int(x) for x in out["random_seed"]) & historical_seeds
    if collision:
        raise RuntimeError("New stage seed collided with historical selection/evaluation seeds: {}".format(sorted(collision)[:10]))
    atomic_csv(os.path.join(RUN_ROOT, "config", "seed_registry.csv"), out)
    return out


def base_context(mod, disease, g, templates, disease_data):
    ctx = dict(templates[int(g)])
    p = disease_data[disease]["prize_df"]
    ctx["bare"] = np.asarray([float(p.loc[eid, "prize"]) for eid in ctx["nodes"]], dtype=np.float64)
    ctx["raw"] = np.asarray([float(p.loc[eid, "mean_abs_SHAP"])
                              if str(p.loc[eid, "mean_abs_SHAP"]) not in ("", "nan", "NaN") else np.nan
                              for eid in ctx["nodes"]], dtype=np.float64)
    ctx["evidence"] = np.isfinite(ctx["raw"])
    rebuilt = np.asarray(ctx["base_costs"] + ctx["edge_penalties"], dtype=np.float64)
    if "frozen_effective_costs" not in ctx or not np.array_equal(
            np.asarray(ctx["frozen_effective_costs"], dtype=np.float64), rebuilt):
        raise RuntimeError("Loaded frozen effective-cost vector differs from base cost plus hub penalty for g={}".format(g))
    ctx["frozen_effective_costs"] = rebuilt
    return ctx


def create_degree_bins(mod, eligible, templates, disease_data):
    active = sorted(set(eligible["Disease"]))
    bins_by_disease, rows = {}, []
    for disease in active:
        ctx = base_context(mod, disease, G_VALUES[0], templates, disease_data)
        evidence_ix = np.flatnonzero(ctx["evidence"])
        bins = mod.degree_bins(evidence_ix, ctx["degree"], min_size=100)
        bins_by_disease[disease] = bins
        for bin_id, item in enumerate(bins, 1):
            for ix in item["indices"]:
                rows.append({"Disease": disease, "degree_bin": bin_id,
                             "degree_min": item["degree_min"], "degree_max": item["degree_max"],
                             "Entrez": ctx["nodes"][int(ix)],
                             "background_degree": int(ctx["degree"][int(ix)]),
                             "bin_size": len(item["indices"]),
                             "tail_handling": "merge trailing remainder into prior bin; retain one small bin if only bin"})
    atomic_csv(os.path.join(RUN_ROOT, "config", "degree_stratified_reward_shuffle_bins.csv"), pd.DataFrame(rows))
    return bins_by_disease


def run_record(mod, disease, b, w, g, control_type, replicate, phase,
               ppi_df, ctx, disease_data, degree_bins, config_hash,
               input_manifest_hash, composite_script_hash, selected_lock_hash,
               candidate_set_hash):
    key = "{}|{}|{}|{}".format(MASTER_SEED, disease, STAGE_SEED_LABELS[phase], replicate)
    seed = int.from_bytes(hashlib.sha256(key.encode("utf-8")).digest()[:4], "little", signed=False)
    mod.CURRENT_PHASE = phase
    mod.CURRENT_RANDOM_SEED = seed
    mod.CURRENT_SEED_KEY = key
    mod.CURRENT_DISEASE = disease
    mod.CURRENT_REPLICATE = int(replicate)
    rec = mod.run_one_control(disease, float(b), int(w), int(g), control_type,
                              int(replicate), ppi_df, ctx, disease_data,
                              degree_bins, config_hash, input_manifest_hash,
                              composite_script_hash)
    rec = dict(rec)
    if rec.get("stage_id") not in (None, phase):
        raise RuntimeError("Existing repeat checkpoint belongs to a different stage")
    rec["stage_id"] = phase
    rec["seed_key"] = key
    rec["seed_derivation"] = "SHA256(key)[0:4] interpreted little-endian unsigned; key excludes parameter_id"
    rec["parameter_noise_pairing_group"] = "{}|{}|{}".format(disease, phase, replicate) if phase == "selection" else ""
    rec["locked_selected_parameters_sha256"] = selected_lock_hash or ""
    rec["selection_candidate_set_sha256"] = candidate_set_hash or ""
    if int(rec.get("random_seed", -1)) != seed:
        raise RuntimeError("Solver runtime did not use the frozen stable seed")
    if control_type == "edge_cost":
        rec["perturbation_version"] = "effective_cost_full_v2"
        rec["cost_perturbation_formula"] = "(base_cost+hub_penalty)*exp(epsilon-sigma^2/2); all real PPI edges; virtual root edges unchanged at w"
        rng = np.random.RandomState(seed)
        z = np.asarray(rng.normal(0.0, SIGMA, size=len(ctx["base_costs"])), dtype=np.float64)
        multiplier = np.exp(z - (SIGMA ** 2) / 2.0)
        frozen_cost = np.asarray(ctx["frozen_effective_costs"], dtype=np.float64)
        perturbed_cost = frozen_cost * multiplier
        vectors_dir = os.path.join(RUN_ROOT, "inputs", "random_vectors", phase, disease)
        os.makedirs(vectors_dir, exist_ok=True)
        vector_path = os.path.join(vectors_dir, "cost_noise_rep_{:03d}.npz".format(int(replicate)))
        if os.path.isfile(vector_path):
            with np.load(vector_path, allow_pickle=False) as saved:
                if not np.array_equal(saved["z"], z) or not np.array_equal(saved["multiplier"], multiplier):
                    raise RuntimeError("Persisted paired cost-noise vector/multiplier differs from regenerated values")
        else:
            tmp_path = vector_path + ".tmp.npz"
            np.savez_compressed(tmp_path, z=z, multiplier=multiplier)
            os.replace(tmp_path, vector_path)
        rec["perturbation_vector_path"] = vector_path
        rec["perturbation_vector_file_sha256"] = sha256(vector_path)
        rec["cost_noise_z_sha256"] = bytes_sha256(z.tobytes(order="C"))
        rec["cost_noise_vector_sha256"] = rec["cost_noise_z_sha256"]
        rec["cost_multiplier_sha256"] = bytes_sha256(multiplier.tobytes(order="C"))
        rec["frozen_effective_real_costs_sha256"] = bytes_sha256(frozen_cost.tobytes(order="C"))
        rec["perturbed_real_costs_sha256"] = bytes_sha256(perturbed_cost.tobytes(order="C"))
        rec["virtual_root_edge_cost"] = int(w)
        rec["virtual_root_cost_unchanged"] = True
        rec["solver_input_real_costs_sha256"] = rec["perturbed_real_costs_sha256"]
        rec["solver_input_root_costs_sha256"] = bytes_sha256(np.full(int(rec["solver_call"]["n_edges"]) - len(z), float(w), dtype=np.float64).tobytes(order="C"))
        rec["base_reward_vector_sha256"] = bytes_sha256(np.asarray(ctx["bare"], dtype=np.float64).tobytes(order="C"))
        rec["altered_reward_vector_sha256"] = rec["base_reward_vector_sha256"]
        rec["rng_replay_contract"] = "RandomState(MT19937, seed).normal(0,0.1,n_real_edges), frozen edge order in source solver arrays; multiply frozen effective cost; roots fixed w"
    elif control_type == "degree_stratified_reward_shuffle":
        rec["perturbation_version"] = "degree_stratified_reward_shuffle_selection_v1"
        rec["reward_permutation_contract"] = "within frozen degree bins over ordered mapped evidence nodes; measured-zero evidence included; no-evidence connectors fixed at zero"
        evidence_ix = np.flatnonzero(np.asarray(ctx["evidence"], dtype=bool))
        permutation_indices = np.arange(len(evidence_ix), dtype=np.int32)
        rng = np.random.RandomState(seed)
        for bin_row in degree_bins[disease]:
            members = np.asarray(bin_row["indices"], dtype=np.int64)
            positions = np.searchsorted(evidence_ix, members)
            permutation = rng.permutation(len(members))
            permutation_indices[positions] = positions[permutation]
        altered = np.asarray(ctx["bare"], dtype=np.float64).copy()
        altered[evidence_ix] = altered[evidence_ix[permutation_indices]]
        vectors_dir = os.path.join(RUN_ROOT, "inputs", "random_vectors", phase, disease)
        os.makedirs(vectors_dir, exist_ok=True)
        vector_path = os.path.join(vectors_dir, "reward_permutation_rep_{:03d}.npz".format(int(replicate)))
        if os.path.isfile(vector_path):
            with np.load(vector_path, allow_pickle=False) as saved:
                if not np.array_equal(saved["permutation_indices"], permutation_indices):
                    raise RuntimeError("Persisted reward permutation differs from regenerated permutation")
        else:
            tmp_path = vector_path + ".tmp.npz"
            np.savez_compressed(tmp_path, permutation_indices=permutation_indices)
            os.replace(tmp_path, vector_path)
        rec["perturbation_vector_path"] = vector_path
        rec["perturbation_vector_file_sha256"] = sha256(vector_path)
        rec["reward_permutation_indices_sha256"] = bytes_sha256(permutation_indices.tobytes(order="C"))
        rec["base_reward_vector_sha256"] = bytes_sha256(np.asarray(ctx["bare"], dtype=np.float64).tobytes(order="C"))
        rec["altered_reward_vector_sha256"] = bytes_sha256(altered.tobytes(order="C"))
        rec["frozen_effective_real_costs_sha256"] = bytes_sha256(np.asarray(ctx["frozen_effective_costs"], dtype=np.float64).tobytes(order="C"))
        rec["solver_input_real_costs_sha256"] = rec["frozen_effective_real_costs_sha256"]
        rec["solver_input_root_costs_sha256"] = bytes_sha256(np.full(int(rec["solver_call"]["n_edges"]) - len(ctx["base_costs"]), float(w), dtype=np.float64).tobytes(order="C"))
        rec["rng_replay_contract"] = "RandomState(MT19937, seed); sequential permutation within saved degree-bin order; evidence-node order in solver arrays; no-evidence connectors fixed zero"
    rec["record_payload_sha256"] = digest_json_record(mod, rec)
    path = mod._control_run_path(disease, parameter_id(mod, b, w, g), control_type, replicate)
    mod.write_json(path, mod._json_clean(rec))
    return rec, path


def checked_record(mod, path, disease, pid, control_type, rep,
                   config_hash, input_manifest_hash, code_hash, stage):
    with open(path, "r", encoding="utf-8-sig") as handle:
        rec = json.load(handle)
    if (rec.get("status") != "complete" or rec.get("Disease") != disease or
            rec.get("parameter_id") != pid or rec.get("control_type") != control_type or
            int(rec.get("replicate", -1)) != rep or rec.get("stage_id") != stage or
            rec.get("run_config_sha256") != config_hash or
            rec.get("input_manifest_sha256") != input_manifest_hash or
            rec.get("execution_script_sha256") != code_hash or
            rec.get("record_payload_sha256") != digest_json_record(mod, rec) or
            rec.get("is_feasible") is not True):
        raise RuntimeError("Repeat checkpoint identity, hash, status or feasibility failure: {}".format(path))
    return rec


def linear_quantile(values, probability):
    vals = np.asarray(values, dtype=np.float64)
    vals = vals[np.isfinite(vals)]
    if not len(vals):
        return None
    return float(np.quantile(vals, probability, interpolation="linear"))


def numeric_distribution(prefix, values):
    vals = np.asarray([float(x) for x in values if x is not None and np.isfinite(float(x))], dtype=np.float64)
    if not len(vals):
        return {prefix + "_n": 0, prefix + "_median": None, prefix + "_q25": None,
                prefix + "_q75": None, prefix + "_q05": None, prefix + "_q95": None}
    return {prefix + "_n": int(len(vals)), prefix + "_median": float(np.median(vals)),
            prefix + "_q25": linear_quantile(vals, .25), prefix + "_q75": linear_quantile(vals, .75),
            prefix + "_q05": linear_quantile(vals, .05), prefix + "_q95": linear_quantile(vals, .95)}


def read_stage_records(mod, disease, b, w, g, phase, control_type, repeats,
                       config_hash, input_manifest_hash, code_hash):
    pid = parameter_id(mod, b, w, g)
    root = "randomizations" if phase == "degree_null" else "perturbations"
    records = []
    for rep in range(1, repeats + 1):
        path = os.path.join(RUN_ROOT, root, phase, disease, pid, control_type,
                            "replicate_{:03d}.json".format(rep))
        if not os.path.isfile(path):
            raise RuntimeError("Missing required repeat {} {}/{}: {}".format(disease, pid, phase, path))
        records.append(checked_record(mod, path, disease, pid, control_type, rep,
                                      config_hash, input_manifest_hash, code_hash, phase))
    return records


def summarize_candidates(mod, eligibility, input_manifest_hash, config_hash, code_hash):
    eligible = eligibility[eligibility["eligible"].astype(str).str.lower().isin(["true", "1"])].copy()
    rows, long = [], []
    for _, c in eligible.iterrows():
        disease, b, w, g, pid = str(c["Disease"]), float(c["b"]), int(c["w"]), int(c["g"]), str(c["parameter_id"])
        records = read_stage_records(mod, disease, b, w, g, "selection", "edge_cost",
                                     SELECTION_REPEATS, config_hash, input_manifest_hash, code_hash)
        if len(records) != SELECTION_REPEATS:
            raise RuntimeError("Selection candidate does not have exactly 20 repeats: {} {}".format(disease, pid))
        for rec in records:
            metrics = rec.get("network_metrics", {})
            long.append({"Disease": disease, "parameter_id": pid, "b": b, "w": w, "g": g,
                         "replicate": int(rec["replicate"]), "random_seed": int(rec["random_seed"]),
                         "seed_key": rec["seed_key"], "cost_noise_vector_sha256": rec.get("cost_noise_vector_sha256"),
                         "edge_order_sha256": rec.get("edge_order_sha256"),
                         "node_jaccard": rec.get("node_jaccard_to_unperturbed_adopted"),
                         "backbone_edge_jaccard": rec.get("backbone_edge_jaccard_to_unperturbed_adopted"),
                         "real_node_count": metrics.get("real_node_count"),
                         "real_skeleton_edge_count": metrics.get("real_skeleton_edge_count"),
                         "component_count": metrics.get("component_count"),
                         "isolated_node_count": metrics.get("isolated_node_count"),
                         "largest_component_fraction": metrics.get("largest_component_fraction"),
                         "retained_SHAP_mass_fraction_full_model": metrics.get("retained_SHAP_mass_fraction_full_model"),
                         "retained_SHAP_mass_fraction_mapped_PPI": metrics.get("retained_SHAP_mass_fraction_mapped_PPI"),
                         "status": rec.get("status"), "is_feasible": rec.get("is_feasible"),
                         "adopted_solution_source": rec.get("adopted_solution_source"),
                         "objective_recomputed": rec.get("objective_recomputed"),
                         "repeat_record_path": mod._control_run_path(disease, pid, "edge_cost", int(rec["replicate"]))})
        node_j = [float(r["node_jaccard_to_unperturbed_adopted"]) for r in records]
        edge_j = [float(r["backbone_edge_jaccard_to_unperturbed_adopted"]) for r in records]
        rows.append({"Disease": disease, "parameter_id": pid, "b": b, "w": w, "g": g,
                     "selection_repeats_complete": len(records),
                     "selection_node_Jaccard_Q25_linear": linear_quantile(node_j, .25),
                     "selection_node_Jaccard_median": float(np.median(np.asarray(node_j, dtype=np.float64))),
                     "selection_node_Jaccard_min": float(np.min(node_j)),
                     "selection_node_Jaccard_max": float(np.max(node_j)),
                     "selection_edge_Jaccard_Q25_linear": linear_quantile(edge_j, .25),
                     "selection_edge_Jaccard_median": float(np.median(np.asarray(edge_j, dtype=np.float64))),
                     "baseline_node_count": int(c["real_node_count"]),
                     "baseline_edge_count": int(c["real_skeleton_edge_count"]),
                     "baseline_full_model_SHAP_retention": float(c["retained_SHAP_mass_fraction_full_model"]),
                     "baseline_mapped_PPI_SHAP_retention": float(c["retained_SHAP_mass_fraction_mapped_PPI"]),
                     "baseline_no_evidence_connector_fraction": float(c["zero_reward_no_evidence_connector_count"]) / float(c["real_node_count"]),
                     "baseline_adopted_solution_source": c["adopted_solution_source"],
                     "selection_node_Jaccard_all20_finite": all(np.isfinite(node_j)),
                     "candidate_eligibility_sha256": sha256(os.path.join(RUN_ROOT, "summaries", "parameter_qualification_all_104.csv"))})
    if not rows:
        raise RuntimeError("No disease has an eligible compact base network")
    # One selection seed and edge-noise vector must be shared by every eligible
    # parameter for the same disease and repetition.
    paired = pd.DataFrame(long)
    for (disease, replicate), group in paired.groupby(["Disease", "replicate"], sort=False):
        if group["parameter_id"].nunique() > 1:
            for col in ["random_seed", "cost_noise_vector_sha256", "edge_order_sha256"]:
                if group[col].nunique(dropna=False) != 1:
                    raise RuntimeError("Paired selection noise differs across parameter candidates: {} rep {} field {}".format(disease, replicate, col))
    summary = pd.DataFrame(rows)
    ranked = []
    for disease, group in summary.groupby("Disease", sort=False):
        ordered = sorted(group.to_dict(orient="records"),
                         key=lambda x: (-float(x["selection_node_Jaccard_Q25_linear"]),
                                        -float(x["selection_node_Jaccard_median"]),
                                        -float(x["baseline_full_model_SHAP_retention"]),
                                        float(x["b"]), int(x["w"]), int(x["g"])))
        for rank, row in enumerate(ordered, 1):
            row["selection_rank"] = rank
            row["selected"] = rank == 1
            row["selection_rule"] = "node Jaccard Q25 desc -> node Jaccard median desc -> baseline full-model SHAP retention desc -> numeric b,w,g asc"
            ranked.append(row)
    ranked_df = pd.DataFrame(ranked).sort_values(["Disease", "selection_rank"]).reset_index(drop=True)
    atomic_csv(os.path.join(RUN_ROOT, "summaries", "selection_candidate_summary.csv"), ranked_df)
    atomic_csv(os.path.join(RUN_ROOT, "summaries", "selection_candidate_20repeats.csv"), pd.DataFrame(long))
    chosen = ranked_df[ranked_df["selected"]].copy()
    final = []
    for _, row in chosen.iterrows():
        base = eligibility[(eligibility["Disease"] == row["Disease"]) &
                           (eligibility["parameter_id"] == row["parameter_id"])].iloc[0]
        final.append({"Disease": row["Disease"], "status": "selected_parameter_locked",
                      "b": float(row["b"]), "w": int(row["w"]), "g": int(row["g"]),
                      "parameter_id": row["parameter_id"],
                      "node_count": int(base["real_node_count"]), "edge_count": int(base["real_skeleton_edge_count"]),
                      "component_count": int(base["component_count"]),
                      "selection_repeat_count": int(row["selection_repeats_complete"]),
                      "selection_node_Jaccard_Q25_linear": float(row["selection_node_Jaccard_Q25_linear"]),
                      "selection_node_Jaccard_median": float(row["selection_node_Jaccard_median"]),
                      "baseline_full_model_SHAP_retention": float(row["baseline_full_model_SHAP_retention"]),
                      "baseline_mapped_PPI_SHAP_retention": float(row["baseline_mapped_PPI_SHAP_retention"]),
                      "selection_edge_Jaccard_Q25_linear": float(row["selection_edge_Jaccard_Q25_linear"]),
                      "selection_edge_Jaccard_median": float(row["selection_edge_Jaccard_median"]),
                      "selection_rank": int(row["selection_rank"]),
                      "selection_reason": row["selection_rule"],
                      "base_selected_nodes_path": base["selected_nodes_path"],
                      "base_selected_nodes_sha256": base["selected_nodes_sha256"],
                      "base_selected_edges_path": base["selected_edges_path"],
                      "base_selected_edges_sha256": base["selected_edges_sha256"],
                      "base_run_manifest_path": base["base_run_manifest_path"],
                      "base_run_manifest_sha256": base["base_run_manifest_sha256"]})
    selected_df = pd.DataFrame(final).sort_values("Disease").reset_index(drop=True)
    selected_path = os.path.join(RUN_ROOT, "summaries", "selected_parameters.csv")
    atomic_csv(selected_path, selected_df)
    selected_sha = sha256(selected_path)
    lock_path = os.path.join(RUN_ROOT, "config", "selected_parameters.lock.json")
    if os.path.isfile(lock_path):
        with open(lock_path, "r", encoding="utf-8-sig") as handle:
            old_lock = json.load(handle)
        if old_lock.get("selected_parameters_sha256") != selected_sha:
            raise RuntimeError("Existing selected parameter lock differs from newly computed selection")
        lock = old_lock
    else:
        lock = {"status": "selection_locked_for_evaluation", "selected_parameters_path": selected_path,
                "selected_parameters_sha256": selected_sha,
                "selected_candidate_summary_sha256": sha256(os.path.join(RUN_ROOT, "summaries", "selection_candidate_summary.csv")),
                "selection_repeat_table_sha256": sha256(os.path.join(RUN_ROOT, "summaries", "selection_candidate_20repeats.csv")),
                "qualification_table_sha256": sha256(os.path.join(RUN_ROOT, "summaries", "parameter_qualification_all_104.csv")),
                "run_config_sha256": config_hash, "input_manifest_sha256": input_manifest_hash,
                "selection_order": "per-disease node-Jaccard Q25 linear desc, median desc, baseline full-model SHAP retention desc, then numeric b,w,g asc",
                "locked_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
        write_json(lock_path, lock)
    return ranked_df, selected_df, lock


def load_locked_parameters(mod, config_hash):
    lock_path = os.path.join(RUN_ROOT, "config", "selected_parameters.lock.json")
    selected_path = os.path.join(RUN_ROOT, "summaries", "selected_parameters.csv")
    if not os.path.isfile(lock_path) or not os.path.isfile(selected_path):
        raise RuntimeError("Selection lock is not present; evaluation cannot begin")
    with open(lock_path, "r", encoding="utf-8-sig") as handle:
        lock = json.load(handle)
    actual = sha256(selected_path)
    if actual != lock.get("selected_parameters_sha256") or lock.get("run_config_sha256") != config_hash:
        raise RuntimeError("Selected parameter lock/config SHA mismatch")
    return mod.read_table(selected_path), lock, actual


def run_stage(mod, stage, selected, eligibility, ppi_df, disease_data,
              bins_by_disease, templates, config_hash, input_hash, code_hash,
              lock_hash, candidate_hash):
    if stage == "selection":
        frame = eligibility[eligibility["eligible"].astype(str).str.lower().isin(["true", "1"])].copy()
        repeat_count, control_type = SELECTION_REPEATS, "edge_cost"
        groups = [(str(r["Disease"]), float(r["b"]), int(r["w"]), int(r["g"]))
                  for _, r in frame.sort_values(["Disease", "b", "w", "g"]).iterrows()]
    elif stage == "evaluation":
        repeat_count, control_type = EVALUATION_REPEATS, "edge_cost"
        groups = [(str(r["Disease"]), float(r["b"]), int(r["w"]), int(r["g"]))
                  for _, r in selected.sort_values("Disease").iterrows()]
    elif stage == "degree_null":
        repeat_count, control_type = DEGREE_NULL_REPEATS, "degree_stratified_reward_shuffle"
        groups = [(str(r["Disease"]), float(r["b"]), int(r["w"]), int(r["g"]))
                  for _, r in selected.sort_values("Disease").iterrows()]
    else:
        raise ValueError("Unknown stage {}".format(stage))
    total = len(groups) * repeat_count
    done = 0
    progress = {"stage": stage, "status": "running", "planned_records": total,
                "completed_records": 0, "last_disease": "", "last_parameter_id": "",
                "last_replicate": 0, "updated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    progress_path = os.path.join(RUN_ROOT, "logs", "progress_{}.json".format(stage))
    for disease, b, w, g in groups:
        ctx = base_context(mod, disease, g, templates, disease_data)
        pid = parameter_id(mod, b, w, g)
        for rep in range(1, repeat_count + 1):
            run_record(mod, disease, b, w, g, control_type, rep, stage,
                       ppi_df, ctx, disease_data, bins_by_disease, config_hash,
                       input_hash, code_hash, lock_hash, candidate_hash)
            done += 1
            progress.update({"completed_records": done, "last_disease": disease,
                             "last_parameter_id": pid, "last_replicate": rep,
                             "updated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")})
            write_json(progress_path, progress)
        log("stage={} disease={} parameter={} repeats={}/{} stage_records={}/{}".format(
            stage, disease, pid, repeat_count, repeat_count, done, total))
    progress.update({"status": "complete", "completed_records": done,
                     "completed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")})
    write_json(progress_path, progress)
    if done != total:
        raise RuntimeError("Stage did not finish its planned record count")


def summarize_repeat_set(mod, records):
    metrics = [r.get("network_metrics", {}) for r in records]
    node_j = [r.get("node_jaccard_to_unperturbed_adopted") for r in records]
    edge_j = [r.get("backbone_edge_jaccard_to_unperturbed_adopted") for r in records]
    names = ["real_node_count", "real_skeleton_edge_count", "component_count", "isolated_node_count",
             "largest_component_fraction", "nonisolated_node_fraction", "connector_fraction",
             "retained_SHAP_mass_fraction_full_model", "retained_SHAP_mass_fraction_mapped_PPI",
             "selected_node_fraction_top1pct_degree", "skeleton_edge_fraction_incident_top1pct_degree"]
    out = {"repeat_count": len(records), "empty_solution_count": sum(int(m.get("real_node_count", 0)) == 0 for m in metrics),
           "empty_solution_fraction": float(sum(int(m.get("real_node_count", 0)) == 0 for m in metrics)) / len(records),
           "single_connected_network_count": sum(int(m.get("real_node_count", 0)) > 0 and int(m.get("component_count", 0)) == 1 for m in metrics),
           "single_connected_network_fraction": float(sum(int(m.get("real_node_count", 0)) > 0 and int(m.get("component_count", 0)) == 1 for m in metrics)) / len(records),
           "within_30_50_count": sum(MIN_NODES <= int(m.get("real_node_count", 0)) <= MAX_NODES for m in metrics),
           "within_30_50_fraction": float(sum(MIN_NODES <= int(m.get("real_node_count", 0)) <= MAX_NODES for m in metrics)) / len(records),
           "at_most_50_count": sum(int(m.get("real_node_count", 0)) <= MAX_NODES for m in metrics),
           "at_most_50_fraction": float(sum(int(m.get("real_node_count", 0)) <= MAX_NODES for m in metrics)) / len(records)}
    for field in names:
        values = [m.get(field) for m in metrics]
        for key, value in numeric_distribution(field, values).items():
            out[key] = value
    for name, values in [("node_Jaccard", node_j), ("backbone_edge_Jaccard", edge_j)]:
        for key, value in numeric_distribution(name, values).items():
            out[key] = value
    return out


def frequency_tables(mod, disease, parameter_id_value, records, base_nodes, base_edges,
                     denominator, stage):
    node_counts, edge_counts = {}, {}
    for rec in records:
        for eid in set(map(str, rec.get("selected_node_Entrez", []))):
            node_counts[eid] = node_counts.get(eid, 0) + 1
        for edge in set(map(str, rec.get("selected_edge_keys", []))):
            edge_counts[edge] = edge_counts.get(edge, 0) + 1
    node_rows = [{"Disease": disease, "parameter_id": parameter_id_value,
                  "stage": stage, "Entrez": eid, "in_unperturbed_base_network": eid in base_nodes,
                  "selected_count": int(count), "frequency_denominator": int(denominator),
                  "selection_frequency": float(count) / float(denominator)}
                 for eid, count in ((eid, node_counts.get(eid, 0))
                                    for eid in sorted(set(node_counts) | set(base_nodes)))]
    edge_rows = [{"Disease": disease, "parameter_id": parameter_id_value,
                  "stage": stage, "edge_key": edge, "in_unperturbed_base_backbone": edge in base_edges,
                  "selected_count": int(count), "frequency_denominator": int(denominator),
                  "selection_frequency": float(count) / float(denominator)}
                 for edge, count in ((edge, edge_counts.get(edge, 0))
                                     for edge in sorted(set(edge_counts) | set(base_edges)))]
    return node_rows, edge_rows


def save_summary_for_stage(mod, stage, selected, config_hash, input_hash, code_hash):
    repeats = EVALUATION_REPEATS if stage == "evaluation" else DEGREE_NULL_REPEATS
    control_type = "edge_cost" if stage == "evaluation" else "degree_stratified_reward_shuffle"
    stage_rows, node_freq, edge_freq = [], [], []
    for _, chosen in selected.sort_values("Disease").iterrows():
        disease = str(chosen["Disease"]); b, w, g = float(chosen["b"]), int(chosen["w"]), int(chosen["g"])
        pid = str(chosen["parameter_id"])
        records = read_stage_records(mod, disease, b, w, g, stage, control_type,
                                     repeats, config_hash, input_hash, code_hash)
        summary = summarize_repeat_set(mod, records)
        summary.update({"Disease": disease, "parameter_id": pid, "b": b, "w": w, "g": g,
                        "stage": stage, "frequency_denominator": repeats,
                        "interpretation": "frozen-input network perturbation" if stage == "evaluation" else "degree-stratified reward null; descriptive, not significance testing"})
        stage_rows.append(summary)
        base_nodes_path = chosen["base_selected_nodes_path"]
        base_edges_path = chosen["base_selected_edges_path"]
        nodes = mod.read_table(base_nodes_path)
        edges = mod.read_table(base_edges_path, "\t")
        bn = set(nodes["Entrez"].astype(str))
        be = {"{}|{}".format(*sorted((str(a), str(b)))) for a, b in zip(edges["Entrez_A"], edges["Entrez_B"])}
        nr, er = frequency_tables(mod, disease, pid, records, bn, be, repeats, stage)
        node_freq.extend(nr); edge_freq.extend(er)
    summary_frame = pd.DataFrame(stage_rows).sort_values("Disease")
    if stage == "evaluation":
        stem = "evaluation"
    else:
        stem = "degree_null"
    atomic_csv(os.path.join(RUN_ROOT, "summaries", stem + "_network_distribution_by_disease.csv"), summary_frame)
    atomic_csv(os.path.join(RUN_ROOT, "summaries", stem + "_node_selection_frequency.csv"), pd.DataFrame(node_freq))
    atomic_csv(os.path.join(RUN_ROOT, "summaries", stem + "_edge_selection_frequency.csv"), pd.DataFrame(edge_freq))
    long_rows = []
    for _, chosen in selected.sort_values("Disease").iterrows():
        disease = str(chosen["Disease"]); b, w, g = float(chosen["b"]), int(chosen["w"]), int(chosen["g"])
        pid = str(chosen["parameter_id"])
        records = read_stage_records(mod, disease, b, w, g, stage, control_type,
                                     repeats, config_hash, input_hash, code_hash)
        for rec in records:
            row = {k: rec.get(k) for k in ["Disease", "parameter_id", "b", "w", "g", "stage_id", "control_type",
                                           "replicate", "random_seed", "seed_key", "status", "termination_status",
                                           "is_feasible", "adopted_solution_source", "objective_recomputed",
                                           "selected_reward", "unselected_reward", "selected_augmented_edge_cost",
                                           "solver_runtime_seconds", "node_jaccard_to_unperturbed_adopted",
                                           "backbone_edge_jaccard_to_unperturbed_adopted"]}
            row["repeat_record_path"] = mod._control_run_path(disease, pid, control_type, int(rec["replicate"]))
            row["record_sha256"] = sha256(row["repeat_record_path"])
            row.update({"network_" + key: val for key, val in rec.get("network_metrics", {}).items()})
            long_rows.append(row)
    atomic_csv(os.path.join(RUN_ROOT, "summaries", stem + "_repeat_metrics.csv"), pd.DataFrame(long_rows))
    return summary_frame


def copy_selected_networks(mod, selected, selection_lock):
    output_rows = []
    delivery_root = os.path.join(RUN_ROOT, "selected_networks")
    os.makedirs(delivery_root, exist_ok=True)
    for _, row in selected.sort_values("Disease").iterrows():
        disease, pid = str(row["Disease"]), str(row["parameter_id"])
        dest = os.path.join(delivery_root, disease)
        os.makedirs(dest, exist_ok=True)
        nodes_src = str(row["base_selected_nodes_path"])
        edges_src = str(row["base_selected_edges_path"])
        nodes_dst = os.path.join(dest, "selected_nodes.csv")
        edges_dst = os.path.join(dest, "selected_forest_edges.csv")
        if sha256(nodes_src) != str(row["base_selected_nodes_sha256"]) or sha256(edges_src) != str(row["base_selected_edges_sha256"]):
            raise RuntimeError("Base network changed before final delivery for {}".format(disease))
        shutil.copy2(nodes_src, nodes_dst)
        edge_frame = mod.read_table(edges_src, "\t")
        atomic_csv(edges_dst, edge_frame)
        manifest = {"Disease": disease, "status": "awaiting_user_review",
                    "b": float(row["b"]), "w": int(row["w"]), "g": int(row["g"]),
                    "parameter_id": pid, "selection_rule": row["selection_reason"],
                    "selected_parameters_sha256": selection_lock["selected_parameters_sha256"],
                    "selected_parameters_lock_path": os.path.join(RUN_ROOT, "config", "selected_parameters.lock.json"),
                    "base_run_id": str(row.get("source_run_id", SCAN_ID)), "base_run_manifest_path": str(row["base_run_manifest_path"]),
                    "base_run_manifest_sha256": str(row["base_run_manifest_sha256"]),
                    "source_selected_nodes_path": nodes_src, "source_selected_nodes_sha256": sha256(nodes_src),
                    "source_selected_edges_path": edges_src, "source_selected_edges_sha256": sha256(edges_src),
                    "outputs": [{"path": nodes_dst, "sha256": sha256(nodes_dst), "size_bytes": os.path.getsize(nodes_dst)},
                                {"path": edges_dst, "sha256": sha256(edges_dst), "size_bytes": os.path.getsize(edges_dst)}],
                    "model_retrained": False, "SHAP_recomputed": False, "formal_network_replaced": False}
        man_path = os.path.join(dest, "selected_parameter_manifest.json")
        write_json(man_path, manifest)
        output_rows.append({"Disease": disease, "status": manifest["status"], "b": row["b"], "w": row["w"], "g": row["g"],
                            "parameter_id": pid, "node_count": int(row["node_count"]), "edge_count": int(row["edge_count"]),
                            "selected_nodes_path": nodes_dst, "selected_nodes_sha256": sha256(nodes_dst),
                            "selected_edges_path": edges_dst, "selected_edges_sha256": sha256(edges_dst),
                            "manifest_path": man_path, "manifest_sha256": sha256(man_path)})
    atomic_csv(os.path.join(RUN_ROOT, "summaries", "selected_network_delivery_registry.csv"), pd.DataFrame(output_rows))


def finalize(mod, selected, lock, config_hash, input_hash, code_hash):
    eval_summary = mod.read_table(os.path.join(RUN_ROOT, "summaries", "evaluation_network_distribution_by_disease.csv"))
    null_summary = mod.read_table(os.path.join(RUN_ROOT, "summaries", "degree_null_network_distribution_by_disease.csv"))
    merged = selected.merge(eval_summary, on=["Disease", "parameter_id", "b", "w", "g"], how="left", suffixes=("", "_evaluation"))
    null_keep = ["Disease", "parameter_id", "empty_solution_fraction", "single_connected_network_fraction",
                 "within_30_50_fraction", "at_most_50_fraction", "node_Jaccard_median", "node_Jaccard_q25", "node_Jaccard_q75",
                 "backbone_edge_Jaccard_median", "backbone_edge_Jaccard_q25", "retained_SHAP_mass_fraction_full_model_median",
                 "retained_SHAP_mass_fraction_mapped_PPI_median"]
    null_sub = null_summary[[c for c in null_keep if c in null_summary.columns]].copy()
    null_sub = null_sub.rename(columns={c: "degree_null_" + c for c in null_sub.columns if c not in ["Disease", "parameter_id"]})
    merged = merged.merge(null_sub, on=["Disease", "parameter_id"], how="left")
    all_rows = []
    chosen_map = {str(r["Disease"]): r for _, r in merged.iterrows()}
    no_eligible = set(mod.read_table(os.path.join(RUN_ROOT, "summaries", "no_eligible_candidates.csv"))["Disease"].astype(str))
    for disease in DISEASES:
        if disease in chosen_map:
            row = chosen_map[disease].to_dict()
            row["status"] = "complete_selected_and_evaluated"
            row["limitations"] = "network is model-contribution-supported and constrained by PPI; stability is conditional on frozen model/reward/PPI and is not external or experimental validation"
            all_rows.append(row)
        else:
            all_rows.append({"Disease": disease, "status": "no_eligible_candidate" if disease in no_eligible else "blocked_input_identity",
                             "b": None, "w": None, "g": None, "parameter_id": None,
                             "limitations": "No candidate entered selection or perturbation evaluation; no parameter-grid expansion was performed."})
    final_frame = pd.DataFrame(all_rows)
    atomic_csv(os.path.join(RUN_ROOT, "summaries", "final_disease_summary.csv"), final_frame)
    copy_selected_networks(mod, selected, lock)

    # Validate every planned stage record count and input identity before declaring complete.
    expected = {"selection": int(sum(1 for _ in open(os.path.join(RUN_ROOT, "summaries", "selection_candidate_summary.csv"), "r", encoding="utf-8-sig")) - 1) * SELECTION_REPEATS,
                "evaluation": len(selected) * EVALUATION_REPEATS,
                "degree_null": len(selected) * DEGREE_NULL_REPEATS}
    actual = {}
    for phase, root_name in [("selection", "perturbations"), ("evaluation", "perturbations"), ("degree_null", "randomizations")]:
        root = os.path.join(RUN_ROOT, root_name, phase)
        actual[phase] = sum(1 for base, _, files in os.walk(root)
                            for name in files if name.startswith("replicate_") and name.endswith(".json")) if os.path.isdir(root) else 0
    if actual != expected:
        raise RuntimeError("Final repeat counts differ: expected {}, actual {}".format(expected, actual))
    output_files = []
    for base, _, files in os.walk(RUN_ROOT):
        for name in files:
            path = os.path.join(base, name)
            if name.endswith((".csv", ".json", ".md", ".py", ".npz")) and os.path.isfile(path):
                output_files.append({"path": path, "sha256": sha256(path), "size_bytes": os.path.getsize(path)})
    run_manifest = {"run_id": RUN_ID, "status": "complete", "updated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                    "selected_disease_count": len(selected), "no_eligible_disease_count": len(no_eligible),
                    "planned_repeat_counts": expected, "observed_repeat_counts": actual,
                    "selected_parameters_sha256": lock["selected_parameters_sha256"],
                    "run_config_sha256": config_hash, "input_manifest_sha256": input_hash,
                    "execution_code_sha256": code_hash, "outputs": output_files,
                    "formal_networks_replaced": False, "drug_analysis_run": False, "formal_figures_run": False, "formal_networks_replaced": False, "drug_analysis_run": False}
    write_json(os.path.join(RUN_ROOT, "summaries", "run_manifest.json"), run_manifest)
    log("finalize complete selected_diseases={} no_eligible={} repeats={}".format(len(selected), len(no_eligible), actual))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["selection", "evaluation", "degree_null", "finalize"], required=True)
    args = parser.parse_args()
    if os.path.exists(os.path.join(RUN_ROOT, "config", "run_config.json")):
        run_config_path = os.path.join(RUN_ROOT, "config", "run_config.json")
    else:
        os.makedirs(RUN_ROOT, exist_ok=True)
        for folder in ["config", "inputs", "solver_audit", "grid", "perturbations", "randomizations",
                       "summaries", "figure_source", "logs", "scripts", "selected_networks"]:
            os.makedirs(os.path.join(RUN_ROOT, folder), exist_ok=True)
    runtime_path, runtime_source_hash, runtime_hash = runtime_copy_and_patch()
    driver_hash = sha256(__file__)
    driver_stem, driver_ext = os.path.splitext(os.path.basename(__file__))
    driver_copy = os.path.join(RUN_ROOT, "scripts", driver_stem + "_" + driver_hash[:12] + driver_ext)
    if not os.path.isfile(driver_copy):
        shutil.copy2(__file__, driver_copy)
    elif sha256(driver_copy) != driver_hash:
        raise RuntimeError("Hash-addressed archived driver differs from executing driver")
    code_hash = bytes_sha256((driver_hash + runtime_hash).encode("ascii"))
    mod = setup_module(load_runtime(runtime_path), runtime_path)
    source_manifest = make_source_manifest(mod, runtime_source_hash, runtime_hash, driver_hash)
    config = expected_config(runtime_source_hash, runtime_hash, driver_hash, source_manifest)
    if os.path.isfile(os.path.join(RUN_ROOT, "config", "run_config.json")):
        with open(os.path.join(RUN_ROOT, "config", "run_config.json"), "r", encoding="utf-8-sig") as handle:
            old = json.load(handle)
        if old != config:
            raise RuntimeError("Frozen run configuration/input/code identity changed; refusing resume")
    else:
        write_json(os.path.join(RUN_ROOT, "config", "run_config.json"), config)
    config_hash = sha256(os.path.join(RUN_ROOT, "config", "run_config.json"))

    ppi_df, node_df, disease_data = mod.load_inputs()
    mod.audit_only(ppi_df, node_df, disease_data)
    if set(disease_data) != set(DISEASES):
        raise RuntimeError("Frozen disease input identity incomplete: {}".format(sorted(set(DISEASES) - set(disease_data))))
    mod.ACTIVE_DISEASES[:] = list(DISEASES)
    metrics = verify_scan_summary(mod)
    eligibility, candidate_files = qualify_base_networks(mod, metrics, ppi_df, disease_data)
    # Reconstruct and verify all 1664 frozen solver arrays before any new PCSF call.
    array_audit = audit_frozen_solver_arrays(mod, ppi_df, disease_data, eligibility)
    if len(array_audit) != 1664 or not all(r.get("status") == "verified" and r.get("effective_cost_formula_verified") and r.get("prize_b_applied_once_verified") for r in array_audit.to_dict("records")):
        raise RuntimeError("Full compact-grid solver-array/effective-cost audit did not pass")
    eligibility = mod.read_table(os.path.join(RUN_ROOT, "summaries", "parameter_qualification_all_104.csv"))
    mod.FROZEN_BASE_RESULT_PATHS = {
        (str(row["Disease"]), str(row["parameter_id"])):
            os.path.dirname(str(row["base_run_manifest_path"]))
        for _, row in eligibility.iterrows()
    }
    if len(mod.FROZEN_BASE_RESULT_PATHS) != 1664:
        raise RuntimeError("Verified frozen-base path registry does not cover the complete 1664-row scan")
    eligible = eligibility[eligibility["eligible"].astype(str).str.lower().isin(["true", "1"])].copy()
    active_candidates = sorted(eligible["Disease"].astype(str).unique().tolist())
    if not set(active_candidates).issubset(set(DISEASES)):
        raise RuntimeError("Eligible base networks contain diseases outside the frozen list: {}".format(active_candidates))
    if not active_candidates:
        raise RuntimeError("No disease has an eligible 30-50-node connected tree; qualification table and no-candidate register were saved")
    seeds = prepare_seed_registry(eligible, active_candidates)
    templates = {g: mod.build_graph(active_candidates[0], g, ppi_df, disease_data) for g in G_VALUES}
    for g in G_VALUES:
        # Freeze the formula-built vector used by patched perturbation calls;
        # build_graph already checks OI's costs against this approved formula.
        templates[g]["frozen_effective_costs"] = np.asarray(
            templates[g]["base_costs"] + templates[g]["edge_penalties"], dtype=np.float64).copy()
    bins_by_disease = create_degree_bins(mod, eligible, templates, disease_data)
    input_manifest_path = os.path.join(RUN_ROOT, "inputs", "input_sha256_manifest.csv")
    input_manifest_hash = sha256(input_manifest_path)
    candidate_hash = sha256(os.path.join(RUN_ROOT, "summaries", "parameter_qualification_all_104.csv"))
    selected, lock, lock_hash = None, None, None
    if args.phase != "selection":
        selected, lock, lock_hash = load_locked_parameters(mod, config_hash)
        if set(selected["Disease"].astype(str)) != set(active_candidates):
            raise RuntimeError("Locked selected disease set differs from current candidate eligibility")
        if lock.get("input_manifest_sha256") != input_manifest_hash:
            raise RuntimeError("Selected lock input hash differs from current disease/PPI input manifest")
    if args.phase == "degree_null":
        evaluation_progress_path = os.path.join(RUN_ROOT, "logs", "progress_evaluation.json")
        evaluation_summary_path = os.path.join(RUN_ROOT, "summaries", "evaluation_network_distribution_by_disease.csv")
        if not os.path.isfile(evaluation_progress_path) or not os.path.isfile(evaluation_summary_path):
            raise RuntimeError("Degree-stratified null phase requires completed cost evaluation artifacts")
        with open(evaluation_progress_path, "r", encoding="utf-8-sig") as handle:
            evaluation_progress = json.load(handle)
        expected_evaluation_records = len(selected) * EVALUATION_REPEATS
        if (evaluation_progress.get("status") != "complete" or
                int(evaluation_progress.get("completed_records", -1)) != expected_evaluation_records or
                int(evaluation_progress.get("planned_records", -1)) != expected_evaluation_records):
            raise RuntimeError("Degree-stratified null phase cannot start before all locked-parameter cost repeats complete")
        evaluation_summary = mod.read_table(evaluation_summary_path)
        expected_evaluation_keys = set(zip(selected["Disease"].astype(str), selected["parameter_id"].astype(str)))
        observed_evaluation_keys = set(zip(evaluation_summary["Disease"].astype(str),
                                          evaluation_summary["parameter_id"].astype(str)))
        if (len(evaluation_summary) != len(selected) or
                evaluation_summary.duplicated(["Disease", "parameter_id"]).any() or
                observed_evaluation_keys != expected_evaluation_keys or
                not evaluation_summary["repeat_count"].astype(int).eq(EVALUATION_REPEATS).all()):
            raise RuntimeError("Cost-evaluation summary is incomplete or does not match the locked parameter table")
    else:
        existing_lock = os.path.join(RUN_ROOT, "config", "selected_parameters.lock.json")
        if os.path.isfile(existing_lock):
            selected, lock, lock_hash = load_locked_parameters(mod, config_hash)
    if args.phase in ("selection", "evaluation", "degree_null"):
        run_stage(mod, args.phase, selected, eligibility, ppi_df, disease_data, bins_by_disease,
                  templates, config_hash, input_manifest_hash, code_hash, lock_hash, candidate_hash)
        if args.phase == "selection":
            _, selected, lock = summarize_candidates(mod, eligibility, input_manifest_hash, config_hash, code_hash)
            lock_hash = lock["selected_parameters_sha256"]
            log("selection locked diseases={} selected_parameters_sha256={}".format(len(selected), lock_hash))
        elif args.phase in ("evaluation", "degree_null"):
            save_summary_for_stage(mod, args.phase, selected, config_hash, input_manifest_hash, code_hash)
    else:
        if selected is None:
            selected, lock, lock_hash = load_locked_parameters(mod, config_hash)
        # Rebuild stage summaries from checkpointed records before finalizing, for safe resume.
        for stage in ["evaluation", "degree_null"]:
            name = "evaluation_network_distribution_by_disease.csv" if stage == "evaluation" else "degree_null_network_distribution_by_disease.csv"
            if not os.path.isfile(os.path.join(RUN_ROOT, "summaries", name)):
                save_summary_for_stage(mod, stage, selected, config_hash, input_manifest_hash, code_hash)
        finalize(mod, selected, lock, config_hash, input_manifest_hash, code_hash)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        os.makedirs(os.path.join(RUN_ROOT, "logs"), exist_ok=True)
        with open(os.path.join(RUN_ROOT, "logs", "fatal_error.log"), "a", encoding="utf-8") as handle:
            handle.write("\n{}\n{}\n".format(time.strftime("%Y-%m-%dT%H:%M:%S%z"), traceback.format_exc()))
        raise
