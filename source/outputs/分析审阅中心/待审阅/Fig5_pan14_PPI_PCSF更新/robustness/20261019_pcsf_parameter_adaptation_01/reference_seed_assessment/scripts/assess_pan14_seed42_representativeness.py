#!/usr/bin/env python3
"""Assess seed42 as a reference network against the ten frozen pan14 PCSF runs.

This is a read-only analysis of existing parameter-grid outputs. It does not
invoke PCSF or alter any network. Outputs are written to a dedicated
reference_seed_assessment directory beneath the completed robustness run.
"""
from __future__ import print_function

import csv
import hashlib
import json
import math
import os
import platform
import statistics
import sys
from collections import defaultdict, deque
from datetime import datetime


ROOT = r"E:\micro_output\分析审阅中心\待审阅\Fig5_pan14_PPI_PCSF更新"
RUN_ID = "20261019_pcsf_parameter_adaptation_01"
RUN = os.path.join(ROOT, "robustness", RUN_ID)
OUT = os.path.join(RUN, "reference_seed_assessment")
SEEDS = list(range(42, 52))
REFERENCE_SEED = 42
PARAMETER_ID = "b10_w3_g3"
EXPECTED_NODE_COUNTS = {42: 82, 43: 60, 44: 41, 45: 203, 46: 76,
                        47: 96, 48: 53, 49: 101, 50: 56, 51: 83}
EXPECTED = {"model_branch": "PTG14 pan14", "features": 4843,
            "ppi_nodes": 12512, "ppi_edges": 83065,
            "ppi_components": 34, "mapped_evidence_nodes": 4175}


def now():
    return datetime.now().astimezone().isoformat(timespec="seconds")


def log(message):
    line = "[{}] {}".format(now(), message)
    print(line, flush=True)
    path = os.path.join(OUT, "logs", "execution.log")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a", encoding="utf-8") as handle:
        handle.write(line + "\n")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path, delimiter=","):
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_rows(path, rows, fieldnames, delimiter=","):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    temp = path + ".tmp"
    with open(temp, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames,
                                delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: csv_value(row.get(key)) for key in fieldnames})
    os.replace(temp, path)


def csv_value(value):
    if value is None:
        return ""
    if isinstance(value, float) and not math.isfinite(value):
        return ""
    if isinstance(value, bool):
        return "True" if value else "False"
    return value


def write_json(path, value):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    temp = path + ".tmp"
    with open(temp, "w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2,
                  sort_keys=True, allow_nan=False)
        handle.write("\n")
    os.replace(temp, path)


def as_int(value, label):
    try:
        return int(str(value).strip())
    except Exception:
        raise ValueError("{} is not an integer: {!r}".format(label, value))


def as_float(value, label):
    if value is None or str(value).strip() == "":
        return None
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("{} is not finite: {!r}".format(label, value))
    return number


def as_bool(value, label):
    normalized = str(value).strip().lower()
    if normalized in ("true", "1", "yes"):
        return True
    if normalized in ("false", "0", "no"):
        return False
    raise ValueError("{} is not boolean: {!r}".format(label, value))


def canonical_entrez(value):
    return str(as_int(value, "Entrez"))


def canonical_edge(a, b):
    left = canonical_entrez(a)
    right = canonical_entrez(b)
    if left == right:
        raise ValueError("Self-loop encountered in selected PPI edge")
    return tuple(sorted((left, right), key=lambda item: int(item)))


def jaccard(left, right):
    left = set(left)
    right = set(right)
    union = left | right
    return float(len(left & right)) / len(union) if union else None


def graph_components(nodes, edges):
    node_set = set(nodes)
    adjacency = dict((node, set()) for node in node_set)
    for edge in edges:
        a, b = edge
        if a not in node_set or b not in node_set:
            raise ValueError("Edge endpoint is outside the assessed node set")
        adjacency[a].add(b)
        adjacency[b].add(a)
    remaining = set(node_set)
    components = []
    while remaining:
        start = min(remaining, key=lambda item: int(item))
        queue = deque([start])
        seen = set()
        while queue:
            node = queue.popleft()
            if node in seen:
                continue
            seen.add(node)
            queue.extend(adjacency[node] - seen)
        components.append(seen)
        remaining.difference_update(seen)
    components.sort(key=lambda comp: (-len(comp), min(comp, key=lambda x: int(x))))
    return components, adjacency


def average(values):
    return statistics.mean(values) if values else None


def median(values):
    return statistics.median(values) if values else None


def rank_desc(rows, field, seed_field="model_seed"):
    ordered = sorted(rows, key=lambda row: (-float(row[field]), int(row[seed_field])))
    return {int(row[seed_field]): position for position, row in enumerate(ordered, 1)}


def output_hash_map(manifest):
    result = {}
    for path, digest in manifest.get("outputs", {}).items():
        result[os.path.normcase(os.path.normpath(path))] = str(digest).lower()
    return result


def require_file(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    return path


def load_model_shap(audit_path, seeds, ppi_node_set):
    rows = read_rows(audit_path)
    if len(rows) != EXPECTED["features"] * len(seeds):
        raise ValueError("Full mapping audit row count does not equal 4,843 x 10")
    by_seed = dict((seed, defaultdict(list)) for seed in seeds)
    feature_counts = dict((seed, set()) for seed in seeds)
    full_max_by_seed = dict((seed, 0.0) for seed in seeds)
    for row in rows:
        seed = as_int(row.get("seed"), "mapping seed")
        if seed not in by_seed:
            raise ValueError("Unexpected model seed in mapping audit: {}".format(seed))
        gene_id = str(row.get("gene_id", "")).strip()
        if not gene_id or gene_id in feature_counts[seed]:
            raise ValueError("Blank or duplicate model feature in seed {}".format(seed))
        feature_counts[seed].add(gene_id)
        raw = as_float(row.get("mean_abs_SHAP"), "mean_abs_SHAP")
        if raw is None or raw < 0:
            raise ValueError("Invalid full-model SHAP value for seed {} {}".format(seed, gene_id))
        full_max_by_seed[seed] = max(full_max_by_seed[seed], raw)
        entrez = str(row.get("chosen_Entrez", "")).strip()
        if not entrez or entrez.lower() in ("na", "nan"):
            continue
        entrez = canonical_entrez(entrez)
        if entrez not in ppi_node_set:
            continue
        by_seed[seed][entrez].append((raw, gene_id))
    for seed in seeds:
        if len(feature_counts[seed]) != EXPECTED["features"]:
            raise ValueError("Seed {} does not contain 4,843 unique model genes".format(seed))
    collapsed = {}
    normalized_by_seed = {}
    for seed in seeds:
        collapsed[seed] = {}
        normalized_by_seed[seed] = {}
        if full_max_by_seed[seed] <= 0:
            raise ValueError("Seed {} full-model SHAP maximum is not positive".format(seed))
        for entrez, values in by_seed[seed].items():
            # Match the accepted pipeline: highest gene-level SHAP per Entrez;
            # gene_id is the deterministic tie-break.
            raw, gene_id = sorted(values, key=lambda item: (-item[0], item[1]))[0]
            collapsed[seed][entrez] = {"raw_mean_abs_SHAP": raw,
                                       "model_gene_id": gene_id}
            normalized_by_seed[seed][entrez] = raw / full_max_by_seed[seed]
        if len(collapsed[seed]) != EXPECTED["mapped_evidence_nodes"]:
            raise ValueError("Seed {} mapped evidence node count is {}, expected {}".format(
                seed, len(collapsed[seed]), EXPECTED["mapped_evidence_nodes"]))
    return collapsed, normalized_by_seed


def load_official_symbols(map_path, accepted_hash):
    actual = sha256(map_path)
    if actual.lower() != accepted_hash.lower():
        raise ValueError("Approved official symbol map hash differs from the frozen source manifest")
    mapping = defaultdict(set)
    for row in read_rows(map_path, delimiter="\t"):
        symbol = str(row.get("SYMBOL", "")).strip()
        entrez = str(row.get("ENTREZID", "")).strip()
        if symbol and entrez:
            mapping[canonical_entrez(entrez)].add(symbol)
    return mapping


def sha_array(array):
    """Match the frozen runner's SHA256 over C-contiguous ndarray bytes."""
    import numpy as np
    contiguous = np.ascontiguousarray(array)
    return hashlib.sha256(contiguous.tobytes(order="C")).hexdigest()


def verify_frozen_run_inputs(config, run_state, frozen_dir):
    """Verify copied inputs and exact solver arrays against frozen manifests."""
    expected_config_hash = run_state.get("config_sha256")
    actual_config_hash = sha256(os.path.join(RUN, "config", "run_config.json"))
    if not expected_config_hash or actual_config_hash.lower() != expected_config_hash.lower():
        raise ValueError("run_config.json SHA256 differs from completed run_state")

    input_hashes = config.get("input_hashes", {})
    if not input_hashes:
        raise ValueError("Frozen run_config has no input_hashes")
    for basename, expected_hash in input_hashes.items():
        if basename == "run_manifest.json":
            path = require_file(os.path.join(ROOT, "networks", "seed_51", basename))
        else:
            path = require_file(os.path.join(frozen_dir, basename))
        actual_hash = sha256(path)
        if actual_hash.lower() != str(expected_hash).lower():
            raise ValueError("Frozen input SHA256 mismatch for {}".format(basename))

    frozen_script_hash = str(config.get("execution_script_sha256", ""))
    if len(frozen_script_hash) != 64:
        raise ValueError("Frozen execution script SHA256 is missing or malformed")
    parameter_manifest_config_hash = hashlib.sha256(
        (actual_config_hash + frozen_script_hash).encode("ascii")).hexdigest()
    # The NPZ is the exact per-seed solver input snapshot. Verify its arrays
    # against the hashes in each run manifest before using the exported network.
    import numpy as np
    return np, parameter_manifest_config_hash


def main():
    for folder in ("config", "inputs", "summaries", "logs", "scripts"):
        os.makedirs(os.path.join(OUT, folder), exist_ok=True)
    log("assessment started: run={}, parameter={}, reference seed={}".format(
        RUN_ID, PARAMETER_ID, REFERENCE_SEED))

    config_path = require_file(os.path.join(RUN, "config", "run_config.json"))
    state_path = require_file(os.path.join(RUN, "config", "run_state.json"))
    config = json.load(open(config_path, "r", encoding="utf-8"))
    state = json.load(open(state_path, "r", encoding="utf-8"))
    if config.get("run_id") != RUN_ID or state.get("status") != "complete":
        raise ValueError("Parameter exploration run identity/status is not complete")
    model_config = config.get("model", {})
    if (model_config.get("branch") != EXPECTED["model_branch"] or
            int(model_config.get("features", -1)) != EXPECTED["features"] or
            model_config.get("model_seeds") != SEEDS or
            model_config.get("labels") != {"Neurological": 0, "Tumor": 1}):
        raise ValueError("Frozen model identity, feature count, seeds, or labels differ")
    ppi_config = config.get("ppi", {})
    if (int(ppi_config.get("nodes", -1)) != EXPECTED["ppi_nodes"] or
            int(ppi_config.get("edges", -1)) != EXPECTED["ppi_edges"] or
            int(ppi_config.get("components", -1)) != EXPECTED["ppi_components"] or
            ppi_config.get("filter") != "number>=2" or
            not ppi_config.get("all_components_retained")):
        raise ValueError("Frozen PPI identity or dimensions differ")

    frozen = os.path.join(RUN, "inputs", "frozen")
    np, expected_parameter_config_hash = verify_frozen_run_inputs(config, state, frozen)
    mapping_audit_path = require_file(os.path.join(frozen, "full_gene_mapping_audit.csv"))
    ppi_nodes_path = require_file(os.path.join(frozen, "ppi_nodes_ge2.csv"))
    ppi_edges_path = require_file(os.path.join(frozen, "oi_input_edges.tsv"))
    source_manifest_path = require_file(os.path.join(frozen, "source_input_sha256_manifest.csv"))
    source_chain_path = require_file(os.path.join(frozen, "seed_source_chain_verification.json"))
    source_chain = json.load(open(source_chain_path, "r", encoding="utf-8"))
    if source_chain.get("status") != "verified":
        raise ValueError("Frozen PTG14 source-chain verification is not complete")
    source_contract_signatures = set()
    for seed in SEEDS:
        seed_source = source_chain.get("seeds", {}).get(str(seed), {})
        if (seed_source.get("status") != "verified" or
                as_int(seed_source.get("input_feature_count"), "source feature count") != EXPECTED["features"] or
                as_int(seed_source.get("full_shap_rows"), "source full SHAP rows") != EXPECTED["features"]):
            raise ValueError("Seed {} is not verified as a full 4,843-feature PTG14 source".format(seed))
        source_contract_signatures.add((seed_source.get("contract_id"),
                                        seed_source.get("model_contract_json_sha256"),
                                        seed_source.get("input_feature_order_sha256")))
        source_manifest = require_file(os.path.join(ROOT, "networks", "seed_{}".format(seed), "run_manifest.json"))
        if sha256(source_manifest).lower() != str(seed_source.get("manifest_sha256", "")).lower():
            raise ValueError("Seed {} original PTG14 model manifest SHA256 differs from source-chain verification".format(seed))
    if len(source_contract_signatures) != 1:
        raise ValueError("The ten model seeds do not share one verified model contract and feature order")
    source_rows = read_rows(source_manifest_path)
    source_by_suffix = {os.path.basename(row["path"]).lower(): row for row in source_rows}
    official_map_row = source_by_suffix.get("orghs_official_symbol_entrez.tsv")
    if official_map_row is None:
        raise ValueError("Approved official symbol mapping resource is not recorded")
    official_map_path = require_file(official_map_row["path"])
    official_symbols = load_official_symbols(official_map_path, official_map_row["sha256"])
    ppi_source_row = source_by_suffix.get("network-1.txt")
    if ppi_source_row is None or sha256(require_file(ppi_source_row["path"])).lower() != str(ppi_source_row["sha256"]).lower():
        raise ValueError("Original mentor PPI SHA256 differs from the frozen source-input manifest")

    ppi_node_rows = read_rows(ppi_nodes_path)
    ppi_nodes = set(canonical_entrez(row["Entrez"]) for row in ppi_node_rows)
    ppi_degree = {canonical_entrez(row["Entrez"]): as_int(row["background_degree"], "PPI background degree")
                  for row in ppi_node_rows}
    if len(ppi_nodes) != EXPECTED["ppi_nodes"]:
        raise ValueError("Frozen PPI node list count/uniqueness mismatch")
    ppi_edge_rows = read_rows(ppi_edges_path, delimiter="\t")
    ppi_edges = set(canonical_edge(row["protein1"].replace("ENTREZ:", ""),
                                   row["protein2"].replace("ENTREZ:", ""))
                    for row in ppi_edge_rows)
    if len(ppi_edge_rows) != EXPECTED["ppi_edges"] or len(ppi_edges) != EXPECTED["ppi_edges"]:
        raise ValueError("Frozen PPI edge count, duplicate, or normalization mismatch")
    if any(a not in ppi_nodes or b not in ppi_nodes for a, b in ppi_edges):
        raise ValueError("Frozen PPI edge endpoint is missing from the node order")

    shap_by_seed, normalized_shap_by_seed = load_model_shap(mapping_audit_path, SEEDS, ppi_nodes)
    node_order_rows = read_rows(os.path.join(frozen, "oi_base_node_order.tsv"), delimiter="\t")
    ordered_ppi_nodes = [canonical_entrez(row["Entrez"]) for row in node_order_rows]
    if len(ordered_ppi_nodes) != EXPECTED["ppi_nodes"] or set(ordered_ppi_nodes) != ppi_nodes:
        raise ValueError("Frozen solver node order is not a bijection to approved PPI nodes")
    node_index_by_entrez = {entrez: index for index, entrez in enumerate(ordered_ppi_nodes)}
    edge_array_rows = read_rows(os.path.join(frozen, "oi_base_edge_array.tsv"), delimiter="\t")
    expected_real_edges = np.asarray([[as_int(row["node_u_index"], "edge u"),
                                       as_int(row["node_v_index"], "edge v")]
                                      for row in edge_array_rows], dtype=np.int64)
    frozen_base_costs = np.asarray([as_float(row["base_cost"], "base cost")
                                    for row in edge_array_rows], dtype=np.float64)
    if expected_real_edges.shape != (EXPECTED["ppi_edges"], 2):
        raise ValueError("Frozen indexed PPI edge array has an unexpected shape")
    if (expected_real_edges.size and
            (int(expected_real_edges.min()) < 0 or int(expected_real_edges.max()) >= EXPECTED["ppi_nodes"])):
        raise ValueError("Frozen indexed PPI edge endpoint index is outside the frozen node order")
    indexed_ppi_edges = set(canonical_edge(ordered_ppi_nodes[int(u)], ordered_ppi_nodes[int(v)])
                             for u, v in expected_real_edges)
    if len(indexed_ppi_edges) != len(expected_real_edges) or indexed_ppi_edges != ppi_edges:
        raise ValueError("Frozen indexed PPI edges do not map exactly to the frozen Entrez PPI edge set")
    degree_by_index = np.asarray([ppi_degree[entrez] for entrez in ordered_ppi_nodes], dtype=np.int64)
    edge_u = expected_real_edges[:, 0]
    edge_v = expected_real_edges[:, 1]
    denominator = ((EXPECTED["ppi_nodes"] - degree_by_index[edge_u] - 1) *
                   (EXPECTED["ppi_nodes"] - degree_by_index[edge_v] - 1) +
                   degree_by_index[edge_u] * degree_by_index[edge_v])
    if not np.allclose(frozen_base_costs, 1.0, rtol=0.0, atol=1e-12):
        raise ValueError("Frozen PPI base edge costs are not all 1")
    expected_effective_costs = frozen_base_costs + (10.0 ** 3 * degree_by_index[edge_u] *
                                                     degree_by_index[edge_v] / denominator)
    if not np.all(np.isfinite(expected_effective_costs)) or np.any(expected_effective_costs <= 0):
        raise ValueError("Reconstructed g=3 effective edge cost is invalid")
    frozen_edge_cost_by_key = {}
    for edge_index, (u, v) in enumerate(expected_real_edges):
        key = canonical_edge(ordered_ppi_nodes[int(u)], ordered_ppi_nodes[int(v)])
        if key in frozen_edge_cost_by_key:
            raise ValueError("Frozen indexed PPI has a duplicated undirected edge")
        frozen_edge_cost_by_key[key] = float(expected_effective_costs[edge_index])
    reference_fixed_arrays = None
    network_dir = os.path.join(RUN, "runs", PARAMETER_ID)
    networks = {}
    input_hash_rows = []
    for path in (config_path, state_path, mapping_audit_path, ppi_nodes_path,
                 ppi_edges_path, source_manifest_path, official_map_path,
                 os.path.join(frozen, "full_seed_gene_weights.csv"),
                 os.path.join(frozen, "oi_base_node_order.tsv"),
                 os.path.join(frozen, "oi_base_edge_array.tsv"),
                 os.path.join(frozen, "oi_input_nodes.csv"),
                 os.path.join(frozen, "solver_objective_audit.csv"),
                 os.path.join(frozen, "solver_quality_audit_manifest.json"),
                 source_chain_path,
                 os.path.join(RUN, "summaries", "seed_pairwise_similarity.csv"),
                 os.path.join(RUN, "summaries", "node_selection_frequency.csv"),
                 os.path.join(RUN, "summaries", "network_metrics_by_run.csv")):
        require_file(path)
        input_hash_rows.append({"path": path, "sha256": sha256(path),
                                "size_bytes": os.path.getsize(path), "role": "frozen_config_or_summary"})
    input_hash_rows.append({"path": ppi_source_row["path"], "sha256": sha256(ppi_source_row["path"]),
                            "size_bytes": os.path.getsize(ppi_source_row["path"]),
                            "role": "original_mentor_PPI"})
    for seed in SEEDS:
        source_manifest = os.path.join(ROOT, "networks", "seed_{}".format(seed), "run_manifest.json")
        input_hash_rows.append({"path": source_manifest, "sha256": sha256(source_manifest),
                                "size_bytes": os.path.getsize(source_manifest),
                                "role": "original_PT G14_model_manifest_seed_{}".format(seed)})

    metric_summary = read_rows(os.path.join(RUN, "summaries", "network_metrics_by_run.csv"))
    metric_rows = [row for row in metric_summary if row.get("parameter_id") == PARAMETER_ID]
    if len(metric_rows) != len(SEEDS):
        raise ValueError("Expected ten rows for the selected parameter in network_metrics_by_run")
    metric_by_seed = {as_int(row["model_seed"], "summary model_seed"): row for row in metric_rows}

    structural_differences = []
    for seed in SEEDS:
        folder = os.path.join(network_dir, "seed_{}".format(seed))
        nodes_path = require_file(os.path.join(folder, "selected_nodes.csv"))
        forest_path = require_file(os.path.join(folder, "selected_forest_edges.tsv"))
        induced_path = require_file(os.path.join(folder, "selected_induced_ppi_edges.tsv"))
        manifest_path = require_file(os.path.join(folder, "run_manifest.json"))
        manifest = json.load(open(manifest_path, "r", encoding="utf-8"))
        if (manifest.get("run_id") != RUN_ID or manifest.get("parameter_id") != PARAMETER_ID or
                as_int(manifest.get("model_seed"), "manifest model_seed") != seed or
                [manifest.get(key) for key in ("b", "w", "g")] != [10, 3, 3] or
                manifest.get("status") != "complete" or
                manifest.get("pruning") != "strong" or
                manifest.get("best_feasible_result_source") != "raw_pcst_fast" or
                manifest.get("dummy_mode") != "terminals" or
                as_int(manifest.get("num_clusters"), "manifest num_clusters") != 1 or
                as_int(manifest.get("root_index"), "manifest root_index") != EXPECTED["ppi_nodes"] or
                as_int(manifest.get("root_edge_count"), "manifest root_edge_count") != EXPECTED["mapped_evidence_nodes"] or
                as_float(manifest.get("root_edge_cost"), "manifest root_edge_cost") is None or
                not math.isclose(as_float(manifest.get("root_edge_cost"), "manifest root_edge_cost"), 3.0,
                                 rel_tol=0.0, abs_tol=1e-12)):
            raise ValueError("Seed {} run identity, parameters, status, or best-solution source differs".format(seed))
        if manifest.get("config_sha256") != expected_parameter_config_hash:
            raise ValueError("Seed {} derived parameter-config hash differs from frozen runner formula".format(seed))
        arrays_path = require_file(os.path.join(folder, "solver_input_arrays.npz"))
        array_file_hash = output_hash_map(manifest).get(os.path.normcase(os.path.normpath(arrays_path)))
        if not array_file_hash or sha256(arrays_path).lower() != array_file_hash:
            raise ValueError("Seed {} solver-input NPZ file hash differs from run manifest".format(seed))
        arrays = np.load(arrays_path, allow_pickle=False)
        array_hash_map = manifest.get("array_sha256", {})
        array_fields = {"real_edges": "real_edges", "scaled_prizes": "real_prizes",
                        "effective_costs": "real_costs", "augmented_edges": "augmented_edges",
                        "augmented_costs": "augmented_costs"}
        for manifest_key, npz_key in array_fields.items():
            if npz_key not in arrays.files or sha_array(arrays[npz_key]).lower() != str(array_hash_map.get(manifest_key, "")).lower():
                raise ValueError("Seed {} {} array does not match its recorded SHA256".format(seed, manifest_key))
        if (not np.array_equal(arrays["real_edges"], expected_real_edges) or
                not np.allclose(arrays["real_costs"], expected_effective_costs, rtol=0.0, atol=1e-12)):
            raise ValueError("Seed {} real PPI arrays differ from the frozen indexed graph/costs".format(seed))
        expected_reward_array = np.zeros(EXPECTED["ppi_nodes"], dtype=np.float64)
        for entrez, normalized in normalized_shap_by_seed[seed].items():
            expected_reward_array[node_index_by_entrez[entrez]] = normalized * 10.0
        if (arrays["real_prizes"].shape != expected_reward_array.shape or
                not np.allclose(arrays["real_prizes"], expected_reward_array, rtol=1e-10, atol=1e-12)):
            raise ValueError("Seed {} solver prizes do not reproduce full-feature SHAP normalization times b=10".format(seed))
        if (not np.array_equal(arrays["augmented_edges"][:len(expected_real_edges)], expected_real_edges) or
                not np.allclose(arrays["augmented_costs"][:len(expected_real_edges)], expected_effective_costs,
                                rtol=0.0, atol=1e-12)):
            raise ValueError("Seed {} augmented arrays do not begin with the frozen real PPI graph".format(seed))
        fixed_arrays = (arrays["real_edges"].copy(), arrays["real_costs"].copy(),
                        arrays["augmented_edges"].copy(), arrays["augmented_costs"].copy())
        if reference_fixed_arrays is None:
            reference_fixed_arrays = fixed_arrays
        elif any(not np.array_equal(left, right) for left, right in zip(reference_fixed_arrays, fixed_arrays)):
            raise ValueError("Seed {} fixed graph/cost/root arrays differ from the other model seeds".format(seed))
        root_index = as_int(manifest["root_index"], "manifest root index")
        augmented_edges = arrays["augmented_edges"]
        augmented_prizes = arrays["augmented_prizes"]
        root_edge_costs = arrays["root_edge_costs"]
        terminals = arrays["terminals"]
        expected_terminals = {node_index_by_entrez[entrez]
                              for entrez in normalized_shap_by_seed[seed]}
        if (root_index != len(arrays["real_prizes"]) or
                len(augmented_prizes) != root_index + 1 or
                not np.allclose(augmented_prizes[:root_index], arrays["real_prizes"], rtol=0.0, atol=0.0) or
                not math.isclose(float(augmented_prizes[root_index]), 0.0, abs_tol=1e-12) or
                len(terminals) != EXPECTED["mapped_evidence_nodes"] or
                len(set(int(x) for x in terminals)) != len(terminals) or
                set(int(x) for x in terminals) != expected_terminals or
                len(root_edge_costs) != len(terminals) or
                len(augmented_edges) != len(expected_real_edges) + len(terminals) or
                len(arrays["augmented_costs"]) != len(augmented_edges) or
                any(not math.isclose(float(cost), 3.0, abs_tol=1e-12) for cost in root_edge_costs) or
                not np.allclose(arrays["augmented_costs"][len(expected_real_edges):], 3.0,
                                rtol=0.0, atol=1e-12)):
            raise ValueError("Seed {} virtual-root index, prize, terminals, or costs are inconsistent".format(seed))
        root_edges = [tuple(sorted((int(a), int(b)))) for a, b in augmented_edges
                      if int(a) == root_index or int(b) == root_index]
        expected_root_edges = [tuple(sorted((int(node), root_index))) for node in terminals]
        if sorted(root_edges) != sorted(expected_root_edges):
            raise ValueError("Seed {} augmented graph root edges differ from recorded terminals".format(seed))
        arrays.close()
        output_hashes = output_hash_map(manifest)
        for output_path in (nodes_path, forest_path, induced_path):
            expected_hash = output_hashes.get(os.path.normcase(os.path.normpath(output_path)))
            if expected_hash is None or sha256(output_path).lower() != expected_hash:
                raise ValueError("Seed {} primary artifact hash mismatch: {}".format(seed, output_path))
        nodes_rows = read_rows(nodes_path)
        ids = [canonical_entrez(row["Entrez"]) for row in nodes_rows]
        if len(ids) != len(set(ids)):
            raise ValueError("Duplicate Entrez in seed {} selected nodes".format(seed))
        node_set = set(ids)
        if not node_set.issubset(ppi_nodes):
            raise ValueError("Seed {} selected node is absent from frozen PPI".format(seed))
        evidence_ids = set()
        node_row_by_id = {}
        for row in nodes_rows:
            eid = canonical_entrez(row["Entrez"])
            row["_entrez"] = eid
            evidence = as_bool(row["model_evidence"], "model_evidence")
            if evidence:
                evidence_ids.add(eid)
                if eid not in shap_by_seed[seed]:
                    raise ValueError("Evidence node {} lacks frozen SHAP for seed {}".format(eid, seed))
                expected_raw = shap_by_seed[seed][eid]["raw_mean_abs_SHAP"]
                expected_normalized = normalized_shap_by_seed[seed][eid]
                if not math.isclose(as_float(row.get("raw_mean_abs_SHAP"), "selected node raw SHAP"),
                                    expected_raw, rel_tol=1e-10, abs_tol=1e-13):
                    raise ValueError("Evidence node raw SHAP differs from frozen full-gene mapping for seed {} {}".format(seed, eid))
                if not math.isclose(as_float(row.get("normalized_prize"), "selected normalized prize"),
                                    expected_normalized, rel_tol=1e-10, abs_tol=1e-13):
                    raise ValueError("Evidence node normalized prize differs from frozen full-gene mapping for seed {} {}".format(seed, eid))
                if not math.isclose(as_float(row.get("solver_prize"), "selected solver prize"),
                                    expected_normalized * 10.0, rel_tol=1e-10, abs_tol=1e-12):
                    raise ValueError("Evidence node solver prize differs from normalized reward times b=10 for seed {} {}".format(seed, eid))
            else:
                if as_float(row.get("raw_mean_abs_SHAP"), "connector raw SHAP") is not None:
                    raise ValueError("Connector {} has a fabricated/nonmissing SHAP value".format(eid))
                if (not math.isclose(as_float(row.get("normalized_prize"), "connector prize"), 0.0, abs_tol=1e-12) or
                        not math.isclose(as_float(row.get("solver_prize"), "connector solver prize"), 0.0, abs_tol=1e-12)):
                    raise ValueError("Connector {} does not have zero normalized and solver reward".format(eid))
            if as_int(row.get("background_degree"), "selected background degree") != ppi_degree[eid]:
                raise ValueError("Seed {} selected node {} background degree differs from frozen PPI metadata".format(seed, eid))
            node_row_by_id[eid] = row

        forest_rows = read_rows(forest_path, delimiter="\t")
        forest_edges = set()
        edge_cost_by_key = {}
        for row in forest_rows:
            edge = canonical_edge(row["Entrez_A"], row["Entrez_B"])
            if edge in forest_edges:
                raise ValueError("Duplicate undirected backbone edge for seed {}".format(seed))
            if edge[0] not in node_set or edge[1] not in node_set or edge not in ppi_edges:
                raise ValueError("Seed {} backbone edge is not an induced frozen PPI edge".format(seed))
            forest_edges.add(edge)
            actual_cost = as_float(row.get("effective_cost"), "effective_cost")
            expected_cost = frozen_edge_cost_by_key.get(edge)
            if expected_cost is None or not math.isclose(actual_cost, expected_cost,
                                                         rel_tol=1e-11, abs_tol=1e-12):
                raise ValueError("Seed {} backbone edge {} cost differs from frozen indexed PPI".format(seed, edge))
            edge_cost_by_key[edge] = actual_cost
        components, adjacency = graph_components(node_set, forest_edges)
        isolates = sum(1 for node in node_set if not adjacency[node])
        induced_rows = read_rows(induced_path, delimiter="\t")
        induced_set = set(canonical_edge(row["Entrez_A"], row["Entrez_B"]) for row in induced_rows)
        if len(induced_rows) != len(induced_set):
            raise ValueError("Seed {} induced PPI edge table contains duplicate undirected edges".format(seed))
        expected_induced = set(edge for edge in ppi_edges if edge[0] in node_set and edge[1] in node_set)
        if induced_set != expected_induced:
            raise ValueError("Seed {} induced PPI table differs from the frozen background induced on selected nodes".format(seed))

        evidence_n = len(evidence_ids)
        connector_n = len(node_set) - evidence_n
        connector_fraction = connector_n / float(len(node_set)) if node_set else None
        selected_raw_mass = sum(shap_by_seed[seed][eid]["raw_mean_abs_SHAP"] for eid in evidence_ids)
        total_raw_mass = sum(item["raw_mean_abs_SHAP"] for item in shap_by_seed[seed].values())
        retained_fraction = selected_raw_mass / total_raw_mass if total_raw_mass else None
        actual_node_count = len(node_set)
        checks = {
            "expected_node_count": EXPECTED_NODE_COUNTS[seed],
            "observed_node_count": actual_node_count,
            "expected_edge_count": EXPECTED_NODE_COUNTS[seed] - 1,
            "observed_edge_count": len(forest_edges),
            "one_component": len(components) == 1,
            "zero_isolates": isolates == 0,
            "best_source_raw": manifest.get("best_feasible_result_source") == "raw_pcst_fast",
        }
        if (actual_node_count != EXPECTED_NODE_COUNTS[seed] or
                len(forest_edges) != actual_node_count - 1 or len(components) != 1 or
                isolates != 0 or not checks["best_source_raw"]):
            structural_differences.append({"model_seed": seed, **checks})

        expected_metrics = metric_by_seed.get(seed, {})
        compares = {
            "n_real_nodes": actual_node_count,
            "n_real_skeleton_edges": len(forest_edges),
            "n_model_evidence_nodes": evidence_n,
            "n_zero_reward_PPI_connectors": connector_n,
            "n_components": len(components),
            "n_isolated_nodes": isolates,
            "largest_component_fraction": len(components[0]) / float(actual_node_count),
            "retained_mapped_SHAP_mass_fraction": retained_fraction,
        }
        for key, actual in compares.items():
            if key not in expected_metrics:
                raise ValueError("Summary metric column missing: {}".format(key))
            expected_value = as_float(expected_metrics[key], "summary " + key)
            if expected_value is None or not math.isclose(float(actual), expected_value,
                                                          rel_tol=1e-9, abs_tol=1e-10):
                raise ValueError("Seed {} independent {}={} differs from summary {}".format(
                    seed, key, actual, expected_value))
        networks[seed] = {"nodes": node_set, "evidence_nodes": evidence_ids,
                          "connectors": node_set - evidence_ids, "edges": forest_edges,
                          "adjacency": adjacency, "components": components,
                          "isolates": isolates, "rows": node_row_by_id,
                          "edge_costs": edge_cost_by_key, "manifest": manifest,
                          "induced_edges": induced_set, "retained_fraction": retained_fraction,
                          "connector_fraction": connector_fraction,
                          "selected_raw_mass": selected_raw_mass,
                          "total_raw_mass": total_raw_mass}
        for path, role in ((nodes_path, "selected_nodes"), (forest_path, "selected_forest_edges"),
                           (induced_path, "selected_induced_ppi_edges"), (manifest_path, "run_manifest")):
            input_hash_rows.append({"path": path, "sha256": sha256(path),
                                    "size_bytes": os.path.getsize(path), "role": "seed_{}_{}".format(seed, role)})

    if structural_differences:
        diff_path = os.path.join(OUT, "logs", "preflight_structural_differences.csv")
        write_rows(diff_path, structural_differences,
                   ["model_seed", "expected_node_count", "observed_node_count",
                    "expected_edge_count", "observed_edge_count", "one_component",
                    "zero_isolates", "best_source_raw"])
        log("preflight stopped: seed network structure differs from stated expectations; see {}".format(diff_path))
        raise ValueError("One or more seed networks differ from the stated structural checks")

    # Independent pairwise set calculations, using backbone edges only.
    node_pairs = []
    evidence_pairs = []
    edge_pairs = []
    node_by_pair = {}
    evidence_by_pair = {}
    edge_by_pair = {}
    for i, seed_a in enumerate(SEEDS):
        for seed_b in SEEDS[i + 1:]:
            a = networks[seed_a]
            b = networks[seed_b]
            node_common = len(a["nodes"] & b["nodes"])
            node_union = len(a["nodes"] | b["nodes"])
            evidence_common = len(a["evidence_nodes"] & b["evidence_nodes"])
            evidence_union = len(a["evidence_nodes"] | b["evidence_nodes"])
            edge_common = len(a["edges"] & b["edges"])
            edge_union = len(a["edges"] | b["edges"])
            key = (seed_a, seed_b)
            node_value = jaccard(a["nodes"], b["nodes"])
            evidence_value = jaccard(a["evidence_nodes"], b["evidence_nodes"])
            edge_value = jaccard(a["edges"], b["edges"])
            base = {"seed_a": seed_a, "seed_b": seed_b}
            node_pairs.append(dict(base, intersection_count=node_common, union_count=node_union,
                                   node_jaccard=node_value))
            evidence_pairs.append(dict(base, intersection_count=evidence_common,
                                       union_count=evidence_union,
                                       evidence_node_jaccard=evidence_value))
            edge_pairs.append(dict(base, intersection_count=edge_common, union_count=edge_union,
                                   backbone_edge_jaccard=edge_value))
            node_by_pair[key] = node_value
            evidence_by_pair[key] = evidence_value
            edge_by_pair[key] = edge_value
    if len(node_pairs) != 45:
        raise ValueError("Expected 45 unique pairwise comparisons")

    # Reconcile the independent pairwise values and node support against the
    # original completed parameter summaries, without using them as the source.
    old_pair_rows = [row for row in read_rows(os.path.join(RUN, "summaries", "seed_pairwise_similarity.csv"))
                     if row.get("parameter_id") == PARAMETER_ID]
    old_pair_by_key = {(as_int(row["seed_a"], "summary seed_a"),
                        as_int(row["seed_b"], "summary seed_b")): row for row in old_pair_rows}
    if len(old_pair_by_key) != 45:
        raise ValueError("Existing seed-pairwise summary is missing or duplicates comparisons")
    for key in node_by_pair:
        old = old_pair_by_key[key]
        for column, actual in (("node_jaccard", node_by_pair[key]),
                               ("evidence_node_jaccard", evidence_by_pair[key]),
                               ("skeleton_edge_jaccard", edge_by_pair[key])):
            reported = as_float(old.get(column), "summary " + column)
            if reported is None or not math.isclose(reported, actual, rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError("Independent pairwise {} differs for {}".format(column, key))

    all_nodes = set().union(*[networks[seed]["nodes"] for seed in SEEDS])
    selected_seeds = {node: [seed for seed in SEEDS if node in networks[seed]["nodes"]]
                      for node in all_nodes}
    evidence_seeds = {node: [seed for seed in SEEDS if node in networks[seed]["evidence_nodes"]]
                      for node in all_nodes}
    nonisolated_seeds = {node: [seed for seed in SEEDS
                                if node in networks[seed]["nodes"] and networks[seed]["adjacency"][node]]
                         for node in all_nodes}
    node_frequency_input = [row for row in read_rows(os.path.join(RUN, "summaries", "node_selection_frequency.csv"))
                           if row.get("parameter_id") == PARAMETER_ID]
    freq_input_by_node = {canonical_entrez(row["Entrez"]): row for row in node_frequency_input}
    if len(freq_input_by_node) != len(all_nodes):
        raise ValueError("Existing node-frequency summary union differs from the ten network node union")

    frequency_rows = []
    symbol_map_status = {}
    for node in sorted(all_nodes, key=lambda x: int(x)):
        shap_values = [shap_by_seed[seed][node]["raw_mean_abs_SHAP"]
                       for seed in SEEDS if node in shap_by_seed[seed]]
        evidence_n = len(evidence_seeds[node])
        seed42_row = networks[REFERENCE_SEED]["rows"].get(node)
        model_symbols = str(seed42_row.get("gene_symbols", "")).strip() if seed42_row else ""
        if model_symbols:
            symbol = model_symbols
            symbol_status = "from_frozen_model_mapping_audit"
        elif len(official_symbols.get(node, set())) == 1:
            symbol = sorted(official_symbols[node])[0]
            symbol_status = "unique_official_symbol"
        else:
            symbol = ""
            symbol_status = "ambiguous_or_missing_official_symbol"
        symbol_map_status[node] = (symbol, symbol_status)
        selected = selected_seeds[node]
        evidence = evidence_seeds[node]
        connector = [seed for seed in SEEDS
                     if node in networks[seed]["nodes"] and
                     node not in networks[seed]["evidence_nodes"]]
        nonisolated = nonisolated_seeds[node]
        summary_row = freq_input_by_node[node]
        if (as_int(summary_row["selected_count_observed"], "summary selected count") != len(selected) or
                as_int(summary_row["backbone_node_count_observed"], "summary nonisolated count") != len(nonisolated)):
            raise ValueError("Independently computed node frequency differs for Entrez {}".format(node))
        existing_median = as_float(summary_row.get("raw_mean_abs_SHAP_median_across_model_seeds"),
                                   "summary SHAP median")
        actual_median = median(shap_values)
        if (actual_median is None) != (existing_median is None) or (
                actual_median is not None and not math.isclose(actual_median, existing_median,
                                                                rel_tol=1e-10, abs_tol=1e-13)):
            raise ValueError("Independent SHAP median differs from existing summary for Entrez {}".format(node))
        frequency_rows.append({"Entrez": node, "Symbol": symbol,
            "symbol_mapping_status": symbol_status,
            "selected_count_of_10": len(selected), "selected_frequency": len(selected) / 10.0,
            "selected_seed_list": "|".join(map(str, selected)),
            "model_evidence_count_of_10": evidence_n,
            "connector_count_of_10": len(connector),
            "connector_seed_list": "|".join(map(str, connector)),
            "evidence_seed_list": "|".join(map(str, evidence)),
            "nonisolated_count_of_10": len(nonisolated),
            "nonisolated_seed_list": "|".join(map(str, nonisolated)),
            "majority_ge6_of10": len(selected) >= 6,
            "core_10_of10": len(selected) == 10,
            "raw_mean_abs_SHAP_median_across_10_models": actual_median,
            "raw_mean_abs_SHAP_min_across_10_models": min(shap_values) if shap_values else None,
            "raw_mean_abs_SHAP_max_across_10_models": max(shap_values) if shap_values else None,
            "model_SHAP_seed_count": len(shap_values),
            "role_across_seeds": ("model_evidence_and_connector" if evidence_n and connector else
                                   "model_evidence_only" if evidence_n else
                                   "zero_reward_PPI_connector_only"),
            "background_degree": as_int(networks[REFERENCE_SEED]["rows"][node]["background_degree"],
                                         "seed42 background degree") if seed42_row else
                                next(as_int(row["background_degree"], "background degree")
                                     for row in ppi_node_rows if canonical_entrez(row["Entrez"]) == node)})

    ref_nodes = networks[REFERENCE_SEED]["nodes"]
    ref_node_rows = networks[REFERENCE_SEED]["rows"]
    support_by_node = {row["Entrez"]: row for row in frequency_rows}
    seed42_node_rows = []
    for node in sorted(ref_nodes, key=lambda x: int(x)):
        source = ref_node_rows[node]
        support = support_by_node[node]
        seed42_node_rows.append({"Entrez": node, "Symbol": support["Symbol"],
            "symbol_mapping_status": support["symbol_mapping_status"],
            "seed42_model_evidence": as_bool(source["model_evidence"], "seed42 model evidence"),
            "seed42_node_role": "model_evidence" if as_bool(source["model_evidence"], "seed42 model evidence") else "zero_reward_PPI_connector",
            "seed42_raw_mean_abs_SHAP": as_float(source.get("raw_mean_abs_SHAP"), "seed42 SHAP"),
            "seed42_normalized_prize": as_float(source.get("normalized_prize"), "seed42 normalized prize"),
            "seed42_solver_prize": as_float(source.get("solver_prize"), "seed42 solver prize"),
            "background_degree": as_int(source["background_degree"], "background degree"),
            "selected_count_of_10": support["selected_count_of_10"],
            "selected_frequency": support["selected_frequency"],
            "selected_seed_list": support["selected_seed_list"],
            "nonisolated_count_of_10": support["nonisolated_count_of_10"],
            "majority_ge6_of10": support["majority_ge6_of10"],
            "core_10_of10": support["core_10_of10"],
            "raw_mean_abs_SHAP_median_across_10_models": support["raw_mean_abs_SHAP_median_across_10_models"],
            "raw_mean_abs_SHAP_min_across_10_models": support["raw_mean_abs_SHAP_min_across_10_models"],
            "raw_mean_abs_SHAP_max_across_10_models": support["raw_mean_abs_SHAP_max_across_10_models"],
            "model_SHAP_seed_count": support["model_SHAP_seed_count"]})

    all_backbone_edges = set().union(*[networks[seed]["edges"] for seed in SEEDS])
    seed42_edge_rows = []
    for edge in sorted(networks[REFERENCE_SEED]["edges"],
                       key=lambda pair: (int(pair[0]), int(pair[1]))):
        a, b = edge
        selected_edge_seeds = [seed for seed in SEEDS if edge in networks[seed]["edges"]]
        co_selected_seeds = [seed for seed in SEEDS
                             if a in networks[seed]["nodes"] and b in networks[seed]["nodes"]]
        costs = [networks[seed]["edge_costs"][edge] for seed in SEEDS
                 if edge in networks[seed]["edge_costs"]]
        if not costs or any(not math.isclose(cost, costs[0], rel_tol=1e-12, abs_tol=1e-12)
                            for cost in costs[1:]):
            raise ValueError("Frozen effective cost is missing or inconsistent for edge {}".format(edge))
        seed42_edge_rows.append({"Entrez_A": a, "Symbol_A": support_by_node[a]["Symbol"],
            "Entrez_B": b, "Symbol_B": support_by_node[b]["Symbol"],
            "backbone_selected_count_of_10": len(selected_edge_seeds),
            "backbone_selected_frequency": len(selected_edge_seeds) / 10.0,
            "selected_seed_list": "|".join(map(str, selected_edge_seeds)),
            "both_endpoints_selected_count_of_10": len(co_selected_seeds),
            "both_endpoints_selected_seed_list": "|".join(map(str, co_selected_seeds)),
            "backbone_frequency_conditional_on_both_endpoints":
                len(selected_edge_seeds) / float(len(co_selected_seeds)) if co_selected_seeds else None,
            "effective_cost": costs[0]})

    # Pairwise representativeness and deterministic tie-breaking by seed.
    repr_rows = []
    for seed in SEEDS:
        other = [s for s in SEEDS if s != seed]
        node_vals = [node_by_pair[(min(seed, s), max(seed, s))] for s in other]
        evidence_vals = [evidence_by_pair[(min(seed, s), max(seed, s))] for s in other]
        edge_vals = [edge_by_pair[(min(seed, s), max(seed, s))] for s in other]
        net = networks[seed]
        repr_rows.append({"model_seed": seed,
            "mean_node_jaccard_to_other_9": average(node_vals),
            "median_node_jaccard_to_other_9": median(node_vals),
            "min_node_jaccard_to_other_9": min(node_vals),
            "max_node_jaccard_to_other_9": max(node_vals),
            "mean_evidence_node_jaccard_to_other_9": average(evidence_vals),
            "mean_backbone_edge_jaccard_to_other_9": average(edge_vals),
            "node_count": len(net["nodes"]),
            "evidence_node_count": len(net["evidence_nodes"]),
            "connector_node_count": len(net["connectors"]),
            "connector_fraction": net["connector_fraction"],
            "node_count_rank_desc_ties_seed_ascending": None,
            "connector_fraction_rank_desc_ties_seed_ascending": None,
            "best_feasible_solution_source": net["manifest"]["best_feasible_result_source"]})
    node_ranks = rank_desc(repr_rows, "node_count")
    connector_ranks = rank_desc(repr_rows, "connector_fraction")
    for row in repr_rows:
        seed = row["model_seed"]
        row["node_count_rank_desc_ties_seed_ascending"] = node_ranks[seed]
        row["connector_fraction_rank_desc_ties_seed_ascending"] = connector_ranks[seed]
    repr_rows.sort(key=lambda row: (-row["mean_node_jaccard_to_other_9"], row["model_seed"]))
    for rank, row in enumerate(repr_rows, 1):
        row["medoid_rank_by_mean_node_jaccard"] = rank
    medoid_seed = repr_rows[0]["model_seed"]

    medoid_key = (min(REFERENCE_SEED, medoid_seed), max(REFERENCE_SEED, medoid_seed))
    if medoid_seed == REFERENCE_SEED:
        medoid_node_j = medoid_evidence_j = medoid_edge_j = 1.0
    else:
        medoid_node_j = node_by_pair[medoid_key]
        medoid_evidence_j = evidence_by_pair[medoid_key]
        medoid_edge_j = edge_by_pair[medoid_key]
    medoid_comparison = [{"reference_seed": REFERENCE_SEED, "medoid_seed_by_node_jaccard": medoid_seed,
        "reference_node_jaccard_rank": next(row["medoid_rank_by_mean_node_jaccard"] for row in repr_rows
                                             if row["model_seed"] == REFERENCE_SEED),
        "reference_mean_node_jaccard": next(row["mean_node_jaccard_to_other_9"] for row in repr_rows
                                             if row["model_seed"] == REFERENCE_SEED),
        "medoid_mean_node_jaccard": repr_rows[0]["mean_node_jaccard_to_other_9"],
        "node_jaccard_reference_vs_medoid": medoid_node_j,
        "evidence_node_jaccard_reference_vs_medoid": medoid_evidence_j,
        "backbone_edge_jaccard_reference_vs_medoid": medoid_edge_j,
        "reference_node_count": len(networks[REFERENCE_SEED]["nodes"]),
        "medoid_node_count": len(networks[medoid_seed]["nodes"]),
        "reference_connector_fraction": networks[REFERENCE_SEED]["connector_fraction"],
        "medoid_connector_fraction": networks[medoid_seed]["connector_fraction"],
        "reference_minus_medoid_node_count": len(networks[REFERENCE_SEED]["nodes"]) - len(networks[medoid_seed]["nodes"]),
        "reference_minus_medoid_connector_fraction": networks[REFERENCE_SEED]["connector_fraction"] - networks[medoid_seed]["connector_fraction"]}]

    # Set-level consensus and two distinct connectivity definitions.
    union_nodes = set(all_nodes)
    majority_nodes = {node for node in all_nodes if len(selected_seeds[node]) >= 6}
    core_nodes = {node for node in all_nodes if len(selected_seeds[node]) == 10}
    seed42_majority_nodes = majority_nodes & networks[REFERENCE_SEED]["nodes"]
    seed42_core_nodes = core_nodes & networks[REFERENCE_SEED]["nodes"]
    consensus_sets = [("union_all_seeds", union_nodes),
                      ("majority_ge6_of10", majority_nodes),
                      ("core_10_of10", core_nodes)]
    ppi_connectivity_edges = ppi_edges
    consensus_rows = []
    for set_name, nodes in consensus_sets:
        evidence_count = sum(1 for node in nodes if evidence_seeds[node])
        connector_count = len(nodes) - evidence_count
        backbone_union = {edge for edge in all_backbone_edges if edge[0] in nodes and edge[1] in nodes}
        ppi_induced = {edge for edge in ppi_connectivity_edges if edge[0] in nodes and edge[1] in nodes}
        for edge_basis, edges in (("seed_backbone_union_restricted_to_node_set", backbone_union),
                                  ("frozen_PPI_induced_subgraph", ppi_induced)):
            if not nodes:
                components = []
                isolates = 0
                largest = 0
            else:
                components, adjacency = graph_components(nodes, edges)
                isolates = sum(1 for node in nodes if not adjacency[node])
                largest = len(components[0])
            consensus_rows.append({"node_set": set_name, "node_count": len(nodes),
                "model_evidence_node_count": evidence_count, "connector_node_count": connector_count,
                "edge_basis": edge_basis, "edge_count": len(edges),
                "component_count": len(components), "isolated_node_count": isolates,
                "largest_component_nodes": largest,
                "largest_component_fraction": largest / float(len(nodes)) if nodes else None,
                "nodes_in_non_singleton_components": sum(len(comp) for comp in components if len(comp) > 1)})

    # The full-background induced graph must be based on PPI edges, never on
    # the displayed backbone edge set.
    majority_ppi_row = next(row for row in consensus_rows
        if row["node_set"] == "majority_ge6_of10" and row["edge_basis"] == "frozen_PPI_induced_subgraph")

    # Cross-check the supplied b10/w3/g3 node frequency table for all nodes.
    if len(frequency_rows) != len(all_nodes):
        raise ValueError("Output frequency table does not cover the full selected-node union")

    # Persist all outputs only after the input and independent reconciliation pass.
    write_rows(os.path.join(OUT, "inputs", "input_sha256_manifest.csv"), input_hash_rows,
               ["path", "sha256", "size_bytes", "role"])
    write_rows(os.path.join(OUT, "summaries", "ten_seed_network_summary.csv"),
               [{"model_seed": seed, "node_count": len(networks[seed]["nodes"]),
                 "model_evidence_node_count": len(networks[seed]["evidence_nodes"]),
                 "zero_reward_connector_node_count": len(networks[seed]["connectors"]),
                 "zero_reward_connector_fraction": networks[seed]["connector_fraction"],
                 "backbone_edge_count": len(networks[seed]["edges"]),
                 "component_count": len(networks[seed]["components"]),
                 "isolated_node_count": networks[seed]["isolates"],
                 "largest_component_node_count": len(networks[seed]["components"][0]),
                 "largest_component_fraction": len(networks[seed]["components"][0]) / float(len(networks[seed]["nodes"])),
                 "retained_mapped_model_SHAP_mass_fraction": networks[seed]["retained_fraction"],
                 "best_feasible_solution_source": networks[seed]["manifest"]["best_feasible_result_source"],
                 "selected_evidence_SHAP_mass": networks[seed]["selected_raw_mass"],
                 "total_mapped_evidence_SHAP_mass": networks[seed]["total_raw_mass"]}
                for seed in SEEDS],
               ["model_seed", "node_count", "model_evidence_node_count", "zero_reward_connector_node_count",
                "zero_reward_connector_fraction", "backbone_edge_count", "component_count",
                "isolated_node_count", "largest_component_node_count", "largest_component_fraction",
                "retained_mapped_model_SHAP_mass_fraction", "best_feasible_solution_source",
                "selected_evidence_SHAP_mass", "total_mapped_evidence_SHAP_mass"])
    write_rows(os.path.join(OUT, "summaries", "pairwise_node_similarity.csv"), node_pairs,
               ["seed_a", "seed_b", "intersection_count", "union_count", "node_jaccard"])
    write_rows(os.path.join(OUT, "summaries", "pairwise_evidence_node_similarity.csv"), evidence_pairs,
               ["seed_a", "seed_b", "intersection_count", "union_count", "evidence_node_jaccard"])
    write_rows(os.path.join(OUT, "summaries", "pairwise_backbone_edge_similarity.csv"), edge_pairs,
               ["seed_a", "seed_b", "intersection_count", "union_count", "backbone_edge_jaccard"])
    write_rows(os.path.join(OUT, "summaries", "seed_representativeness.csv"), repr_rows,
               ["model_seed", "mean_node_jaccard_to_other_9", "median_node_jaccard_to_other_9",
                "min_node_jaccard_to_other_9", "max_node_jaccard_to_other_9",
                "mean_evidence_node_jaccard_to_other_9", "mean_backbone_edge_jaccard_to_other_9",
                "node_count", "evidence_node_count", "connector_node_count", "connector_fraction",
                "node_count_rank_desc_ties_seed_ascending", "connector_fraction_rank_desc_ties_seed_ascending",
                "medoid_rank_by_mean_node_jaccard", "best_feasible_solution_source"])
    write_rows(os.path.join(OUT, "summaries", "seed42_vs_medoid_comparison.csv"), medoid_comparison,
               list(medoid_comparison[0].keys()))
    node_fields = ["Entrez", "Symbol", "symbol_mapping_status", "seed42_model_evidence", "seed42_node_role",
        "seed42_raw_mean_abs_SHAP", "seed42_normalized_prize", "seed42_solver_prize", "background_degree",
        "selected_count_of_10", "selected_frequency", "selected_seed_list", "nonisolated_count_of_10",
        "majority_ge6_of10", "core_10_of10", "raw_mean_abs_SHAP_median_across_10_models",
        "raw_mean_abs_SHAP_min_across_10_models", "raw_mean_abs_SHAP_max_across_10_models",
        "model_SHAP_seed_count"]
    write_rows(os.path.join(OUT, "summaries", "seed42_nodes_with_seed_support.csv"), seed42_node_rows, node_fields)
    write_rows(os.path.join(OUT, "summaries", "seed42_backbone_edges_with_seed_support.csv"), seed42_edge_rows,
               ["Entrez_A", "Symbol_A", "Entrez_B", "Symbol_B", "backbone_selected_count_of_10",
                "backbone_selected_frequency", "selected_seed_list", "both_endpoints_selected_count_of_10",
                "both_endpoints_selected_seed_list", "backbone_frequency_conditional_on_both_endpoints",
                "effective_cost"])
    write_rows(os.path.join(OUT, "summaries", "all_seed_node_frequency.csv"), frequency_rows,
               ["Entrez", "Symbol", "symbol_mapping_status", "selected_count_of_10", "selected_frequency",
                "selected_seed_list", "model_evidence_count_of_10", "connector_count_of_10",
                "connector_seed_list", "evidence_seed_list",
                "nonisolated_count_of_10", "nonisolated_seed_list", "majority_ge6_of10", "core_10_of10",
                "raw_mean_abs_SHAP_median_across_10_models", "raw_mean_abs_SHAP_min_across_10_models",
                "raw_mean_abs_SHAP_max_across_10_models", "model_SHAP_seed_count", "role_across_seeds",
                "background_degree"])
    write_rows(os.path.join(OUT, "summaries", "consensus_set_connectivity.csv"), consensus_rows,
               ["node_set", "node_count", "model_evidence_node_count", "connector_node_count", "edge_basis",
                "edge_count", "component_count", "isolated_node_count", "largest_component_nodes",
                "largest_component_fraction", "nodes_in_non_singleton_components"])

    outputs = []
    for subdir in ("inputs", "summaries"):
        for name in sorted(os.listdir(os.path.join(OUT, subdir))):
            path = os.path.join(OUT, subdir, name)
            if os.path.isfile(path):
                outputs.append({"path": path, "sha256": sha256(path), "size_bytes": os.path.getsize(path)})
    assessment = {
        "status": "complete", "created_local": now(), "run_id": RUN_ID,
        "parameter_id": PARAMETER_ID, "b": 10, "w": 3, "g": 3,
        "reference_seed": REFERENCE_SEED, "model_seeds": SEEDS,
        "model_identity": model_config, "ppi_identity": ppi_config,
        "solver": config.get("solver", {}),
        "input_checks": {"ten_networks_match_expected_node_counts": True,
            "all_networks_single_connected_backbone_component": True,
            "all_backbone_edges_equal_nodes_minus_one": True,
            "all_best_solution_sources_raw_pcst_fast": True,
            "selected_induced_ppi_reconciled_to_frozen_background": True,
            "pairwise_results_reconciled_to_existing_summary": True,
            "node_frequency_and_shap_medians_reconciled_to_existing_summary": True,
            "network_metrics_reconciled_to_existing_summary": True},
        "representativeness": {"medoid_seed_by_node_jaccard": medoid_seed,
            "seed42_rank": next(row["medoid_rank_by_mean_node_jaccard"] for row in repr_rows
                                 if row["model_seed"] == REFERENCE_SEED),
            "seed42_mean_node_jaccard": next(row["mean_node_jaccard_to_other_9"] for row in repr_rows
                                               if row["model_seed"] == REFERENCE_SEED),
            "majority_node_count": len(majority_nodes), "majority_background_ppi_largest_component":
                majority_ppi_row["largest_component_nodes"],
            "majority_model_evidence_count": sum(1 for node in majority_nodes if evidence_seeds[node]),
            "majority_connector_count": sum(1 for node in majority_nodes if not evidence_seeds[node]),
            "seed42_majority_node_count": len(seed42_majority_nodes),
            "seed42_majority_model_evidence_count": sum(1 for node in seed42_majority_nodes if evidence_seeds[node]),
            "seed42_majority_connector_count": sum(1 for node in seed42_majority_nodes if not evidence_seeds[node]),
            "core_node_count": len(core_nodes),
            "core_model_evidence_count": sum(1 for node in core_nodes if evidence_seeds[node]),
            "core_connector_count": sum(1 for node in core_nodes if not evidence_seeds[node]),
            "seed42_core_node_count": len(seed42_core_nodes),
            "seed42_core_model_evidence_count": sum(1 for node in seed42_core_nodes if evidence_seeds[node]),
            "seed42_core_connector_count": sum(1 for node in seed42_core_nodes if not evidence_seeds[node])},
        "outputs": outputs,
        "script_sha256": sha256(os.path.abspath(__file__)),
        "software": {"python": sys.version, "platform": platform.platform(),
                     "dependencies": "Python standard library only"},
        "interpretation_limit": "This compares model-seed network selection; it is not a causal network or a decomposition of seed-to-seed variation mechanisms."
    }
    write_json(os.path.join(OUT, "assessment_manifest.json"), assessment)
    log("assessment complete: seed42 rank {}, medoid seed {}, majority nodes {}, PPI-induced majority LCC {}".format(
        assessment["representativeness"]["seed42_rank"], medoid_seed,
        len(majority_nodes), majority_ppi_row["largest_component_nodes"]))
    print(json.dumps({"status": "complete", "output_dir": OUT,
                      "seed42_rank": assessment["representativeness"]["seed42_rank"],
                      "seed42_mean_node_jaccard": assessment["representativeness"]["seed42_mean_node_jaccard"],
                      "medoid_seed": medoid_seed, "majority_nodes": len(majority_nodes),
                      "core_nodes": len(core_nodes), "core_model_evidence": assessment["representativeness"]["core_model_evidence_count"],
                      "core_connectors": assessment["representativeness"]["core_connector_count"],
                      "seed42_majority_nodes": len(seed42_majority_nodes),
                      "seed42_majority_evidence": assessment["representativeness"]["seed42_majority_model_evidence_count"],
                      "seed42_majority_connectors": assessment["representativeness"]["seed42_majority_connector_count"],
                      "seed42_core_nodes": len(seed42_core_nodes),
                      "seed42_core_evidence": assessment["representativeness"]["seed42_core_model_evidence_count"],
                      "seed42_core_connectors": assessment["representativeness"]["seed42_core_connector_count"],
                      "majority_ppi_lcc": majority_ppi_row["largest_component_nodes"]},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        try:
            log("assessment failed: {}: {}".format(type(error).__name__, error))
        except Exception:
            pass
        raise
