from __future__ import annotations
import csv, hashlib, json, math, os, sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

SHAP_ROOT = Path(r"E:\micro_output\分析审阅中心\待审阅\Fig6_疾病vsCtrl_全基因CauFinder")
OUT_ROOT = Path(r"E:\micro_output\分析审阅中心\待审阅\Fig6_疾病vsCtrl_SHAP_PCSF模块")
PPI_PATH = Path(r"D:\data\microglia_new\input\drug_repurposing_38core\mentor_20260923\Network-1.txt")
CONFIG_PATH = SHAP_ROOT / "config" / "all_disease_training_configurations.json"
OFFICIAL_PATH = OUT_ROOT / "data" / "orgHs_official_symbol_entrez.tsv"
ALIAS_PATH = OUT_ROOT / "data" / "orgHs_alias_entrez.tsv"
EXPECTED_GENES = 17645
STAMP = datetime.now().astimezone().isoformat(timespec="seconds")

for sub in ("config", "data", "tables", "networks", "scripts", "logs"):
    (OUT_ROOT / sub).mkdir(parents=True, exist_ok=True)

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def mtime(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime).astimezone().isoformat(timespec="seconds")

def read_tsv_map(path: Path, key_col: str) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key, eid = row.get(key_col, ""), row.get("ENTREZID", "")
            if key and eid:
                out[key].add(str(eid))
    return out

official = read_tsv_map(OFFICIAL_PATH, "SYMBOL")
alias = read_tsv_map(ALIAS_PATH, "ALIAS")
official_by_entrez: dict[str, set[str]] = defaultdict(set)
for sym, ids in official.items():
    for eid in ids:
        official_by_entrez[eid].add(sym)

def map_symbol(symbol: str) -> tuple[str, str, str, str]:
    ids = official.get(symbol, set())
    if len(ids) == 1:
        eid = next(iter(ids)); return "official_unique", eid, "official_symbol", eid
    if len(ids) > 1:
        return "official_ambiguous", "", "official_symbol", "|".join(sorted(ids))
    ids = alias.get(symbol, set())
    if len(ids) == 1:
        eid = next(iter(ids)); return "alias_unique", eid, "alias", eid
    if len(ids) > 1:
        return "alias_ambiguous", "", "alias", "|".join(sorted(ids))
    return "unmapped", "", "none", ""

manifest_rows: list[dict[str, object]] = []
def add_manifest(kind: str, disease: str, path: Path, row_count: object, validation: str) -> None:
    manifest_rows.append({"input_type":kind,"Disease":disease,"path":str(path),"sha256":sha256(path),
                          "size_bytes":path.stat().st_size,"modified_local":mtime(path),
                          "row_count":row_count,"validation":validation})

for kind, path in (("frozen_disease_config", CONFIG_PATH), ("mentor_ppi", PPI_PATH),
                   ("official_symbol_map", OFFICIAL_PATH), ("alias_map", ALIAS_PATH),
                   ("mapping_resource_manifest", OUT_ROOT / "data" / "orgHs_mapping_resource_manifest.json")):
    if not path.is_file():
        raise FileNotFoundError(path)
    add_manifest(kind, "", path, "", "present")
with CONFIG_PATH.open("r", encoding="utf-8-sig") as fh:
    frozen = json.load(fh)
diseases = frozen["current_diseases"]
if len(diseases) != 16 or len(set(diseases)) != 16:
    raise ValueError(f"冻结配置中的疾病列表不是16个唯一病种: {diseases}")

# Preserve every original PPI row and audit the declared source-count field.
raw_ppi: list[dict[str, str]] = []
source_mismatches = 0
invalid_numbers = 0
with PPI_PATH.open("r", encoding="utf-8-sig", newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    required = {"Protein_A_Entrez_ID","Protein_B_Entrez_ID","data_source.s.","number"}
    if not required.issubset(reader.fieldnames or []):
        raise ValueError(f"PPI字段不匹配: {reader.fieldnames}")
    for rowid, row in enumerate(reader, 1):
        sources_raw = row["data_source.s."] or ""
        source_items = [x.strip() for x in sources_raw.split(",") if x.strip()]
        declared = row["number"].strip()
        try:
            nval = int(float(declared))
        except Exception:
            nval = None
            invalid_numbers += 1
        mismatch = nval is not None and len(set(source_items)) != nval
        source_mismatches += int(mismatch)
        if nval is None:
            disposition = "invalid_number"
        elif nval < 2:
            disposition = "excluded_number_lt_2"
        elif str(row["Protein_A_Entrez_ID"]).strip() == str(row["Protein_B_Entrez_ID"]).strip():
            disposition = "excluded_self_loop"
        else:
            disposition = "eligible_for_dedup"
        raw_ppi.append({"record_id":str(rowid),"A":str(row["Protein_A_Entrez_ID"]).strip(),
                        "B":str(row["Protein_B_Entrez_ID"]).strip(),"source_text":sources_raw,
                        "source_unique_count":str(len(set(source_items))),"number_raw":declared,
                        "number_numeric":"" if nval is None else str(nval),
                        "source_count_matches_number":"" if nval is None else str(not mismatch).lower(),
                        "disposition":disposition})

eligible = [r for r in raw_ppi if r["disposition"] == "eligible_for_dedup"]
edges: dict[tuple[str,str], list[dict[str,str]]] = defaultdict(list)
for r in eligible:
    a,b=sorted((r["A"],r["B"]))
    edges[(a,b)].append(r)

edge_rows=[]
for (a,b), recs in sorted(edges.items(), key=lambda z:(int(z[0][0]) if z[0][0].isdigit() else z[0][0], int(z[0][1]) if z[0][1].isdigit() else z[0][1])):
    source_union=sorted({s.strip() for r in recs for s in r["source_text"].split(",") if s.strip()})
    nums=sorted({r["number_numeric"] for r in recs}, key=lambda x:int(x))
    edge_rows.append({"Protein_A_Entrez_ID":a,"Protein_B_Entrez_ID":b,"cost":1,
                      "data_source_union":";".join(source_union),
                      "raw_source_texts":" || ".join(sorted({r["source_text"] for r in recs})),
                      "number_values":";".join(nums),"raw_record_count":len(recs),
                      "raw_record_ids":";".join(r["record_id"] for r in recs)})

nodes=set()
degree=defaultdict(int)
parent={}
def find(x):
    parent.setdefault(x,x)
    if parent[x]!=x: parent[x]=find(parent[x])
    return parent[x]
def union(a,b):
    ra,rb=find(a),find(b)
    if ra!=rb: parent[rb]=ra
for e in edge_rows:
    a,b=e["Protein_A_Entrez_ID"],e["Protein_B_Entrez_ID"]
    nodes.update((a,b)); degree[a]+=1; degree[b]+=1; union(a,b)
components=defaultdict(list)
for node in nodes: components[find(node)].append(node)
component_ids={}
comp_order=sorted((sorted(v,key=lambda x:(not x.isdigit(),int(x) if x.isdigit() else x)) for v in components.values()),
                  key=lambda v:(-len(v),v[0]))
for i,comp in enumerate(comp_order,1):
    for node in comp: component_ids[node]=f"PPI_C{i:05d}"

with (OUT_ROOT/"tables"/"ppi_edges_ge2.tsv").open("w",encoding="utf-8",newline="") as fh:
    fields=["Protein_A_Entrez_ID","Protein_B_Entrez_ID","cost","data_source_union","raw_source_texts","number_values","raw_record_count","raw_record_ids"]
    w=csv.DictWriter(fh,fields,delimiter="\t");w.writeheader();w.writerows(edge_rows)
with (OUT_ROOT/"tables"/"ppi_filter_audit.csv").open("w",encoding="utf-8-sig",newline="") as fh:
    fields=["record_type","record_id","A","B","source_text","source_unique_count","number_raw","number_numeric","source_count_matches_number","disposition","metric","value"]
    w=csv.DictWriter(fh,fields);w.writeheader()
    for r in raw_ppi: w.writerow({"record_type":"raw_record",**r})
    summaries={"raw_record_count":len(raw_ppi),"number_ge2_record_count":sum(r["number_numeric"] and int(r["number_numeric"])>=2 for r in raw_ppi),
               "source_count_mismatch_count":source_mismatches,"invalid_number_count":invalid_numbers,
               "self_loop_ge2_record_count":sum(r["disposition"]=="excluded_self_loop" for r in raw_ppi),
               "eligible_nonself_raw_record_count":len(eligible),"deduplicated_edge_count":len(edge_rows),
               "duplicate_excess_record_count":len(eligible)-len(edge_rows),"node_count":len(nodes),"component_count":len(comp_order)}
    for k,v in summaries.items(): w.writerow({"record_type":"summary","metric":k,"value":v})
with (OUT_ROOT/"tables"/"ppi_node_audit.csv").open("w",encoding="utf-8-sig",newline="") as fh:
    fields=["Entrez","official_symbols","background_degree","component_id"]
    w=csv.DictWriter(fh,fields);w.writeheader()
    for eid in sorted(nodes,key=lambda x:(not x.isdigit(),int(x) if x.isdigit() else x)):
        w.writerow({"Entrez":eid,"official_symbols":"|".join(sorted(official_by_entrez.get(eid,set()))),
                    "background_degree":degree[eid],"component_id":component_ids[eid]})

input_summary=[]
mapping_audit_path=OUT_ROOT/"tables"/"gene_mapping_audit.csv"
map_fields=["Disease","gene_set","gene_id","mean_abs_SHAP","mean_abs_SHAP_normalized","rank","mapping_status","chosen_Entrez","mapping_source","candidate_Entrez_ids","in_qualified_PPI","PPI_background_degree","PPI_component_id"]
with mapping_audit_path.open("w",encoding="utf-8-sig",newline="") as maf:
    mw=csv.DictWriter(maf,map_fields);mw.writeheader()
    for disease in diseases:
        tpath=SHAP_ROOT/"tables"/f"{disease}_SHAP_target_full_ranking.csv"
        cpath=SHAP_ROOT/"tables"/f"{disease}_SHAP_Ctrl_full_ranking.csv"
        mpath=SHAP_ROOT/"models"/disease/"model_config.json"
        rpath=SHAP_ROOT/"models"/disease/"runtime.json"
        disease_issues=[]
        for p in (tpath,cpath,mpath):
            if not p.is_file(): disease_issues.append(f"missing:{p.name}")
        if disease_issues:
            input_summary.append({"Disease":disease,"status":"blocked_missing_input","issues":";".join(disease_issues)})
            continue
        model=json.loads(mpath.read_text(encoding="utf-8-sig"))
        if model.get("Disease")!=disease or model.get("positive_class")!=disease or model.get("negative_class")!="Ctrl":
            disease_issues.append("model_identity_or_class_mismatch")
        if model.get("input_object_sha256")!=frozen.get("input_object_sha256") or model.get("n_model_genes")!=EXPECTED_GENES:
            disease_issues.append("model_frozen_input_or_feature_count_mismatch")
        target_rows=[];ctrl_rows=[]
        for p,kind,dest,expected_group,expected_analysis in ((tpath,"target",target_rows,f"disease:{disease}","target_disease_train"),
                                                             (cpath,"Ctrl",ctrl_rows,"disease:Ctrl","Ctrl_train")):
            with p.open("r",encoding="utf-8-sig",newline="") as fh:
                reader=csv.DictReader(fh); required={"target_disease","contribution_group","n_cells","gene_id","mean_abs_SHAP","mean_abs_SHAP_normalized","mean_signed_SHAP_descriptive_only","rank","include_in_model"}
                if not required.issubset(reader.fieldnames or []):
                    disease_issues.append(f"{kind}:missing_fields"); continue
                dest.extend(reader)
            if len(dest)!=EXPECTED_GENES: disease_issues.append(f"{kind}:row_count_{len(dest)}")
            seen=set(); ranks=[]; vals=[]
            for row in dest:
                if row["target_disease"]!=disease: disease_issues.append(f"{kind}:Disease_identity_mismatch"); break
                if row["contribution_group"]!=expected_group: disease_issues.append(f"{kind}:contribution_group_mismatch"); break
                if row.get("analysis_group") and row["analysis_group"]!=expected_analysis: disease_issues.append(f"{kind}:analysis_group_mismatch"); break
                g=row["gene_id"]
                if not g or g in seen: disease_issues.append(f"{kind}:duplicate_or_empty_gene_id"); break
                seen.add(g)
                try:
                    v=float(row["mean_abs_SHAP"]); rk=int(row["rank"])
                    if not math.isfinite(v) or v<0: raise ValueError
                except Exception:
                    disease_issues.append(f"{kind}:invalid_mean_abs_or_rank"); break
                ranks.append(rk); vals.append(v)
                if row["include_in_model"].lower() not in ("true","1","yes"): disease_issues.append(f"{kind}:gene_not_in_model_manifest"); break
            if len(ranks)==EXPECTED_GENES and sorted(ranks)!=list(range(1,EXPECTED_GENES+1)):
                disease_issues.append(f"{kind}:rank_not_complete_consecutive")
        target_by_gene={r["gene_id"]:r for r in target_rows}
        ctrl_by_gene={r["gene_id"]:r for r in ctrl_rows}
        if set(target_by_gene)!=set(ctrl_by_gene): disease_issues.append("target_ctrl_gene_set_mismatch")
        if target_rows and ctrl_rows:
            denominator=max(float(r["mean_abs_SHAP"]) for r in target_rows)
            if denominator<=0: disease_issues.append("zero_max_SHAP_no_reward_input")
        else: denominator=0.0
        mapped_target=defaultdict(list); mapped_ctrl=defaultdict(list)
        for gene_set,rows,dest in (("target",target_rows,mapped_target),("Ctrl",ctrl_rows,mapped_ctrl)):
            for row in rows:
                status,eid,source,candidates=map_symbol(row["gene_id"])
                inp=eid in nodes if eid else False
                mw.writerow({"Disease":disease,"gene_set":gene_set,"gene_id":row["gene_id"],"mean_abs_SHAP":row["mean_abs_SHAP"],
                             "mean_abs_SHAP_normalized":row["mean_abs_SHAP_normalized"],"rank":row["rank"],"mapping_status":status,
                             "chosen_Entrez":eid,"mapping_source":source,"candidate_Entrez_ids":candidates,"in_qualified_PPI":str(inp).lower(),
                             "PPI_background_degree":degree.get(eid,""),"PPI_component_id":component_ids.get(eid,"")})
                if eid and inp:
                    rec=dict(row);rec.update(mapping_status=status,mapping_source=source,entrez=eid)
                    dest[eid].append(rec)
        selection_fields=["Disease","gene_id","mean_abs_SHAP","mean_abs_SHAP_normalized","mean_signed_SHAP_descriptive_only","rank","mapping_status","Entrez","mapping_source","in_qualified_PPI","PPI_background_degree","PPI_component_id","prize","selection_audit_status"]
        with (OUT_ROOT/"tables"/f"{disease}_all_gene_selection_audit.csv").open("w",encoding="utf-8-sig",newline="") as af:
            aw=csv.DictWriter(af,selection_fields);aw.writeheader()
            for row in target_rows:
                status,eid,source,candidates=map_symbol(row["gene_id"]); inppi=eid in nodes if eid else False
                raw=float(row["mean_abs_SHAP"]); prize=raw/denominator if denominator>0 and raw>0 else 0.0
                aw.writerow({"Disease":disease,"gene_id":row["gene_id"],"mean_abs_SHAP":row["mean_abs_SHAP"],
                             "mean_abs_SHAP_normalized":row["mean_abs_SHAP_normalized"],"mean_signed_SHAP_descriptive_only":row["mean_signed_SHAP_descriptive_only"],
                             "rank":row["rank"],"mapping_status":status,"Entrez":eid,"mapping_source":source,
                             "in_qualified_PPI":str(inppi).lower(),"PPI_background_degree":degree.get(eid,""),
                             "PPI_component_id":component_ids.get(eid,""),"prize":prize,
                             "selection_audit_status":"eligible_reward" if inppi and raw>0 else "mapped_zero_reward" if inppi else "not_in_PPI_or_unmapped"})
        aggregated=[]
        for eid in sorted(nodes,key=lambda x:(not x.isdigit(),int(x) if x.isdigit() else x)):
            trg=mapped_target.get(eid,[]); ctrl=mapped_ctrl.get(eid,[])
            maxrow=max(trg,key=lambda r:(float(r["mean_abs_SHAP"]),-int(r["rank"]))) if trg else None
            raw=float(maxrow["mean_abs_SHAP"]) if maxrow else None
            state="no_SHAP_evidence" if maxrow is None else "measured_zero_reward" if raw==0 else "positive_model_reward"
            aggregated.append({"Disease":disease,"Entrez":eid,"Symbol":"|".join(sorted({r["gene_id"] for r in trg})) or "|".join(sorted(official_by_entrez.get(eid,set()))),
                               "source_model_symbols":"|".join(sorted({r["gene_id"] for r in trg})),"mean_abs_SHAP":raw if raw is not None else "",
                               "mean_abs_SHAP_normalized":max((float(r["mean_abs_SHAP_normalized"]) for r in trg),default=""),
                               "prize":(raw/denominator if raw is not None and denominator>0 and raw>0 else 0.0),
                               "rank_for_max_SHAP":maxrow["rank"] if maxrow else "","ctrl_mean_abs_SHAP_annotation":max((float(r["mean_abs_SHAP"]) for r in ctrl),default=""),
                               "model_evidence_status":state,"background_degree":degree[eid],"background_component_id":component_ids[eid],
                               "positive_prize_terminal":str(raw is not None and raw>0 and denominator>0).lower()})
        prize_fields=list(aggregated[0].keys()) if aggregated else ["Disease","Entrez","Symbol","source_model_symbols","mean_abs_SHAP","mean_abs_SHAP_normalized","prize","rank_for_max_SHAP","ctrl_mean_abs_SHAP_annotation","model_evidence_status","background_degree","background_component_id","positive_prize_terminal"]
        with (OUT_ROOT/"tables"/f"{disease}_gene_prizes.csv").open("w",encoding="utf-8-sig",newline="") as pf:
            pw=csv.DictWriter(pf,prize_fields);pw.writeheader();pw.writerows(aggregated)
        total_shap=sum(float(r["mean_abs_SHAP"]) for r in target_rows)
        retained_shap=sum(float(r["mean_abs_SHAP"]) for eid,recs in mapped_target.items() if eid in nodes for r in [max(recs,key=lambda x:float(x["mean_abs_SHAP"]))])
        mapped_count=sum(1 for r in target_rows if (lambda z:z[0] in ("official_unique","alias_unique") and z[1] in nodes)(map_symbol(r["gene_id"])))
        input_summary.append({"Disease":disease,"status":"ready_for_PCSF_input" if not disease_issues else "blocked_input_mismatch",
                              "issue_count":len(disease_issues),"issues":";".join(sorted(set(disease_issues))),
                              "model_genes":len(target_rows),"target_ctrl_gene_set_equal":str(set(target_by_gene)==set(ctrl_by_gene)).lower(),
                              "max_target_SHAP":denominator,"target_SHAP_total":total_shap,"mapped_model_genes_in_PPI":mapped_count,
                              "mapped_gene_fraction":mapped_count/len(target_rows) if target_rows else "",
                              "retained_raw_SHAP_by_Entrez_max":retained_shap,
                              "retained_SHAP_mass_fraction":retained_shap/total_shap if total_shap>0 else "",
                              "unique_reward_Entrez_nodes":len(mapped_target),"qualified_PPI_nodes":len(nodes),"positive_reward_terminals":sum(bool(r["positive_prize_terminal"]=="true") for r in aggregated),
                              "gene_prize_file":str(OUT_ROOT/"tables"/f"{disease}_gene_prizes.csv"),"gene_audit_file":str(OUT_ROOT/"tables"/f"{disease}_all_gene_selection_audit.csv")})
        for kind,p,rows in (("SHAP_target",tpath,target_rows),("SHAP_Ctrl",cpath,ctrl_rows)):
            add_manifest(kind,disease,p,len(rows),";".join(disease_issues) if disease_issues else "validated_identity_count_values_rank_features")
        add_manifest("disease_model_config",disease,mpath,"",";".join(disease_issues) if disease_issues else "validated Disease/classes/object/features")
        if rpath.is_file(): add_manifest("disease_runtime",disease,rpath,"","present; cross-checked but not sole completion criterion")

with (OUT_ROOT/"tables"/"input_disease_audit_summary.csv").open("w",encoding="utf-8-sig",newline="") as fh:
    fields=sorted({k for r in input_summary for k in r})
    w=csv.DictWriter(fh,fields);w.writeheader();w.writerows(input_summary)
with (OUT_ROOT/"tables"/"input_manifest.csv").open("w",encoding="utf-8-sig",newline="") as fh:
    fields=["input_type","Disease","path","sha256","size_bytes","modified_local","row_count","validation"]
    w=csv.DictWriter(fh,fields);w.writeheader();w.writerows(manifest_rows)
summary={"prepared_at":STAMP,"disease_count":len(diseases),"diseases":diseases,"ppi_summary":summaries,
         "mapping_resource":"org.Hs.eg.db 3.19.1 (local installed Bioconductor resource)",
         "pcsf_status":"blocked_missing_dependency_authorization","pcsf_solver_run":False,
         "module_analysis_run":False,"parameter_search":False,"robustness_validation":False}
(OUT_ROOT/"logs"/"input_preparation_summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
print(json.dumps(summary,ensure_ascii=False))
print("disease_statuses",[(r["Disease"],r["status"],r.get("issue_count")) for r in input_summary])
