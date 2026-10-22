#!/usr/bin/env python3
"""Assemble the frozen 104-point Fig6 disease PCSF base grid by verified provenance."""
import csv, hashlib, json, math, os, shutil
import pandas as pd

ROOT = r'E:\micro_output\分析审阅中心\待审阅\Fig6_疾病vsCtrl_SHAP_PCSF模块'
RUN = os.path.join(ROOT, 'robustness', '20261021_disease_compact_le50_parameter_selection_01')
BASE = os.path.join(RUN, 'base_scan')
DISEASES = ['AD','ALS','ALS_FTD','DLB','DLBD_PD','EOAD','FTD','FTD_MND','GBM','HD','LCBM','LGG','MB','MS','PD','PSP']
B = [0.25,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0]
W = [3,6]
G = [2,3,4,5]
RUN_56 = os.path.join(ROOT,'robustness','20261020_disease_size_target_scan_01')
RUN_B2 = os.path.join(ROOT,'robustness','20261020_disease_parameter_adaptation_08')
RUN_B2G4 = os.path.join(RUN,'base_scan_b2g4_addendum')


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def read_json(path):
    with open(path,'r',encoding='utf-8-sig') as f: return json.load(f)


def read_csv(path, **kwargs):
    # Opening the Unicode path in Python avoids pandas 0.23's Windows path issue.
    with open(path, 'r', encoding='utf-8-sig', newline='') as handle:
        return pd.read_csv(handle, **kwargs)


def write_csv(frame, path):
    with open(path, 'w', encoding='utf-8-sig', newline='') as handle:
        frame.to_csv(handle, index=False)


def parameter_id(b, w, g):
    b_text = format(float(b), '.15g').replace('-', 'm').replace('.', 'p')
    return 'b{}_w{}_g{}'.format(b_text, int(w), int(g))


def add_record(rows, disease, pid, origin, status, manifest_path, config_paths, expected_params):
    if not os.path.isfile(manifest_path): raise FileNotFoundError(manifest_path)
    manifest=read_json(manifest_path)
    if manifest.get('status')!='complete' or manifest.get('Disease')!=disease or manifest.get('parameter_id')!=pid:
        raise RuntimeError('Base manifest status or identity mismatch: '+manifest_path)
    b,w,g=map(float,(manifest['b'],manifest['w'],manifest['g']))
    if ((disease, int(w), int(g)) != (expected_params[0], int(expected_params[2]), int(expected_params[3]))
            or not math.isclose(b, float(expected_params[1]), rel_tol=0.0, abs_tol=1e-12)):
        raise RuntimeError('Manifest parameter mismatch: '+manifest_path)
    config_path=None
    for candidate in config_paths:
        if os.path.isfile(candidate) and sha(candidate)==manifest.get('run_config_sha256'):
            config_path=candidate;break
    if config_path is None: raise RuntimeError('No matching source run_config for '+manifest_path)
    config=read_json(config_path)
    configured_script=config.get('execution_script_sha256')
    if configured_script != manifest.get('execution_script_sha256'):
        raise RuntimeError('Source config/script SHA mismatch: '+manifest_path)
    script_path=config.get('execution_script_path')
    if script_path and (not os.path.isfile(script_path) or sha(script_path)!=configured_script):
        raise RuntimeError('Source execution script missing or SHA mismatch: '+str(script_path))
    outputs=manifest.get('outputs',[])
    output_map={os.path.basename(x['path']):x for x in outputs}
    need=['adopted_nodes.csv','adopted_forest_edges.tsv','solver_arrays.npz','candidate_objective_comparison.csv','raw_solver_vertices.csv','raw_solver_edges.csv']
    if not set(need).issubset(output_map): raise RuntimeError('Required base outputs absent from manifest: '+manifest_path)
    for name,item in output_map.items():
        path=item['path']
        if not os.path.isfile(path) or sha(path)!=item.get('sha256') or os.path.getsize(path)!=int(item.get('size_bytes',-1)):
            raise RuntimeError('Base output SHA/size mismatch: '+path)
    metrics=dict(manifest.get('metrics',{}))
    if metrics.get('Disease')!=disease or metrics.get('parameter_id')!=pid:
        raise RuntimeError('Metrics identity mismatch: '+manifest_path)
    if not math.isclose(float(metrics.get('b')),b,rel_tol=0.0,abs_tol=1e-12) or int(metrics.get('w'))!=int(w) or int(metrics.get('g'))!=int(g):
        raise RuntimeError('Metrics parameters mismatch: '+manifest_path)
    row=metrics
    row.update({'result_origin':origin,'status':status,'source_run_id':config.get('run_id',''),
                'source_run_config_path':config_path,'source_run_config_sha256':sha(config_path),
                'source_execution_script_path':script_path or '',
                'source_execution_script_sha256':configured_script,
                'base_run_manifest_path':manifest_path,'base_run_manifest_sha256':sha(manifest_path),
                'base_selected_nodes_path':output_map['adopted_nodes.csv']['path'],
                'base_selected_nodes_sha256':output_map['adopted_nodes.csv']['sha256'],
                'base_selected_edges_path':output_map['adopted_forest_edges.tsv']['path'],
                'base_selected_edges_sha256':output_map['adopted_forest_edges.tsv']['sha256'],
                'solver_arrays_path':output_map['solver_arrays.npz']['path'],
                'solver_arrays_sha256':output_map['solver_arrays.npz']['sha256'],
                'input_manifest_sha256':manifest.get('input_manifest_sha256',''),
                'source_run_status':status})
    rows.append(row)


def main():
    os.makedirs(os.path.join(BASE,'summaries'),exist_ok=True)
    low=os.path.join(BASE,'summaries','parameter_network_metrics.csv')
    low_status=os.path.join(BASE,'summaries','base_run_status.csv')
    low_archive=os.path.join(BASE,'summaries','lowb_only_per_disease_parameter_network_metrics.csv')
    low_status_archive=os.path.join(BASE,'summaries','lowb_only_base_run_status.csv')
    # Preserve the original 640-row low-b summaries before the assembled tables
    # are written to the same canonical filenames. On resume, never recopy the
    # expanded 1664-row tables over these low-b-only snapshots.
    for source,dest in [(low,low_archive),(low_status,low_status_archive)]:
        if not os.path.isfile(dest):
            if not os.path.isfile(source): raise FileNotFoundError('Missing completed low-b base scan: '+source)
            with open(source,'r',encoding='utf-8-sig',newline='') as src, open(dest,'w',encoding='utf-8-sig',newline='') as dst:
                shutil.copyfileobj(src,dst)
    scan_reg=read_csv(os.path.join(RUN_56,'summaries','network_result_registry.csv'),dtype=str,keep_default_na=False)
    scan_status=read_csv(os.path.join(RUN_56,'summaries','base_run_status.csv'),dtype=str,keep_default_na=False)
    b2_status=read_csv(os.path.join(RUN_B2,'summaries','base_run_status.csv'),dtype=str,keep_default_na=False)
    b2g4_metrics=read_csv(os.path.join(RUN_B2G4,'summaries','parameter_network_metrics.csv'),dtype=str,keep_default_na=False)
    b2g4_status=read_csv(os.path.join(RUN_B2G4,'summaries','base_run_status.csv'),dtype=str,keep_default_na=False)
    low_status_df=read_csv(low_status_archive,dtype=str,keep_default_na=False)
    if len(scan_reg)!=896 or len(scan_status)!=896 or len(low_status_df)!=640 or len(b2g4_metrics)!=32 or len(b2g4_status)!=32:
        raise RuntimeError('Expected exact source counts 896 + 96 + 32 + 640; one source grid is incomplete')
    for label, metrics_df, status_df in [('low-b', read_csv(os.path.join(BASE,'summaries','lowb_only_per_disease_parameter_network_metrics.csv'),dtype=str,keep_default_na=False), low_status_df),
                                         ('b2-g4', b2g4_metrics, b2g4_status)]:
        metric_keys=set(zip(metrics_df['Disease'].astype(str),metrics_df['parameter_id'].astype(str)))
        status_keys=set(zip(status_df['Disease'].astype(str),status_df['parameter_id'].astype(str)))
        if (metrics_df.duplicated(['Disease','parameter_id']).any() or
                status_df.duplicated(['Disease','parameter_id']).any() or metric_keys!=status_keys or
                set(status_df['status'])!={'complete'}):
            raise RuntimeError('{} status rows do not exactly match complete metrics keys'.format(label))
    if scan_status['status'].isin(['complete','complete_reused_prior_base']).sum()!=896:
        raise RuntimeError('Historical 56-point base scan contains non-complete status')
    b2=b2_status[b2_status['parameter_id'].str.match(r'^b2_w[36]_g[235]$')]
    if len(b2)!=96 or set(b2['status'])!={'complete'}:
        raise RuntimeError('The 96 exact b=2 prior base results for g=2,3,5 are incomplete')
    config_candidates=[os.path.join(RUN_56,'config','run_config.json'),os.path.join(RUN_B2,'config','run_config.json'),os.path.join(BASE,'config','run_config.json'),os.path.join(RUN_B2G4,'config','run_config.json')]
    rows=[]
    # Existing b=3..6 network outputs: use the accepted result registry and exact source paths.
    for _,r in scan_reg.iterrows():
        disease=str(r['Disease']);pid=str(r['parameter_id']);b=float(r['b']);w=int(r['w']);g=int(r['g'])
        if b not in [3.0,3.5,4.0,4.5,5.0,5.5,6.0]: raise RuntimeError('Unexpected b in accepted 56 grid')
        st=scan_status[(scan_status.Disease==disease)&(scan_status.parameter_id==pid)]
        if len(st)!=1: raise RuntimeError('Duplicate/missing 56-grid status '+disease+' '+pid)
        add_record(rows,disease,pid,'reused_accepted_56_grid',str(st.iloc[0]['status']),str(r['manifest_path']),config_candidates,(disease,b,w,g))
    # Exact b=2 results from the earlier complete grid; these are base outputs only.
    for _,r in b2.iterrows():
        disease=str(r['Disease']);pid=str(r['parameter_id']);w=int(pid.split('_w')[1].split('_')[0]);g=int(pid.split('_g')[1])
        b2dir=os.path.join(RUN_B2,'grid',disease,pid)
        add_record(rows,disease,pid,'reused_verified_b2_base_grid','complete',os.path.join(b2dir,'run_manifest.json'),config_candidates,(disease,2.0,w,g))
    # b=2,g=4 was absent from the older grid (which instead contained g=6); use the exact in-grid addendum.
    if b2g4_metrics.duplicated(['Disease','parameter_id']).any() or b2g4_status.duplicated(['Disease','parameter_id']).any():
        raise RuntimeError('b=2,g=4 addendum metrics/status contains duplicate keys')
    b2g4_status_by={(str(r.Disease),str(r.parameter_id)):str(r.status) for r in b2g4_status.itertuples()}
    for _,r in b2g4_metrics.iterrows():
        disease=str(r['Disease']);pid=str(r['parameter_id']);b=float(r['b']);w=int(r['w']);g=int(r['g'])
        key=(disease,pid)
        if key not in b2g4_status_by or pid!=parameter_id(b,w,g) or not math.isclose(b,2.0,rel_tol=0.0,abs_tol=1e-12) or g!=4:
            raise RuntimeError('b=2,g=4 addendum identity/status mismatch: {} {}'.format(disease,pid))
        manifest_path=os.path.join(RUN_B2G4,'grid',disease,pid,'run_manifest.json')
        add_record(rows,disease,pid,'newly_computed_missing_b2_g4',b2g4_status_by[key],manifest_path,config_candidates,(disease,2.0,w,4))
    # New low-b base networks, excluding b=2 which has a complete hash-verified source.
    low_metrics=read_csv(low_archive,dtype=str,keep_default_na=False)
    low_status_by={(str(r.Disease),str(r.parameter_id)):str(r.status) for r in low_status_df.itertuples()}
    if len(low_metrics)!=640: raise RuntimeError('New low-b metrics must contain exactly 640 rows')
    if low_metrics.duplicated(['Disease','parameter_id']).any() or low_status_df.duplicated(['Disease','parameter_id']).any():
        raise RuntimeError('New low-b metrics or status has duplicate disease/parameter keys')
    low_registry_rows=[]
    for _,r in low_metrics.iterrows():
        disease=str(r['Disease']);pid=str(r['parameter_id']);w=int(r['w']);g=int(r['g']);b=float(r['b'])
        key=(disease,pid)
        if key not in low_status_by: raise RuntimeError('New low-b result absent from status table')
        expected_pid=parameter_id(b,w,g)
        if pid != expected_pid or b not in [0.25,0.5,1.0,1.5,2.5]:
            raise RuntimeError('New low-b parameter identity mismatch: {} {}'.format(disease,pid))
        manifest_path=os.path.join(BASE,'grid',disease,pid,'run_manifest.json')
        add_record(rows,disease,pid,'newly_computed_missing_low_b',low_status_by[key],manifest_path,config_candidates,(disease,b,w,g))
        low_registry_rows.append({'Disease':disease,'parameter_id':pid,'b':b,'w':w,'g':g,
            'manifest_path':manifest_path,'manifest_sha256':sha(manifest_path),
            'selected_nodes_path':os.path.join(BASE,'grid',disease,pid,'adopted_nodes.csv'),
            'selected_edges_path':os.path.join(BASE,'grid',disease,pid,'adopted_forest_edges.tsv'),
            'solver_arrays_path':os.path.join(BASE,'grid',disease,pid,'solver_arrays.npz'),
            'result_origin':'newly_computed_missing_low_b','status':low_status_by[key]})
    write_csv(pd.DataFrame(low_registry_rows),os.path.join(BASE,'summaries','lowb_reconstructed_result_registry.csv'))
    df=pd.DataFrame(rows)
    expected={(d,b,w,g) for d in DISEASES for b in B for w in W for g in G}
    observed={(str(r.Disease),float(r.b),int(r.w),int(r.g)) for r in df.itertuples()}
    if len(df)!=1664 or observed!=expected or df.duplicated(['Disease','parameter_id']).any():
        raise RuntimeError('Assembled grid does not equal the unique frozen 16x104 candidate set')
    input_hashes=set(df['input_manifest_sha256'].astype(str))
    if len(input_hashes)!=1 or '' in input_hashes:
        raise RuntimeError('Reused and newly computed grids do not share one exact input-manifest identity')
    df=df.sort_values(['Disease','b','w','g']).reset_index(drop=True)
    write_csv(df,os.path.join(BASE,'summaries','per_disease_parameter_network_metrics.csv'))
    status=df[['Disease','parameter_id','status','solver_runtime_seconds','result_origin','base_run_manifest_path','base_run_manifest_sha256']].copy()
    write_csv(status,os.path.join(BASE,'summaries','base_run_status.csv'))
    write_csv(df[['Disease','parameter_id','b','w','g','result_origin','source_run_id','source_run_config_path','source_run_config_sha256','base_run_manifest_path','base_run_manifest_sha256','base_selected_nodes_path','base_selected_nodes_sha256','base_selected_edges_path','base_selected_edges_sha256','solver_arrays_path','solver_arrays_sha256','input_manifest_sha256','status']],os.path.join(BASE,'summaries','base_result_registry.csv'))
    # Replace only the base-grid summary descriptor inside this new run; original sources remain untouched.
    descriptor={'run_id':'20261021_disease_compact_le50_parameter_selection_01_base_scan','status':'complete_verified_104_point_base_grid',
                'parameter_count_per_disease':104,'disease_count':16,'rows':1664,'b':B,'w':W,'g':G,
                'reused_accepted_3_to_6_rows':896,'reused_verified_b2_rows':96,'newly_computed_missing_b2_g4_rows':32,
                'newly_computed_missing_low_b_rows':640,
                'input_manifest_sha256':next(iter(input_hashes)),
                'qualification_rule':'30<=real nodes<=50; exactly one connected true PCSF backbone tree; no truncation or manual node/edge changes'}
    with open(os.path.join(BASE,'summaries','size_target_summary_manifest.json'),'w',encoding='utf-8') as f: json.dump(descriptor,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({'rows':len(df),'unique_input_manifest_sha256':list(input_hashes),'origins':df.result_origin.value_counts().to_dict(),'statuses':df.status.value_counts().to_dict()},ensure_ascii=False,indent=2))

if __name__=='__main__': main()
