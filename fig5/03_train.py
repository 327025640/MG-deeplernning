"""PTG01–04并集：全冻结训练池、10个独立种子、固定解释集；不保存逐细胞SHAP。"""
import argparse
import gc
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path('E:/micro_output/分析审阅中心/待审阅/PTG02_CD_分支疾病与调控候选')
OUT = ROOT / 'ptg14_training_v2'
ADAPTER = ROOT / 'data/pan14_disease/caufinder_adapter'
HELPER = Path('E:/micro_output/fig3/script/ptg02_p14d_07_network_fulltrain_stability_20260914.py')
PRIOR = ROOT / 'data/pan14_disease/resources/NicheNet_human.csv'
SEEDS = list(range(42,52))
COLORS = {'PTG01':'#2166AC','PTG02':'#B2182B','PTG03':'#1B9E77','PTG04':'#7570B3'}
CONFIG = dict(epochs=39, batch_size=128, n_alpha=5, n_beta=10,
              seeds=SEEDS, background_seed=4242, background_n=100,
              training='all frozen train cells', explanation='all frozen train cells',
              input='PTG01-04 union; detected >=10 train cells and nonconstant',
              additivity_check=False)
os.environ.setdefault('OMP_NUM_THREADS','8')
sys.path.insert(0,str(ADAPTER))
import numpy as np
import pandas as pd
import networkx as nx

def sha(path, algorithm='sha256'):
    h=hashlib.new(algorithm)
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''): h.update(b)
    return h.hexdigest()

def dump(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf8')
    tmp.replace(path)

def helper():
    spec=importlib.util.spec_from_file_location('frozen_network_helpers',HELPER)
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod

def inputs():
    record=json.loads((OUT/'data/export_complete.json').read_text(encoding='utf8'))
    cells=pd.read_csv(OUT/'data/cells.csv')
    genes=pd.read_csv(OUT/'data/genes.csv')
    if cells.CellID.duplicated().any() or genes.gene_id.duplicated().any(): raise ValueError('重复ID')
    if len(cells)!=record['n_cells'] or len(genes)!=record['n_genes']: raise ValueError('导出维度变化')
    matrix=OUT/'data/expression_cells_x_genes.f32'
    if matrix.stat().st_size!=len(cells)*len(genes)*4: raise ValueError('矩阵大小不符')
    for path,key in [(matrix,'matrix_md5'),(OUT/'data/cells.csv','cells_md5'),(OUT/'data/genes.csv','genes_md5')]:
        if sha(path,'md5')!=record[key]: raise ValueError('导出矩阵或行列清单指纹不符：'+str(path))
    tr=np.flatnonzero(cells.split.eq('train').to_numpy())
    if len(tr)!=203080: raise ValueError('训练池改变')
    if cells.groupby('SampleID_Global')['split'].nunique().max()!=1: raise ValueError('样本跨split')
    sources=[matrix,OUT/'data/export_complete.json',OUT/'data/cells.csv',OUT/'data/genes.csv',
             HELPER,PRIOR,Path(__file__)] + sorted((ADAPTER/'CauFinder').glob('*.py'))
    contract=dict(config=CONFIG,fingerprints={str(p):sha(p) for p in sources},
                  matrix_md5=record['matrix_md5'])
    ident=hashlib.sha256(json.dumps(contract,sort_keys=True).encode()).hexdigest()
    locked=OUT/'data/contract.json'
    if locked.exists() and json.loads(locked.read_text(encoding='utf8'))['id']!=ident:
        raise RuntimeError('配置、输入或源码改变，不能混入之前种子；另行审查')
    dump(locked,dict(id=ident,**contract))
    X=np.memmap(matrix,dtype='<f4',mode='r',shape=(len(cells),len(genes)))
    return cells,genes,tr,X,ident

def require_memory(byte_count):
    h=helper(); status=h.memory_status()
    # 复用已审查助手中的Windows GlobalMemoryStatusEx，避免新增依赖。
    available_gib=status.get('available_ram_gib')
    if available_gib is None or available_gib*2**30<byte_count:
        raise MemoryError(f'需要至少约{byte_count/1024**3:.1f} GiB可用内存；当前{status}')

def build_network():
    cells,genes,tr,X,ident=inputs()
    target=OUT/'data/network_complete.json'
    if target.exists():
        old=json.loads(target.read_text(encoding='utf8'))
        if old['contract']!=ident or old['edges_sha256']!=sha(OUT/'data/network_edges.csv'): raise ValueError('网络缓存改变')
        return
    h=helper(); prior=pd.read_csv(PRIOR)
    _,_,nodes=h.make_prior_graph(prior,pd.Index(genes.gene_id))
    idx=pd.Index(genes.gene_id).get_indexer(nodes)
    require_memory(len(tr)*len(nodes)*4+len(nodes)**2*80+2*1024**3)
    # 使用已审查的全细胞Spearman实现，保留并列秩和零表达。
    arr=np.asarray(X[np.ix_(tr,idx)],dtype=np.float32)
    work=OUT/'data/network_rank_work.f64'
    if work.exists(): raise RuntimeError('发现中断的秩缓存；检查无在运行任务后由Luna清理该文件')
    import shutil
    if shutil.disk_usage(OUT).free<arr.size*8+2*1024**3: raise OSError('网络秩临时文件磁盘不足')
    corr,audit=h.exact_spearman_by_rank_matrix(arr,nodes,work)
    if audit['zero_variance_genes']: raise ValueError('网络有常量基因，与输入门槛不符')
    graph,edges,stats=h.reconstruct_from_correlation(prior,pd.Index(genes.gene_id),nodes,corr)
    if stats['final_edges_with_undefined_weight']: raise ValueError('网络存在未定义权重')
    edges.to_csv(OUT/'data/network_edges.csv',index=False)
    dump(target,dict(contract=ident,audit=audit,stats=stats,edges_sha256=sha(OUT/'data/network_edges.csv')))

def train(seed):
    import torch
    import shap
    import anndata as ad
    from CauFinder.caufinder_main import CausalFinder,ShapModel
    from CauFinder.utils import set_seed,get_influence_score
    from sklearn.metrics import roc_auc_score,average_precision_score,balanced_accuracy_score
    torch.set_num_threads(8)
    cells,genes,tr,X,ident=inputs()
    run=OUT/f'runs/seed_{seed}'; run.mkdir(parents=True,exist_ok=True)
    complete=run/'complete.json'
    if complete.exists():
        info=json.loads(complete.read_text(encoding='utf8'))
        if info['contract']!=ident: raise ValueError('种子来源不一致')
        for name,fingerprint in info['outputs'].items():
            if sha(run/name)!=fingerprint: raise ValueError('已完成种子产物改变')
        print(f'seed {seed}已完成，跳过',flush=True); return
    require_memory(len(tr)*len(genes)*16+3*1024**3)
    set_seed(seed)
    train_x=np.asarray(X[tr,:],dtype=np.float32)
    labels=cells.iloc[tr].Group3.map({'Neurological':0.,'Tumor':1.}).to_numpy()
    if not np.isfinite(labels).all() or not np.isfinite(train_x).all(): raise ValueError('输入异常')
    adata=ad.AnnData(train_x,obs=pd.DataFrame({'labels':labels},index=cells.iloc[tr].CellID),var=pd.DataFrame(index=genes.gene_id))
    model=CausalFinder(adata=adata,n_latent=25,n_causal=5,n_hidden=128,
        n_layers_encoder=0,n_layers_decoder=0,n_layers_dpd=0,
        dropout_rate_encoder=0.,dropout_rate_decoder=0.,dropout_rate_dpd=0.,
        use_batch_norm='none',use_batch_norm_dpd=True,pdp_linear=True)
    ckpt=run/'model.pt'
    if ckpt.exists():
        saved=torch.load(ckpt,map_location='cpu',weights_only=False)
        if saved['contract']!=ident or saved['seed']!=seed: raise ValueError('checkpoint配置不符')
        model.module.load_state_dict(saved['state_dict'])
    else:
        # 不沿用旧rolling state_dict恢复优化器；完整epoch结束的本轮模型才可复用。
        os.environ.pop('P14D_CKPT',None)
        start=time.time()
        model.train(max_epochs=39,stage_training=True,train_size=1.,validation_size=None,
                    batch_size=128,drop_last=False,n_alpha=5,n_beta=10)
        history=pd.DataFrame(model.history)
        if len(history)!=39 or not np.isfinite(history.to_numpy(dtype=float)).all(): raise ValueError('训练loss不完整/非有限')
        history.to_csv(run/'training_history.csv',index=False)
        tmp=run/'model.partial.pt'
        torch.save(dict(contract=ident,seed=seed,state_dict=model.module.state_dict(),seconds=time.time()-start),tmp)
        tmp.replace(ckpt)
    device=torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model.module.to(device).eval(); wrapped=ShapModel(model.module,'prob').to(device).eval()
    # 每个模型共享同一背景和解释细胞，直接固定背景索引，不依赖训练消耗的随机数。
    set_seed(4242)
    bg=np.random.RandomState(4242).permutation(len(tr))[:100]
    pd.DataFrame({'CellID':cells.iloc[tr[bg]].CellID.to_numpy()}).to_csv(run/'background.csv',index=False)
    background=torch.as_tensor(train_x[bg],dtype=torch.float32,device=device)
    explainer=shap.DeepExplainer(wrapped,background)
    absolute=np.zeros(len(genes)); signed=np.zeros(len(genes)); errors=[]; outputs=[]
    expected=float(np.asarray(explainer.expected_value).reshape(-1)[0])
    for start in range(0,len(tr),128):
        v=torch.as_tensor(train_x[start:start+128],dtype=torch.float32,device=device)
        s=np.asarray(explainer.shap_values(v,check_additivity=False))
        if s.ndim==3 and s.shape[-1]==1: s=s[:,:,0]
        if s.shape!=tuple(v.shape) or not np.isfinite(s).all(): raise ValueError('SHAP维度/数值异常')
        absolute+=np.abs(s).sum(axis=0,dtype=np.float64); signed+=s.sum(axis=0,dtype=np.float64)
        with torch.no_grad(): p=wrapped(v).detach().cpu().numpy().reshape(-1)
        errors.extend(np.abs(p-expected-s.sum(axis=1)).tolist()); outputs.extend(p.tolist())
        if start%12800==0: print(f'seed {seed} SHAP {start}/{len(tr)}',flush=True)
    a=absolute/len(tr)
    if a.sum()<=0: raise ValueError('所有SHAP权重为0')
    weights=genes.copy(); weights['mean_abs_SHAP']=a; weights['weight']=a/a.sum(); weights['mean_signed_SHAP']=signed/len(tr)
    weights.to_csv(run/'weights.csv',index=False)
    dump(run/'shap_diagnostics.json',dict(expected=expected,mean_abs_residual=float(np.mean(errors)),
          median_abs_residual=float(np.median(errors)),p95_abs_residual=float(np.quantile(errors,.95)),
          max_abs_residual=float(np.max(errors)),n_explained=len(tr),background_n=100,
          explanation='all train cells; no claim that disabled additivity check guarantees rank accuracy'))
    metrics=[]
    for split in ['train','val','holdout']:
        ids=np.flatnonzero(cells.split.eq(split).to_numpy()); predictions=[]
        if split=='train': predictions=outputs
        else:
            with torch.no_grad():
                for start in range(0,len(ids),128):
                    z=torch.as_tensor(np.asarray(X[ids[start:start+128]],dtype=np.float32),device=device)
                    predictions.extend(wrapped(z).cpu().numpy().reshape(-1).tolist())
        y=cells.iloc[ids].Group3.map({'Neurological':0,'Tumor':1}).to_numpy()
        p=np.asarray(predictions)
        if not np.isfinite(p).all(): raise ValueError('预测非有限')
        metrics.append(dict(split=split,n=len(ids),AUROC=roc_auc_score(y,p),AUPRC=average_precision_score(y,p),
                            balanced_accuracy=balanced_accuracy_score(y,p>=.5),min_p=float(p.min()),max_p=float(p.max())))
    pd.DataFrame(metrics).to_csv(run/'metrics.csv',index=False)
    del explainer,train_x,adata; gc.collect()
    edge=pd.read_csv(OUT/'data/network_edges.csv')
    graph=nx.from_pandas_edgelist(edge,'source','target',edge_attr='weight',create_using=nx.DiGraph)
    scores=get_influence_score(graph,weights.set_index('gene_id'),lam=1.,power=1.)
    h=helper(); drivers,events,error=h.optimize_candidates(h.load_driver_function(),graph,scores)
    dump(run/'solver.json',events)
    if error or any(e.get('cvxpy_status')!='optimal' for e in events): raise RuntimeError(f'优化器未成功：{error}')
    if drivers.empty: raise RuntimeError('优化器返回空表')
    drivers.rename_axis('gene_id').to_csv(run/'candidates.csv')
    files=['weights.csv','candidates.csv','metrics.csv','shap_diagnostics.json','solver.json','background.csv','model.pt','training_history.csv']
    dump(complete,dict(contract=ident,seed=seed,outputs={name:sha(run/name) for name in files}))

def summarize():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    _,genes,_,_,ident=inputs()
    frames=[]; selected=[]; diagnostics=[]
    for seed in SEEDS:
        run=OUT/f'runs/seed_{seed}'
        info=json.loads((run/'complete.json').read_text(encoding='utf8'))
        if info['contract']!=ident: raise ValueError('不同契约混入')
        for name,fp in info['outputs'].items():
            if sha(run/name)!=fp: raise ValueError('运行产物指纹改变')
        w=pd.read_csv(run/'weights.csv'); w['seed']=seed; frames.append(w)
        d=pd.read_csv(run/'candidates.csv')
        selected.extend(dict(gene_id=g,seed=seed) for g in d.loc[d.is_CauFVS_driver.astype(str).str.lower().isin(['true','1']),'gene_id'])
        diagnostics.append(dict(seed=seed,**json.loads((run/'shap_diagnostics.json').read_text(encoding='utf8'))))
    long=pd.concat(frames,ignore_index=True)
    if not long.groupby('gene_id').size().eq(10).all(): raise ValueError('不足10次')
    freq=pd.DataFrame(selected).groupby('gene_id').size().rename('selected_runs')
    summary=long.groupby('gene_id').weight.agg(['mean','median','std','min','max']).reset_index().merge(genes,on='gene_id',validate='one_to_one')
    summary['selected_runs']=summary.gene_id.map(freq).fillna(0).astype(int)
    summary['selection_frequency']=summary.selected_runs/10
    summary=summary.sort_values(['selected_runs','median','gene_id'],ascending=[False,False,True])
    for folder in ['tables','figures','source_data']: (OUT/folder).mkdir(exist_ok=True)
    long.to_csv(OUT/'source_data/all_genes_weights_10runs.csv.gz',index=False)
    summary.to_csv(OUT/'tables/gene_weight_and_selection_summary.csv',index=False)
    pd.DataFrame(diagnostics).to_csv(OUT/'tables/shap_diagnostics_all_runs.csv',index=False)
    # 主图显示所有至少一次入选的候选；超过30个分独立页，不丢弃候选。
    candidates=summary.loc[summary.selected_runs>0]
    if candidates.empty: raise ValueError('没有候选')
    plt.rcParams.update({'font.family':'Arial','pdf.fonttype':42,'svg.fonttype':'none','font.size':8})
    rng=np.random.default_rng(42)
    for page,start in enumerate(range(0,len(candidates),30),1):
        sub=candidates.iloc[start:start+30]; names=sub.gene_id.tolist()
        values=[long.loc[long.gene_id.eq(g),'weight'].to_numpy() for g in names]
        fig,ax=plt.subplots(figsize=(max(7.2,len(names)*.32),3.6))
        b=ax.boxplot(values,patch_artist=True,showfliers=False,widths=.6)
        for patch,p in zip(b['boxes'],sub.PTG): patch.set_facecolor(COLORS[p])
        for i,v in enumerate(values,1): ax.scatter(i+rng.uniform(-.12,.12,len(v)),v,s=6,c='#333333',zorder=3)
        ax.set_xticks(range(1,len(names)+1),names,rotation=60,ha='right',color='black')
        ax.set_ylabel('Normalized global SHAP weight'); ax.set_ylim(bottom=0)
        ax.grid(axis='y',alpha=.2); ax.set_axisbelow(True)
        ax.spines[['top','right']].set_visible(False)
        ax.legend(handles=[Patch(facecolor=c,label=p) for p,c in COLORS.items()],ncol=4,frameon=False,fontsize=7)
        fig.tight_layout()
        for ext in ['png','pdf','svg']: fig.savefig(OUT/f'figures/PTG14_SHAP_10runs_boxplot_{page:02d}.{ext}',dpi=600)
        plt.close(fig)
        sub.to_csv(OUT/f'source_data/plot_order_{page:02d}.csv',index=False)
        grid=long.pivot(index='gene_id',columns='seed',values='weight').loc[names,SEEDS]
        fig,ax=plt.subplots(figsize=(5.2,max(3,len(names)*.19)))
        im=ax.imshow(grid.to_numpy(),aspect='auto',cmap='Blues',vmin=0)
        ax.set_yticks(range(len(names)),names,color='black')
        ax.set_xticks(range(10),SEEDS); ax.set_xlabel('Training seed')
        fig.colorbar(im,ax=ax,label='Normalized global SHAP weight')
        fig.tight_layout()
        for ext in ['png','pdf','svg']: fig.savefig(OUT/f'figures/PTG14_SHAP_10runs_heatmap_{page:02d}.{ext}',dpi=600)
        plt.close(fig)
        grid.to_csv(OUT/f'source_data/heatmap_{page:02d}.csv')
    report='\n## PTG01–04并集：10次训练结果\n\n'
    report+=f'输入基因{len(genes)}；训练及解释细胞203080；10次独立初始化，共同背景100细胞。\n\n'
    report+='全局权重为训练细胞mean(abs(SHAP))后在全部输入基因内归一化；固定网络上逐模型重新优化候选。\n\n'
    report+=f'至少一次入选候选{len(candidates)}个，10次全部入选{int((candidates.selected_runs==10).sum())}个。\n\n'
    report+='保留原39 epochs、MC 5/10和关闭加性检查适配；逐运行误差见shap_diagnostics_all_runs.csv。高误差不被自动视为解释可靠。\n'
    report+='这是初始化稳定性，不是独立生物学重复或已验证因果。当前输入未增加表达层批次校正。\n'
    report+='\n前10项（基因、PTG、中位权重、入选次数）：\n'
    for r in candidates.head(10).itertuples(): report+=f'\n- {r.gene_id} / {r.PTG} / {r.median:.6g} / {r.selected_runs}/10'
    path=ROOT/'解读.md'; existing=path.read_text(encoding='utf-8-sig') if path.exists() else ''
    marker='\n## PTG01–04并集：10次训练结果'
    end='\n<!-- END_PTG14_TRAINING_V2 -->'
    if marker in existing:
        before,rest=existing.split(marker,1)
        if end not in rest: raise RuntimeError('旧章节缺少结束标记，保留原文请人工合并')
        after=rest.split(end,1)[1]
        existing=before+after
    path.write_text(existing+report+end+'\n',encoding='utf8')

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('stage',choices=['network','seed','summary','all']); p.add_argument('--seed',type=int,choices=SEEDS)
    args=p.parse_args()
    if args.stage=='network': build_network()
    elif args.stage=='seed':
        if args.seed is None: p.error('seed阶段需要--seed')
        build_network(); train(args.seed)
    elif args.stage=='summary': summarize()
    else:
        build_network()
        for seed in SEEDS:
            subprocess.run([sys.executable,str(Path(__file__).resolve()),'seed','--seed',str(seed)],check=True)
        summarize()
