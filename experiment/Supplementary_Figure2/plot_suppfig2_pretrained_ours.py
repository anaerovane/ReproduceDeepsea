"""Infer local DeepSEA checkpoints on the gkm held-out rows and draw separate panels."""
from __future__ import annotations
import csv, hashlib, json, os, sys
from pathlib import Path
import numpy as np
import torch
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
from huggingface_hub import hf_hub_download

HERE = Path(__file__).resolve().parent
FUXIAN = Path(os.environ.get("DEEPSEA_FUXIAN_ROOT", str(HERE.parents[1]))).resolve()
OUT = HERE / "results"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(FUXIAN))
from deepsea_models import load_pretrained_model, load_our_model
from model_assets import resolve_ours_checkpoint, resolve_pretrained_predict_checkpoint

TRAINED_FEATURES = OUT / "gkm_auc_results_completed.csv"
TEST_X = OUT / "test_forward_agct_1000x4.bin"
TEST_Y = OUT / "test_forward_labels_rowmajor.bin"
HF_REPO = os.environ.get("DEEPSEA_HF_REPO", "aer0vane/reproduce_deepsea")
HF_RESULTS = "experiment/Supplementary_Figure2/results"
N_TEST, WIDTH, N_FEATURES = 227_512, 1000, 919
BATCH = 256
BASES = np.array(list("AGCT"))
CHECKPOINTS = {
    "pretrained": resolve_pretrained_predict_checkpoint(FUXIAN),
    "ours": resolve_ours_checkpoint(FUXIAN),
}


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def resolve_hf_result(filename: str) -> Path:
    """Use the project cache when available, otherwise fetch the verified Hub artifact."""
    return Path(hf_hub_download(
        repo_id=HF_REPO,
        filename=f"{HF_RESULTS}/{filename}",
        repo_type="model",
    ))


def load_test_index_map(features: list[int]) -> dict[int, np.ndarray]:
    """Load the compact, flat bundle; fall back to an in-progress run layout."""
    if (OUT / "test_indices.npz").exists():
        index_bundle = OUT / "test_indices.npz"
    else:
        try:
            index_bundle = resolve_hf_result("test_indices.npz")
        except Exception:
            index_bundle = None
    if index_bundle is not None:
        with np.load(index_bundle) as bundle:
            return {i: bundle[f"feature_{i}"].astype(np.int32) for i in features}


    result = {}
    for i in features:
        retry = OUT / "retry_incomplete" / "test_indices" / f"feature_{i}.npy"
        path = retry if retry.exists() else OUT / "test_indices" / f"feature_{i}.npy"
        result[i] = np.load(path).astype(np.int32)
    return result


def load_results():
    with TRAINED_FEATURES.open() as f: all_rows=list(csv.DictReader(f))
    rows=[r for r in all_rows if r.get('status')=='ok']
    feat=[int(r['feature']) for r in rows]
    idx_by_feature=load_test_index_map(feat)
    for r in rows:
        i=int(r['feature'])
        if len(idx_by_feature[i]) != int(r['n_test']):
            raise RuntimeError(f'Feature {i}: current test indices no longer match the AUC row')
    return rows,feat,idx_by_feature


def infer(model_name, features, union, xmap, device):
    cache=OUT/("pretrained_same_test_predictions_agct.npy" if model_name=='pretrained' else "ours_same_test_predictions.npy")
    cache_for_read=cache
    index_cache=OUT/'deepsea_inference_test_indices.npy'
    if not index_cache.exists():
        try:
            index_cache=resolve_hf_result(index_cache.name)
        except Exception:
            pass
    if not cache.exists():
        try:
            cache_for_read=resolve_hf_result(cache.name)
        except Exception:
            pass
    if cache_for_read.exists() and index_cache.exists():
        old_idx=np.load(index_cache,mmap_mode='r')
        p=np.load(cache_for_read,mmap_mode='r')
        if np.array_equal(old_idx,union) and p.shape==(len(union),len(features)):
            print(f'{model_name}: reusing verified prediction cache {p.shape}',flush=True)
            return p
    if model_name=='pretrained': model=load_pretrained_model(device)
    else: model=load_our_model(device)
    if model_name=='pretrained':

        for seq in (model.features,model.classifier):
            for j,m in enumerate(seq):
                if isinstance(m,torch.nn.ReLU): seq[j]=torch.nn.Threshold(0,1e-6)
    dest=OUT/f'.{model_name}_pred_tmp.npy'
    pred=np.lib.format.open_memmap(dest,mode='w+',dtype=np.float32,shape=(len(union),len(features)))
    col=np.asarray(features,dtype=np.int64)
    model.eval()
    with torch.inference_mode():
        for start in range(0,len(union),BATCH):
            ids=union[start:start+BATCH]

            raw=np.asarray(xmap[ids],dtype=np.uint8)
            agct=torch.from_numpy(raw.transpose(0,2,1).copy()).to(device=device,dtype=torch.float32)
            rc=raw[:,::-1,:][:,:,[3,2,1,0]]
            rc=torch.from_numpy(rc.transpose(0,2,1).copy()).to(device=device,dtype=torch.float32)

            both=torch.cat([agct.unsqueeze(2),rc.unsqueeze(2)],dim=0)
            values=model(both)
            b=len(ids)
            avg=(values[:b]+values[b:])*0.5
            pred[start:start+b]=avg[:,col].detach().cpu().numpy()
            if start==0 or start+BATCH>=len(union) or (start//BATCH)%100==0:
                pred.flush()
                print(f'{model_name}: {min(start+BATCH,len(union))}/{len(union)} sequences',flush=True)
    pred.flush(); del model
    final=cache
    dest.replace(final)
    return np.load(final,mmap_mode='r')


def model_aucs(pred, rows, features, union, ymap, idx_by_feature):
    position={int(v):i for i,v in enumerate(union)}
    fcol={f:i for i,f in enumerate(features)}
    auc={}
    for r in rows:
        f=int(r['feature']); idx=idx_by_feature[f]
        y=np.asarray(ymap[idx,f],dtype=np.uint8)
        ix=np.fromiter((position[int(v)] for v in idx),dtype=np.int64,count=len(idx))
        auc[f]=float(roc_auc_score(y,np.asarray(pred[ix,fcol[f]],dtype=np.float64)))
    return auc


def draw(model_name, rows, deepsea_auc):
    data=[]
    for r in rows:
        f=int(r['feature'])
        data.append((f,r['predictor'],deepsea_auc[f],float(r['auc_gkm1000']),float(r['auc_gkm300'])))
    xa=np.array([x[4] for x in data]); ya=np.array([x[2] for x in data])
    d1=np.array([x[2] for x in data]); d2=np.array([x[3] for x in data]); d3=np.array([x[4] for x in data])
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
    fig,(a,b)=plt.subplots(1,2,figsize=(12.5,4.7),gridspec_kw={'width_ratios':[1.05,1]})
    a.scatter(xa,ya,s=11,alpha=.62,color='#3d86b8',edgecolors='none',rasterized=True)
    a.plot([.5,1],[.5,1],color='#aaa',lw=.8)
    a.set(xlim=(.48,1.02),ylim=(.48,1.02),xlabel='gkm-SVM (300bp)',ylabel=f'DeepSEA {model_name} (1000bp)')
    a.set_xticks(np.arange(.5,1.01,.1));a.set_yticks(np.arange(.5,1.01,.1))
    groups=[d1,d2,d3]; pos=[1,2,3]
    vp=b.violinplot(groups,positions=pos,widths=.72,showextrema=False)
    for body in vp['bodies']:body.set_facecolor('white');body.set_edgecolor('#bcbcbc');body.set_linewidth(.8);body.set_alpha(1)
    rng=np.random.default_rng(3547)
    sc=None
    for p,v in zip(pos,groups):
        sc=b.scatter(p+rng.uniform(-.17,.17,len(v)),v,c=v,cmap='Reds',vmin=.5,vmax=1,s=7,alpha=.56,edgecolors='none',rasterized=True)
    b.set(xlim=(.45,3.55),ylim=(.48,1.02),ylabel='AUC',xticks=pos,xticklabels=['DeepSEA\n(1000bp)','gkm-SVM\n(1000bp)','gkm-SVM\n(300bp)'])
    b.set_yticks(np.arange(.5,1.01,.1))
    cb=fig.colorbar(sc,ax=b,fraction=.048,pad=.06);cb.set_ticks([.5,.6,.7,.8,.9,1]);cb.set_label('AUC',rotation=0,labelpad=9)
    fig.tight_layout(w_pad=2)
    for ext in ('png','pdf'):fig.savefig(OUT/f'gkm_vs_deepsea_{model_name}.{ext}',dpi=300 if ext=='png' else None,bbox_inches='tight')
    plt.close(fig)
    return data


def main():
    torch.set_num_threads(4)
    device=torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    rows,features,idxmap=load_results()
    union=np.unique(np.concatenate(list(idxmap.values()))).astype(np.int32)
    xmap=np.memmap(TEST_X,dtype=np.uint8,mode='r',shape=(N_TEST,WIDTH,4))
    ymap=np.memmap(TEST_Y,dtype=np.uint8,mode='r',shape=(N_TEST,N_FEATURES))
    np.save(OUT/'deepsea_inference_test_indices.npy',union)
    auc_rows=[]
    metadata={'device':str(device),'n_features':len(features),'n_unique_test_sequences':len(union),'strand_averaging':True,'input_channel_order':'AGCT',
              'test_indices_source':'per-feature saved held-out indices used for gkm AUCs','checkpoints':{}}
    for name,path in CHECKPOINTS.items():
        print(f'Running {name} on {len(union)} shared held-out sequences ({len(features)} TFs)',flush=True)
        pred=infer(name,features,union,xmap,device)
        aucs=model_aucs(pred,rows,features,union,ymap,idxmap)
        gkm_by_feature={int(r['feature']):r for r in rows}
        for f,v in aucs.items():
            g=gkm_by_feature[f]
            auc_rows.append({'model':name,'feature':f,'predictor':g['predictor'],'auc_deepsea':v,
                             'auc_gkm1000':float(g['auc_gkm1000']),'auc_gkm300':float(g['auc_gkm300']),
                             'n_test':int(g['n_test']),'n_test_positive':int(g['n_test_pos'])})
        draw(name,rows,aucs)
        metadata['checkpoints'][name]={'path':str(path),'sha256':sha256(path)}
        del pred
    outcsv=OUT/'pretrained_ours_same_test_auc.csv'
    with outcsv.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['model','feature','predictor','auc_deepsea','auc_gkm1000','auc_gkm300','n_test','n_test_positive']);w.writeheader();w.writerows(auc_rows)
    metadata['output_figures']=['gkm_vs_deepsea_pretrained.png','gkm_vs_deepsea_ours.png']
    (OUT/'pretrained_ours_same_test_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('Wrote both model figures and paired AUC table',flush=True)

if __name__=='__main__':main()
