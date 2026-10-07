#!/usr/bin/env python3
"""Real ref/alt inference; Figure 3 chromatin-only replication on supplied SNP sets.

No HGMD substitution, conservation features, or copied baseline scores.
Default uses ALL rows in supplementary tables 5/6. --random-cap optionally
limits random negatives, retaining every positive and distance-group negative.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
import xgboost as xgb

HERE = Path(__file__).resolve().parent
EXPERIMENT_DIR = HERE.parent
ROOT = EXPERIMENT_DIR.parent.parent
OUTPUT_DIR = EXPERIMENT_DIR / 'outputs'
FIGURE_DIR = EXPERIMENT_DIR / 'figures'
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT))
from deepsea_models import load_our_model, load_pretrained_model
from model_assets import resolve_ours_checkpoint, resolve_pretrained_predict_checkpoint

GROUPS = ['Random negative SNP', '31kbp negative SNP', '6.3kbp negative SNP',
          '710bp negative SNP', '360bp negative SNP']


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def write_metrics_provenance(out):
    source_manifest = out/'baseline_source_manifest.json'
    provenance = {
        'metric_source': 'DeepSEA AUC recalculated from OOF probabilities; CADD, GWAVA, and FunSeq2 AUC recalculated from their per-variant scores',
        'plot_source': 'two separate pretrained and ours figures generated from dynamically computed AUCs; all methods use common complete cases per task/group',
        'sampling': f"random_negative_cap={json.loads((out/'manifest.json').read_text()).get('random_negative_cap')}; positives and distance groups retained",
        'inputs': {
            'DeepSEA pretrained sequence model': sha256(resolve_pretrained_predict_checkpoint(ROOT)),
            'Hugging Face: aer0vane/reproduce_deepsea/training_checkpoints/best_model_FINAL_EPOCH53.pth': sha256(resolve_ours_checkpoint(ROOT)),
            'newdata/41592_2015_BFnmeth3547_MOESM649_ESM.csv': sha256(ROOT/'newdata/41592_2015_BFnmeth3547_MOESM649_ESM.csv'),
            'newdata/41592_2015_BFnmeth3547_MOESM650_ESM.csv': sha256(ROOT/'newdata/41592_2015_BFnmeth3547_MOESM650_ESM.csv'),
        },
        'reference_genome_sha256': {
            p.name:sha256(p) for p in sorted((ROOT/'reference_genome').glob('chr*.fa'))
            if p.name not in {'chrM.fa'}
        },
        'code_sha256': {name:sha256(HERE/name) for name in [
            'figure3_rerun.py','deepsea_models.py','model_assets.py',
            'retrieve_gwava_scores.py','retrieve_funseq2_scores.py']},
        'oof_sha256': {p.name: sha256(p) for p in sorted(out.glob('*_oof.tsv.gz'))},
        'cadd_score_sha256': {p.name: sha256(p) for p in sorted(out.glob('*_CADD_scores.tsv.gz'))},
        'gwava_score_sha256': sha256(out/'gwava_GRCh37_scores.tsv.gz'),
        'funseq2_score_sha256': sha256(out/'funseq2_v2.1.0_scores.tsv.gz'),
        'baseline_source_manifest_sha256': sha256(source_manifest) if source_manifest.exists() else None,
        'auc_summary_sha256': sha256(out/'auc_summary.csv'),
        'comparison_auc_sha256': sha256(out/'pretrained_ours_baselines_auc.csv'),
        'figure_sha256': {p.name: sha256(p) for p in sorted(FIGURE_DIR.glob('figure3_*_separate.*'))},
    }
    (out/'metrics_provenance.json').write_text(json.dumps(provenance, indent=2))


def read_tables(cap):
    tables = {}
    for task, number, positive in [('eqtl', 649, 'eQTL'), ('gwas', 650, 'GWAS Catalog')]:
        path = ROOT / 'newdata' / f'41592_2015_BFnmeth3547_MOESM{number}_ESM.csv'
        df = pd.read_csv(path, comment='#')
        df['positive'] = df.label.eq(positive)
        if cap:
            random = df[df.label.eq(GROUPS[0])]
            df = pd.concat([df[~df.label.eq(GROUPS[0])], random.sample(min(cap, len(random)), random_state=42)])
        tables[task] = df.reset_index(drop=True)
    return tables


def nearest_distances(df):
    result = np.full(len(df), np.nan)
    for chrom, sub in df.groupby('chr'):
        positions = np.sort(sub.loc[sub.positive, 'pos'].unique())
        if not len(positions):
            continue
        query = sub.pos.to_numpy()
        ix = np.searchsorted(positions, query)
        distance = np.minimum(np.abs(query - positions[np.clip(ix, 0, len(positions)-1)]),
                              np.abs(query - positions[np.clip(ix-1, 0, len(positions)-1)]))
        result[sub.index] = distance
    return result


def prepare(tables, out):

    keycols = ['chr', 'pos', 'ref', 'alt']
    variants = pd.concat([d[keycols] for d in tables.values()]).drop_duplicates().reset_index(drop=True)
    variants['variant_id'] = np.arange(len(variants))
    seqpath = out / 'sequences.npy'
    if (out / 'prepared.json').exists():
        return pd.read_csv(out / 'variants.tsv', sep='\t'), np.load(seqpath, mmap_mode='r')
    bases = np.lib.format.open_memmap(seqpath, mode='w+', dtype=np.uint8, shape=(len(variants), 1000))
    lookup = np.full(256, 4, dtype=np.uint8)
    lookup[np.frombuffer(b'AGCT', np.uint8)] = np.arange(4)
    variants['valid'] = False
    variants['exclusion'] = ''
    for chrom, sub in variants.groupby('chr', sort=False):
        path = ROOT / 'reference_genome' / f'{chrom}.fa'
        if not path.exists():
            variants.loc[sub.index, 'exclusion'] = 'missing_chromosome'
            continue
        with open(path, 'rb') as f:
            genome = b''.join(line.strip() for line in f if not line.startswith(b'>')).upper()
        for row in sub.itertuples():
            pos = int(row.pos)-1
            if len(row.ref) != 1 or len(row.alt) != 1 or row.ref not in 'AGCT' or row.alt not in 'AGCT':
                reason = 'not_biallelic_snv'
            elif pos < 499 or pos+501 > len(genome):
                reason = 'window_out_of_bounds'
            elif genome[pos:pos+1] != row.ref.encode():
                reason = 'reference_mismatch'
            else:
                encoded = lookup[np.frombuffer(genome[pos-499:pos+501], dtype=np.uint8)]
                if np.any(encoded == 4):
                    reason = 'ambiguous_window'
                else:
                    bases[row.Index] = encoded
                    variants.loc[row.Index, 'valid'] = True
                    reason = ''
            variants.loc[row.Index, 'exclusion'] = reason
        print(f'Prepared {chrom}: {len(sub):,} variants', flush=True)
    bases.flush()
    variants.to_csv(out / 'variants.tsv', sep='\t', index=False)
    (out / 'prepared.json').write_text(json.dumps({'total':len(variants), 'valid':int(variants.valid.sum()),
                                                'excluded':variants.exclusion.value_counts().to_dict()}, indent=2))
    return variants, bases


def effects(model, encoded, alternate, device):

    ref = torch.nn.functional.one_hot(torch.as_tensor(encoded.astype(np.int64), device=device), 4)
    ref = ref.permute(0, 2, 1).float().unsqueeze(2).contiguous()
    alt = ref.clone()
    alt[:, :, 0, 499] = 0
    alt[torch.arange(len(alt), device=device), torch.as_tensor(alternate, device=device), 0, 499] = 1
    abs_effect, odds_effect = [], []
    for reverse in [False, True]:
        r, a = (ref.flip([1,3]), alt.flip([1,3])) if reverse else (ref, alt)
        predictions = model(torch.cat([r, a]))
        pr, pa = predictions[:len(ref)], predictions[len(ref):]
        abs_effect.append((pr-pa).abs())
        odds_effect.append(torch.logit(pr.clamp(1e-6,1-1e-6))-torch.logit(pa.clamp(1e-6,1-1e-6)))
    return torch.cat([sum(abs_effect)/2, sum(odds_effect)/2], dim=1).cpu().numpy()


def infer(variants, bases, out, name, batch_size, device):
    path = out / f'{name}_effects.npy'
    progress = out / f'{name}_progress.json'
    valid = np.flatnonzero(variants.valid.to_numpy())
    start = json.loads(progress.read_text())['completed'] if progress.exists() else 0
    features = np.lib.format.open_memmap(path, mode='r+' if path.exists() else 'w+',
                                        dtype=np.float32, shape=(len(variants),1838))
    if start == len(valid):
        return features
    model = (load_pretrained_model if name == 'pretrained' else load_our_model)(device)
    if name == 'pretrained':

        for container in [model.features, model.classifier]:
            for i, module in enumerate(container):
                if isinstance(module, torch.nn.ReLU):
                    container[i] = torch.nn.Threshold(0, 1e-6)
    alt = variants.alt.map(dict(zip('AGCT', range(4)))).to_numpy()
    begin = time.time()
    with torch.inference_mode():
        for offset in range(start, len(valid), batch_size):
            indices = valid[offset:offset+batch_size]
            features[indices] = effects(model, bases[indices], alt[indices].astype(np.int64), device)
            end = min(offset+batch_size, len(valid))
            if offset == start or end == len(valid) or (offset//batch_size)%100 == 0:
                features.flush()
                progress.write_text(json.dumps({'completed':end,'total':len(valid)}))
                speed = (end-start)/max(time.time()-begin, 1e-9)
                print(f'{name}: {end:,}/{len(valid):,}; {speed:.1f} variants/s; ETA {(len(valid)-end)/speed/60:.1f} min', flush=True)
    del model
    torch.cuda.empty_cache()
    return features


def evaluate(df, features, out, task, model_name, groups=GROUPS,
             train_group='Random negative SNP', reg_alpha=0, reg_lambda=10, iterations=100):
    labels = df.positive.astype(int).to_numpy()
    base = df.positive | df.label.eq(train_group)
    oof = np.full(len(df), np.nan)
    for fold in range(10):
        train = (df.fold != fold) & base
        test = df.fold == fold
        if not test.any():
            continue
        fold_path = out / f'{task}_{model_name}_fold{fold}.npz'
        indices_all = np.flatnonzero(test.to_numpy())
        if fold_path.exists():
            cached = np.load(fold_path)
            if not np.array_equal(cached['test_indices'],indices_all):
                raise RuntimeError('Fold cache does not match cohort')
            oof[indices_all] = cached['probabilities']
            print(f'{task}/{model_name}: restored fold {fold}',flush=True)
            continue
        scaler = StandardScaler()
        x = scaler.fit_transform(np.asarray(features[df.loc[train,'variant_id']], dtype=np.float32))
        y = labels[train]
        weights = np.where(y == 1, len(y)/(2*(y == 1).sum()), len(y)/(2*(y == 0).sum()))


        penalty_scale = float(weights.sum())
        class Progress(xgb.callback.TrainingCallback):
            def after_iteration(self, model, epoch, evals_log):
                if (epoch+1)%20 == 0 or epoch+1 == iterations:
                    print(f'{task}/{model_name}/fold{fold}: classifier iteration {epoch+1}/{iterations}',flush=True)
                return False
        clf = xgb.XGBClassifier(booster='gblinear', updater='coord_descent', feature_selector='cyclic',
                                device='cuda' if torch.cuda.is_available() else 'cpu',
                                reg_alpha=reg_alpha/penalty_scale, reg_lambda=reg_lambda/penalty_scale,
                                n_estimators=iterations, learning_rate=0.1, callbacks=[Progress()],
                                n_jobs=8, random_state=42, eval_metric='logloss')
        clf.fit(x, y, sample_weight=weights)
        del x
        for start in range(0, int(test.sum()), 10000):
            indices = np.flatnonzero(test.to_numpy())[start:start+10000]
            x = scaler.transform(np.asarray(features[df.iloc[indices].variant_id], dtype=np.float32))
            oof[indices] = clf.predict_proba(x)[:,1]
        np.savez(fold_path,test_indices=indices_all,probabilities=oof[indices_all],
                 scaler_mean=scaler.mean_,scaler_scale=scaler.scale_,
                 reg_alpha=reg_alpha/penalty_scale,reg_lambda=reg_lambda/penalty_scale)
        clf.save_model(out / f'{task}_{model_name}_fold{fold}_classifier.json')
        del clf, x
        print(f'{task}/{model_name}: evaluated fold {fold}', flush=True)
    saved = df[['chr','pos','ref','alt','label','fold','variant_id','positive','distance_bp']].copy()
    saved['oof_probability'] = oof
    saved.to_csv(out / f'{task}_{model_name}_oof.tsv.gz', sep='\t', index=False)
    rows = []
    for group in groups:
        selection = df.positive | df.label.eq(group)
        neg = df.label.eq(group)
        rows.append(dict(task=task, model=model_name, group=group,
                         mean_distance_bp=float(df.loc[neg,'distance_bp'].mean()),
                         n_positive=int(df.positive.sum()), n_negative=int(neg.sum()),
                         auc=float(roc_auc_score(labels[selection],oof[selection]))))
    return rows


def plot_pretrained_and_ours_separately(summary, out):
    """Create independent pretrained/ours figures against real method scores."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    task_info = {
        'eqtl': ('GRASP eQTL (noncoding)', 649),
        'gwas': ('GWAS Catalog (noncoding)', 650),
    }
    comparison_rows = []
    baseline_results = {}
    common_cohorts = {}
    positive_counts = {}
    gwava_scores = pd.read_csv(out/'gwava_GRCh37_scores.tsv.gz', sep='\t')
    funseq_path = out/'funseq2_v2.1.0_scores.tsv.gz'
    raw_funseq = pd.read_csv(funseq_path, sep='\t', dtype={'chr':str,'ref':str,'alt_set':str})
    if not raw_funseq.empty:
        raw_funseq['funseq2_score'] = pd.to_numeric(raw_funseq.funseq2_score, errors='coerce')
        raw_funseq = raw_funseq.dropna(subset=['funseq2_score']).drop_duplicates()


        exact = raw_funseq[raw_funseq.alt_set.ne('.')]
        exact_key = ['chr','pos','ref','alt_set']
        conflict = exact.groupby(exact_key).funseq2_score.nunique()
        bad = set(conflict[conflict.gt(1)].index.tolist())
        if bad:
            exact = exact[~exact.set_index(exact_key).index.isin(bad)]
        generic = raw_funseq[raw_funseq.alt_set.eq('.')]
        generic_conflict = generic.groupby(['chr','pos','ref']).funseq2_score.nunique()
        bad_generic = set(generic_conflict[generic_conflict.gt(1)].index.tolist())
        if bad_generic:
            generic = generic[~generic.set_index(['chr','pos','ref']).index.isin(bad_generic)]
        raw_funseq = pd.concat([exact, generic], ignore_index=True)
    funseq_available = not raw_funseq.empty
    keys = ['chr','pos','ref','alt','label']


    cohorts = {}
    for task in task_info:
        full = pd.read_csv(out/f'{task}_pretrained_oof.tsv.gz', sep='\t')
        random_neg = full[full.label.eq(GROUPS[0])].sample(frac=.05, random_state=42)
        retained = full[full.positive | ~full.label.eq(GROUPS[0])]
        cohort = pd.concat([retained, random_neg], ignore_index=True)
        cohorts[task] = cohort

        scores = pd.read_csv(out/f'{task}_CADD_scores.tsv.gz', sep='\t')
        scores = cohort[keys+['positive','distance_bp']].merge(
            scores[keys+['CADD_PHRED']], on=keys, how='left', validate='one_to_one')
        gw = cohort[keys+['positive','distance_bp']].merge(
            gwava_scores, on=['chr','pos'], how='left', validate='many_to_one')
        fun = None
        if funseq_available:
            exact_scores = raw_funseq[raw_funseq.alt_set.ne('.')][['chr','pos','ref','alt_set','funseq2_score']].rename(columns={'alt_set':'alt'})
            generic_scores = raw_funseq[raw_funseq.alt_set.eq('.')][['chr','pos','ref','funseq2_score']].drop_duplicates(['chr','pos','ref'])
            fun = cohort[keys+['positive','distance_bp']].merge(exact_scores, on=['chr','pos','ref','alt'], how='left', validate='many_to_one')
            fallback = cohort[['chr','pos','ref']].merge(generic_scores,on=['chr','pos','ref'],how='left',validate='many_to_one').funseq2_score
            fun['funseq2_score'] = fun.funseq2_score.where(fun.funseq2_score.notna(), fallback)

        common = cohort[keys+['positive','distance_bp']].merge(
            scores[keys+['CADD_PHRED']], on=keys, how='left', validate='one_to_one')
        gw_cols = ['GWAVA_Unmatched','GWAVA_TSS','GWAVA_Region']
        common = common.merge(gwava_scores[['chr','pos']+gw_cols], on=['chr','pos'],
                              how='left', validate='many_to_one')
        score_cols = ['CADD_PHRED'] + gw_cols
        if fun is not None:
            common = common.merge(fun[keys+['funseq2_score']], on=keys, how='left', validate='one_to_one')
            score_cols.append('funseq2_score')
        common = common.dropna(subset=score_cols).copy()
        common_cohorts[task] = common
        positive_counts[task] = int(common.positive.sum())
        score_frames = {'CADD (PHRED; version unverified)':(common,'CADD_PHRED')}
        score_frames.update({name:(common,name) for name in gw_cols})
        if fun is not None:
            score_frames['FunSeq2 v2.1.0'] = (common,'funseq2_score')

        baseline_results[task] = {}
        for method,(frame,column) in score_frames.items():
            baseline_results[task][method] = {}
            frame[column] = pd.to_numeric(frame[column], errors='coerce')
            for group in GROUPS:
                selected = frame[frame.positive | frame.label.eq(group)].dropna(subset=[column])
                neg = frame[frame.label.eq(group)].dropna(subset=[column])
                row = {'task':task,'model':method,'group':group,
                       'mean_distance_bp':float(neg.distance_bp.mean()),
                       'n_positive':int(selected.positive.sum()),'n_negative':int(len(neg)),
                       'auc':float(roc_auc_score(selected.positive.astype(int),selected[column])),
                       'scored_positive':int(selected.positive.sum()),'scored_negative':int(len(neg))}
                baseline_results[task][method][group] = row
                comparison_rows.append(row)

    deepsea_results = {}
    for deepsea_name in ['pretrained','ours']:
        deepsea_results[deepsea_name] = {}
        for task in task_info:
            pred = pd.read_csv(out/f'{task}_{deepsea_name}_oof.tsv.gz', sep='\t')
            frame = common_cohorts[task][keys+['positive','distance_bp']].merge(
                pred[keys+['oof_probability']], on=keys, how='left', validate='one_to_one')
            deepsea_results[deepsea_name][task] = {}
            for group in GROUPS:
                selected = frame[frame.positive | frame.label.eq(group)].dropna(subset=['oof_probability'])
                neg = frame[frame.label.eq(group)].dropna(subset=['oof_probability'])
                row = {'task':task,'model':f'DeepSEA {deepsea_name}','group':group,
                       'mean_distance_bp':float(neg.distance_bp.mean()),
                       'n_positive':int(selected.positive.sum()),'n_negative':int(len(neg)),
                       'auc':float(roc_auc_score(selected.positive.astype(int),selected.oof_probability)),
                       'scored_positive':int(selected.positive.sum()),'scored_negative':int(len(neg))}
                deepsea_results[deepsea_name][task][group] = row
                comparison_rows.append(row)

    metric_values = np.asarray([row['auc'] for row in comparison_rows], dtype=float)
    scale_pad = max((float(metric_values.max())-float(metric_values.min()))*.10, .01)
    y_low = max(.5, float(np.floor((float(metric_values.min())-scale_pad)/.05)*.05))
    y_high = min(1., float(np.ceil((float(metric_values.max())+scale_pad)/.05)*.05))

    for deepsea_name in ['pretrained','ours']:

        fig, axes = plt.subplots(1,2,figsize=(10.3,3.35),gridspec_kw={'wspace':.34})
        for ax, task in zip(axes,task_info):
            title, _ = task_info[task]
            model_rows = [deepsea_results[deepsea_name][task][group] for group in GROUPS]
            cadd_rows = [baseline_results[task]['CADD (PHRED; version unverified)'][group] for group in GROUPS]
            x = np.arange(len(GROUPS))
            model_auc = np.asarray([row['auc'] for row in model_rows],dtype=float)
            cadd_auc = np.asarray([row['auc'] for row in cadd_rows],dtype=float)
            ax.plot(x,model_auc,'o-',lw=1.45,ms=3.1,color='#13a9d2',label=f'DeepSEA {deepsea_name}')
            ax.plot(x,cadd_auc,'o-',lw=1.25,ms=2.8,color='#8a8a8a',label='CADD')
            for gwava_name, color, linestyle in [
                    ('GWAVA_Unmatched','#579d42','-'),('GWAVA_TSS','#579d42','--'),
                    ('GWAVA_Region','#579d42',':')]:
                scores = np.asarray([baseline_results[task][gwava_name][group]['auc'] for group in GROUPS],dtype=float)
                ax.plot(x,scores,'o'+linestyle,lw=1.15,ms=2.6,color=color,label=gwava_name.replace('_',' '))
            if funseq_available:
                scores = np.asarray([baseline_results[task]['FunSeq2 v2.1.0'][group]['auc'] for group in GROUPS],dtype=float)
                ax.plot(x,scores,'o-',lw=1.15,ms=2.6,color='#444444',label='FunSeq2 v2.1.0')
            labels = ['All']+[f"{row['mean_distance_bp']:,.0f}" for row in cadd_rows[1:]]
            n_positive=positive_counts[task]
            ax.set_xticks(x,labels)
            ax.set_title(f'{title}\n($n$ = {n_positive:,})',fontsize=9,pad=4)
            ax.set_xlabel('Negative SNP group (mean distance, bp)',fontsize=8,labelpad=3)
            ax.set_ylabel('AUC',fontsize=8)
            ax.set_ylim(y_low,y_high)
            ax.tick_params(axis='both',labelsize=7.5,length=2,pad=2)
            ax.spines[['top','right']].set_visible(False)
        handles, legend_labels = axes[0].get_legend_handles_labels()
        fig.legend(handles,legend_labels,frameon=False,fontsize=7.4,loc='center left',
                   bbox_to_anchor=(.805,.52),borderaxespad=0)
        fig.subplots_adjust(left=.075,right=.785,bottom=.22,top=.94,wspace=.34)
        for extension in ['png','pdf']:
            FIGURE_DIR.mkdir(parents=True, exist_ok=True)
            fig.savefig(FIGURE_DIR/f'figure3_{deepsea_name}_separate.{extension}',dpi=240,bbox_inches='tight')
        plt.close(fig)

    pd.DataFrame(comparison_rows).to_csv(out/'pretrained_ours_baselines_auc.csv',index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--random-cap',type=int, default=0)
    parser.add_argument('--batch-size',type=int, default=128)
    parser.add_argument('--output',type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    torch.set_num_threads(8)
    torch.backends.cudnn.benchmark = True
    tables = read_tables(args.random_cap)
    ours_checkpoint = resolve_ours_checkpoint(ROOT)
    pretrained_checkpoint = resolve_pretrained_predict_checkpoint(ROOT)
    input_paths = [pretrained_checkpoint, ours_checkpoint] + list((ROOT/'newdata').glob('*649*csv')) + list((ROOT/'newdata').glob('*650*csv'))
    input_hashes = {
        ("Hugging Face: aer0vane/reproduce_deepsea/training_checkpoints/best_model_FINAL_EPOCH53.pth" if p == ours_checkpoint
         else "DeepSEA pretrained sequence model" if p == pretrained_checkpoint
         else str(p.relative_to(ROOT))): sha256(p)
        for p in input_paths
    }
    manifest = {'random_negative_cap':args.random_cap, 'method':'1838 ref/alt chromatin effects; strand effects averaged; gblinear; 10 supplied spatial folds',
                'regularization_convention':'legacy summed-loss parameters divided by training weight sum for modern gblinear',
                'classifier_solver':'GPU coordinate descent when CUDA available; cyclic features',
                'missing':['HGMD 2014.4 original variants','four conservation features'],
                'baseline_methods':['CADD PHRED (version unverified)','GWAVA v1.0 Unmatched/TSS/Region','FunSeq2 v2.1.0'],
                'source':'https://www.nature.com/articles/nmeth.3547',
                'inputs': input_hashes,
                'reference_genome_sha256':{p.name:sha256(p) for p in sorted((ROOT/'reference_genome').glob('chr*.fa')) if p.name!='chrM.fa'}}
    manifest_path = args.output/'manifest.json'
    if manifest_path.exists() and json.loads(manifest_path.read_text()) != manifest:
        raise RuntimeError('Output cache belongs to different inputs/options; use a new --output directory')
    manifest_path.write_text(json.dumps(manifest,indent=2))
    variants, bases = prepare(tables,args.output)
    for task, df in tables.items():
        df['distance_bp'] = nearest_distances(df)
        df = df.merge(variants, on=['chr','pos','ref','alt'],how='left',validate='many_to_one')
        print(task, 'exclusions by group:', df.groupby('label').valid.agg(['size','sum']).to_dict(),flush=True)

        df = df[df.valid].sort_values('positive',ascending=False,kind='stable').drop_duplicates(['chr','pos','label'])
        tables[task] = df.reset_index(drop=True)
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    rows = []
    for model_name in ['pretrained','ours']:
        features = infer(variants,bases,args.output,model_name,args.batch_size,device)
        for task, df in tables.items():


            result = evaluate(df,features,args.output,task,model_name)
            (args.output/f'{task}_{model_name}_summary.json').write_text(json.dumps(result,indent=2))
            rows.extend(result)
    summary = pd.DataFrame(rows)
    summary.to_csv(args.output/'auc_summary.csv',index=False)
    plot_pretrained_and_ours_separately(summary,args.output)
    write_metrics_provenance(args.output)
    print(summary.to_string(index=False),flush=True)


if __name__ == '__main__':
    main()
