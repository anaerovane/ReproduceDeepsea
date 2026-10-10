import json
import os
import random
import time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from transformers import AutoModelForSequenceClassification, AutoTokenizer
ROOT = Path(os.environ.get('FIG4_DIR', Path(__file__).resolve().parent.parent))
DATA = ROOT / 'data'
PROJECT_ROOT = Path(__file__).resolve().parents[3]
GENOME = Path(os.environ.get('GENOME_DIR', PROJECT_ROOT / 'reference_genome'))
MODEL = Path(os.environ['NT_MODEL_PATH']) if os.environ.get('NT_MODEL_PATH') else None
POSITIVES = DATA / 'figure4d_evaluation_positives.tsv.gz'
NEGATIVES = DATA / 'figure4d_negative_variants.tsv.gz'
OUT = ROOT / 'outputs'
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261010
CONTEXT = 6000
MAX_TOKENS = 1100
BATCH = int(os.environ.get('NT_BATCH_SIZE', '8'))
EPOCHS = int(os.environ.get('NT_EPOCHS', '1'))
COHORTS = ('eqtl', 'meqtl', 'clinvar')

def genome_sequence(chrom):
    path = GENOME / f'chr{chrom}.fa'
    with path.open('rt') as handle:
        return ''.join((line.strip() for line in handle if not line.startswith('>')))

def make_alt_sequence(row, genome):
    zero = int(row.pos) - 1
    left = CONTEXT // 2
    start = zero - left
    end = start + CONTEXT
    if start < 0 or end > len(genome):
        return None
    window = genome[start:end].upper()
    center = zero - start
    ref, alt = (str(row.ref).upper(), str(row.alt).upper())
    if len(ref) != 1 or len(alt) != 1 or ref not in 'ACGT' or (alt not in 'ACGT'):
        return None
    if window[center] != ref or any((base not in 'ACGT' for base in window)):
        return None
    return window[:center] + alt + window[center + 1:]

def load_cohort(cohort):
    pos = pd.read_csv(POSITIVES, sep='\t', dtype={'chrom': str, 'ref': str, 'alt': str})
    neg = pd.read_csv(NEGATIVES, sep='\t', dtype={'chrom': str, 'ref': str, 'alt': str})
    frame = pd.concat([pos[pos.cohort.eq(cohort)], neg[neg.cohort.eq(cohort)]], ignore_index=True)
    frame['positive'] = frame.positive.astype(bool)
    frame['chrom'] = frame.chrom.astype(str)
    frame['sequence'] = None
    for chrom, indices in frame.groupby('chrom', sort=False).groups.items():
        genome = genome_sequence(chrom)
        frame.loc[indices, 'sequence'] = [make_alt_sequence(row, genome) for row in frame.loc[indices].itertuples()]
        print(f'{cohort}: chr{chrom} sequences={len(indices):,}', flush=True)
    frame = frame.dropna(subset=['sequence']).reset_index(drop=True)
    counts = frame.positive.value_counts().to_dict()
    if counts.get(True, 0) < 100 or counts.get(False, 0) < 100:
        raise RuntimeError(f'Insufficient balanced, ref-matched data for {cohort}: {counts}')
    return frame

def encode(tokenizer, seqs):
    batch = tokenizer(seqs, return_tensors='pt', padding=True, truncation=True, max_length=MAX_TOKENS)
    return {k: v.cuda() for k, v in batch.items()}

def run():
    if not torch.cuda.is_available():
        raise RuntimeError('CUDA is required for this run')
    if MODEL is None or not MODEL.is_dir():
        raise RuntimeError('Set NT_MODEL_PATH to the local pinned NT-v2 50M checkpoint directory')
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    torch.cuda.manual_seed_all(SEED)
    torch.set_num_threads(4)
    tokenizer = AutoTokenizer.from_pretrained(str(MODEL), trust_remote_code=True)
    all_scores = []
    manifest = {'model': 'InstaDeepAI/nucleotide-transformer-v2-50m-multi-species', 'model_path': str(MODEL), 'input_context_bp': CONTEXT, 'max_tokens': MAX_TOKENS, 'cv': '5-fold chromosome GroupKFold; chromosome groups never overlap between train and test', 'seed': SEED, 'batch_size': BATCH, 'epochs_per_fold': EPOCHS, 'cohorts': {}, 'device': torch.cuda.get_device_name(0), 'status': 'running'}
    (OUT / 'nt50m_fig4d_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    for cohort in COHORTS:
        frame = load_cohort(cohort)
        group_splitter = GroupKFold(n_splits=5)
        labels = frame.positive.astype(int).to_numpy()
        groups = frame.chrom.to_numpy()
        cohort_scores = np.full(len(frame), np.nan, dtype=np.float32)
        fold_records = []
        for fold, (train_idx, test_idx) in enumerate(group_splitter.split(frame, labels, groups), 1):
            started = time.time()
            print(f'{cohort} fold={fold}/5 train={len(train_idx):,} test={len(test_idx):,}', flush=True)
            model = AutoModelForSequenceClassification.from_pretrained(str(MODEL), trust_remote_code=True, num_labels=2, ignore_mismatched_sizes=True).cuda()
            optimizer = torch.optim.AdamW(model.parameters(), lr=1e-05, weight_decay=0.01)
            model.train()
            rng = np.random.default_rng(SEED + fold)
            order = rng.permutation(train_idx)
            batch_losses = []
            for step, offset in enumerate(range(0, len(order), BATCH), 1):
                idx = order[offset:offset + BATCH]
                tokens = encode(tokenizer, frame.iloc[idx].sequence.tolist())
                y = torch.as_tensor(labels[idx], dtype=torch.long, device='cuda')
                optimizer.zero_grad(set_to_none=True)
                result = model(**tokens, labels=y)
                result.loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()
                batch_losses.append(float(result.loss.detach().cpu()))
                if step % 100 == 0 or step * BATCH >= len(order):
                    torch.cuda.synchronize()
                    elapsed = time.time() - started
                    print(f'{cohort} fold={fold}/5 step={step}/{(len(order) + BATCH - 1) // BATCH} loss={np.mean(batch_losses[-100:]):.5f} elapsed={elapsed / 60:.1f}min', flush=True)
            model.eval()
            predictions = []
            with torch.inference_mode():
                for offset in range(0, len(test_idx), BATCH):
                    idx = test_idx[offset:offset + BATCH]
                    tokens = encode(tokenizer, frame.iloc[idx].sequence.tolist())
                    predictions.extend(model(**tokens).logits.softmax(-1)[:, 1].cpu().numpy().tolist())
            cohort_scores[test_idx] = predictions
            fold_auc = roc_auc_score(labels[test_idx], cohort_scores[test_idx])
            fold_record = {'fold': fold, 'train_count': len(train_idx), 'test_count': len(test_idx), 'train_chromosomes': sorted(set(groups[train_idx]), key=lambda x: (not str(x).isdigit(), int(x) if str(x).isdigit() else str(x))), 'test_chromosomes': sorted(set(groups[test_idx]), key=lambda x: (not str(x).isdigit(), int(x) if str(x).isdigit() else str(x))), 'mean_train_loss': float(np.mean(batch_losses)), 'test_auc': float(fold_auc), 'elapsed_seconds': time.time() - started}
            fold_records.append(fold_record)
            print(f'{cohort} fold={fold}/5 auc={fold_auc:.6f} elapsed={time.time() - started:.1f}s', flush=True)
            del model, optimizer
            torch.cuda.empty_cache()
            partial = frame[['cohort', 'label', 'chrom', 'pos', 'ref', 'alt', 'positive']].copy()
            partial['fold'] = 0
            partial.loc[test_idx, 'fold'] = fold
            partial['nt_finetuned_score'] = cohort_scores
            partial = partial.loc[np.isfinite(cohort_scores)]
            partial.to_csv(OUT / f'{cohort}_nt50m_oof.partial.tsv.gz', sep='\t', index=False, compression='gzip')
        scored = frame[['cohort', 'label', 'chrom', 'pos', 'ref', 'alt', 'positive']].copy()
        scored['nt_finetuned_score'] = cohort_scores
        scored.to_csv(OUT / f'{cohort}_nt50m_oof.tsv.gz', sep='\t', index=False, compression='gzip')
        auc = roc_auc_score(labels, cohort_scores)
        manifest['cohorts'][cohort] = {'positive_count': int(labels.sum()), 'negative_count': int(len(labels) - labels.sum()), 'ref_matched_count': int(len(frame)), 'oof_auc': float(auc), 'folds': fold_records}
        (OUT / 'nt50m_fig4d_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
        all_scores.append(scored)
        print(f'{cohort} OOF_AUC={auc:.6f}', flush=True)
    pd.concat(all_scores, ignore_index=True).to_csv(OUT / 'nt50m_fig4d_oof.tsv.gz', sep='\t', index=False, compression='gzip')
    manifest['status'] = 'complete'
    (OUT / 'nt50m_fig4d_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
if __name__ == '__main__':
    run()
