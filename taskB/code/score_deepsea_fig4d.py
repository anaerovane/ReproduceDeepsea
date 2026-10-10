import json
import os
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch import nn
ROOT = Path(os.environ.get('FIG4_DIR', Path(__file__).resolve().parent.parent))
PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))
from model_assets import resolve_ours_checkpoint, resolve_pretrained_predict_checkpoint
DATA = ROOT / 'data'
OUT = ROOT / 'outputs'
GENOME = Path(os.environ.get('GENOME_DIR', PROJECT_ROOT / 'reference_genome'))
PRETRAINED = os.environ.get('DEEPSEA_PRETRAINED')
OURS = os.environ.get('DEEPSEA_OURS')
CONTEXT = 1000
BASE_INDEX = {base: i for i, base in enumerate('AGCT')}

class DeepSEA(nn.Module):

    def __init__(self):
        super().__init__()
        self.features = nn.Sequential(nn.Conv2d(4, 320, (1, 8)), nn.ReLU(), nn.MaxPool2d((1, 4)), nn.Dropout(0.2), nn.Conv2d(320, 480, (1, 8)), nn.ReLU(), nn.MaxPool2d((1, 4)), nn.Dropout(0.2), nn.Conv2d(480, 960, (1, 8)), nn.ReLU(), nn.Dropout(0.5))
        self.classifier = nn.Sequential(nn.Linear(50880, 925), nn.ReLU(), nn.Linear(925, 919), nn.Sigmoid())

    def forward(self, x):
        return self.classifier(self.features(x).reshape(x.size(0), -1))

def load_model(name):
    model = DeepSEA()
    if name == 'pretrained':
        path = Path(PRETRAINED) if PRETRAINED else resolve_pretrained_predict_checkpoint(PROJECT_ROOT)
    else:
        path = Path(OURS) if OURS else resolve_ours_checkpoint(PROJECT_ROOT)
    state = torch.load(path, map_location='cpu', weights_only=False)
    state = state.get('model_state_dict', state) if isinstance(state, dict) else state
    if name == 'pretrained':
        key_map = {'2.0.': 'features.0.', '2.4.': 'features.4.', '2.8.': 'features.8.', '2.12.1.': 'classifier.0.', '2.14.1.': 'classifier.2.'}
        state = {next((new + key[len(old):] for old, new in key_map.items() if key.startswith(old)), key): value for key, value in state.items()}
        for container in (model.features, model.classifier):
            for i, layer in enumerate(container):
                if isinstance(layer, nn.ReLU):
                    container[i] = nn.Threshold(0, 1e-06)
    model.load_state_dict(state)
    model.cuda().eval()
    return model

def load_sequence(chrom):
    with (GENOME / f'chr{chrom}.fa').open('rt') as handle:
        return ''.join((line.strip() for line in handle if not line.startswith('>')))

def encode_variant(row, genome):
    pos = int(row.pos) - 1
    start = pos - CONTEXT // 2
    end = start + CONTEXT
    if start < 0 or end > len(genome):
        return None
    seq = genome[start:end].upper()
    center = pos - start
    ref, alt = (str(row.ref).upper(), str(row.alt).upper())
    if len(ref) != 1 or len(alt) != 1 or ref not in BASE_INDEX or (alt not in BASE_INDEX) or (seq[center] != ref):
        return None
    codes = np.fromiter((BASE_INDEX.get(base, 4) for base in seq), dtype=np.uint8)
    if np.any(codes == 4):
        return None
    codes[center] = BASE_INDEX[alt]
    return codes

def score_variants(model, frame, model_name, batch_size=128):
    alt_codes, valid_rows = ([], [])
    for chrom, indices in frame.groupby('chrom', sort=False).groups.items():
        genome = load_sequence(chrom)
        for row in frame.loc[indices].itertuples():
            codes = encode_variant(row, genome)
            if codes is not None:
                alt_codes.append(codes)
                valid_rows.append(row.Index)
        print(f'{model_name}: chr{chrom} valid={len(valid_rows):,}', flush=True)
        del genome
    if not valid_rows:
        raise RuntimeError('No valid reference-matched variants')
    alt_codes = np.stack(alt_codes)
    scores = np.full(len(frame), np.nan, dtype=np.float32)
    started = time.time()
    with torch.inference_mode():
        for offset in range(0, len(valid_rows), batch_size):
            idx = np.asarray(valid_rows[offset:offset + batch_size])
            alt = torch.as_tensor(alt_codes[offset:offset + batch_size].astype(np.int64), device='cuda')
            ref = alt.clone()
            ref[:, 500] = torch.as_tensor([BASE_INDEX[str(frame.iloc[i].ref).upper()] for i in idx], dtype=torch.long, device='cuda')
            ref_hot = torch.nn.functional.one_hot(ref, 4).permute(0, 2, 1).float().unsqueeze(2)
            alt_hot = torch.nn.functional.one_hot(alt, 4).permute(0, 2, 1).float().unsqueeze(2)
            all_inputs = torch.cat([ref_hot, alt_hot, ref_hot.flip([1, 3]), alt_hot.flip([1, 3])], dim=0)
            pred = model(all_inputs)
            n = len(idx)
            effect = ((pred[:n] - pred[n:2 * n]).abs() + (pred[2 * n:3 * n] - pred[3 * n:]).abs()) / 2
            scores[idx] = effect.sum(dim=1).float().cpu().numpy()
            if offset == 0 or (offset // batch_size + 1) % 100 == 0 or offset + batch_size >= len(valid_rows):
                torch.cuda.synchronize()
                speed = (offset + len(idx)) / max(time.time() - started, 1e-08)
                print(f'{model_name}: variants={offset + len(idx):,}/{len(valid_rows):,}; {speed:.1f}/s', flush=True)
    frame_out = frame.copy()
    frame_out[f'deepsea_{model_name}_sum_abs_delta'] = scores
    return frame_out

def main():
    if not torch.cuda.is_available():
        raise RuntimeError('CUDA is required for DeepSEA inference')
    torch.set_num_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = [DATA / 'figure4d_evaluation_positives.tsv.gz', DATA / 'figure4d_negative_variants.tsv.gz']
    present = [path for path in inputs if path.is_file()]
    frames = [pd.read_csv(path, sep='\t', dtype={'chrom': str, 'ref': str, 'alt': str}) for path in present]
    frame = pd.concat(frames, ignore_index=True).drop_duplicates(['cohort', 'chrom', 'pos', 'ref', 'alt'])
    for name in ('pretrained', 'ours'):
        target = OUT / f'deepsea_{name}_scores.tsv.gz'
        if target.is_file():
            previous = pd.read_csv(target, sep='\t', dtype={'chrom': str, 'ref': str, 'alt': str})
            have = set(zip(previous.chrom, previous.pos, previous.ref, previous.alt))
            missing = frame[[tuple(x) not in have for x in frame[['chrom', 'pos', 'ref', 'alt']].itertuples(index=False, name=None)]]
        else:
            previous = pd.DataFrame()
            missing = frame
        if len(missing):
            model = load_model(name)
            scored = score_variants(model, missing.reset_index(drop=True), name)
            del model
            torch.cuda.empty_cache()
            combined = pd.concat([previous, scored], ignore_index=True).drop_duplicates(['cohort', 'chrom', 'pos', 'ref', 'alt'], keep='last')
            combined.to_csv(target, sep='\t', index=False, compression='gzip')
        scored_count = len(pd.read_csv(target, sep='\t'))
        print(f'{name}: cached={scored_count:,}', flush=True)
if __name__ == '__main__':
    main()
