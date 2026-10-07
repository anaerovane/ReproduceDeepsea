"""
================================================================================
实验 #2: 上下文长度消融 (Context Length Ablation)
================================================================================
测试输入窗口大小 (1000bp / 500bp / 200bp) 对AUC的影响。
两个模型都跑：pretrained 和 ours。

做法：取中心子序列，两端填0补到1000bp，模型不变。
================================================================================
"""
import os
import sys
import numpy as np
import scipy.io
import torch
import torch.nn as nn
from sklearn.metrics import roc_curve, auc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, WORK)
from model_assets import resolve_ours_checkpoint, resolve_pretrained_predict_checkpoint


TEST_DATA_PATH = f'{WORK}/deepsea_train/test.mat'
PREDICT_MODEL_PATH = str(resolve_pretrained_predict_checkpoint(WORK))
OUR_MODEL_PATH = str(resolve_ours_checkpoint(WORK))
CONTEXT_LENGTHS = [200, 500, 1000]
FULL_LENGTH = 1000


class DeepSEA_Ours(nn.Module):
    def __init__(self):
        super().__init__()
        self.features = nn.Sequential(
            nn.Conv2d(4, 320, (1, 8)), nn.ReLU(), nn.MaxPool2d((1, 4)), nn.Dropout(0.2),
            nn.Conv2d(320, 480, (1, 8)), nn.ReLU(), nn.MaxPool2d((1, 4)), nn.Dropout(0.2),
            nn.Conv2d(480, 960, (1, 8)), nn.ReLU(), nn.Dropout(0.5),
        )
        self.classifier = nn.Sequential(
            nn.Linear(50880, 925), nn.ReLU(), nn.Linear(925, 919), nn.Sigmoid()
        )
    def forward(self, x):
        x = self.features(x)
        x = x.view(x.size(0), -1)
        return self.classifier(x)


class LambdaBase(nn.Sequential):
    def __init__(self, fn, *args):
        super().__init__(*args)
        self.lambda_func = fn
    def forward_prepare(self, input):
        output = []
        for module in self._modules.values():
            output.append(module(input))
        return output if output else input

class Lambda(LambdaBase):
    def forward(self, input):
        return self.lambda_func(self.forward_prepare(input))

class ReCodeAlphabet(nn.Module):
    def forward(self, input):
        input_reordered = [input[:, i, ...] for i in [0, 2, 1, 3]]
        return torch.stack(input_reordered, dim=1)

class ConcatenateRC(nn.Module):
    def forward(self, input):
        input_rc = input.clone()
        for idim in [1, 3]:
            idxs = list(range(input_rc.size(idim) - 1, -1, -1))
            input_rc = input_rc.flip(dims=[idim])
        return torch.cat([input, input_rc], dim=0)

class AverageRC(nn.Module):
    def forward(self, input):
        n = input.shape[0] // 2
        return input[:n] / 2 + input[n:] / 2

def build_backbone():
    return nn.Sequential(
        nn.Conv2d(4, 320, (1, 8), (1, 1)),
        nn.Threshold(0, 1e-06),
        nn.MaxPool2d((1, 4), (1, 4)),
        nn.Dropout(0.2),
        nn.Conv2d(320, 480, (1, 8), (1, 1)),
        nn.Threshold(0, 1e-06),
        nn.MaxPool2d((1, 4), (1, 4)),
        nn.Dropout(0.2),
        nn.Conv2d(480, 960, (1, 8), (1, 1)),
        nn.Threshold(0, 1e-06),
        nn.Dropout(0.5),
        Lambda(lambda x: x.view(x.size(0), -1)),
        nn.Sequential(Lambda(lambda x: x.view(1, -1) if len(x.size()) == 1 else x), nn.Linear(50880, 925)),
        nn.Threshold(0, 1e-06),
        nn.Sequential(Lambda(lambda x: x.view(1, -1) if len(x.size()) == 1 else x), nn.Linear(925, 919)),
        nn.Sigmoid(),
    )

def build_predict_model():
    backbone = build_backbone()
    return nn.Sequential(ReCodeAlphabet(), ConcatenateRC(), backbone, AverageRC())


def extract_center_pad(seq_data, target_len):
    """从中心提取target_len长度，两端填0补到FULL_LENGTH"""
    n, c, l = seq_data.shape
    if l == FULL_LENGTH:
        start = (FULL_LENGTH - target_len) // 2
        out = np.zeros_like(seq_data)
        out[:, :, start:start+target_len] = seq_data[:, :, start:start+target_len]
        return out
    else:
        raise ValueError(f"Expected length {FULL_LENGTH}, got {l}")

def compute_per_feature_auc(labels, preds):
    """计算每个特征的AUC，返回数组"""
    n_features = labels.shape[1]
    aucs = []
    for i in range(n_features):
        yt = labels[:, i]
        ys = preds[:, i]
        if yt.sum() > 0 and yt.sum() < len(yt):
            fpr, tpr, _ = roc_curve(yt, ys)
            aucs.append(auc(fpr, tpr))
        else:
            aucs.append(np.nan)
    return np.array(aucs)


def main():
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Device: {device}")


    print("Loading test data...")
    mat = scipy.io.loadmat(TEST_DATA_PATH)
    testxdata = mat['testxdata'].astype(np.float32)
    testlabels = mat['testdata']


    half = testxdata.shape[0] // 2
    testxdata_fwd = testxdata[:half]
    testlabels_fwd = testlabels[:half]
    print(f"  Sequences: {testxdata_fwd.shape}, Labels: {testlabels_fwd.shape}")

    results = {}


    print("\n" + "="*60)
    print("OUR MODEL")
    print("="*60)
    model_ours = DeepSEA_Ours()
    ckpt = torch.load(OUR_MODEL_PATH, map_location=device, weights_only=False)
    sd = ckpt['model_state_dict'] if isinstance(ckpt, dict) and 'model_state_dict' in ckpt else ckpt
    model_ours.load_state_dict(sd)
    model_ours.to(device)
    model_ours.eval()

    for cl in CONTEXT_LENGTHS:
        print(f"\n  Context length: {cl}bp")

        if cl == 1000:
            data = testxdata_fwd.copy()
        else:
            data = extract_center_pad(testxdata_fwd, cl)


        preds = []
        bs = 256
        with torch.no_grad():
            for i in range(0, half, bs):
                batch = data[i:i+bs]
                batch_tensor = torch.FloatTensor(batch).unsqueeze(2).to(device)
                pred_fwd = model_ours(batch_tensor)
                batch_rc = batch_tensor.flip(dims=[1, 3])
                pred_rc = model_ours(batch_rc)
                pred = ((pred_fwd + pred_rc) / 2).cpu().numpy()
                preds.append(pred)
        preds = np.vstack(preds)

        aucs = compute_per_feature_auc(testlabels_fwd, preds)
        valid_aucs = aucs[~np.isnan(aucs)]
        print(f"    AUC: median={np.median(valid_aucs):.4f}, mean={np.mean(valid_aucs):.4f}, n={len(valid_aucs)}")
        results[('ours', cl)] = aucs


    print("\n" + "="*60)
    print("PRETRAINED MODEL")
    print("="*60)
    model_pre = build_predict_model()
    sd = torch.load(PREDICT_MODEL_PATH, map_location=device)
    model_pre.load_state_dict(sd)
    model_pre.to(device)
    model_pre.eval()


    testxdata_acgt = testxdata_fwd[:, [0, 2, 1, 3], :].copy()

    for cl in CONTEXT_LENGTHS:
        print(f"\n  Context length: {cl}bp")
        if cl == 1000:
            data = testxdata_acgt.copy()
        else:
            data = extract_center_pad(testxdata_acgt, cl)


        preds = []
        bs = 256
        with torch.no_grad():
            for i in range(0, half, bs):
                batch = data[i:i+bs]
                batch_tensor = torch.FloatTensor(batch).unsqueeze(2).to(device)
                pred = model_pre(batch_tensor).cpu().numpy()
                preds.append(pred)
        preds = np.vstack(preds)

        aucs = compute_per_feature_auc(testlabels_fwd, preds)
        valid_aucs = aucs[~np.isnan(aucs)]
        print(f"    AUC: median={np.median(valid_aucs):.4f}, mean={np.mean(valid_aucs):.4f}, n={len(valid_aucs)}")
        results[('pretrained', cl)] = aucs


    print("\n" + "="*60)
    print("STATISTICAL TESTS")
    print("="*60)
    for model_name in ['ours', 'pretrained']:
        aucs_1000 = results[(model_name, 1000)]
        aucs_200 = results[(model_name, 200)]
        aucs_500 = results[(model_name, 500)]


        mask = ~(np.isnan(aucs_1000) | np.isnan(aucs_200))
        t_stat, p_val = stats.ttest_rel(aucs_1000[mask], aucs_200[mask])
        print(f"  {model_name}: 1000bp vs 200bp → t={t_stat:.2f}, p={p_val:.2e}")

        mask = ~(np.isnan(aucs_1000) | np.isnan(aucs_500))
        t_stat, p_val = stats.ttest_rel(aucs_1000[mask], aucs_500[mask])
        print(f"  {model_name}: 1000bp vs 500bp → t={t_stat:.2f}, p={p_val:.2e}")

        mask = ~(np.isnan(aucs_500) | np.isnan(aucs_200))
        t_stat, p_val = stats.ttest_rel(aucs_500[mask], aucs_200[mask])
        print(f"  {model_name}: 500bp vs 200bp → t={t_stat:.2f}, p={p_val:.2e}")


    print("\nPlotting...")

    with open(f'{WORK}/predictor_names.txt') as f:
        predictor_names = [l.strip() for l in f if l.strip()]

    def classify(name):
        if 'DNase' in name: return 'DNase-seq'
        if any(h in name for h in ['H3K','H4K','H2A','H2B','Histone']): return 'Histone'
        return 'TF'

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    colors = {'ours': '#D73027', 'pretrained': '#2166AC'}
    labels_plot = {'ours': 'Ours (epoch 53)', 'pretrained': 'Pretrained'}

    for idx, (model_name, ax) in enumerate(zip(['ours', 'pretrained'], axes)):
        data_for_box = []
        positions = []
        tick_labels = []

        for ci, cl in enumerate(CONTEXT_LENGTHS):
            aucs = results[(model_name, cl)]

            cats = {'TF': [], 'DNase-seq': [], 'Histone': []}
            for i, name in enumerate(predictor_names):
                cat = classify(name)
                if not np.isnan(aucs[i]):
                    cats[cat].append(aucs[i])

            offset = ci * 4
            for j, (cat, color) in enumerate(zip(['TF', 'DNase-seq', 'Histone'],
                                                   ['#4575B4', '#F4622E', '#74ADD1'])):
                positions.append(offset + j)
                data_for_box.append(cats[cat])

        bp = ax.boxplot(data_for_box, positions=positions, widths=0.6,
                        patch_artist=True, showfliers=False)


        color_cycle = ['#4575B4', '#F4622E', '#74ADD1'] * 3
        for patch, c in zip(bp['boxes'], color_cycle):
            patch.set_facecolor(c)
            patch.set_alpha(0.7)


        ax.set_xticks([1, 5, 9])
        ax.set_xticklabels(['200bp', '500bp', '1000bp'])
        ax.set_title(labels_plot[model_name], fontsize=11, fontweight='bold',
                     color=colors[model_name])
        ax.set_ylabel('Per-feature AUC')
        ax.set_ylim(0.4, 1.0)
        ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.3)


        aucs_1000 = results[(model_name, 1000)]
        aucs_200 = results[(model_name, 200)]
        mask = ~(np.isnan(aucs_1000) | np.isnan(aucs_200))
        _, p_val = stats.ttest_rel(aucs_1000[mask], aucs_200[mask])
        ax.text(0.5, 0.95, f'1000bp vs 200bp: p={p_val:.1e}',
                transform=ax.transAxes, ha='center', va='top', fontsize=8,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))


    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#4575B4', alpha=0.7, label='TF binding'),
                       Patch(facecolor='#F4622E', alpha=0.7, label='DNase-seq'),
                       Patch(facecolor='#74ADD1', alpha=0.7, label='Histone marks')]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3, fontsize=9,
               bbox_to_anchor=(0.5, -0.02))

    plt.suptitle('Experiment #2: Context Length Ablation', fontsize=12, fontweight='bold')
    plt.tight_layout(rect=[0, 0.04, 1, 0.96])

    out_path = f'{HERE}/experiment2_context_ablation.png'
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


    np.savez(f'{HERE}/experiment2_auc_data.npz',
             ours_200=results[('ours', 200)],
             ours_500=results[('ours', 500)],
             ours_1000=results[('ours', 1000)],
             pretrained_200=results[('pretrained', 200)],
             pretrained_500=results[('pretrained', 500)],
             pretrained_1000=results[('pretrained', 1000)])

    print("\nDone!")

if __name__ == '__main__':
    main()
