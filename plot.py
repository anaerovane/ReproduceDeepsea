"""
用真实复现数据画图：
  Figure 2a: ROC (test_labels + test_predictions)
  Figure 2b: 散点图 (dgf_ref_preds + dgf_alt_preds)
  Figure 2c: Accuracy vs Margin (dgf_ref_preds + dgf_alt_preds)

Usage:
  python3 plot.py pretrained
  python3 plot.py ours
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import matplotlib.gridspec as gridspec
from sklearn.metrics import roc_curve, auc, accuracy_score
from huggingface_hub import hf_hub_download
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.size': 8,
    'axes.labelsize': 9,
    'axes.titlesize': 9,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 6,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'font.family': 'sans-serif',
    'axes.edgecolor': 'black',
    'axes.linewidth': 0.8,
})

WORK = os.path.dirname(os.path.abspath(__file__))
OUT  = os.path.join(WORK, 'experiment', 'Figure2')
HF_REPO_ID = 'aer0vane/reproduce_deepsea'


def input_path(filename):
    """Resolve shared metadata, source inputs, or generated Figure 2 caches."""
    if filename == 'predictor_names.txt':
        local_path = os.path.join(WORK, filename)
        hub_path = filename
    elif filename == '41592_2015_BFnmeth3547_MOESM647_ESM.csv':
        local_path = os.path.join(WORK, 'experiment', 'Figure2', 'inputs', filename)
        hub_path = f'experiment/Figure2/inputs/{filename}'
    else:
        local_path = os.path.join(WORK, 'experiment', 'Figure2', 'outputs', filename)
        hub_path = f'experiment/Figure2/outputs/{filename}'
    if os.path.isfile(local_path):
        return local_path
    print(f"  Fetching {filename} from Hugging Face cache...")
    return hf_hub_download(repo_id=HF_REPO_ID, filename=hub_path, repo_type='model')


model_name = sys.argv[1] if len(sys.argv) > 1 else 'pretrained'
if model_name not in ('pretrained', 'ours'):
    print(f"Usage: python3 {sys.argv[0]} [pretrained|ours]")
    sys.exit(1)

os.makedirs(OUT, exist_ok=True)
print(f"=== Generating Figure 2 for: {model_name} ===")


print("Loading data...")


labels = np.load(input_path('test_labels.npy'))
pred_file = input_path(f'test_predictions_{model_name}.npy')
preds = np.load(pred_file)
with open(input_path('predictor_names.txt')) as f:
    predictor_names = [l.strip() for l in f if l.strip()]
print(f"  model={model_name}, test_labels={labels.shape}, test_preds={preds.shape}")


dgf_ref_path = input_path(f'dgf_ref_preds_{model_name}.npy')
dgf_alt_path = input_path(f'dgf_alt_preds_{model_name}.npy')
dgf_ref = np.load(dgf_ref_path)
dgf_alt = np.load(dgf_alt_path)
print(f"  dgf_ref={dgf_ref.shape}, dgf_alt={dgf_alt.shape}")


dgf_csv = pd.read_csv(input_path('41592_2015_BFnmeth3547_MOESM647_ESM.csv'), skiprows=1)

variants = dgf_csv[['CHR', 'POS', 'REF', 'ALT']].drop_duplicates().reset_index(drop=True)
variant_to_idx = {}
for i, row in variants.iterrows():
    key = (row['CHR'], row['POS'], row['REF'], row['ALT'])
    variant_to_idx[key] = i

dgf_csv['true_alt_biased'] = dgf_csv['Alt reads'] > dgf_csv['Ref reads']


grouped_mean = dgf_csv.groupby(['CHR', 'POS', 'REF', 'ALT'])['true_alt_biased'].mean()
true_dir_experimental = np.array([grouped_mean[(r['CHR'], r['POS'], r['REF'], r['ALT'])] > 0.5 for _, r in variants.iterrows()])
print(f"  Experimental labels: alt-biased={true_dir_experimental.sum()}, ref-biased={(~true_dir_experimental).sum()}")


def classify(name):
    if 'DNase' in name:
        return 'DNase-seq'
    if any(h in name for h in ['H3K','H4K','H2A','H2B','Histone']):
        return 'Histone marks'
    return 'TF binding'


print("\n=== Figure 2 ===")

fig = plt.figure(figsize=(13, 6.5))
fig.subplots_adjust(hspace=0.4, wspace=0.35, top=0.92, bottom=0.08, left=0.08, right=0.98)

gs = gridspec.GridSpec(2, 3, figure=fig, width_ratios=[1,1,1], height_ratios=[1,1])


axes_a = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])]

cats_a = [
    ('TF binding',    'Transcription factors'),
    ('DNase-seq',     'DNase I-hypersensitive sites'),
    ('Histone marks', 'Histone marks'),
]


model_color = '#2166AC' if model_name == 'pretrained' else '#D73027'
model_label = 'Pretrained' if model_name == 'pretrained' else 'Ours (epoch 53)'

for ax, (cat, title) in zip(axes_a, cats_a):
    idxs = [i for i, n in enumerate(predictor_names) if classify(n) == cat]
    aucs = []
    for i in idxs:
        yt = labels[:, i]
        ys = preds[:, i]
        if yt.sum() > 0 and yt.sum() < len(yt):
            fpr, tpr, _ = roc_curve(yt, ys)
            aucs.append(auc(fpr, tpr))
            ax.plot(fpr, tpr, color=model_color, alpha=0.15, lw=0.4)


    if aucs:
        mean_auc = np.mean(aucs)
        ax.text(0.6, 0.15, f'{model_label}\nmean AUC={mean_auc:.3f}',
                fontsize=7, color=model_color, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    ax.plot([0,1],[0,1],'k--',lw=0.6,alpha=0.5)
    ax.set_xlim(0,1); ax.set_ylim(0,1)
    ax.set_title(title, fontsize=9, fontweight='normal')
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    ax.set_aspect('equal')
    ax.tick_params(direction='out', length=2)

axes_a[0].text(-0.22, 1.18, 'a', fontsize=13, fontweight='bold', transform=axes_a[0].transAxes)


ax_b = fig.add_subplot(gs[1, 0:2])

dnase_idx = [i for i, n in enumerate(predictor_names) if 'DNase' in n]
if not dnase_idx:
    dnase_idx = [0]

p_ref = dgf_ref[:, dnase_idx[0]]
p_alt = dgf_alt[:, dnase_idx[0]]

effect = p_alt - p_ref
alt_biased = effect > 0.05
ref_biased = effect < -0.05
neutral = ~(alt_biased | ref_biased)

print(f"  1b: total={len(p_ref)}, alt_biased={alt_biased.sum()}, ref_biased={ref_biased.sum()}, neutral={neutral.sum()}")

np.random.seed(42)
if neutral.sum() > 15000:
    ni = np.random.choice(np.where(neutral)[0], 15000, replace=False)
else:
    ni = np.where(neutral)[0]
ax_b.scatter(p_ref[ni], p_alt[ni], s=1.0, c='#BDBDBD', alpha=0.35, label='Neutral')

if alt_biased.sum() > 8000:
    ai = np.random.choice(np.where(alt_biased)[0], 8000, replace=False)
else:
    ai = np.where(alt_biased)[0]
ax_b.scatter(p_ref[ai], p_alt[ai], s=2.0, c='#D73027', alpha=0.7, label=f'Alt-biased (n={alt_biased.sum()})')

if ref_biased.sum() > 8000:
    ri = np.random.choice(np.where(ref_biased)[0], 8000, replace=False)
else:
    ri = np.where(ref_biased)[0]
ax_b.scatter(p_ref[ri], p_alt[ri], s=2.0, c='#4575B4', alpha=0.7, label=f'Ref-biased (n={ref_biased.sum()})')

margin = 0.07
ax_b.plot([0, 1-margin], [margin, 1], 'k-', lw=0.8, alpha=0.7)
ax_b.plot([margin, 1], [0, 1-margin], 'k-', lw=0.8, alpha=0.7)

ax_b.plot([0,1],[0,1],'k--',lw=0.6,alpha=0.5)
ax_b.set_xlim(0,1); ax_b.set_ylim(0,1)
ax_b.set_xlabel('P$_{DHS}$ (reference)')
ax_b.set_ylabel('P$_{DHS}$ (alternative)')
ax_b.set_aspect('equal')
ax_b.legend(loc='upper left', fontsize=6, markerscale=2, framealpha=0.9)
ax_b.tick_params(direction='out', length=2)
ax_b.text(-0.10, 1.12, 'b', fontsize=13, fontweight='bold', transform=ax_b.transAxes)


ax_c = fig.add_subplot(gs[1, 2])


c_color = '#2166AC' if model_name == 'pretrained' else '#D73027'
c_label = 'Pretrained' if model_name == 'pretrained' else 'Ours (epoch 53)'

margins = np.arange(0, 0.41, 0.02)

for idx, d_idx in enumerate(dnase_idx[:min(20, len(dnase_idx))]):
    p_ref_d = dgf_ref[:, d_idx]
    p_alt_d = dgf_alt[:, d_idx]
    margin_vals = np.abs(p_ref_d - p_alt_d)

    accs = []
    for m in margins:
        mask = margin_vals >= m
        if mask.sum() < 10:
            accs.append(np.nan)
            continue
        pred_dir = (p_alt_d[mask] > p_ref_d[mask]).astype(int)
        true_dir = true_dir_experimental[mask]
        acc = (pred_dir == true_dir).mean()
        accs.append(acc)

    accs = np.array(accs)
    valid = ~np.isnan(accs)
    if valid.sum() > 3:
        ax_c.plot(margins[valid], accs[valid], color=c_color, alpha=0.2, lw=0.5)

mean_accs = []
for m in margins:
    accs = []
    for d_idx in dnase_idx:
        p_ref_d = dgf_ref[:, d_idx]
        p_alt_d = dgf_alt[:, d_idx]
        margin_vals = np.abs(p_ref_d - p_alt_d)
        mask = margin_vals >= m
        if mask.sum() < 10:
            continue
        pred_dir = (p_alt_d[mask] > p_ref_d[mask]).astype(int)
        true_dir = true_dir_experimental[mask]
        acc = (pred_dir == true_dir).mean()
        accs.append(acc)
    mean_accs.append(np.mean(accs) if accs else np.nan)

mean_accs = np.array(mean_accs)
valid = ~np.isnan(mean_accs)
ax_c.plot(margins[valid], mean_accs[valid], '-', color=c_color, lw=2.0,
          label=f'{c_label} mean acc')


if valid.sum() > 0:
    final_acc = mean_accs[valid][-1]
    ax_c.text(0.35, 0.55 if model_name == 'pretrained' else 0.45,
              f'{c_label}: {final_acc:.1%}', fontsize=7, color=c_color, fontweight='bold')

ax_c.set_xlabel('Margin')
ax_c.set_ylabel('Accuracy')
ax_c.tick_params(direction='out', length=2)
ax_c.set_ylim(0.48, 1.02)
ax_c.set_xlim(-0.01, 0.42)
ax_c.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax_c.text(-0.18, 1.12, 'c', fontsize=13, fontweight='bold', transform=ax_c.transAxes)

suffix = f'_{model_name}'
plt.savefig(f'{OUT}/figure2{suffix}.png', dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"  Saved {OUT}/figure2{suffix}.png")

print("\n完成！")
