"""
Supplementary Figure 1: overall AUC box plots for pretrained and ours
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pathlib import Path

WORK = Path(__file__).resolve().parent


data = np.load(WORK / 'supplementary_figure1_auc_data.npz')

contexts = [200, 500, 1000]

for model_name, key_prefix, color, title in [
    ('pretrained', 'pretrained', '#2166AC', 'Pretrained'),
    ('ours', 'ours', '#D73027', 'Ours (epoch 53)')
]:
    fig, ax = plt.subplots(figsize=(4, 4))

    box_data = []
    for cl in contexts:
        aucs = data[f'{key_prefix}_{cl}']
        box_data.append(aucs[~np.isnan(aucs)])

    bp = ax.boxplot(box_data, positions=[0, 1, 2], widths=0.5,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color='black', linewidth=1.2))

    for patch in bp['boxes']:
        patch.set_facecolor(color)
        patch.set_alpha(0.7)


    for i, bd in enumerate(box_data):
        med = np.median(bd)
        ax.text(i, med + 0.008, f'{med:.3f}', ha='center', va='bottom', fontsize=8)

    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['200bp', '500bp', '1000bp'], fontsize=10)
    ax.set_ylabel('Per-feature AUC', fontsize=10)
    ax.set_xlabel('Context length', fontsize=10)
    ax.set_ylim(0.5, 1.0)
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(length=3)
    ax.set_title(title, fontsize=12, fontweight='bold', color=color)

    plt.tight_layout()
    out = WORK / f'supplementary_figure1_{model_name}.png'
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}")
