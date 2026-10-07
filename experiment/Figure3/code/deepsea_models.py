#!/usr/bin/env python3
"""DeepSEA architectures and checkpoint loaders used by the Figure 3 rerun."""

from pathlib import Path
import sys

import torch
import torch.nn as nn

HERE = Path(__file__).resolve().parent
WORK_DIR = HERE.parent.parent.parent
sys.path.insert(0, str(WORK_DIR))
from model_assets import resolve_ours_checkpoint, resolve_pretrained_predict_checkpoint

MODEL_OURS = WORK_DIR / 'training' / 'checkpoints' / 'best_model_FINAL_EPOCH53.pth'
MODEL_PRETRAINED = resolve_pretrained_predict_checkpoint(WORK_DIR)


class DeepSEA(nn.Module):
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
        return self.classifier(x.reshape(x.size(0), -1))


def load_our_model(device='cuda'):
    model = DeepSEA()
    checkpoint_path = resolve_ours_checkpoint(WORK_DIR)
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    state = checkpoint.get('model_state_dict', checkpoint) if isinstance(checkpoint, dict) else checkpoint
    model.load_state_dict(state)
    return model.to(device).eval()


def load_pretrained_model(device='cuda'):
    model = DeepSEA()
    checkpoint = torch.load(MODEL_PRETRAINED, map_location=device, weights_only=False)
    state = checkpoint.get('model_state_dict', checkpoint) if isinstance(checkpoint, dict) else checkpoint
    key_map = {
        '2.0.': 'features.0.', '2.4.': 'features.4.', '2.8.': 'features.8.',
        '2.12.1.': 'classifier.0.', '2.14.1.': 'classifier.2.',
    }
    remapped = {}
    for key, value in state.items():
        new_key = next((new + key[len(old):] for old, new in key_map.items()
                        if key.startswith(old)), key)
        remapped[new_key] = value
    model.load_state_dict(remapped)
    return model.to(device).eval()
