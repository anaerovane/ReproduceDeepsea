# Figure 3 outputs on Hugging Face

大型 Figure 3 推理结果、缓存和评估表托管在 Hugging Face：

- [Browse files](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Figure3/outputs)
- [GitHub source and figures](https://github.com/anaerovane/ReproduceDeepsea/tree/main/experiment/Figure3)

恢复到 fuxian3 项目根目录：

```bash
hf download aer0vane/reproduce_deepsea --include "experiment/Figure3/outputs/**" --local-dir .
```

目录约 28 GB，包含序列缓存、变异表、pretrained/ours 效应矩阵、10 折预测与分类器缓存、基线分数和结果核验文件。模型权重不在此目录重复存放。
