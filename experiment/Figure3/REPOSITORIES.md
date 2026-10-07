# Figure 3 文件位置

- GitHub： [代码、图片、方法说明和轻量结果摘要](https://github.com/anaerovane/ReproduceDeepsea/tree/main/experiment/Figure3)
- Hugging Face： [推理矩阵、逐样本预测和折缓存](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Figure3/outputs)

Hugging Face outputs 约 29.4 GB。恢复到 fuxian3 项目根目录：

```bash
hf download aer0vane/reproduce_deepsea --include "experiment/Figure3/outputs/**" --local-dir .
```

模型权重继续从已有 Hugging Face 权重路径按需下载；GitHub 不存放模型权重，也不存放大型推理缓存。
