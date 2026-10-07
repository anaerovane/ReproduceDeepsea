# Supple3 核查与运行说明

本次使用本地两个权重重新完成四个位点、每个位点 3000 个单碱基突变的推理。位点、等位基因、特征索引、权重路径和序列设置集中在 `experiment6_config.json`；绘图在使用每个特征前，会核对 `predictor_names.txt` 中对应的精确名称。
`prepare_sequences.py` 从本地 GRCh37/hg19 FASTA 生成四条 1000bp 输入，并核对参考碱基及歧义字符。新生成的四个 one-hot 数组与最终推理使用的旧输入逐元素完全一致，中心参考碱基为 G/T/A/C。
重构后的推理入口成功严格加载两个 checkpoint 并完成全部 24,000 个突变预测。八组效应均由保存的原始概率按 log2-odds 公式复算，最大误差约 `6.2e-7`（浮点舍入）。
旧 `mayexperiment/experiment6_all_data.npz` 中 ours 数组与最终结果完全一致，但 pretrained 数组最大差异为 6.52，属于旧的混合缓存，已删除；旧单序列图和旧脚本也已删除。
这证明结果来自所提供权重的推理；不独立证明权重训练过程或参考基因组版本符合原论文。

修复并固化：移除旧的单序列测试集例子与旧混合缓存；统一变异与模型参数配置；GATA1/GATA-1 名称按明确特征索引和预期名称校验；pretrained 输入通道在 AGCT 突变生成后按配置转换；序列输入及中心参考碱基均校验。对同名重复特征使用配置中记录的索引，不按效应大小挑选。

颜色为负值蓝色、零白色、正值黄色。同一 TF 两模型与全长/放大图共用对称色标，范围来自实际数值。
四行分别为 CEBPB/HepG2、GATA1/K562、FOXA1/HepG2、FOXA2/HepG2；每行上 1000bp，下 400–600bp，底部为真实序列。热图行顺序 AGCT，参考碱基效应为零。
这些是模型预测的 in silico 效应，不是实验测量；数值与条纹未调整为参考图片。未核对原论文具体 replicate 和效应定义，不能声称逐像素复现。

运行（在项目根目录；先准备序列）：

```bash
python deepsea/fuxian3/experiment/Supple3/prepare_sequences.py
python deepsea/fuxian3/experiment/Supple3/experiment6_saturation_full.py
# 仅从最终效应缓存重绘
python deepsea/fuxian3/experiment/Supple3/experiment6_plot.py
```

数据：experiment6_all_data.npz 保存效应；experiment6_raw_predictions.npz 保存原始和突变预测；两张图各输出 PNG 和 PDF。
