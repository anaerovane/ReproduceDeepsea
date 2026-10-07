# Figure 3 数据核查与真实推理复跑

## 旧实验的结论

旧版复跑脚本已清理，不作为结果来源。当前主图代码入口为 `code/figure3_rerun.py`；队列定义、标签和空间交叉验证 fold 来自论文补充表 5/6。本地模型重新推理生成主图预测，不把补充表中的旧预测列当作本次模型结果。

这里的来源判断基于本地文件名、表内说明、样本数量和内容；已保存本地 SHA256。尝试与出版方下载文件做逐字节核验时，远端访问失败，未完成来源认证（见 `outputs/publisher_source_verification.json`）。不能把本次模型确实运行过，进一步推断为所有既有文件的历史下载来源均已独立验证。

但“有真实数据”不等于“正确复现 Figure 3”。旧实验存在以下问题：

- `extract_sequences` 读取 alt，却没有构建 alt 序列；模型只预测参考序列的 919 个染色质概率，分类器不是基于 ref/alt 变异效应。
- 只使用随机负样本，并按最多正样本数量的 5 倍抽样；没有评估论文的距离负样本组。
- 使用随机分层 5 折及树模型 XGBoost；论文使用空间 10 折和正则化线性 logistic 模型。
- 旧版曾把 ClinVar 当作 HGMD proxy；ClinVar 不是 HGMD 原始队列，旧 proxy 文件已清理。
- 旧 NPZ 只保存标签和预测，不含样本 ID、输入哈希、fold 或运行日志，无法仅凭该文件独立核实具体运行来源。未发现硬编码不等于已验证旧结果的完整来源。

旧 NPZ 的 pooled AUC 与当前空间折、负样本组和评估口径不一致，因此不作为 Figure 3 结果引用；当前结果仅以逐样本 OOF 文件重新计算的 AUC 汇总为准。

## 新实验的方法

论文：[Zhou & Troyanskaya, Nature Methods 2015](https://www.nature.com/articles/nmeth.3547)。详细方法还可检查本地 `../../nmeth.3547.pdf`。

1. GRASP/GWAS 使用补充表所有正样本和所有距离组/随机负样本，默认不抽样。合并两个任务的 SNP 后共 1,855,704 个唯一 chr/pos/ref/alt 组合。
2. 用 hg19 FASTA 取 1,000 bp 窗口，变异位于第 500 个碱基（数组索引 499）；严格检查 ref 和边界。无法验证 ref、有 N、非单碱基替换等样本排除并记录。
3. pretrained 使用 DeepSEA 发布的 sequence-prediction 权重和原始 Threshold 激活，`model_assets.py` 按需下载；ours 使用 Hugging Face 上 epoch-53 checkpoint（若 `training/checkpoints/` 有本地训练权重则优先用本地文件，否则自动取回）。
4. AGCT 编码。分别预测 ref/alt 和两个链，计算 919 个绝对概率差及 919 个 log-odds 差，然后平均两个链的效应。不能先平均预测再计算非线性效应。
5. eQTL/GWAS 沿用原始 fold 0–9。每折仅用其他 fold 的正样本和 All 随机负样本训练，再预测该 fold 的正样本和所有距离负样本。标准化只 fit 训练集；样本权重平衡正负类别。
6. XGBoost `gblinear`，eta 0.1，原始 alpha 0、lambda 10、100 轮。使用 GPU 确定性 coordinate descent；新版的 alpha/lambda 必须除以该折训练权重总和，才能保持旧版对梯度总和的正则强度。每折保存分类器、标准化参数及预测缓存。输出每个距离组与正样本的 pooled out-of-fold AUC。
7. 距离用负样本到同染色体最近正样本的实际 bp 距离计算，不硬编码横坐标。补充表的组名、Figure 3 的名义距离和按当前完整正样本集计算的最近距离有差异。例如 eQTL 的 31kbp/6.3kbp/710bp/360bp 组实际均值为约 14,781/3,838/1,251/401 bp；图中直接报告计算值，保留原始组名于 CSV。这些点不能声称与论文横坐标完全相同。

**局限：**缺少原文的四项进化保守性特征。DeepSEA 两个 checkpoint 的曲线属于 chromatin-effects-only 复跑，不能承诺得到原文完整 DeepSEA 曲线。图中的 CADD PHRED、三种 GWAVA 曲线和 FunSeq2 v2.1.0 曲线由对应变异分数重新计算；CADD 分数的版本没有可复核记录。HGMD 面板因原始队列不可得而省略。具体来源、哈希及 AUC 口径见 `FIGURE3_DATA_PROVENANCE.md`。

**原始补充表的标签重叠：**核查发现 eQTL 表有 16,377 行负样本与正样本的 chr/pos/ref/alt 完全相同，GWAS 表有 1,924 行。对应预测概率和 fold 也相同。本次保留作者 CSV 的原始标签，两个 checkpoint 使用完全相同的队列和空间 fold；同坐标样本不会跨 fold，但重复的正负标签会产生评价噪声。因此本次结果是“原始补充表标签上的模型比较”，不是独立重建并严格排除正负交集后的负样本实验。论文补充表参考 AUC 也按同样的原始标签计算，不能混用清理后的队列与未清理的参考值。详细核查见 `source_cohort_checks.json`。若另做清理队列实验，应移走所有分类器/OOF 缓存后重新训练，不能仅删掉测试行并沿用本次分类器。

## HGMD 队列不可得

HGMD Professional 2014.4 原始 2,977 个 regulatory SNV 坐标不可得，主图只展示 GRASP eQTL 与 GWAS Catalog。ClinVar proxy 不等价于 HGMD，相关旧脚本、结果和本地下载已清理；详情见 `HGMD_unavailable.md`。

## 代码与数据

- GitHub 代码、图片和说明：[ReproduceDeepsea / Figure3](https://github.com/anaerovane/ReproduceDeepsea/tree/main/experiment/Figure3)。
- Hugging Face 大型结果文件：[Figure3 outputs](https://huggingface.co/aer0vane/reproduce_deepsea/tree/main/experiment/Figure3/outputs)。需要本地数据时，在 fuxian3 根目录运行：`hf download aer0vane/reproduce_deepsea --include "experiment/Figure3/outputs/**" --local-dir .`。

## 运行与输出

目录一眼看懂：`code/` 放脚本，`figures/` 放最终图片，`outputs/` 放输入缓存、逐样本预测、模型折缓存和核验结果；本目录根部只放说明文档。

复跑需要 fuxian3 根目录下的论文补充表 5/6（`newdata/`）和 hg19 染色体 FASTA（`reference_genome/`）；下载位置与上游链接见 [项目下载指南](../../DOWNLOAD.md)。大型缓存和基线分数可从上面的 Hugging Face 地址恢复。

在 `Figure3/` 根目录执行主复跑（依赖 numpy/pandas/torch/sklearn/xgboost/matplotlib）：

```bash
python -u code/figure3_rerun.py > outputs/figure3_run.log 2>&1
```

数据来源、分组、完整案例限制和版本限制见 `FIGURE3_DATA_PROVENANCE.md`。

- `outputs/manifest.json`：原始补充表和 checkpoint 的 SHA256、是否抽样、方法及缺失数据。
- `outputs/variants.tsv`、`prepared.json`：唯一变异 ID、参考序列检查及排除原因。
- `outputs/*_effects.npy`：两个模型真实预测得到的 1,838 维效应；`*_progress.json` 支持中断恢复。
- `outputs/*_oof.tsv.gz`：逐样本坐标、标签、fold、变异 ID、距离及交叉验证预测。
- `outputs/auc_summary.csv`：完整队列的 DeepSEA OOF 分组 AUC；`outputs/pretrained_ours_baselines_auc.csv`：两张图所用共同完整案例上的 DeepSEA/CADD/GWAVA/FunSeq2 AUC、样本数和平均距离。
- `figures/figure3_pretrained_separate.png/pdf` 与 `figure3_ours_separate.png/pdf`：两张独立的 pretrained、ours 方法比较图。All 随机负例组不设数量上限；距离组也使用补充表提供的全部有效阴性样本。所有比较方法在每个 task/group 上使用相同的完整案例。

缓存只适用于生成它的输入和代码。修改算法后应指定新输出目录或先移走旧缓存；不可直接把旧预测当新算法的结果。 本次整理只移动文件并修正路径引用，没有重新训练或评估；`outputs/metrics_provenance.json` 中的代码哈希仍对应生成这些结果时的版本。

## 当前可用队列结果与核验

主图所有 AUC、样本数和实际平均距离均由逐样本 OOF 预测及输入队列计算，不在绘图代码中手工填写。当前数值见 `outputs/pretrained_ours_baselines_auc.csv`，来源哈希及计算说明见 `outputs/metrics_provenance.json` 和 `FIGURE3_DATA_PROVENANCE.md`。

eQTL/GWAS 的 40 个分类器折均已完成。逐样本 OOF 文件重算的 AUC 与保存汇总一致，同坐标始终同 fold。两个模型全部 1,855,581 个有效 SNP 的 1,838 维效应已扫描，全部有限，919 个绝对效应全部处于 [0, 1]。另已核对 pretrained 与本地原始 variant 模型在测试序列上的输出一致；相同 ref/alt 的效应在浮点误差内为零。

HGMD 数据可得性结论见 `HGMD_unavailable.md`。主图只含 eQTL/GWAS 两项。核验记录位于 `outputs/metrics_provenance.json`、`eqtl_verification.json`、`gwas_verification.json` 和 `feature_verification.json`。
