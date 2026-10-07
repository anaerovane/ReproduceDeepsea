# Figure 3 复跑：数据来源与数字计算

本文档是 Figure 3 主图的来源说明。当前主图只有两张：pretrained 与 ours；每张各有 GRASP eQTL、GWAS Catalog 两个面板。HGMD 原始队列没有取得，因此不画 HGMD 面板；ClinVar proxy 已从主目录移除，不替代 HGMD。

## 图和结果表

- `figures/figure3_pretrained_separate.png` / `.pdf`
- `figures/figure3_ours_separate.png` / `.pdf`
- `outputs/pretrained_ours_baselines_auc.csv`：两张图共同使用的 AUC、样本数和实际平均距离，每个队列、模型、负例组一行。
- `outputs/metrics_provenance.json`：输入与结果 SHA256。
- `outputs/manifest.json`、`baseline_source_manifest.json`：模型和输入队列、基线来源、版本及匹配规则。

## 输入和来源

| 用途 | 输入/来源 | 处理和可核对记录 |
|---|---|---|
| eQTL 队列 | 论文 Supplementary Table 5：`newdata/41592_2015_BFnmeth3547_MOESM649_ESM.csv`；[论文页面](https://www.nature.com/articles/nmeth.3547)；[Supplementary Table 5 CSV](https://media-springernature-com.libproxy1.nus.edu.sg/original/springer-static/esm/art%3A10.1038%2Fnmeth.3547/MediaObjects/41592_2015_BFnmeth3547_MOESM649_ESM.csv) | 采用表内 GRASP eQTL 正例、随机 1000 Genomes SNP 负例及按距离分组的负例。原表标签、fold、坐标和 ref/alt 均保留。SHA256 在 `manifest.json`。 |
| GWAS 队列 | 论文 Supplementary Table 6：`newdata/41592_2015_BFnmeth3547_MOESM650_ESM.csv`；[论文页面](https://www.nature.com/articles/nmeth.3547)；[Supplementary Table 6 CSV](https://media-springernature-com.libproxy1.nus.edu.sg/original/springer-static/esm/art%3A10.1038%2Fnmeth.3547/MediaObjects/41592_2015_BFnmeth3547_MOESM650_ESM.csv) | 采用表内 GWAS Catalog 正例和对应负例组。处理与 eQTL 相同。SHA256 在 `manifest.json`。 |
| 参考基因组 | `reference_genome/chr1.fa` … `chr22.fa`、`chrX.fa`、`chrY.fa` | 以 hg19/GRCh37 坐标截取 1,000 bp 上下文并核对 ref 碱基。所有使用 fasta 的 SHA256 记录于 `manifest.json` 和 `metrics_provenance.json`。 |
| pretrained DeepSEA | DeepSEA sequence-prediction 权重，由 `model_assets.py` 从 [Hugging Face](https://huggingface.co/aer0vane/reproduce_deepsea/blob/main/models/deepsea_predict.pth) 或 [Zenodo](https://zenodo.org/records/1466993) 取回 | 输入变异两侧各 500 bp，核对 hg19/GRCh37 ref；对 ref/alt、正反链计算 919 个输出的效应。文件 SHA256 在 `manifest.json`。 |
| ours DeepSEA | Hugging Face `training_checkpoints/best_model_FINAL_EPOCH53.pth` | 使用相同的变异、序列、参考检查与空间 fold；文件 SHA256 在 `manifest.json`。 |
| CADD | `eqtl_CADD_scores.tsv.gz`、`gwas_CADD_scores.tsv.gz` | 图和 AUC 使用文件中的 `CADD_PHRED` 每变异分数。旧文件名曾写 `v1.0`，但本项目没有保留可核对的下载/请求记录，**版本不能确认**，所以图例只标“CADD PHRED，版本未核实”。分数表 SHA256 在 `metrics_provenance.json`。 |
| GWAVA | Sanger [GWAVA v1.0](https://ftp.sanger.ac.uk/pub/resources/software/gwava/v1.0/)；本次可访问的 [GRCh37 索引镜像](https://zenodo.org/records/3956168) | 使用镜像 BED 的三个分数字段，来源字段次序映射为 Region、TSS、Unmatched。BED 无 allele 字段，按 GRCh37 chr/position 匹配；同一位置若有冲突分数则排除。原始 694 MB 文件和索引已清理，下载文件与索引的哈希、队列分数表哈希保存在 `baseline_source_manifest.json`。 |
| FunSeq2 | [FunSeq2 building data context](https://info.gersteinlab.org/Funseq2#Building_data_context)；原版 hg19 v2.1.0 全基因组评分轨道 `http://archive.gersteinlab.org/funseq2/hg19_wg_score.tsv.gz` 及同名 `.tbi` | 下载完整索引轨道后按 GRCh37 chr/position/ref/alt 提取，优先用等位基因精确匹配；只有没有 allele-specific 记录时才回退到 alt `.` 的通用分数。全量轨道已清理。轨道/索引哈希、字节数和提取表哈希在 `baseline_source_manifest.json`。 |

## 图中数字怎样产生

1. DeepSEA 两个 checkpoint 的每变异效应输入由实际 ref/alt 模型推理产生；两条 DeepSEA 曲线的分类器采用补充表给定的 10 个空间 fold。fold 外样本的预测写入 `*_oof.tsv.gz`，AUC 从该逐样本 OOF 概率计算。
2. 每个任务、每个负例组单独计算 ROC AUC：该组负例为 0、该任务正例为 1。CADD、GWAVA 和 FunSeq2 使用对应变异分数；所有 AUC 均由数据和代码实时计算，绘图脚本没有手填 AUC 数组。
3. 为保证方法之间可比，每个 task/group 先取 DeepSEA、CADD、三条 GWAVA、FunSeq2 都有有效分数的共同完整样本，再对 pretrained 与 ours 使用这同一批样本。
4. 当前运行清单记录 `random_negative_cap=0`，即没有对补充表中的随机负例设置数量上限；正例和距离负例组也保留。原文 caption 只说随机抽取 1000 Genomes 负例，没有说明从更大池中抽取的比例。横坐标距离根据保留负例到同染色体最近正例的 bp 距离计算，显示组内实际均值；“All” 显示为 All。
5. 面板标题的 `n` 和结果 CSV 中的样本数是完成全部方法共同完整案例筛选后的正例数/负例数，不是原补充表的原始行数。

## 结果边界

- 缺少原文四项保守性特征，因此 DeepSEA 曲线只包含 chromatin effects，不能称作完整特征版 DeepSEA。
- HGMD Professional 原始 2,977 个 regulatory SNV 坐标不可得；ClinVar 不是等价替代，未放进主图。
- Supplementary Tables 5/6 本身存在正负样本完全相同的 chr/pos/ref/alt 行。本复跑保留原表标签，且在方法间使用相同队列；这会带来标签重叠噪声，详见本目录 `README.md` 中的原始补充表标签重叠说明。
- 当前 “All” 组不设置随机负例数量上限；结果范围仍受论文补充表所提供的候选队列限制，不能称为对完整 1000 Genomes SNP 集的重建。
- CADD 原始分数表的来源版本缺乏可复核请求记录；对外报告时应按图例写“CADD PHRED（版本未核实）”，不要引用旧文件名推断版本。

## 正确代码入口

- `code/figure3_rerun.py`：真实 ref/alt 推理、10-fold OOF 评估、动态 AUC 表与两张图。
- `code/deepsea_models.py`：两种 checkpoint 的 DeepSEA 结构与加载器。
- `code/retrieve_gwava_scores.py`、`code/retrieve_funseq2_scores.py`：从带索引的来源轨道提取队列位点分数。

运行主复跑：`python -u code/figure3_rerun.py`。现有缓存位于 `outputs/`，按 manifest 和输入哈希使用；变更输入/算法时应指定新的 `--output` 目录。基线轨道的原始大文件已清理，因此从头重新提取 GWAVA/FunSeq2 分数时，需要按 `outputs/baseline_source_manifest.json` 的来源重新获取对应轨道和 tabix 索引。
