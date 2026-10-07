# HGMD 原始队列当前无法取得

## 结论

截至 2026-10-06，本项目没有取得可用于复现 DeepSEA Figure 3 的 HGMD 原始正例坐标。原论文所述队列为 **HGMD Professional 2014.4 的 2,977 个 regulatory SNV**。论文补充材料没有提供逐变异坐标表；本机也没有 HGMD Professional 授权或数据文件。

HGMD 官网说明，Public 版仅供注册的学术/非营利机构用户访问，且不能复制、保存或再分发；Professional 数据需授权。已检查的公开论文和补充材料未提供可确认等同于原图 2,977 个变异的队列。因此本实验不能声称复现 HGMD 面板。

另有 HGMD-PUBLIC/Ensembl 历史记录及其他人工整理的调控变异集线索，但它们未被证实与 HGMD Professional 2014.4 的 2,977 个 SNV 一致，不作为 HGMD 数据使用。ClinVar 正例也不是 HGMD 替代的等价数据，不能放进主图并标成 HGMD。

## 本项目 Figure 3 的实际范围

主复跑只展示目前能取得原始队列及论文补充预测表的两项：

- GRASP noncoding eQTL
- GWAS Catalog noncoding SNPs

主结果分成两张独立图，保存在 `figures/`：`figure3_pretrained_separate.png/.pdf` 和 `figure3_ours_separate.png/.pdf`。每张图含 GRASP eQTL 与 GWAS Catalog 两个可得面板，并比较 DeepSEA checkpoint、CADD PHRED、GWAVA unmatched/TSS/region 和 FunSeq2 v2.1.0。CADD 来源版本没有保留可复核的请求记录，图例明确标注“版本未核实”。由于缺少四项保守性特征，DeepSEA 曲线是 chromatin-effects-only，不能称为完整复现原图的 DeepSEA 方法表现。

历史上做过 ClinVar proxy 探索，但它不是 HGMD 等价队列，现已从 Figure3 文件夹删除其脚本、图片、结果和下载数据，避免与目标图混淆。

## 已核实的限制

- CADD PHRED eQTL/GWAS 每变异分数表保存在 `outputs/eqtl_CADD_scores.tsv.gz` 和 `outputs/gwas_CADD_scores.tsv.gz`；少量正例无 CADD 分数，主图比较时与其他方法共同取完整案例。该分数的 CADD 版本在本项目没有可复核记录。
- eQTL/GWAS 原始补充表的负标签中有与正例坐标完全相同的行；现有复跑沿用原表标签，结果需按该限制解释。
- 主图不呈现 HGMD 或 ClinVar 面板；CADD、GWAVA 和 FunSeq2 只作为可得 eQTL/GWAS 面板中的比较方法，不替代 HGMD 队列。

## 参考

- DeepSEA 原论文：https://www.nature.com/articles/nmeth.3547
- HGMD 官网与访问限制：https://www.hgmd.cf.ac.uk/
- 历史 HGMD-PUBLIC 来源审计：`HGMD_source_audit.md`
