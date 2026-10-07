# HGMD Figure 3 source audit

## Target cohort

The DeepSEA paper describes the positive set as single nucleotide substitutions annotated as regulatory mutations in **HGMD Professional 2014.4**. Figure 3 reports `n = 2,977`. The paper's supplementary PDF contains figures and methods, but no per-variant HGMD coordinate table. The public Supplementary Tables 5 and 6 are the GRASP/eQTL and GWAS predictions, not HGMD positives.

## Historical public source lead (not an obtainable validated cohort)

Ensembl's March 2015 GRCh37 update states that it incorporated the latest public HGMD data, version 2014.4. Archived Ensembl Variation release 80 metadata identifies source 8 as `HGMD-PUBLIC`, version `20144` / December 2014, and an `HGMD-PUBLIC variants` set (set 29). This is only a lead to historical public annotations; it has not yielded a validated downloadable cohort for this project and is **not proof of equivalence to the paper's HGMD Professional 2014.4 cohort**.

## Why this does not yet reproduce the 2,977 positives

- HGMD Professional 2014.4 is the version explicitly named in the paper.
- The accessible Ensembl archive is HGMD-PUBLIC. The archived records expose genomic variation and consequence information, but the published material does not establish that this public set retains the professional release's `DM` classification and original regulatory mutation category.
- The archived Ensembl BioMart endpoint for that period timed out while exporting HGMD records. The historical FTP database dump is available, but its whole `variation_feature` table is several GB compressed and would still require checking whether the necessary HGMD category labels are present.
- A later paper also reports 2,977 HGMD pathogenic SNPs, but its underlying raw data are listed as available upon request; its count alone does not establish that its variants are the same as DeepSEA's.

## Current conclusion

The project has **not obtained a usable HGMD cohort**. I have not recovered or validated the exact 2,977 HGMD Professional regulatory SNVs, so it would be misleading to replace them with Ensembl HGMD-PUBLIC records, ClinVar, or another variant set and label it as the original HGMD cohort. Exact reconstruction likely requires the original licensed HGMD Professional release or the authors' cohort file.

## Sources

- DeepSEA paper: https://www.nature.com/articles/nmeth.3547
- Ensembl March 2015 GRCh37 update: https://www.ensembl.info/2015/03/31/first-update-of-the-ensembl-grch37-site/
- Ensembl release 80 source metadata: https://ftp.ensembl.org/pub/grch37/release-80/mysql/homo_sapiens_variation_80_37/source.txt.gz
- Ensembl release 80 variation-set metadata: https://ftp.ensembl.org/pub/grch37/release-80/mysql/homo_sapiens_variation_80_37/variation_set.txt.gz
- HGMD access model: https://pmc.ncbi.nlm.nih.gov/articles/PMC3898141/
