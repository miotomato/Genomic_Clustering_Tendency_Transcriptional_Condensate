cobinding_TFBS_annotation.py
- Identifies TF binding sites (TFBS) with the highest Cluster Propensity (CP).
- Identifies co-binding TFBS and adds genomic annotations.
- Results shown in Figures 5a and 5b (COAD) and Figures 6a and 6b (BRCA).

DCI_ov_cobinding_TFBS.py
DCI_ov_cobinding_TFBS_plot.py
- Calculates differential chromatin interactions (DCI) using public Hi-C data across different subsets of co-binding TFBS.
- Plot results shown in Figure 5c.

ATAC_ov_cobinding_TFBS.py
ATAC_ov_cobinding_TFBS_plot.py
- Identifies TFBS overlapping with ATAC-seq peaks associated with patient survival.
- Results shown in Figure 5d.

ATAC_ov_cobinding_TFBS_target_gene.py
- Identifies target genes of TFBS-overlapping ATAC-seq peaks associated with patient survival.
- Demonstration plot shown in Figure 6c, with results for different subsets shown in Figure 6d.

openness_survival_analysis.py
- Performs patient survival analysis using chromatin accessibility data.
- Example shown in Figure 5e, with the corresponding IGV snapshot shown in Figure 5f.