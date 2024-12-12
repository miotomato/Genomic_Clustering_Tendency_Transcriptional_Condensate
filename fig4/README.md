clustered_TFBS_merge.py
- Identifies clustered TF binding sites (C-TFBS) and non-clustered TF binding sites (NC-TFBS) by merging neighboring TFBS within a defined distance.
- Demonstration plot shown in Figure 4a.

ATAC_ov_TFBS.py
- Analyzes the overlap between ATAC-seq peaks and TF binding sites (TFBS) to assess chromatin accessibility.

ATAC_RP.py
find_overlap_keep_info_NOT_sep_strand_revised.py
get-regulatory-potential-on-genes_peak_level.py
- Calculates the regulatory potential (RP) of ATAC-seq signals at TFBS.

summarize_ATAC_HiC_plot.py
- Summarizes chromatin accessibility (ATAC) signals and regulatory potential (RP) at C-TFBS and NC-TFBS.
- Performs Fisher's exact test to compare enrichment at super-enhancers (SEs), chromatin accessibility, differential ATAC signals, and chromatin interactions between C-TFBS and NC-TFBS.
- Demonstration plot shown in Figure 4a, example plots shown in Figures 4b (BRCA) and 4c (COAD).

pancancer_mutation_plot.py
- Performs a t-test to compare mutation rates between C-TFBS and NC-TFBS.
- Results shown in Figures 4d and 4e.