TFBS_vs_TFMS.py
- Calculates the Cluster Propensity (CP) for transcription factor binding sites (TFBS) using transcription factor motif sites (TFMS) as controls.
- CP quantifies the clustering tendency of TFBS in the genome relative to TFMS using the Kolmogorov-Smirnov (KS) test statistic.
- Figure 2a demonstrates the clustering of TFBS compared to random TFMS distributions.

combine_statistics_results.py
- Summarizes the statistical results of TFBS CP calculations, including enrichment analyses and correlations with super-enhancers.

summary_TFBS_CP_SE_plot.py
- Computes the TFBS Cluster Propensity (CP) and enrichment at super-enhancers (SEs) for each cell type.
- Generates heatmaps for visualizing TFBS CP and SE enrichment, shown in Figures 2b and 2c.
- Calculates the correlation between TFBS CP and SE enrichment and produces correlation plots, shown in Figure 2d.

TFBS_zscore_vs_CP.py
- Computes the Z-score for TFBS CP to normalize clustering tendencies across cell types.
- Plots the relationship between TFBS Z-scores and TFBS CP for each cell type shown in Figure 2e, highlighting TF binding clustering patterns.

