TFMS_vs_nonBlackList_random.py
- Calculates the Cluster Propensity (CP) for each transcription factor motif site (TFMS) using non-blacklist random regions as controls.
- CP quantifies the clustering tendency of TFMS in the genome using the Kolmogorov-Smirnov (KS) test statistic.
- Figure 1a demonstrates the clustering concept, while Figures 1b and 1c show specific examples of clustered TFMS associated with super-enhancers (SEs).

TFMS_enrich_at_SE.py
- Performs enrichment analysis of TF motif sites (TFMS) at super-enhancers (SEs) using Fisher's exact test.
- Compares TFMS between clustered and non-clustered categories to assess their association with SEs.

TFMS_CP_SE_plot.py
- Generates plots of TFMS Cluster Propensity (CP), enrichment at SEs, and log2 odds ratios.
- Visualizes the relationship between TFMS clustering and super-enhancer enrichment, aiding in the interpretation of transcription factor behavior in the genome, shown in Figure 1d