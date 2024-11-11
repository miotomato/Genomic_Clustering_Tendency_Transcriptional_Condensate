# Transciptional_condensates



### f6_TCGA_atacseq_survival

TCGA survival top/bottom 50% patients -> differential ATAC score -> BART


### f7_TF_condensates_test

-- Hypothesis test --
TF clusters/condensates are associated with SE
check if the clustered TFs, 
either from motif matching data, or chip-seq data in a cell type,
have any association/correlation with IDR and SE

NOTES: focus on cell-type specific, or differential ATAC-seq regions


### f8_TF_condensates_V2

NOTES: use TFBS data from all cell types 

cluster potential (CP) of 
- TF motif sites (TFMS) 
- ChIP-seq TF binding sites (TFBS) with motif
- TFBS without motif
 
Correlation among
- TF cluster potential
- TF protein AICAP score
- TF enrichment at SEs

Clinical outcomes by
- ATAC-seq signals at cancer-specific TFs
- ATAC-seq signals at SEs with cancer-specific TFs  




### f9_TF_condensates_V3

NOTES: only look at distance to down-stream region

cluster potential (CP) of 
- TF motif sites (TFMS) vs. random
- TFBS vs. TFMS

only look at distance to down-stream region in bedtools closest for
- TFMS, TFBS, and 
- randomly selected TFMS as control for TFBS




### f10_TF_condensates_top10k_TFBS


NOTES: 
only use top10K peaks from each dataset
only look at distance to down-stream region



### f11_TF_condensates_KS_test

NOTES: 
- only look at distance to down-stream region
- use KS test for TFMS/TFBS CP



### f12_KS_test_with_rename

- only look at distance to down-stream region
- use KS test for TFMS/TFBS CP
- Jurkat and CUTLL1 as T-ALL
- RBPJ as NOTCH1











