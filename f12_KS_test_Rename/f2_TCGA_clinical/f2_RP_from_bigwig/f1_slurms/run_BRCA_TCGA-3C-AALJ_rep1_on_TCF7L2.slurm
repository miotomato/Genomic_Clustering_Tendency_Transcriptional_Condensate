#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurms/slurm_BRCA_TCGA-3C-AALJ_rep1_on_TCF7L2.out

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//MCF-7/MCF-7_TCF7L2_percentile_T.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/BRCA/BRCA_TCGA-3C-AALJ_rep1.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/BRCA/BRCA_TCGA-3C-AALJ_rep1_halflife_10000_on_MCF-7_TCF7L2_percentile_T.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//MCF-7/MCF-7_TCF7L2_percentile_T.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/BRCA/BRCA_TCGA-3C-AALJ_rep1.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/BRCA/BRCA_TCGA-3C-AALJ_rep1_halflife_10000_on_MCF-7_TCF7L2_percentile_T.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//MCF-7/MCF-7_TCF7L2_percentile_T_ExtendMerge.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/BRCA/BRCA_TCGA-3C-AALJ_rep1.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/BRCA/BRCA_TCGA-3C-AALJ_rep1_halflife_10000_on_MCF-7_TCF7L2_percentile_T_ExtendMerge.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//MCF-7/MCF-7_TCF7L2_percentile_T_ExtendMerge.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/BRCA/BRCA_TCGA-3C-AALJ_rep1.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/BRCA/BRCA_TCGA-3C-AALJ_rep1_halflife_10000_on_MCF-7_TCF7L2_percentile_T_ExtendMerge.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//MCF-7/MCF-7_TCF7L2_percentile_C.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/BRCA/BRCA_TCGA-3C-AALJ_rep1.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/BRCA/BRCA_TCGA-3C-AALJ_rep1_halflife_10000_on_MCF-7_TCF7L2_percentile_C.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//MCF-7/MCF-7_TCF7L2_percentile_C.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/BRCA/BRCA_TCGA-3C-AALJ_rep1.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/BRCA/BRCA_TCGA-3C-AALJ_rep1_halflife_10000_on_MCF-7_TCF7L2_percentile_C.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//MCF-7/MCF-7_TCF7L2_percentile_T.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/BRCA/BRCA_TCGA-3C-AALJ_rep1.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/BRCA/BRCA_TCGA-3C-AALJ_rep1_halflife_5000_on_MCF-7_TCF7L2_percentile_T.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//MCF-7/MCF-7_TCF7L2_percentile_T.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/BRCA/BRCA_TCGA-3C-AALJ_rep1.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/BRCA/BRCA_TCGA-3C-AALJ_rep1_halflife_5000_on_MCF-7_TCF7L2_percentile_T.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//MCF-7/MCF-7_TCF7L2_percentile_T_ExtendMerge.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/BRCA/BRCA_TCGA-3C-AALJ_rep1.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/BRCA/BRCA_TCGA-3C-AALJ_rep1_halflife_5000_on_MCF-7_TCF7L2_percentile_T_ExtendMerge.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//MCF-7/MCF-7_TCF7L2_percentile_T_ExtendMerge.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/BRCA/BRCA_TCGA-3C-AALJ_rep1.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/BRCA/BRCA_TCGA-3C-AALJ_rep1_halflife_5000_on_MCF-7_TCF7L2_percentile_T_ExtendMerge.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//MCF-7/MCF-7_TCF7L2_percentile_C.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/BRCA/BRCA_TCGA-3C-AALJ_rep1.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/BRCA/BRCA_TCGA-3C-AALJ_rep1_halflife_5000_on_MCF-7_TCF7L2_percentile_C.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//MCF-7/MCF-7_TCF7L2_percentile_C.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/BRCA/BRCA_TCGA-3C-AALJ_rep1.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/BRCA/BRCA_TCGA-3C-AALJ_rep1_halflife_5000_on_MCF-7_TCF7L2_percentile_C.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//MCF-7/MCF-7_TCF7L2_percentile_T.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/BRCA/BRCA_TCGA-3C-AALJ_rep1.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/BRCA/BRCA_TCGA-3C-AALJ_rep1_halflife_20000_on_MCF-7_TCF7L2_percentile_T.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//MCF-7/MCF-7_TCF7L2_percentile_T.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/BRCA/BRCA_TCGA-3C-AALJ_rep1.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/BRCA/BRCA_TCGA-3C-AALJ_rep1_halflife_20000_on_MCF-7_TCF7L2_percentile_T.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//MCF-7/MCF-7_TCF7L2_percentile_T_ExtendMerge.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/BRCA/BRCA_TCGA-3C-AALJ_rep1.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/BRCA/BRCA_TCGA-3C-AALJ_rep1_halflife_20000_on_MCF-7_TCF7L2_percentile_T_ExtendMerge.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//MCF-7/MCF-7_TCF7L2_percentile_T_ExtendMerge.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/BRCA/BRCA_TCGA-3C-AALJ_rep1.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/BRCA/BRCA_TCGA-3C-AALJ_rep1_halflife_20000_on_MCF-7_TCF7L2_percentile_T_ExtendMerge.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//MCF-7/MCF-7_TCF7L2_percentile_C.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/BRCA/BRCA_TCGA-3C-AALJ_rep1.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/BRCA/BRCA_TCGA-3C-AALJ_rep1_halflife_20000_on_MCF-7_TCF7L2_percentile_C.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//MCF-7/MCF-7_TCF7L2_percentile_C.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/BRCA/BRCA_TCGA-3C-AALJ_rep1.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/BRCA/BRCA_TCGA-3C-AALJ_rep1_halflife_20000_on_MCF-7_TCF7L2_percentile_C.merge.SE_overlapped.csv

