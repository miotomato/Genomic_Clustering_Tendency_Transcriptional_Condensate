#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurms/slurm_PRAD_TCGA-YL-A9WK_rep1_on_HNF4G.out

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//LNCaP/LNCaP_HNF4G_percentile_T.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/PRAD/PRAD_TCGA-YL-A9WK_rep1.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/PRAD/PRAD_TCGA-YL-A9WK_rep1_halflife_10000_on_LNCaP_HNF4G_percentile_T.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//LNCaP/LNCaP_HNF4G_percentile_T.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/PRAD/PRAD_TCGA-YL-A9WK_rep1.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/PRAD/PRAD_TCGA-YL-A9WK_rep1_halflife_10000_on_LNCaP_HNF4G_percentile_T.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//LNCaP/LNCaP_HNF4G_percentile_T_ExtendMerge.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/PRAD/PRAD_TCGA-YL-A9WK_rep1.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/PRAD/PRAD_TCGA-YL-A9WK_rep1_halflife_10000_on_LNCaP_HNF4G_percentile_T_ExtendMerge.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//LNCaP/LNCaP_HNF4G_percentile_T_ExtendMerge.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/PRAD/PRAD_TCGA-YL-A9WK_rep1.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/PRAD/PRAD_TCGA-YL-A9WK_rep1_halflife_10000_on_LNCaP_HNF4G_percentile_T_ExtendMerge.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//LNCaP/LNCaP_HNF4G_percentile_C.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/PRAD/PRAD_TCGA-YL-A9WK_rep1.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/PRAD/PRAD_TCGA-YL-A9WK_rep1_halflife_10000_on_LNCaP_HNF4G_percentile_C.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//LNCaP/LNCaP_HNF4G_percentile_C.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/PRAD/PRAD_TCGA-YL-A9WK_rep1.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/PRAD/PRAD_TCGA-YL-A9WK_rep1_halflife_10000_on_LNCaP_HNF4G_percentile_C.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//LNCaP/LNCaP_HNF4G_percentile_T.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/PRAD/PRAD_TCGA-YL-A9WK_rep1.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/PRAD/PRAD_TCGA-YL-A9WK_rep1_halflife_5000_on_LNCaP_HNF4G_percentile_T.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//LNCaP/LNCaP_HNF4G_percentile_T.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/PRAD/PRAD_TCGA-YL-A9WK_rep1.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/PRAD/PRAD_TCGA-YL-A9WK_rep1_halflife_5000_on_LNCaP_HNF4G_percentile_T.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//LNCaP/LNCaP_HNF4G_percentile_T_ExtendMerge.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/PRAD/PRAD_TCGA-YL-A9WK_rep1.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/PRAD/PRAD_TCGA-YL-A9WK_rep1_halflife_5000_on_LNCaP_HNF4G_percentile_T_ExtendMerge.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//LNCaP/LNCaP_HNF4G_percentile_T_ExtendMerge.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/PRAD/PRAD_TCGA-YL-A9WK_rep1.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/PRAD/PRAD_TCGA-YL-A9WK_rep1_halflife_5000_on_LNCaP_HNF4G_percentile_T_ExtendMerge.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//LNCaP/LNCaP_HNF4G_percentile_C.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/PRAD/PRAD_TCGA-YL-A9WK_rep1.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/PRAD/PRAD_TCGA-YL-A9WK_rep1_halflife_5000_on_LNCaP_HNF4G_percentile_C.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//LNCaP/LNCaP_HNF4G_percentile_C.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/PRAD/PRAD_TCGA-YL-A9WK_rep1.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/PRAD/PRAD_TCGA-YL-A9WK_rep1_halflife_5000_on_LNCaP_HNF4G_percentile_C.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//LNCaP/LNCaP_HNF4G_percentile_T.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/PRAD/PRAD_TCGA-YL-A9WK_rep1.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/PRAD/PRAD_TCGA-YL-A9WK_rep1_halflife_20000_on_LNCaP_HNF4G_percentile_T.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//LNCaP/LNCaP_HNF4G_percentile_T.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/PRAD/PRAD_TCGA-YL-A9WK_rep1.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/PRAD/PRAD_TCGA-YL-A9WK_rep1_halflife_20000_on_LNCaP_HNF4G_percentile_T.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//LNCaP/LNCaP_HNF4G_percentile_T_ExtendMerge.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/PRAD/PRAD_TCGA-YL-A9WK_rep1.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/PRAD/PRAD_TCGA-YL-A9WK_rep1_halflife_20000_on_LNCaP_HNF4G_percentile_T_ExtendMerge.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//LNCaP/LNCaP_HNF4G_percentile_T_ExtendMerge.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/PRAD/PRAD_TCGA-YL-A9WK_rep1.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/PRAD/PRAD_TCGA-YL-A9WK_rep1_halflife_20000_on_LNCaP_HNF4G_percentile_T_ExtendMerge.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//LNCaP/LNCaP_HNF4G_percentile_C.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/PRAD/PRAD_TCGA-YL-A9WK_rep1.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/PRAD/PRAD_TCGA-YL-A9WK_rep1_halflife_20000_on_LNCaP_HNF4G_percentile_C.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//LNCaP/LNCaP_HNF4G_percentile_C.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/PRAD/PRAD_TCGA-YL-A9WK_rep1.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/PRAD/PRAD_TCGA-YL-A9WK_rep1_halflife_20000_on_LNCaP_HNF4G_percentile_C.merge.SE_overlapped.csv

