#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurms/slurm_CESC_TCGA-EA-A3HT_rep5_on_ZNF143.out

time python get_RP_from_bedGraph.py \
-a ../../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_ZNF143_percentile_T.merge.bed \
-b ../../../data/TCGA/bigwig/CESC/CESC_TCGA-EA-A3HT_rep5.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/CESC/CESC_TCGA-EA-A3HT_rep5_halflife_10000_on_HeLa_ZNF143_percentile_T.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_ZNF143_percentile_T.merge.SE_overlapped.bed \
-b ../../../data/TCGA/bigwig/CESC/CESC_TCGA-EA-A3HT_rep5.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/CESC/CESC_TCGA-EA-A3HT_rep5_halflife_10000_on_HeLa_ZNF143_percentile_T.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_ZNF143_percentile_T_ExtendMerge.merge.bed \
-b ../../../data/TCGA/bigwig/CESC/CESC_TCGA-EA-A3HT_rep5.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/CESC/CESC_TCGA-EA-A3HT_rep5_halflife_10000_on_HeLa_ZNF143_percentile_T_ExtendMerge.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_ZNF143_percentile_T_ExtendMerge.merge.SE_overlapped.bed \
-b ../../../data/TCGA/bigwig/CESC/CESC_TCGA-EA-A3HT_rep5.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/CESC/CESC_TCGA-EA-A3HT_rep5_halflife_10000_on_HeLa_ZNF143_percentile_T_ExtendMerge.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_ZNF143_percentile_C.merge.bed \
-b ../../../data/TCGA/bigwig/CESC/CESC_TCGA-EA-A3HT_rep5.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/CESC/CESC_TCGA-EA-A3HT_rep5_halflife_10000_on_HeLa_ZNF143_percentile_C.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_ZNF143_percentile_C.merge.SE_overlapped.bed \
-b ../../../data/TCGA/bigwig/CESC/CESC_TCGA-EA-A3HT_rep5.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/CESC/CESC_TCGA-EA-A3HT_rep5_halflife_10000_on_HeLa_ZNF143_percentile_C.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_ZNF143_percentile_T.merge.bed \
-b ../../../data/TCGA/bigwig/CESC/CESC_TCGA-EA-A3HT_rep5.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/CESC/CESC_TCGA-EA-A3HT_rep5_halflife_5000_on_HeLa_ZNF143_percentile_T.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_ZNF143_percentile_T.merge.SE_overlapped.bed \
-b ../../../data/TCGA/bigwig/CESC/CESC_TCGA-EA-A3HT_rep5.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/CESC/CESC_TCGA-EA-A3HT_rep5_halflife_5000_on_HeLa_ZNF143_percentile_T.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_ZNF143_percentile_T_ExtendMerge.merge.bed \
-b ../../../data/TCGA/bigwig/CESC/CESC_TCGA-EA-A3HT_rep5.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/CESC/CESC_TCGA-EA-A3HT_rep5_halflife_5000_on_HeLa_ZNF143_percentile_T_ExtendMerge.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_ZNF143_percentile_T_ExtendMerge.merge.SE_overlapped.bed \
-b ../../../data/TCGA/bigwig/CESC/CESC_TCGA-EA-A3HT_rep5.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/CESC/CESC_TCGA-EA-A3HT_rep5_halflife_5000_on_HeLa_ZNF143_percentile_T_ExtendMerge.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_ZNF143_percentile_C.merge.bed \
-b ../../../data/TCGA/bigwig/CESC/CESC_TCGA-EA-A3HT_rep5.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/CESC/CESC_TCGA-EA-A3HT_rep5_halflife_5000_on_HeLa_ZNF143_percentile_C.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_ZNF143_percentile_C.merge.SE_overlapped.bed \
-b ../../../data/TCGA/bigwig/CESC/CESC_TCGA-EA-A3HT_rep5.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/CESC/CESC_TCGA-EA-A3HT_rep5_halflife_5000_on_HeLa_ZNF143_percentile_C.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_ZNF143_percentile_T.merge.bed \
-b ../../../data/TCGA/bigwig/CESC/CESC_TCGA-EA-A3HT_rep5.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/CESC/CESC_TCGA-EA-A3HT_rep5_halflife_20000_on_HeLa_ZNF143_percentile_T.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_ZNF143_percentile_T.merge.SE_overlapped.bed \
-b ../../../data/TCGA/bigwig/CESC/CESC_TCGA-EA-A3HT_rep5.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/CESC/CESC_TCGA-EA-A3HT_rep5_halflife_20000_on_HeLa_ZNF143_percentile_T.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_ZNF143_percentile_T_ExtendMerge.merge.bed \
-b ../../../data/TCGA/bigwig/CESC/CESC_TCGA-EA-A3HT_rep5.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/CESC/CESC_TCGA-EA-A3HT_rep5_halflife_20000_on_HeLa_ZNF143_percentile_T_ExtendMerge.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_ZNF143_percentile_T_ExtendMerge.merge.SE_overlapped.bed \
-b ../../../data/TCGA/bigwig/CESC/CESC_TCGA-EA-A3HT_rep5.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/CESC/CESC_TCGA-EA-A3HT_rep5_halflife_20000_on_HeLa_ZNF143_percentile_T_ExtendMerge.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_ZNF143_percentile_C.merge.bed \
-b ../../../data/TCGA/bigwig/CESC/CESC_TCGA-EA-A3HT_rep5.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/CESC/CESC_TCGA-EA-A3HT_rep5_halflife_20000_on_HeLa_ZNF143_percentile_C.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_ZNF143_percentile_C.merge.SE_overlapped.bed \
-b ../../../data/TCGA/bigwig/CESC/CESC_TCGA-EA-A3HT_rep5.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/CESC/CESC_TCGA-EA-A3HT_rep5_halflife_20000_on_HeLa_ZNF143_percentile_C.merge.SE_overlapped.csv

