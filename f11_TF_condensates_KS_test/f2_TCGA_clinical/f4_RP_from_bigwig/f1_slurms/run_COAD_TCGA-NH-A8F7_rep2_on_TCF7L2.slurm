#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurms/slurm_COAD_TCGA-NH-A8F7_rep2_on_TCF7L2.out

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged//HCT-116/HCT-116_TCF7L2_percentile_T.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/COAD/COAD_TCGA-NH-A8F7_rep2.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/COAD/COAD_TCGA-NH-A8F7_rep2_halflife_10000_on_HCT-116_TCF7L2_percentile_T.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged//HCT-116/HCT-116_TCF7L2_percentile_T.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/COAD/COAD_TCGA-NH-A8F7_rep2.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/COAD/COAD_TCGA-NH-A8F7_rep2_halflife_10000_on_HCT-116_TCF7L2_percentile_T.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged//HCT-116/HCT-116_TCF7L2_percentile_T_ExtendMerge.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/COAD/COAD_TCGA-NH-A8F7_rep2.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/COAD/COAD_TCGA-NH-A8F7_rep2_halflife_10000_on_HCT-116_TCF7L2_percentile_T_ExtendMerge.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged//HCT-116/HCT-116_TCF7L2_percentile_T_ExtendMerge.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/COAD/COAD_TCGA-NH-A8F7_rep2.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/COAD/COAD_TCGA-NH-A8F7_rep2_halflife_10000_on_HCT-116_TCF7L2_percentile_T_ExtendMerge.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged//HCT-116/HCT-116_TCF7L2_percentile_C.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/COAD/COAD_TCGA-NH-A8F7_rep2.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/COAD/COAD_TCGA-NH-A8F7_rep2_halflife_10000_on_HCT-116_TCF7L2_percentile_C.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged//HCT-116/HCT-116_TCF7L2_percentile_C.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/COAD/COAD_TCGA-NH-A8F7_rep2.bdg \
-s hg38 -u 10000 -o f1_get_RP_from_bedGraph/COAD/COAD_TCGA-NH-A8F7_rep2_halflife_10000_on_HCT-116_TCF7L2_percentile_C.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged//HCT-116/HCT-116_TCF7L2_percentile_T.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/COAD/COAD_TCGA-NH-A8F7_rep2.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/COAD/COAD_TCGA-NH-A8F7_rep2_halflife_5000_on_HCT-116_TCF7L2_percentile_T.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged//HCT-116/HCT-116_TCF7L2_percentile_T.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/COAD/COAD_TCGA-NH-A8F7_rep2.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/COAD/COAD_TCGA-NH-A8F7_rep2_halflife_5000_on_HCT-116_TCF7L2_percentile_T.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged//HCT-116/HCT-116_TCF7L2_percentile_T_ExtendMerge.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/COAD/COAD_TCGA-NH-A8F7_rep2.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/COAD/COAD_TCGA-NH-A8F7_rep2_halflife_5000_on_HCT-116_TCF7L2_percentile_T_ExtendMerge.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged//HCT-116/HCT-116_TCF7L2_percentile_T_ExtendMerge.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/COAD/COAD_TCGA-NH-A8F7_rep2.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/COAD/COAD_TCGA-NH-A8F7_rep2_halflife_5000_on_HCT-116_TCF7L2_percentile_T_ExtendMerge.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged//HCT-116/HCT-116_TCF7L2_percentile_C.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/COAD/COAD_TCGA-NH-A8F7_rep2.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/COAD/COAD_TCGA-NH-A8F7_rep2_halflife_5000_on_HCT-116_TCF7L2_percentile_C.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged//HCT-116/HCT-116_TCF7L2_percentile_C.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/COAD/COAD_TCGA-NH-A8F7_rep2.bdg \
-s hg38 -u 5000 -o f1_get_RP_from_bedGraph/COAD/COAD_TCGA-NH-A8F7_rep2_halflife_5000_on_HCT-116_TCF7L2_percentile_C.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged//HCT-116/HCT-116_TCF7L2_percentile_T.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/COAD/COAD_TCGA-NH-A8F7_rep2.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/COAD/COAD_TCGA-NH-A8F7_rep2_halflife_20000_on_HCT-116_TCF7L2_percentile_T.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged//HCT-116/HCT-116_TCF7L2_percentile_T.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/COAD/COAD_TCGA-NH-A8F7_rep2.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/COAD/COAD_TCGA-NH-A8F7_rep2_halflife_20000_on_HCT-116_TCF7L2_percentile_T.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged//HCT-116/HCT-116_TCF7L2_percentile_T_ExtendMerge.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/COAD/COAD_TCGA-NH-A8F7_rep2.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/COAD/COAD_TCGA-NH-A8F7_rep2_halflife_20000_on_HCT-116_TCF7L2_percentile_T_ExtendMerge.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged//HCT-116/HCT-116_TCF7L2_percentile_T_ExtendMerge.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/COAD/COAD_TCGA-NH-A8F7_rep2.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/COAD/COAD_TCGA-NH-A8F7_rep2_halflife_20000_on_HCT-116_TCF7L2_percentile_T_ExtendMerge.merge.SE_overlapped.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged//HCT-116/HCT-116_TCF7L2_percentile_C.merge.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/COAD/COAD_TCGA-NH-A8F7_rep2.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/COAD/COAD_TCGA-NH-A8F7_rep2_halflife_20000_on_HCT-116_TCF7L2_percentile_C.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged//HCT-116/HCT-116_TCF7L2_percentile_C.merge.SE_overlapped.bed \
-b ../../../f9_TF_condensates_V3/data/TCGA/bigwig/COAD/COAD_TCGA-NH-A8F7_rep2.bdg \
-s hg38 -u 20000 -o f1_get_RP_from_bedGraph/COAD/COAD_TCGA-NH-A8F7_rep2_halflife_20000_on_HCT-116_TCF7L2_percentile_C.merge.SE_overlapped.csv

