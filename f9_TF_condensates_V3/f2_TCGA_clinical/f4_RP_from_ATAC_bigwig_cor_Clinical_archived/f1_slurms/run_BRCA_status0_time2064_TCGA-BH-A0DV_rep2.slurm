#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurms/slurm_BRCA_status0_time2064_TCGA-BH-A0DV_rep2.out

time python get_RP_from_bedGraph.py \
-a ../../data/TCGA/tcga_atac.bed \
-b ../../data/TCGA/bigwig/BRCA_status0_time2064_TCGA-BH-A0DV_rep2.bdg \
-s hg38 -o f1_get_RP_from_bedGraph/BRCA_status0_time2064_TCGA-BH-A0DV_rep2_on_tcga_atac.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merge/HCT-116.merge.bed \
-b ../../data/TCGA/bigwig/BRCA_status0_time2064_TCGA-BH-A0DV_rep2.bdg \
-s hg38 -o f1_get_RP_from_bedGraph/BRCA_status0_time2064_TCGA-BH-A0DV_rep2_on_HCT-116.merge.csv

time python get_RP_from_bedGraph.py \
-a ../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merge/MCF-7.merge.bed \
-b ../../data/TCGA/bigwig/BRCA_status0_time2064_TCGA-BH-A0DV_rep2.bdg \
-s hg38 -o f1_get_RP_from_bedGraph/BRCA_status0_time2064_TCGA-BH-A0DV_rep2_on_MCF-7.merge.csv

