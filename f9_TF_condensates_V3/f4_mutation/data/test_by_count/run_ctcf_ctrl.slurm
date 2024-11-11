#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o ctcf_constitutive_mutatioinCount.out

time bedtools intersect \
-a /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_constitutive_0.8.bed \
-b pancancer_mutation.bed \
-c > ctcf_constitutive_mutatioinCount.bed
