#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 48:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_TFMS_jarspar_by_pvalue.log

module load gcc/7.1.0
module load bedtools/2.26.0
	 
# file_dir=/nv/vol190/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/
file_dir=../../data/motif_fimo_jaspar_by_pvalue/
outdir=data_TFMS_jarspar_by_pvalue
mkdir $outdir

time for ii in $file_dir/*.bed
do
  prename=$(basename $ii .bed)
  echo $prename
  bedtools closest -a $ii -b $ii -D ref -fd -io -t first > $outdir/${prename}.tsv
done

# bedtools closest -a $notch_motif -b $notch_motif -d -io -t first > $outdir/notch_motif_closest.tsv

