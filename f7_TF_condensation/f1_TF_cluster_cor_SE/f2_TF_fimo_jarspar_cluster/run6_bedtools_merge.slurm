#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 48:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_run6.log


  
file_dir=f5_motif_within_300bp/
outdir=f6_bedtools_merge
mkdir $outdir

for ii in $file_dir/*.bed
do
  prename=$(basename $ii .bed)
  echo $ii $prename
  bedtools merge -d 300 -i $ii > $outdir/${prename}.bed
done
