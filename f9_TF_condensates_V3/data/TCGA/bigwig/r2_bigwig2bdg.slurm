#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 48:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_r2.log

# Run program


time 
for cancertype in BRCA COAD CESC LIHC PRAD 
do
  for ii in ${cancertype}/*bw
  do
    prename=$(basename $ii .insertions.bw)
    echo $prename $ii ${cancertype}/${prename}.bdg
    bigWigToBedGraph $ii ${cancertype}/${prename}.bdg
    wait
  done
done