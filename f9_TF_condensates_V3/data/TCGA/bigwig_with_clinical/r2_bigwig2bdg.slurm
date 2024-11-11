#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 48:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o r2.log

# Run program


time for ii in *.bw
do
   prename=$(basename $ii .insertions.bw)
   echo $prename
   bigWigToBedGraph $ii ${prename}.bdg
   wait
done