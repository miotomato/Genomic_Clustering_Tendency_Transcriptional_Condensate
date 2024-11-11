export PATH="/nv/vol190/zanglab/zw5j/env/anaconda2/bin:$PATH"
export PATH="/nv/vol190/zanglab/zw5j/env/hicpro/installation/HiC-Pro_2.10.0/bin/utils/:$PATH"
# module load R
module load bowtie2
module load samtools
module load gcc/7.1.0  openmpi/3.1.4
module load R/3.5.3


HiC-Pro -i /nv/vol190/zanglab/zw5j/projects_data/tf_condensates/GSE104333_hct116_hic -o hct116_rad21_auxin -c confit_test_latest_Mboi.txt -p 


