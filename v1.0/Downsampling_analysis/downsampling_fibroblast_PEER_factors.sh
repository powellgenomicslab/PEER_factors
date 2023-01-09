## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -pe smp 1
#$ -r yes
#$ -l mem_requested=10G
#$ -N identify_PEER_factors_fibroblast_cluster1
#$ -o stdout_PEER_fibroblast_cluster1
#$ -e stderr_PEER_fibroblast_cluster1
#$ -t 1-30

hostname

cd $SGE_O_WORKDIR

i=${SGE_TASK_ID};


Rscript=/share/ScratchGeneral/angxue/software/R-4.0.5/bin/Rscript

$Rscript --vanilla downsampling_fibroblast_cluster1_PEER_factors.R ${i}

