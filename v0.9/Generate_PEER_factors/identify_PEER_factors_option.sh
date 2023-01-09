## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -r yes
#$ -l mem_requested=10G,tmp_requested=10G,tmpfree=10G
#$ -N PEER_factors_option_nr
#$ -t 1-14


hostname

cd $SGE_O_WORKDIR

i=${SGE_TASK_ID};

Rscript=/share/ScratchGeneral/angxue/software/R-4.0.5/bin/Rscript

dir="/share/ScratchGeneral/angxue/proj/vQTL/PEER_factors/scale_center_test/no_outliers/"

$Rscript --vanilla ${dir}identify_PEER_factors.R ${i} nr 




