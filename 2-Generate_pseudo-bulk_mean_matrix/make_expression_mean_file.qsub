## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -pe smp 1
#$ -r yes
#$ -l mem_requested=20G
#$ -N generate_mean_mx
#$ -t 1-14

hostname

cd $SGE_O_WORKDIR

i=${SGE_TASK_ID};

# touch ${i}.txt

#module load briglo/R/3.6.0

Rscript=/share/ScratchGeneral/angxue/software/R-4.0.5/bin/Rscript

touch make_expression_mean_file_${i}.log

> make_expression_mean_file_${i}.log

echo ${JOB_ID} >> make_expression_mean_file_${i}.log

$Rscript --vanilla make_expression_mean_file.R ${i} >> make_expression_mean_file_${i}.log  

