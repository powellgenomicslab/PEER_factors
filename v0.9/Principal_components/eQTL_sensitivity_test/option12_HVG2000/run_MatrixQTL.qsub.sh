## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -pe smp 1
#$ -r yes
#$ -l mem_requested=20G
#$ -N run_MatrixQTL_ct_name_PF_num
#$ -o ./ct_name/PF_num/stdout_run_MatrixQTL_ct_name_PF_num
#$ -e ./ct_name/PF_num/stderr_run_MatrixQTL_ct_name_PF_num
#$ -t 1-22

cd $SGE_O_WORKDIR

i=${SGE_TASK_ID};

# touch ${i}.txt

#module load briglo/R/3.6.0

Rscript=/share/ScratchGeneral/angxue/software/R-4.0.5/bin/Rscript

# mkdir -p ct_name

touch ./ct_name/PF_num/chr${i}_run_MatrixQTL_ct_name.log
> ./ct_name/PF_num/chr${i}_run_MatrixQTL_ct_name.log

$Rscript --vanilla run_MatrixQTL.R ${i} PF_num ct_name >> ./ct_name/PF_num/chr${i}_run_MatrixQTL_ct_name.log 

