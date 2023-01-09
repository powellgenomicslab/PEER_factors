## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -pe smp 1
#$ -r yes
#$ -l mem_requested=20G
#$ -N summarize_nr_signal
#$ -o ./stdout_summarize_nr_signal
#$ -e ./stderr_summarize_nr_signal
#$ -t 1-14

cd $SGE_O_WORKDIR

i=${SGE_TASK_ID};



Rscript=/share/ScratchGeneral/angxue/software/R-4.0.5/bin/Rscript


$Rscript --vanilla summarize_number_signals.R ${i} 
