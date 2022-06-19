## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -pe smp 1
#$ -r yes
#$ -l mem_requested=10G
#$ -N summarize_nr_signal_ct_name
#$ -o ./stdout_summarize_nr_signal_ct_name
#$ -e ./stderr_summarize_nr_signal_ct_name
#$ -t 1-14

cd $SGE_O_WORKDIR

i=${SGE_TASK_ID};



Rscript=/share/ScratchGeneral/angxue/software/R-4.0.5/bin/Rscript


$Rscript --vanilla summarize_number_signals.R ${i} 
