## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -pe smp 1
#$ -r yes
#$ -l mem_requested=10G
#$ -N merge_results
#$ -o ./stdout_merge_results
#$ -e ./stderr_merge_results
#$ -t 1-14

cd $SGE_O_WORKDIR

k=${SGE_TASK_ID};

cd /share/ScratchGeneral/angxue/proj/vQTL/MatrixQTL/PF_test/raw/

new_names=(NA B_IN B_MEM CD4_NC CD4_ET CD4_SOX4 CD8_NC CD8_ET CD8_S100B DC Mono_C Mono_NC NK_R NK Plasma)

echo ${new_names[k]}

for i in {0..50}
do

echo ${i}

touch ./${new_names[k]}/${i}/${new_names[k]}_eQTL_all.txt
> ./${new_names[k]}/${i}/${new_names[k]}_eQTL_all.txt

head -n1 ./${new_names[k]}/${i}/${new_names[k]}_eQTL_chr1_FDR_005.txt >> ./${new_names[k]}/${i}/${new_names[k]}_eQTL_all.txt

for j in {1..22}
do
tail -n +2 ./${new_names[k]}/${i}/${new_names[k]}_eQTL_chr${j}_FDR_005.txt >> ./${new_names[k]}/${i}/${new_names[k]}_eQTL_all.txt

done

# wc -l ./${i}/CD4_NC_eQTL_all.txt

done




