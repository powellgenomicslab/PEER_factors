##

new_names=(B_IN B_MEM CD4_NC CD4_ET CD4_SOX4 CD8_NC CD8_ET CD8_S100B DC Mono_C Mono_NC NK_R NK Plasma)

# new_names2=(BimmNaive Bmem CD4all CD4effCM CD4TGFbStim CD8all CD8eff CD8unknown DC MonoC MonoNC NKact NKmat Plasma)

for k in {0..13}
do

for i in {0..50}
do

mkdir -p ./${new_names[k]}/${i}/
echo ${new_names[k]}
echo ${i} 

sed 's/PF_num/'${i}'/g;s/ct_name/'${new_names[k]}'/g' run_MatrixQTL.qsub.sh > run_MatrixQTL_${new_names[k]}_PF${i}.qsub.sh


qsub run_MatrixQTL_${new_names[k]}_PF${i}.qsub.sh

rm -f run_MatrixQTL_${new_names[k]}_PF${i}.qsub.sh

done

done

##
