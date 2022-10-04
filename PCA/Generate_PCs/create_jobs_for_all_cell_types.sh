##

# new_names=(B_IN B_MEM CD4_NC CD4_ET CD4_SOX4 CD8_NC CD8_ET CD8_S100B DC Mono_C Mono_NC NK_R NK Plasma)

# run_MatrixQTL.R

for i in {1..13}
do
mkdir -p option_${i}
sed 's/nr/'${i}'/g' identify_PCs_option.qsub.sh > ./option_${i}/identify_PCs_option_${i}.qsub.sh

cd ./option_${i}
qsub identify_PCs_option_${i}.qsub.sh
cd ../

done


##
