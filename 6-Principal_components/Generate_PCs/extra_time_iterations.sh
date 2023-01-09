##

# new_names=(B_IN B_MEM CD4_NC CD4_ET CD4_SOX4 CD8_NC CD8_ET CD8_S100B DC Mono_C Mono_NC NK_R NK Plasma)

# run_MatrixQTL.R

for i in {1..13}
do

cd ./option_${i}

grep mins *.o* > time_option_${i}.txt
grep Converged *.o* > iterations_${i}.txt

echo option_${i}

cd ../

done


##
