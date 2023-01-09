####
for i in {1..13}
do
mkdir -p option_${i}
sed 's/nr/'${i}'/g' identify_PEER_factors_option.sh > ./option_${i}/identify_PEER_factors_option_${i}.sh

cd ./option_${i}
sh identify_PEER_factors_option_${i}.sh
cd ../

done


####
