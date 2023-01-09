####

for i in {1..6}
do

sed 's/cluster1/cluster'${i}'/g' downsampling_fibroblast_PEER_factors.R > downsampling_fibroblast_cluster${i}_PEER_factors.R


sed 's/cluster1/cluster'${i}'/g' downsampling_fibroblast_PEER_factors.sh > downsampling_fibroblast_cluster${i}_PEER_factors.sh

qsub downsampling_fibroblast_cluster${i}_PEER_factors.sh

done




####
