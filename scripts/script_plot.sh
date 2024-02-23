# !/bin/bash


# pwd=/Volumes/Expansion/Mac_data/posidonius/from_baobab
pwd=/Volumes/Expansion/baobab/posidonius

# # Plot Millholland
# json_dir=${pwd}/json_millholland
# data_dir=${pwd}/millholland
# cases=("test_ef_disk1g_outerplanet") # "test_ef_disk01g_outerplanet")


# # Plot Sub-stellar comparison
# json_dir=${pwd}/json_files/json_substellar
# data_dir=${pwd}/substellar
# cases=( "t1_earth")
# # cases=( "t1_earth" "t1_intermediate" "t1_mercury" "t1_small")

# Plot Sub-stellar comparison
json_dir=${pwd}/json_files/json_substellar
data_dir=${pwd}/substellar/andrade
cases=(
"t1_earth_andrade_0001_h01" 
"t1_earth_andrade_0001"
"t1_earth_andrade_001"
"t1_earth_andrade_01")
# data_dir=${pwd}/substellar/sunbergcooper
# cases=( 
# "t1_earth_sunbergcooper_0001_h01"
# "t1_earth_sunbergcooper_0001"
# "t1_earth_sunbergcooper_001"
# "t1_earth_sunbergcooper_01")

# # Plot Sub-stellar comparison
# json_dir=${pwd}/json_files/json_mmr_test
# data_dir=${pwd}/mmr_test
# cases=("t1_AS_notides")
# # cases=("t1_AS-t1a-mass-l-sma_res_allearth" "t1_AS-t1a-mass-l_res_allearth" "t1_AS-t1a-mass_res_allearth" "t1_AS_allearth" "t1_AS_allearth_O20" "t1_AS_gr" "t1_AS_gr_allearth" "t1_AS_gr_allearth_O20" "t1_AS_notides")

# # Plot Sub-stellar tests WHFast
# json_dir=~/scratch/posidonius/json_files/json_bolmont_ttv/AS-ci/start_0_1yrs/t1_ttv_AS_allearth
# data_dir=~/scratch/posidonius/bolmont_ttv/AS-ci/start_0_1yrs/t1_ttv_AS_allearth/
# cases=("test_tstp0.001_his_0.01")


# Get the length of the list
cases_length=${#cases[@]}
# Initialize a counter
counter=0
echo "Start ... "
# Loop over the list using a while loop
while [ $counter -lt $cases_length ]; do
    case="${cases[$counter]}"
    echo "\n \t\t Processing ${case} \n"
    # python scripts/all_explore_history.py ${json_dir}/${case}.json ${data_dir}/${case}/${case}_history.bin
    python scripts/all_explore_history.py ${json_dir}/${case}.json ${data_dir}/${case}_history.bin
    ((counter++))
done
