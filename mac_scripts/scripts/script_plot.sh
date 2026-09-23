# !/bin/bash

# pwd=/Volumes/Expansion/Mac_data/posidonius/from_baobab
# pwd=/Volumes/Expansion/baobab/posidonius

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
# json_dir=${pwd}/json_files/json_substellar
# data_dir=${pwd}/substellar/andrade
# cases=(
# "t1_big_sulfr_andrade_ias15_h1"
# "t1_inter_sulfr_andrade_ias15_h1"
# "t1_earth_sulfr_andrade_ias15_h1"
# "t1_earth_andrade_ias15_h001"
# "t1_earth_andrade_ias15_h01"
# "t1_earth_andrade_ias15_h1"
# "t1_earth_andrade_ias15_h1_b"
# "t1_earth_andrade_ias15_h1_bc"
# "t1_earth_andrade_ias15_h1_bcd"
# "t1_earth_andrade_ias15_h1_bcde"
# "t1_earth_andrade_ias15_h1_bcdef"
# "t1_earth_andrade_ias15_h1_bcdefg"
#
# "t1_earth_andrade_whfast_h01_b"
# "t1_earth_andrade_whfast_h01_bc"
# "t1_earth_andrade_whfast_h01_bcd"
# "t1_earth_andrade_whfast_h01_bcde"
# "t1_earth_andrade_whfast_h01_bcdef"
# "t1_earth_andrade_whfast_h01_bcdefg"
# "t1_earth_andrade_whfast_h01"
# "t1_big_sulfr_andrade_whfast_h01"
# "t1_inter_sulfr_andrade_whfast_h01"
# "t1_earth"
# "t1_earth"
# "t1_earth_andrade_whfast_h01_hist007"
# )
# data_dir=${pwd}/substellar
# cases=(
# "t1_earth"
# )


# pwd=/Volumes/Expansion/baobab/posidonius

# json_dir=${pwd}/json_files/json_substellar
# data_dir=${pwd}/substellar/andrade/add_modes/tstep/AS
# cases=(
#     "t1_earth_addmodes"
#     # "t1_earth_nomodes"
#     # "t1_earth_addmodes_tstep"
#     # "t1_earth_nomodes_tstep"
# )
# # Get the length of the list
# cases_length=${#cases[@]}
# # Initialize a counter
# counter=0
# echo "Start ... "
# # Loop over the list using a while loop
# while [ $counter -lt $cases_length ]; do
#     case="${cases[$counter]}"
#     echo "\n \t\t Processing ${case} \n"
#     # python scripts/all_explore_history.py ${json_dir}/${case}.json ${data_dir}/${case}/${case}_history.bin
#     python scripts/all_explore_history.py ${json_dir}/${case}.json ${data_dir}/${case}_history.bin
#     ((counter++))
# done



# json_dir=${pwd}/json_files/json_substellar
# data_dir=${pwd}/substellar/andrade/add_modes/tstep/AS-k2corr
# cases=(
#     "t1_earthcorr_addmodes"
#     # "t1_noMMR_earthcorr_addmodes"
#     # "t1_breakMMR_earthcorr_addmodes"
#     # "t1_earthcorr_nomodes"
#     "t1_earthcorr_addmodes_tstep"
#     # "t1_earthcorr_nomodes_tstep"
# )

# # Get the length of the list
# cases_length=${#cases[@]}
# # Initialize a counter
# counter=0
# echo "Start ... "
# # Loop over the list using a while loop
# while [ $counter -lt $cases_length ]; do
#     case="${cases[$counter]}"
#     echo "\n \t\t Processing ${case} \n"
#     # python scripts/all_explore_history.py ${json_dir}/${case}.json ${data_dir}/${case}/${case}_history.bin
#     python scripts/all_explore_history.py ${json_dir}/${case}.json ${data_dir}/${case}_history.bin
#     ((counter++))
# done







# json_dir=${pwd}/json_files/json_substellar
# data_dir=${pwd}/substellar/andrade/add_modes/toygg
# cases=(
#     # "t1_earth_addmodes"
#     "t1_earthcorr_addmodes"
#     # "t1_earthx10_addmodes"
#     # "t1_earthcorrx10_addmodes"
# )
# # Get the length of the list
# cases_length=${#cases[@]}
# # Initialize a counter
# counter=0
# echo "Start ... "
# # Loop over the list using a while loop
# while [ $counter -lt $cases_length ]; do
#     case="${cases[$counter]}"
#     echo "\n \t\t Processing ${case} \n"
#     # python scripts/all_explore_history.py ${json_dir}/${case}.json ${data_dir}/${case}/${case}_history.bin
#     python scripts/all_explore_history.py ${json_dir}/${case}.json ${data_dir}/${case}_history.bin
#     ((counter++))
# done






pwd=/Volumes/Expansion/ygg/posidonius
json_dir=${pwd}/json_files/json_substellar/thesis
data_dir=${pwd}/substellar/thesis
cases=(
# "t1_kaula_ctl_earth10pc"
# "t1_kaula_ctl_earth"
"t1_ctl_earth10pc"
"t1_ctl_earth"
"t1_maxwell_biger"
"t1_maxwell_inter"
"t1_maxwell_earth"
"t1_sundbergcooper_biger"
"t1_sundbergcooper_inter"
"t1_sundbergcooper_earth"
"t1_andrade_biger"
"t1_andrade_inter"
"t1_andrade_earth"
"t1_notides"
)
pwd=/Volumes/Expansion/ygg/posidonius
json_dir=${pwd}/json_files/json_substellar/thesis
data_dir=${pwd}/substellar/thesis/1p
cases=(
"t1_andrade_earth_b"
"t1_andrade_earth_bc"
"t1_andrade_earth_bcd"
"t1_andrade_earth_bcde"
"t1_andrade_earth_bcdef"
"t1_andrade_earth_bcdefg"

)
# pwd=/Volumes/Expansion/ygg/posidonius
# json_dir=${pwd}/json_files/json_highrotation/thesis
# data_dir=${pwd}/highrotation/thesis
# cases=(
# "t1-26_ctl_earth10pc"
# "t1-26_ctl_earth"
# "t1-26_maxwell_biger"
# "t1-26_maxwell_inter"
# "t1-26_maxwell_earth"
# "t1-26_sundbergcooper_biger"
# "t1-26_sundbergcooper_inter"
# "t1-26_sundbergcooper_earth"
# "t1-26_andrade_biger"
# "t1-26_andrade_inter"
# "t1-26_andrade_earth"
# )

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
