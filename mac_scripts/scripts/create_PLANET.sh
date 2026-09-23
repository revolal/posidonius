# !/bin/bash

# dir=/media/arevol/Expansion/Dynaclim_data/thesis
# #/Volumes/Expansion/ygg/posidonius
# # json_dir=${dir}/json_files/json_substellar/thesis
# # data_dir=${dir}/substellar/thesis
# json_dir=${dir}
# data_dir=${dir}
# cases=(
#     "t1_andrade_biger"
#     "t1_andrade_inter"
#     "t1_andrade_earth"
#     "t1_notides"
# )
# # "t1_kaula_ctl_earth10pc"
# # "t1_kaula_ctl_earth"
# # "t1_ctl_earth10pc"
# # "t1_ctl_earth"
# # "t1_maxwell_biger"
# # "t1_maxwell_inter"
# # "t1_maxwell_earth"
# # "t1_sundbergcooper_biger"
# # "t1_sundbergcooper_inter"
# # "t1_sundbergcooper_earth"
# # "t1_andrade_biger"
# # "t1_andrade_inter"
# # "t1_andrade_earth"
# # "t1_notides"
# # )

dir=/home/arevol/Desktop/Dynaclim/posidonius
#/Volumes/Expansion/ygg/posidonius
# json_dir=${dir}/json_files/json_substellar/thesis
# data_dir=${dir}/substellar/thesis
json_dir=${dir}
data_dir=${dir}
cases=(
    "k2_earthtides"
)

# # Get the length of the list
cases_length=${#cases[@]}

# Initialize a counter
counter=0
echo "Start ... create PLANET "
# Loop over the list using a while loop
while [ $counter -lt $cases_length ]; do
    case="${cases[$counter]}"
    echo "> Processing ${case} "
    python create_PLANETi.py ${json_dir}/${case}.json ${data_dir}/${case}_history.bin
    ((counter++))
done


counter=0
echo "Start ... create Timing"
# Loop over the list using a while loop
while [ $counter -lt $cases_length ]; do
    case="${cases[$counter]}"
    echo "> Processing ${case} "
    python create_Timings.py ${data_dir}/${case}_history_
    ((counter++))
done