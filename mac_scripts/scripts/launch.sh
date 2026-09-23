pwd=/Volumes/Expansion/ygg/posidonius
json_dir=${pwd}/json_files/json_highrotation/thesis
data_dir=${pwd}/highrotation/thesis

case=t1-ctl10pc
# case=t1-26_andrade_inter

posidonius start ${json_dir}/${case}.json ${data_dir}/${case}.bin ${data_dir}/${case}_history.bin