pwd=/Volumes/Expansion/ygg/posidonius
json_dir=${pwd}/json_files/json_highrotation/thesis
data_dir=${pwd}/highrotation/thesis

# case=t1-ctl10pc
case=t1-26_andrade_inter

posidonius resume --historic-snapshot-period 3652.5 --recovery-snapshot-period 36525.0 $data_dir/${case}.bin $data_dir/${case}_history.bin