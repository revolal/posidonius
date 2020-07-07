#!/bin/sh
python setup.py install --user
cargo install --path . --force
python cases/trappist1_multi.py target/trappist1_multi.json
posidonius start target/trappist1_multi.json target/trappist1_multi.bin target/trappist1_multi_history.bin 


# #!/bin/bash
# #

# #SBATCH --job-name=MULTI_T1_POS
# #SBATCH --output=output.txt
# #
# #SBATCH --partition=s2
# #SBATCH --ntasks=1
# path=$PWD
# NAME=${path##*/}
# data_dir=/hpcstorage/revolal/Posidonius/Multi-Trappist/$NAME
# posidonius start --silent trappist1_multi.json $data_dir/trappist1_multi.bin $data_dir/trappist1_history_multi.bin
            