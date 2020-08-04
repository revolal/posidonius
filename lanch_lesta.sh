#!/bin/bash
#

#SBATCH --job-name=TTV_T1_POS
#SBATCH --output=output.txt
#
#SBATCH --partition=p5
#SBATCH --ntasks=1
path=$PWD
NAME=${path##*/}
data_dir=/hpcstorage/revolal/Posidonius/T1e
posidonius start --silent $data_dir/trappist_1e.json $data_dir/trappist_1e.bin $data_dir/trappist_1e_history.bin

