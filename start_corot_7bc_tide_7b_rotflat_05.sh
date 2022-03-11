#!/bin/bash
#
#SBATCH --job-name=corot_7bc_tide_7b_rotflat_05
#SBATCH --output=corot_7bc_tide_7b_rotflat_05.txt
#
#SBATCH --partition=public-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --dependency=singleton
#SBATCH --time=96:00:00

data_dir=~/scratch/corot_7bc_tides
RESULTS=$data_dir/corot_7bc_tide_7b_tides_05.bin
if [ ! -f $RESULTS ]
then
    posidonius start -s corot_7bc_tide_7b_tides_05.json $data_dir/corot_7bc_tide_7b_tides_05.bin $data_dir/corot_7bc_tide_7b_tides_05_history.bin
else
  RESULT=$(grep -i -E "Simulation completed|PANIC" corot_7bc_tide_7b_tides_05.txt)
  if [[ $? -eq 0 ]]; then
    echo "Cancel rest of the jobs!"
    scancel --jobname $SLURM_JOB_NAME
  else
    echo "Do NOT cancel anything, the simulation has not finished"
    posidonius resume -s $data_dir/corot_7bc_tide_7b_tides_05.bin $data_dir/corot_7bc_tide_7b_tides_05_history.bin 
  fi
fi
