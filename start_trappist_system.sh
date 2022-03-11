#!/bin/bash
#
#SBATCH --job-name=trappist_system
#SBATCH --output=trappist_system.txt
#
#SBATCH --partition=public-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --dependency=singleton
#SBATCH --time=96:00:00

data_dir=~/scratch/trappist
RESULTS=$data_dir/trappist_system.bin
if [ ! -f $RESULTS ]
then
    posidonius start -s trappist_system.json $data_dir/trappist_system.bin $data_dir/trappist_system_history.bin
else
  RESULT=$(grep -i -E "Simulation completed|PANIC" trappist_system.txt)
  if [[ $? -eq 0 ]]; then
    echo "Cancel rest of the jobs!"
    scancel --jobname $SLURM_JOB_NAME
  else
    echo "Do NOT cancel anything, the simulation has not finished"
    posidonius resume -s $data_dir/trappist_system.bin $data_dir/trappist_system_history.bin
  fi
fi
