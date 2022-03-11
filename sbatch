#!/bin/bash
#
#SBATCH --job-name=corot_7bc_tide_7b_rotflat_05
#SBATCH --output=corot_7bc_tide_7b_rotflat_05.txt
#
#SBATCH --partition=private-astro-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --dependency=singleton
#SBATCH --time=96:00:00

data_dir=/home/revolal/scratch/corot_7bc
RESULTS=$data_dir/corot_7bc_tide_7b_rotflat_05.bin
if [ ! -f $RESULTS ]
then
    posidonius start -s corot_7bc_tide_7b_rotflat_05.json 
else
  RESULT=$(grep -i -E "Simulation completed|PANIC" corot_7bc_tide_7b_rotflat_05.txt)
  if [[ $? -eq 0 ]]; then
    echo "Cancel rest of the jobs!"
    scancel --jobname $SLURM_JOB_NAME
  else
    echo "Do NOT cancel anything, the simulation has not finished"
    posidonius resume -s 
  fi
fi

