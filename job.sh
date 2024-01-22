#!/bin/bash 
#SBATCH -J test
#SBATCH -D /data/scratch2/mpstar
#SBATCH -o out.out 
#SBATCH --partition=nowick
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=30-00:00:00
#SBATCH --mail-type=end 
#SBATCH --mail-user=mpstar@zedat.fu-berlin.de
hostname
date
eval "$(/home/mpstar/anaconda3/bin/conda shell.bash hook)"
conda activate r_env
Rscript ./test.R
