#!/bin/bash
#SBATCH --job-name=PythonTest
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --partition=wzhctdnormal
#SBATCH -o log/%j.loop
#SBATCH -e log/%j.loop

. /public/software/apps/anaconda3/5.2.0/etc/profile.d/conda.sh

echo "Maxwell"

module load apps/anaconda3/5.2.0
conda activate test

python  Full_spectrum_main.py --id 1 --Albedo 1 --Ntheta 5 --Nwave 200

