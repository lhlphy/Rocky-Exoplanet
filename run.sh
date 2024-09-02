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

python  lib/Full_spectrum_main.py --id 1 --Ntheta 20 --Nwave 5
# python plot_lib.py
# python demo_vertify.py
# python Tmap_2D_plot.py


