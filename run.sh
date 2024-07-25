#!/bin/bash
#SBATCH --job-name=PythonTest
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=1
#SBATCH --partition=wzhcnormal
#SBATCH -o log/%j.loop
#SBATCH -e log/%j.loop

. /public/software/apps/anaconda3/5.2.0/etc/profile.d/conda.sh

echo "Maxwell"

module load apps/anaconda3/5.2.0
conda activate test

# python  Full_spectrum_main.py --id 3 --Ntheta 15 --Nwave 200
python plot_lib.py
# python demo_vertify.py
# python Tmap_2D_plot.py


