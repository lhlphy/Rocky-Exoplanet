#!/bin/bash
#SBATCH --job-name=PythonTest
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --partition=wzhctdnormal
#SBATCH -o log/%j.loop
#SBATCH -e log/%j.loop

. /public/software/apps/anaconda3/5.2.0/etc/profile.d/conda.sh

echo "Maxwell" 

module load apps/anaconda3/5.2.0
conda activate test

python  lib/Full_spectrum_main.py --id 10 --Ntheta 19 --Nwave 1201 --LB 0.3 --UB 12.3
# python  lib/plot_lib.py
# python  lib/bond_albedo_calculator.py
# python demo_vertify.py
# python Tmap_2D_plot.py
# python lib/lava_data.py
# python telescope_measure/test.py


