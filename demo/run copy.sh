#!/bin/bash
#SBATCH --job-name=PythonTest
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --partition=wzhctdnormal
#SBATCH -o log/%j.loop
#SBATCH -e log/%j.loop

. /public/software/apps/anaconda3/5.2.0/etc/profile.d/conda.sh

echo "Maxwell" 

module load apps/anaconda3/5.2.0
conda activate test
export lavatype=zero

PARAMS=$(cat params.txt)
python lib/Full_spectrum_main.py $PARAMS

# python  lib/Full_spectrum_main.py --id 1 --Ntheta 1 --Nwave 1201 --LB 0.3 --UB 12.3 --mode TR --lavatype zero --Nsubpro 190 --heat_redist No
# python  lib/Full_spectrum_main.py --id 2 --Ntheta 1 --Nwave 1201 --LB 0.3 --UB 12.3 --mode TR --lavatype zero --Nsubpro 190 --heat_redist Full
# python  lib/Full_spectrum_main.py --id 3 --Ntheta 1 --Nwave 1201 --LB 0.3 --UB 12.3 --mode TR --lavatype low --Nsubpro 190 --heat_redist No
# python  lib/Full_spectrum_main.py --id 4 --Ntheta 1 --Nwave 1201 --LB 0.3 --UB 12.3 --mode TR --lavatype high --Nsubpro 190 --heat_redist No
# python  lib/Full_spectrum_main.py --id 5 --Ntheta 1 --Nwave 1201 --LB 0.3 --UB 12.3 --mode TR --lavatype high_OH --Nsubpro 190 --heat_redist No
# python  lib/Full_spectrum_main.py --id 6 --Ntheta 30 --Nwave 1201 --LB 0.3 --UB 12.3 --mode PC --lavatype high --Nsubpro 190 --heat_redist No
# python  lib/Full_spectrum_main.py --id 7 --Ntheta 30 --Nwave 1201 --LB 0.3 --UB 12.3 --mode PC --lavatype low --Nsubpro 190 --heat_redist No
# python  lib/Full_spectrum_main.py --id 8 --Ntheta 1 --Nwave 1201 --LB 0.3 --UB 12.3 --mode TR --lavatype zero --Nsubpro 62 --heat_redist Full
# python  lib/Full_spectrum_main.py --id 9 --Ntheta 30 --Nwave 200 --LB 2.3 --UB 5.2 --mode PC --lavatype zero --Nsubpro 190 --heat_redist No
# python  lib/Full_spectrum_main.py --id 5 --Ntheta 60 --Nwave 5 --LB 2.99 --UB 3.01 --mode PC --lavatype one --Nsubpro 124 --heat_redist No --roughness 0
# python  lib/Full_spectrum_main.py --id 8 --Ntheta 60 --Nwave 5 --LB 2.99 --UB 3.01 --mode PC --lavatype one --Nsubpro 250 --heat_redist No --roughness 0
# python  lib/Full_spectrum_main.py --id 9 --Ntheta 15 --Nwave 5 --LB 2.99 --UB 3.01 --mode PC --lavatype one --Nsubpro 125 --heat_redist No --roughness 1000
# python  lib/Full_spectrum_main.py --id 11 --Ntheta 1 --Nwave 500 --LB 0.3 --UB 10 --mode TR --lavatype high --Nsubpro 250 --heat_redist No 

# python  lib/Full_spectrum_main.py --id 11 --Ntheta 60 --Nwave 5 --LB 2.99 --UB 3.01 --mode PC --lavatype one --Nsubpro 125 --heat_redist No --roughness 0
# python  lib/Full_spectrum_main.py --id 12 --Ntheta 60 --Nwave 5 --LB 2.99 --UB 3.01 --mode PC --lavatype one --Nsubpro 125 --heat_redist No --roughness 1000
# python  lib/Full_spectrum_main.py --id 11 --Ntheta 60 --Nwave 5 --LB 2.99 --UB 3.01 --mode PC --lavatype one --Nsubpro 125 --heat_redist No --roughness 0
# python  lib/Full_spectrum_main.py --id 12 --Ntheta 60 --Nwave 5 --LB 2.99 --UB 3.01 --mode PC --lavatype one --Nsubpro 125 --heat_redist No --roughness 1000

# python  lib/Full_spectrum_main.py --id 13 --Ntheta 60 --Nwave 100 --LB 0.1 --UB 5 --mode PC --lavatype low --Nsubpro 125 --heat_redist No 
# python  lib/Full_spectrum_main.py --id 14 --Ntheta 60 --Nwave 100 --LB 0.1 --UB 5 --mode PC --lavatype high --Nsubpro 125 --heat_redist No 

# python lib/transit_cal.py --name R6

# python  lib/plot_lib.py
# python  lib/bond_albedo_calculator.py
# python demo_vertify.py
# python Tmap_2D_plot.py
# python lib/lava_data.py
# python telescope_measure/test.py

# python lib/plot_paper.py
echo "DONE"


