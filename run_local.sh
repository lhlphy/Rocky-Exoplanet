#!/bin/bash
# test script for wsl, test project on local machine
# this script is for wsl, not for windows
echo "Maxwell" 
export lavatype=zero

python3 lib/Full_spectrum_main.py --id 12 --Ntheta 6 --Nwave 5 --LB 2.99 --UB 3.01 --mode PC --lavatype one --Nsubpro 5 --heat_redist No --roughness 1000

# python lib/transit_cal.py

# python  lib/plot_lib.py
# python  lib/bond_albedo_calculator.py
# python demo_vertify.py
# python Tmap_2D_plot.py
# python lib/lava_data.py
# python telescope_measure/test.py

# python lib/plot_paper.py
echo "DONE"