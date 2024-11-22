#!/bin/bash
#SBATCH --job-name=PythonTest
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=wzhcnormal
#SBATCH -o log/%j.loop
#SBATCH -e log/%j.loop

. /public/software/apps/anaconda3/5.2.0/etc/profile.d/conda.sh

echo "Maxwell" 

module load apps/anaconda3/5.2.0
conda activate test

# env
id=5

# python  lib/Full_spectrum_main.py --id $id --Ntheta 1 --Nwave 1201 --LB 0.3 --UB 12.3 --mode TR --lavatype zero --Nsubpro 190 --heat_redist No
# python  lib/Full_spectrum_main.py --id $id --Ntheta 1 --Nwave 1201 --LB 0.3 --UB 12.3 --mode TR --lavatype zero --Nsubpro 190 --heat_redist Full
# python  lib/Full_spectrum_main.py --id $id --Ntheta 1 --Nwave 1201 --LB 0.3 --UB 12.3 --mode TR --lavatype low --Nsubpro 190 --heat_redist No
# python  lib/Full_spectrum_main.py --id $id --Ntheta 1 --Nwave 1201 --LB 0.3 --UB 12.3 --mode TR --lavatype high --Nsubpro 190 --heat_redist No
# python  lib/Full_spectrum_main.py --id $id --Ntheta 1 --Nwave 1201 --LB 0.3 --UB 12.3 --mode TR --lavatype high_OH --Nsubpro 190 --heat_redist No
# python  lib/Full_spectrum_main.py --id $id --Ntheta 30 --Nwave 1201 --LB 0.3 --UB 12.3 --mode PC --lavatype high --Nsubpro 190 --heat_redist No
# python  lib/Full_spectrum_main.py --id $id --Ntheta 30 --Nwave 1201 --LB 0.3 --UB 12.3 --mode PC --lavatype low --Nsubpro 190 --heat_redist No
# python  lib/Full_spectrum_main.py --id $id --Ntheta 1 --Nwave 1201 --LB 0.3 --UB 12.3 --mode TR --lavatype zero --Nsubpro 62 --heat_redist Full

# 获取当前脚本的路径
SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_NAME="run.sh"

echo "<$SCRIPT_NAME>"

# 目标目录
DEST_DIR="$SCRIPT_PATH/temp/R$id"

# 确保目标目录存在
mkdir -p "$DEST_DIR"

# 复制当前脚本到目标目录
cat "$SCRIPT_PATH/$SCRIPT_NAME" > "$DEST_DIR/$SCRIPT_NAME"

echo "Script copied to $DEST_DIR/$SCRIPT_NAME"


# python  lib/plot_lib.py
# python  lib/bond_albedo_calculator.py
# python demo_vertify.py
# python Tmap_2D_plot.py
# python lib/lava_data.py
# python telescope_measure/test.py
echo "DONE"


