import numpy as np
from main_function import *
from parameter_list import  Temperature


# 理论模型与我的模型的相互验证
vertify_radiation(np.array([4000,5000])* 1e-9, Temperature = Temperature , id=1, Ntheta = 10, NWavelength = 2)