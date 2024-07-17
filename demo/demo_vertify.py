import numpy as np
from main_function import *
from parameter_list import Albedo, Temperature

vertify_radiation(np.array([5000,5000])* 1e-9, Temperature = Temperature, Albedo=0 , id=9, Ntheta = 10, NWavelength = 1)