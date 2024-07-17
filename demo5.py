import numpy as np
import main_function as mf
import time
import warnings  
from scipy.integrate._quadpack_py import IntegrationWarning  
import argparse
import os
from parameter_list import Albedo, Temperature

import matplotlib.pyplot as plt

  
# 忽略 IntegrationWarning  
warnings.filterwarnings('ignore', category=IntegrationWarning)  


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--id', default="1" , type=str)
    parser.add_argument('--Albedo', default = Albedo , type=float)
    parser.add_argument('--Ntheta', default= 5 , type=int)
    parser.add_argument('--Nwave', default = 1 , type=int)
    args = parser.parse_args()


    os.makedirs(f'temp/R{args.id}/variables', exist_ok=True)
    os.makedirs(f'temp/R{args.id}/Results', exist_ok=True)

    t5 = time.time()
    wave_bound = np.array([5000 , 5000])*1e-9
    mf.thermal_spectrum(wave_bound, Temperature, Albedo=args.Albedo, id=args.id, Ntheta = args.Ntheta, NWavelength= args.Nwave)

    t6 = time.time()
    print("Total Time = ", t6 - t5, "s, Processing ALL DONE!")





