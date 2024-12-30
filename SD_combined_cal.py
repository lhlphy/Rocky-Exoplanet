import numpy as np
import argparse
import warnings
from scipy.integrate._quadpack_py import IntegrationWarning 
import os
import sys

# Add the lib directory to the system path
sys.path.append(os.path.join(os.path.dirname(__file__), 'lib'))
from function_library import Tmap

# # 忽略 IntegrationWarning  
# warnings.filterwarnings('ignore', category=IntegrationWarning) 

class SD_combined_cal:
    def __init__(self):
        warnings.filterwarnings('ignore', category=IntegrationWarning) 
        warnings.filterwarnings("ignore")
        parser = argparse.ArgumentParser()
        parser.add_argument('--target', default='none' , type=str)
        parser.add_argument('--id', default=0 , type=int)
        # parser.add_argument('--coarse', default=0 , type=float)
        # parser.add_argument('--specular', default=0.5 , type=float)
        # parser.add_argument('--diffuse', default=0.5 , type=float)
        # parser.add_argument('--Albedo', default=0.5 , type=float)
        parser.add_argument('--Ntheta', default=5 , type=int)
        parser.add_argument('--Nwave', default=1 , type=int)
        parser.add_argument('--LB', default=1, type=float)  # lower bound (unit: um)
        parser.add_argument('--UB', default=1, type=float)  # upper bound (unit: um)
        parser.add_argument('--mode', default='PC', type=str) # Calculation target "PC":Phase curve ; "TR":Transit
        parser.add_argument('--lavatype', default= 'low', type = str) # albedo model type of lava "low","high","mode1","zero"
        parser.add_argument('--Nsubpro', default = 1, type = int)  # 求解时分为多少个子进程
        parser.add_argument('--heat_redist', default = 'No', type = str) # 热再分配 "No","Full","Yes"
        parser.add_argument('--roughness', default = 0, type = float)  # roughness
        parser.add_argument('--FRnormal', default=0.1, type=float) # Fresnel reflection coefficient when normal incidence
        parser.add_argument('--polarization', default = "None", type = str)  # polarization "None", "P", "S"
        args = parser.parse_args()
        Wavelength_bound = np.array([args.LB, args.UB]) *1e-6   # Wavelength range (m)

        # 设置环境变量
        os.environ['mode'] = args.mode
        os.environ['lavatype'] = args.lavatype
        os.environ['heat_redist'] = args.heat_redist
        os.environ['roughness'] = str(args.roughness)
        os.environ['polarization'] = args.polarization
        os.environ['FRnormal'] = str(args.FRnormal)
        os.makedirs('log', exist_ok= True)
        os.makedirs('temp', exist_ok=True)
        os.makedirs(f'temp/R{args.id}/Results', exist_ok=True)
        self.id = args.id
        
        self.TMAP0 = Tmap(0, args.id)
        
    def data_load(self, Fresnel_folder, Lambert_folder):
        
        pass

    



if __name__ == '__main__':

    TMAP0 = Tmap(0, id)
    