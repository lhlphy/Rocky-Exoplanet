import numpy as np
import argparse
import warnings
from scipy.integrate._quadpack_py import IntegrationWarning 
import os

# # 忽略 IntegrationWarning  
# warnings.filterwarnings('ignore', category=IntegrationWarning) 

if __name__ == '__main__':
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
    args = parser.parse_args()
    Wavelength_bound = np.array([args.LB, args.UB]) *1e-6   # Wavelength range (m)
    parser.add_argument('--Model', default = "None", type = str ) 

    # 设置环境变量
    os.environ['mode'] = args.mode
    os.environ['lavatype'] = args.lavatype
    os.environ['heat_redist'] = args.heat_redist
    os.environ['roughness'] = str(args.roughness)
    os.environ['Model'] = args.Model
    os.makedirs('log', exist_ok= True)
    os.makedirs('temp', exist_ok=True)
    
    # # 存储变量，传递
    # os.makedirs('log', exist_ok= True)
    # with open('log/temp_vars.txt', 'w') as f:
    #     f.write(f"{args.mode}\n")
    #     f.write(f"{args.lavatype}\n")
    #     f.write(f"{args.heat_redist}\n")

    from Full_spectrum import Full_spectrum   # have to import after log/temp_vars.txt is created
    Full_spectrum(Wavelength_bound, args)
    from transit_cal import Transit_cal  # do the flux correction immediately, generate R{id}copy folder
    Transit_cal(f'R{args.id}')
    

        
    