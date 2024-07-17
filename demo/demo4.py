import numpy as np
import main_function as mf
import time
import warnings  
from scipy.integrate._quadpack_py import IntegrationWarning  
import argparse
import os

import matplotlib.pyplot as plt

  
# 忽略 IntegrationWarning  
warnings.filterwarnings('ignore', category=IntegrationWarning)  


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--id', default="0" , type=str)
    parser.add_argument('--Albedo', default=0.1 , type=float)
    args = parser.parse_args()


    os.makedirs(f'temp/R{args.id}/variables', exist_ok=True)
    os.makedirs(f'temp/R{args.id}/Results', exist_ok=True)

    t5 = time.time()
    # 创建一个空字典来存储变量  
    variables = {}  
    
    cal_size = [18]
    Theta_list = np.linspace(0, 2*np.pi, cal_size[0]+1)
    Theta_list = Theta_list[:-1]
    #Theta_list = [np.pi/4, np.pi/3]

    for j, Theta in enumerate(Theta_list):
        t1 = time.time()
        star_flux = mf.Cal_star_flux(Theta)
        
        var_name = "Tmap_"  + str(int(Theta*180/np.pi))
        variables[var_name] = mf.Tmap(Theta, args.id, star_flux, Albedo=args.Albedo)

        t2 = time.time()
        print("Theta = ", Theta, "Time = ", t2 - t1, "s, Processing Done!")

    t6 = time.time()
    print("Total Time = ", t6 - t5, "s, Processing ALL DONE!")

    # save "variables" to temp/ folder
    #np.save(f'temp/{args.id}/variables/variables.npy', variables)
    np.save(f'temp/R{args.id}/variables/Theta_list.npy', Theta_list)
    np.save(f'temp/R{args.id}/variables/star_flux.npy', star_flux)

    # plot the results
    # result_ploter(TOT_Intensity, "TOT_Intensity", Theta_list, Coarse, args.id)
    # result_ploter(TOT_Diffuse, "TOT_Diffuse", Theta_list, Coarse, args.id)
    # result_ploter(TOT_Specular, "TOT_Specular", Theta_list, Coarse, args.id)




