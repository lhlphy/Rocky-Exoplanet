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
    parser.add_argument('--coarse', default=0 , type=float)
    parser.add_argument('--specular', default=0.5 , type=float)
    parser.add_argument('--diffuse', default=0.5 , type=float)
    args = parser.parse_args()


    os.makedirs(f'temp/{args.id}/variables', exist_ok=True)
    os.makedirs(f'temp/{args.id}/Results', exist_ok=True)

    t5 = time.time()
    # 创建一个空字典来存储变量  
    variables = {}  
    
    cal_size = [18]
    TOT_Intensity = np.zeros(cal_size)
    TOT_Diffuse = np.zeros(cal_size)
    TOT_Specular = np.zeros(cal_size)
    # # 定义一个字符串作为变量名  
    # variable_name = "my_variable"  
    # # 在字典中为变量赋值  
    # variables[variable_name] = 42  
    # # 现在，可以通过字典访问该变量  
    # print(variables[variable_name])  # 输出: 42
    # Coarse_list = np.linspace(0 , np.pi/3 , cal_size[0])
    Theta_list = np.linspace(0, 2*np.pi, cal_size[0]+1)
    Theta_list = Theta_list[:-1]
    #Theta_list = [np.pi/4, np.pi/3]
    Coarse = args.coarse * np.pi/180

    for j, Theta in enumerate(Theta_list):
        t1 = time.time()

        var_name = "Intensity_" + str(int(Coarse*180/np.pi)) + "_" + str(int(Theta*180/np.pi))
        var_name2 = "Diffuse_ratio_" + str(int(Coarse*180/np.pi)) + "_" + str(int(Theta*180/np.pi))

        variables[var_name], variables[var_name2] = mf.global_intensity(Theta, Coarse = Coarse, SPE_REF=args.specular,DIF_REF=args.diffuse, id = args.id, Model = 'Gaussian_wave')
        TOT_Intensity[j] = variables[var_name].sum()
        TOT_Diffuse[j] = (variables[var_name2] * variables[var_name]).sum()
        TOT_Specular[j] = TOT_Intensity[j] - TOT_Diffuse[j]

        t2 = time.time()
        print("Coarse = ", Coarse, "Theta = ", Theta, "Time = ", t2 - t1, "s, Processing Done!")

    t6 = time.time()
    print("Total Time = ", t6 - t5, "s, Processing ALL DONE!")

    star_flux = mf.Cal_star_flux(Theta_list)

    # save "variables" to temp/ folder
    np.save(f'temp/{args.id}/variables/variables.npy', variables)
    np.save(f'temp/{args.id}/variables/TOT_Intensity.npy', TOT_Intensity)
    np.save(f'temp/{args.id}/variables/TOT_Diffuse.npy', TOT_Diffuse)
    np.save(f'temp/{args.id}/variables/TOT_Specular.npy', TOT_Specular)
    np.save(f'temp/{args.id}/variables/Theta_list.npy', Theta_list)
    np.save(f'temp/{args.id}/variables/star_flux.npy', star_flux)

    # plot the results
    # result_ploter(TOT_Intensity, "TOT_Intensity", Theta_list, Coarse, args.id)
    # result_ploter(TOT_Diffuse, "TOT_Diffuse", Theta_list, Coarse, args.id)
    # result_ploter(TOT_Specular, "TOT_Specular", Theta_list, Coarse, args.id)

    vars = [TOT_Intensity, TOT_Diffuse, TOT_Specular]/star_flux
    name = ["Total reflect", "Diffuse", "Specular"]
    mf.multi_result_plotter(vars, name, Theta_list, Coarse, args.id)



