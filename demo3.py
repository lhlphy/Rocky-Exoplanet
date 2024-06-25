import numpy as np
import main_function as mf
import time
import warnings  
from scipy.integrate._quadpack_py import IntegrationWarning  
import matplotlib.pyplot as plt
import argparse

  
# 忽略 IntegrationWarning  
warnings.filterwarnings('ignore', category=IntegrationWarning)  
  
# 如果你想在之后的代码中恢复警告，可以使用以下代码：  
# warnings.filterwarnings('default', category=IntegrationWarning)


    # 按列绘制TOT_Intensity,x轴为Theta_list,绘制在同一副图中
def result_ploter(Var, name, cal_size, Theta_list, Coarse_list):
    plt.figure()
    for i in range(cal_size[0]):
        plt.plot(Theta_list, Var[i,:], label = "Coarse = " + str(int(Coarse_list[i]*180/np.pi)))
    plt.xlabel('Theta')
    plt.ylabel(name)
    plt.legend()

    plt.savefig('temp/Results/'+name+'.png')
    plt.close()
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--id', default="0" , type=str)
    parser.add_argument('--coarse', default=0 , type=float)
    args = parser.parse_args()


    t5 = time.time()
    # 创建一个空字典来存储变量  
    variables = {}  
    
    cal_size = (3, 24)
    TOT_Intensity = np.zeros(cal_size)
    TOT_Diffuse = np.zeros(cal_size)
    TOT_Specular = np.zeros(cal_size)
    # # 定义一个字符串作为变量名  
    # variable_name = "my_variable"  
    
    # # 在字典中为变量赋值  
    # variables[variable_name] = 42  
    
    # # 现在，可以通过字典访问该变量  
    # print(variables[variable_name])  # 输出: 42
    Coarse_list = np.linspace(0 , np.pi/3 , cal_size[0])
    Theta_list = np.linspace(0, 2*np.pi - np.pi/12, cal_size[1])
    Coarse = args.coarse

    for j, Theta in enumerate(Theta_list):
        t1 = time.time()

        var_name = "Intensity_" + str(int(Coarse*180/np.pi)) + "_" + str(int(Theta*180/np.pi))
        var_name2 = "Diffuse_ratio_" + str(int(Coarse*180/np.pi)) + "_" + str(int(Theta*180/np.pi))

        variables[var_name], variables[var_name2] = mf.global_intensity(Theta, Coarse = Coarse)
        TOT_Intensity[i,j] = variables[var_name].sum()
        TOT_Diffuse[i,j] = (variables[var_name2] * variables[var_name]).sum()
        TOT_Specular[i,j] = TOT_Intensity[i,j] - TOT_Diffuse[i,j]

        t2 = time.time()
        print("Coarse = ", Coarse, "Theta = ", Theta, "Time = ", t2 - t1, "s, Processing Done!")

    t6 = time.time()
    print("Total Time = ", t6 - t5, "s, Processing ALL DONE!")

    # save "variables" to temp/ folder
    np.save('temp/variables/variables.npy', variables)
    np.save('temp/variables/TOT_Intensity.npy', TOT_Intensity)
    np.save('temp/variables/TOT_Diffuse.npy', TOT_Diffuse)
    np.save('temp/variables/TOT_Specular.npy', TOT_Specular)
    np.save('temp/variables/Coarse_list.npy', Coarse_list)
    np.save('temp/variables/Theta_list.npy', Theta_list)

    # plot the results
    result_ploter(TOT_Intensity, "TOT_Intensity", cal_size, Theta_list, Coarse_list)
    result_ploter(TOT_Diffuse, "TOT_Diffuse", cal_size, Theta_list, Coarse_list)
    result_ploter(TOT_Specular, "TOT_Specular", cal_size, Theta_list, Coarse_list)

