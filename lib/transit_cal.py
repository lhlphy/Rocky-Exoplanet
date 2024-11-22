import os
os.environ['roughness'] = '0'
import numpy as np
from parameter_list import PPs, APs
from function_library import Cal_star_flux, B

# 将程序结果进行修正，并考虑transit（原程序未考虑transit）
# 原程序计算contrast ratio 使用的是各个轨道位置的star flux, 这显然是不正确的，应该使用恒星的全光通量 Int(Pi*Rs^2*B(Ts, wave), lam1, lam2), 这是个常数
# 所以该程序通过还原绝对光通量，然后除以恒定的恒星全光通量，从而得到正确的contrast ratio
if __name__ == "__main__":
    import argparse 
    parser = argparse.ArgumentParser()
    parser.add_argument('--name', type=str, default='R6')  # read in the name of folder to be processed
    args = parser.parse_args()
    name = args.name
    
    I_specular = np.load(f'temp/{name}/variables/I_specular.npy')
    I_intensity = np.load(f'temp/{name}/variables/I_intensity.npy')
    I_diffuse = np.load(f'temp/{name}/variables/I_diffuse.npy')
    Theta_list = np.load(f'temp/{name}/variables/Theta.npy')
    Wave_list = np.load(f'temp/{name}/variables/wave_list.npy')
    Thermal = np.load(f'temp/{name}/variables/Thermal.npy')
    Star_flux = np.load(f'temp/{name}/variables/Star_flux.npy')

    SF = np.zeros(I_specular.shape) # 存储不同orbital_angle and wavelength下的恒星光通量

    for i, theta in enumerate(Theta_list):
        for j, wave in enumerate(Wave_list):
            SF[j, i] = Cal_star_flux(theta, wave)
            
    # 计算恒星全光通量
    SF_full = SF[:,np.size(SF,1)//2]  # only consider the half of the star
    SF_full = np.array([SF_full]).transpose()
    # print(SF_full)

    # 计算出行星的绝对光通量
    IS = I_specular * Star_flux 
    ID = I_diffuse * Star_flux
    II = I_intensity * Star_flux
    Thermal = Thermal * Star_flux

    # 使用 np.nan_to_num 将 NaN 替换为 0  
    IS = np.nan_to_num(IS, nan=0.0)  
    ID = np.nan_to_num(ID, nan=0.0)
    II = np.nan_to_num(II, nan=0.0)
    Thermal = np.nan_to_num(Thermal, nan=0.0)

    # 将行星的绝对光通量除以恒星全光通量，得到正确的contrast ratio
    IS2 = IS / SF_full 
    ID2 = ID / SF_full
    II2 = II  / SF_full
    Thermal = (Thermal + SF - SF_full) / SF_full
    # thermal同时包含了thermal emission 和transit的修正项， 前者为正值或0，后者为负值
    # 而且在我们考虑的earth-like planet的情况下，应该有在transit附近 thermal emission << transit修正

    os.makedirs(f'temp/{name}copy/variables', exist_ok=True)
    os.makedirs(f'temp/{name}copy/Results', exist_ok=True)

    np.save(f'temp/{name}copy/variables/I_specular.npy', IS2)
    np.save(f'temp/{name}copy/variables/I_diffuse.npy', ID2)
    np.save(f'temp/{name}copy/variables/I_intensity.npy', II2)
    np.save(f'temp/{name}copy/variables/Theta.npy', Theta_list)
    np.save(f'temp/{name}copy/variables/wave_list.npy', Wave_list)
    np.save(f'temp/{name}copy/variables/Star_flux.npy', SF)
    np.save(f'temp/{name}copy/variables/Thermal.npy', Thermal)

