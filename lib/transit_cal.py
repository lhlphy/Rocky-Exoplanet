import numpy as np
from parameter_list import PPs, APs
from function_library import Cal_star_flux, B
import os

# 将程序结果进行修正，并考虑transit
name = 'R9'
I_specular = np.load(f'temp/{name}/variables/I_specular.npy')
I_intensity = np.load(f'temp/{name}/variables/I_intensity.npy')
I_diffuse = np.load(f'temp/{name}/variables/I_diffuse.npy')
Theta_list = np.load(f'temp/{name}/variables/Theta.npy')
Wave_list = np.load(f'temp/{name}/variables/wave_list.npy')
Thermal = np.load(f'temp/{name}/variables/Thermal.npy')
Star_flux = np.load(f'temp/{name}/variables/Star_flux.npy')

SF = np.zeros(I_specular.shape)

for i, theta in enumerate(Theta_list):
    for j, wave in enumerate(Wave_list):
        SF[j, i] = Cal_star_flux(theta, wave)
        
SF_full = SF[:,np.size(SF,1)//2]  # only consider the half of the star
SF_full = np.array([SF_full]).transpose()
# print(SF_full)

IS = I_specular * Star_flux 
ID = I_diffuse * Star_flux
II = I_intensity * Star_flux
Thermal = Thermal * Star_flux

# 使用 np.nan_to_num 将 NaN 替换为 0  
IS = np.nan_to_num(IS, nan=0.0)  
ID = np.nan_to_num(ID, nan=0.0)
II = np.nan_to_num(II, nan=0.0)
Thermal = np.nan_to_num(Thermal, nan=0.0)

IS2 = IS / SF_full
ID2 = ID / SF_full
II2 = II  / SF_full
Thermal = (Thermal + SF - SF_full) / SF_full

os.makedirs(f'temp/{name}copy/variables', exist_ok=True)
os.makedirs(f'temp/{name}copy/Results', exist_ok=True)

np.save(f'temp/{name}copy/variables/I_specular.npy', IS2)
np.save(f'temp/{name}copy/variables/I_diffuse.npy', ID2)
np.save(f'temp/{name}copy/variables/I_intensity.npy', II2)
np.save(f'temp/{name}copy/variables/Theta.npy', Theta_list)
np.save(f'temp/{name}copy/variables/wave_list.npy', Wave_list)
np.save(f'temp/{name}copy/variables/Star_flux.npy', SF)
np.save(f'temp/{name}copy/variables/Thermal.npy', Thermal)

