# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 02:23:09 2024

@author: dell
"""
import json
import numpy as np

rho_J = 1.326
R_J = 69911*1000
R_E = 6357*1000
M_J = 1.8982 * 1e27
G = 6.67 *1e-11
CG = G * M_J/(4*np.pi**2)
h = 6.626e-34  # Planck's constant
c = 3.0e8  # Speed of light
k = 1.38e-23  # Boltzmann constant
A = h * c/ k
lam = 0.82e-6

def B(lam,T):
    return 1/(np.exp(A/lam/T)-1)

def comp(x):
    Teq = x['Teq']
    Tday = x['Tday']
    Ts = x['Teff']
    Rp = x['Rp'] * R_J
    Per = x['Period'] * 24*3600
    # gs = 10**x['log10g_s']
    # TR = (Teq/Ts)**4
    a = (Ts/Teq)**2
    
    return Rp / a  
# 读取 JSON 数据
with open(r'C:\Users\dell\Desktop\Exo\exoplanet_table.json', 'r') as file:
    exoplanets = json.load(file)
    
exo = exoplanets['data']    
exo = [planet for planet in exo if  type(planet['Rp']) is float and type(planet['Period']) is float and type(planet['Tday']) is float and type(planet['log10g_s']) is float]
# 按 R_p/a 排序
# exo = [planet for planet in exo if rho_J * planet['Mp']/planet['Rp']**3 > 3.0 ]
exo = [planet for planet in exo if planet['Rp']* R_J < 1.6 * R_E ]

sorted_exo = sorted(exo, key= comp, reverse=True)

# 输出前十个对象
top_exoplanets = sorted_exo[:10]

# 打印结果
for planet in top_exoplanets:
    print(planet['Planet_Name'],' R_p/a ', ' T:',planet['Tday'],' K-band SNR:',planet['SNR_Transmission_K_mag'])

