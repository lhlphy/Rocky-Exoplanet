import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.integrate import quad
import pandas as pd 
# import matplotlib.pyplot as plt
# Constants List

### orbital parameters
AU = 149_597_870.7  # km, 1 Astronomical Unit 149597870.7
R1 = 695_500 * 0.19* 0.1 # km, radius of the Star
R2 = 69_911 * 0.116      # km, radius of the Planet
e = 0  # Eccentricity of the Earth's orbit
a = AU * 0.006 # km, semi-major axis of the Earth's orbit / 10 
# m1 = 1.989e30   # kg, mass of the Sun
# m2 = 5.972e24   # kg, mass of the Earth
# m1 = m1 / m2    # Normalized mass of the Sun
# m2 = 1          # Normalized mass of the Earth
#Theta = np.pi/3 # The orbit angle (Old version)

### Observation parameters
# Define the direction of the camera (observation vector)
camera = np.array([1, 0, 0])
#normalize the camera vector
camera = camera / np.linalg.norm(camera)


### Accuracy control parameters
SIZE = [720, 1440]  # Size of the meshgrid
# Create meshgrid for the planet
phiP_list = np.linspace(-np.pi / 2, np.pi / 2, SIZE[0])
thetaP_list = np.linspace(0, 2 * np.pi, SIZE[1])


### Thermal and optical parameters
Temperature = 3036  # K, temperature of the Star
Wavelength = 5e-6  # m, wavelength of the light  ### Attention: wavelength and wavelengh
# SPE_REF_g = 555  # Specular reflection coefficient  1.if SPE_REF <0 , using the experiment data   2.if SPE_REF >1, using the Fresnel equation model
# DIF_REF_g = 0.1 # Diffuse reflection coefficient
Coarse_g = 0  # Coarseness of the surface
# Albedo = 0  # Albedo of the Earth
Sigma = 5.67e-8  # W/m^2/K^4, Stefan-Boltzmann constant

Wind_speed = 10   # wind speed in m/s (only available for the Gaussian wave model)
Obs_array = np.array([6e-7, 8e-7, 1e-6, 1.5*1e-6, 2*1e-6, 3*1e-6, 5*1e-6])

# 读取CSV文件  
file_path = 'Teide.csv'  # 请确保将路径更改为正确的文件路径  
data = pd.read_csv(file_path) 
# data = data.dropna()  # 删除所有包含NaN的行
Data = np.zeros(data.shape)
Data[:,0] = pd.to_numeric(data['1187 K'], errors='coerce')  
Data[:,1] = pd.to_numeric(data['Unnamed: 9'], errors='coerce')
Data[:,2] = pd.to_numeric(data['1413 K'], errors='coerce')
Data[:,3] = pd.to_numeric(data['Unnamed: 1'], errors='coerce')
Data[:,4] = pd.to_numeric(data['1673 K'], errors='coerce')
Data[:,5] = pd.to_numeric(data['Unnamed: 3'], errors='coerce')
Data[:,6] = pd.to_numeric(data['1963 K'], errors='coerce')
Data[:,7] = pd.to_numeric(data['Unnamed: 5'], errors='coerce')
Data[:,8] = pd.to_numeric(data['2164 K'], errors='coerce')
Data[:,9] = pd.to_numeric(data['Unnamed: 7'], errors='coerce')
spl= [0] * (Data.shape[1]//2)
for i in range(Data.shape[1]//2):
    x = Data[~np.isnan(Data[:,2*i]),2*i]
    x = x + np.random.rand(x.size) * 1e-6
    spl[i] = interp1d(x , Data[~np.isnan(Data[:,2*i+1]),2*i+1], kind='linear', fill_value="extrapolate")

def Albedo(lam, T):
    T_list = np.array([1187, 1413, 1673, 1963, 2164])

    E1 = np.zeros(len(T_list))
    for i in range(len(T_list)):
        if lam * 1e6 > 16:
            E1[i] = spl[i](16)
        else:
            E1[i] = spl[i](lam * 1e6)

    spl2 = interp1d(T_list, E1, kind='linear', fill_value="extrapolate")
    E2 = spl2(T)
    return  1 - E2

def A_Specular(lam, T):
    # data_lam = np.linspace(0.1, 10) *1e-6
    # data_A = 0.5 * np.ones(data_lam.size)
    # spl = interp1d(data_lam, data_A, kind='cubic')
    # return spl(lam)
    return 0

def A_diffuse(lam, T):
    # data_lam = np.linspace(0.1, 10) *1e-6
    # data_A = 0.5 * np.ones(data_lam.size)
    # spl = interp1d(data_lam, data_A, kind='cubic')
    # return spl(lam)
    return 1 - Albedo(lam, T)


# w = np.linspace(1, 10,100) * 1e-6
# A = np.zeros(w.size)
# for i, wave in enumerate(w):
#     A[i] = Albedo(wave, 00)
# plt.plot(w, A)
# plt.show()

