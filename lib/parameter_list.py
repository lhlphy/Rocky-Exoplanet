import numpy as np
import pandas as pd 
from lava_data import LA
import os

# Constants List
AU = 149_597_870.7  # km, 1 Astronomical Unit 149597870.7
Sigma_const = 5.67e-8  # W/m^2/K^4, Stefan-Boltzmann constant

# parser = argparse.ArgumentParser()
print('parameter_list.py')
# read mode = "PC" ? "TR"
# with open('log/temp_vars.txt', 'r') as f:
#     lines = f.readlines()
#     mode = lines[0].strip()  # 去除换行符
# mode = os.getenv('mode')

class Accuracy_parameters:
    ### Accuracy control parameters
    def __init__(self):
        self.SIZE = [181, 361]  # Size of the meshgrid
        # Create meshgrid for the planet
        self.phiP_list = np.linspace(-np.pi / 2, np.pi / 2, self.SIZE[0])
        self.thetaP_list = np.linspace(0, 2 * np.pi, self.SIZE[1])
        self.Obs_array = np.array([3]) * 1e-6  # The wavelength of the observation array
        
        if float(os.getenv('roughness')) < 1e-3: # if roughness = 0, use the Specular_Only model
            self.Model =  'Specular_Only'
        else:
            self.Model =  'Gaussian_wave' 
        
        Mode = os.getenv('mode') # get mode from environment variable
        if Mode == 'PC':
            self.mode = 'Phase curve'   # Phase curve & Transit
        elif Mode == "TR":
            self.mode = 'Transit'
        else:
            self.mode = Mode


class Planet_parameters:
    def __init__(self, Nline):
        data_base = pd.read_csv('PS.csv', header = 96)
        row_data = data_base.iloc[Nline]
        ### orbital parameters
        self.Rs = row_data['st_rad'] * 696340  # km, radius of the Star
        self.Rp = row_data['pl_rade'] * 6371.4  # km, radius of the Planet
        self.eccentricity = 0 # row_data['pl_orbeccen'] # Eccentricity of the planet's orbit
        self.semi_axis = row_data['pl_orbsmax'] * AU  # km, semi-major axis of the planet's orbit
        ### Thermal and optical parameters
        self.Stellar_T = row_data['st_teff'] # K, temperature of the Star
        # self.pl_eqT = row_data['pl_eqt']  # K, fully redistribution, planet equilibrium Temperature [K] (from database)
        self.pl_eqT = self.Stellar_T * np.sqrt(self.Rs / 2 /self.semi_axis)  # from theoretical calculation
        
        self.Coarse_g = 0  # Coarseness of the surface
        self.Wind_speed = 10   # wind speed in m/s (only available for the Gaussian wave model)
        self.roughness = float(os.getenv('roughness'))
        
        ### Observation parameters
        camera = np.array([1, 0, 0]) # Define the direction of the camera (observation vector)
        self.camera = camera / np.linalg.norm(camera) # normalize the camera vector
        
    def Albedo(self, lam, T = 0):
        if lam < LA.Wmin *1e-6 or lam > LA.Wmax *1e-6:
            raise ValueError('Albedo exceed the bounds! ')
        return  LA.A_interp(lam*1e6)
    
    def A_diffuse(self, lam, T = 0):
        return self.Albedo(lam, T)
    
    def A_Specular(self, lam, T = 0):
        return self.Albedo(lam, T)
    

PPs = Planet_parameters(733 - 98)  
# K2-141 b : 4264 - 98
# 55 Cnc e : 215 - 98
# TOI-2445 b: 34287 - 98
# GJ-367 b: 733 - 98
APs = Accuracy_parameters()






