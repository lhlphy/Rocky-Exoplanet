import numpy as np
import pandas as pd 
from lava_data import LA
import os
from scipy.interpolate import interp1d

# Constants List
AU = 149_597_870.7  # km, 1 Astronomical Unit 149597870.7
Sigma_const = 5.67e-8  # W/m^2/K^4, Stefan-Boltzmann constant

# parser = argparse.ArgumentParser()
print('parameter_list.py')
# if environment variable 'roughness' is not set, use default value 0
if 'roughness' not in os.environ:
    os.environ['roughness'] = '0'
if 'FRnormal' not in os.environ:
    os.environ['FRnormal'] = '0'
if 'lavatype' not in os.environ:
    os.environ['lavatype'] = 'zero'
class Accuracy_parameters:
    ### Accuracy control parameters
    def __init__(self):
        self.SIZE = [37, 74]  # Size of the meshgrid
        # Create meshgrid for the planet
        self.phiP_list = np.linspace(-np.pi / 2, np.pi / 2, self.SIZE[0])
        self.thetaP_list = np.linspace(0, 2 * np.pi, self.SIZE[1])
        self.Obs_array = np.array([3]) * 1e-6  # The wavelength of the observation array
        self.Polarization = os.getenv('Polarization', default="None")
        self.OpticalFrame = os.getenv('OpticalFrame', default="Demo")

        #####################################################################################################
        if os.getenv('Model') != 'None':
            self.Model = os.getenv('Model')
        else:
            if float(os.getenv('roughness')) < 1e-3: # if roughness = 0, use the Specular_Only model
                self.Model =  'Specular_Only'
            elif float(os.getenv('roughness')) >= 1e-3 and float(os.getenv('roughness')) <= 90 :
                self.Model =  'Gaussian_wave' 
            else:
                self.Model = 'Lambert_Only'
            # self.Model = 'Lambert'
        ######################################################################################################
        Mode = os.getenv('mode') # get mode from environment variable
        if Mode == 'PC':
            self.mode = 'Phase curve'   # Phase curve & Transit
            # main_function.py -> BRDF()->line 40: controled by this parameter (whether consider the block effect of star)
        elif Mode == "TR":
            self.mode = 'Transit'
        else:
            self.mode = Mode

APs = Accuracy_parameters()

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
        if APs.OpticalFrame == 'Demo':
            self.std_FR = float(os.getenv('FRnormal'))  # Fresnel reflection coefficient when normal incidence
        elif APs.OpticalFrame == 'Full_cal':
            def std_FR(self, Lam):
                return self.Albedo(Lam)
        else:
            raise ValueError("Version ERROR: OpticalFrame must be either 'Demo' or 'Full_cal', other options are under construction")
        
        ### Observation parameters
        camera = np.array([1, 0, 0]) # Define the direction of the camera (observation vector)
        self.camera = camera / np.linalg.norm(camera) # normalize the camera vector
        self.T_liq = 1201.16  # K, liquid temperature, 作为完全融化区域的温度阈值
        
        self.An_list = np.linspace(0, 1, 100)
        self.A_Mean_list = self.A_Mean_list_Cal(self.An_list)
        self.A_Mean_interp = interp1d(self.An_list, self.A_Mean_list, kind = 'linear')
       
            
    def Albedo(self, lam, T = 0):
        if  lam > LA.Wmax *1e-6:
            raise ValueError('Albedo exceed the high bounds! ')
        return  LA.A_interp(lam*1e6)
    
    def A_diffuse(self, lam, T = 0):
        return self.Albedo(lam, T)
    
    def A_Specular(self, lam, T = 0):
        return self.Albedo(lam, T)
    
    def Fresnel(self, Lam = -1, I_angle = 0, Polarization = 'default', A_normal = 0):
        if Polarization == 'default':  # get the polarization from the environment variable
            Polarization = APs.Polarization
            
        SINI = np.sin(I_angle)
        COSI = np.cos(I_angle)
        if Lam < 0:
            A = A_normal
        else:
            A = self.Albedo(Lam)
            
        n = 2/(1- np.sqrt(A)) -1
        Co1 = np.sqrt(n**2 - SINI**2)
        if Polarization == 'S':
            Rs = ((COSI - Co1) / (COSI + Co1)) **2
            return Rs
        elif Polarization == 'P':
            Rp = ((Co1 - n**2 *COSI)/ (Co1 + n**2 *COSI))**2
            return Rp
        elif Polarization == 'None':

            Rs = ((COSI - Co1) / (COSI + Co1)) **2
            Rp = ((Co1 - n**2 *COSI)/ (Co1 + n**2 *COSI))**2
            return (Rs+Rp)/2
        else:
            raise ValueError("Polarization must be either 'S', 'P', or 'None'")
        
    def A_Mean_list_Cal(self, An_list):  # calculate the mean albedo for each An
        from scipy.integrate import quad
        A_Mean_list = np.zeros(len(An_list))
        Polarization = APs.Polarization
        
        for i, An in enumerate(An_list):
            A_Mean_list[i] = quad(lambda theta: self.Fresnel(-1, theta, Polarization, A_normal= An) *np.sin(theta), 0, np.pi/2)[0]
            
        return A_Mean_list
        
    def A_Mean(self, An = -1, Wavelength = -1):
        if (Wavelength == -1) and (An == -1):
            raise ValueError('lam or An must be specified')
        elif (Wavelength != -1) and (An != -1):
            raise ValueError('lam and An cannot be both specified')
        
        if Wavelength > 0:
            A = self.Albedo(Wavelength)
        elif An > 1e-6:
            A = An
        elif An < 1e-6:
            return 0
        else:
            raise ValueError('A_Mean: lam or An must have a proper value')
        
        return self.A_Mean_interp(A)
        

PPs = Planet_parameters(4170) 
# K2-141 b : 4264 - 98  /4170
# 55 Cnc e : 215 - 98
# TOI-2445 b: 34287 - 98
# GJ-367 b: 733 - 98
# Kepler-808 b: 30420 ## this planet has best specular/thermal*Melt_area ratio in 1 um (from exo_rank.py)
# 10781  Kepler-1320 b


if __name__ == '__main__':
    os.environ['roughness'] = '0'
    print("This is a module that contains all the parameters.")
    print('Rs:', PPs.Rs)
    print('Rp:', PPs.Rp)
    print('semi_axis:', PPs.semi_axis)
    print('Stellar_T:', PPs.Stellar_T)
    print('pl_eqT:', PPs.pl_eqT)



