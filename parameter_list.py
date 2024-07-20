import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad
# Constants List

### orbital parameters
AU = 149_597_870.7  # km, 1 Astronomical Unit 149597870.7
R1 = 695_500 * 0.19 # km, radius of the Star
R2 = 69_911 * 0.116      # km, radius of the Planet
e = 0  # Eccentricity of the Earth's orbit
a = 0.006 * AU  # km, semi-major axis of the Earth's orbit / 10 
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
SIZE = [18, 36]  # Size of the meshgrid
# Create meshgrid for the planet
phiP_list = np.linspace(-np.pi / 2, np.pi / 2, SIZE[0])
thetaP_list = np.linspace(0, 2 * np.pi, SIZE[1])


### Thermal and optical parameters
Temperature = 3036  # K, temperature of the Star
Wavelengh = 5e-6  # m, wavelength of the light  ### Attention: wavelength and wavelengh
# SPE_REF_g = 555  # Specular reflection coefficient  1.if SPE_REF <0 , using the experiment data   2.if SPE_REF >1, using the Fresnel equation model
# DIF_REF_g = 0.1 # Diffuse reflection coefficient
Coarse_g = 0  # Coarseness of the surface
N1 = 1  # Refractive index of the Vacuum
N2 = 1.333  # Refractive index of the water
# Albedo = 0  # Albedo of the Earth
Sigma = 5.67e-8  # W/m^2/K^4, Stefan-Boltzmann constant

Wind_speed = 10   # wind speed in m/s (only available for the Gaussian wave model)

def Albedo(lam):
    # data_lam = np.linspace(0.1, 10) *1e-6
    # data_A = 0.5 * np.ones(data_lam.size)
    # spl = interp1d(data_lam, data_A, kind='cubic')
    # return spl(lam)
    return 0

def A_Specular(lam):
    # data_lam = np.linspace(0.1, 10) *1e-6
    # data_A = 0.5 * np.ones(data_lam.size)
    # spl = interp1d(data_lam, data_A, kind='cubic')
    # return spl(lam)
    return 1

def A_diffuse(lam):
    # data_lam = np.linspace(0.1, 10) *1e-6
    # data_A = 0.5 * np.ones(data_lam.size)
    # spl = interp1d(data_lam, data_A, kind='cubic')
    # return spl(lam)
    return 1


