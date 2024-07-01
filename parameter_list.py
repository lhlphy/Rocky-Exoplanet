import numpy as np
# Constants List

### orbital parameters
AU = 149_597_870.7  # km, 1 Astronomical Unit 149597870.7
R1 = 696300 # km, radius of the Star
R2 = 6378      # km, radius of the Planet
e = 0  # Eccentricity of the Earth's orbit
a = AU  # km, semi-major axis of the Earth's orbit / 10
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
SIZE = [90, 180]  # Size of the meshgrid
# Create meshgrid for the planet
phiP_list = np.linspace(-np.pi / 2, np.pi / 2, SIZE[0])
thetaP_list = np.linspace(0, 2 * np.pi, SIZE[1])


### Thermal and optical parameters
Temperature = 5800  # K, temperature of the Star
Wavelengh = 6e-7  # m, wavelength of the light
SPE_REF_g = 0.5  # Specular reflection coefficient  #if SPE_REF == -1, using the experiment data
DIF_REF_g = 0.5  # Diffuse reflection coefficient
Coarse_g = 0.5  # Coarseness of the surface



