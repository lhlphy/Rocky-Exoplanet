import numpy as np
# Constants List
AU = 149_597_870.7  # km, 1 Astronomical Unit
R1 = 1_392_000  # km, radius of the Star
R2 = 4880      # km, radius of the Planet
e = 0.205630   # Eccentricity of the Earth's orbit
a = 0.02 * AU  # km, semi-major axis of the Earth's orbit / 10
# m1 = 1.989e30   # kg, mass of the Sun
# m2 = 5.972e24   # kg, mass of the Earth
# m1 = m1 / m2    # Normalized mass of the Sun
# m2 = 1          # Normalized mass of the Earth

# Define the direction of the camera (observation vector)
camera = np.array([1, 0, 0])
#normalize the camera vector
camera = camera / np.linalg.norm(camera)
SIZE = [10, 30]  # Size of the meshgrid
# Create meshgrid for the planet
phiP_list = np.linspace(-np.pi / 2, np.pi / 2, SIZE[0])
thetaP_list = np.linspace(0, 2 * np.pi, SIZE[1])

#Theta = np.pi/3 # The orbit angle


