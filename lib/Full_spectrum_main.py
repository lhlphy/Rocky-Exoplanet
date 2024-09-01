import numpy as np
from lib.Full_spectrum import Full_spectrum
import argparse
import warnings
from scipy.integrate._quadpack_py import IntegrationWarning 


# # 忽略 IntegrationWarning  
# warnings.filterwarnings('ignore', category=IntegrationWarning) 

if __name__ == '__main__':
    warnings.filterwarnings('ignore', category=IntegrationWarning) 
    warnings.filterwarnings("ignore")
    parser = argparse.ArgumentParser()
    parser.add_argument('--id', default=0 , type=int)
    # parser.add_argument('--coarse', default=0 , type=float)
    # parser.add_argument('--specular', default=0.5 , type=float)
    # parser.add_argument('--diffuse', default=0.5 , type=float)
    # parser.add_argument('--Albedo', default=0.5 , type=float)
    parser.add_argument('--Ntheta', default=5 , type=int)
    parser.add_argument('--Nwave', default=1 , type=int)
    args = parser.parse_args()
    Wavelength_bound = np.array([2.77, 2.78]) *1e-6   # Wavelength range (m)

    Full_spectrum(Wavelength_bound, args)
    