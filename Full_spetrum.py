import numpy as np
from parameter_list import *
import main_function as mf
from scipy.integrate import dblquad
from scipy import interpolate
import matplotlib.pyplot as plt
import os
import time


def Full_spetrum(wavelength_bound, args = None, Temperature = Temperature, Albedo = Albedo, id = 0, Ntheta = 5, Nwave = 1):

    if args != None:
        id = args.id
        Albedo = args.Albedo
        Ntheta = args.Ntheta
        Nwave = args.Nwave
    
    # Create the folder to store the results
    os.makedirs(f'temp/R{id}/variables', exist_ok=True)
    os.makedirs(f'temp/R{id}/Results', exist_ok=True)

    t0 = time.time()
    # Calculate the thermal spectrum
    mf.thermal_spectrum(wavelength_bound, Temperature, Albedo= Albedo, id= id, Ntheta = Ntheta, NWavelength= Nwave)
    # Load the results
    thermal_ratio = np.load(f'temp/R{id}/Results/Ratio.npy')
    Theta_list  = np.load(f'temp/R{id}/Results/Theta.npy')

    t1 = time.time()
    print("Total Time = ", t1 - t0, "s, Processing ALL DONE!")
