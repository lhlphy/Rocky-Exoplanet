import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os

def IDS_plot(name, Obs_wave):
    # load data
    I_specular = np.load(f'temp/{name}/variables/I_specular.npy')
    I_intensity = np.load(f'temp/{name}/variables/I_intensity.npy')
    I_diffuse = np.load(f'temp/{name}/variables/I_diffuse.npy')
    Theta_list = np.load(f'temp/{name}/variables/Theta.npy')
    Wave_list = np.load(f'temp/{name}/variables/wave_list.npy')
    Thermal = np.load(f'temp/{name}/variables/Thermal.npy')

    if type(Obs_wave) == float:
        Obs_wave = np.array([Obs_wave])

    # interpolate
    spls = interp1d(Wave_list, I_specular.T, kind='cubic')
    spli = interp1d(Wave_list, I_intensity.T, kind='cubic')
    spld = interp1d(Wave_list, I_diffuse.T, kind='cubic')
    splt = interp1d(Wave_list, Thermal, kind='cubic')
    I_s = spls(Obs_wave).T
    I_i = spli(Obs_wave).T
    I_d = spld(Obs_wave).T
    I_t = splt(Obs_wave).T

    theta = np.linspace(0, 2*np.pi, 200)
    spls = interp1d(Theta_list, I_s, kind='cubic')
    spli = interp1d(Theta_list, I_i, kind='cubic')
    spld = interp1d(Theta_list, I_d, kind='cubic')
    splt = interp1d(Theta_list, I_t, kind='cubic')
    I_s = spls(theta)
    I_i = spli(theta)
    I_d = spld(theta)
    I_t = splt(theta)

    # plot
    plt.rcParams['font.size'] = 12
    plt.figure(figsize=(8, 5))
    plt.plot(theta, I_s * 1e6, label='Specular')
    plt.plot(theta, I_d * 1e6, label='Diffuse')
    plt.plot(theta, I_t * 1e6, label='Thermal')

    plt.xlabel('Orbital angle (rad)')
    plt.ylabel('Contrast ratio (ppm)')
    plt.show()
    os.makedirs(f'temp/P0', exist_ok= True)
    plt.savefig(f'temp/P0/compare.png')
    plt.close()



if __name__ =='main':
    IDS_plot('LHS 3884 re', np.array([0.5]) * 1e-6)