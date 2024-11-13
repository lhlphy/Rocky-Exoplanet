import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def data_loader(name, Obs_wave):
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
    # print(Wave_list.shape, I_specular.T.shape)
    I_specular = np.nan_to_num(I_specular, nan=0) 
    I_intensity = np.nan_to_num(I_intensity, nan=0) 
    I_diffuse = np.nan_to_num(I_diffuse, nan=0) 
    spls = interp1d(Wave_list, I_specular.T, kind='linear', fill_value='extrapolate')
    spli = interp1d(Wave_list, I_intensity.T, kind='linear', fill_value='extrapolate')
    spld = interp1d(Wave_list, I_diffuse.T, kind='linear', fill_value='extrapolate')
    splt = interp1d(Wave_list, Thermal.T, kind='linear', fill_value='extrapolate')
    I_s = spls(Obs_wave).T
    I_i = spli(Obs_wave).T
    I_d = spld(Obs_wave).T
    I_t = splt(Obs_wave).T
    
    
    if np.size(Theta_list) == 2:
        print('Wavelength: ', Obs_wave[0]*1e6, 'um')
        print('Secondary eclipse depth difference: ',(I_d-I_s) *1e6,' ppm')
        return 

    theta = np.linspace(0, 2*np.pi, 200)
    Theta_list[np.size(Theta_list)//2] += 1e-10
    spls = interp1d(Theta_list, I_s, kind='linear', fill_value='extrapolate')
    spli = interp1d(Theta_list, I_i, kind='linear', fill_value='extrapolate')
    spld = interp1d(Theta_list, I_d, kind='linear', fill_value='extrapolate')
    splt = interp1d(Theta_list, I_t, kind='linear', fill_value='extrapolate')
    I_s = spls(theta)
    I_i = spli(theta)
    I_d = spld(theta)
    I_t = splt(theta)
    return I_s, I_i, I_d, I_t, theta
    
    
def specular_diffuse_plot(name_specular, name_diffuse, Obs_wave):
    Is_specular, Ii_specular, Id_specular, It_specular, theta = data_loader(name_specular, Obs_wave)
    Is_diffuse, Ii_diffuse, Id_diffuse, It_diffuse, theta = data_loader(name_diffuse, Obs_wave)
    
    i = 0
    theta = theta/(2 *np.pi)

    fig, ax = plt.subplots(figsize=(9,6))
    ax.plot(theta, Is_specular[i,:], label='specular')
    ax.plot(theta, Id_diffuse[i,:], label='diffuse')
    plt.legend()
    plt.savefig(f"temp/{name_specular}/specular_diffuse_{Obs_wave[0]*1e6}.png")
    plt.savefig(f"temp/{name_diffuse}/specular_diffuse_{Obs_wave[0]*1e6}.png")
    plt.close()
    
    
if __name__ == "__main__":
    specular_diffuse_plot("Spectrum_1", "Spectrum_3", np.array([3]) * 1e-6)
