import numpy as np
from parameter_list import PPs, APs
import main_function as mf
from function_library import B, sym_complete, decorator_timer, reflection_plot
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os
 
 
def Fresnel_coefficient(Lam, Mtheta=-1, Mphi=-1, Polarization="default"):
    """
    calculate the reflection coefficient for specular reflector
    OpticalFrame: 
        "Demo": no input parameter needed
        "Full_cal": Lam is the only needed parameter
    """
    if Polarization == "default":
        Polarization = APs.Polarization

    # calculate the reflection coefficent in different wavelength and different location
    if Mtheta == -1 and Mphi == -1:
        Mtheta, Mphi = np.meshgrid(APs.thetaP_list, APs.phiP_list)
        
    # calculate the cosine of the angle between the normal vector and the incident ray
    COSI = np.cos(Mtheta) * np.cos(Mphi)
    SINI = np.sqrt(1 - COSI**2)
        
    if APs.OpticalFrame == "Demo":
        An = PPs.std_FR
    elif APs.OpticalFrame == "Full_cal":
        An = PPs.Albedo(Lam)
    elif APs.OpticalFrame == "Non_Fresnel":  # Non_Fresnel : Do not use Fresnel equation, just return Albedo
        An = PPs.Albedo(Lam)
        return An
    else:
        raise ValueError("Version ERROR: OpticalFrame must be either 'Demo', 'Full_cal' or 'Non_Fresnel', other options are under construction")
        
    n = 2/(1- np.sqrt(An)) -1 
    Co1 = np.sqrt(n**2 - SINI**2)
    Rs = ((COSI - Co1) / (COSI + Co1)) **2
    Rp = ((Co1 - n**2 *COSI)/ (Co1 + n**2 *COSI))**2
    R = (Rs + Rp) / 2
    if Polarization == "None":  # natural light, no polarization
        return R
    elif Polarization == "S":  # Observer is S polarization 
        Rs_O = (Rp *np.sin(Mphi)**2 + Rs * np.cos(Mphi)**2 )/2
        return Rs_O
    elif Polarization == "P": # Observer is P polarization
        Rp_O = (Rp *np.cos(Mphi)**2 + Rs * np.sin(Mphi)**2 )/2
        return Rp_O
    else:
        raise ValueError("Polarization must be either 'S', 'P', or 'None'")
        
    
@decorator_timer('Full_spectrum')
def Full_spectrum(wavelength_bound, args = None, id = 0, Ntheta = 5, Nwave = 1):
    """
    calculate thermal spectrum and reflection spectrum
    """

    if args != None:
        id = args.id
        Ntheta = args.Ntheta
        Nwave = args.Nwave
    
    # Create the folder to store the results
    os.makedirs(f'temp/R{id}/variables', exist_ok=True)
    os.makedirs(f'temp/R{id}/Results', exist_ok=True)

    """ Calculate the thermal spectrum """ 
    mf.thermal_spectrum(wavelength_bound, id= id, Ntheta = Ntheta, NWavelength= Nwave, Nsubpro= args.Nsubpro)
    # Load the results
    thermal_ratio = np.load(f'temp/R{id}/variables/Thermal.npy')
    Theta_list  = np.load(f'temp/R{id}/variables/Theta.npy')
    Star_flux  = np.load(f'temp/R{id}/variables/Star_flux.npy')


    """ Calculate the reflected and diffused light   """
    if Ntheta == 1:  # Create Theta_list
        Theta_list2 = np.array([np.pi])
    else:
        Theta_list2 = np.linspace(0, np.pi, Ntheta)  # 0-pi 与 pi-2pi 重复
    Wave_list = np.linspace(wavelength_bound[0], wavelength_bound[1], Nwave)

    # 与 thermal_ratio 一致
    I_intensity = np.zeros([Nwave, Ntheta])
    I_diffuse = I_intensity.copy()
    I_specular = I_intensity.copy()

    # Create vector to store the geometry results, and using the gemotry result to generate the intensity I_xxxxxx
    Tmap = np.load(f'temp/R{id}/plots/Tmap0.npy')

    # calculate the reflection coefficent in different wavelength and different location
    # different wavelength and different surface temperature result in different albedo
    # consider_temperature = False   # this key determined wheather consider the influence of T, False will save the time #######
    # if consider_temperature:  # consider the influence of Temperature
    #     # need a matrix to save the albedo in each location and each wavelength
    #     AD_matrix = np.zeros([Nwave, APs.SIZE[0], APs.SIZE[1]])
    #     AS_matrix = np.zeros([Nwave, APs.SIZE[0], APs.SIZE[1]])
        
    #     for j, wave in enumerate(Wave_list): 
    #         for n, phi in enumerate(APs.phiP_list):
    #             for m, theta in enumerate(APs.thetaP_list):
    #                 AD_matrix[j, n, m] = PPs.A_diffuse( wave, Tmap[n, m])
    #                 AS_matrix[j, n, m] = PPs.A_Specular( wave, Tmap[n, m])
    #         # 补上planck函数, global_intensity给出的结果不包含普朗克函数和反射率
    #         AD_matrix[j,:,:] = AD_matrix[j,:,:] * B(wave, PPs.Stellar_T)   
    #         AS_matrix[j,:,:] = AS_matrix[j,:,:] * B(wave, PPs.Stellar_T)
            
    # else:   # do not consider the influence of Temperature
    #     # only need to save albedo in each wavelength
    #     AD_matrix = np.zeros(Nwave)
    #     AS_matrix = np.zeros(Nwave)
    #     for j, wave in enumerate(Wave_list):
    #         AD_matrix[j] = PPs.A_diffuse(wave)
    #         AS_matrix[j] = PPs.A_Specular(wave)
    #         # 补上planck函数, global_intensity给出的结果不包含普朗克函数和反射率
    #         AD_matrix[j] = AD_matrix[j] * B(wave, PPs.Stellar_T)   
    #         AS_matrix[j] = AS_matrix[j] * B(wave, PPs.Stellar_T)
            
    # calculate the reflection using global_intensity
    for i, Theta in enumerate(Theta_list2):
        I, D, S = mf.global_intensity(Theta, id, Model= APs.Model)
        # I: total intensity; D:diffuse ; S:specular reflection
        # both under the condition of albedo = 1
        if APs.Model == 'SD_combined':
            TMAP0 = np.load(f'temp/R{id}/plots/Tmap0.npy')
            D[TMAP0 > PPs.T_liq] = 0
            S[TMAP0 < PPs.T_liq] = 0
            I = D + S

        for j, wave in enumerate(Wave_list):
            # D1 = D * AD_matrix[j]
            # S1 = S * AS_matrix[j]
            # consider Fresnel_coefficient
            Fresnel_matrix = Fresnel_coefficient(wave)* B(wave, PPs.Stellar_T)
            D1 = D * Fresnel_matrix  # PPs.std_FR * B(wave, PPs.Stellar_T)
            S1 = S * Fresnel_matrix

            # consider the reflection coefficent, D1,S1 is diffuse and specular reflection map in different location of planet
            # I1 = D1 + S1  
            I_diffuse[j, i] = D1.sum()
            I_specular[j, i] = S1.sum()
            I_intensity[j, i] = I_diffuse[j, i] + I_specular[j, i]
            # sum up intensity in different location get the total intensity
        reflection_plot(S1, Theta, id)  # plot the specular reflection map

    # symmetry
    I_intensity = sym_complete(I_intensity,1) / Star_flux
    I_diffuse = sym_complete(I_diffuse,1) / Star_flux
    I_specular = sym_complete(I_specular,1) / Star_flux
    Theta_list2 = sym_complete(Theta_list2)
    Theta_list2[Ntheta: 2*Ntheta] = 2*np.pi - Theta_list2[Ntheta: 2*Ntheta]

    if (Theta_list != Theta_list2).any():  # Check the consistency of the Theta_list [size: Ntheta*2]
        print(Theta_list)
        print(Theta_list2)
        raise ValueError("Theta_list is not consistent!")
    
    # Save the results
    np.save(f'temp/R{id}/variables/I_intensity.npy', I_intensity)
    np.save(f'temp/R{id}/variables/I_diffuse.npy', I_diffuse)
    np.save(f'temp/R{id}/variables/I_specular.npy', I_specular)
    np.save(f'temp/R{id}/variables/wave_list.npy', Wave_list)
    #np.save(f'temp/R{id}/variables/Theta_list.npy', Theta_list2)


    """Combine the thermal and optical results, generate the full spectrum"""

    FS = I_intensity + thermal_ratio   # Full Spectrum  [size: Nwave * (2*Ntheta)]  xaixs: Wave_list, yaixs: Theta_list

    FS_plotter(FS, Wave_list, Theta_list, Nwave, Ntheta, id)


def FS_plotter(FS, Wave_list, Theta_list, Nwave, Ntheta, id = 0 , Obs_wavelength = APs.Obs_array):
    """
    Plot the Contrast_ratio and Full_spectrum and phase_curve
    FS: Full Spectrum (matrix: Nwave* (2*Ntheta))
    Obs_wavelength: wavelength be observed (array)
    """
    plt.rcParams['font.size'] = 12
    ## plot Full spectrum
    fig, ax = plt.subplots(figsize=(8,4))
    #Xwave = np.linspace(Wave_list[0], Wave_list[-1], 1000)
    ## plot Contrast ratio
    for i in range(Ntheta):
        #spl = interp1d(Wave_list, FS[:,i], kind='cubic')
        ax.plot(Wave_list* 1e6, FS[:,i] *1e6, label = f'{round(Theta_list[i],3)}')
    ax.set_xlabel('Wavelength ($\mu$m)')
    ax.set_ylabel('Contrast ratio (ppm)')
    ax.set_title(f'Full Spectrum')
    legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    legend.set_title("Orbital Phase")
    fig.subplots_adjust(right=0.8)
    plt.show() 
    plt.savefig(f'temp/R{id}/Results/Constrat_ratio.png')
    plt.close()

    # plot Full Spectrum
    fig, ax = plt.subplots(figsize = (8, 5))
    Star_flux  = np.load(f'temp/R{id}/variables/Star_flux.npy')
    for i in range(Ntheta):
        #spl = interp1d(Wave_list, FS[:,i], kind='cubic')
        ax.plot(Wave_list* 1e6, FS[:,i]* Star_flux[:,i] , label = f'{round(Theta_list[i],3)}')
    ax.set_xlabel('Wavelength ($\mu$m)')
    ax.set_ylabel('$\mathrm{Spectrum\; of\; Planet \; (W \cdot sr^{-1}\cdot nm^{-1})}$')
    ax.set_title(f'Full Spectrum')
    legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    legend.set_title("Orbital Phase")
    fig.subplots_adjust(right=0.8)
    plt.show()
    plt.savefig(f'temp/R{id}/Results/Full_spectrum.png')
    plt.close()

    ## plot phase curve
    if type(Obs_wavelength) != np.ndarray:
        Obs_wavelength = np.array([Obs_wavelength])

    spl = interp1d(Wave_list, FS.T, kind = 'cubic')
    y = spl(Obs_wavelength)
    y = y.T

    fig , ax = plt.subplots(figsize=(8, 4))
    for i, Owave in enumerate(Obs_wavelength):
        ax.plot(Theta_list, y[i,:] *1e6, label = f'{round(Owave* 1e6,1)} $\mu$m')
    ax.set_xlabel('Orbital phase')
    ax.set_ylabel('Contrast ratio (ppm)')
    ax.set_title(f'Phase curve')
    legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    legend.set_title("Wavelength")
    fig.subplots_adjust(right=0.8)
    plt.show()
    plt.savefig(f'temp/R{id}/Results/Phase_curve.png')
    plt.close()




    






